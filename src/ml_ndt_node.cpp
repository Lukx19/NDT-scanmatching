#include <ml_ndt_node.h>

// ******************************* CONSTRUCTORS ********************

MlNdt::MlNdt(ros::NodeHandle &n, ros::NodeHandle &n_private)
    : nh_(n), nh_private_(n_private), msg_sync_(ImuSyncPolicy(10)), seq_(0),
      is_initialized(false) {

  initParameters();
  matcher_.setMaxRange(max_range_);
  matcher_.setResolution(resolution_);
  matcher_.setLayers(layers_);

  new_odom_pub_ =
      nh_.advertise<nav_msgs::Odometry>(new_odom_topic_, 100, false);

  odom_sub_.subscribe(nh_, odom_topic_, 10);
  laser_sub_.subscribe(nh_, laser_topic_, 10);

  // sync messages using approximate alghorithm
  msg_sync_.connectInput(odom_sub_, laser_sub_);
  msg_sync_.registerCallback(boost::bind(&MlNdt::data_cb, this, _1, _2));
  ROS_INFO("ML-NDT: Node is initialized.");
}

// ******************** PUBLIC FUNSTIONS ******************************

void MlNdt::start() { ros::spin(); }

// ******************** PRIVATE FUNCTIONS *****************************

void MlNdt::data_cb(const nav_msgs::Odometry::ConstPtr &odom,
                    const sensor_msgs::LaserScan::ConstPtr &laser) {
  ROS_INFO_STREAM("ML-NDT: Laser msg and odometry received.");
  sensor_msgs::PointCloud2 laser_pcl_msg;
  pcl::PointCloud<pcl::PointXYZ> laser_pcl;
  pcl::PointCloud<pcl::PointXYZ> laser_pcl_base;
  // project laser message to point cloud class in laser frame_id
  projector_.projectLaser(*laser, laser_pcl_msg);
  pcl::fromROSMsg(laser_pcl_msg, laser_pcl);
  // transform point cloud from laser frame_id -> robot base frame
  //TODO: change to time from laser msg. Resolve problems with tf timming error
  tf::StampedTransform trans;
  try {
    tf_list_.waitForTransform(robot_base_frame_, laser->header.frame_id,
                              ros::Time(0), ros::Duration(10.0));
    tf_list_.lookupTransform(robot_base_frame_, laser->header.frame_id,
                             ros::Time(0), trans);
    pcl_ros::transformPointCloud(laser_pcl, laser_pcl_base,trans);
  } catch (tf::TransformException &e) {
    ROS_ERROR_STREAM("ML_NDT: error in transforming laser point cloud from "
                     "laser_frame to robot base "
                     << e.what());
  }

  ROS_INFO_STREAM(laser_pcl_base.size());

  if(laser_pcl_base.size() == 0)
    return;

  // transform robot odometry too odometry frame
  tf::Pose pose_tf;
  Eigen::Vector3d pose;
  tf_list_.waitForTransform(odom_frame_, odom->header.frame_id,
                            odom->header.stamp, ros::Duration(5.0));

  tf::poseMsgToTF(odom->pose.pose, pose_tf);
  pose << pose_tf.getOrigin().getX(), pose_tf.getOrigin().getY(),
      tf::getYaw(odom->pose.pose.orientation);
  ROS_INFO("ML-NDT: Messages transformed.");


  pose_t calc_trans;
  if(mode_ == MLNDT2D){
    if(!getTransfMlNdt(calc_trans,pose,laser_pcl_base)){
      return;
    }
  }else if(mode_ == PCL2D){
    if(!getTransfPclNdt(calc_trans,pose,laser_pcl_base))
      return;
  }
  ROS_INFO_STREAM("NDT res:"<<calc_trans.transpose());

    // tf_broadcast_.sendTransform(tf::StampedTransform(
    //     transform, ros::Time::now(), new_odom_frame_, odom_frame_));
    ++seq_;
    old_odom_ = pose;
    //ROS_INFO("ML-NDT: New odom sent");
}

void MlNdt::initParameters() {
  // find tf_prefix if exists add it to tf_prefix_ variable
  tf_prefix_ = "";
  std::string tf_prefix_path;
  if (nh_private_.searchParam("tf_prefix", tf_prefix_path)) {
    nh_private_.getParam(tf_prefix_path, tf_prefix_);
  }

  robot_base_frame_ =
      nh_private_.param<std::string>("robot_base_frame_id", "base_link");

  odom_frame_ = nh_private_.param<std::string>("odom_farme_id", "odom");

  new_odom_frame_ =
      nh_private_.param<std::string>("new_odom_farme_id", "odom_clac");

  odom_topic_ = nh_private_.param<std::string>("odom_topic", "/odom");

  new_odom_topic_ =
      nh_private_.param<std::string>("calculated_odom_topic", "/odom_calc");

  laser_topic_ = nh_private_.param<std::string>("laser_topic", "/laser");

  resolution_ =
      static_cast<size_t>(nh_private_.param<int>("resolution", 8));

  max_range_ =
      static_cast<float>(nh_private_.param<double>("maximal_laser_range", 4.0));

  layers_ = static_cast<size_t>(nh_private_.param<int>("number_of_layers", 4));

  std::string mode = nh_private_.param<std::string>("scanmatcher_mode", "MLNDT2D");
  if(mode == "MLNDT2D"){
    mode_ = MLNDT2D;
  }else if( mode == "PCL2D"){
    mode_ = PCL2D;
  }

  if (tf_prefix_ != "") {
    robot_base_frame_ = tf_prefix_ + "/" + robot_base_frame_;
    odom_frame_ = tf_prefix_ + "/" + odom_frame_;
    new_odom_frame_ = tf_prefix_ + "/" + new_odom_frame_;
  }
}

bool MlNdt::getTransfMlNdt(pose_t & transform,
                           const pose_t & odom_pose,
                           const pcl_t & points)
{
  if (!is_initialized) {
    matcher_.initialize(odom_pose, points);
    is_initialized = true;
    ROS_INFO("ML-NDT: First laser scan inserted to structure");
    transform = pose_t::Zero();
    return false;
  } else {
    // calculate new odometry based on scan matching
    matcher_.calculate(odom_pose, points);
    transform = matcher_.getTransformation();
    ROS_INFO("ML-NDT: New odom calculated %f  %f  %f", transform[0], transform[1],
             transform[2]);
  }
    return true;
}

bool MlNdt::getTransfPclNdt(pose_t & transform,
                          const pose_t & odom_pose,
                          const pcl_t & points)
{
  if(!is_initialized){
      old_scan_ = points.makeShared();
      pcl_matcher_.setInputTarget (old_scan_);
      is_initialized = true;
      return false;
    }else{

      eigt::transform2d_t<double> odom_trans =
                                  eigt::transBtwPoses<double>(old_odom_,odom_pose);
      double distance =  eigt::getDisplacement(odom_trans);
      double angle = eigt::getAngle(odom_trans);
      if(distance < 0.3 || angle <0.2){
         return false;
       }
      pcl_t::Ptr new_pcl = points.makeShared();
      pcl_matcher_.setInputSource (new_pcl);
      // Set initial alignment estimate found using robot odometry.
        Eigen::Matrix<double, 4, 4>init_guess = eigt::convertTransform<double>(odom_trans); 
       // Calculating required rigid transform to align the input cloud to the target cloud.
       pcl_t::Ptr output_cloud (new pcl_t());
       pcl_matcher_.align (*output_cloud,init_guess.cast<float>());

       ROS_INFO_STREAM("Normal Distributions Transform has converged:" << pcl_matcher_.hasConverged ()
            << " score: " << pcl_matcher_.getFitnessScore ());
       transform = eigt::getPoseFromTransform<double>(
                      eigt::convertTransform<double>(
                        pcl_matcher_.getFinalTransformation().cast<double>()));
       old_scan_ = std::move(new_pcl);
       pcl_matcher_.setInputTarget (old_scan_);
    }
    return true;
}


int main(int argc, char **argv) {

  ros::init(argc, argv, "ml_ndt");
  ros::NodeHandle n;
  ros::NodeHandle n_private("~");
  MlNdt ndt(n, n_private);
  ndt.start();
  return 0;
}
