#include <ml_ndt_scanmatching/ml_ndt_node.h>

// ******************************* CONSTRUCTORS ********************

MlNdt::MlNdt(ros::NodeHandle &n, ros::NodeHandle &n_private)
    : nh_(n), nh_private_(n_private),is_ready_(false), msg_sync_(ImuSyncPolicy(10)), seq_(0),
      is_initialized(false) {

  initParameters();
  matcher_.setMaxRange(max_range_);
  matcher_.setResolution(resolution_);
  matcher_.setLayers(layers_);

  Eigen::Vector2f step;
  step<< 0.5,0.5;
  pcl_matcher_.setGridStep(step);
  pcl_matcher_.setOptimizationStepSize(0.05);
  pcl_matcher_.setOptimizationStepSize(0.5);
  pose_pub_ =
      nh_.advertise<geometry_msgs::PoseWithCovarianceStamped>(pose_pub_topic_, 100, false);

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
  ROS_INFO_STREAM("ML-NDT: Laser points received: "<< laser_pcl_base.size());
  if(laser_pcl_base.size() == 0)
    return;

  // transform robot odometry too odometry frame
  tf::Pose pose_tf;
  Eigen::Vector3d pose;
  //tf_list_.waitForTransform(odom_frame_, odom->header.frame_id,
   //                         odom->header.stamp, ros::Duration(5.0));

  tf::poseMsgToTF(odom->pose.pose, pose_tf);
  pose << pose_tf.getOrigin().getX(), pose_tf.getOrigin().getY(),
      tf::getYaw(odom->pose.pose.orientation);
  ROS_INFO("ML-NDT: Messages transformed.");
  // prepare initial data
  if(!is_ready_){
    old_odom_ = pose;
    old_position_ = pose;
    is_ready_ = true;
  }
  //publishTFTransform(eigt::transBtwPoses(old_position_,pose));
  // calculate transform and return new position
  eigt::transform2d_t<double> transformation;
  Eigen::Matrix3d covar(Eigen::Matrix3d::Identity());
  if(mode_ == MLNDT2D){
    if(!getTransfMlNdt(transformation,pose,laser_pcl_base)){
      return;
    }
  }else if(mode_ == PCL2D){
    
    if(!getTransfPclNdt(transformation,covar,pose,laser_pcl_base))
      return;
  }else if(mode_ == D2DNDT){
    if(!getTransfMlNdt3d(transformation,pose,laser_pcl_base))
      return;
  }
  eigt::pose2d_t<double> calc_pose = eigt::transformPose(old_odom_,transformation);
  ROS_INFO_STREAM("NDT res trans: \n"<<eigt::getPoseFromTransform(transformation).transpose());
  ROS_INFO_STREAM("NDT res pose:"<< calc_pose.transpose());
  //publishTFTransform(eigt::transBtwPoses(calc_pose,pose));
  

  

  geometry_msgs::PoseWithCovarianceStamped p_msg;
  p_msg.header.seq = seq_;
  p_msg.header.stamp = laser->header.stamp;
  p_msg.header.frame_id = pose_frame_;
  p_msg.pose.pose.position.x = calc_pose(0);
  p_msg.pose.pose.position.y = calc_pose(1);
  p_msg.pose.pose.position.z = 0;
  tf::Quaternion orient;
  orient.setRPY(0,0,calc_pose(2));
  p_msg.pose.pose.orientation.x = orient.getX();
  p_msg.pose.pose.orientation.y = orient.getY();
  p_msg.pose.pose.orientation.z = orient.getZ();
  p_msg.pose.pose.orientation.w = orient.getW();
  p_msg.pose.covariance = convertCovariance(covar);
  pose_pub_.publish(p_msg);
  
  old_position_ = calc_pose;
  old_odom_ = pose;
  ROS_INFO("ML-NDT: DATA sent");
    ++seq_;
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

  pose_frame_ =
      nh_private_.param<std::string>("new_pose_farme_id", "odom_clac");

  odom_topic_ = nh_private_.param<std::string>("odom_topic", "/odom");

  pose_pub_topic_ =
      nh_private_.param<std::string>("calculated_pose_topic", "/pose_calc");

  laser_topic_ = nh_private_.param<std::string>("laser_topic", "/laser");

  resolution_ =
      static_cast<size_t>(nh_private_.param<int>("resolution", 8));

  max_range_ =
      static_cast<float>(nh_private_.param<double>("maximal_laser_range", 4.0));
  
  min_angle_ = nh_private_.param<double>("min_rotated_angle", 0);

  min_displacement_ = nh_private_.param<double>("min_traveled_distance", 0);

  layers_ = static_cast<size_t>(nh_private_.param<int>("number_of_layers", 4));

  std::string mode = nh_private_.param<std::string>("scanmatcher_mode", "MLNDT2D");
  if(mode == "MLNDT2D"){
    mode_ = MLNDT2D;
  }else if( mode == "PCL2D"){
    mode_ = PCL2D;
  } else if (mode == "D2DNDT"){
    mode_ = D2DNDT;
  }

  if (tf_prefix_ != "") {
    robot_base_frame_ = tf_prefix_ + "/" + robot_base_frame_;
    odom_frame_ = tf_prefix_ + "/" + odom_frame_;
    pose_frame_ = tf_prefix_ + "/" + pose_frame_;
  }
}

bool MlNdt::getTransfMlNdt(eigt::transform2d_t<double> & trans,
                           const pose_t & odom_pose,
                           const pcl_t & points)
{
  if (!is_initialized) {
    old_scan_ = points.makeShared();
    is_initialized = true;
    return false;
  } else {
    eigt::transform2d_t<double> odom_trans =
                                eigt::transBtwPoses(old_odom_,odom_pose);
    double distance =  std::abs(eigt::getDisplacement(odom_trans));
    double angle = std::abs(eigt::getAngleDiffrence(old_odom_,odom_pose));
    std::cout << "old position" << old_position_ << "\n\nnew odom" << odom_pose
              << "\n angle: " << angle << "\n\n"
              << "displacement: " << distance << "\n";
    std::cout << "init guess tranform \n" << eigt::getPoseFromTransform(odom_trans)
              << "\n\n";
    if(distance < min_displacement_ && angle <min_angle_){
       return false;
     }
    // calculate new odometry based on scan matching
    pcl_t::Ptr source = points.makeShared();
    if(!matcher_.matchBiber(*source,*old_scan_,trans,odom_trans))
      return false;
    old_scan_ = std::move(source);
  }
  return true;
}

bool MlNdt::getTransfPclNdt(eigt::transform2d_t<double> & trans,
                            Eigen::Matrix3d & covar,
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
                                  eigt::transBtwPoses(old_odom_,odom_pose);
      //std::cout<<odom_trans.matrix()<<"\n";
      //odom_trans.setIdentity();
      double distance =  std::abs(eigt::getDisplacement(odom_trans));
      double angle = std::abs(eigt::getAngleDiffrence(odom_pose,old_odom_));
      std::cout<<"old position"<<old_position_<<"\n\nnew odom"<<odom_pose<<"\n angle: "<<angle<<"\n\n";
      std::cout<<"tranform \n"<<eigt::getPoseFromTransform(odom_trans)<<"\n\n";
      if(distance < min_displacement_ && angle < min_angle_){
         return false;
       }
      pcl_t::Ptr new_pcl = points.makeShared();
      pcl_matcher_.setInputSource (new_pcl);
      // Set initial alignment estimate found using robot odometry.
        Eigen::Matrix<double, 4, 4>init_guess = eigt::convertFromTransform(odom_trans); 
       // Calculating required rigid transform to align the input cloud to the target cloud.
       pcl_t::Ptr output_cloud (new pcl_t());
       init_guess.setIdentity();
       pcl_matcher_.align (*output_cloud,init_guess.cast<float>());

       ROS_INFO_STREAM("Normal Distributions Transform has converged:" << pcl_matcher_.hasConverged ()
            << " score: " << pcl_matcher_.getFitnessScore ());
       trans = eigt::convertToTransform<double>(
                        pcl_matcher_.getFinalTransformation().cast<double>()).inverse();
       covar = pcl_matcher_.getCovariance();
       old_scan_ = std::move(new_pcl);
       pcl_matcher_.setInputTarget (old_scan_);
    }
    return true;
}

bool MlNdt::getTransfMlNdt3d(eigt::transform2d_t<double> & trans,
                          const pose_t & odom_pose,
                          const pcl_t & points){
   if(!is_initialized){
      old_scan_ = points.makeShared();
      is_initialized = true;
      return false;
    }else{
      eigt::transform2d_t<double> odom_trans =
                                  eigt::transBtwPoses(old_odom_,odom_pose);
      double distance =  std::abs(eigt::getDisplacement(odom_trans));
      double angle = std::abs(eigt::getAngleDiffrence(old_odom_,odom_pose));
      std::cout<<"old position"<<old_position_<<"\n\nnew odom"<<odom_pose<<"\n angle: "<<angle<<"\n\n";
      std::cout<<"tranform \n"<<eigt::getPoseFromTransform(odom_trans)<<"\n\n";
      if(distance < min_displacement_ && angle <min_angle_){
         return false;
       }
      Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T(eigt::convertFromTransform(odom_trans));
      pcl_t::Ptr source = points.makeShared();
      //... load information into cloud1 and cloud2 ...
      bool ret = d2d_matcher_.match(*(old_scan_),*source,T);
      std::cout<<ret<"\n";
      std::cout<<T.matrix();
      trans = eigt::convertToTransform(T.matrix());
      //we now have T, which minimizes |cloud1 - cloud2*T| 
      //lslgeneric::transformPointCloudInPlace<pcl::PointXYZ>(T,cloud2);
      old_scan_ = std::move(source);
    }
  return true;
}

void MlNdt::publishTFTransform(const eigt::transform2d_t<double> & tf_transform){
  tf::Transform t_msg;
  tf::Quaternion orientation;
  double angle = eigt::getAngle(tf_transform);
  orientation.setRPY(0, 0,angle);
  t_msg.setOrigin(
      tf::Vector3(tf_transform.translation()(0), tf_transform.translation()(1), 0));
  t_msg.setRotation(orientation);
  tf_broadcast_.sendTransform(tf::StampedTransform(
         t_msg, ros::Time::now(), pose_frame_, odom_frame_));
}

boost::array<double,36> MlNdt::convertCovariance(const Eigen::Matrix3d & covar) const{
  boost::array<double,36> out_cov;
  for(size_t i =0; i< 36;++i){
    out_cov[i] = 0;
  }
  out_cov[0] = covar(0,0);
  out_cov[1] = covar(0,1);
  out_cov[6] = covar(1,0);
  out_cov[7] = covar(1,1);
  out_cov[5] = covar(0,2);
  out_cov[11] = covar(1,2);
  out_cov[30] = covar(2,0);
  out_cov[31] = covar(2,1);
  out_cov[35] = covar(2,2);
  return out_cov;
}

int main(int argc, char **argv) {

  ros::init(argc, argv, "ml_ndt");
  ros::NodeHandle n;
  ros::NodeHandle n_private("~");
  MlNdt ndt(n, n_private);
  ndt.start();
  return 0;
}
