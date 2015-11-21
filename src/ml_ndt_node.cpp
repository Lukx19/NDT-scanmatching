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

  odom_sub_.subscribe(nh_, odom_topic_, 1000);
  laser_sub_.subscribe(nh_, laser_topic_, 1000);

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
  Eigen::Vector3f pose;
  tf_list_.waitForTransform(odom_frame_, odom->header.frame_id,
                            odom->header.stamp, ros::Duration(5.0));

  tf::poseMsgToTF(odom->pose.pose, pose_tf);
  pose << pose_tf.getOrigin().getX(), pose_tf.getOrigin().getY(),
      tf::getYaw(odom->pose.pose.orientation);
  ROS_INFO("ML-NDT: Messages transformed.");

  if (!is_initialized) {
    matcher_.initialize(pose, laser_pcl_base);
    is_initialized = true;
    ROS_INFO("ML-NDT: First laser scan inserted to structure");
  } else {
    // calculate new odometry based on scan matching
    matcher_.calculate(pose, laser_pcl_base);
    pose_t new_pose = matcher_.getTransformation();
    ROS_INFO("ML-NDT: New odom calculated %f  %f  %f", new_pose[0], new_pose[1],
             new_pose[2]);
    // fill in new odometry msg
    nav_msgs::Odometry msg;
    msg.header.frame_id = new_odom_frame_;
    msg.header.stamp = ros::Time::now();
    msg.header.seq = seq_;
    msg.twist = odom->twist;
    tf::Quaternion orientation;
    orientation.setRPY(0, 0, new_pose[2]);
    msg.pose.pose.orientation.x = orientation.getX();
    msg.pose.pose.orientation.y = orientation.getY();
    msg.pose.pose.orientation.z = orientation.getZ();
    msg.pose.pose.orientation.w = orientation.getW();
    msg.pose.pose.position.x = new_pose[0];
    msg.pose.pose.position.y = new_pose[1];
    msg.pose.pose.position.z = 0;
    msg.pose.covariance = odom->pose.covariance;
    new_odom_pub_.publish(msg);
    tf::Transform transform;
    transform.setOrigin(
        tf::Vector3(pose[0] - new_pose[0], pose[1] - new_pose[1], pose[2]));
    transform.setRotation(pose_tf.getRotation() - orientation);
    tf_broadcast_.sendTransform(tf::StampedTransform(
        transform, ros::Time::now(), new_odom_frame_, odom_frame_));
    ++seq_;
    ROS_INFO("ML-NDT: New odom sent");
  }
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
      static_cast<float>(nh_private_.param<double>("resolution", 0.5));

  max_range_ =
      static_cast<float>(nh_private_.param<double>("maximal_laser_range", 4.0));

  layers_ = static_cast<size_t>(nh_private_.param<int>("number_of_layers", 4));

  if (tf_prefix_ != "") {
    robot_base_frame_ = tf_prefix_ + "/" + robot_base_frame_;
    odom_frame_ = tf_prefix_ + "/" + odom_frame_;
    new_odom_frame_ = tf_prefix_ + "/" + new_odom_frame_;
  }
}

int main(int argc, char **argv) {

  ros::init(argc, argv, "ml_ndt");
  ros::NodeHandle n;
  ros::NodeHandle n_private("~");
  MlNdt ndt(n, n_private);
  ndt.start();
  return 0;
}
