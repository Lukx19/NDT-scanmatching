#include <ml_ndt_node.h>

MlNdt::MlNdt(ros::NodeHandle &n, ros::NodeHandle &n_private)
    : nh_(n), nh_private_(n_private) {

  // find tf_prefix if exists add it to tf_prefix_ variable
  tf_prefix_ = "";
  std::string tf_prefix_path;
  if (nh_.searchParam("tf_prefix", tf_prefix_path)) {
    nh_.getParam(tf_prefix_path, tf_prefix_);
  }

  robot_base_frame_ =
      nh_private_.param<std::string>("robot_base_frame", "base_link");
  world_farme_ = nh_private_.param<std::string>("world_farme", "odom");
  odom_topic_ = nh_private_.param<std::string>("odom_topic", "/odom");
  laser_topic_ = nh_private_.param<std::string>("laser_topic", "/laser");

  if (tf_prefix_ != "") {
    robot_base_frame_ = tf_prefix_ + "/" + robot_base_frame_;
  }

  // subscribers to acceleromer and gyroscope
  message_filters::Subscriber<nav_msgs::Odometry> odom_sub(nh_, odom_topic_,
                                                           1000);
  message_filters::Subscriber<sensor_msgs::LaserScan> laser_sub(
      nh_, laser_topic_, 1000);

  // sync messages using approximate alghorithm
  constexpr int sync_delay = 10;
  typedef message_filters::sync_policies::ApproximateTime<
      nav_msgs::Odometry, sensor_msgs::LaserScan> ImuSyncPolicy;
  message_filters::Synchronizer<ImuSyncPolicy> msg_sync(
      ImuSyncPolicy(sync_delay), odom_sub, laser_sub);
  msg_sync.registerCallback(boost::bind(&MlNdt::data_cb, this, _1, _2));
}

void MlNdt::data_cb(const nav_msgs::Odometry::ConstPtr &odom,
                    const sensor_msgs::LaserScan::ConstPtr &laser) {

  sensor_msgs::PointCloud2 laser_msg;
  sensor_msgs::PointCloud2 laser_base;
  pcl::PointCloud<pcl::PointXYZ> laser_base_points;

  projector_.projectLaser(*laser, laser_msg);
  tf_list_.waitForTransform(robot_base_frame_, laser->header.frame_id,
                            laser->header.stamp, ros::Duration(5.0));
  pcl_ros::transformPointCloud(robot_base_frame_, laser_msg, laser_base,
                               tf_list_);
  pcl::fromROSMsg(laser_base, laser_base_points);
  // std::cout<<laser_base_points.getMatrixXfMap()<<std::endl;

  tf::Pose pose_tf;
  Eigen::Vector3f pose;
  tf_list_.waitForTransform(world_farme_, odom->header.frame_id,
                            odom->header.stamp, ros::Duration(5.0));

  tf::poseMsgToTF(odom->pose.pose, pose_tf);
  pose << pose_tf.getOrigin().getX(), pose_tf.getOrigin().getY(),
      tf::getYaw(odom->pose.pose.orientation);
}

int main(int argc, char **argv) {

  ros::init(argc, argv, "ml_ndt");
  ros::NodeHandle n;
  ros::NodeHandle n_private("~");
  MlNdt ndt(n, n_private);
  ros::spin();
  return 0;
}
