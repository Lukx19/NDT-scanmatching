#ifndef ML_NDT
#define ML_NDT

#include <string>
#include <iostream>
#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <nav_msgs/Odometry.h>
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <message_filters/synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <laser_geometry/laser_geometry.h>
#include <sensor_msgs/PointCloud2.h>

#include <tf/transform_datatypes.h>
#include <tf/transform_broadcaster.h>

#include <pcl_conversions/pcl_conversions.h>
#include <pcl_ros/transforms.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <Eigen/Dense>
#include <scanmatcher.h>

class MlNdt {
public:
  typedef message_filters::sync_policies::ApproximateTime<
      nav_msgs::Odometry, sensor_msgs::LaserScan> ImuSyncPolicy;
  typedef message_filters::Subscriber<nav_msgs::Odometry> odom_sub_t;
  typedef message_filters::Subscriber<sensor_msgs::LaserScan> laser_sub_t;
  typedef Eigen::Vector3f pose_t;
  MlNdt(ros::NodeHandle &n, ros::NodeHandle &n_private);
  void start();
private:
  ros::NodeHandle nh_;
  ros::NodeHandle nh_private_;
  Scanmatcher matcher_;

  // parameters from launch file
  std::string odom_frame_;
  std::string new_odom_frame_;
  std::string robot_base_frame_;
  std::string tf_prefix_;
  std::string odom_topic_;
  std::string new_odom_topic_;
  std::string laser_topic_;
  float max_range_;
  float resolution_;
  size_t layers_;

  laser_geometry::LaserProjection projector_;
  tf::TransformListener tf_list_;
  tf::TransformBroadcaster tf_broadcast_;
  ros::Publisher new_odom_pub_;
  odom_sub_t  odom_sub_;
  laser_sub_t laser_sub_;
  message_filters::Synchronizer<ImuSyncPolicy> msg_sync_;
  uint seq_;
  bool is_initialized;


  void data_cb(const nav_msgs::Odometry::ConstPtr &odom,
               const sensor_msgs::LaserScan::ConstPtr &laser);
  void initParameters();
};

#endif
