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

#include <pcl_conversions/pcl_conversions.h>
#include <pcl_ros/transforms.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>


#include <Eigen/Dense>


class MlNdt{
public:
    MlNdt(ros::NodeHandle & n,ros::NodeHandle & n_private);
      
private:
    ros::NodeHandle nh_;
    ros::NodeHandle nh_private_;
    std::string world_farme_;
    std::string robot_base_frame_;
    std::string tf_prefix_;
    std::string odom_topic_;
    std::string laser_topic_;
    laser_geometry::LaserProjection projector_;
    tf::TransformListener tf_list_;
    void init();
    void data_cb(const nav_msgs::Odometry::ConstPtr & odom, const sensor_msgs::LaserScan::ConstPtr& laser);
    void pointCloud2Data_cb();
    
};



#endif
