#include <gtest/gtest.h>
#include <ros/ros.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/TransformStamped.h>
#include <sensor_msgs/LaserScan.h>
#include <Eigen/Dense>
#include <fstream>
#include <eigen_tools.h>

double EPSILON = 0.001;
ros::Subscriber sub_node_transforms;
ros::Publisher pub_odom;
ros::Publisher pub_laser_scan;
std::vector<sensor_msgs::LaserScan> laser_scans;
std::vector<nav_msgs::Odometry> real_odom;
std::vector<geometry_msgs::TransformStamped::ConstPtr> trans_calc;

void prepareLaserScans(){
    laser_scans.clear();
    std::fstream laser_msg;
    laser_msg.open("data/data.scans",std::ios_base::in);
    if(!laser_msg){
        std::cout<<"File data/data.scans NOT FOUND"<<std::endl;
    }
    float num_f;
    unsigned long num_l;
    std::string line;
    ros::Time time;
    unsigned int seq =0;
    while(std::getline(laser_msg,line)){
        std::stringstream str(line);
        sensor_msgs::LaserScan msg;
        msg.header.seq = seq;
        str >> num_l;
        msg.header.stamp = time.fromNSec(num_l);
        msg.header.frame_id = "base_link";
        str.ignore();
        str >> num_f;
        msg.angle_min = num_f;
        msg.angle_max = 2.26456475258f;
        msg.range_min = 0.0230000000447f;
        msg.range_max = 60.0f;
        msg.time_increment =  1.73611115315e-05f;
        str.ignore();
        str >> num_f;
        msg.angle_increment = num_f;
        str.ignore();
        str >> num_l;
        for(size_t i=0; i< num_l ; ++i){
            str.ignore();
            str>>num_f;
            msg.ranges.push_back(num_f);
        }
        laser_scans.push_back(msg);
        ++seq;
    }
    laser_msg.close();
}

void preparePoseData(){
    real_odom.clear();
    std::fstream pose_msg;
    pose_msg.open("data/data.poses",std::ios_base::in);
    if(!pose_msg){
        std::cout<<"File data/data.poses NOT FOUND"<<std::endl;
    }
    float num_f;
    unsigned long num_l;
    std::string line;
    ros::Time time;
    unsigned int seq =0;
    while(std::getline(pose_msg,line)){
        std::stringstream str(line);
        nav_msgs::Odometry msg;
        msg.header.seq = seq;
        str>>num_l;
        str.ignore();
        msg.header.stamp = time.fromNSec(num_l);
        msg.header.frame_id = "odom";
        str >> num_f;
        str.ignore();
        msg.pose.pose.position.x =num_f;
        str >> num_f;
        str.ignore();
        msg.pose.pose.position.y =num_f;
        msg.pose.pose.position.z = 0;
        str >> num_f;
        Eigen::Quaternionf rot(Eigen::AngleAxisf(num_f,Eigen::Vector3f::UnitZ()));
        msg.pose.pose.orientation.x = rot.x();
        msg.pose.pose.orientation.y = rot.y();
        msg.pose.pose.orientation.z = rot.z();
        msg.pose.pose.orientation.w = rot.w();
        real_odom.push_back(msg);
        ++seq;      
    }
    pose_msg.close();
}

void transformsCallback(const geometry_msgs::TransformStamped::ConstPtr & trans ){
    trans_calc.push_back(trans);
}

eigt::pose2d_t<double> getPoseFromOdomMsg(const nav_msgs::Odometry & msg){
    eigt::pose2d_t<double> pose;
    auto rot_msg = msg.pose.pose.orientation;
    Eigen::Quaterniond rot(rot_msg.x,rot_msg.y,rot_msg.z,rot_msg.w); 
    pose << msg.pose.pose.position.x,msg.pose.pose.position.y,rot.toRotationMatrix().eulerAngles(0,1,2)(2);
    return pose;
}

eigt::pose2d_t<double> getTransFromTransMsg(const geometry_msgs::TransformStamped & msg){
    eigt::pose2d_t<double> trans;
    auto rot_msg = msg.transform.rotation;
    Eigen::Quaterniond rot(rot_msg.x,rot_msg.y,rot_msg.z,rot_msg.w); 
    trans << msg.transform.translation.x,msg.transform.translation.y,rot.toRotationMatrix().eulerAngles(0,1,2)(2);
    return trans;
}


TEST(MlNdtNode, twoLaserScansTransformationCheck)
{
  trans_calc.clear();
  pub_odom.publish(real_odom[1]);
  pub_laser_scan.publish(laser_scans[1]);
  pub_odom.publish(real_odom[2]);
  ros::spinOnce();
  // spin event loop until receiving calculated odom
  while(trans_calc.size() <= 0 && ros::ok()){
    ros::spinOnce();
  }
  auto real_trans = eigt::getPoseFromTransform(eigt::transBtwPoses(
      getPoseFromOdomMsg(real_odom[1]), getPoseFromOdomMsg(real_odom[1])));
  auto calc_trans = getTransFromTransMsg(*trans_calc[0]);
  EXPECT_LT((real_trans - calc_trans).cwiseAbs().sum(),EPSILON);
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "ml_ndt_tester");
  preparePoseData();
  prepareLaserScans();
  ros::NodeHandle n;
  ros::NodeHandle n_private("~");
  sub_node_transforms = n.subscribe("transform",1000, transformsCallback);
  pub_odom = n.advertise<nav_msgs::Odometry>("odom",1000);
  pub_laser_scan = n.advertise<sensor_msgs::LaserScan>("laser",1000);
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
