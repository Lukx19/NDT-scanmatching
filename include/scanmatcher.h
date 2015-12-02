#ifndef SCANMATCHER
#define SCANMATCHER

#include <exception>
#include <memory>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <nav_msgs/Odometry.h>
#include <tf/transform_datatypes.h>

#include <Eigen/Dense>
#include <Eigen/Geometry> 

#include <layer.h>
#include <field.h>

#ifndef DEBUG
  #define DEBUG(out); std::cout<<out<<"\n";
#endif

class Field;
class Layer;

class Scanmatcher {
public:
  typedef size_t Id_t;
  typedef Eigen::Vector3f pose_t;
  typedef Eigen::Vector2f point_t;
  typedef Eigen::Transform<float,2, Eigen::TransformTraits::Affine> transform_t;
  typedef pcl::PointCloud<pcl::PointXYZ> pcl_t;
  typedef std::vector<point_t> points2_t;
  typedef std::vector<Eigen::Vector3f> points3_t;
  enum scan_type { S180, S360 };
  Scanmatcher(float max_range=4.0, size_t resolution=8, size_t layers=4,
              scan_type scan=scan_type::S180)
      : pose_(pose_t::Zero()), max_range_(max_range), resolution_(resolution), layers_count_(layers),
        initialized_(false), s_type_(scan), points_(points2_t(360)) {}

  void initialize(pose_t &pose, points2_t &points);
  void initialize(pose_t &pose, pcl_t &points);
  bool calculate(pose_t &pose, points2_t &points);
  bool calculate(pose_t &pose, pcl_t &points);

  pose_t calculate(pose_t &prev_pose, points2_t &first_scan, pose_t &curr_pose,
                   points2_t &second_scan);

  pose_t getTransformation() const;
  tf::Transform getTFTransform() const;
  pose_t getPose() const;
  nav_msgs::Odometry getOdom() const;

  point_t getPoint(Id_t id) const;

  void setResolution(const size_t res);
  void setLayers(const size_t layers);
  void setMaxRange(const float range);

private:
  pose_t pose_;
  pose_t last_odom_;
  transform_t transform_;
 // transformation beetween last odometry pose and calculated pose
  //transform_t odom_pose_trans_;
  float max_range_;
  size_t resolution_;
  size_t layers_count_;
  bool initialized_;
  scan_type s_type_;
  points2_t points_;
  std::vector<Layer> layer_;

  const float PI_F = 3.14159265358979f;

  points2_t projectPointsTo2D(pcl_t &points);
  void createLayers();

  void initializeNdt(pose_t &pose, points2_t &points);
  void initializeNdt(pose_t &pose, points2_t &&points);

  void updateLayers(const pose_t & calc_pose,const pose_t & odom,const points2_t & points);

  bool calculateNdt(pose_t &pose, points2_t &points);
  bool calculateNdt(pose_t &pose, points2_t &&points);
  // calculates transform between last stored pose and new pose as a parameter
  transform_t calcTransformation(const pose_t &first_pose,const pose_t &second_pose) const;
  pose_t transformPose(const pose_t &pose,const transform_t &trans) const;
  float getAngleFromTransform(const transform_t & trans) const;

};

#endif
