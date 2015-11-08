#ifndef SCANMATCHER
#define SCANMATCHER

#include <exception>

#include <Eigen/Dense>
#include "layer.h"
#include "field.h"
class Field;
class Layer;

class Scanmatcher {
public:
  typedef size_t Id_t;
  typedef Eigen::Vector3f pose_t;
  typedef Eigen::Vector2f point_t;
  typedef std::vector<point_t> points2_t;
  typedef std::vector<Eigen::Vector3f> points3_t;
  enum scan_type { S180, S360 };
  Scanmatcher(size_t max_range, size_t resolution, size_t layers,
              scan_type scan)
      : pose_(pose_t::Zero()), transform_(pose_t::Zero()),
        max_range_(max_range), resolution_(resolution), layers_count_(layers),
        initialized_(false), s_type_(scan) {}

  void initialize(pose_t &pose, points2_t &points);
  void initialize(pose_t &pose, points3_t &points);
  void calculate(pose_t &pose, points2_t &points);
  void calculate(pose_t &pose, points3_t &points);

  pose_t calculate(pose_t &prev_pose, points2_t &first_scan, pose_t &curr_pose,
                   points2_t &second_scan);

  pose_t getTransformation() const;
  point_t getPoint(Id_t id) const;

private:
  pose_t pose_;
  pose_t transform_;
  size_t max_range_;
  size_t resolution_;
  size_t layers_count_;
  bool initialized_;
  scan_type s_type_;
  points2_t points_;
  std::vector<Layer> layer_;

  const float PI_F = 3.14159265358979f;

  points2_t projectPointsTo2D(points3_t &points);

  bool calculateNdt(pose_t &pose, points2_t &points);
  bool calculateNdt(pose_t &pose, points2_t &&points);

  void initializeNdt(pose_t &pose, points2_t &points);
  void initializeNdt(pose_t &pose, points2_t &&points);
  // calculates transform between last stored pose and new pose as a parameter
  pose_t calcTransformation(pose_t &first_pose, pose_t &second_pose);
};

#endif
