#include <scanmatcher.h>

void Scanmatcher::initialize(pose_t &pose, points2_t &points) {
  initializeNdt(pose, points);
}
void Scanmatcher::initialize(pose_t &pose, pcl_t &points) {
  initializeNdt(pose, projectPointsTo2D(points));
}
void Scanmatcher::calculate(pose_t &pose, points2_t &points) {
  calculateNdt(pose, points);
}
void Scanmatcher::calculate(pose_t &pose, pcl_t &points) {
  calculateNdt(pose, projectPointsTo2D(points));
}

Scanmatcher::pose_t Scanmatcher::calculate(pose_t &prev_pose,
                                           points2_t &first_scan,
                                           pose_t &curr_pose,
                                           points2_t &second_scan) {
  points2_t old_points = std::move(points_);
  pose_t old_pose = std::move(pose_);
  initializeNdt(prev_pose, first_scan);
  calculateNdt(curr_pose, second_scan);
  initializeNdt(old_pose, old_points);
  return transform_;
}

Scanmatcher::pose_t Scanmatcher::getTransformation() const{
  return transform_;
}

void Scanmatcher::setResolution(const size_t res){
  resolution_ = res;
}
void Scanmatcher::setLayers(const size_t layers){
  layers_count_ = layers;  
}

void Scanmatcher::setMaxRange(const float range){
  max_range_ = range;
}

// ****************** PRIVATE FUNCTIONS *********************************
Scanmatcher::points2_t Scanmatcher::projectPointsTo2D(pcl_t &points) {
  points2_t points2d;
  point_t one_point;
  for (auto pt : points) {
    one_point << pt.x,pt.y;
    points2d.emplace_back(one_point);
  }
  return points2d;
}
void Scanmatcher::initializeNdt(pose_t &pose, points2_t &points) {
  // initialize layers
  layer_.clear();
  points_ = std::move(points);
  layer_.push_back(Layer(&points_, resolution_, max_range_));
  for (size_t i = 1; i < layers_count_; ++i) {
    layer_.push_back(Layer(
        &points_, static_cast<size_t>(std::pow(resolution_, i)), max_range_));
  }

  pose_ = std::move(pose);

  initialized_ = true;
}

void Scanmatcher::initializeNdt(pose_t &pose, points2_t &&points) {
  initializeNdt(pose, points);
}

bool Scanmatcher::calculateNdt(pose_t &pose, points2_t &points) {
  if (!initialized_)
    return false;
  pose_t transform = calcTransformation(pose_, pose);
  for (auto &l : layer_) {
    if (l.calculateNdt(transform, points)) {
      transform = std::move(l.getTransformation());
    }
  }
  transform_ = std::move(transform);
  initializeNdt(pose, points);
  return true;
}

bool Scanmatcher::calculateNdt(pose_t &pose, points2_t &&points) {
  return calculateNdt(pose, points);
}

Scanmatcher::pose_t Scanmatcher::calcTransformation(pose_t &first_pose,
                                                    pose_t &second_pose) {
  pose_t transform;
  transform << second_pose(0) - first_pose(0), second_pose(1) - first_pose(1);
  float theta_first = first_pose(2) + PI_F;
  float theta_second = second_pose(2) + PI_F;
  transform(2) = theta_second - theta_first;
  return std::move(transform);
}


