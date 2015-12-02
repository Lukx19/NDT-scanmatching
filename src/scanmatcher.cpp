#include <scanmatcher.h>

void Scanmatcher::initialize(pose_t &pose, points2_t &points) {
  initializeNdt(pose, points);
}
void Scanmatcher::initialize(pose_t &pose, pcl_t &points) {
  initializeNdt(pose, projectPointsTo2D(points));
}
bool Scanmatcher::calculate(pose_t &pose, points2_t &points) {
  return calculateNdt(pose, points);
}
bool Scanmatcher::calculate(pose_t &pose, pcl_t &points) {
  return calculateNdt(pose, projectPointsTo2D(points));
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
  return getTransformation();
}

Scanmatcher::pose_t Scanmatcher::getTransformation() const{
  pose_t pose;
  float angle = getAngleFromTransform(transform_);
  pose<<transform_.translation(),angle;
  return pose;
}

tf::Transform Scanmatcher::getTFTransform()const{
  tf::Transform msg;
  transform_t t = calcTransformation(pose_,last_odom_);
  tf::Quaternion orientation;
  float angle = getAngleFromTransform(t);
  orientation.setRPY(0, 0,angle);
  msg.setOrigin(
      tf::Vector3(t.translation()(0), t.translation()(1), 0));
  msg.setRotation(orientation);
  return msg;
}

Scanmatcher::pose_t Scanmatcher::getPose()const{
  return pose_;
}

nav_msgs::Odometry Scanmatcher::getOdom()const{
  nav_msgs::Odometry odom;
  pose_t p = getPose();
  odom.pose.pose.position.x=p(0);
  odom.pose.pose.position.y=p(1);
  odom.pose.pose.position.z=0;

  tf::Quaternion orientation;
  orientation.setRPY(0, 0, p(2));
  odom.pose.pose.orientation.x = orientation.getX();
  odom.pose.pose.orientation.y = orientation.getY();
  odom.pose.pose.orientation.z = orientation.getZ();
  odom.pose.pose.orientation.w = orientation.getW();
  return odom;
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

void Scanmatcher::createLayers(){
  for (size_t i = 0; i < layers_count_; ++i) {
    size_t power = static_cast<size_t>(pow(2,i));
    layer_.push_back(Layer(&points_,resolution_ * power, max_range_));
  }
}

void Scanmatcher::initializeNdt(pose_t &pose,points2_t & points) {
  // initialize layers
  layer_.clear();
  //DEBUG("ML-NDT: Layer created with"<<points_.size()<<"pts"); 
  points_.clear();

  for(auto & pt: points)
    points_.push_back(pt);
  createLayers();
  pose_ = pose;
  last_odom_ = pose;

  initialized_ = true;
}

void Scanmatcher::initializeNdt(pose_t &pose, points2_t && points) {
  layer_.clear();
  //DEBUG("ML-NDT: Layer created with2"<<points_.size()<<"pts"); 
  points_.clear();

  for(auto && pt: points)
    points_.push_back(std::move(pt));
  createLayers();
  pose_ = pose;
  last_odom_ = pose;

  initialized_ = true;

}

void Scanmatcher::updateLayers(const pose_t & calc_pose,const pose_t & odom,const points2_t & points){
  layer_.clear();
  points_.clear();
  for(auto && pt: points)
    points_.push_back(std::move(pt));
  createLayers();
  pose_ = calc_pose;
  last_odom_ = odom;
}

bool Scanmatcher::calculateNdt(pose_t & current_pose_odom, points2_t &points) {
  if (!initialized_)
    return false;
 //transform_t transform = calcTransformation(last_odom_,current_pose_odom);
  transform_t transform = transform_t::Identity();
  //float angle = getAngleFromTransform(transform);
  //DEBUG("tranform in odom \n"<<transform.matrix());
  //if(std::abs(angle)< 0.261799388)
  // return false;
  for (auto &l : layer_) {
    if (l.calculateNdt(transform, points)) {
      transform = std::move(l.getTransformation());
    }
  }
  DEBUG("tranform calculated \n"<<transform.matrix())
  transform_ = std::move(transform);
  
  pose_t transformed_pose =transformPose(pose_,transform);
  DEBUG("calculated transformed pose \n"<<transformed_pose);
  DEBUG("odometry pose \n"<<current_pose_odom);
  updateLayers(transformed_pose,current_pose_odom, points);
  return true;
}

bool Scanmatcher::calculateNdt(pose_t &pose, points2_t &&points) {
  return calculateNdt(pose, points);
}

Scanmatcher::transform_t Scanmatcher::calcTransformation(const pose_t &from,
                                                    const pose_t & to)const {
  // p0 ---> t --> p1
  transform_t p0;
  transform_t p1;
  transform_t t;
  p0 = p0.fromPositionOrientationScale(from.head<2>(), 
                                  Eigen::Rotation2Df(from(2)), 
                                  point_t::Ones());
  p1 = p1.fromPositionOrientationScale(to.head<2>(),
                                  Eigen::Rotation2Df(to(2)),
                                  point_t::Ones());
  t =  p1 * p0.inverse();
  return std::move(t);
}

Scanmatcher::pose_t Scanmatcher::transformPose(const pose_t & pose,const transform_t &trans) const{
  transform_t t;
  pose_t p;
  t = t.fromPositionOrientationScale(pose.head<2>(),
                                  Eigen::Rotation2Df(pose(2)),
                                 point_t::Ones());
  t=trans*t;
  float angle = getAngleFromTransform(t); 
  p<<t.translation()(0),t.translation()(1),angle;
  return p;
}

float Scanmatcher::getAngleFromTransform(const transform_t & trans) const{
  return std::atan2(trans.rotation()(1,0),trans.rotation()(0,0)); 
}

