#include <ml_ndt_scanmatching/other/scanmatcher.h>

void Scanmatcher::initialize(const pose_t &pose, const points2_t &points) {
  initializeNdt(pose, points);
}

void Scanmatcher::initialize(const pose_t &pose, const pcl_t &points) {
  initializeNdt(pose, projectPointsTo2D(points));
}

bool Scanmatcher::calculate(const pose_t &pose, const points2_t &points) {
  return calculateNdt(pose, points);
}

bool Scanmatcher::calculate(const pose_t &pose,const pcl_t &points) {
  return calculateNdt(pose, projectPointsTo2D(points));
}

Scanmatcher::pose_t Scanmatcher::calculate(const pose_t &prev_pose,
                                           const points2_t &first_scan,
                                           const pose_t &curr_pose,
                                           const points2_t &second_scan) {
  points2_t old_points = std::move(points_);
  pose_t old_pose = std::move(pose_);
  initializeNdt(prev_pose, first_scan);
  calculateNdt(curr_pose, second_scan);
  initializeNdt(old_pose, old_points);
  return getTransformation();
}

bool Scanmatcher::match(const pcl_t &source, const pcl_t &target,
                        transform_t &trans, const transform_t & init_guess)
{
  std::vector<Layer> layers;
  points2_t target_pts = projectPointsTo2D(target);
  transform_t offset;
  //create layers
  offset.setIdentity();
  //offset(1,2)-=1.0;
  for (size_t i = 0; i < layers_count_; ++i) {
    size_t power = static_cast<size_t>(pow(2,i));
    layers.emplace_back(Layer(&target_pts,resolution_ * power, max_range_,offset));
  }

  // start calculation
  DEBUG("start calc");
  //trans = init_guess;
  trans.setIdentity();
  points2_t source_pts = projectPointsTo2D(source);
  for (auto &l : layers) {
     if (l.calculateNdt(trans, source_pts)) {
       trans = l.getTransformation();
     }else{
      return false;
     }
  }
  return true;
}

bool Scanmatcher::matchBiber(const pcl_t &source, const pcl_t &target,
                        transform_t &trans, const transform_t & init_guess)
{
  std::vector<Layer> layers;
  points2_t target_pts = projectPointsTo2D(target);
  transform_t offset;
  offset.setIdentity();
  //create layers
  offset.setIdentity();
  //offset(1,2)-=1.0;
  double shift = max_range_ / resolution_;
  std::array<transform_t,4> offsets;
  offsets[1] = Eigen::Translation<double ,2>(0,shift);
  offsets[2] = Eigen::Translation<double ,2>(shift,0);
  offsets[3] = Eigen::Translation<double ,2>(shift,shift);
  for (size_t i = 0; i < 4; ++i) {
    layers.emplace_back(Layer(&target_pts,resolution_, max_range_,offset));
  }

  // calculate optimalisation data (hessian,...)
  DEBUG("start calc");
  points2_t source_pts = projectPointsTo2D(source);
  bool converged = false;
  eigt::pose2d_t<double> p , old_p , delta_p, temp;
  p = eigt::getPoseFromTransform(init_guess);
  old_p = p;
  size_t iter =0;
  OptimalisationData opt_data; 
  while(!converged){
    for (auto &l : layers) {
       opt_data += l.calcOptData(eigt::getTransFromPose(p),source_pts);
    }
    opt_data.makeToSPD2();
    delta_p = opt_data.hessian.ldlt().solve(-opt_data.gradient);
    temp = p+delta_p;
    DEBUG("score after Newton: "<< scoreLayers(eigt::getTransFromPose(temp),source_pts));
    double step = 0.05;
    delta_p = delta_p * step;
    p+= delta_p;
    DEBUG("score after LS: "<< scoreLayers(eigt::getTransFromPose(p),source_pts));
    if((p - old_p).cwiseAbs().sum() < EPSILON)
      converged = true;
    if(iter >= MAX_ITER)
      converged = true;
    old_p = p; 
    ++iter;
  }

  trans = eigt::getTransFromPose(p);


  return true;
}

Scanmatcher::pose_t Scanmatcher::getTransformation() const{
  pose_t pose;
  double angle = eigt::getAngle(transform_);
  pose<<transform_.translation(),angle;
  return pose;
}

tf::Transform Scanmatcher::getTFTransform()const{
  tf::Transform msg;
  transform_t t = eigt::transBtwPoses(pose_,last_odom_);
  tf::Quaternion orientation;
  double angle = eigt::getAngle(t);
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

ml_ndt_scanmatching::NDTMapMsg Scanmatcher::getLayerData(size_t layer_id) const{
  if(layer_id > layer_.size())
    throw new std::invalid_argument("Showing data for this layer id is not possible. Check call of Scanmatcher::getLayerData(layer_id)");
  ml_ndt_scanmatching::NDTMapMsg msg;
  msg = layer_[layer_id].getLayerData();
  msg.x_cen = pose_(0);
  msg.y_cen = pose_(1);

  return msg;
}

void Scanmatcher::setResolution(const size_t res){
  resolution_ = res;
}
void Scanmatcher::setLayers(const size_t layers){
  layers_count_ = layers;
}

void Scanmatcher::setMaxRange(const double range){
  max_range_ = range;
}

// ****************** PRIVATE FUNCTIONS *********************************
Scanmatcher::points2_t Scanmatcher::projectPointsTo2D(const pcl_t &points) {
  points2_t points2d;
  point_t one_point;
  for (auto pt : points) {
    one_point << pt.x,pt.y;
    points2d.emplace_back(one_point);
  }
  return points2d;
}

void Scanmatcher::createLayers(){
transform_t offset;
offset.setIdentity();
//offset(1,2)-=1.0;
  for (size_t i = 0; i < layers_count_; ++i) {
    size_t power = static_cast<size_t>(pow(2,i));
    DEBUG("Creating Layer from create Layer");
    layer_.emplace_back(Layer(&points_,resolution_ * power, max_range_,offset));
  }
}

void Scanmatcher::initializeNdt(const pose_t &pose,const points2_t & points) {
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

void Scanmatcher::initializeNdt(const pose_t &pose, points2_t && points) {
  layer_.clear();
  DEBUG("ML-NDT: Layer created with2"<<points_.size()<<"pts");
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

bool Scanmatcher::calculateNdt(const pose_t & current_pose_odom, const points2_t &points) {
  if (!initialized_)
    return false;
  transform_t transform = eigt::transBtwPoses(last_odom_,current_pose_odom);
  DEBUG("odom trans:\n"<<transform.matrix());
  double angle = eigt::getAngle(transform);
  double distance = eigt::getDisplacement(transform);
  DEBUG("angle:  "<<angle);
  DEBUG("distance: "<<distance);
  transform.setIdentity();
  if(std::abs(angle)< 0.10 && distance < 0.5)
   return false;
   for (auto &l : layer_) {
     if (l.calculateNdt(transform, points)) {
       transform = l.getTransformation();
     }
  }
  transform_ = std::move(transform);
  DEBUG("tranform calculated \n"<<transform_.matrix());
  pose_t transformed_pose =eigt::transformPose(pose_,transform_);
  DEBUG("calculated transformed pose \n"<<transformed_pose);
  DEBUG("odometry pose \n"<<current_pose_odom);
  DEBUG("Diff odom<->calc_odom:\n"<<current_pose_odom - transformed_pose);
  updateLayers(transformed_pose,current_pose_odom, points);
  return true;
}

double Scanmatcher::scoreLayers(const transform_t & trans, const points2_t & cloud_in)const{
  double score  = 0;
  for (auto &l : layer_) {
     score += l.scoreLayer(trans,cloud_in);
  }
  return score;
}
