#include <layer.h>
// ********************* PUBLIC FUNCTIONS *******************************
Layer::point_t Layer::getPoint(const Id_t id) const { 
  DEBUG("points size:"<<points_->size());
  return points_->at(id); }

bool Layer::calculateNdt(transform_t &transf, points_t &points) {

  //printLaserPoints(*points_);
  //printLaserPoints(points);
  transform_t p = transf;
  //p.setIdentity();
  transform_t old_p = transf;
  transform_t delta;
  point_t difference;
  double best_score=0;
  size_t iter=0;
  bool converged = false;
  //DEBUG("Solving layer "<<size_);
  while (!converged) {
    //DEBUG("Starting new iteration");
    old_p = p;
    pose_t g = pose_t::Zero();
    hessian_t hessian = hessian.Zero();
    Eigen::Matrix<double, 2, 3> jacobian;

    Eigen::Rotation2D<double> rot(0);
    rot = rot.fromRotationMatrix(p.rotation());
    double si = std::sin(rot.angle());
    double co = std::cos(rot.angle());
    double total_score = 0;

    for (const auto & point : points) {
      point_t trans_point = p*point;
      double x = trans_point(0);
      double y = trans_point(1);
      Field field;
      if(!getPointField(trans_point,field))
        continue;
      covar_t inv_covar =field.getInvCovar();
      //DEBUG(inv_variace);
      difference = trans_point - field.getMean();
      double point_score = scorePoint(field,trans_point);
      //DEBUG(difference.dot(field.calcInvertedVariance() * difference)* -0.5F);
      //DEBUG(point_score);
      total_score += point_score;
      // prepare for calculating H * delta_p = -g where delta_p is change in
      // transformation H is hessian and g is gradient
      // incrementaly summing parts of  g for each point
      jacobian << 1, 0, -x * si - y * co, 
                  0, 1,  x * co - y * si;
      point_t hess_derivative;
      hess_derivative << -x * co + y * si,
                         -x * si - y * co;
      g += pointGradient(difference,inv_covar,point_score,jacobian);
      hessian+=pointHessian(difference,inv_covar,point_score,jacobian,hess_derivative);
      // create hessian for each point
      
    }
    if(total_score  > -0.0001){
     DEBUG("Score is too small finished.");
     break;
    }

    //DEBUG("score:"<<total_score);
    pose_t delta_p;
    hessian = makeToSPD(hessian,g);

    Eigen::SelfAdjointEigenSolver<hessian_t> saes(hessian);
    //DEBUG("eigenvalues: "<<saes.eigenvalues().transpose());
    double min_eigenval = saes.eigenvalues().minCoeff();
    if(min_eigenval<0 || saes.eigenvalues().sum() * 0 != 0){
      DEBUG("EIGEN values are negative");
      break;
    }
    Eigen::LDLT<hessian_t> ldtl(hessian);
    if (ldtl.info() != Eigen::Success){
      delta_p.setZero();
      ROS_ERROR_STREAM("Error in computation of hessian");
      DEBUG("Error in computation of hessian");
      //throw std::logic_error("Error in computation of hessian");
    }else {
      delta_p = ldtl.solve(-g);
      if (ldtl.info() != Eigen::Success){
        ROS_ERROR_STREAM("Error in computation of diffrence in transformation");
        DEBUG("Error in computation of diffrence in transformation");
        delta_p.setZero();
        //throw std::logic_error(
        //    "Error in computation of diffrence in transformation");
      }
    }
    delta = delta.fromPositionOrientationScale(delta_p.head<2>(),
                                              Eigen::Rotation2Dd(delta_p(2)),
                                              point_t::Ones());//delta_p;
      p=delta * p;
      int step = lineSearchMT(p,g,delta_p,points);
      DEBUG("step: "<<step);
   if(iter >= MAX_ITER || delta_p.array().abs().sum() < 0.001)
     converged = true;
    ++iter;
  }

  transform_ = p;
  return true;
}

Layer::transform_t Layer::getTransformation(){
  return transform_;
}

//*********************************** PRIVATE FUNCTIONS *********************

void Layer::initializeFields(points_t * points) {
 for (size_t row = 0; row < size_; ++row) {
    field_line_t line;
    for (size_t col = 0; col < size_; ++col) {
      line.push_back(Field(points));
    }
    fields_.push_back(std::move(line));
  }
  //DEBUG("ML_NDT:Fields created for layer"<<size_<<"created "<<fields_.size());
  size_t id = 0;

  for (const auto & point : (*points_)) {
    if (isInBoundries(point)) {
      std::pair<size_t, size_t> coord = getFieldCoordintes(point);
      // DEBUG("ML_NDT: point:"<<point<< "coordinates"<<coord.first<<" "<<coord.second);
      fields_[coord.first][coord.second].addPoint(id);
    }
    ++id;
  }

  // for (size_t row = 0; row < size_; ++row) {
  //   field_line_t line;
  //   for (size_t col = 0; col < size_; ++col) {
  //     if(fields_[row][col].getPoints()>MIN_POINTS_IN_FIELD)
  //       std::cout<<"D";
  //     else
  //       std::cout<<".";
  //   }
  //   DEBUG("");
  // }

}

void Layer::initializeParams(){
  double lfc1,lfc2,lfd3;
  double integral, outlier_ratio, support_size;
  integral = 0.1;
  outlier_ratio = 0.35;
  support_size = size_;
  lfc1 = (1-outlier_ratio)/integral;
  lfc2 = outlier_ratio/pow(support_size,2);
  lfd3 = -log(lfc2);
  LFD1 = -log( lfc1 + lfc2 ) - lfd3;
  LFD2 = -log((-log( lfc1 * exp( -0.5 ) + lfc2 ) - lfd3 ) / LFD1);
}

bool Layer::isInBoundries(const point_t &point) const {
  if (std::abs(point(0)) < max_range_ && std::abs(point(1)) < max_range_)
    return true;
  else
    return false;
}

Layer::hessian_t Layer::makeToSPD(const hessian_t &hess,const Eigen::Vector3d & s_gradient)const {
  //Eigen::LDLT<Eigen::Matrix3f> ldtl(hess);
  Eigen::SelfAdjointEigenSolver<hessian_t> saes(hess);
  Eigen::Vector3d eigenvalues=saes.eigenvalues();
  double min = eigenvalues.minCoeff();
  double max = eigenvalues.maxCoeff();
  if(min < 0){
    hessian_t eigenvectors = saes.eigenvectors();
    double additive = s_gradient.norm();
    if(min + additive <= 0){
      additive =  0.001 *max - min;
    }
    Eigen::Vector3d addition;
    addition<<additive,additive,additive;
    eigenvalues += addition;
    return eigenvectors*eigenvalues.asDiagonal()*eigenvectors.transpose();
  }
  //DEBUG("ML_NDT: eigen_valuses of hessian:"<<saes.eigenvalues());
  return hess;
}

Layer::hessian_t Layer::makeToSPD2(const hessian_t & hessian,const Eigen::Vector3d & gradient)const {

  Eigen::EigenSolver<hessian_t> solver;
  solver.compute (hessian, false);
  double min_eigenval = solver.eigenvalues().real().minCoeff();
  if(min_eigenval<0){
    double lambda = 1.1*min_eigenval -1.0;
    hessian_t addition;
    addition = -lambda * hessian_t::Identity(); 
    return hessian + addition;
  }
  return hessian;
}

Layer::pose_t Layer::pointGradient(const point_t & difference,
                                   const covar_t & inv_covar,
                                   double score,
                                   const Eigen::Matrix<double, 2, 3> & jacobian)const{
  pose_t g = pose_t::Zero();
  for (long r = 0; r < 3; ++r) {
    g(r) = (difference.dot(inv_covar * jacobian.col(r)))*score;
  }
  return g;
}

Layer::hessian_t Layer::pointHessian(const point_t & difference,
                                     const covar_t & inv_covar,
                                     double score,
                                     const Eigen::Matrix<double, 2, 3> & jacobian,
                                     const point_t & hessian_derivative)const{
  hessian_t hessian; 
  for (long r = 0; r < 3; ++r) {
    for (long c = 0; c < 3; ++c) {
      double jacc_r = difference.dot(inv_covar * jacobian.col(r));
      double jacc_c =-1*(difference.dot(inv_covar * jacobian.col(c)));
      point_t hess;
      if (r == 2 && c == 2)
        hess = hessian_derivative;
      else
        hess << 0, 0;
      hessian(r, c) = score *
                      (jacc_r * jacc_c + difference.transpose() * inv_covar * hess +
                      jacobian.col(c).transpose() * inv_covar * jacobian.col(r));
    }
  }
  return hessian;
}
/*
  Brief: Calculates coordinates in layer's grid from point in parameter
*/
std::pair<size_t, size_t>
Layer::getFieldCoordintes(const point_t &point) const {
  // how many  meters of point locations are in one field
  double resolution = 2*max_range_ / size_;
  // move space of points from [-range,range] to [0,2range] in x axis
  size_t fieldx = 
      static_cast<size_t>(std::floor((max_range_-point(0)) / resolution));
  // flips space of points in y axis from [0,range] to [range,0]
  size_t fieldy =
      static_cast<size_t>(std::floor((max_range_ - point(1)) / resolution));
  return std::make_pair(fieldx, fieldy);
}


void Layer::printLaserPoints(const points_t & points)const{
  DEBUG("points: \n");
  for(auto & pt : points){
    std::cout<<pt.transpose() << ";";
  }
  DEBUG("\n");
}

bool Layer::getPointField(const point_t & pt,Field & field)const {
  if(!isInBoundries(pt))
    return false;
  std::pair<size_t, size_t> coords = getFieldCoordintes(pt);
  if(fields_[coords.first][coords.second].getPoints() < MIN_POINTS_IN_FIELD)
    return false;
  field = fields_[coords.first][coords.second];
  return true;
}

double Layer::scorePoint(const Field & field, const point_t & trans_point)const{
  point_t difference = trans_point - field.getMean();
  return LFD1*std::exp((difference.dot(field.getInvCovar() * difference) * -0.5F*LFD2));
}

double Layer::scoreLayer(const transform_t & trans, const points_t & cloud_in) const {
   double total_score =0;
   for (const auto & point : cloud_in) {
      point_t trans_point = trans*point;
      Field field;
      if(!getPointField(trans_point,field))
        continue;
      total_score += scorePoint(field,trans_point);
    }
    return total_score;
}
