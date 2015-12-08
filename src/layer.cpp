#include <layer.h>
// ********************* PUBLIC FUNCTIONS *******************************
Layer::point_t Layer::getPoint(const Id_t id) const { 
  DEBUG("points size:"<<points_->size());
  return points_->at(id); }

bool Layer::calculateNdt(transform_t &transf, points_t &points) {

  float lfc1,lfc2,lfd3, lfd1,lfd2;
  float integral, outlier_ratio, support_size;
  integral = 0.1f;
  outlier_ratio = 0.35f;
  support_size = size_;
  lfc1 = (1-outlier_ratio)/integral;
  lfc2 = outlier_ratio/powf(support_size,2);
  lfd3 = -logf(lfc2);
  lfd1 = -logf( lfc1 + lfc2 ) - lfd3;
  lfd2 = -logf((-logf( lfc1 * expf( -0.5f ) + lfc2 ) - lfd3 ) / lfd1);



  //printLaserPoints(*points_);
  //printLaserPoints(points);
  transform_t p = transf;
  //p.setIdentity();
  transform_t old_p = transf;
  transform_t delta;
  point_t difference;
  float best_score=0;
  size_t iter=0;
  bool converged = false;
  //DEBUG("Solving layer "<<size_);
  while (!converged) {
    //DEBUG("Starting new iteration");
    old_p = p;
    pose_t g = pose_t::Zero();
    Eigen::Matrix3f hessian = hessian.Zero();
    Eigen::Matrix<float, 2, 3> jacobian;

    Eigen::Rotation2D<float> rot(0);
    rot = rot.fromRotationMatrix(p.rotation());
    float si = std::sin(rot.angle());
    float co = std::cos(rot.angle());
    double total_score = 0;

    for (const auto & point : points) {
      point_t trans_point = p*point;
      float x = trans_point(0);
      float y = trans_point(1);
      Field field;
      if(!getPointField(trans_point,field))
        continue;
      Eigen::Matrix2f inv_variace =field.calcInvertedVariance();
      //DEBUG(inv_variace);
      difference = trans_point - field.calcMean();
      double point_score =lfd1*std::exp((double)(difference.dot(field.calcInvertedVariance() * difference) * -0.5*lfd2));
      //DEBUG(difference.dot(field.calcInvertedVariance() * difference)* -0.5F);
      //DEBUG(point_score);
      total_score += point_score;
      // prepare for calculating H * delta_p = -g where delta_p is change in
      // transformation H is hessian and g is gradient
      // incrementaly summing parts of  g for each point
      jacobian << 1, 0, -x * si - y * co, 
                  0, 1,  x * co - y * si;
      for (long r = 0; r < 3; ++r) {
        g(r) += (difference.dot(inv_variace * jacobian.col(r)))*point_score;
      }
      // create hessian for each point
      for (long r = 0; r < 3; ++r) {
        for (long c = 0; c < 3; ++c) {
          float jacc_r =
              difference.dot(inv_variace * jacobian.col(r));
          float jacc_c =
              -1*(difference.dot(inv_variace * jacobian.col(c)));
          point_t hess;
          if (r == 2 && c == 2)
            hess << -x * co + y * si,
                    -x * si - y * co;
          else
            hess << 0, 0;
          hessian(r, c) +=
              point_score *
              (jacc_r * jacc_c + difference.transpose() * inv_variace * hess +
               jacobian.col(c).transpose() * inv_variace * jacobian.col(r));
        }
      }
    }
    if(total_score  > -0.0001){
     DEBUG("Score is too small finished.");
     break;
    }

    //DEBUG("score:"<<total_score);
    pose_t delta_p;
    hessian = makeToSPD(hessian,g);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> saes(hessian);
    //DEBUG("eigenvalues: "<<saes.eigenvalues().transpose());
    float min_eigenval = saes.eigenvalues().minCoeff();
    if(min_eigenval<0 || saes.eigenvalues().sum() * 0 != 0){
      DEBUG("EIGEN values are negative");
      break;
    }
    Eigen::LDLT<Eigen::Matrix3f> ldtl(hessian);
    if (ldtl.info() != Eigen::Success){
      delta_p << 0,0,0;
      ROS_ERROR_STREAM("Error in computation of hessian");
      DEBUG("Error in computation of hessian");
      //throw std::logic_error("Error in computation of hessian");
    }else {
      delta_p = ldtl.solve(-g);
      if (ldtl.info() != Eigen::Success){
        ROS_ERROR_STREAM("Error in computation of diffrence in transformation");
        DEBUG("Error in computation of diffrence in transformation");
        delta_p <<0,0,0;
        //throw std::logic_error(
        //    "Error in computation of diffrence in transformation");
      }
    }
    delta = delta.fromPositionOrientationScale(delta_p.head<2>(),
                                              Eigen::Rotation2Df(delta_p(2)),
                                              point_t::Ones());//delta_p;
      p=delta * p;
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

bool Layer::isInBoundries(const point_t &point) const {
  if (std::abs(point(0)) < max_range_ && std::abs(point(1)) < max_range_)
    return true;
  else
    return false;
}

Eigen::Matrix3f Layer::makeToSPD(const Eigen::Matrix3f &hess,const Eigen::Vector3f & s_gradient)const {
  //Eigen::LDLT<Eigen::Matrix3f> ldtl(hess);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> saes(hess);
  Eigen::Vector3f eigenvalues=saes.eigenvalues();
  float min = eigenvalues.minCoeff();
  float max = eigenvalues.maxCoeff();
  if(min < 0){
    Eigen::Matrix3f eigenvectors = saes.eigenvectors();
    float additive = s_gradient.norm();
    if(min + additive <= 0){
      additive =  0.001f *max - min;
    }
    Eigen::Vector3f addition;
    addition<<additive,additive,additive;
    eigenvalues += addition;
    return eigenvectors*eigenvalues.asDiagonal()*eigenvectors.transpose();
  }
  //DEBUG("ML_NDT: eigen_valuses of hessian:"<<saes.eigenvalues());
  return hess;
}

Eigen::Matrix3f Layer::makeToSPD2(const Eigen::Matrix3f & hessian,const Eigen::Vector3f & gradient)const {

  Eigen::EigenSolver<Eigen::Matrix3f> solver;
  solver.compute (hessian, false);
  float min_eigenval = solver.eigenvalues().real().minCoeff();
  if(min_eigenval<0){
    float lambda = 1.1F*min_eigenval -1.0F;
    Eigen::Matrix3f addition;
    addition = -lambda * Eigen::Matrix3f::Identity(); 
    return hessian + addition;
  }
  return hessian;
}


/*
  Brief: Calculates coordinates in layer's grid from point in parameter
*/
std::pair<size_t, size_t>
Layer::getFieldCoordintes(const point_t &point) const {
  // how many  meters of point locations are in one field
  float resolution = 2*max_range_ / size_;
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