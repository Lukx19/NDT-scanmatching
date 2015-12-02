#include <layer.h>
// ********************* PUBLIC FUNCTIONS *******************************
Layer::point_t Layer::getPoint(const Id_t id) const { 
  DEBUG("points size:"<<points_->size());
  return points_->at(id); }

bool Layer::calculateNdt(transform_t &transf, points_t &points) {

  //printLaserPoints(*points_);
  //printLaserPoints(points);

  transform_t p = transf;
  transform_t delta;
  size_t max_iteration = 5;
  //DEBUG("Solving layer "<<size_);
  float error = 0;
  for (size_t iter = 0; iter < MAX_ITER; ++iter) {
    //DEBUG("Starting new iteration");
    pose_t g = pose_t::Zero();
    Eigen::Matrix3f hessian = hessian.Zero();
    Eigen::Matrix<float, 2, 3> jacobian;
    float si = p.rotation()(0,1);
    float co = p.rotation()(0,0);
    float error = 0;

    for (const auto & point : points) {
      point_t trans_point = p*point;//transformPoint(point, p);
      //DEBUG("Point: "<<point<<" transformed to:"<<trans_point);
      if(!isInBoundries(trans_point))
        continue;
      std::pair<size_t, size_t> new_coords =
          getFieldCoordintes(trans_point);
      if(fields_[new_coords.first][new_coords.second].getPoints() < MIN_POINTS_IN_FIELD)
        continue;
      //DEBUG("Point: "<<point<<" transformed to:"<<trans_point);
      Eigen::Matrix2f inv_variace =
          fields_[new_coords.first][new_coords.second].calcInvertedVariance();
      //DEBUG("Variance: "<<fields_[new_coords.first][new_coords.second].calcVariance());
      //DEBUG("Variance inv: "<<inv_variace);
      point_t difference =
          trans_point -
          fields_[new_coords.first][new_coords.second].calcMean();
      float score_part =
          expf((difference.transpose() * inv_variace * difference)(0) * -0.5);
      error += score_part;
      // prepare for calculating H * delta_p = -g where delta_p is change in
      // transformation H is hessian and g is gradient
      // incrementaly summing parts of  g for each point
      jacobian << 1, 0, -difference(0) * si - difference(1) * co, 
                  0, 1,  difference(0) * co - difference(1) * si;
      for (long r = 0; r < 3; ++r) {
        g(r) += (difference.transpose() * inv_variace * jacobian.col(r))(0) *
                score_part;
      }
      // create hessian for each point
      for (long r = 0; r < 3; ++r) {
        for (long c = 0; c < 3; ++c) {
          float jacc_r =
              (difference.transpose() * inv_variace * jacobian.col(r))(0);
          float jacc_c =
              -1*(difference.transpose() * inv_variace * jacobian.col(c))(0);
          point_t hess;
          if (r == 2 && c == 2)
            hess << -difference(0) * co + difference(1) * si,
                    -difference(0) * si - difference(1) * co;
          else
            hess << 0, 0;
          hessian(r, c) +=
              score_part *
              (jacc_r * jacc_c + difference.transpose() * inv_variace * hess +
               jacobian.col(c).transpose() * inv_variace * jacobian.col(r));
        }
      }
    }
    DEBUG("score:"<<error);
    pose_t delta_p;
    hessian = makeToSPD(hessian,g);

    //Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> saes(hessian);
    //DEBUG("eigenvalues: "<<saes.eigenvalues().transpose());

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
    DEBUG("delta"<<delta_p);
    p =p*delta.fromPositionOrientationScale(delta_p.head<2>(),
                                            Eigen::Rotation2Df(delta_p(2)),
                                            point_t::Ones());//delta_p;
    //DEBUG("Calculation of iteration done"<<delta_p);
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
}

bool Layer::isInBoundries(const point_t &point) const {
  if (std::abs(point(0)) < max_range_ && std::abs(point(1)) < max_range_)
    return true;
  else
    return false;
}

// Layer::point_t Layer::transformPoint(const Id_t id,
//                                      const pose_t &transform) const {
//   return std::move(transformPoint(points_->at(id), transform));
// }

// Layer::point_t Layer::transformPoint(const point_t &point,
//                                      const pose_t &transform) const {
//   float si = sinf(transform_(2));
//   float co = cosf(transform_(2));
//   Eigen::Matrix2f rotate;
//   rotate << co, -si, si, co;
//   return rotate * point + transform.segment(0, 2);
// }

Eigen::Matrix3f Layer::makeToSPD(const Eigen::Matrix3f &hess,const Eigen::Vector3f & s_gradient)const {
  Eigen::Matrix3f increase = Eigen::Matrix3f::Zero();
  //Eigen::LDLT<Eigen::Matrix3f> ldtl(hess);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> saes(hess);
  Eigen::Vector3f eigenvalues=saes.eigenvalues();
  float min = eigenvalues.minCoeff();
  float max = eigenvalues.maxCoeff();
  if(min < 0){
    Eigen::Matrix3f eigenvectors = saes.eigenvectors();
    float additive = s_gradient.norm();
    if(min + additive <= 0){
      additive =  0.001*max - min;
    }
    Eigen::Vector3f addition;
    addition<<additive,additive,additive;
    eigenvalues += addition;
    return eigenvectors*eigenvalues.asDiagonal()*eigenvectors.transpose();
  }
  //DEBUG("ML_NDT: eigen_valuses of hessian:"<<saes.eigenvalues());
  return hess;
  // float k = 0;
  // while (!ldtl.isPositive()) {
  //   float min_eigen = saes.eigenvalues().minCoeff();
  //   increase += Eigen::Matrix3f::Identity() * (powf(k, 2) * min_eigen);
  //   ++k;
  //   ldtl.compute(mat + increase);
  //   saes.compute(mat + increase);
  // }
  // return std::move(mat + increase);
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
      static_cast<size_t>(std::floor((point(0)+max_range_) / resolution));
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