#include <layer.h>

Layer::point_t Layer::getPoint(const Id_t id) const { return points_->at(id); }

bool Layer::calculateNdt(pose_t &transf, points_t &points) {

  pose_t p = transf;
  size_t max_iteration = 4;

  for (size_t iter = 0; iter < max_iteration; ++iter) {

    pose_t g = pose_t::Zero();
    Eigen::Matrix3f hessian = hessian.Zero();
    Eigen::Matrix<float, 2, 3> jacobian;
    float si = sinf(p(2));
    float co = cosf(p(2));

    for (auto point : points) {
      point_t transform_point = transformPoint(point, p);
      std::pair<size_t, size_t> new_coords =
          getFieldCoordintes(transform_point);

      Eigen::Matrix2f inv_variace =
          fields_[new_coords.first][new_coords.second].calcVariance().inverse();
      point_t diffrence =
          transform_point -
          fields_[new_coords.first][new_coords.second].calcMean();
      float score_part =
          expf((-diffrence.transpose() * inv_variace * diffrence)(0) / 2);
      // prepare for calculating H * delta_p = -g where delta_p is change in
      // transformation
      // incrementaly summing parts of  g for each point
      jacobian << 1, 0, -point(0) * si - point(1) * co, 0, 1,
          point(0) * co - point(1) * si;
      for (long r = 0; r < 3; ++r) {
        g(r) += (diffrence.transpose() * inv_variace * jacobian.col(r))(0) *
                score_part;
      }

      // create hessian
      for (long r = 0; r < 3; ++r) {
        for (long c = 0; c < 3; ++c) {
          float jacc_r =
              (diffrence.transpose() * inv_variace * jacobian.col(r))(0);
          float jacc_c =
              (diffrence.transpose() * inv_variace * jacobian.col(c))(0);
          point_t hess;
          if (r == 2 && c == 2)
            hess << -point(0) * co + point(1) * si,
                -point(0) * si - point(1) * co;
          else
            hess << 0, 0;
          hessian(r, c) +=
              score_part *
              (jacc_r * jacc_c + diffrence.transpose() * inv_variace * hess +
               jacobian.col(c).transpose() * inv_variace * jacobian.col(r));
        }
      }
    }

    hessian = makeToSPD(hessian);
    Eigen::LDLT<Eigen::Matrix3f> ldtl(hessian);
    if (ldtl.info() != Eigen::Success)
      throw std::logic_error("Error in computation of hessian");
    pose_t delta_p = ldtl.solve(-g);
    if (ldtl.info() != Eigen::Success)
      throw std::logic_error(
          "Error in computation of diffrence in transformation");
    p += delta_p;
  }

  transform_ = p;
  // if(hessian.)
  return true;
}

Layer::pose_t Layer::getTransformation(){
  return transform_;
}

//*********************************** PRIVATE FUNCTIONS *********************

void Layer::initializeFields() {
  for (size_t row = 0; row < size_; ++row) {
    field_line_t line;
    for (size_t col = 0; col < size_; ++col) {
      line.push_back(Field(this));
    }
    fields_.push_back(std::move(line));
  }

  size_t id = 0;

  for (const auto point : (*points_)) {
    if (isInBoundries(point)) {
      std::pair<size_t, size_t> coord = getFieldCoordintes(point);
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

Layer::point_t Layer::transformPoint(const Id_t id,
                                     const pose_t &transform) const {
  return std::move(transformPoint(points_->at(id), transform));
}

Layer::point_t Layer::transformPoint(const point_t &point,
                                     const pose_t &transform) const {
  float si = sinf(transform_(2));
  float co = cosf(transform_(2));
  Eigen::Matrix2f rotate;
  rotate << co, -si, si, co;
  return rotate * point + transform.segment(0, 2);
}

Eigen::Matrix3f Layer::makeToSPD(Eigen::Matrix3f &mat) {
  Eigen::Matrix3f increase = Eigen::Matrix3f::Zero();
  Eigen::LDLT<Eigen::Matrix3f> ldtl(mat);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> saes(mat,
                                                      Eigen::EigenvaluesOnly);
  float k = 0;
  while (!ldtl.isPositive()) {
    float min_eigen = saes.eigenvalues().minCoeff();
    increase += Eigen::Matrix3f::Identity() * (powf(k, 2) * min_eigen);
    ++k;
    ldtl.compute(mat + increase);
    saes.compute(mat + increase);
  }
  return std::move(mat + increase);
}
std::pair<size_t, size_t>
Layer::getFieldCoordintes(const point_t &point) const {
  // how many  meters of point locations are in one field
  float resolution = max_range_ / size_;
  // move space of points from [-range,range] to [0,2range] in x axis
  size_t fieldx =
      static_cast<size_t>(std::floor((point(0) + max_range_) / resolution));
  // flips space of points in y axis from [0,range] to [range,0]
  size_t fieldy =
      static_cast<size_t>(std::floor((max_range_ - point(1)) / resolution));
  return std::make_pair(fieldx, fieldy);
}


