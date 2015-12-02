#include <field.h>

void Field::addPoint(Id_t id) { points_ids_.emplace_back(id); }

Field::point_t Field::calcMean() const {
  if (points_ids_.size() < 3)
    return point_t::Zero();
  Eigen::Vector2f mean = Eigen::Vector2f::Zero();
  for (const auto & id : points_ids_) {
    mean += points_->at(id);
  }
  return mean / points_ids_.size();
}

Field::var_t Field::calcVariance() const {
  if (points_ids_.size() < 3)
    return var_t::Zero();

  var_t var = Eigen::Matrix2f::Zero();
  point_t mean = calcMean();
  for (const auto & id : points_ids_) {
    point_t pt = points_->at(id);
    var += (pt - mean) * (pt - mean).transpose();
  }
  return var / (points_ids_.size() - 1);
}

Field::var_t Field::calcInvertedVariance()const{
  const var_t var = calcVariance();
  if(var.determinant() < 0.0001){
    return pinv<var_t>(var,1.e-06f);
  }else{
    return var.inverse();
  }
}

size_t Field::getPoints() const{
  return points_ids_.size();
}

// ******************** PRIVATE FUNCTIONS ***************



