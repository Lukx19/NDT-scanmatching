#include <field.h>

void Field::addPoint(Id_t id) { points_ids_.emplace_back(id); }

Field::point_t Field::calcMean() const {
  if (points_ids_.size() < 3)
    return point_t::Zero();

  Eigen::Vector2f mean = Eigen::Vector2f::Zero();
  for (auto id : points_ids_) {
    mean += layer_->getPoint(id);
  }
  return mean / points_ids_.size();
}

Field::var_t Field::calcVariance() const {
  if (points_ids_.size() < 3)
    return var_t::Zero();

  Eigen::Matrix2f var = Eigen::Matrix2f::Zero();
  point_t mean = calcMean();
  for (auto id : points_ids_) {
    point_t pt = layer_->getPoint(id);
    var += (pt - mean) * (pt - mean).transpose();
  }
  return var / (points_ids_.size() - 1);
}
