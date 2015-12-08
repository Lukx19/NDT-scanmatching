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
  Eigen::MatrixXf mp;
  point_t mean = calcMean();
  mp.resize(points_ids_.size(),2);
  for (size_t i=0; i<points_ids_.size();++i){
    point_t pt = points_->at(points_ids_[i]);
    mp(i,0) = pt(0) -mean(0);
    mp(i,1) = pt(1) -mean(1); 
  }
  return (mp.transpose()*mp) / (points_ids_.size()-1);
  
  // var_t var = Eigen::Matrix2f::Zero();
  // point_t mean = calcMean();
  // var_t sxx = var_t::Zero();
  // point_t sx = point_t::Zero();
  // for (const auto & id : points_ids_) {
  //   point_t pt = points_->at(id);
  //   sx+=pt;
  //   sxx+=pt*pt.transpose();
  //   //var += (pt - mean) * (pt - mean).transpose();
  // }
  //DEBUG(var / (points_ids_.size() - 1));
  //DEBUG(((sxx - 2 *(sx* mean.transpose()))/points_ids_.size()) + mean * mean.transpose());
  //return var / (points_ids_.size() - 1);
  //return ((sxx - 2 *(sx* mean.transpose()))/points_ids_.size()) + mean * mean.transpose();
}

Field::var_t Field::calcInvertedVariance()const{
  var_t covar = calcVariance();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> solver (covar);
  Eigen::Matrix2f evecs = solver.eigenvectors().real();
  Eigen::Vector2f evals = solver.eigenvalues().real();
  float max_eval = evals.maxCoeff();
  for(size_t i=0;i<2;++i){
    if(evals(i)< max_eval* EVAL_FACTOR)
      evals(i) = max_eval / EVAL_FACTOR;
  }
  return evecs*(evals.asDiagonal().inverse())*(evecs.transpose());
  // if (solver.eigenvalues ()[0] < 0.0001F * solver.eigenvalues ()[1])
  // {
  //   Eigen::Matrix2f v = solver.eigenvalues ().asDiagonal ();
  //   Eigen::Matrix2f q = solver.eigenvectors ();
  //   // set minimum smallest eigenvalue:
  //   v (0,0) = v (1,1) * 0.0001F;
  //   var = q * v * q.transpose();
  // }
  // return var.inverse();

  // if(var.determinant() < 0.0001){
  //   return pinv<var_t>(var,1.e-06f);
  // }else{
  //   return var.inverse();
  // }
}

size_t Field::getPoints() const{
  return points_ids_.size();
}

// ******************** PRIVATE FUNCTIONS ***************



