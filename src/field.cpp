#include <field.h>

void Field::addPoint(Id_t id) { points_ids_.emplace_back(id); }

void Field::calcNormDistribution(){
  prepNormDist();
}

Field::point_t Field::getMean() const{
  if(is_calc_){
    return mean_;
  }
  else{
    return point_t::Zero();
  }
}

Field::var_t Field::getCovar() const{
    if(is_calc_){
    return covar_;
  }
  else{
    return var_t::Zero();
  }
}

Field::var_t Field::getInvCovar() const{
  if(is_calc_){
    return inv_covar_;
  }
  else{
    return var_t::Zero();
  }
}

size_t Field::getPoints() const{
  return points_ids_.size();
}

bool Field::isReady() const{
  return is_calc_;
}

std::string Field::toString() const{
  //return std::to_string(points_ids_.size())+" ";//+":"+mean_[1];
  if(points_ids_.size() < MIN_PTR_EVAL)
    return ". ";
  else
    return "O ";
}

// ******************** PRIVATE FUNCTIONS ***************
void Field::prepNormDist(){
  mean_ = calcMean();
  covar_ = calcCovariance(mean_);
  inv_covar_ = calcInvertedCovariance(covar_);
  is_calc_ = true;
}

Field::point_t Field::calcMean() const {
  if (points_ids_.size() < MIN_PTR_EVAL)
    return point_t::Zero();
  Eigen::Vector2d mean = Eigen::Vector2d::Zero();
  for (const auto & id : points_ids_) {
    mean += points_->at(id);
  }
  //DEBUG("mean: "<<mean / points_ids_.size());
  return mean / points_ids_.size();
}

Field::var_t Field::calcCovariance(const point_t & mean) const {
  if (points_ids_.size() < MIN_PTR_EVAL)
    return var_t::Zero();
  Eigen::MatrixXd mp;
  mp.resize(points_ids_.size(),2);
  for (size_t i=0; i<points_ids_.size();++i){
    point_t pt = points_->at(points_ids_[i]);
    mp(i,0) = pt(0) -mean(0);
    mp(i,1) = pt(1) -mean(1); 
  }
  //DEBUG("var: "<<(mp.transpose()*mp) / (points_ids_.size()-1))
  return (mp.transpose()*mp) / (points_ids_.size()-1);
}

Field::var_t Field::calcInvertedCovariance(const var_t & covar)const{
  if (points_ids_.size() < MIN_PTR_EVAL)
    return var_t::Zero();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver (covar);
  Eigen::Matrix2d evecs = solver.eigenvectors().real();
  Eigen::Vector2d evals = solver.eigenvalues().real();
  float max_eval = evals.maxCoeff();
  for(size_t i=0;i<2;++i){
    if(evals(i)< max_eval* EVAL_FACTOR)
      evals(i) = max_eval / EVAL_FACTOR;
  }
  //DEBUG("covar: "<<evecs*(evals.asDiagonal().inverse())*(evecs.transpose()));
  return evecs*(evals.asDiagonal().inverse())*(evecs.transpose());
}

