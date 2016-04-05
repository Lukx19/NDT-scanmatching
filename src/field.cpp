#include <ml_ndt_scanmatching/other/field.h>

void Field::addPoint(Id_t id) { points_ids_.emplace_back(id); }

void Field::calcNormDistribution(){
  prepNormDist();
}

ml_ndt_scanmatching::NDTCellMsg Field::getCellData()const
{
  ml_ndt_scanmatching::NDTCellMsg msg;
  if(isReady()){
    msg.mean_x = mean_(0);
    msg.mean_y = mean_(1);
    msg.mean_z = 0.0;
    msg.occupancy = 100;
    msg.N = static_cast<long>(points_ids_.size());
    // filling from our covar 2x2 to 3d ndt msg covar 3x3
    for(long row = 0; row < 2;++row){
      for(long col = 0;col<2;++col){
        msg.cov_matrix.push_back(covar_(row,col));
      }
      msg.cov_matrix.push_back(0.0);
    }
    for(size_t i =0; i< 3; ++i){
      msg.cov_matrix.push_back(0.0);
    }
  }else{
    msg.mean_x = 0.0;
    msg.mean_y = 0.0;
    msg.mean_z = 0.0;
    msg.occupancy = 0.0;
    msg.N = static_cast<long>(points_ids_.size());
    for(size_t i = 0; i< 9;++i){
      msg.cov_matrix.push_back(0.0);
    }
    msg.cov_matrix[0] = 1;
    msg.cov_matrix[4] = 1;
    msg.cov_matrix[8] = 1; 
  }
  
  return msg;
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
  if(!is_calc_)
    return ". ";
  else
    return "O ";
}

// ******************** PRIVATE FUNCTIONS ***************
void Field::prepNormDist(){
  mean_ = calcMean();
  auto res = calcInvertedCovariance(calcCovariance(mean_));
  covar_=std::get<0>(res);
  inv_covar_ = std::get<1>(res);
  is_calc_ = std::get<2>(res);
  if(!is_calc_){
    mean_ = point_t::Zero();
  }
}

Field::point_t Field::calcMean() const {
  if (points_ids_.size() < MIN_PTR_EVAL)
    return point_t::Zero();
  Eigen::Vector2d mean = Eigen::Vector2d::Zero();
  for (const auto & id : points_ids_) {
    mean += points_->at(id);
  }
  DEBUG_FIE("mean: "<<(mean / points_ids_.size()).transpose());
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
  DEBUG_FIE("covar original: \n"<<(mp.transpose()*mp) / (points_ids_.size()-1))
  return (mp.transpose()*mp) / (points_ids_.size()-1);
}

std::tuple<Field::var_t,Field::var_t,bool>
Field::calcInvertedCovariance(const var_t & covar)const{
  std::tuple<Field::var_t,Field::var_t,bool> ret;
  if (points_ids_.size() < MIN_PTR_EVAL){
    return std::make_tuple(var_t::Zero(),var_t::Zero(),false);
  }
  std::pair<Field::var_t,Field::var_t> pair;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver (covar);
  Eigen::Matrix2d evecs = solver.eigenvectors().real();
  Eigen::Vector2d evals = solver.eigenvalues().real();
  DEBUG_FIE("eigenvalues befor scaling: "<<evals.transpose());
  DEBUG_FIE("determinant before scaling: "<<covar.determinant());
  if(evals(0) <= 0 || evals(1) <= 0){
    return std::make_tuple(var_t::Zero(),var_t::Zero(),false);
  }
  double max_eval = evals.maxCoeff();
  bool recalc = false;
  for(long i=0;i<2;++i){
    if(evals(i) * EVAL_FACTOR < max_eval){
      evals(i) = max_eval / EVAL_FACTOR;
      recalc = true;
    }
  }
  DEBUG_FIE("inv_covar: \n "<<evecs*(evals.asDiagonal().inverse())*(evecs.transpose()));
  DEBUG_FIE("covar after scale: \n"<<evecs*(evals.asDiagonal())*(evecs.transpose()));

  Field::var_t new_covar = covar;
  if(recalc){
   new_covar = evecs*(evals.asDiagonal())*(evecs.transpose());
 }
 DEBUG_FIE("covar determinat: "<<new_covar.determinant()<<"\n");
  return std::make_tuple(new_covar,evecs*(evals.asDiagonal().inverse())*(evecs.transpose()),true);
}

