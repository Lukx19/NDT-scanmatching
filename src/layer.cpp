#include <layer.h>
// OPTIMALISATIONDATA **************************************************

void OptimalisationData::setZero()
{
  hessian.setZero();
  gradient.setZero();
  score = 0.0;
  is_valid = false;
}

OptimalisationData OptimalisationData::operator+=(const OptimalisationData & other)
{
  OptimalisationData data;
  data.hessian += other.hessian;
  data.gradient += other.gradient;
  score += other.score;
  return *this;
}

void OptimalisationData::makeToSPD()
{
  // Eigen::LDLT<Eigen::Matrix3f> ldtl(hess);
  Eigen::SelfAdjointEigenSolver<hessian_t> saes(hessian);
  Eigen::Vector3d eigenvalues = saes.eigenvalues();
  double min = eigenvalues.minCoeff();
  double max = eigenvalues.maxCoeff();
  if (min < 0) {
    hessian_t eigenvectors = saes.eigenvectors();
    double additive = gradient.norm();
    if (min + additive <= 0) {
      additive = 0.001 * max - min;
    }
    Eigen::Vector3d addition;
    addition << additive, additive, additive;
    eigenvalues += addition;
    hessian =
        eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
  }
  // DEBUG_LAY("ML_NDT: eigen_valuses of hessian:"<<saes.eigenvalues());
}

void OptimalisationData::makeToSPD2()
{
  Eigen::EigenSolver<hessian_t> solver;
  solver.compute(hessian, false);
  double min_eigenval = solver.eigenvalues().real().minCoeff();
  if (min_eigenval < 0) {
    double lambda = 1.1 * min_eigenval - 1.0;
    hessian_t addition;
    addition = -lambda * hessian_t::Identity();
    hessian += addition;
  }
}

// ********************* PUBLIC FUNCTIONS *******************************
Layer::point_t Layer::getPoint(const Id_t id) const { 
  DEBUG_LAY("points size:"<<points_->size());
  return points_->at(id); }

ml_ndt_scanmatching::NDTMapMsg Layer::getLayerData() const
{
  ml_ndt_scanmatching::NDTMapMsg msg;
  msg.x_size = 2 * max_range_;
  msg.y_size = 2 * max_range_;
  msg.z_size = 0.0;
  msg.x_cell_size = 2* max_range_ / size_;
  msg.y_cell_size = 2* max_range_ / size_;
  msg.z_cell_size = 0.0;
  msg.x_cen = 0;
  msg.y_cen = 0;
  msg.z_cen = 0;
  for(auto & grid_line : fields_){
    for(auto & field : grid_line){
      msg.cells.emplace_back(field.getCellData());
    }
  }
  return msg;
}
bool Layer::calculateNdt(const transform_t &transf, const points_t &points) {

  transform_t p = transf;
  //p.setIdentity();
  transform_t old_p = transf;
  transform_t best_p = transform_t::Identity();
  transform_t delta;

  point_t difference;
  long best_score=std::numeric_limits<long>::max();
  size_t iter=0;
  bool converged = false;
  //DEBUG_LAY("Solving layer "<<size_);
  DEBUG_LAY("Score at the beggining: "<<-scoreLayer(p,points));
  while (!converged) {
    DEBUG_LAY("Starting new iteration "<<iter <<" best score: "<<best_score);
    DEBUG_LAY("calc_trans before opt: "<<eigt::getPoseFromTransform(p).transpose());
    old_p = p;
    OptimalisationData opt_data = calcOptData(p,points);
    DEBUG_LAY("score:"<<-opt_data.score);
    //DEBUG_LAY("gradient: "<<opt_data.gradient.transpose())
    //DEBUG_LAY("hessian before SPD: \n"<<opt_data.hessian);
    pose_t delta_p;
    long total_score_l;
    opt_data.makeToSPD2();
   // DEBUG_LAY("hessian after SPD: \n"<<hessian);
    Eigen::SelfAdjointEigenSolver<hessian_t> saes(opt_data.hessian);
    double min_eigenval = saes.eigenvalues().minCoeff();
    if(min_eigenval<0){
      DEBUG_LAY("EIGEN values are negative - matrix not SPD");
      break;
    }
    //delta_p = hessian.inverse() *(-g); 
    delta_p = opt_data.hessian.ldlt().solve(-opt_data.gradient);
    delta = eigt::getTransFromPose(delta_p);
    total_score_l = -static_cast<long>(scoreLayer(delta * p,points));
    //long total_score_l = static_cast<long>(total_score);
    DEBUG_LAY("Score after Newton: "<<total_score_l);
    DEBUG_LAY("delta transform: "<<delta_p);
    
    
    //DEBUG_LAY("calc_trans after Newton: "<<eigt::getPoseFromTransform(p).transpose());
    //DEBUG_LAY("linsearch: "<<lineSearchMT(p,g,delta_p,points));
    double step = lineSearchMT(p,opt_data.gradient,delta_p,points);
    //delta_p = step * delta_p;
    DEBUG_LAY("Linear search step: "<<step);
    delta = eigt::getTransFromPose(delta_p);
    total_score_l = -static_cast<long>(scoreLayer(delta * p,points));
    DEBUG_LAY("Score after LS: "<<total_score_l)
    if(total_score_l < best_score && total_score_l > 0)
    {
        p = delta * p;
        best_p = p;
        best_score = total_score_l;
    }
    DEBUG_LAY("delta transform: "<<delta_p);
    DEBUG_LAY("calc_trans after Lin. step: "<<eigt::getPoseFromTransform(p).transpose());
    
   double change_in_trans = (eigt::getPoseFromTransform(old_p) - eigt::getPoseFromTransform(p)).norm();
   DEBUG_LAY("change in pose: "<<change_in_trans);
   if(iter >= MAX_ITER || change_in_trans < INC_CHANGE)
     converged = true;
    
    ++iter;
    DEBUG_LAY("\n");
  } // end of converge while

  transform_ = best_p;
  return true;
}

Layer::transform_t Layer::getTransformation(){
  return transform_;
}

OptimalisationData Layer::calcOptData(const transform_t &transf, const points_t  &points) const{
  OptimalisationData opt_data;


  Eigen::Matrix<double, 2, 3> jacobian;
  point_t hess_derivative;

  double angle = eigt::getAngle(transf);
  double si = std::sin(angle);
  double co = std::cos(angle);

  point_t difference;

  for (const auto & point : points) {
    point_t trans_point = transf*point;
    double x = trans_point(0);
    double y = trans_point(1);
    Field field;
    if(!getPointField(trans_point,field))
      continue;
    //DEBUG_LAY(field.getPoints());
    covar_t inv_covar =field.getInvCovar();
    //DEBUG_LAY(field.getInvCovar()<<"\n");
    difference = trans_point - field.getMean();
    double point_score = scorePoint(field,trans_point);
    //DEBUG_LAY(difference.dot(field.calcInvertedVariance() * difference)* -0.5F);
    //DEBUG_LAY(point_score);
    opt_data.score += point_score;
    // prepare for calculating H * delta_p = -g where delta_p is change in
    // transformation H is hessian and g is gradient
    // incrementaly summing parts of  g for each point
    jacobian << 1, 0, -x * si - y * co, 
                0, 1,  x * co - y * si;
    hess_derivative << -x * co + y * si,
                       -x * si - y * co;
    opt_data.gradient += pointGradient(difference,inv_covar,point_score,jacobian);
    opt_data.hessian += pointHessian(difference,inv_covar,point_score,jacobian,hess_derivative);
  }
  return opt_data;
}


std::string Layer::toString() const{
  std::string visual;
  for (size_t row = 0; row < size_; ++row) {
    for (size_t col = 0; col < size_; ++col) {
      if(row == size_/2 && col == size_/2)
        visual+="^ ";
      else
        visual+=fields_[row][col].toString();
    }
    visual+="\n";
  }
  return visual;

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
  //DEBUG_LAY("ML_NDT:Fields created for layer"<<size_<<"created "<<fields_.size());
  size_t id = 0;

  for (const auto & point : (*points_)) {
    if (isInBoundries(point)) {
      std::pair<size_t, size_t> coord = getFieldCoordintes(point);
      //DEBUG_LAY("ML_NDT: point:"<<point<< "coordinates"<<coord.first<<" "<<coord.second);
      fields_[coord.first][coord.second].addPoint(id);
    }
    ++id;
  }

  for(auto & row:  fields_ ){
    for(auto & field : row){
      field.calcNormDistribution();
    }
  }
  DEBUG_LAY("Layer initialized\n"+toString());
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
  DEBUG_LAY("LFD1 = "<<LFD1);
  DEBUG_LAY("LFD2 = "<<LFD2);
}

bool Layer::isInBoundries(const point_t & point) const {
  point_t pt = offset_inv_ * point;
 // DEBUG_LAY("point: "<<point<<" transform: "<<offset_.affine()<<" result: "<<pt);
  //point_t pt = point;
  if (std::abs(pt(0)) < max_range_ && std::abs(pt(1)) < max_range_)
    return true;
  else
    return false;
}



Layer::pose_t Layer::pointGradient(
    const point_t &difference, const covar_t &inv_covar, double score,
    const Eigen::Matrix<double, 2, 3> &jacobian) const
{
  pose_t g = pose_t::Zero();
  point_t diff_cvi = difference.transpose() * inv_covar;
  for (long r = 0; r < 3; ++r) {
    double diff_cvi_jacc = diff_cvi.dot(jacobian.col(r));
    g(r) = (diff_cvi_jacc * (-score));
  }
  return g;
}

Layer::hessian_t
Layer::pointHessian(const point_t &difference, const covar_t &inv_covar,
                    double score, const Eigen::Matrix<double, 2, 3> &jacobian,
                    const point_t &hessian_derivative) const
{
  hessian_t hessian;
  point_t diff_cvi = difference.transpose() * inv_covar;
  for (long i = 0; i < 3; ++i) {
    for (long j = 0; j < 3; ++j) {
      double jacc_i = diff_cvi.dot(jacobian.col(i));
      double jacc_j = (-1 * diff_cvi).dot(jacobian.col(j));
      point_t hess;
      if (i == 2 && j == 2)
        hess = hessian_derivative;
      else
        hess << 0, 0;
      hessian(i, j) = -score * (jacc_i * jacc_j + diff_cvi.dot(hess) +
                                (jacobian.col(j).transpose() * inv_covar)
                                    .dot(jacobian.col(i)));
    }
  }
  return hessian;
}
/*
  Brief: Calculates coordinates in layer's grid from point in parameter
*/
std::pair<size_t, size_t>
Layer::getFieldCoordintes(const point_t &point) const {
  point_t pt = offset_inv_ * point;
  // how many  meters of point locations are in one field
  double resolution = 2*max_range_ / size_;
  // move space of points from [-range,range] to [0,2range] in x axis
  size_t fieldx = 
      static_cast<size_t>(std::floor(((max_range_)- pt(0)) / resolution));
  // flips space of points in y axis from [0,range] to [range,0]
  size_t fieldy =
      static_cast<size_t>(std::floor(((max_range_) - pt(1)) / resolution));
  return std::make_pair(fieldx, fieldy);
}


void Layer::printLaserPoints(const points_t & points)const{
  DEBUG_LAY("points: \n");
  for(auto & pt : points){
    std::cout<<pt.transpose() << ";";
  }
  DEBUG_LAY("\n");
}

bool Layer::getPointField(const point_t & pt,Field & field)const {
  if(!isInBoundries(pt))
    return false;
  std::pair<size_t, size_t> coords = getFieldCoordintes(pt);
  if(!fields_[coords.first][coords.second].isReady())
    return false;
  field = fields_[coords.first][coords.second];
  return true;
}
// return score in form -score for purpose of minimalization probelm
double Layer::scorePoint(const Field & field, const point_t & trans_point)const{
  point_t difference = trans_point - field.getMean();
  return LFD1*std::exp((difference.dot(field.getInvCovar() * difference) * -0.5F*LFD2));
  //return -std::exp((difference.dot(field.getInvCovar() * difference) * -0.5F));
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
