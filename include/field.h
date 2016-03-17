
#ifndef NDT_FIELD
#define NDT_FIELD
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <layer.h>
#include <ros/ros.h>
#include <ml_ndt_scanmatching/NDTCellMsg.h>
#define EIGEN_ALIGN_TO_BOUNDARY(n) __attribute__((aligned(n)))

#ifndef DEBUG
  #define DEBUG(out); std::cout<<out<<"\n";
#endif

//EIGEN_ALIGN_TO_BOUNDARY(8) 
class Layer;

class Field{
  public:
    typedef std::size_t Id_t;
    typedef Eigen::Vector2d point_t;
    typedef std::vector<point_t> points2_t;
    typedef Eigen::Matrix2d var_t;
    Field()
      :mean_(point_t::Zero()),
      covar_(var_t::Zero()),
      inv_covar_(var_t::Zero()),
      is_calc_(false)
      {}

    Field (points2_t * points)
      :points_(points),
      mean_(point_t::Zero()),
      covar_(var_t::Zero()),
      inv_covar_(var_t::Zero()),
      is_calc_(false)
      {}
    //Field(Scanmatcher & ndt,std::size_t x0,std::size_t y0,std::size_t x1,std::size_t y1):
    //  ndt_(ndt),x0_(x0),y0_(y0),x1_(x1),y1_(y1){}
    void addPoint(Id_t id);
    void calcNormDistribution();
    ml_ndt_scanmatching::NDTCellMsg getCellData()const;
    point_t getMean() const;
    var_t getCovar() const;
    var_t getInvCovar() const;
    size_t getPoints() const;
    bool isReady() const;
    std::string toString() const;


  private:
    points2_t * points_;
    point_t mean_;
    var_t covar_;
    var_t inv_covar_;
    bool is_calc_;
    std::vector<Id_t> points_ids_;

    double EVAL_FACTOR = 100;
    size_t MIN_PTR_EVAL = 4;

    void prepNormDist();
    point_t calcMean() const;
    var_t calcCovariance(const point_t & mean) const;
    std::tuple<var_t,var_t,bool>
      calcInvertedCovariance(const var_t & covar)const;

    template<typename T>
    T pinv(const T & mat, double tolerance = 1.e-06f)const;

};

template<typename T>
T Field::pinv(const T & mat, double tolerance)const
{
  Eigen::MatrixXd matX(mat.rows(),mat.cols());
  for(long row = 0; row < mat.rows();++row)
    for(long col = 0; col< mat.cols();++col)
      matX(row,col) = mat(row,col);


  Eigen::JacobiSVD<Eigen::MatrixXd> svd_mat(matX, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd u = svd_mat.matrixU();
  Eigen::MatrixXd v = svd_mat.matrixV();
  Eigen::VectorXd s = svd_mat.singularValues();
  double max = 0;
  for(long i =0; i < s.rows();++i){
    if(std::abs(s(i)) > max)
      max = std::abs(s(i));
  }

  for(long i =0; i < s.rows();++i){
    if(std::abs(s(i)) > max * tolerance){
      s(i) = 1/(s(i));
    }else{
      s(i) = 0;
    }
  }
  Eigen::MatrixXd res = v * s.asDiagonal() * u.transpose();
  T res_correct_type;
  for(long row = 0; row < res.rows();++row)
    for(long col = 0; col< res.cols();++col)
      res_correct_type(row,col) = res(row,col);
  return res_correct_type;
}

#endif
