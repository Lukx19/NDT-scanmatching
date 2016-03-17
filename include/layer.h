#ifndef NDT_LAYER
#define NDT_LAYER

#include <math.h>
#include <ros/ros.h>
#include <field.h>
#include <scanmatcher.h>
#include <ml_ndt_scanmatching/NDTMapMsg.h>

#ifndef DEBUG
  #define DEBUG(out); std::cout<<out<<"\n";
#endif
class Scanmatcher;
class Field;

class Layer {
private:
  typedef size_t Id_t;
  typedef Eigen::Vector3d pose_t;
  typedef Eigen::Vector2d point_t;
  typedef Eigen::Matrix2d covar_t;
  typedef Eigen::Matrix3d hessian_t;
  typedef Eigen::Transform<double,2, Eigen::TransformTraits::Affine> transform_t;
  typedef std::vector<point_t> points_t;
  typedef Eigen::Matrix<Eigen::Vector2d, Eigen::Dynamic, Eigen::Dynamic>
      dyn_matrix_t;
  typedef std::vector<Field> field_line_t;
  typedef std::vector<field_line_t> field_grid_t;
public:
  Layer(points_t * points, size_t size, double max_range, const transform_t & offset)
      : points_(points), size_(size),
        max_range_(max_range),
        offset_(offset),
        offset_inv_(offset.inverse()){
    initializeFields(points);
    initializeParams();
  }
  // void addPoint(Id_t id,point_t pt);
  point_t getPoint(const Id_t id) const;
  ml_ndt_scanmatching::NDTMapMsg getLayerData() const;

  bool calculateNdt(const transform_t &transf, const points_t  &points);
  transform_t getTransformation();
  std::string toString()const;

private:
  // Scanmatcher & ndt_;
  points_t *points_;
  transform_t transform_;
  size_t size_; // number of fields in layer. layer = size_ x size_ fields
  double max_range_; // max range in meters -> size_ = 2 *max_range meters
  transform_t offset_;
  transform_t offset_inv_;
  field_grid_t fields_;
  const size_t MAX_ITER = 50;
  const double INC_CHANGE = 0.001;
  const size_t MIN_POINTS_IN_FIELD = 4;
  double LFD1,LFD2;

  void initializeFields(points_t * points);
  void initializeParams();
  bool isInBoundries(const point_t &point) const;
  //point_t transformPoint(const Id_t id, const transform_t &transform) const;
  //point_t transformPoint(const point_t &point, const transform_t &transform) const;
  std::pair<size_t, size_t> getFieldCoordintes(const point_t &pt) const;
  Eigen::Matrix3d makeToSPD(const Eigen::Matrix3d &hess,const Eigen::Vector3d & s_gradient)const;
  Eigen::Matrix3d makeToSPD2(const Eigen::Matrix3d &hess,const Eigen::Vector3d & s_gradient)const;

  pose_t pointGradient(const point_t &difference, const covar_t &inv_covar,
                       double score,
                       const Eigen::Matrix<double, 2, 3> &jacobian) const;

  hessian_t pointHessian(const point_t &difference, const covar_t &inv_covar,
                         double score,
                         const Eigen::Matrix<double, 2, 3> &jacobian,
                         const point_t &hessian_derivative) const;

  void printLaserPoints(const points_t &points) const;
  bool getPointField(const point_t & pt, Field & field)const;
  double scorePoint(const Field & field, const point_t & transformed_pt)const;
  double scoreLayer(const transform_t & trans, const points_t & cloud_in) const;

  //perform line search to find the best descent rate (Mohre&Thuente)
  //adapted from NOX based on implementation:
  double lineSearchMT(const transform_t & trans,
                      const pose_t & gradient,
                      const pose_t & increment,
                      const points_t & cloud_in) const;
   struct MoreThuente
    {
        static double min(double a, double b);
        static double max(double a, double b);
        static double absmax(double a, double b, double c);
        static int cstep(double& stx, double& fx, double& dx,
                         double& sty, double& fy, double& dy,
                         double& stp, double& fp, double& dp,
                         bool& brackt, double stmin, double stmax);
    }; //end MoreThuente
};

#endif
