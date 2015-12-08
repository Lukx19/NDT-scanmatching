#ifndef NDT_LAYER
#define NDT_LAYER

#include <math.h>
#include <ros/ros.h>
#include "field.h"
#include "scanmatcher.h"

#ifndef DEBUG
  #define DEBUG(out); std::cout<<out<<"\n";
#endif
class Scanmatcher;
class Field;

class Layer {
private:
  typedef size_t Id_t;
  typedef Eigen::Vector3f pose_t;
  typedef Eigen::Vector2f point_t;
  typedef Eigen::Transform<float,2, Eigen::TransformTraits::Affine> transform_t;
  typedef std::vector<point_t> points_t;
  typedef Eigen::Matrix<Eigen::Vector2f, Eigen::Dynamic, Eigen::Dynamic>
      dyn_matrix_t;
  typedef std::vector<Field> field_line_t;
  typedef std::vector<field_line_t> field_grid_t;
public:
  Layer(points_t * points, size_t size, float max_range)
      : points_(points), size_(size),
        max_range_(max_range) {
    initializeFields(points);
    initializeParams();
  }
  // void addPoint(Id_t id,point_t pt);
  point_t getPoint(const Id_t id) const;
  //dyn_matrix_t getMeanVectors();
  //dyn_matrix_t getVarianceMatrices();

  bool calculateNdt(transform_t &transf, points_t &points);
  transform_t getTransformation();

private:
  // Scanmatcher & ndt_;
  points_t *points_;
  transform_t transform_;
  size_t size_;
  float max_range_;
  field_grid_t fields_;
  const size_t MAX_ITER = 10;
  const size_t MIN_POINTS_IN_FIELD = 3;
  float LFD1,LFD2;

  void initializeFields(points_t * points);
  void initializeParams();
  bool isInBoundries(const point_t &point) const;
  //point_t transformPoint(const Id_t id, const transform_t &transform) const;
  //point_t transformPoint(const point_t &point, const transform_t &transform) const;
  std::pair<size_t, size_t> getFieldCoordintes(const point_t &pt) const;
  Eigen::Matrix3f makeToSPD(const Eigen::Matrix3f &hess,const Eigen::Vector3f & s_gradient)const;
  Eigen::Matrix3f makeToSPD2(const Eigen::Matrix3f &hess,const Eigen::Vector3f & s_gradient)const;
  void printLaserPoints(const points_t & points)const;
  bool getPointField(const point_t & pt, Field & field)const;
  float scorePoint(const Field & field, const point_t & transformed_pt)const;
  float scoreLayer(const transform_t & trans, const points_t & cloud_in) const;
};

#endif
