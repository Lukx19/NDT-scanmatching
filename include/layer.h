#ifndef NDT_LAYER
#define NDT_LAYER

#include <math.h>
#include "field.h"
#include "scanmatcher.h"

class Scanmatcher;
class Field;

class Layer {
public:
  typedef size_t Id_t;
  typedef Eigen::Vector3f pose_t;
  typedef Eigen::Vector2f point_t;
  typedef std::vector<point_t> points_t;
  typedef Eigen::Matrix<Eigen::Vector2f, Eigen::Dynamic, Eigen::Dynamic>
      dyn_matrix_t;
  typedef std::vector<Field> field_line_t;
  typedef std::vector<field_line_t> field_grid_t;
  Layer(points_t *points, size_t size, float max_range)
      : points_(points), transform_(pose_t::Zero()), size_(size),
        max_range_(max_range) {
    initializeFields();
  }
  // void addPoint(Id_t id,point_t pt);
  point_t getPoint(const Id_t id) const;
  dyn_matrix_t getMeanVectors();
  dyn_matrix_t getVarianceMatrices();

  bool calculateNdt(pose_t &transf, points_t &points);
  pose_t getTransformation();

private:
  // Scanmatcher & ndt_;
  points_t *points_;
  pose_t transform_;
  size_t size_;
  float max_range_;
  field_grid_t fields_;

  void initializeFields();
  bool isInBoundries(const point_t &point) const;
  point_t transformPoint(const Id_t id, const pose_t &transform) const;
  point_t transformPoint(const point_t &point, const pose_t &transform) const;
  std::pair<size_t, size_t> getFieldCoordintes(const point_t &pt) const;
  Eigen::Matrix3f makeToSPD(Eigen::Matrix3f &mat);
};

#endif
