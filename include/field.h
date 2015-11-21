#ifndef NDT_FIELD
#define NDT_FIELD
#include <vector>
#include <Eigen/Dense>
#include "layer.h"
#define EIGEN_ALIGN_TO_BOUNDARY(n) __attribute__((aligned(n)))
class Layer;

class EIGEN_ALIGN_TO_BOUNDARY(8) Field{
  public:
    typedef std::size_t Id_t;
    typedef Eigen::Vector2f point_t;
    typedef Eigen::Matrix2f var_t;
    Field (Layer * layer)
      :layer_(layer),
      mean_(point_t::Zero()),
      variance_(var_t::Zero())
      {}
    //Field(Scanmatcher & ndt,std::size_t x0,std::size_t y0,std::size_t x1,std::size_t y1):
    //  ndt_(ndt),x0_(x0),y0_(y0),x1_(x1),y1_(y1){}
    void addPoint(Id_t id);
    point_t calcMean() const;
    var_t calcVariance() const;
  private:
    const Layer * layer_;
    point_t mean_;
    var_t variance_;

    //Scanmatcher & ndt_;
    //std::size_t x0_,y0_,x1_,y1_;

    std::vector<Id_t> points_ids_;

};


#endif
