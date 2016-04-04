#ifndef PCL_MLNDT2D
#define PCL_MLNDT2D

#include <ros/ros.h>
#include <pcl/registration/registration.h>
#include <pcl/registration/impl/ndt_2d.hpp>
#include <Eigen/Dense>

namespace pcl
{
template <typename PointSource, typename PointTarget>
class MultiLevelNdt2D : public Registration<PointSource, PointTarget>
{
  typedef typename Registration<PointSource, PointTarget>::PointCloudSource
      PointCloudSource;
  typedef typename PointCloudSource::Ptr PointCloudSourcePtr;
  typedef typename PointCloudSource::ConstPtr PointCloudSourceConstPtr;

  typedef typename Registration<PointSource, PointTarget>::PointCloudTarget
      PointCloudTarget;

  typedef PointIndices::Ptr PointIndicesPtr;
  typedef PointIndices::ConstPtr PointIndicesConstPtr;

public:
  typedef boost::shared_ptr<MultiLevelNdt2D<PointSource, PointTarget>> Ptr;
  typedef boost::shared_ptr<const MultiLevelNdt2D<PointSource, PointTarget>>
      ConstPtr;

  /** \brief Empty constructor. */
  MultiLevelNdt2D()
    : Registration<PointSource, PointTarget>()
    , grid_centre_(0, 0)
    , grid_extent_(20, 20)
    , newton_lambda_(1, 1, 1)
    , covariance_(Eigen::Matrix3d::Identity())
    , inform_matrix_(Eigen::Matrix3d::Identity())
  {
    reg_name_ = "NormalDistributionsTransform2D";
    Eigen::Matrix2d step;
    step << 0.5, 0.5;
    grid_step_perlevel_.push_back(step);
    grid_step_perlevel_.push_back(2 * step);
    grid_step_perlevel_.push_back(3 * step);
    grid_step_perlevel_.push_back(4 * step);
  }

  /** \brief Empty destructor */
  virtual ~MultiLevelNdt2D()
  {
  }

  /** \brief centre of the ndt grid (target coordinate system)
    * \param centre value to set
    */
  virtual void setGridCentre(const Eigen::Vector2d& centre)
  {
    grid_centre_ = centre;
  }

  /** \brief Grid spacing (step) of the NDT grid
    * \param[in] steps vector of values to set for each level
    */
  virtual void setGridStep(const std::vector<Eigen::Vector2d>& steps)
  {
    grid_step_perlevel_ = steps;
  }

  /** \brief NDT Grid extent (in either direction from the grid centre)
    * \param[in] extent value to set
    */
  virtual void setGridExtent(const Eigen::Vector2d& extent)
  {
    grid_extent_ = extent;
  }

  /** \brief NDT Newton optimisation step size parameter
    * \param[in] lambda step size: 1 is simple newton optimisation, smaller
   * values may improve convergence
    */
  virtual void setOptimizationStepSize(const double& lambda)
  {
    newton_lambda_ = Eigen::Vector3d(lambda, lambda, lambda);
  }

  /** \brief NDT Newton optimisation step size parameter
    * \param[in] lambda step size: (1,1,1) is simple newton optimisation,
    * smaller values may improve convergence, or elements may be set to
    * zero to prevent optimisation over some parameters
    *
    * This overload allows control of updates to the individual (x, y,
    * theta) free parameters in the optimisation. If, for example, theta is
    * believed to be close to the correct value a small value of lambda[2]
    * should be used.
    */
  virtual void setOptimizationStepSize(const Eigen::Vector3d& lambda)
  {
    newton_lambda_ = lambda;
  }

  virtual Eigen::Matrix3f getCovariance() const
  {
    return covariance_;
  }
  virtual Eigen::Matrix3f getInformMatrix() const
  {
    return inform_matrix_;
  }

protected:
  using Registration<PointSource, PointTarget>::reg_name_;
  using Registration<PointSource, PointTarget>::target_;
  using Registration<PointSource, PointTarget>::converged_;
  using Registration<PointSource, PointTarget>::nr_iterations_;
  using Registration<PointSource, PointTarget>::max_iterations_;
  using Registration<PointSource, PointTarget>::transformation_epsilon_;
  using Registration<PointSource, PointTarget>::transformation_;
  using Registration<PointSource, PointTarget>::previous_transformation_;
  using Registration<PointSource, PointTarget>::final_transformation_;
  using Registration<PointSource, PointTarget>::update_visualizer_;
  using Registration<PointSource, PointTarget>::indices_;

  Eigen::Vector2d grid_centre_;
  Eigen::Vector2d grid_extent_;
  Eigen::Vector3d newton_lambda_;
  Eigen::Matrix3d covariance_;
  Eigen::Matrix3d inform_matrix_;
  std::vector<Eigen::Vector2d> grid_step_perlevel_;

  /** \brief Rigid transformation computation method with initial guess.
    * \param[out] output the transformed input point cloud dataset using the
   * rigid transformation found
    * \param[in] guess the initial guess of the transformation to compute
    */
  virtual void computeTransformation(PointCloudSource& output,
                                     const Eigen::Matrix4f& guess);
  virtual void
  computeSingleGrid(const PointCloudSource & source_cloud,
                    ndt2d::ValueAndDerivatives<3, double>& result_data,
                    Eigen::Matrix4f& trans,
                    const ndt2d::NDTSingleGrid<PointTarget>& target_ndt);

  virtual Eigen::Matrix3d makeToSPD(const Eigen::Matrix3d& hessian,
                                    const Eigen::Vector3d& gradient) const;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


template <typename PointSource, typename PointTarget>
void MultiLevelNdt2D<PointSource,PointTarget>
      ::computeTransformation(PointCloudSource& output,const Eigen::Matrix4f& guess)
{
  std::vector<ndt2d::NDTSingleGrid<PointTarget>> grids;
  for(auto & level:grid_step_perlevel_){
    grids.emplace_back(ndt2d::NDTSingleGrid<PointTarget>(target_,grid_centre_,grid_extent_,level));
  }
  // calculate transformation
  ndt2d::ValueAndDerivatives<3, double> res;
  Eigen::Matrix4f trans = guess;
  for(auto it = grids.rbegin(); it != grids.rend();++it){
    computeSingleGrid(output,res,trans,*it);
  }
  // save covariance of meassurement
  Eigen::FullPivLU<Eigen::Matrix<double, 3, 3>> dec(res.hessian);
  Eigen::Matrix<double, 3, 3> invH;
  if (dec.isInvertible()) {
    ROS_DEBUG_STREAM("nice, we have invertible Hessian\n" << res.hessian
                                                          << "\n");
    covariance_ = dec.inverse();
  } else {
    ROS_ERROR_STREAM("Hessian is not invertible:\n" << res.hessian << "\n");
    covariance_.setIdentity();
  }
  inform_matrix_ = res.hessian;
  final_transformation_ = trans;
  // trasform original source cloud to new position based on calculated transform
  PointCloudSource intm_cloud
  pcl::transformPointCloud(output,intm_cloud,trans);
  output = intm_cloud;
}

template <typename PointSource, typename PointTarget>
void MultiLevelNdt2D<PointSource, PointTarget>::computeSingleGrid(
    const PointCloudSource& output,
    ndt2d::ValueAndDerivatives<3, double>& result_data, 
    Eigen::Matrix4f& trans,
    const ndt2d::NDTSingleGrid<PointTarget>& target_ndt)
{
  PointCloudSource intm_cloud = output;

  if (trans != Eigen::Matrix4f::Identity()) {
    transformation_ = trans;
    transformPointCloud(output, intm_cloud, transformation_);
  }

  Eigen::Matrix4f& transformation = transformation_;

  // transform from transformation matrix to transformation vector [dx,dy,theta]
  const Eigen::Matrix3f initial_rot(transformation.block<3, 3>(0, 0));
  const Eigen::Vector3f rot_x(initial_rot * Eigen::Vector3f::UnitX());
  const double z_rotation = std::atan2(rot_x[1], rot_x[0]);

  Eigen::Vector3d xytheta_transformation(transformation(0, 3),
                                         transformation(1, 3), z_rotation);

  ndt2d::ValueAndDerivatives<3, double> score;
  while (!converged_) {
    // prepare data for hessian and jacobian calculation
    const double cos_theta = std::cos(xytheta_transformation[2]);
    const double sin_theta = std::sin(xytheta_transformation[2]);
    previous_transformation_ = transformation;

    // calculate hessian, score and jaccobian
    score = ndt2d::ValueAndDerivatives<3, double>::Zero();
    for (size_t i = 0; i < intm_cloud.size(); i++)
      score += target_ndt.test(intm_cloud[i], cos_theta, sin_theta);

    ROS_DEBUG("[pcl::NormalDistributionsTransform2D::computeTransformation] "
              "NDT score %f (x=%f,y=%f,r=%f)\n",
              float(score.value), xytheta_transformation[0],
              xytheta_transformation[1], xytheta_transformation[2]);

    if (score.value <= 0.000001 && score.value > -0.000001) {
      score.hessian = makeToSPD(score.hessian, score.grad);
      Eigen::Vector3d delta_transformation;
      delta_transformation = score.hessian.ldlt().solve(-score.grad);
      // TODO replace by more sophisticated linear search method
      Eigen::Vector3d new_transformation =
          xytheta_transformation +
          newton_lambda_.cwiseProduct(delta_transformation);

      xytheta_transformation = new_transformation;

      // update transformation matrix from x, y, theta:
      transformation.block<3, 3>(0, 0).matrix() = Eigen::Matrix3f(
          Eigen::AngleAxisf(static_cast<float>(xytheta_transformation[2]),
                            Eigen::Vector3f::UnitZ()));
      transformation.block<3, 1>(0, 3).matrix() =
          Eigen::Vector3f(static_cast<float>(xytheta_transformation[0]),
                          static_cast<float>(xytheta_transformation[1]), 0.0f);
    } else {
      ROS_ERROR("[pcl::NormalDistributionsTransform2D::computeTransformation] "
                "no overlap of scans");
    }

    transformPointCloud(output, intm_cloud, transformation);
    nr_iterations_++;
    if (update_visualizer_ != 0)
      update_visualizer_(output, *indices_, *target_, *indices_);

    if (nr_iterations_ > max_iterations_ ||
        (transformation - previous_transformation_).norm() <
            transformation_epsilon_)
      converged_ = true;

  }  // end of convergence loop
  //output = intm_cloud;
  result_data = score;
}

template <typename PointSource, typename PointTarget>
Eigen::Matrix3d MultiLevelNdt2D<PointSource, PointTarget>::makeToSPD(
    const Eigen::Matrix3d& hessian, const Eigen::Vector3d& gradient) const
{
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> saes(hessian);
  Eigen::Vector3d eigenvalues = saes.eigenvalues();
  double min = eigenvalues.minCoeff();
  double max = eigenvalues.maxCoeff();
  if (min < 0) {
    Eigen::Matrix3d eigenvectors = saes.eigenvectors();
    double additive = gradient.norm();
    if (min + additive <= 0) {
      additive = 0.001 * max - min;
    }
    Eigen::Vector3d addition;
    addition << additive, additive, additive;
    eigenvalues += addition;
    return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
  } else {
    return hessian;
  }
}

}  // end of pcl namespace
#endif
