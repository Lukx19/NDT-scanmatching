#ifndef PCL_NDT2D_EXTENDED
#define PCL_NDT2D_EXTENDED

#include <ros/ros.h>
#include <pcl/registration/ndt_2d.h>
#include <Eigen/Dense>
namespace pcl
{
template <typename PointSource, typename PointTarget>
class NormalDistributionsTransform2DEx
    : public NormalDistributionsTransform2D<PointSource, PointTarget>
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
  NormalDistributionsTransform2DEx()
    : NormalDistributionsTransform2D<PointSource, PointTarget>()
  {
  }
  virtual ~NormalDistributionsTransform2DEx()
  {
  }
  virtual Eigen::Matrix3d getCovariance() const;
  virtual Eigen::Matrix3d getInformMatrix() const;

protected:
  using NormalDistributionsTransform2D<PointSource,PointTarget>::converged_;
  using NormalDistributionsTransform2D<PointSource,PointTarget>::nr_iterations_;
  using NormalDistributionsTransform2D<PointSource,PointTarget>::max_iterations_;
  using NormalDistributionsTransform2D<PointSource,PointTarget>::transformation_epsilon_;
  using NormalDistributionsTransform2D<PointSource,PointTarget>::previous_transformation_;
  using NormalDistributionsTransform2D<PointSource,PointTarget>::final_transformation_;
  using NormalDistributionsTransform2D<PointSource,PointTarget>::transformation_;
  using NormalDistributionsTransform2D<PointSource,PointTarget>::update_visualizer_;
  using NormalDistributionsTransform2D<PointSource,PointTarget>::indices_;
  using NormalDistributionsTransform2D<PointSource,PointTarget>::target_;
  using NormalDistributionsTransform2D<PointSource,PointTarget>::grid_centre_;
  using NormalDistributionsTransform2D<PointSource,PointTarget>::grid_extent_;
  using NormalDistributionsTransform2D<PointSource,PointTarget>::grid_step_;
  using NormalDistributionsTransform2D<PointSource,PointTarget>::newton_lambda_;
  Eigen::Matrix3d covariance_;
  Eigen::Matrix3d inform_matrix_;

  virtual void computeTransformation(PointCloudSource &output,
                                     const Eigen::Matrix4f &guess);
  virtual Eigen::Matrix3d makeToSPD(const Eigen::Matrix3d & hessian, const Eigen::Vector3d & gradient) const;
};

template <typename PointSource, typename PointTarget>
Eigen::Matrix3d
NormalDistributionsTransform2DEx<PointSource, PointTarget>::getCovariance() const
{
  return covariance_;
}

template <typename PointSource, typename PointTarget>
Eigen::Matrix3d
NormalDistributionsTransform2DEx<PointSource, PointTarget>::getInformMatrix() const
{
  return inform_matrix_;
}

template <typename PointSource, typename PointTarget>
void NormalDistributionsTransform2DEx<PointSource,PointTarget>
      ::computeTransformation(PointCloudSource &output,const Eigen::Matrix4f &guess)
{
  PointCloudSource intm_cloud = output;

  if (guess != Eigen::Matrix4f::Identity()) {
    transformation_ = guess;
    transformPointCloud(output, intm_cloud, transformation_);
  }

  // build Normal Distribution Transform of target cloud:
  ndt2d::NDT2D<PointTarget> target_ndt(target_, grid_centre_, grid_extent_,
                                       grid_step_);
  Eigen::Matrix4f &transformation = transformation_;

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
  
  // save covariance of meassurement
  Eigen::FullPivLU<Eigen::Matrix<double, 3, 3>> dec(score.hessian);
  Eigen::Matrix<double, 3, 3> invH;
  if (dec.isInvertible()) {
    ROS_DEBUG_STREAM("nice, we have invertible Hessian\n" << score.hessian
                                                          << "\n");
    covariance_ = dec.inverse();
  } else {
    ROS_ERROR_STREAM("Hessian is not invertible:\n" << score.hessian << "\n");
    covariance_.setIdentity();
  }
  inform_matrix_ = score.hessian;
  final_transformation_ = transformation;
  output = intm_cloud;
}


template <typename PointSource, typename PointTarget>
Eigen::Matrix3d NormalDistributionsTransform2DEx<PointSource,PointTarget>
                  ::makeToSPD(const Eigen::Matrix3d & hessian, const Eigen::Vector3d & gradient) const
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

}  // end of namespace pcl

#endif
