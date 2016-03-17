#ifndef EIGEN_TOOLS
#define EIGEN_TOOLS

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace eigt
{
template <typename T>
using transform2d_t = Eigen::Transform<T, 2, Eigen::TransformTraits::Affine>;

template <typename T>
using point2d_t = Eigen::Matrix<T, 2, 1>;

template <typename T>
using pose2d_t = Eigen::Matrix<T, 3, 1>;
// ************************* DECLARATION *************************************
template <typename T>
transform2d_t<T> transBtwPoses(const pose2d_t<T> &from, const pose2d_t<T> &to);

template <typename T>
pose2d_t<T> transformPose(const pose2d_t<T> &pose,
                          const transform2d_t<T> &trans);

template <typename T>
pose2d_t<T> getPoseFromTransform(const transform2d_t<T> &trans);

template <typename T>
T getAngle(const transform2d_t<T> &trans);

template <typename T>
T getDisplacement(const transform2d_t<T> &trans);

template <typename T>
T getAngleDiffrence(const pose2d_t<T> &from, const pose2d_t<T> &to);

template <typename T>
transform2d_t<T>
convertTransform(const Eigen::Matrix<T, 4, 4> &trans);

template <typename T>
Eigen::Matrix<T, 4, 4> convertTransform(
    const transform2d_t<T> &trans);
}
// ************************* IMPLEMENTATION****************************
template <typename T>
eigt::transform2d_t<T> eigt::transBtwPoses(const pose2d_t<T> &from,
                                           const pose2d_t<T> &to)
{
  transform2d_t<T> t;
  t.setIdentity();
  T angle = getAngleDiffrence(from, to);
  Eigen::Rotation2D<T> rot(angle);
  t.matrix().block(0,0,2,2) = rot.toRotationMatrix();
  t.matrix().block(0,2,2,1) = to.head(2) - rot.toRotationMatrix() * from.head(2);
  return t;
}

template <typename T>
eigt::pose2d_t<T> eigt::transformPose(const pose2d_t<T> &pose,
                                      const transform2d_t<T> &trans)
{
  transform2d_t<T> t;
  t = trans *
      t.fromPositionOrientationScale(
          pose.head(2), Eigen::Rotation2D<T>(pose(2)), point2d_t<T>::Ones());
  return getPoseFromTransform<T>(t);
}

template <typename T>
eigt::pose2d_t<T> eigt::getPoseFromTransform(const transform2d_t<T> &trans)
{
  pose2d_t<T> pose;
  pose << trans(0, 2), trans(1, 2), std::atan2(trans(1, 0), trans(0, 0));
  return pose;
}

template <typename T>
T eigt::getAngle(const transform2d_t<T> &trans)
{
  return std::atan2(trans.rotation()(1, 0), trans.rotation()(0, 0));
}

template <typename T>
T eigt::getDisplacement(const transform2d_t<T> &trans)
{
  return std::sqrt(std::pow(trans.matrix()(0, 2), 2.0) +
                   std::pow(trans.matrix()(1, 2), 2.0));
}

template <typename T>
T eigt::getAngleDiffrence(const pose2d_t<T> &from, const pose2d_t<T> &to)
{
  return std::atan2(std::sin(to(2) - from(2)), std::cos(to(2) - from(2)));
}

template <typename T>
eigt::transform2d_t<T> eigt::convertTransform(const Eigen::Matrix<T, 4, 4> &trans)
{
  transform2d_t<T> new_trans;
  new_trans.matrix().block(0,0,2,2) = trans.block(0,0,2,2);
  new_trans.matrix().block(0,2,2,0) = trans.block(0,3,2,1);
  return new_trans;
}

template <typename T>
Eigen::Matrix<T, 4, 4> eigt::convertTransform(const transform2d_t<T> &trans)
{
  Eigen::Matrix<T, 4, 4> new_trans;
  pose2d_t<T> pose = getPoseFromTransform<T>(trans);
  Eigen::AngleAxis<T> init_rotation(pose(2), Eigen::Matrix<T,3,1>::UnitZ());
  Eigen::Translation<T, 3> init_translation(pose(0), pose(1), 0);
  return (init_translation * init_rotation).matrix();
}

#endif
