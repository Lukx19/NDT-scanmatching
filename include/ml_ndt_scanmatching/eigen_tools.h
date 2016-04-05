#ifndef EIGEN_TOOLS
#define EIGEN_TOOLS

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592
#endif

namespace eigt
{
template <typename T>
using transform2d_t = Eigen::Transform<T, 2, Eigen::TransformTraits::Affine>;

template <typename T>
using point2d_t = Eigen::Matrix<T, 2, 1>;

template <typename T>
using pose2d_t = Eigen::Matrix<T, 3, 1>;
// ************************* DECLARATION *************************************

/// Calculates transformation between inputed poses. 
template <typename T>
transform2d_t<T> transBtwPoses(const pose2d_t<T> &from, const pose2d_t<T> &to);

// Transforms pose in global frame to its new posistion based on transform
template <typename T>
pose2d_t<T> transformPose(const pose2d_t<T> &pose,
                          const transform2d_t<T> &trans);

// maps from transformation matrix 3x3 to vector 
// [delta_x,delta_y, detla_angle]
template <typename T>
pose2d_t<T> getPoseFromTransform(const transform2d_t<T> &trans);

// returns angle which is rotational part of transformation
template <typename T>
T getAngle(const transform2d_t<T> &trans);

// returns square distance from translation part of transform
template <typename T>
T getDisplacement(const transform2d_t<T> &trans);

// returns absolute angular diffrence between poses. Selects shorter angle.
template <typename T>
T getAngleDiffrence(const pose2d_t<T> &from, const pose2d_t<T> &to);

// Maps 4x4 transformation matrix to 3x3 transformation matrix
template <typename T>
transform2d_t<T> convertToTransform(const Eigen::Matrix<T, 4, 4> &trans);

// MAPS from 3x3 transformation matrix to 4x4 Transformation matrix
template <typename T>
Eigen::Matrix<T, 4, 4> convertFromTransform(const transform2d_t<T> &trans);

// Maps transformation encoded in pose vector [delta_x,delta_y,delta_angle] to
// transformation matrix
template <typename T>
eigt::transform2d_t<T> getTransFromPose(const pose2d_t<T> &trans);

template <typename T>
T normalizeAngle(T angle);
}
// ************************* IMPLEMENTATION****************************
template <typename T>
eigt::transform2d_t<T> eigt::transBtwPoses(const pose2d_t<T> &from,
                                           const pose2d_t<T> &to)
{
  transform2d_t<T> t;
  Eigen::Rotation2D<T> rot_from(from(2));
  Eigen::Rotation2D<T> rot_to(to(2));
  t.setIdentity();
  t.matrix().block(0, 0, 2, 2) =
      (rot_from.toRotationMatrix().transpose() * rot_to.toRotationMatrix());
  t.matrix().block(0, 2, 2, 1) =
      rot_from.toRotationMatrix() * (to.head(2) - from.head(2));
  return t;
}

template <typename T>
eigt::pose2d_t<T> eigt::transformPose(const pose2d_t<T> &pose,
                                      const transform2d_t<T> &trans)
{
  Eigen::Rotation2D<T> rot(pose(2));
  pose2d_t<T> res_pose;
  pose2d_t<T> inc = getPoseFromTransform(trans);
  res_pose.head(2) = pose.head(2) + rot.toRotationMatrix().transpose() * inc.head(2);
  res_pose(2) = normalizeAngle(pose(2) + inc(2));
  return res_pose;
}

template <typename T>
eigt::pose2d_t<T> eigt::getPoseFromTransform(const transform2d_t<T> &trans)
{
  pose2d_t<T> pose;
  auto translation = trans.translation();
  pose(0) = translation(0);
  pose(1) = translation(1);
  pose(2) = eigt::getAngle(trans);
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
eigt::transform2d_t<T>
eigt::convertToTransform(const Eigen::Matrix<T, 4, 4> &trans)
{
  transform2d_t<T> new_trans;
  new_trans.setIdentity();
  new_trans.matrix().block(0, 0, 2, 2) = trans.block(0, 0, 2, 2);
  new_trans.matrix().block(0, 2, 2, 1) = trans.block(0, 3, 2, 1);
  return new_trans;
}

template <typename T>
Eigen::Matrix<T, 4, 4> eigt::convertFromTransform(const transform2d_t<T> &trans)
{
  Eigen::Matrix<T, 4, 4> new_trans;
  pose2d_t<T> pose = getPoseFromTransform<T>(trans);
  Eigen::AngleAxis<T> init_rotation(pose(2), Eigen::Matrix<T, 3, 1>::UnitZ());
  Eigen::Translation<T, 3> init_translation(pose(0), pose(1), 0);
  return (init_translation * init_rotation).matrix();
}

template <typename T>
eigt::transform2d_t<T> eigt::getTransFromPose(const pose2d_t<T> &trans)
{
  transform2d_t<T> trans_mtx;
  Eigen::Rotation2D<T> rot(trans(2));
  Eigen::Translation<T, 2> transl(trans.head(2));
  return transl * rot;
}

template <typename T>
T eigt::normalizeAngle(T angle)
{
  return atan2(sin(angle), cos(angle));
  //return angle - 2 * M_PI *std::floor((angle + M_PI) / (2 * M_PI));
}

#endif
