#pragma once

#include "rio_config.h"

#include <math.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace rio {

using Vec3  = Eigen::Matrix<float, 3, 1>;
using Mat3  = Eigen::Matrix<float, 3, 3>;
using Quat  = Eigen::Quaternionf;

using Mat12 = Eigen::Matrix<float, 12, 12>;

using Vec21 = Eigen::Matrix<float, 21, 1>;
using Mat21 = Eigen::Matrix<float, 21, 21>;
using Row21 = Eigen::Matrix<float, 1, 21>;
using Mat21x12 = Eigen::Matrix<float, 21, 12>;

using Vec23 = Eigen::Matrix<float, 23, 1>;
using Mat23 = Eigen::Matrix<float, 23, 23>;
using Row23 = Eigen::Matrix<float, 1, 23>;

inline float clampf(float x, float lo, float hi) {
  return (x < lo) ? lo : (x > hi) ? hi : x;
}

inline Mat3 skew(const Vec3& v) {
  Mat3 S;
  S <<   0.0f, -v.z(),  v.y(),
        v.z(),  0.0f, -v.x(),
       -v.y(),  v.x(),  0.0f;
  return S;
}

inline Quat quatExpSmall(const Vec3& dtheta) {
  const float a2 = dtheta.squaredNorm();
  if (a2 < 1e-12f) {
    return Quat(1.0f, 0.5f*dtheta.x(), 0.5f*dtheta.y(), 0.5f*dtheta.z());
  }
  const float a = sqrtf(a2);
  const float half_a = 0.5f * a;
  const float s = sinf(half_a) / a;
  return Quat(cosf(half_a), dtheta.x()*s, dtheta.y()*s, dtheta.z()*s);
}

} // namespace rio