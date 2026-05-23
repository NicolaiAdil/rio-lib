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

// Hypsometric formula for differential altitude given two pressure samples
// taken at (approximately) the same temperature. Uses local temperature T,
// so the result tracks the actual atmosphere better than ISA-with-T0=288.15K
// would over a short interval.
//
//   Δh = (R · T / g) · ln(p_prev / p_curr)
//
// where T is in Kelvin. p_prev and p_curr in Pa (units cancel).
inline float differentialAltitude(float p_prev_pa, float p_curr_pa,
                                  float T_kelvin) {
  constexpr float R_spec = 287.05f;
  constexpr float g_std  = 9.80665f;
  return (R_spec * T_kelvin / g_std) * logf(p_prev_pa / p_curr_pa);
}

// NASA Earth Atmosphere Model (Troposphere; valid up to ~11 km) —
// absolute pressure (Pa) → altitude (m) using fixed ISA constants
// rather than the locally measured temperature. Eq. 12 of Girod et
// al., "A robust baro-radar-inertial odometry m-estimator",
// arXiv:2408.05764. Used by the BRIO barometric factor: anchor on
// the first reading, then the residual measures state z against
// (z_baro_i − z_baro_anchor), independent of the sensor's reported
// temperature.
inline float pressureToAltitude(float p_pa) {
  return (288.08f * powf(p_pa / 101290.0f, 1.0f / 5.256f)
          - 273.1f - 15.04f) / -0.00649f;
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