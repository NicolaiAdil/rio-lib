// src/rio_eskf.cpp
#include "rio/rio_eskf.h"

namespace rio {

RioEskf::RioEskf() = default;

void RioEskf::setParams(const Params& p) {
  params_ = p;
  params_.q_RB.normalize();
  params_set_ = true;
}

bool RioEskf::paramsSet() const { return params_set_; }

void RioEskf::reset(const NominalState& x0, const float* P0_diag_15, float t0) {
  if (!params_set_) return;

  x_ = x0;
  x_.q_WB.normalize();

  P_.setZero();
  if (P0_diag_15) {
    for (int i = 0; i < 15; ++i) P_(i, i) = P0_diag_15[i];
  }

  initialized_ = true;
  t_last_ = t0;
}

bool RioEskf::isInitialized() const { return initialized_; }
float RioEskf::lastTime() const { return t_last_; }

const NominalState& RioEskf::state() const { return x_; }
const Mat15& RioEskf::covariance() const { return P_; }

void RioEskf::propagate(const ImuSample& s) {
  if (!params_set_ || !initialized_) return;

  float dt = s.t - t_last_;
  if (dt < params_.min_dt) { t_last_ = s.t; return; }
  dt = clampf(dt, params_.min_dt, params_.max_dt);

#if defined(RIO_EIGEN_NO_MALLOC)
  Eigen::internal::set_is_malloc_allowed(false);
#endif

  // Unbias
  const Vec3 w = s.gyr - x_.b_g;
  const Vec3 a = s.acc - x_.b_a;

  // Nominal propagation
  {
    const Quat dq = quatExpSmall(w * dt);
    x_.q_WB = (x_.q_WB * dq).normalized();

    const Vec3 a_W = x_.q_WB * a + params_.g_W;

    x_.p_WB += x_.v_WB * dt + 0.5f * a_W * (dt * dt);
    x_.v_WB += a_W * dt;

    // Optional first-order bias decay
    if (params_.tau_ba > 0.0f) x_.b_a *= expf(-dt / params_.tau_ba);
    if (params_.tau_bg > 0.0f) x_.b_g *= expf(-dt / params_.tau_bg);
  }

  // Covariance propagation: P = F P F^T + Qd
  Mat15 F = Mat15::Identity();

  const Mat3 R_WB = x_.q_WB.toRotationMatrix();

  // dp += dv*dt
  F.block<3,3>(0, 3).diagonal().array() += dt;

  // dv += -R*[a]_x * dtheta * dt    (dtheta is 9..11)
  F.block<3,3>(3, 9) = -(R_WB * skew(a)) * dt;

  // dv += -R * dba * dt             (dba is 6..8)
  F.block<3,3>(3, 6) = -R_WB * dt;

  // dtheta += -dbg*dt               (dbg is 12..14)
  F.block<3,3>(9, 12).diagonal().array() += -dt;

  // Optional bias first-order terms (discrete approx)
  if (params_.tau_ba > 0.0f) {
    const float a_ba = -dt / params_.tau_ba;
    F.block<3,3>(6, 6).diagonal().array() += a_ba;
  }
  if (params_.tau_bg > 0.0f) {
    const float a_bg = -dt / params_.tau_bg;
    F.block<3,3>(12, 12).diagonal().array() += a_bg;
  }

  // Discrete noise (dt^2 approximation)
  const float sg2  = params_.sigma_gyr * params_.sigma_gyr;
  const float sa2  = params_.sigma_acc * params_.sigma_acc;
  const float sbg2 = params_.sigma_bg  * params_.sigma_bg;
  const float sba2 = params_.sigma_ba  * params_.sigma_ba;
  const float dt2 = dt * dt;

  Mat15 Qd = Mat15::Zero();
  Qd.block<3,3>(9, 9).diagonal().array()   += sg2  * dt2; // gyro noise -> dtheta
  Qd.block<3,3>(3, 3).diagonal().array()   += sa2  * dt2; // accel noise -> dv
  Qd.block<3,3>(6, 6).diagonal().array()   += sba2 * dt2; // accel-bias RW
  Qd.block<3,3>(12,12).diagonal().array()  += sbg2 * dt2; // gyro-bias RW

  Mat15 FP;
  FP.noalias() = F * P_;
  P_.noalias() = FP * F.transpose();
  P_ += Qd;

  t_last_ = s.t;

#if defined(RIO_EIGEN_NO_MALLOC)
  Eigen::internal::set_is_malloc_allowed(true);
#endif
}

void RioEskf::updateDoppler(const RadarDoppler* meas, size_t n) {
  if (!params_set_ || !initialized_ || !meas || n == 0) return;

#if defined(RIO_EIGEN_NO_MALLOC)
  Eigen::internal::set_is_malloc_allowed(false);
#endif

  const Mat3 C_BW = x_.q_WB.conjugate().toRotationMatrix();
  const Mat3 C_RB = params_.q_RB.toRotationMatrix();
  const Mat3 M = C_RB * C_BW; // dv_W -> dv_R

  const Vec3 v_B = C_BW * x_.v_WB;
  const Vec3 v_R = C_RB * v_B;

  // dv_R ≈ -C_RB*[v_B]_x * dtheta  (dtheta is 9..11)
  const Mat3 J_theta = -C_RB * skew(v_B);

  for (size_t i = 0; i < n; ++i) {
    Vec3 u = meas[i].u_R;
    const float un = u.norm();
    if (un < 1e-6f) continue;
    u /= un;

    const float sigma = (meas[i].sigma > 1e-6f) ? meas[i].sigma : 1e-2f;
    const float R = sigma * sigma;

    const float vr_hat = u.dot(v_R);
    const float r = meas[i].vr - vr_hat;

    Row15 H = Row15::Zero();
    H.block<1,3>(0, 3).noalias() = u.transpose() * M;       // dv at 3..5
    H.block<1,3>(0, 9).noalias() = u.transpose() * J_theta; // dtheta at 9..11

    scalarUpdate_(H, r, R);
  }

#if defined(RIO_EIGEN_NO_MALLOC)
  Eigen::internal::set_is_malloc_allowed(true);
#endif
}

void RioEskf::updateDopplerWithOmega(const RadarDoppler* meas, size_t n, const Vec3& omega_B) {
  if (!params_set_ || !initialized_ || !meas || n == 0) return;

#if defined(RIO_EIGEN_NO_MALLOC)
  Eigen::internal::set_is_malloc_allowed(false);
#endif

  const Mat3 C_BW = x_.q_WB.conjugate().toRotationMatrix();
  const Mat3 C_RB = params_.q_RB.toRotationMatrix();
  const Mat3 M = C_RB * C_BW;

  const Vec3 v_B = C_BW * x_.v_WB;
  const Vec3 v_radar_B = v_B + omega_B.cross(params_.p_BR_B);
  const Vec3 v_R = C_RB * v_radar_B;

  const Mat3 J_theta = -C_RB * skew(v_radar_B);

  for (size_t i = 0; i < n; ++i) {
    Vec3 u = meas[i].u_R;
    const float un = u.norm();
    if (un < 1e-6f) continue;
    u /= un;

    const float sigma = (meas[i].sigma > 1e-6f) ? meas[i].sigma : 1e-2f;
    const float R = sigma * sigma;

    const float vr_hat = u.dot(v_R);
    const float r = meas[i].vr - vr_hat;

    Row15 H = Row15::Zero();
    H.block<1,3>(0, 3).noalias() = u.transpose() * M;
    H.block<1,3>(0, 9).noalias() = u.transpose() * J_theta;

    scalarUpdate_(H, r, R);
  }

#if defined(RIO_EIGEN_NO_MALLOC)
  Eigen::internal::set_is_malloc_allowed(true);
#endif
}

void RioEskf::updateRadarVelocity3(const RadarVel3& z) {
  if (!params_set_ || !initialized_) return;

#if defined(RIO_EIGEN_NO_MALLOC)
  Eigen::internal::set_is_malloc_allowed(false);
#endif

  const Mat3 C_BW = x_.q_WB.conjugate().toRotationMatrix();
  const Mat3 C_RB = params_.q_RB.toRotationMatrix();
  const Mat3 M = C_RB * C_BW;

  const Vec3 v_B = C_BW * x_.v_WB;
  const Vec3 v_R = C_RB * v_B;

  const Vec3 r = z.v_R - v_R;
  const Mat3 J_theta = -C_RB * skew(v_B);

  for (int axis = 0; axis < 3; ++axis) {
    const float sigma = (z.sigma(axis) > 1e-6f) ? z.sigma(axis) : 1e-2f;
    const float R = sigma * sigma;

    Row15 H = Row15::Zero();
    H.block<1,3>(0, 3) = M.row(axis);       // dv
    H.block<1,3>(0, 9) = J_theta.row(axis); // dtheta

    scalarUpdate_(H, r(axis), R);
  }

#if defined(RIO_EIGEN_NO_MALLOC)
  Eigen::internal::set_is_malloc_allowed(true);
#endif
}

void RioEskf::scalarUpdate_(const Row15& H, float residual, float R) {
  // Scalar update (stable, symmetrized).
  const float S = (H * P_ * H.transpose())(0,0) + R;
  if (!(S > 1e-12f)) return;
  const float invS = 1.0f / S;

  const Vec15 PHt = P_ * H.transpose(); // 15x1
  const Vec15 K   = PHt * invS;         // 15x1

  inject_(K * residual);

  const Row15 HP = H * P_;

  Mat15 Pnew = P_;
  Pnew.noalias() -= K * HP;
  Pnew.noalias() -= PHt * K.transpose();
  Pnew.noalias() += (PHt * PHt.transpose()) * invS;

  P_ = 0.5f * (Pnew + Pnew.transpose());
}

void RioEskf::inject_(const Vec15& dx) {
  // dx order: [dp dv dba dtheta dbg]
  const Vec3 dp  = dx.segment<3>(0);
  const Vec3 dv  = dx.segment<3>(3);
  const Vec3 dba = dx.segment<3>(6);
  const Vec3 dth = dx.segment<3>(9);
  const Vec3 dbg = dx.segment<3>(12);

  x_.p_WB += dp;
  x_.v_WB += dv;
  x_.b_a  += dba;
  x_.b_g  += dbg;

  const Quat dq = quatExpSmall(dth);
  x_.q_WB = (x_.q_WB * dq).normalized();

  // Covariance reset for SO(3): affects dtheta block at 9..11
  const Mat3 Gth = Mat3::Identity() - 0.5f * skew(dth);

  Mat15 G = Mat15::Identity();
  G.block<3,3>(9, 9) = Gth;

  Mat15 GP;
  GP.noalias() = G * P_;
  P_.noalias() = GP * G.transpose();
  P_ = 0.5f * (P_ + P_.transpose());
}

} // namespace rio