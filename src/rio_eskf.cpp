// src/rio_eskf.cpp 
#include "rio/rio_eskf.h"

namespace rio {

RioEskf::RioEskf() = default;

void RioEskf::setParams(const Params& p) {
  params_ = p;
  params_.q_IR.normalize();

  // Build continuous-time process noise covariance Q (12x12 diagonal).
  // Column ordering matches generateE: [acc, ba, gyr, bg]
  const float sa2  = params_.sigma_acc * params_.sigma_acc;
  const float sba2 = params_.sigma_ba  * params_.sigma_ba;
  const float sg2  = params_.sigma_gyr * params_.sigma_gyr;
  const float sbg2 = params_.sigma_bg  * params_.sigma_bg;

  Q_ = Mat12::Zero();
  Q_.block<3,3>(0, 0).diagonal().setConstant(sa2);
  Q_.block<3,3>(3, 3).diagonal().setConstant(sba2);
  Q_.block<3,3>(6, 6).diagonal().setConstant(sg2);
  Q_.block<3,3>(9, 9).diagonal().setConstant(sbg2);

  params_set_ = true;
}

bool RioEskf::paramsSet() const { return params_set_; }

void RioEskf::reset(const NominalState& x0, const float* P0_diag_21, float t0) {
  if (!params_set_) return;

  x_ = x0;
  x_.q_WI.normalize();

  P_hat_.setZero();
  if (P0_diag_21) {
    for (int i = 0; i < 21; ++i) P_hat_(i, i) = P0_diag_21[i];
  }
  P_hat_prior_ = P_hat_;

  delta_x_hat_.setZero();
  delta_x_hat_prior_.setZero();

  initialized_ = true;
  t_last_ = t0;
}

bool RioEskf::initAttitudeFromGravity(const Vec3& f_b, const float* P0_diag,
                                       float t0, float g_tol) {
  if (!params_set_) return false;

  const float g_mag = params_.g_W.norm();        // expected ~9.81
  const float fn    = f_b.norm();

  // Reject if |f_b| is too far from expected gravity magnitude
  if (fn < (g_mag - g_tol) || fn > (g_mag + g_tol)) return false;

  // Normalised gravity direction in body frame
  const Vec3 gb = -f_b / std::max(1e-6f, fn);

  // Roll and pitch from gravity (yaw unobservable without compass)
  const float roll  = std::atan2(gb.y(), gb.z());
  const float pitch = std::atan2(-gb.x(),
                        std::sqrt(gb.y() * gb.y() + gb.z() * gb.z()));
  const float yaw   = 0.0f;

  // Build quaternion from Euler angles (ZYX / RPY convention)
  const float cr = std::cos(roll  * 0.5f), sr = std::sin(roll  * 0.5f);
  const float cp = std::cos(pitch * 0.5f), sp = std::sin(pitch * 0.5f);
  const float cy = std::cos(yaw   * 0.5f), sy = std::sin(yaw   * 0.5f);

  Quat q_WI;
  q_WI.w() = cr * cp * cy + sr * sp * sy;
  q_WI.x() = sr * cp * cy - cr * sp * sy;
  q_WI.y() = cr * sp * cy + sr * cp * sy;
  q_WI.z() = cr * cp * sy - sr * sp * cy;
  q_WI.normalize();

  // Preserve current extrinsics in the reset state
  NominalState x0 = x_;
  x0.q_WI = q_WI;

  reset(x0, P0_diag, t0);
  return true;
}

bool RioEskf::isInitialized() const { return initialized_; }
float RioEskf::lastTime() const { return t_last_; }

const NominalState& RioEskf::getState() const { return x_; }
const Mat21& RioEskf::getCovariance() const { return P_hat_; }

void RioEskf::predict(const ImuSample& s, float dt) {
  if (!params_set_ || !initialized_) return;

  if (dt < params_.min_dt || dt > params_.max_dt) { t_last_ = s.t; return; }
  // dt = clampf(dt, params_.min_dt, params_.max_dt);

#if defined(RIO_EIGEN_NO_MALLOC)
  Eigen::internal::set_is_malloc_allowed(false);
#endif

  // Unbias
  const Vec3 w_nom = s.gyr - x_.b_g;
  const Vec3 f_nom = s.acc - x_.b_a;

  const Mat21 A = generateA(f_nom, w_nom);
  const Mat21 Ad = Mat21::Identity() + A * dt;

  const Mat21x12 E = generateE();
  // Q is the PSD matrix (Eq. 129 and 130 in Trawny05b), make rate-independent
  Mat12 Qc = Q_;
  Qc.block<3,3>(0, 0) /= dt;   // acc noise:  σ_a² / dt
  Qc.block<3,3>(3, 3) *= dt;   // acc bias RW: σ_ba² * dt
  Qc.block<3,3>(6, 6) /= dt;   // gyro noise: σ_g² / dt
  Qc.block<3,3>(9, 9) *= dt;   // gyro bias RW: σ_bg² * dt

  const Mat21 Qd = E * Qc * E.transpose() * dt;

  // Predictor update (Fossen 2nd, eq. 14.206, 14.207)
  delta_x_hat_prior_ = Ad * delta_x_hat_; // This is always zero, could comment it out.
  P_hat_prior_ = Ad * P_hat_ * Ad.transpose() + Qd;

  t_last_ = s.t;

#if defined(RIO_EIGEN_NO_MALLOC)
  Eigen::internal::set_is_malloc_allowed(true);
#endif
}

void RioEskf::insPropagation(const ImuSample& s, float dt) {
  const auto f_nom = s.acc - x_.b_a;
  const auto w_nom = s.gyr - x_.b_g;
  // Rotation body -> NED from current quaternion
  const Mat3 R_WI = x_.q_WI.toRotationMatrix();

  // Linear acceleration in NED
  const Vec3 a_W = R_WI * f_nom + params_.g_W;

  // Position and velocity update
  x_.p_WI += x_.v_WI * dt + 0.5f * a_W * dt * dt;
  x_.v_WI += a_W * dt;

  // Attitude update: R_next = R_WI * exp([w_nom]_x * dt)
  const Quat dq = quatExpSmall(w_nom * dt);
  x_.q_WI = (x_.q_WI * dq).normalized();

  // Biases and extrinsics are unchanged (no propagation)
}

CorrectionResult RioEskf::correct(const RadarDoppler* meas, size_t n, const ImuSample& s) {
  CorrectionResult res;
  res.n_total = n;

  if (!params_set_ || !initialized_ || !meas || n == 0) return res;

  const float R_meas = params_.sigma_vr * params_.sigma_vr;
  const float gate_thresh = params_.gate_nsigma * params_.gate_nsigma;

  bool any_update = false;

  for (size_t i = 0; i < n; ++i) {
    Vec3 mu_r = meas[i].u_R;
    const float un = mu_r.norm();
    if (un < 1e-6f) { res.n_skipped++; continue; }
    mu_r /= un;

    // Compute H (1x21) and h (predicted vr)
    auto w_nom = s.gyr - x_.b_g;
    const Row21 H = computeRadarH_(mu_r, w_nom);
    const float h = computeRadarh_(mu_r, w_nom);

    // Residual: e = z - h(x)
    const float e = params_.vr_sign * meas[i].vr - h;

    // Gating
    if (params_.gating_enable) {
      const float S = (H * P_hat_prior_ * H.transpose())(0, 0) + R_meas;
      if (S <= 0.0f) { res.n_skipped++; continue; }
      const float gamma = e * e / S;
      if (gamma > gate_thresh) { res.n_rejected++; continue; }
    }

    // Scalar Kalman correction (Joseph form)
    scalarCorrect_(H, e, R_meas);

    // Inject into nominal state
    updateStateEstimate(delta_x_hat_);

    res.n_accepted++;
    any_update = true;

    // Carry posterior into prior for next measurement in this batch
    P_hat_prior_       = P_hat_;
    delta_x_hat_prior_ = delta_x_hat_;
  }

  if (!any_update) {
    P_hat_       = P_hat_prior_;
    delta_x_hat_ = delta_x_hat_prior_;
  }

  return res;
}

void RioEskf::scalarCorrect_(const Row21& H, float residual, float R) {
  // S = H P H^T + R
  const float S = (H * P_hat_prior_ * H.transpose())(0, 0) + R;
  if (!(S > 1e-12f)) return;
  const float invS = 1.0f / S;

  // K = P_prior * H^T * S^{-1}
  const Vec21 PHt = P_hat_prior_ * H.transpose();
  const Vec21 K   = PHt * invS;

  // Error-state update: δx = K * e
  delta_x_hat_ = K * residual;

  // Covariance update (Joseph form):
  // IKH = I - K H
  // P = IKH * P_prior * IKH^T + K * R * K^T
  const Mat21 IKH = Mat21::Identity() - K * H;
  P_hat_ = IKH * P_hat_prior_ * IKH.transpose() + (K * K.transpose()) * R;
  P_hat_ = 0.5f * (P_hat_ + P_hat_.transpose());
}

Row21 RioEskf::computeRadarH_(const Vec3& mu_r, const Vec3& w_nom) const {
  const Mat3 R_WI = x_.q_WI.toRotationMatrix();
  const Mat3 R_IW = R_WI.transpose();
  const Mat3 R_IR = x_.q_IR.toRotationMatrix();
  const Mat3 R_RI = R_IR.transpose();
  const Vec3& p_IR = x_.p_IR;

  const Vec3 v_I = R_IW * x_.v_WI;

  Row21 H = Row21::Zero();

  // d e / d v_W = -μ^T R_RI R_IW
  H.block<1,3>(0, 3) = -(mu_r.transpose() * (R_RI * R_IW));

  // d e / d δθ = -μ^T R_RI [R_IW v_W]_x = -μ^T R_RI [v_I]_x
  H.block<1,3>(0, 9) = -(mu_r.transpose() * (R_RI * skew(v_I)));

  // d e / d b_g = -μ^T R_RI [p_IR]_x
  H.block<1,3>(0, 12) = -(mu_r.transpose() * (R_RI * skew(p_IR)));

  // d e / d p_IR = -μ^T R_RI [w_nom]_x
  H.block<1,3>(0, 15) = -(mu_r.transpose() * R_RI) * skew(w_nom);

  // d e / d δθ_IR = -μ^T [v_R_nom]_x
  const Vec3 s_I = v_I + skew(w_nom) * p_IR;
  const Vec3 v_R_nom = R_RI * s_I;
  H.block<1,3>(0, 18) = -(mu_r.transpose() * skew(v_R_nom));

  return H;
}

float RioEskf::computeRadarh_(const Vec3& mu_r, const Vec3& w_nom) const {
  const Mat3 R_WI = x_.q_WI.toRotationMatrix();
  const Mat3 R_IW = R_WI.transpose();
  const Mat3 R_IR = x_.q_IR.toRotationMatrix();
  const Mat3 R_RI = R_IR.transpose();
  const Vec3& p_IR = x_.p_IR;

  const Vec3 v_I  = R_IW * x_.v_WI;
  const Vec3 spin = w_nom.cross(p_IR);
  const Vec3 v_R  = R_RI * (v_I + spin);
  return -(mu_r.dot(v_R));
}

void RioEskf::updateStateEstimate(const Vec21& delta_x) {
  // Inject error-state correction into nominal state
  x_.p_WI += delta_x.segment<3>(0);
  x_.v_WI += delta_x.segment<3>(3);
  x_.b_a  += delta_x.segment<3>(6);
  x_.b_g  += delta_x.segment<3>(12);

  const Vec3 dth = delta_x.segment<3>(9);
  const Quat dq = quatExpSmall(dth);
  x_.q_WI = (x_.q_WI * dq).normalized();

  // Extrinsic corrections
  x_.p_IR += delta_x.segment<3>(15);
  const Vec3 dth_IR = delta_x.segment<3>(18);
  const Quat dq_IR = quatExpSmall(dth_IR);
  x_.q_IR = (x_.q_IR * dq_IR).normalized();
}

void RioEskf::setExtrinsics(const Vec3& p_IR, const Quat& q_IR) {
  if (!initialized_) return;
  x_.p_IR = p_IR;
  x_.q_IR = q_IR.normalized();
  // Zero the extrinsic slots of the accumulated error-state so that the new
  // nominal values are the reference point and no stale correction is injected.
  delta_x_hat_.segment<3>(15).setZero();
  delta_x_hat_.segment<3>(18).setZero();
  delta_x_hat_prior_.segment<3>(15).setZero();
  delta_x_hat_prior_.segment<3>(18).setZero();
}

void RioEskf::advancePriorToPosterior() {
  P_hat_       = P_hat_prior_;
  delta_x_hat_ = delta_x_hat_prior_;
}

const Mat21 RioEskf::generateA(Vec3 f_nom, Vec3 w_nom) const {
  Mat21 A = Mat21::Zero();

  // Top-left 15x15 block
  // Row 0: [O3, I3, O3, O3, O3]
  A.block<3,3>(0, 3) = Mat3::Identity();

  // Row 1: [O3, O3, -R_nb, -R_nb @ skew(f_nom), O3]
  const Mat3 R_nb = x_.q_WI.toRotationMatrix();
  A.block<3,3>(3, 6)  = -R_nb;
  A.block<3,3>(3, 9)  = -R_nb * skew(f_nom);

  // Row 2: [O3, O3, -(1/T_acc)*I3, O3, O3]
  if (params_.tau_ba > 0.0f)
    A.block<3,3>(6, 6) = -(1.0f / params_.tau_ba) * Mat3::Identity();

  // Row 3: [O3, O3, O3, -skew(w_nom), -I3]
  A.block<3,3>(9, 9)  = -skew(w_nom);
  A.block<3,3>(9, 12) = -Mat3::Identity();

  // Row 4: [O3, O3, O3, O3, -(1/T_ars)*I3]
  if (params_.tau_bg > 0.0f)
    A.block<3,3>(12, 12) = -(1.0f / params_.tau_bg) * Mat3::Identity();

  // Rows/cols 15:21 remain zero (no dynamics for extrinsics)
  return A;
}

const Mat21x12 RioEskf::generateE() const {
  Mat21x12 E = Mat21x12::Zero();

  const Mat3 R_nb = x_.q_WI.toRotationMatrix();

  // Top-left 15x12 block
  // Row 0 (dp):  [O3, O3, O3, O3]
  // Row 1 (dv):  [-R_nb, O3, O3, O3]
  E.block<3,3>(3, 0) = -R_nb;

  // Row 2 (dba): [O3, I3, O3, O3]
  E.block<3,3>(6, 3) = Mat3::Identity();

  // Row 3 (dth): [O3, O3, -I3, O3]
  E.block<3,3>(9, 6) = -Mat3::Identity();

  // Row 4 (dbg): [O3, O3, O3, I3]
  E.block<3,3>(12, 9) = Mat3::Identity();

  // Rows 15:21 remain zero (no noise on extrinsics)
  return E;
}

} // namespace rio