#include "rio/rio_eskf.h"

namespace rio {

RioEskf::RioEskf() = default;

void RioEskf::setParams(const Params& p) {
  params_ = p;
  params_.q_IR.normalize();

  // Q column ordering matches generateE: [acc, ba, gyr, bg].
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

  const float g_mag = params_.g_W.norm();
  const float fn    = f_b.norm();

  if (fn < (g_mag - g_tol) || fn > (g_mag + g_tol)) return false;

  const Vec3 gb = f_b / std::max(1e-6f, fn);

  const float roll  = std::atan2(gb.y(), gb.z());
  const float pitch = std::atan2(-gb.x(),
                        std::sqrt(gb.y() * gb.y() + gb.z() * gb.z()));
  const float yaw   = 0.0f;

  const float cr = std::cos(roll  * 0.5f), sr = std::sin(roll  * 0.5f);
  const float cp = std::cos(pitch * 0.5f), sp = std::sin(pitch * 0.5f);
  const float cy = std::cos(yaw   * 0.5f), sy = std::sin(yaw   * 0.5f);

  Quat q_WI;
  q_WI.w() = cr * cp * cy + sr * sp * sy;
  q_WI.x() = sr * cp * cy - cr * sp * sy;
  q_WI.y() = cr * sp * cy + sr * cp * sy;
  q_WI.z() = cr * cp * sy - sr * sp * cy;
  q_WI.normalize();

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

#if defined(RIO_EIGEN_NO_MALLOC)
  Eigen::internal::set_is_malloc_allowed(false);
#endif

  const Vec3 w_nom = s.gyr - x_.b_g;
  const Vec3 f_nom = s.acc - x_.b_a;

  const Mat21 A = generateA(f_nom, w_nom);
  const Mat21 Ad = Mat21::Identity() + A * dt;

  const Mat21x12 E = generateE();
  Mat12 Qc = Q_;
  Qc.block<3,3>(0, 0) /= dt;
  Qc.block<3,3>(3, 3) *= dt;
  Qc.block<3,3>(6, 6) /= dt;
  Qc.block<3,3>(9, 9) *= dt;

  const Mat21 Qd = E * Qc * E.transpose() * dt;

  delta_x_hat_prior_ = Ad * delta_x_hat_;
  P_hat_prior_ = Ad * P_hat_ * Ad.transpose() + Qd;

  t_last_ = s.t;

#if defined(RIO_EIGEN_NO_MALLOC)
  Eigen::internal::set_is_malloc_allowed(true);
#endif
}

void RioEskf::insPropagation(const ImuSample& s, float dt) {
  const auto f_nom = s.acc - x_.b_a;
  const auto w_nom = s.gyr - x_.b_g;
  const Mat3 R_WI = x_.q_WI.toRotationMatrix();

  const Vec3 a_W = R_WI * f_nom + params_.g_W;

  x_.p_WI += x_.v_WI * dt + 0.5f * a_W * dt * dt;
  x_.v_WI += a_W * dt;

  const Quat dq = quatExpSmall(w_nom * dt);
  x_.q_WI = (x_.q_WI * dq).normalized();
}

MeasurementUpdate RioEskf::correct(Measurement& m,
                                    const MeasurementContext& ctx) {
  MeasurementUpdate u;

  if (!params_set_ || !initialized_) {
    u.status = MeasurementUpdate::Skipped;
    return u;
  }

  // Expose P_prior so measurements can compute second-order terms.
  MeasurementContext ctx_local = ctx;
  ctx_local.P_prior = &P_hat_prior_;

  Row21 H;
  float e = 0.f, R = 0.f;
  const Measurement::Eval ev = m.evaluate(x_, ctx_local, H, e, R);

  if (ev == Measurement::Eval::NotReady) {
    u.status = MeasurementUpdate::NotReady;
    return u;
  }
  if (ev == Measurement::Eval::Skip) {
    u.status = MeasurementUpdate::Skipped;
    return u;
  }

  const float B = m.computeB(x_, ctx_local);

  // Record S unconditionally so callers can log rejects/skips too.
  const float S = (H * P_hat_prior_ * H.transpose())(0, 0) + R + B;
  u.residual = e;
  u.S        = S;
  u.B        = B;

  // !(S>0) and !isfinite(e) catch NaN — a NaN residual would otherwise
  // sneak past the chi-square gate and poison the state.
  if (!(S > 0.0f) || !isfinite(e)) {
    u.status = MeasurementUpdate::Rejected;
    return u;
  }

  if (m.gatingEnabled()) {
    const float gate_thresh = m.gateNSigma() * m.gateNSigma();
    if (!(e * e / S <= gate_thresh)) {
      u.status = MeasurementUpdate::Rejected;
      return u;
    }
  }

  // Joseph update with effective noise R+B (consistent with the gating S).
  scalarCorrect_(H, e, R + B);
  updateStateEstimate(delta_x_hat_);

  // Carry posterior into prior so the next correct() (e.g. next radar
  // point in the batch) sees the latest covariance.
  P_hat_prior_       = P_hat_;
  delta_x_hat_prior_ = delta_x_hat_;

  m.onAccepted(x_);

  u.status = MeasurementUpdate::Accepted;
  return u;
}

void RioEskf::scalarCorrect_(const Row21& H, float residual, float R) {
  const float S = (H * P_hat_prior_ * H.transpose())(0, 0) + R;
  if (!(S > 1e-12f)) return;
  const float invS = 1.0f / S;

  const Vec21 PHt = P_hat_prior_ * H.transpose();
  const Vec21 K   = PHt * invS;

  delta_x_hat_ = K * residual;

  const Mat21 IKH = Mat21::Identity() - K * H;
  P_hat_ = IKH * P_hat_prior_ * IKH.transpose() + (K * K.transpose()) * R;
  P_hat_ = 0.5f * (P_hat_ + P_hat_.transpose());
}

void RioEskf::updateStateEstimate(const Vec21& delta_x) {
  x_.p_WI += delta_x.segment<3>(0);
  x_.v_WI += delta_x.segment<3>(3);
  x_.b_a  += delta_x.segment<3>(6);
  x_.b_g  += delta_x.segment<3>(12);

  const Vec3 dth = delta_x.segment<3>(9);
  const Quat dq = quatExpSmall(dth);
  x_.q_WI = (x_.q_WI * dq).normalized();

  x_.p_IR += delta_x.segment<3>(15);
  const Vec3 dth_IR = delta_x.segment<3>(18);
  const Quat dq_IR = quatExpSmall(dth_IR);
  x_.q_IR = (x_.q_IR * dq_IR).normalized();
}

void perturbErrorState(NominalState& x, Vec3& w_nom, int i, float eps) {
  if (i < 3) {
    x.p_WI[i] += eps;
  } else if (i < 6) {
    x.v_WI[i - 3] += eps;
  } else if (i < 9) {
    x.b_a[i - 6] += eps;
  } else if (i < 12) {
    Vec3 dtheta = Vec3::Zero();
    dtheta[i - 9] = eps;
    x.q_WI = (x.q_WI * quatExpSmall(dtheta)).normalized();
  } else if (i < 15) {
    x.b_g[i - 12] += eps;
    w_nom[i - 12] -= eps;
  } else if (i < 18) {
    x.p_IR[i - 15] += eps;
  } else {
    Vec3 dtheta = Vec3::Zero();
    dtheta[i - 18] = eps;
    x.q_IR = (x.q_IR * quatExpSmall(dtheta)).normalized();
  }
}

void RioEskf::advancePriorToPosterior() {
  P_hat_       = P_hat_prior_;
  delta_x_hat_ = delta_x_hat_prior_;
}

const Mat21 RioEskf::generateA(Vec3 f_nom, Vec3 w_nom) const {
  Mat21 A = Mat21::Zero();

  A.block<3,3>(0, 3) = Mat3::Identity();

  const Mat3 R_nb = x_.q_WI.toRotationMatrix();
  A.block<3,3>(3, 6)  = -R_nb;
  A.block<3,3>(3, 9)  = -R_nb * skew(f_nom);

  if (params_.tau_ba > 0.0f)
    A.block<3,3>(6, 6) = -(1.0f / params_.tau_ba) * Mat3::Identity();

  A.block<3,3>(9, 9)  = -skew(w_nom);
  A.block<3,3>(9, 12) = -Mat3::Identity();

  if (params_.tau_bg > 0.0f)
    A.block<3,3>(12, 12) = -(1.0f / params_.tau_bg) * Mat3::Identity();

  return A;
}

const Mat21x12 RioEskf::generateE() const {
  Mat21x12 E = Mat21x12::Zero();

  const Mat3 R_nb = x_.q_WI.toRotationMatrix();

  E.block<3,3>(3, 0) = -R_nb;
  E.block<3,3>(6, 3) = Mat3::Identity();
  E.block<3,3>(9, 6) = -Mat3::Identity();
  E.block<3,3>(12, 9) = Mat3::Identity();

  return E;
}

} // namespace rio
