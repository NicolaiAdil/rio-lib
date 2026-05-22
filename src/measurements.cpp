// Per-modality scalar measurement implementations. Math kernels lifted
// verbatim from the previous monolithic rio_eskf.cpp so behavior is
// preserved bit-for-bit.

#include "rio/measurements.h"
#include "rio/rio_eskf.h"   // for NominalState, ImuSample, BarometerSample

#include <cmath>

namespace rio {
namespace {

// Radar Jacobian (1x21). Verbatim from old RioEskf::computeRadarH_.
Row21 computeRadarH(const NominalState& x, const Vec3& mu_r, const Vec3& w_nom) {
  const Mat3 R_WI = x.q_WI.toRotationMatrix();
  const Mat3 R_IW = R_WI.transpose();
  const Mat3 R_IR = x.q_IR.toRotationMatrix();
  const Mat3 R_RI = R_IR.transpose();
  const Vec3& p_IR = x.p_IR;

  const Vec3 v_I = R_IW * x.v_WI;

  Row21 H = Row21::Zero();

  H.block<1,3>(0, 3)  = -(mu_r.transpose() * (R_RI * R_IW));
  H.block<1,3>(0, 9)  = -(mu_r.transpose() * (R_RI * skew(v_I)));
  H.block<1,3>(0, 12) = -(mu_r.transpose() * (R_RI * skew(p_IR)));
  H.block<1,3>(0, 15) = -(mu_r.transpose() * R_RI) * skew(w_nom);

  const Vec3 s_I     = v_I + skew(w_nom) * p_IR;
  const Vec3 v_R_nom = R_RI * s_I;
  H.block<1,3>(0, 18) = -(mu_r.transpose() * skew(v_R_nom));

  return H;
}

// Predicted radial velocity. Verbatim from old RioEskf::computeRadarh_.
float computeRadarh(const NominalState& x, const Vec3& mu_r, const Vec3& w_nom) {
  const Mat3 R_WI = x.q_WI.toRotationMatrix();
  const Mat3 R_IW = R_WI.transpose();
  const Mat3 R_IR = x.q_IR.toRotationMatrix();
  const Mat3 R_RI = R_IR.transpose();
  const Vec3& p_IR = x.p_IR;

  const Vec3 v_I  = R_IW * x.v_WI;
  const Vec3 spin = w_nom.cross(p_IR);
  const Vec3 v_R  = R_RI * (v_I + spin);
  return -(mu_r.dot(v_R));
}

// Perturb the nominal state by ε along error-state direction i ∈ [0,21).
// Mirrors the error-state layout used by the ESKF covariance:
//   0:3   δp      (position, additive)
//   3:6   δv      (velocity, additive)
//   6:9   δb_a    (accel bias, additive)
//   9:12  δθ_WI   (attitude, right-multiplicative quaternion perturbation)
//   12:15 δb_g    (gyro bias, additive — also shifts w_nom = gyr − b_g)
//   15:18 δp_IR   (radar lever arm, additive)
//   18:21 δθ_IR   (radar attitude offset, right-multiplicative)
// w_nom is updated in place so the perturbed Jacobian sees a consistent
// (gyr − b_g) — the ESKF caches w_nom once per evaluate() and passes it
// into computeRadarH().
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
    w_nom[i - 12] -= eps;       // w_nom = gyr − b_g
  } else if (i < 18) {
    x.p_IR[i - 15] += eps;
  } else {
    Vec3 dtheta = Vec3::Zero();
    dtheta[i - 18] = eps;
    x.q_IR = (x.q_IR * quatExpSmall(dtheta)).normalized();
  }
}

// Numerical Hessian H_xx (21×21) of the scalar radar Doppler measurement
// w.r.t. the error state, computed by central finite differences of the
// analytical Jacobian along each of the 21 error-state directions. Exact
// to within FD truncation O(ε²) and round-off; avoids the tedium of
// deriving and verifying eight separate analytical Hessian blocks.
// Cost per call: 42 evaluations of computeRadarH. Only invoked when
// Params::underweight is true, so non-underweighted radar updates pay nothing.
Mat21 computeRadarHxxNumerical(const NominalState& x,
                               const Vec3& mu_r, const Vec3& w_nom) {
  // ε balances FD truncation (∝ ε²) against round-off (∝ 1/ε). For
  // single-precision Eigen floats with O(1) state components, ε ≈ 1e-3
  // sits near the round-off / truncation knee.
  constexpr float eps = 1e-3f;
  constexpr float inv_2eps = 1.0f / (2.0f * eps);

  Mat21 Hxx = Mat21::Zero();
  for (int i = 0; i < 21; ++i) {
    NominalState x_p = x; Vec3 w_p = w_nom; perturbErrorState(x_p, w_p, i,  eps);
    NominalState x_m = x; Vec3 w_m = w_nom; perturbErrorState(x_m, w_m, i, -eps);
    const Row21 dH = (computeRadarH(x_p, mu_r, w_p)
                    - computeRadarH(x_m, mu_r, w_m)) * inv_2eps;
    Hxx.col(i) = dH.transpose();
  }
  // True Hessian is symmetric; FD picks up small asymmetric round-off
  // that would bias B = ½ tr(H_xx P H_xx P) if left in place.
  Hxx = 0.5f * (Hxx + Hxx.transpose());
  return Hxx;
}

}  // namespace

// ── RadarDopplerMeasurement ───────────────────────────────────────────────────

ScalarMeasurement::Eval
RadarDopplerMeasurement::evaluate(const NominalState& x,
                                   const MeasurementContext& ctx,
                                   Row21& H, float& e, float& R) {
  if (!ctx.last_imu) return Eval::Skip;

  Vec3 mu_r = u_R_;
  const float un = mu_r.norm();
  // !(un >= …) also rejects NaN.
  if (!(un >= 1e-6f)) return Eval::Skip;
  mu_r /= un;

  const Vec3 w_nom = ctx.last_imu->gyr - x.b_g;

  H = computeRadarH(x, mu_r, w_nom);
  const float h = computeRadarh(x, mu_r, w_nom);
  e = p_.vr_sign * vr_ - h;
  R = p_.sigma_vr * p_.sigma_vr;
  return Eval::Apply;
}

// Second-order underweighting B = ½ tr(H_xx P H_xx P), with H_xx the
// numerical Hessian of h(x). Returns 0 when underweighting is disabled
// in Params or when the ESKF didn't provide a prior covariance.
float
RadarDopplerMeasurement::computeB(const NominalState& x,
                                  const MeasurementContext& ctx) const {
  if (!p_.underweight || !ctx.last_imu || !ctx.P_prior) return 0.0f;

  Vec3 mu_r = u_R_;
  const float un = mu_r.norm();
  if (!(un >= 1e-6f)) return 0.0f;
  mu_r /= un;

  const Vec3 w_nom = ctx.last_imu->gyr - x.b_g;

  const Mat21 Hxx = computeRadarHxxNumerical(x, mu_r, w_nom);
  const Mat21& P  = *ctx.P_prior;

  // B = ½ tr(M M)  where M = H_xx · P. Implemented as sum_{i,j} M[i,j]·M[j,i]
  // so we only need one matmul + one Frobenius-style sum.
  const Mat21 M = Hxx * P;
  float B = 0.0f;
  for (int i = 0; i < 21; ++i)
    for (int j = 0; j < 21; ++j)
      B += M(i, j) * M(j, i);
  B *= 0.5f;
  // Guard: H_xx · P is finite by construction, but a non-positive-definite
  // P (which shouldn't happen with UDU updates but could in pathological
  // cases) could in principle yield a negative B that would shrink S below
  // R. Clip to zero — underweighting should only inflate W, never shrink.
  return (B > 0.0f) ? B : 0.0f;
}

// ── BarometerDiffMeasurement ──────────────────────────────────────────────────

void BarometerDiffMeasurement::setSample(const BarometerSample& s) {
  has_pending_     = true;
  pending_p_pa_    = s.pressure_pa;
  pending_temp_c_  = s.temp_c;
}

void BarometerDiffMeasurement::resetAnchor() {
  has_anchor_ = false;
  p_prev_     = 0.f;
  z_prev_     = 0.f;
}

ScalarMeasurement::Eval
BarometerDiffMeasurement::evaluate(const NominalState& x,
                                    const MeasurementContext& /*ctx*/,
                                    Row21& H, float& e, float& R) {
  // Reset per-call telemetry.
  last_dz_meas_     = 0.f;
  last_dz_pred_     = 0.f;
  just_initialized_ = false;

  if (!has_pending_) return Eval::NotReady;
  // Consume pending sample.
  has_pending_ = false;

  if (!(pending_p_pa_ > 1.0f) || !std::isfinite(pending_p_pa_) ||
      !std::isfinite(pending_temp_c_)) {
    return Eval::Skip;
  }

  // Anchor on first valid reading; signal NotReady so caller knows there
  // was no update this call (matches old "initialized=true, accepted=false").
  if (!has_anchor_) {
    has_anchor_       = true;
    p_prev_           = pending_p_pa_;
    z_prev_           = x.p_WI.z();
    just_initialized_ = true;
    return Eval::NotReady;
  }

  // Δz from pressure (hypsometric, with local temperature).
  const float T_kelvin = pending_temp_c_ + 273.15f;
  const float dz_baro  = differentialAltitude(p_prev_, pending_p_pa_, T_kelvin);
  const float dz_meas  = p_.z_sign * dz_baro;

  // Predicted Δz from current state. z_prev_ is a frozen snapshot of state
  // z at last accept (or anchor init) — H only sees current state z.
  const float dz_pred = x.p_WI.z() - z_prev_;

  e = dz_meas - dz_pred;

  H = Row21::Zero();
  H(0, 2) = 1.0f;

  R = p_.sigma_dz * p_.sigma_dz;

  last_dz_meas_ = dz_meas;
  last_dz_pred_ = dz_pred;
  return Eval::Apply;
}

void BarometerDiffMeasurement::onAccepted(const NominalState& x_post) {
  // Differential mode (default): re-anchor to the pressure consumed in the
  // last evaluate() and the posterior state z. evaluate() clears
  // has_pending_ but leaves pending_p_pa_ holding the consumed value, so we
  // can read it here.
  // Absolute mode: leave the anchor at its boot-time value so each
  // measurement constrains x.z against the same reference.
  if (!p_.reset_anchor_on_accept) return;
  p_prev_ = pending_p_pa_;
  z_prev_ = x_post.p_WI.z();
}

}  // namespace rio
