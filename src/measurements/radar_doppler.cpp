#include "rio/measurements/radar_doppler.h"
#include "rio/rio_eskf.h"

#include <cmath>

namespace rio {
namespace {

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

// Central-difference numerical Hessian of computeRadarH along each of the
// 21 error-state axes; symmetrized to remove FD round-off asymmetry that
// would otherwise bias B = ½ tr(H_xx P H_xx P).
Mat21 computeRadarHxxNumerical(const NominalState& x,
                               const Vec3& mu_r, const Vec3& w_nom) {
  // ε ≈ 1e-3 sits near the FD truncation / round-off knee for float Eigen.
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
  Hxx = 0.5f * (Hxx + Hxx.transpose());
  return Hxx;
}

}  // namespace

Measurement::Eval
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

  // B = ½ tr(M M) with M = H_xx · P, computed as Σ M[i,j]·M[j,i].
  const Mat21 M = Hxx * P;
  float B = 0.0f;
  for (int i = 0; i < 21; ++i)
    for (int j = 0; j < 21; ++j)
      B += M(i, j) * M(j, i);
  B *= 0.5f;
  // Underweighting must only inflate S; clip pathological negatives.
  return (B > 0.0f) ? B : 0.0f;
}

}  // namespace rio
