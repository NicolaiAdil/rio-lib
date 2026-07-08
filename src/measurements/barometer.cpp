#include "rio/measurements/barometer.h"
#include "rio/rio_eskf.h"

#include <cmath>

namespace rio {
namespace {

// NASA Earth Atmosphere Model (Troposphere, valid to ~11 km): absolute
// pressure (Pa) → altitude (m). Eq. 12 of Girod et al., arXiv:2408.05764.
float pressureToAltitude(float p_pa) {
  return (288.08f * powf(p_pa / 101290.0f, 1.0f / 5.256f)
          - 273.1f - 15.04f) / -0.00649f;
}

}  // namespace

void BarometerMeasurement::setSample(const BarometerSample& s) {
  has_pending_     = true;
  pending_p_pa_    = s.pressure_pa;
  pending_temp_c_  = s.temp_c;
}

void BarometerMeasurement::resetAnchor() {
  has_anchor_     = false;
  z_baro_anchor_  = 0.f;
  z_state_anchor_ = 0.f;
}

Measurement::Eval
BarometerMeasurement::evaluate(const NominalState& x,
                                const MeasurementContext& /*ctx*/,
                                Row21& H, float& e, float& R) {
  last_dz_meas_     = 0.f;
  last_dz_pred_     = 0.f;
  just_initialized_ = false;

  if (!has_pending_) return Eval::NotReady;
  has_pending_ = false;

  if (!(pending_p_pa_ > 1.0f) || !std::isfinite(pending_p_pa_) ||
      !std::isfinite(pending_temp_c_)) {
    return Eval::Skip;
  }

  // First valid reading sets the anchor; no update this call.
  if (!has_anchor_) {
    has_anchor_       = true;
    z_baro_anchor_    = pressureToAltitude(pending_p_pa_);
    z_state_anchor_   = x.p_WI.z();
    just_initialized_ = true;
    return Eval::NotReady;
  }

  // Residual = sign-of-paper-r_B inverted (BRIO Eq. 11); Jacobian matches.
  const float z_baro_i = pressureToAltitude(pending_p_pa_);
  pending_z_baro_      = z_baro_i;
  const float dz_meas  = p_.z_sign * (z_baro_i - z_baro_anchor_);
  const float dz_pred  = x.p_WI.z() - z_state_anchor_;

  e = dz_meas - dz_pred;

  H = Row21::Zero();
  H(0, 2) = 1.0f;

  R = p_.sigma_dz * p_.sigma_dz;

  last_dz_meas_ = dz_meas;
  last_dz_pred_ = dz_pred;
  return Eval::Apply;
}

void BarometerMeasurement::onAccepted(const NominalState& x_post) {
  // Differential: re-anchor every accept. Absolute (BRIO default): never.
  if (!p_.reset_anchor_on_accept) return;
  z_baro_anchor_  = pending_z_baro_;
  z_state_anchor_ = x_post.p_WI.z();
}

}  // namespace rio
