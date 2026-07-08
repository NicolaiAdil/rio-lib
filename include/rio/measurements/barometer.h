#pragma once

#include "../measurement.h"

namespace rio {

// Barometric altimetry factor — BRIO Eq. 11 + 12 (Girod et al.,
// arXiv:2408.05764). Anchors (z_baro_anchor, z_state_anchor) on the
// first valid reading. Params::reset_anchor_on_accept selects mode:
//   true  → differential (re-anchor every accept; constrains Δz).
//   false → absolute, BRIO default (constrains x.z against z_p^0).
class BarometerMeasurement : public Measurement {
public:
  struct Params {
    float sigma_dz              = 0.3f;
    float z_sign                = 1.0f;
    float gate_nsigma           = 5.0f;
    bool  gating                = true;
    bool  reset_anchor_on_accept = true;
  };

  explicit BarometerMeasurement(const Params& p) noexcept : p_(p) {}

  void setSample(const BarometerSample& s);
  void resetAnchor();

  Eval evaluate(const NominalState& x, const MeasurementContext& ctx,
                Row21& H, float& e, float& R) override;
  void onAccepted(const NominalState& x_post) override;

  bool  gatingEnabled() const override { return p_.gating; }
  float gateNSigma()    const override { return p_.gate_nsigma; }

  float lastDzMeas()      const { return last_dz_meas_; }
  float lastDzPred()      const { return last_dz_pred_; }
  bool  justInitialized() const { return just_initialized_; }
  bool  hasAnchor()       const { return has_anchor_; }

private:
  Params p_;

  bool   has_pending_     = false;
  float  pending_p_pa_    = 0.f;
  float  pending_temp_c_  = 0.f;

  // Anchor: absolute altitude (NASA model) and state z at anchor time.
  bool   has_anchor_      = false;
  float  z_baro_anchor_   = 0.f;
  float  z_state_anchor_  = 0.f;

  // Cached pressureToAltitude(pending_p_pa_) so onAccepted can re-anchor.
  float  pending_z_baro_  = 0.f;

  float  last_dz_meas_    = 0.f;
  float  last_dz_pred_    = 0.f;
  bool   just_initialized_ = false;
};

}  // namespace rio
