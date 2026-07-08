#pragma once

#include "../measurement.h"

namespace rio {

// Per-point radar Doppler factor. Constructed cheaply each iteration;
// stateless across calls.
class RadarDopplerMeasurement : public Measurement {
public:
  struct Params {
    float sigma_vr    = 0.038f;
    float vr_sign     = 1.0f;
    float gate_nsigma = 5.0f;
    bool  gating      = true;
    bool  underweight = false;  // enable 2nd-order term in computeB()
  };

  RadarDopplerMeasurement(const Params& p, Vec3 u_R, float vr) noexcept
      : p_(p), u_R_(u_R), vr_(vr) {}

  Eval evaluate(const NominalState& x, const MeasurementContext& ctx,
                Row21& H, float& e, float& R) override;
  float computeB(const NominalState& x,
                 const MeasurementContext& ctx) const override;

  bool  gatingEnabled() const override { return p_.gating; }
  float gateNSigma()    const override { return p_.gate_nsigma; }

private:
  Params p_;
  Vec3   u_R_;
  float  vr_;
};

}  // namespace rio
