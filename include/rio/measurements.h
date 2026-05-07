#pragma once

#include <stdint.h>
#include <math.h>

#include "rio_types.h"

namespace rio {

// Forward decl — defined in rio_eskf.h.
struct ImuSample;
struct BarometerSample;
struct NominalState;

// Auxiliary data for measurement evaluation. Extend as new modalities need
// more context (e.g. lever-arm offsets, current biases). Pointers are
// non-owning and must remain valid through applyScalar().
struct MeasurementContext {
  const ImuSample* last_imu = nullptr;
};

// Per-call diagnostic returned by RioEskf::applyScalar().
struct ScalarUpdate {
  enum Status : uint8_t {
    Accepted = 0,
    Rejected = 1,   // failed gating (chi^2 too large)
    Skipped  = 2,   // bad data, S<=0, or filter not initialized
    NotReady = 3,   // measurement reported it has nothing to update
  };
  Status status   = NotReady;
  float  residual = NAN;
  float  S        = NAN;
};

// Abstract scalar measurement (1-D update). Subclasses supply h, H, R for
// their modality. ESKF handles gating and the scalar Joseph update.
class ScalarMeasurement {
public:
  enum class Eval : uint8_t {
    Apply,    // H, e, R filled — ESKF should run gating + Joseph update
    Skip,     // bad data (mapped to ScalarUpdate::Skipped)
    NotReady, // nothing to apply this call (e.g. anchor not set yet)
  };

  virtual ~ScalarMeasurement() = default;

  // Compute Jacobian H (1x21), residual e = z - h(x), and noise R.
  // Return Eval::Apply only if all three are valid.
  virtual Eval evaluate(const NominalState& x,
                        const MeasurementContext& ctx,
                        Row21& H, float& e, float& R) = 0;

  // Optional post-accept hook. Called only when the update is accepted.
  // Use for differential anchors, residual statistics, etc.
  virtual void onAccepted(const NominalState& /*x_post*/) {}

  virtual bool  gatingEnabled() const { return true; }
  virtual float gateNSigma()    const { return 5.0f; }
};

// Per-point radar Doppler scalar measurement. One instance per radar point.
// Constructed cheaply each iteration; holds no anchor / accumulation state.
class RadarDopplerMeasurement : public ScalarMeasurement {
public:
  struct Params {
    float sigma_vr   = 0.038f;
    float vr_sign    = 1.0f;
    float gate_nsigma = 5.0f;
    bool  gating     = true;
  };

  RadarDopplerMeasurement(const Params& p, Vec3 u_R, float vr) noexcept
      : p_(p), u_R_(u_R), vr_(vr) {}

  Eval evaluate(const NominalState& x, const MeasurementContext& ctx,
                Row21& H, float& e, float& R) override;

  bool  gatingEnabled() const override { return p_.gating; }
  float gateNSigma()    const override { return p_.gate_nsigma; }

private:
  Params p_;
  Vec3   u_R_;
  float  vr_;
};

// Differential barometer (Δz) scalar measurement. Persistent instance —
// holds the (p_prev, z_prev) anchor across calls and re-anchors in
// onAccepted().
class BarometerDiffMeasurement : public ScalarMeasurement {
public:
  struct Params {
    float sigma_dz    = 0.3f;
    float z_sign      = 1.0f;
    float gate_nsigma = 5.0f;
    bool  gating      = true;
  };

  explicit BarometerDiffMeasurement(const Params& p) noexcept : p_(p) {}

  // Stage the next sample. evaluate() consumes it; subsequent calls without
  // a fresh setSample() return false (NotReady).
  void setSample(const BarometerSample& s);
  void resetAnchor();

  Eval evaluate(const NominalState& x, const MeasurementContext& ctx,
                Row21& H, float& e, float& R) override;
  void onAccepted(const NominalState& x_post) override;

  bool  gatingEnabled() const override { return p_.gating; }
  float gateNSigma()    const override { return p_.gate_nsigma; }

  // Sensor-specific telemetry from the last evaluate() call.
  float lastDzMeas()      const { return last_dz_meas_; }
  float lastDzPred()      const { return last_dz_pred_; }
  bool  justInitialized() const { return just_initialized_; }
  bool  hasAnchor()       const { return has_anchor_; }

private:
  Params p_;

  // Pending (caller-provided) sample.
  bool   has_pending_     = false;
  float  pending_p_pa_    = 0.f;
  float  pending_temp_c_  = 0.f;

  // Differential anchor.
  bool   has_anchor_      = false;
  float  p_prev_          = 0.f;
  float  z_prev_          = 0.f;

  // Telemetry from last evaluate().
  float  last_dz_meas_    = 0.f;
  float  last_dz_pred_    = 0.f;
  bool   just_initialized_ = false;
};

}  // namespace rio
