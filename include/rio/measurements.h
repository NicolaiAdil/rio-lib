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
  // Prior error-state covariance. Populated by the ESKF before calling
  // evaluate()/computeB() so measurements that want to apply second-order
  // (underweighting) corrections can read P without a side channel.
  // Null means "not provided" — measurements must fall back gracefully.
  const Mat21*     P_prior  = nullptr;
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
  // Second-order underweighting term added to the residual covariance:
  //   S = H P H^T + R + B   (Navigation Filter Best Practices §5.2.3)
  // 0 for measurements that don't apply underweighting. NaN until evaluated.
  float  B        = NAN;
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

  // Optional second-order underweighting term per NavFilter Best Practices §5.2.3:
  //   B = ½ tr(H_xx P H_xx P)
  // where H_xx is the Hessian of h(x). Added to the residual covariance:
  //   S = H P H^T + R + B
  // Default 0 means "no underweighting" — only override for measurements
  // where the nonlinear second-order terms in h are significant compared
  // to R (e.g. radar Doppler at large attitude uncertainty).
  // Called by the ESKF after evaluate() returns Apply, with ctx.P_prior set.
  virtual float computeB(const NominalState& /*x*/,
                         const MeasurementContext& /*ctx*/) const {
    return 0.0f;
  }

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
    float sigma_vr    = 0.038f;
    float vr_sign     = 1.0f;
    float gate_nsigma = 5.0f;
    bool  gating      = true;
    // When true, computeB() returns the second-order Hessian-based
    // underweighting term B = ½ tr(H_xx P H_xx P) per §5.2.3 of
    // NASA Navigation Filter Best Practices. When false, computeB
    // returns 0 (and the cheap base implementation path is taken).
    bool  underweight = false;
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

// Barometer scalar measurement. Persistent instance — holds the
// (p_prev, z_prev) anchor across calls. Two modes, selected via
// Params::reset_anchor_on_accept:
//   true  (default): differential. Anchor is reset to (pending_p, x_post.z)
//                    after every accepted update — measurement constrains
//                    Δz between consecutive samples.
//   false:           absolute. Anchor is set once on the first valid reading
//                    (typically right after attitude init, when x.z ≈ 0)
//                    and never moves — measurement constrains x.z against the
//                    boot-time pressure reference.
class BarometerDiffMeasurement : public ScalarMeasurement {
public:
  struct Params {
    float sigma_dz              = 0.3f;
    float z_sign                = 1.0f;
    float gate_nsigma           = 5.0f;
    bool  gating                = true;
    bool  reset_anchor_on_accept = true;  // true=differential, false=absolute
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
