#pragma once

#include <stdint.h>
#include <math.h>

#include "rio_types.h"

namespace rio {

struct ImuSample;
struct BarometerSample;
struct NominalState;

// Shared evaluation context. Pointers are non-owning and must outlive correct().
struct MeasurementContext {
  const ImuSample* last_imu = nullptr;
  const Mat21*     P_prior  = nullptr;  // set by the ESKF before evaluate()
};

struct MeasurementUpdate {
  enum Status : uint8_t {
    Accepted = 0,
    Rejected = 1,   // failed gating (chi^2 too large)
    Skipped  = 2,   // bad data, S<=0, or filter not initialized
    NotReady = 3,   // measurement reported it has nothing to update
  };
  Status status   = NotReady;
  float  residual = NAN;
  float  S        = NAN;
  float  B        = NAN;   // 2nd-order term added to S (0 unless overridden)
};

// Base class for scalar (1-D) measurements. ESKF runs gating + Joseph update.
class Measurement {
public:
  enum class Eval : uint8_t {
    Apply,    // H, e, R filled — ESKF runs gating + Joseph
    Skip,     // bad data (→ MeasurementUpdate::Skipped)
    NotReady, // nothing to apply this call (e.g. anchor not set yet)
  };

  virtual ~Measurement() = default;

  // Fill Jacobian H (1×21), residual e = z − h(x), and noise R.
  virtual Eval evaluate(const NominalState& x,
                        const MeasurementContext& ctx,
                        Row21& H, float& e, float& R) = 0;

  // Optional B = ½ tr(H_xx P H_xx P) added to S (NavFilter Best Practices
  // §5.2.3). See RadarDopplerMeasurement for an example.
  virtual float computeB(const NominalState& /*x*/,
                         const MeasurementContext& /*ctx*/) const {
    return 0.0f;
  }

  // Optional post-accept hook (e.g. re-anchoring for differential factors).
  virtual void onAccepted(const NominalState& /*x_post*/) {}

  virtual bool  gatingEnabled() const { return true; }
  virtual float gateNSigma()    const { return 5.0f; }
};

}  // namespace rio
