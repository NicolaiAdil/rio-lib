#pragma once

#include <stddef.h>
#include <stdint.h>
#include "rio_types.h"

namespace rio {

struct ImuSample {
  float t;
  Vec3  acc;
  Vec3  gyr;
};

struct RadarDoppler {
  float t;
  Vec3  u_R;
  float vr;
  float sigma;
};

struct BarometerSample {
  float t;
  float pressure_pa;
  float temp_c;
};

// Process-model and dynamics parameters only. Sensor sigmas / gating /
// signs live on per-modality measurement classes (see measurements.h).
struct Params {
  Vec3 g_W = Vec3(0, 0, -9.80665f);

  float sigma_gyr = 0.01f;
  float sigma_acc = 0.10f;
  float sigma_bg  = 0.001f;
  float sigma_ba  = 0.010f;

  float tau_bg = -1.0f;
  float tau_ba = -1.0f;

  Quat q_IR = Quat::Identity();
  Vec3 p_IR = Vec3::Zero();

  float min_dt = 1e-4f;
  float max_dt = 0.05f;
};

struct NominalState {
  Vec3 p_WI = Vec3::Zero();
  Vec3 v_WI = Vec3::Zero();
  Vec3 b_a  = Vec3::Zero();
  Quat q_WI = Quat::Identity();
  Vec3 b_g  = Vec3::Zero();
  Vec3 p_IR = Vec3::Zero();
  Quat q_IR = Quat::Identity();
};

}  // namespace rio

// Bring in the measurement interface (ScalarMeasurement, ScalarUpdate,
// MeasurementContext) so callers only need rio_eskf.h.
#include "measurements.h"

namespace rio {

// Error-state ordering:
// δx = [ δp^n, δv^n, δb_acc^b, δθ_nb, δb_ars^b, δp_IR, δθ_IR ]^T
class RioEskf {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  RioEskf();

  void setParams(const Params& p);
  bool paramsSet() const;

  void reset(const NominalState& x0, const float* P0_diag_21, float t0);

  /// Initialize attitude from a gravity-aligned accelerometer reading.
  /// @param f_b       Specific-force measurement in body frame (m/s²).
  /// @param P0_diag   Pointer to 21-element initial covariance diagonal.
  /// @param t0        Timestamp to seed the filter with.
  /// @param g_tol     Allowed deviation of |f_b| from |g_W| (m/s²).
  /// @return true if |f_b| was close enough to gravity and the filter was reset.
  bool initAttitudeFromGravity(const Vec3& f_b, const float* P0_diag, float t0,
                               float g_tol = 0.8f);

  bool isInitialized() const;
  float lastTime() const;

  const NominalState& getState() const;
  const Mat21& getCovariance() const;

  void predict(const ImuSample& s, float dt);
  void insPropagation(const ImuSample& s, float dt);

  /// Apply one scalar measurement update. Generic over modality —
  /// subclasses of ScalarMeasurement supply h, H, R via evaluate().
  /// On reject/skip/not-ready, P_hat_ and δx are NOT touched; if the
  /// caller wants P_hat_ snapped to the predicted prior (for diagnostics
  /// after a batch with no accepts), call advancePriorToPosterior().
  ScalarUpdate applyScalar(ScalarMeasurement& m,
                           const MeasurementContext& ctx = {});

  void updateStateEstimate(const Vec21& delta_x);
  void advancePriorToPosterior();

private:
  void scalarCorrect_(const Row21& H, float residual, float R);

  const Mat21 generateA(Vec3 f_nom, Vec3 w_nom) const;
  const Mat21x12 generateE() const;

  Params params_{};
  NominalState x_{};
  Mat12 Q_{Mat12::Zero()};

  Vec21 delta_x_hat_{Vec21::Zero()};
  Vec21 delta_x_hat_prior_{Vec21::Zero()};
  Mat21 P_hat_{Mat21::Zero()};
  Mat21 P_hat_prior_{Mat21::Zero()};

  bool params_set_{false};
  bool initialized_{false};
  float t_last_{0.0f};
};

} // namespace rio
