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

// Process-model and dynamics parameters. Sensor sigmas / gating / signs
// live on per-modality measurement classes.
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

// Apply a scalar perturbation eps along error-state axis i ∈ [0,21) to x.
// w_nom (cached gyr − b_g) is updated in lockstep when i lands in the
// gyro-bias block so downstream Jacobians see a consistent (gyr − b_g).
// Per-axis counterpart to RioEskf::updateStateEstimate — both encode the
// 21-dim error-state layout documented above RioEskf.
void perturbErrorState(NominalState& x, Vec3& w_nom, int i, float eps);

}  // namespace rio

#include "measurement.h"

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
  /// @param g_tol Allowed deviation of |f_b| from |g_W| (m/s²).
  bool initAttitudeFromGravity(const Vec3& f_b, const float* P0_diag, float t0,
                               float g_tol = 0.8f);

  bool isInitialized() const;
  float lastTime() const;

  const NominalState& getState() const;
  const Mat21& getCovariance() const;

  void predict(const ImuSample& s, float dt);
  void insPropagation(const ImuSample& s, float dt);

  /// Apply one measurement update. Dispatch is by runtime type of `m`
  /// (virtual evaluate/computeB/onAccepted). If gating/skip/not-ready
  /// fires, P and δx are unchanged — call advancePriorToPosterior() to
  /// snap P_hat_ to the predicted prior.
  MeasurementUpdate correct(Measurement& m,
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
