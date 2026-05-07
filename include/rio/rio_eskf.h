#pragma once

#include <stddef.h>
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

  float sigma_vr = 0.1f;
  bool  gating_enable = true;
  float gate_nsigma = 3.0f;
  float vr_sign = -1.0f;

  // Differential barometer aiding (z only).
  // sigma_baro_dz is the std-dev (m) of one Δz measurement (combines noise
  // of two pressure samples and short-term local pressure disturbances).
  float sigma_baro_dz       = 0.3f;
  bool  baro_gating_enable  = true;
  float baro_gate_nsigma    = 5.0f;
  // Sign convention: world-frame z is "up" if g_W.z() < 0 (default), and
  // increasing altitude increases p_WI.z(). Set to -1.0f if z is "down".
  float baro_z_sign         = 1.0f;
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

struct CorrectionResult {
  size_t n_total    = 0;   // measurements received
  size_t n_accepted = 0;   // passed gating and used for update
  size_t n_rejected = 0;   // failed gating (chi^2 too large)
  size_t n_skipped  = 0;   // skipped (zero-norm direction, S<=0, etc.)
};

struct BaroCorrectionResult {
  bool  initialized = false;  // anchor was just set on this call (no update)
  bool  accepted    = false;
  bool  rejected    = false;  // gating
  bool  skipped     = false;  // bad data / not initialized
  float dz_meas     = 0.0f;   // Δz from pressure (m)
  float dz_pred     = 0.0f;   // Δz from state (m)
  float residual    = 0.0f;
};

// Error-state ordering:
// δx = [ δp^n, δv^n, δb_acc^b, δθ_nb, δb_ars^b ]^T
class RioEskf {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  RioEskf();

  void setParams(const Params& p);
  bool paramsSet() const;

  void reset(const NominalState& x0, const float* P0_diag_21, float t0);

  /// Attempt to initialize attitude from a gravity-aligned accelerometer reading.
  /// @param f_b       Specific-force measurement in body frame (m/s²).
  /// @param P0_diag   Pointer to 21-element initial covariance diagonal.
  /// @param t0        Timestamp to seed the filter with.
  /// @param g_tol     Allowed deviation of |f_b| from |g_W| (m/s²).  Default 0.8.
  /// @return true if |f_b| was close enough to gravity and the filter was reset.
  bool initAttitudeFromGravity(const Vec3& f_b, const float* P0_diag, float t0,
                               float g_tol = 0.8f);

  bool isInitialized() const;
  float lastTime() const;

  const NominalState& getState() const;
  const Mat21& getCovariance() const;

  void predict(const ImuSample& s, float dt);
  void insPropagation(const ImuSample& s, float dt);
  CorrectionResult correct(const RadarDoppler* meas, size_t n, const ImuSample& s);
  BaroCorrectionResult correctBarometer(const BarometerSample& s);
  void resetBarometer();
  void updateStateEstimate(const Vec21& delta_x);
  void advancePriorToPosterior();

private:
  void scalarCorrect_(const Row21& H, float residual, float R);
  Row21 computeRadarH_(const Vec3& mu_r, const Vec3& w_nom) const;
  float computeRadarh_(const Vec3& mu_r, const Vec3& w_nom) const;

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

  // Differential barometer anchor: pressure and state z at last accepted
  // (or initializing) barometer reading.
  bool  baro_has_prev_{false};
  float baro_p_prev_{0.0f};
  float baro_z_prev_{0.0f};
};

} // namespace rio