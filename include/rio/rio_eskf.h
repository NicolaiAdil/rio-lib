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
  Vec3  u_R;
  float vr;
  float sigma;
};

struct RadarVel3 {
  Vec3 v_R;
  Vec3 sigma;
};

struct Params {
  Vec3 g_W = Vec3(0, 0, -9.80665f);

  float sigma_gyr = 0.01f;
  float sigma_acc = 0.10f;
  float sigma_bg  = 0.001f;
  float sigma_ba  = 0.010f;

  float tau_bg = -1.0f;
  float tau_ba = -1.0f;

  Quat q_RB = Quat::Identity();
  Vec3 p_BR_B = Vec3::Zero();

  float min_dt = 1e-4f;
  float max_dt = 0.05f;
};

struct NominalState {
  Vec3 p_WB = Vec3::Zero();
  Vec3 v_WB = Vec3::Zero();
  Quat q_WB = Quat::Identity();
  Vec3 b_a  = Vec3::Zero();
  Vec3 b_g  = Vec3::Zero();
};

// Error-state ordering:
// δx = [ δp^n, δv^n, δb_acc^b, δθ_nb, δb_ars^b ]^T
class RioEskf {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  RioEskf();

  void setParams(const Params& p);
  bool paramsSet() const;

  void reset(const NominalState& x0, const float* P0_diag_15, float t0);

  bool isInitialized() const;
  float lastTime() const;

  const NominalState& state() const;
  const Mat15& covariance() const;

  void propagate(const ImuSample& s);
  void updateDoppler(const RadarDoppler* meas, size_t n);
  void updateDopplerWithOmega(const RadarDoppler* meas, size_t n, const Vec3& omega_B);
  void updateRadarVelocity3(const RadarVel3& z);

private:
  void scalarUpdate_(const Row15& H, float residual, float R);
  void inject_(const Vec15& dx);

  Params params_{};
  NominalState x_{};
  Mat15 P_{Mat15::Zero()};

  bool params_set_{false};
  bool initialized_{false};
  float t_last_{0.0f};
};

} // namespace rio