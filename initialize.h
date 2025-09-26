#ifndef INITIALIZE_H
#define INITIALIZE_H

#include "Config.h"
#include "Metric.h"
#include "kt.h"

namespace Initialize {
  State get_q0(const Config& config, const Metric& g);

  void make_initial_gaussian_state(State& q, const Config& config, const Metric& g);
  void make_initial_ep_shock_state(State& q, const Config& config, const Metric& g);
  void make_initial_v_shock_state(State& q, const Config& config, const Metric& g);
  void make_initial_v_gaussian_state(State& q, const Config& config, const Metric& g);
  void make_initial_ep_v_gaussian_state(State& q, const Config& config, const Metric& g);

  Grid get_fermi_dirac(const Grid& x, double x0, double fL, double fR, double Delta);
  Grid get_gaussian(const Grid& x, double x0, double fmax, double fmin, double w);
  Grid get_negative_gaussian(const Grid& x, double x0, double f0, double fscale, double w);
  Grid get_temp(const Grid& ep, double ep_coeff);
  Grid get_shear_scalar(const Grid& T, double shear_coeff);
  Grid get_gaussian_dTdx(double A_ep, const Grid& ep, const Grid& x, double Delta, double ep_coeff);
  Grid get_gamma(const Grid& v);
};

#endif