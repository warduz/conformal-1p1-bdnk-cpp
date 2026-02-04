#ifndef EVOLVER_H
#define EVOLVER_H

#include "preamble.h"
#include "Config.h"
#include "initialize.h"

struct Evolver {
  Evolver(const String& config_file_name);

  Config config;
  Metric g {};

  std::filesystem::path directory;
  const Grid& x;
  const Grid& ts;

  int Nx;
  double hx;
  int Nt {};
  int Niter {};
  double ht;
  double t0;
  double tf;
  double t;

  IGrid Km2;
  IGrid Km1;
  IGrid K;
  IGrid Kp1;
  IGrid Kp2;

  State q;
  State qold;

  State qpKp1 {};
  State qmKp1 {};
  State qpKm1 {};
  State qmKm1 {};

  State HL {};
  State HR {};
  State flux_M {};
  State flux_P {};
  State Zq {};
  State src {};

  void run();
  void update_src(int id = -1);
  void do_rk_stage(int stage, int id = -1);
  void shift_q(int id = -1);
  void calculate_left_and_right_fluxes(int id = -1);
  // void calculate_residue(int id = -1);
  // void calculate_residue_and_update_q(int stage, int id = -1);
  void update_q(int stage, int id = -1);
  void set_boundary_condition_on_q(int id = -1);

  void Hflux(State& H, const State& qM, const State qP, int id = -1);
  State get_flux(const State& q, int id = -1);
  void make_flux(State& f, const State& q, int id = -1);
  State get_src(const State& q, int id = -1);
  void make_src(State& src, const State& q, int id = -1);

  FourVector get_C(const State& q, int id = -1) const;
  FourVectorConstRef get_C_ref(const State& q, int id = -1) const;
  FourVector get_u(const State& q, const Grid& T, int id = -1) const;
  Grid get_temp(const State& q, int id = -1) const;
  Grid get_Dmn(int m, int n, const FourVector& umu, int id = -1) const;
  Grid get_shear_scalar(const Grid& T) const;
  Grid get_ep(const Grid& T) const;
  Grid get_P(const Grid& ep) const;
  Grid get_P(double ep_coeff, const Grid& T) const;
  Grid get_v(const FourVector& umu) const;
  Grid get_tau_ep(const Grid& T, double a1, double etaovs) const;
  Grid get_tau_Q(const Grid& T, double a2, double etaovs) const;

  Grid get_Tmn(int m, int n, const State& q, const Grid& T, const FourVector& umu, const Grid& ep, int id = -1);
  Grid get_T11(const State& q, const Grid& T, const FourVector& umu, const Grid& ep, int id = -1);
  Grid get_Tmn_ideal(int m, int n, const State& q, const Grid& T, const FourVector& umu, int id = -1) const;
  Grid get_Tmn_ideal(int m, int n, const State& q, const Grid& T, const FourVector& umu, const Grid& ep, int id = -1) const;

  Grid get_X_mn(int m, int n, const State& q, const Grid& T, const FourVector& umu, const Grid& ep, int id = -1);
  Grid get_X_0n(int n, const State& q, const Grid& T, const FourVector& umu, const Grid& ep, int id = -1);

  Grid get_Ymnab(int m, int n, int a, int b, const State& q, const Grid& T, const FourVector& umu, int id = -1) const;
  Grid get_Ymnab(int m, int n, int a, int b, const State& q, const Grid& T, const FourVector& umu, const Grid& ep, const Grid& P, const Grid& tau_ep, const Grid& tau_P, const Grid& tau_Q, const Grid& shear_scalar, int id = -1) const;
  Grid get_Ymnab_inv(int m, int n, int a, int b, const State& q, const Grid& T, const FourVector& umu, const std::array<Grid, 4>& Y0n0b, int id = -1) const;
  Grid get_Hmnab(int m, int n, int a, int b, const State& q, const Grid& T, const FourVector& umu, int id = -1) const;
  Grid get_Hmnab(
    int m, int n, int a, int b, const Grid& T, const FourVector& umu,
    const Grid& ep, const Grid& P, const Grid& tau_ep, const Grid& tau_P, const Grid& tau_Q, const Grid& shear_scalar, int id=-1
  ) const;

  Grid get_Kn_u(double etaovs, const Grid& T, const Grid& v, const FourVector& umu, const Grid& X_x0, const Grid& X_xx, int id = -1) const;
  Grid get_Kn_T(double etaovs, const Grid& T, const Grid& v, const FourVector& umu, const Grid& X_x0, const Grid& X_xx, int id = -1) const;
  Grid get_Anorm(const Grid& T, const FourVector& umu, const Grid& X_00, const Grid& X_0x, const Grid& X_x0, const Grid& X_xx, const Grid& ep, const Grid& P, int id = -1) const;
  Grid get_Qnorm_sqrt(const Grid& T, const FourVector& umu, const Grid& X_00, const Grid& X_0x, const Grid& X_x0, const Grid& X_xx, const Grid& ep, const Grid& P, int id = -1) const;
};

#endif