#include "initialize.h"

namespace Initialize {
  State get_q0(const Config& config, const Metric& g) {
    State q {};

    switch (config.problem) {
      case 0: make_initial_gaussian_state(q, config, g); break;
      case 1: make_initial_ep_shock_state(q, config, g); break;
      case 2: make_initial_v_shock_state(q, config, g); break;
      case 3: make_initial_v_gaussian_state(q, config, g); break;
      case 4: make_initial_ep_v_gaussian_state(q, config, g); break;
      default: throw std::invalid_argument("Invalid problem. Expected: 0, 1, 2, 3, 4, 5.");
    }

    return q;
  }


  void make_initial_gaussian_state(State& q, const Config& config, const Metric& g) {
    double epR { config.epR };
    double epL { config.epL };
    double vL { config.vL };
    double ep_coeff { config.ep_coeff };
    double shear_coeff { config.shear_coeff };
    double theta_kt { config.theta_kt };
    double hx { config.hx };
    double a2 { config.a2 };
    double etaovs { config.etaovs };
    double x0 { config.xM };
    int N { config.Ncells };
    double Delta { config.Delta };
    const Grid& x { config.xs };

    Grid zero { Grid::Zero(N) };

    Grid ep { get_gaussian(x, x0, epL, epR, Delta) };
    Grid T { get_temp(ep, ep_coeff) };

    Grid v { zero };
    Grid gamma { get_gamma(v) };
    FourVector umu { gamma, gamma*v };
    FourVector u_mu { zero, zero };
    g.lower_index(umu, u_mu);

    Grid dTdx { get_gaussian_dTdx(epL, ep, x, Delta, ep_coeff) };
    // Grid dTdx_num { KT::general_ddx(T, theta_kt, hx) };

    Grid dxu_t { KT::general_ddx(u_mu[0], theta_kt, hx) };
    Grid dxu_x { KT::general_ddx(u_mu[1], theta_kt, hx) };

    Grid C_t { T*u_mu[0] };
    Grid C_x { T*u_mu[1] };

    Grid X_xx { zero };
    Grid X_xt { - dTdx };

    Grid Ttt { ep };
    Grid Ttx { zero };

    q = { Ttt, Ttx, C_t, C_x, X_xx, X_xt };

    #pragma omp parallel for
    for (int i=0; i<q.size(); ++i) {
      const double right_ghost { q[i](N-3) };
      const double left_ghost { q[i](2) };
      q[i].coeffRef(N-2) = right_ghost;
      q[i].coeffRef(N-1) = right_ghost;
      q[i].coeffRef(1) = left_ghost;
      q[i].coeffRef(0) = left_ghost;
    }

  }


  void make_initial_ep_shock_state(State& q, const Config& config, const Metric& g) {
    double epR { config.epR };
    double epL { config.epL };
    double vL { config.vL };
    double vR { config.vR };
    double ep_coeff { config.ep_coeff };
    double shear_coeff { config.shear_coeff };
    double theta_kt { config.theta_kt };
    double hx { config.hx };
    double a2 { config.a2 };
    double etaovs { config.etaovs };
    double x0 { config.xM };
    int N { config.Ncells };
    double Delta { config.Delta };
    const Grid& x { config.xs };

    Grid zero { Grid::Zero(N) };
    Grid ones { Grid::Ones(N) };

    Grid ep { get_fermi_dirac(x, x0, epL, epR, Delta) };
    Grid T { get_temp(ep, ep_coeff) };

    Grid v { vL * ones };
    Grid gamma { get_gamma(v) };
    FourVector umu { gamma, gamma*v };
    FourVector u_mu { zero, zero };
    g.lower_index(umu, u_mu);

    Grid dxep { KT::general_ddx(ep, theta_kt, hx) };
    Grid dxu_t { KT::general_ddx(u_mu[0], theta_kt, hx) };
    Grid dxu_x { KT::general_ddx(u_mu[1], theta_kt, hx) };

    Grid dtu_x { 1./umu[0] * ( (u_mu[0]*u_mu[0] * dxu_x - umu[1]*umu[1] * dxu_x - umu[1] * dxep/(4.*ep))*u_mu[1]/(3.+ 2.*u_mu[1]*u_mu[1]) - umu[1]*dxu_x-dxep/(4.*ep) ) };
    Grid dtu_t { - umu[1] * dtu_x / umu[0] };

    Grid C_t { T*u_mu[0] };
    Grid C_x { T*u_mu[1] };

    Grid X_xx { KT::general_ddx(C_x, theta_kt, hx) };
    Grid X_xt { KT::general_ddx(C_t, theta_kt, hx) };

    Grid sigmatt { - 2./3.*(umu[0]*umu[0] - 1.) * dtu_t + u_mu[0] * u_mu[1] * dxu_t - (umu[0]*umu[0] - 1.) * dxu_x/3 };
    Grid sigmatx { u_mu[0]*u_mu[1]*dtu_t/6. + dtu_x * (umu[0]*umu[0] - 1.)/2. - dxu_x * u_mu[0] * u_mu[1] /6. - umu[0]*umu[0] * dxu_t/2. };

    Grid Ttt { ep*(4*(umu[0]*umu[0])/3. -1./3.) - 2. * get_shear_scalar(T, shear_coeff) * sigmatt };
    Grid Ttx { 4./3.*ep*umu[0]*umu[1] - 2. * get_shear_scalar(T, shear_coeff) * sigmatx };

    q = { Ttt, Ttx, C_t, C_x, X_xx, X_xt };

    #pragma omp parallel for
    for (int i=0; i<q.size(); ++i) {
      const double right_ghost { q[i](N-3) };
      const double left_ghost { q[i](2) };
      q[i].coeffRef(N-2) = right_ghost;
      q[i].coeffRef(N-1) = right_ghost;
      q[i].coeffRef(1) = left_ghost;
      q[i].coeffRef(0) = left_ghost;
    }


  }


  void make_initial_v_shock_state(State& q, const Config& config, const Metric& g) {
    double epR { config.epR };
    double epL { config.epL };
    double vL { config.vL };
    double vR { config.vR };
    double ep_coeff { config.ep_coeff };
    double shear_coeff { config.shear_coeff };
    double theta_kt { config.theta_kt };
    double hx { config.hx };
    double a2 { config.a2 };
    double etaovs { config.etaovs };
    double x0 { config.xM };
    int N { config.Ncells };
    double Delta { config.Delta };
    const Grid& x { config.xs };

    Grid zero { Grid::Zero(N) };
    Grid ones { Grid::Ones(N) };

    Grid ep { epL*ones };
    Grid T { get_temp(ep, ep_coeff) };

    Grid v { get_fermi_dirac(x, x0, vL, vR, Delta)};
    Grid gamma { get_gamma(v) };
    FourVector umu { gamma, gamma*v };
    FourVector u_mu { zero, zero };
    g.lower_index(umu, u_mu);

    Grid dxep { zero };
    Grid dxu_t { KT::general_ddx(u_mu[0], theta_kt, hx) };
    Grid dxu_x { KT::general_ddx(u_mu[1], theta_kt, hx) };

    Grid dtu_x { 1./umu[0] * ( (u_mu[0]*u_mu[0] * dxu_x - umu[1]*umu[1] * dxu_x - umu[1] * dxep/(4.*ep))*u_mu[1]/(3.+ 2.*u_mu[1]*u_mu[1]) - umu[1]*dxu_x-dxep/(4.*ep) ) };
    Grid dtu_t { - umu[1] * dtu_x / umu[0] };

    Grid C_t { T*u_mu[0] };
    Grid C_x { T*u_mu[1] };

    Grid X_xx { KT::general_ddx(C_x, theta_kt, hx) };
    Grid X_xt { KT::general_ddx(C_t, theta_kt, hx) };

    Grid sigmatt { - 2./3.*(umu[0]*umu[0] - 1.) * dtu_t + u_mu[0] * u_mu[1] * dxu_t - (umu[0]*umu[0] - 1.) * dxu_x/3 };
    Grid sigmatx { u_mu[0]*u_mu[1]*dtu_t/6. + dtu_x * (umu[0]*umu[0] - 1.)/2. - dxu_x * u_mu[0] * u_mu[1] /6. - umu[0]*umu[0] * dxu_t/2. };

    Grid Ttt { ep*(4*(umu[0]*umu[0])/3. -1./3.) - 2. * get_shear_scalar(T, shear_coeff) * sigmatt };
    Grid Ttx { 4./3.*ep*umu[0]*umu[1] - 2. * get_shear_scalar(T, shear_coeff) * sigmatx };

    q = { Ttt, Ttx, C_t, C_x, X_xx, X_xt };

    #pragma omp parallel for
    for (int i=0; i<q.size(); ++i) {
      const double right_ghost { q[i](N-3) };
      const double left_ghost { q[i](2) };
      q[i].coeffRef(N-2) = right_ghost;
      q[i].coeffRef(N-1) = right_ghost;
      q[i].coeffRef(1) = left_ghost;
      q[i].coeffRef(0) = left_ghost;
    }

  }

  
  void make_initial_v_gaussian_state(State& q, const Config& config, const Metric& g) {
    double epR { config.epR };
    double epL { config.epL };
    double vL { config.vL };
    double vR { config.vR };
    double ep_coeff { config.ep_coeff };
    double shear_coeff { config.shear_coeff };
    double theta_kt { config.theta_kt };
    double hx { config.hx };
    double a2 { config.a2 };
    double etaovs { config.etaovs };
    double x0 { config.xM };
    int N { config.Ncells };
    double Delta { config.Delta };
    const Grid& x { config.xs };

    Grid zero { Grid::Zero(N) };
    Grid ones { Grid::Ones(N) };

    Grid ep { epL * ones };
    Grid T { get_temp(ep, ep_coeff) };

    Grid v { get_negative_gaussian(x, x0, vL, vR, Delta) };
    Grid gamma { get_gamma(v) };
    FourVector umu { gamma, gamma*v };
    FourVector u_mu { zero, zero };
    g.lower_index(umu, u_mu);

    Grid dxep { zero };
    Grid dxu_t { KT::general_ddx(u_mu[0], theta_kt, hx) };
    Grid dxu_x { KT::general_ddx(u_mu[1], theta_kt, hx) };
    Grid dtu_x { 1./umu[0] * ( (u_mu[0]*u_mu[0] * dxu_x - umu[1]*umu[1] * dxu_x - umu[1] * dxep/(4.*ep))*u_mu[1]/(3.+ 2.*u_mu[1]*u_mu[1]) - umu[1]*dxu_x-dxep/(4.*ep) ) };
    Grid dtu_t { - umu[1] * dtu_x / umu[0] };

    Grid C_t { T*u_mu[0] };
    Grid C_x { T*u_mu[1] };

    Grid X_xx { KT::general_ddx(C_x, theta_kt, hx) };
    Grid X_xt { KT::general_ddx(C_t, theta_kt, hx) };

    Grid sigmatt { - 2./3.*(umu[0]*umu[0] - 1.) * dtu_t + u_mu[0] * u_mu[1] * dxu_t - (umu[0]*umu[0] - 1.) * dxu_x/3 };
    Grid sigmatx { u_mu[0]*u_mu[1]*dtu_t/6. + dtu_x * (umu[0]*umu[0] - 1.)/2. - dxu_x * u_mu[0] * u_mu[1] /6. - umu[0]*umu[0] * dxu_t/2. };

    Grid Ttt { ep*(4*(umu[0]*umu[0])/3. -1./3.) - 2. * get_shear_scalar(T, shear_coeff) * sigmatt };
    Grid Ttx { 4./3.*ep*umu[0]*umu[1] - 2. * get_shear_scalar(T, shear_coeff) * sigmatx };

    q = { Ttt, Ttx, C_t, C_x, X_xx, X_xt };

    #pragma omp parallel for
    for (int i=0; i<q.size(); ++i) {
      const double right_ghost { q[i](N-3) };
      const double left_ghost { q[i](2) };
      q[i].coeffRef(N-2) = right_ghost;
      q[i].coeffRef(N-1) = right_ghost;
      q[i].coeffRef(1) = left_ghost;
      q[i].coeffRef(0) = left_ghost;
    }


  }
  

  void make_initial_ep_v_gaussian_state(State& q, const Config& config, const Metric& g) {
    double epR { config.epR };
    double epL { config.epL };
    double vL { config.vL };
    double vR { config.vR };
    double ep_coeff { config.ep_coeff };
    double shear_coeff { config.shear_coeff };
    double theta_kt { config.theta_kt };
    double hx { config.hx };
    double a2 { config.a2 };
    double etaovs { config.etaovs };
    double x0 { config.xM };
    int N { config.Ncells };
    double Delta { config.Delta };
    const Grid& x { config.xs };

    Grid zero { Grid::Zero(N) };

    Grid ep { get_gaussian(x, x0, epL, epR, Delta) };
    Grid T { get_temp(ep, ep_coeff) };

    Grid v { get_gaussian(x, x0, vL, vR, Delta) };
    Grid gamma { get_gamma(v) };
    FourVector umu { gamma, gamma*v };
    FourVector u_mu { zero, zero };
    g.lower_index(umu, u_mu);

    // Grid dTdx { get_gaussian_dTdx(epL, ep, x, Delta, ep_coeff) };

    Grid dxep { KT::general_ddx(ep, theta_kt, hx) };
    Grid dxu_t { KT::general_ddx(u_mu[0], theta_kt, hx) };
    Grid dxu_x { KT::general_ddx(u_mu[1], theta_kt, hx) };
    Grid dtu_x { 1./umu[0] * ( (u_mu[0]*u_mu[0] * dxu_x - umu[1]*umu[1] * dxu_x - umu[1] * dxep/(4.*ep))*u_mu[1]/(3.+ 2.*u_mu[1]*u_mu[1]) - umu[1]*dxu_x-dxep/(4.*ep) ) };
    Grid dtu_t { - umu[1] * dtu_x / umu[0] };

    Grid C_t { T*u_mu[0] };
    Grid C_x { T*u_mu[1] };

    Grid X_xx { KT::general_ddx(C_x, theta_kt, hx) };
    Grid X_xt { KT::general_ddx(C_t, theta_kt, hx) };

    Grid sigmatt { - 2./3.*(umu[0]*umu[0] - 1.) * dtu_t + u_mu[0] * u_mu[1] * dxu_t - (umu[0]*umu[0] - 1.) * dxu_x/3 };
    Grid sigmatx { u_mu[0]*u_mu[1]*dtu_t/6. + dtu_x * (umu[0]*umu[0] - 1.)/2. - dxu_x * u_mu[0] * u_mu[1] /6. - umu[0]*umu[0] * dxu_t/2. };

    Grid Ttt { ep*(4*(umu[0]*umu[0])/3. -1./3.) - 2. * get_shear_scalar(T, shear_coeff) * sigmatt };
    Grid Ttx { 4./3.*ep*umu[0]*umu[1] - 2. * get_shear_scalar(T, shear_coeff) * sigmatx };

    q = { Ttt, Ttx, C_t, C_x, X_xx, X_xt };

    #pragma omp parallel for
    for (int i=0; i<q.size(); ++i) {
      const double right_ghost { q[i](N-3) };
      const double left_ghost { q[i](2) };
      q[i].coeffRef(N-2) = right_ghost;
      q[i].coeffRef(N-1) = right_ghost;
      q[i].coeffRef(1) = left_ghost;
      q[i].coeffRef(0) = left_ghost;
    }

  }


  Grid get_fermi_dirac(const Grid& x, double x0, double fL, double fR, double Delta) {
    return fR + (fL-fR)/(1. + ( (x-x0)/Delta ).exp() );
  }


  Grid get_gaussian(const Grid& x, double x0, double fmax, double fmin, double w) {
    return fmax * (-(x-x0).square()/(w*w)).exp() + fmin;
  }


  Grid get_negative_gaussian(const Grid& x, double x0, double f0, double fscale, double w) {
    return f0 - fscale * (-(x-x0).square()/(w*w)).exp();
  }


  Grid get_temp(const Grid& ep, double ep_coeff) {
    return (ep/ep_coeff).pow(1./4.);
  }


  Grid get_shear_scalar(const Grid& T, double shear_coeff) {
    return shear_coeff * (T*T*T);
  }


  Grid get_gaussian_dTdx(double A_ep, const Grid& ep, const Grid& x, double Delta, double ep_coeff) {
    return - ( A_ep *  x * (-x*x/(Delta*Delta)).exp() ) / (2. * (Delta*Delta) * ep_coeff * (ep/ep_coeff).pow(3./4.));
  }


  Grid get_gamma(const Grid& v) {
    return 1. / (1. - v.square()).sqrt();
  }
};