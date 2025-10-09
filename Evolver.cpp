#include "Evolver.h"

Evolver::Evolver(const String& config_file_name)
  : config { config_file_name },
    Nx { config.Ncells }, hx { config.hx }, x(config.xs), ts(config.ts),
    ht { config.ht }, t0 { config.t0 }, tf { ts[ts.size()-1] }, t { t0 } {
      std::cout << "<------------------------------------------------------------------------------->\n";
      std::cout << "1D CONFORMAL BDNK SOLVER\n\n";

      std::cout << "CURRENT CONFIGURATION:\n\n";

      std::cout <<
        "Problem: " << config.problem_type << '\n' <<
        "(epL,epR) = (" << config.epL << "," << config.epR << ")" << '\n' <<
        "(vL,vR) = (" << config.vL << "," << config.vR << ")" << '\n' <<
        "a1 = " << config.a1 << '\n' << 
        "a2 = " << config.a2 << '\n' <<
        "c+ = " << config.cplus << '\n' << 
        "etaovs = " << config.etaovs << '\n' <<
        "Delta = " << config.Delta << '\n' <<
        "hx = " << config.hx << '\n' <<
        "ht = " << config.ht << '\n' <<
        "Ncells = " << config.Ncells << '\n' <<
        "Niter = " << config.Nt << '\n' <<
        "grid = [" << config.xL << "," << config.xR << "]";

      std::cout << "\n\n";

      // std::cout << "START OF SOLVER\n";

      Km2 = IGrid::LinSpaced(Nx-4, 0, Nx-4);
      Km1 = IGrid::LinSpaced(Nx-4, 1, Nx-3);
      K   = IGrid::LinSpaced(Nx-4, 2, Nx-2);
      Kp1 = IGrid::LinSpaced(Nx-4, 3, Nx-1);
      Kp2 = IGrid::LinSpaced(Nx-4, 4, Nx-0);

      q = Initialize::get_q0(config, g);
      qold = q;

      directory = config.save_dir;

      if (!std::filesystem::exists(directory)) 
        std::filesystem::create_directories(directory);
      
      // Initialize src, flux_M, flux_P
      #pragma omp parallel for
      for (int i=0; i<src.size(); ++i) {
        src[i].setZero(Nx);
        flux_M[i].setZero(Nx-4);
        flux_P[i].setZero(Nx-4);
      }

      run();

      std::cout << "END OF PROGRAM\n";      
      std::cout << "<------------------------------------------------------------------------------->\n";
    }


void Evolver::run() {
  /*****************************************************/
  auto start = std::chrono::high_resolution_clock::now();
  /*****************************************************/

  String unique {""};

  if (!config.save_in_unique_dir) unique = "_" + config.unique;

  std::ofstream t_file(config.save_dir + "/t" + unique + ".txt");
  std::ofstream all_ts_file(config.save_dir + "/all_ts" + unique + ".txt"); all_ts_file << ts;
  std::ofstream x_file(config.save_dir + "/x" + unique + ".txt"); x_file << x;

  std::ofstream T00_file(config.save_dir + "/T00" + unique + ".txt");
  std::ofstream T01_file(config.save_dir + "/T0x" + unique + ".txt");
  std::ofstream C_0_file(config.save_dir + "/C_0" + unique + ".txt");
  std::ofstream C_1_file(config.save_dir + "/C_1" + unique + ".txt");
  std::ofstream X_xx_file(config.save_dir + "/X_xx" + unique + ".txt");
  std::ofstream X_x0_file(config.save_dir + "/X_x0" + unique + ".txt");
  std::ofstream Txx_file(config.save_dir + "/Txx" + unique + ".txt");
  std::ofstream X_00_file(config.save_dir + "/X_00" + unique + ".txt");
  std::ofstream X_0x_file(config.save_dir + "/X_0x" + unique + ".txt");
  std::ofstream T_file(config.save_dir + "/T" + unique + ".txt");
  std::ofstream depdx_file(config.save_dir + "/depdx" + unique + ".txt");
  std::ofstream ep_file(config.save_dir + "/ep" + unique + ".txt");
  std::ofstream v_file(config.save_dir + "/v" + unique + ".txt");
  std::ofstream dvdx_file(config.save_dir + "/dvdx" + unique + ".txt");
  std::ofstream P_file(config.save_dir + "/P" + unique + ".txt");
  std::ofstream Kn_u_file(config.save_dir + "/Kn_u" + unique + ".txt");
  std::ofstream Kn_T_file(config.save_dir + "/Kn_T" + unique + ".txt");
  std::ofstream Anorm_file(config.save_dir + "/Anorm" + unique + ".txt");
  std::ofstream Qnorm_sqrt_file(config.save_dir + "/Qnorm_sqrt" + unique + ".txt");

  for (int nt=0; nt<=ts.size(); ++nt) {
    /*************************************************************/
    auto start_in_loop = std::chrono::high_resolution_clock::now();
    /*************************************************************/

    if (nt == ts.size()) break;

    if (nt % config.tmod == 0) {
      // std::cout << "SAVING INITIAL DATA. i = " << nt << "\n";
      const Grid curr_T { get_temp(q) };
      const Grid curr_ep { get_ep(curr_T) };
      const Grid curr_depdx { KT::general_ddx(curr_ep, config.theta_kt, config.hx) };
      const FourVector curr_umu { get_u(q, curr_T) };
      const Grid curr_T11 { get_T11(q, curr_T, curr_umu, curr_ep) };
      const Grid curr_X_00 { get_X_0n(0, q, curr_T, curr_umu, curr_ep) };
      const Grid curr_X_0x { get_X_0n(1, q, curr_T, curr_umu, curr_ep) };
      const Grid curr_v { get_v(curr_umu) };
      const Grid curr_dvdx { KT::general_ddx(curr_v, config.theta_kt, config.hx) };
      const Grid curr_P { get_P(curr_ep) };
      const Grid curr_Kn_u { get_Kn_u(config.etaovs, curr_T, curr_v, curr_umu, q[5], q[4]) };
      const Grid curr_Kn_T { get_Kn_T(config.etaovs, curr_T, curr_v, curr_umu, q[5], q[4]) };
      const Grid curr_Anorm { get_Anorm(curr_T, curr_umu, curr_X_00, curr_X_0x, q[5], q[4], curr_ep, curr_P) };
      const Grid curr_Qnorm_sqrt { get_Qnorm_sqrt(curr_T, curr_umu, curr_X_00, curr_X_0x, q[5], q[4], curr_ep, curr_P) };
      
      for (int nx=0; nx<Nx; ++nx) {
        T00_file << q[0](nx) << '\t';
        T01_file << q[1](nx) << '\t';
        C_0_file << q[2](nx) << '\t';
        C_1_file << q[3](nx) << '\t';
        X_xx_file << q[4](nx) << '\t';
        X_x0_file << q[5](nx) << '\t';
        Txx_file << curr_T11(nx) << '\t';
        X_00_file << curr_X_00(nx) << '\t';
        X_0x_file << curr_X_0x(nx) << '\t';
        T_file << curr_ep(nx) << '\t';
        depdx_file << curr_depdx(nx) << '\t';
        ep_file << curr_ep(nx) << '\t';
        v_file << curr_v(nx) << '\t';
        dvdx_file << curr_dvdx(nx) << '\t';
        P_file << curr_P(nx) << '\t';
        Kn_u_file << curr_Kn_u(nx) << '\t';
        Kn_T_file << curr_Kn_T(nx) << '\t';
        Anorm_file << curr_Anorm(nx) << '\t';
        Qnorm_sqrt_file << curr_Qnorm_sqrt(nx) << '\t';
      }

      T00_file << '\n';
      T01_file << '\n';
      t_file << ts[nt] << '\n';
      C_0_file << '\n';
      C_1_file << '\n';
      X_xx_file << '\n';
      X_x0_file << '\n';
      Txx_file << '\n';
      X_00_file << '\n';
      X_0x_file << '\n';
      T_file << '\n';
      depdx_file << '\n';
      ep_file << '\n';
      v_file <<'\n';
      dvdx_file << '\n';
      P_file << '\n';
      Kn_u_file << '\n';
      Kn_T_file << '\n';
      Anorm_file << '\n';
      Qnorm_sqrt_file << '\n';
    }

    do_rk_stage(1);
    // std::cout << "Done with RK Stage 1\n";
    do_rk_stage(2);
    // std::cout << "Done with RK Stage 2\n";

    qold = q;

    /***********************************************************/
    auto end_in_loop = std::chrono::high_resolution_clock::now();
    auto duration_in_loop = std::chrono::duration_cast<std::chrono::microseconds>(end_in_loop - start_in_loop);
    if (nt % config.tmod == 0) printf("Iteration %d: %.4f seconds\n", nt, duration_in_loop.count()/std::pow(10, 6));
    /**************************************************************************************************************************************/
  }

  /***************************************************/
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  printf("Time taken to run program: %.6f seconds\n", duration.count()/std::pow(10, 6));
  /****************************************************************************/
}


void Evolver::do_rk_stage(int stage, int id) {
  update_src(id);
  // std::cout << "DONE WITH update_src\n\n";
  shift_q(id);
  // std::cout << "DONE WITH shift_q\n\n";
  calculate_left_and_right_fluxes(id);
  // std::cout << "DONE WITH calculate_left_and_right_fluxes\n\n";
  // calculate_residue(id);
  // // std::cout << "DONE WITH calculate_residue\n\n";
  // calculate_residue_and_update_q(stage, id);
  update_q(stage, id);
  // std::cout << "DONE WITH update_q\n\n";
  // set_boundary_condition_on_q(id);
  // std::cout << "DONE WITH set_boundary_condition_on_q\n\n";
}


void Evolver::update_src(int id) {
  make_src(src, q);
}


void Evolver::shift_q(int id) {
  #pragma omp parallel for
  for (int i=0; i<q.size(); ++i) {
    qpKp1[i] = KT::qp(q[i](K)  , q[i](Kp1), q[i](Kp2), hx, config.theta_kt);
    qmKp1[i] = KT::qm(q[i](Km1), q[i](K)  , q[i](Kp1), hx, config.theta_kt);
    qpKm1[i] = KT::qp(q[i](Km1), q[i](K)  , q[i](Kp1), hx, config.theta_kt);
    qmKm1[i] = KT::qm(q[i](Km2), q[i](Km1), q[i](K)  , hx, config.theta_kt);
  }
}


void Evolver::calculate_left_and_right_fluxes(int id) {
  Hflux(HL, qmKm1, qpKm1, id); // Km1
  Hflux(HR, qmKp1, qpKp1, id); // Kp1
}


// void Evolver::calculate_residue(int id) {
//   // State src { get_src(q) };
//   #pragma omp parallel for
//   for (int i=0; i<6; ++i) {
//     Zq[i] = - (HR[i] - HL[i]) / hx + src[i](K);
//   }
// }


void Evolver::update_q(int stage, int id) {
  // #pragma omp parallel for
  for (int i=0; i<q.size(); ++i) {
    Zq[i] = - (HR[i] - HL[i]) / hx + src[i](K);

    switch (stage) {
      case 1: q[i](K) = q[i](K) + ht * Zq[i]; break;
      case 2: q[i](K) = 0.5*q[i](K) + 0.5*qold[i](K) + 0.5*ht*Zq[i]; break;
      default:
        throw std::runtime_error("Invalid RK stage.");
    }

    q[i](Nx-2) = q[i](Nx-3);
    q[i](Nx-1) = q[i](Nx-2);  
    q[i](1) = q[i](2);
    q[i](0) = q[i](1);
  }
}


// void Evolver::set_boundary_condition_on_q(int id) {
//   for (int i=0; i<q.size(); ++i) {
//     q[i](Nx-2) = q[i](Nx-3);
//     q[i](Nx-1) = q[i](Nx-2);  
//     q[i](1) = q[i](2);
//     q[i](0) = q[i](1);
//   }
// }


void Evolver::Hflux(State& H, const State& qM, const State qP, int id) {
  const double a { 1. }; // make a=1 always, no matter which cplus. Why??

  make_flux(flux_M, qM);
  make_flux(flux_P, qP);

  #pragma omp parallel for
  for (int i=0; i<H.size(); ++i) {
    H[i] = 0.5 * (flux_M[i] + flux_P[i]) - 0.5 * a * (qP[i] - qM[i]);
  }
}


State Evolver::get_flux(const State& q, int id) {
  const int current_size { static_cast<int>(q[0].size()) };
  State f {};

  f[2].setZero(current_size); // Zero
  f[3].setZero(current_size); // Zero
  make_flux(f, q, id);

  return f;
}


void Evolver::make_flux(State& f, const State& q, int id) {
  const Grid T { get_temp(q) };
  const FourVector umu { get_u(q, T) };
  const Grid ep { get_ep(T) };

  f[0] = q[1];
  f[1] = get_Tmn(1, 1, q, T, umu, ep);
  f[4] = - get_X_mn(0, 1, q, T, umu, ep);
  f[5] = - get_X_mn(0, 0, q, T, umu, ep);
}


State Evolver::get_src(const State& q, int id) {
  const int current_size { static_cast<int>(q[0].size()) };
  State src {};

  src[0].setZero(current_size); // Zero
  src[1].setZero(current_size); // Zero
  make_src(src, q, id);
  src[4].setZero(current_size); // Zero
  src[5].setZero(current_size); // Zero

  return src;
}


void Evolver::make_src(State& src, const State& q, int id) {
  const Grid T { get_temp(q) };
  const FourVector umu { get_u(q, T) };
  const Grid ep { get_ep(T) };

  src[2] = get_X_mn(0, 0, q, T, umu, ep, id);
  src[3] = get_X_mn(0, 1, q, T, umu, ep);
}


FourVector Evolver::get_C(const State& q, int id) const {
  return { q[2], q[3] };
}


FourVectorConstRef Evolver::get_C_ref(const State& q, int id) const {
  return { std::cref(q[2]), std::cref(q[3]) };
}


FourVector Evolver::get_u(const State& q, const Grid& T, int id) const {
  FourVectorConstRef C_mu { get_C_ref(q) };
  return { g.inv(0, 0) * C_mu[0].get() / T, g.inv(1, 1) * C_mu[1].get() / T };
}


Grid Evolver::get_temp(const State& q, int id) const {
  Grid CmuSquare { Grid::Zero(q[0].size()) };
  FourVectorConstRef C_mu { get_C_ref(q) };
  for (int i=0; i<C_mu.size(); ++i) {
    CmuSquare += g.inv(i, i) * C_mu[i].get() * C_mu[i].get();
  }
  return CmuSquare.abs().sqrt();
}


Grid Evolver::get_Dmn(int m, int n, const FourVector& umu, int id) const {
  return umu[m] * umu[n] + g(m, n);
}


Grid Evolver::get_shear_scalar(const Grid& T) const {
  return config.shear_coeff * (T*T*T);
}


Grid Evolver::get_ep(const Grid& T) const {
  return config.ep_coeff * (T*T*T*T);
}


Grid Evolver::get_P(const Grid& ep) const {
  return ep / 3.0;
}


Grid Evolver::get_P(double ep_coeff, const Grid& T) const {
  return ep_coeff * (T*T*T*T) / 3.0;
}


Grid Evolver::get_v(const FourVector& umu) const {
  return umu[1] / umu[0];
}


Grid Evolver::get_tau_ep(const Grid& T, double a1, double etaovs) const {
  return a1 * etaovs / T;
}


Grid Evolver::get_tau_Q(const Grid& T, double a2, double etaovs) const {
  return a2 * etaovs / T;
}


Grid Evolver::get_Tmn(int m, int n, const State& q, const Grid& T, const FourVector& umu, const Grid& ep, int id) {
  switch (2*m + n) {
		case 0: return q[0];
		case 1: return q[1];
		case 2: return q[1];
		case 3: return get_T11(q, T, umu, ep, id);
		default: throw std::out_of_range("In Evolver::get_Tmn. Invalid m or n.");
	}
}


Grid Evolver::get_T11(const State& q, const Grid& T, const FourVector& umu, const Grid& ep, int id) {
  const Grid T11_ideal { get_Tmn_ideal(1, 1, q, T, umu, ep) };
  Grid Hmnab_01 { get_Hmnab(1, 1, 0, 1, q, T, umu) };

  std::array<Grid, 4> X_mn {
    get_X_mn(0, 0, q, T, umu, ep), get_X_mn(0, 1, q, T, umu, ep),
    get_X_mn(1, 0, q, T, umu, ep), get_X_mn(1, 1, q, T, umu, ep)
  };

  std::array<Grid, 4> Hmnab {
    get_Hmnab(1, 1, 0, 0, q, T, umu), Hmnab_01,
    Hmnab_01, get_Hmnab(1, 1, 1, 1, q, T, umu)
  };

  std::array<Grid, 2> Fmna { Grid::Zero(T.size()), Grid::Zero(T.size()) };

  for (auto a: {0, 1}) {
    for (auto r: {0, 1}) {
      for (auto s: {0, 1}) {
        Fmna[a] += 2. * Hmnab[2*a+r] * umu[s] * g.inv(r, s) / (T*T);
      }
    }
  }
  
  Grid T11 { Grid::Zero(T.size()) };

  T11 += T11_ideal;

  for (auto a: {0, 1}) {
    for (auto r: {0, 1}) {
      T11 += (Hmnab[2*a+r] * X_mn[2*a+r]) / (T*T) + Fmna[a] * umu[r] * X_mn[2*a+r];
    }
  }

  return T11;
}


Grid Evolver::get_Tmn_ideal(int m, int n, const State& q, const Grid& T, const FourVector& umu, int id) const {
  Grid ep { get_ep(T) };
  return ep*umu[m]*umu[n] + get_P(ep)*get_Dmn(m, n, umu);
}


Grid Evolver::get_Tmn_ideal(int m, int n, const State& q, const Grid& T, const FourVector& umu, const Grid& ep, int id) const {
  return ep*umu[m]*umu[n] + get_P(ep)*get_Dmn(m, n, umu);
}


Grid Evolver::get_X_mn(int m, int n, const State& q, const Grid& T, const FourVector& umu, const Grid& ep, int id) {
  switch (2*m + n) {
		case 0: return get_X_0n(0, q, T, umu, ep, id);
		case 1: return get_X_0n(1, q, T, umu, ep, id);
		case 2: return q[5];
		case 3: return q[4];
		default: throw std::out_of_range("In Evolver::get_X_mn. Invalid m or n.");
	}
}


Grid Evolver::get_X_0n(int n, const State& q, const Grid& T, const FourVector& umu, const Grid& ep, int id) {
  const std::array<GridConstRef, 2> T0n { q[0], q[1] };
  const std::array<GridConstRef, 2> X_1n { q[5], q[4] };
  const std::array<Grid, 2> T0n_ideal { get_Tmn_ideal(0, 0, q, T, umu, ep), get_Tmn_ideal(0, 1, q, T, umu, ep) };

  const Grid P { get_P(ep) };
  const Grid tau_ep { get_tau_ep(T, config.a1, config.etaovs) };
  const Grid tau_P { tau_ep / 3. };
  const Grid tau_Q { get_tau_Q(T, config.a2, config.etaovs) };
  const Grid shear_scalar { get_shear_scalar(T) };

  const std::array<Grid, 4> Y0n1b {
    get_Ymnab(0, 0, 1, 0, q, T, umu, ep, P, tau_ep, tau_P, tau_Q, shear_scalar), get_Ymnab(0, 0, 1, 1, q, T, umu, ep, P, tau_ep, tau_P, tau_Q, shear_scalar),
    get_Ymnab(0, 1, 1, 0, q, T, umu, ep, P, tau_ep, tau_P, tau_Q, shear_scalar), get_Ymnab(0, 1, 1, 1, q, T, umu, ep, P, tau_ep, tau_P, tau_Q, shear_scalar)
  };

  const std::array<Grid, 4> Y0n0b {
    get_Ymnab(0, 0, 0, 0, q, T, umu, ep, P, tau_ep, tau_P, tau_Q, shear_scalar), get_Ymnab(0, 0, 0, 1, q, T, umu, ep, P, tau_ep, tau_P, tau_Q, shear_scalar),
    get_Ymnab(0, 1, 0, 0, q, T, umu, ep, P, tau_ep, tau_P, tau_Q, shear_scalar), get_Ymnab(0, 1, 0, 1, q, T, umu, ep, P, tau_ep, tau_P, tau_Q, shear_scalar)
  };

  const std::array<Grid, 2> Y0n0b_inv {
    get_Ymnab_inv(0, n, 0, 0, q, T, umu, Y0n0b), get_Ymnab_inv(0, n, 0, 1, q, T, umu, Y0n0b)
  };


  Grid X_0n { Grid::Zero(T.size()) };

  for (auto i: {0, 1})
    X_0n += Y0n0b_inv[i] * (T*T * (T0n[i].get() - T0n_ideal[i]) - (Y0n1b[2*i+0]*X_1n[0].get() + Y0n1b[2*i+1]*X_1n[1].get()));

  return X_0n;
}


Grid Evolver::get_Ymnab(int m, int n, int a, int b, const State& q, const Grid& T, const FourVector& umu, int id) const {
  Grid Fmna { Grid::Zero(T.size()) };

  const std::array<Grid, 2> Hmnab {
    get_Hmnab(m, n, a, 0, q, T, umu), get_Hmnab(m, n, a, 1, q, T, umu)
  };

  for (auto r: {0, 1})
  for (auto l: {0, 1})
    Fmna += 2. * Hmnab[r] * umu[l] * g.inv(l, r) / (T*T);

  return Hmnab[b] + T*T * Fmna * umu[b];
}



Grid Evolver::get_Ymnab(int m, int n, int a, int b, const State& q, const Grid& T, const FourVector& umu,
  const Grid& ep, const Grid& P, const Grid& tau_ep, const Grid& tau_P, const Grid& tau_Q, const Grid& shear_scalar, int id) const {
  Grid Fmna { Grid::Zero(T.size()) };

  const std::array<Grid, 2> Hmnab {
    get_Hmnab(m, n, a, 0, T, umu, ep, P, tau_ep, tau_P, tau_Q, shear_scalar), get_Hmnab(m, n, a, 1, T, umu, ep, P, tau_ep, tau_P, tau_Q, shear_scalar)
  };

  for (auto r: {0, 1})
  for (auto l: {0, 1})
    Fmna += 2. * Hmnab[r] * umu[l] * g.inv(l, r) / (T*T);

  return Hmnab[b] + T*T * Fmna * umu[b];
}


Grid Evolver::get_Ymnab_inv(int m, int n, int a, int b, const State& q, const Grid& T, const FourVector& umu, const std::array<Grid, 4>& Y0n0b, int id) const {
  // (Y0101 * Y0000 - Y0100 * Y0001 )
  const Grid factor { 1. / (Y0n0b[2*1+1]*Y0n0b[2*0+0] - Y0n0b[2*1+0]*Y0n0b[2*0+1]) };

  switch (2*n+b) {
    case 0: return factor * Y0n0b[2*1+1];
    case 1: return - factor * Y0n0b[2*0+1];
    case 2: return - factor * Y0n0b[2*1+0];
    case 3: return factor * Y0n0b[2*0+0];
    default: return Grid::Zero(T.size());
  }
}


Grid Evolver::get_Hmnab(int m, int n, int a, int b, const State& q, const Grid& T, const FourVector& umu, int id) const {
  const Grid ep { get_ep(T) };
  const Grid P { get_P(ep) };
  const Grid tau_ep { get_tau_ep(T, config.a1, config.etaovs) };
  const Grid tau_P { tau_ep / 3. };
  const Grid tau_Q { get_tau_Q(T, config.a2, config.etaovs) };
  const double cs { config.cs };
  const double bulk_scalar { config.bulk_scalar };
  const Grid shear_scalar { get_shear_scalar(T) };

  const std::array<Grid, 4> Dmn {
    get_Dmn(0, 0, umu), get_Dmn(0, 1, umu),
    get_Dmn(1, 0, umu), get_Dmn(1, 1, umu)
  };

  return (
      tau_ep * (ep + P) * T / (cs*cs) * umu[m] * umu[n] * (umu[a] * umu[b] + cs*cs * Dmn[2*a+b])
    + tau_P  * (ep + P) * T / (cs*cs) * Dmn[2*m+n]      * (umu[a] * umu[b] + cs*cs * Dmn[2*a+b])
    - 2. * shear_scalar * T * ((Dmn[2*m+a]*Dmn[2*n+b] + Dmn[2*m+b]*Dmn[2*n+a])/2. - Dmn[2*a+b]*Dmn[2*m+n]/3. )
    + tau_Q  * (ep + P) * T * (umu[m]*umu[b]*Dmn[2*a+n] + umu[n]*umu[b]*Dmn[2*m+a] + umu[m]*umu[a]*Dmn[2*b+n] + umu[n]*umu[a]*Dmn[2*m+b])
  );
}


Grid Evolver::get_Hmnab(
    int m, int n, int a, int b, const Grid& T, const FourVector& umu,
    const Grid& ep, const Grid& P, const Grid& tau_ep, const Grid& tau_P, const Grid& tau_Q, const Grid& shear_scalar, int id
  ) const {
    const double cs { config.cs };
    const double bulk_scalar { config.bulk_scalar };

    const std::array<Grid, 4> Dmn {
      get_Dmn(0, 0, umu), get_Dmn(0, 1, umu),
      get_Dmn(1, 0, umu), get_Dmn(1, 1, umu)
    };

    return (
        tau_ep * (ep + P) * T / (cs*cs) * umu[m] * umu[n] * (umu[a] * umu[b] + cs*cs * Dmn[2*a+b])
      + tau_P  * (ep + P) * T / (cs*cs) * Dmn[2*m+n]      * (umu[a] * umu[b] + cs*cs * Dmn[2*a+b])
      - 2. * shear_scalar * T * ((Dmn[2*m+a]*Dmn[2*n+b] + Dmn[2*m+b]*Dmn[2*n+a])/2. - Dmn[2*a+b]*Dmn[2*m+n]/3. )
      + tau_Q  * (ep + P) * T * (umu[m]*umu[b]*Dmn[2*a+n] + umu[n]*umu[b]*Dmn[2*m+a] + umu[m]*umu[a]*Dmn[2*b+n] + umu[n]*umu[a]*Dmn[2*m+b])
    );
}


Grid Evolver::get_Kn_u(double etaovs, const Grid& T, const Grid& v, const FourVector& umu, const Grid& X_x0, const Grid& X_xx, int id) const {
  return umu[0]*umu[0] * ( etaovs * X_xx / (T*T) - v * etaovs * X_x0 / (T*T));
}


Grid Evolver::get_Kn_T(double etaovs, const Grid& T, const Grid& v, const FourVector& umu, const Grid& X_x0, const Grid& X_xx, int id) const {
  return umu[0] * ( etaovs * X_x0 / (T*T) - v * etaovs * X_xx / (T*T));
}


Grid Evolver::get_Anorm(const Grid& T, const FourVector& umu, const Grid& X_00, const Grid& X_0x, const Grid& X_x0, const Grid& X_xx, const Grid& ep, const Grid& P, int id) const {
  Grid A { Grid::Zero(T.size()) };

  const std::array<GridConstRef, 4> X_mn_ref { X_00, X_0x, X_x0, X_xx };
  const Grid tau_ep { get_tau_ep(T, config.a1, config.etaovs) };

  for (auto j: { 0, 1 }) {
    for (auto k: { 0, 1 }) {
      A += 1./T * 4. * ep * tau_ep * (-umu[j]*umu[k] + get_Dmn(j, k, umu)/3.) * X_mn_ref[2*j+k].get();
    }
  }

  return A/(ep+P);
}


Grid Evolver::get_Qnorm_sqrt(const Grid& T, const FourVector& umu, const Grid& X_00, const Grid& X_0x, const Grid& X_x0, const Grid& X_xx, const Grid& ep, const Grid& P, int id) const {
  Grid Q0 { Grid::Zero(T.size()) };
  Grid Q1 { Grid::Zero(T.size()) };

  const std::array<GridConstRef, 4> X_mn_ref { X_00, X_0x, X_x0, X_xx };
  const Grid tau_Q { get_tau_Q(T, config.a2, config.etaovs) };

  for (auto j: { 0, 1 }) {
    for (auto k: { 0, 1 }) {
      Q0 += 1./T * tau_Q * (ep + P) * (umu[j] * get_Dmn(0, k, umu) - umu[k] * get_Dmn(0, j, umu) ) * X_mn_ref[2*j+k].get();
      Q1 += 1./T * tau_Q * (ep + P) * (umu[j] * get_Dmn(1, k, umu) - umu[k] * get_Dmn(1, j, umu) ) * X_mn_ref[2*j+k].get();
    }
  }

  return ( (Q0*Q0*g(0,0) + Q1*Q1*g(1,1)) / (ep + P).square() ).sqrt();
}