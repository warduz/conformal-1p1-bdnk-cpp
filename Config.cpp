#include "Config.h"

Config::Config(const String& args)  {
  // std::cout << "INSIDE CONFIG\n";
  try {
    auto tbl { toml::parse_file(args) };

    String _save_dir { tbl["params"]["save_dir"].value_or("") };
    bool _save_in_unique_dir { tbl["params"]["save_in_unique_dir"].value_or(true) };
    int _problem { tbl["params"]["problem"].value_or(0) };
    bool _x0 { tbl["params"]["x0"].value_or(true) };
    double _cCFL { tbl["params"]["cCFL"].value_or(0.5) };
    double _hx { tbl["params"]["hx"].value_or(0.1) };
    double _ht { tbl["params"]["ht"].value_or(0.5) };
    bool _choose_ht { tbl["params"]["choose_ht"].value_or(false) };
    double _Nx { tbl["params"]["Nx"].value_or(75.) };
    double _theta_kt { tbl["params"]["theta_kt"].value_or(1.0) };
    double _Neta { tbl["params"]["Neta"].value_or(1.0) };
    double _ep_coeff { tbl["params"]["ep_coeff"].value_or(10.0) };
    bool _choose_cplus_a1 { tbl["params"]["choose_cplus_a1"].value_or(false) };
    int _frame { tbl["params"]["frame"].value_or(1) };
    double _epL { tbl["params"]["epL"].value_or(0.48) };
    double _epR { tbl["params"]["epR"].value_or(0.12) };
    double _vL { tbl["params"]["vL"].value_or(0.) };
    double _vR { tbl["params"]["vR"].value_or(0.) };
    double _Delta { tbl["params"]["Delta"].value_or(1.0) };

    // std::cout <<
    //   _save_in_unique_dir << '\n' <<
    //   _problem << '\n' <<
    //   _x0 << '\n' <<
    //   _cCFL << '\n' <<
    //   _hx << '\n' <<
    //   _ht << '\n' <<
    //   _choose_ht << '\n' <<
    //   _Nx << '\n' <<
    //   _theta_kt << '\n' <<
    //   _Neta << '\n' <<
    //   _ep_coeff << '\n' <<
    //   _choose_cplus_a1 << '\n' <<
    //   _frame << '\n';

    int _Nt { static_cast<int>(_Nx/_ht) };
    int _tmod { static_cast<int>(1./_ht) };

    set_problem(_problem);
    set_numerical_params(_x0, _cCFL, _hx, _ht, _choose_ht, _Nx, _Nt, _tmod, _theta_kt);
    set_hydro_params(_Neta, _ep_coeff);
    set_initial_condition(_problem, _epL, _epR, _vL, _vR, _Delta);
    set_spacetime();
    set_save_dir(_save_dir, _save_in_unique_dir);
  } catch (const toml::parse_error& err) {
    std::cerr << "Error loading config file: " << err.description() << "\n";
    exit(1);
  }
}

void Config::set_problem(int i) {
  if (i < 0 || i > 4) throw std::invalid_argument("Invalid problem while setting config. Expected 0 <= i <= 4. Got i = " + std::to_string(i));

  std::array<String, 5> problem_list { "gaussian", "ep-shock", "v-shock", "v-gaussian", "ep-v-gaussian" };

  problem = i;
  problem_type = problem_list[i];
}

void Config::set_numerical_params(double _x0, double _cCFL, double _hx, double _ht, bool _choose_ht, double _Nx, double _Nt, int _tmod, double _theta_kt) {
  if (_ht > _cCFL*_hx) throw std::invalid_argument("CFL Condition not met. Expected ht <= cCFL*hx.");
  std::cout << "Nx = " << _Nx << '\n';
  if (_Nx <= 0) throw std::invalid_argument("Invalid Nx. Expected Nx > 0. Got Nx = " + std::to_string(_Nx));
  if (_theta_kt < 1.0 || _theta_kt > 2.0 ) throw std::invalid_argument("Invalid theta_kt. Expected 1.0 <= theta_kt <= 2.0. Got theta_kt = " + Util::round(_theta_kt));
  
  x0 = _x0;
  cCFL = _cCFL;
  hx = _hx;
  
  if   (_choose_ht) ht = _ht;
  else              ht = _cCFL*hx;

  Nx = _Nx;
  Nt = _Nt;
  tmod = _tmod;
  theta_kt = _theta_kt;
}

void Config::set_hydro_params(double _Neta, double _ep_coeff) {
  if (Neta <= 0) throw std::invalid_argument("Invalid Neta. Expected Neta > 0. Got Neta = " + std::to_string(Neta));
  if (ep_coeff <= 0) throw std::invalid_argument("Invalid ep_coeff. Expected ep_coeff > 0. Got ep_coeff = " + std::to_string(Neta));

  Neta = _Neta;
  ep_coeff = _ep_coeff;

  bulk_scalar = 0.;
  etaovs = Neta * 1./(4. * M_PI);
  shear_coeff =  4./3. * pow(ep_coeff, 0.5) * etaovs;
  cs = 1.0 / std::sqrt(3.0);
}

void Config::set_initial_condition(int problem, double _epL, double _epR, double _vL, double _vR, double _Delta) {
  epL = _epL;
  epR = _epR;
  vL = _vL;
  vR = _vR;
  Delta = _Delta;
}

void Config::set_spacetime() {
  t0 = 0.;
  ts = t0 + ht * Eigen::ArrayXd::LinSpaced(Nt+1, 0, Nt);

  int N { static_cast<int>(std::floor(Nx/hx)) };
  Grid xpos { Grid::LinSpaced(N, hx, N * hx) };
  Grid xneg { -xpos.reverse() };

  if (x0) {
    xs = Grid(xpos.size() + xneg.size() + 1);
    xs << xneg, 0.0, xpos;
  } else {
    xs = Grid(xpos.size() + xneg.size());
    xs << xneg, xpos;
  }
  
  Ncells = static_cast<int>(xs.size());
  xL = xs[0];
  xR = xs[Ncells-1];
  xM = xs[Ncells/2];
}

void Config::set_save_dir(const String& _save_dir, bool _save_in_unique_dir) {
  unique = 
        "bdnk_"
      + problem_type
      + "_(epL,epR)=(" + Util::trim(epL) + "," + Util::trim(epR) + ")"
      + "_(vL,vR)=(" + Util::trim(vL) + "," + Util::trim(vR) + ")"
      + "_a1=" + Util::round(a1) + "_a2=" + Util::round(a2)
      + "_etaovs=" + Util::round(etaovs) + "_Delta=" + Util::trim(Delta)
      + "_hx=" + Util::trim(hx) + "_ht=" + Util::trim(ht)
      + "_Ncells=" + std::to_string(Ncells) + "_grid=[" + Util::round(xL) + "," + Util::round(xR) + "]"
      + "_c+=" + Util::round(cplus);

  save_in_unique_dir = _save_in_unique_dir;

  if (_save_in_unique_dir) {
    if (_save_dir == "" || _save_dir == "./") {
      save_dir = unique;
    } else {
      save_dir = _save_dir + "/" + unique;
    }
  } else {
    save_dir = _save_dir;
  }

  if (!std::filesystem::exists(save_dir)) {
    std::filesystem::create_directories(save_dir);
  }
}