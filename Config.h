#ifndef CONFIG_H
#define CONFIG_H

#include "preamble.h"
#include "util.h"
#include "toml.hpp"

class Config {
  public:
    Config(const String& args);

    String save_dir;
    bool save_in_unique_dir;
    String unique;

    int problem;
    String problem_type;

    bool x0;
    double cCFL;
    double hx;
    double ht;
    int Nx;
    int Nt;
    int tmod;
    double theta_kt;
    
    Grid xs;
    Grid ts;
    double xL;
    double xR;
    int Ncells;
    double xM;
    double t0;

    double Neta;
    double ep_coeff;
    double shear_coeff;
    double bulk_scalar;
    double etaovs;
    double cs;

    double cplus;
    double a1;
    double a2;

    double epL;
    double epR;
    double vL;
    double vR;
    double Delta;
  
  private:
    void set_problem(int i);
    void set_numerical_params(double _x0, double _cCFL, double _hx, double _ht, bool _choose_ht, double _Nx, double _Nt, int _tmod, double _theta_kt);
    void set_hydro_params(double _Neta, double _ep_coeff);
    void set_initial_condition(int problem, double _epL, double _epR, double _vL, double _vR, double _Delta);
    void set_spacetime();
    void set_hydro_frame(bool _choose_cplus_a1, int _frame);
    void set_save_dir(const String& _save_dir, bool _save_in_unique_dir);
};

#endif