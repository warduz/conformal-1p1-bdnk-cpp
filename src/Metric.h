#ifndef SPACETIME_H
#define SPACETIME_H

#include "preamble.h"

class Metric {
  private:
    double g[2][2];
  
  public:
    Metric();
    double operator()(int a, int b) const;
    double inv(int a, int b) const;
    void raise_index(FourVector& vmu, const FourVector& v_mu) const;
    void lower_index(const FourVector& vmu, FourVector& v_mu) const;
};

#endif