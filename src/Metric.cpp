#include "Metric.h"

Metric::Metric() {
  for (int a=0; a<2; ++a)
  for (int b=0; b<2; ++b)
    g[a][b] = (a == b) ? (a == 0 ? -1 : 1) : 0;
}


double Metric::operator()(int a, int b) const {
  return g[a][b];
}


double Metric::inv(int a, int b) const {
  return g[a][b];
}


void Metric::raise_index(FourVector& vmu, const FourVector& v_mu) const {
  for (int a=0; a<2; ++a)
  for (int b=0; b<2; ++b)
    vmu[a] += g[a][b] * v_mu[b];
}


void Metric::lower_index(const FourVector& vmu, FourVector& v_mu) const {
  for (int a=0; a<2; ++a)
  for (int b=0; b<2; ++b)
    v_mu[a] += inv(a, b) * vmu[b];
}