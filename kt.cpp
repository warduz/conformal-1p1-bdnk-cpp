#include "kt.h"

namespace KT {
  Grid minmod(const Grid& a, const Grid& b) {
    return ((a.sign() + b.sign()) * (a.abs()).min(b.abs())) * 0.5;
  }


  Grid phi_minmod3(const Grid& r, double theta) {
    return minmod3(theta*r, 0.5*(1.0+r), theta*Grid::Ones(r.size()));
  }


  Grid minmod3(const Grid& a, const Grid& b, const Grid& c) {
    return minmod(a, minmod(b, c));
  }


  Grid qx(const Grid& qL, const Grid& qM, const Grid& qR, double theta_num, double hx) {
    double tol { 1e-12 };

    Grid dqL { (qM - qL) / hx };
    Grid dqR { (qR - qM) / hx };
    Grid r { dqL / (dqR + tol) };
    
    return dqR * phi_minmod3(r, theta_num);
  }


  Grid deriv_x(const Grid& qL, const Grid& qM, const Grid& qR, double theta_num, double hx) {
    return minmod3(
      theta_num * (qM - qL) / hx,
      0.5       * (qR - qL) / hx,
      theta_num * (qR - qM) / hx
    );
  }


  Grid general_ddx(const Grid& f, double theta_num, double hx) {
    int original_size { (int)(f.size()) };
    
    Grid extended_f(original_size + 2);

    extended_f(0) = f(0);
    extended_f.segment(1, original_size) = f;
    extended_f(original_size + 1) = f(original_size - 1);

    int size { original_size + 2 };

    IGrid Km1 { IGrid::LinSpaced(size-2, 0, size-2) };
    IGrid K   { IGrid::LinSpaced(size-2, 1, size-1) };
    IGrid Kp1 { IGrid::LinSpaced(size-2, 2, size-0) };

    return qx(extended_f(Km1), extended_f(K), extended_f(Kp1), theta_num, hx);
  }


  Grid qp(const Grid& qL, const Grid& qM, const Grid& qR, double hx, double theta_num, int id) {
    return qM - 0.5 * hx * deriv_x(qL, qM, qR, theta_num, hx);
  }


  Grid qm(const Grid& qL, const Grid& qM, const Grid& qR, double hx, double theta_num, int id) {
    return qM + 0.5 * hx * deriv_x(qL, qM, qR, theta_num, hx);
  }
};