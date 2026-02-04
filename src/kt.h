#ifndef KT_H
#define KT_H

#include "preamble.h"

namespace KT {
  Grid minmod(const Grid& a, const Grid& b);
  Grid minmod3(const Grid& a, const Grid& b, const Grid& c);
  Grid qx(const Grid& qL, const Grid& qM, const Grid& qR, double theta_num, double hx);
  Grid deriv_x(const Grid& qL, const Grid& qM, const Grid& qR, double theta_num, double hx);
  Grid general_ddx(const Grid& f, double theta_num, double hx);
  Grid qp(const Grid& qL, const Grid& qM, const Grid& qR, double hx, double theta_num, int id = 0);
  Grid qm(const Grid& qL, const Grid& qM, const Grid& qR, double hx, double theta_num, int id = 0);
};

#endif