#ifndef UTIL_H
#define UTIL_H

#include "preamble.h"

namespace Util {
  inline String trim(double value) {
    std::ostringstream oss;
    oss << std::defaultfloat << std::setprecision(15) << value;
    String str { oss.str() };

    size_t dot_pos = str.find('.');
    if (dot_pos != std::string::npos) {
      while (!str.empty() && str.back() == '0') str.pop_back();
      if (!str.empty() && str.back() == '.') str.pop_back();
    }

    return str;
  }

  inline String round(double value, int decimals = 2) {
    double rounded = std::round(value * std::pow(10, decimals)) / std::pow(10, decimals);

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(decimals) << rounded;
    return oss.str();
}
};

#endif