#ifndef PREAMBLE_H
#define PREAMBLE_H

#include <array>
#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <omp.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

using Grid = Eigen::ArrayXd;
using Vector = std::vector<double>;
using IGrid = Eigen::ArrayXi;
using FourVector = std::array<Grid, 2>;
using State = std::array<Grid, 6>;
using String = std::string;

using GridRef = std::reference_wrapper<Grid>;
using GridConstRef = std::reference_wrapper<const Grid>;
using FourVectorRef = std::array<GridRef, 2>;
using FourVectorConstRef = std::array<GridConstRef, 2>;
using StateRef = std::array<GridRef, 6>;
using StateConstRef = std::array<GridConstRef, 6>;

#endif