#ifndef CUBIC_HERMITE_SPLINE_HPP
#define CUBIC_HERMITE_SPLINE_HPP

#include <array>
#include <iostream>

#include "utils.hpp"

template<size_t N>
using Point = std::array<double, N>;

template<size_t n_points, size_t dimensions>
class CubicHermiteSpline {
public:
  explicit CubicHermiteSpline(const std::array<Point<dimensions>, 4>& p) noexcept : p(p) {
    const double delta_u = 1./(n_points - 1);
    double u = 0;
    for (int i = 0; i <= n_points - 1; i++) {
      const auto u2 = u * u;
      const auto u3 = u2 * u;
      coefficients[i][0] = 2 * u3 - 3 * u2 + 1;
      coefficients[i][1] = -2 * u3 + 3 * u2;
      coefficients[i][2] = u3 - 2 * u2 + u;
      coefficients[i][3] = u3 - u2;
      for (int j = 0; j < dimensions; j++) {
        points[i][j] = 0;
        for (int k = 0; k < 4; k++) {
          points[i][j] += coefficients[i][k] * p[k][j];
        }
      }
      u += delta_u;
    }
  }

  void print() const {
    for (const auto& p : points) {
      for (int j = 0; j < dimensions; j++) {
        std::cout << p[j] << ", ";
      }
      std::cout << '\n';
    }
  }
private:
  std::array<std::array<double, 4>, n_points> coefficients;
  std::array<Point<dimensions>, n_points> points;
  std::array<Point<dimensions>, 4> p;
};

#endif  // CUBIC_HERMITE_SPLINE_HPP
