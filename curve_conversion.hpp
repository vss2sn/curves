#ifndef CURVE_CONVERSION_HPP
#define CURVE_CONVERSION_HPP

#include <array>

#include "bezier_curve.hpp"
#include "cubic_hermite_spline.hpp"
#include "utils.hpp"

template <size_t dimensions>
std::array<Point<dimensions>, 4> convertHermitePointsToBezierPoints(const std::array<Point<dimensions>, 4>& p) {

  constexpr std::array<std::array<float, 4>, 4> conversion_matrix = {
    std::array<float, 4>{1, 0, 0, 0 },
    std::array<float, 4>{1, 0, 1./3, 0},
    std::array<float, 4>{0, 1, 0, -1./3},
    std::array<float, 4>{0, 1, 0, 0}
  };

  std::array<Point<dimensions>, 4> ans;
  for (int row = 0; row < conversion_matrix.size(); row++) {
    for (int dimension = 0; dimension < dimensions; dimension++) {
      ans[row][dimension] = 0;
      for (int col = 0; col < conversion_matrix.size(); col++) {
        ans[row][dimension] += conversion_matrix[row][col] * p[col][dimension];
      }
    }
  }
  return ans;
}

template <size_t dimensions>
std::array<Point<dimensions>, 4> convertBezierPointsToHermitePoints(const std::array<Point<dimensions>, 4>& p) {

  constexpr std::array<std::array<float, 4>, 4> conversion_matrix = {
    std::array<float, 4>{1, 0, 0, 0 },
    std::array<float, 4>{0, 0, 0, 1},
    std::array<float, 4>{-3, 3, 0, 0},
    std::array<float, 4>{0, 0, -3, 3}
  };

  std::array<Point<dimensions>, 4> ans;
  for (int row = 0; row < conversion_matrix.size(); row++) {
    for (int dimension = 0; dimension < dimensions; dimension++) {
      ans[row][dimension] = 0;
      for (int col = 0; col < conversion_matrix.size(); col++) {
        ans[row][dimension] += conversion_matrix[row][col] * p[col][dimension];
      }
    }
  }
  return ans;
}

#endif  // CURVE_CONVERSION_HPP
