#ifndef CURVE_CONVERSION_HPP
#define CURVE_CONVERSION_HPP

#include <array>

#include "bezier_curve.hpp"
#include "cubic_hermite_spline.hpp"
#include "hermite_splines.hpp"
#include "utils.hpp"

template <size_t dimensions>
std::array<Point<dimensions>, 4> convertHermitePointsToBezierPoints(const std::array<Point<dimensions>, 4>& p) {

  constexpr std::array<std::array<double, 4>, 4> conversion_matrix = {
    std::array<double, 4>{1, 0, 0, 0 },
    std::array<double, 4>{1, 0, 1./3, 0},
    std::array<double, 4>{0, 1, 0, -1./3},
    std::array<double, 4>{0, 1, 0, 0}
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

  constexpr std::array<std::array<double, 4>, 4> conversion_matrix = {
    std::array<double, 4>{1, 0, 0, 0 },
    std::array<double, 4>{0, 0, 0, 1},
    std::array<double, 4>{-3, 3, 0, 0},
    std::array<double, 4>{0, 0, -3, 3}
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

template<size_t degree, size_t n_points, size_t dimensions, size_t N>
void create_coversion_matrix(
  std::array<std::array<int, degree+1>, degree+1>& conversion_matrix,
  const BezierCurve<degree, n_points, dimensions>& b
)
{
// std::cout << N << '\n';
 if constexpr (N <= degree/2 && N >= 0) {
    std::array<std::array<int, degree+1>, degree+1-N> temp = b.template get_nth_derivative_to_current_weights_relation<N>();
    // debug_print(temp);
    conversion_matrix[2*N] = temp[0];
    conversion_matrix[2*N+1] = temp[temp.size()-1];
    if (N != 0) {
      for (auto& ele : conversion_matrix[2*N]) {
        ele *= (degree-N + 1);
      }
      for (auto& ele : conversion_matrix[2*N+1]) {
        ele *= (degree-N + 1);
      }
    }
    // debug_print(conversion_matrix);
    create_coversion_matrix<degree, n_points, dimensions, N+1>(conversion_matrix, b);
  } else if constexpr (N%2 == 0 && N == degree/2) {
    auto temp = b.template get_nth_derivative_to_current_weights_relation<N>();
    conversion_matrix[N] = temp[0];
    for (auto& ele : conversion_matrix[N]) {
      ele *= (degree-N);
    }
    // debug_print(conversion_matrix);
  }
  // debug_print(conversion_matrix);
  return;
}

template<size_t degree, size_t n_points, size_t dimensions>
HermiteSplines<degree, n_points, dimensions> B2H(BezierCurve<degree, n_points, dimensions>& b) {
  auto weights = b.get_weights();
  std::array<std::array<int, degree+1>, degree+1> conversion_matrix;
  for(auto& row : conversion_matrix) {
    for (auto& ele : row) {
      ele = 0;
    }
  }
  constexpr size_t N = 0;
  create_coversion_matrix<degree, n_points, dimensions, N>(conversion_matrix, b);
  // debug_print(conversion_matrix);
  std::array<Point<dimensions>, degree+1> p;
  for (int row = 0; row < conversion_matrix.size(); row++) {
    for (int dimension = 0; dimension < dimensions; dimension++) {
      p[row][dimension] = 0;
      for (int col = 0; col < conversion_matrix.size(); col++) {
        p[row][dimension] += conversion_matrix[row][col] * weights[col][dimension];
      }
    }
  }
  return HermiteSplines<degree, n_points, dimensions>(p);
}
#endif  // CURVE_CONVERSION_HPP
