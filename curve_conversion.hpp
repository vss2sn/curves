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
void create_Bezier_to_Hermite_conversion_matrix_impl(std::array<std::array<double, degree+1>, degree+1>& conversion_matrix) {
  // std::cout << N << '\n';
  if constexpr (N < (degree+1)/2 && N >= 0) {
    // std::cout << N << '\n';
    const auto temp = BezierCurve<degree, n_points, dimensions>::template get_nth_derivative_to_current_weights_relation<N>();
    conversion_matrix[2*N] = temp[0];
    conversion_matrix[2*N+1] = temp[temp.size()-1];
    if (N != 0) {
      size_t multiplicative_factor = 1;
      for (int i = degree; i > degree-N && degree-N > 0; i--) {
        multiplicative_factor *= i;
      }
      for (auto& ele : conversion_matrix[2*N]) {
        ele *= multiplicative_factor;
      }
      for (auto& ele : conversion_matrix[2*N+1]) {
        ele *= multiplicative_factor;
      }
    }
    create_Bezier_to_Hermite_conversion_matrix_impl<degree, n_points, dimensions, N+1>(conversion_matrix);
  } else if constexpr (degree%2 == 0 && N == degree/2) {
    // std::cout << N << '\n';
    const auto temp = BezierCurve<degree, n_points, dimensions>::template get_nth_derivative_to_current_weights_relation<N>();
    conversion_matrix[2 * N] = temp[0];
    // if (N != 0) {
      size_t multiplicative_factor = 1;
      for (int i = degree; i > degree-N && degree-N > 0; i--) {
        multiplicative_factor *= i;
      }
      for (auto& ele : conversion_matrix[2*N]) {
        ele *= multiplicative_factor;
      }
    // }
  }
  return;
}

template<size_t degree, size_t n_points, size_t dimensions>
std::array<std::array<double, degree+1>, degree+1> create_Bezier_to_Hermite_conversion_matrix() {
  std::array<std::array<double, degree+1>, degree+1> conversion_matrix;
  for(auto& row : conversion_matrix) {
    for (auto& ele : row) {
      ele = 0;
    }
  }
  create_Bezier_to_Hermite_conversion_matrix_impl<degree, n_points, dimensions, 0>(conversion_matrix);
  // debug_print(conversion_matrix);
  return conversion_matrix;
}

template<size_t degree, size_t n_points, size_t dimensions>
HermiteSplines<degree, n_points, dimensions> B2H(BezierCurve<degree, n_points, dimensions>& b) {
  const auto weights = b.get_weights();
  const auto conversion_matrix = create_Bezier_to_Hermite_conversion_matrix<degree, n_points, dimensions>();
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

template<size_t degree, size_t n_points, size_t dimensions>
BezierCurve<degree, n_points, dimensions> H2B(HermiteSplines<degree, n_points, dimensions>& hs) {
  const auto p = hs.get_p();
  const auto conversion_matrix = create_Bezier_to_Hermite_conversion_matrix<degree, n_points, dimensions>();
  const auto inverted_conversion_matrix = inverse_using_LU_decomp<degree+1>(conversion_matrix);
  std::array<Point<dimensions>, degree+1> weights;
  for (int row = 0; row < inverted_conversion_matrix.size(); row++) {
    for (int dimension = 0; dimension < dimensions; dimension++) {
      weights[row][dimension] = 0;
      for (int col = 0; col < inverted_conversion_matrix.size(); col++) {
        weights[row][dimension] += inverted_conversion_matrix[row][col] * p[col][dimension];
      }
    }
  }
  return BezierCurve<degree, n_points, dimensions>(weights);
}


#endif  // CURVE_CONVERSION_HPP
