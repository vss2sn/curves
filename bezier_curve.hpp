#ifndef BEZIER_CURVE_HPP
#define BEZIER_CURVE_HPP

#include <array>
#include <iostream>

#include "binomial.hpp"
#include "utils.hpp"

template <size_t degree, size_t n_points, size_t dimensions>
class BezierCurve {
public:
  constexpr explicit BezierCurve(const std::array<Point<dimensions>, degree +1>& weights) : weights(weights) {
    constexpr std::array<int, degree + 1> binomial_coeffs = find_all_binomial_coefficients<degree>();
    for (int i = 0; i < n_points; i++) {
      const auto [t_values, one_minus_t_values] = BinomialParamterValues<degree>(i * 1./(n_points - 1));
      for (int j = 0; j < degree + 1; j++) {
        coefficients[j] = binomial_coeffs[j] * t_values[j] * one_minus_t_values[j];
      }
      for (int k = 0; k < dimensions; k++) {
        points[i][k] = 0;
        for (int j = 0; j < degree + 1; j++) {
          points[i][k] += coefficients[j] * weights[j][k];
        }
      }
    }
  }

  constexpr void print() const {
    for (int i = 0; i < n_points; i++) {
      for (int j = 0; j < dimensions; j++) {
        std::cout <<  points[i][j] << ", ";
      }
      std::cout << '\n';
    }
  }

  constexpr BezierCurve<degree-1, n_points, dimensions> get_derivative() const {
    if constexpr (degree > 0) {
      std::array<Point<dimensions>, degree> derivative_weights;
      for (int i = 0; i < degree; i++) {
        for (int d = 0; d < dimensions; d++)
        derivative_weights[i][d] = degree * (weights[i+1][d] - weights[i][d]);
      }
      return BezierCurve<degree-1, n_points, dimensions>(derivative_weights);
    } else {
      return BezierCurve<0, n_points, dimensions>();
    }
  }

  template<size_t N>
  constexpr static std::array<std::array<double, degree+1>, degree+1-N> get_nth_derivative_to_current_weights_relation() {
    if constexpr (N >= 1) {
      // assert(degree+1-N+1 > 0);  // Debug assert
      std::array<std::array<double, degree+1-N+1>, degree+1-N> relation_matrix;
      for (int row = 0; row <degree+1-N; row++) {
        for (int col = 0; col < degree+1-N+1; col++) {
          if (row == col) {
            relation_matrix[row][col] = -1;
          }
          else if (row + 1 == col){
            relation_matrix[row][col] = 1;
          } else {
            relation_matrix[row][col] = 0;
          }
        }
      }
      return multiply_two_matrices(relation_matrix, get_nth_derivative_to_current_weights_relation<N-1>());
    } else {
      std::array<std::array<double, degree+1>, degree+1> relation_matrix;
      for (int row = 0; row <degree+1; row++) {
        for (int col = 0; col < degree+1; col++) {
          if (row  == col) {
            relation_matrix[row][col] = 1;
          } else {
            relation_matrix[row][col] = 0;
          }
        }
      }
      return relation_matrix;
    }
  }

  constexpr std::array<Point<dimensions>, degree + 1> get_weights() const {
    return weights;
  }
private:
  const std::array<Point<dimensions>, degree + 1> weights;
  std::array<Point<dimensions>, n_points> points;
  std::array<double, degree+1> coefficients;
};

#endif  // BEZIER_CURVE_HPP
