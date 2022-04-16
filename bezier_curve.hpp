#ifndef BEZIER_CURVE_HPP
#define BEZIER_CURVE_HPP

#include <cassert>

#include "binomial.hpp"
#include "utils.hpp"

template <size_t degree, size_t n_points, size_t dimensions>
class BezierCurve {
public:
  explicit BezierCurve(const std::array<Point<dimensions>, degree +1>& weights) noexcept : weights(weights) {
    std::array<int, degree + 1> binomial_coeffs = find_all_binomial_coefficients<degree>();
    for (int i = 0; i < n_points; i++) {
      BinomialParamterValues<degree> p( i * 1./(n_points - 1));
      for (int j = 0; j < degree + 1; j++) {
        coefficients[j] = binomial_coeffs[j] * p.t_values[j] * p.one_minus_t_values[j];
      }
      for (int k = 0; k < dimensions; k++) {
        points[i][k] = 0;
        for (int j = 0; j < degree + 1; j++) {
          points[i][k] += coefficients[j] * weights[j][k];
        }
      }
    }
  }

  void print() const {
    for (int i = 0; i < n_points; i++) {
      for (int j = 0; j < dimensions; j++) {
        std::cout <<  points[i][j] << ", ";
      }
      std::cout << '\n';
    }
  }

  BezierCurve<degree-1, n_points, dimensions> get_derivative() {
    if constexpr (degree > 0) {
      std::array<double, degree> derivative_weights;
      for (int i = 0; i < degree; i++) {
        derivative_weights[i] = degree * (weights[i+1] - weights[i]);
      }
      return BezierCurve<degree-1, n_points, dimensions>(derivative_weights);
    } else {
      return BezierCurve<0, n_points, dimensions>();
    }
  }

  template<size_t N>
  static std::array<std::array<double, degree+1>, degree+1-N> get_nth_derivative_to_current_weights_relation() {
    if constexpr (N >= 1) {
      assert(degree+1-N+1 > 0);
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

  std::array<Point<dimensions>, degree + 1> get_weights() const {
    return weights;
  }
private:
  std::array<Point<dimensions>, degree + 1> weights;
  std::array<Point<dimensions>, n_points> points;
  std::array<double, degree+1> coefficients;
};

#endif  // BEZIER_CURVE_HPP
