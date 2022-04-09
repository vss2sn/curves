#ifndef BEZIER_CURVE_HPP
#define BEZIER_CURVE_HPP

#include "binomial.hpp"
#include "utils.hpp"

template <size_t degree, size_t n_points, size_t dimensions>
class BezierCurve {
public:
  explicit BezierCurve(const std::array<Point<dimensions>, degree +1>& weights) noexcept : weights(weights) {
    std::array<int, degree + 1> binomial_coeffs = find_all_binomial_coefficients<degree>();
    for (int i = 0; i < n_points; i++) {
      BinomialParamterValues<degree> p( i * 1./(n_points - 1));
      for (int k = 0; k < dimensions; k++) {
        points[i][k] = 0;
        for (int j = 0; j < degree + 1; j++) {
          points[i][k] += binomial_coeffs[j] * p.t_values[j] * p.one_minus_t_values[j] * weights[j][k];
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

private:
  std::array<Point<dimensions>, degree + 1> weights;
  std::array<Point<dimensions>, n_points> points;
  std::array<double, n_points> coefficients;
};

#endif  // BEZIER_CURVE_HPP
