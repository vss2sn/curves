#ifndef BEZIER_CURVE_HPP
#define BEZIER_CURVE_HPP

#include "binomial.hpp"
#include "utils.hpp"

template<size_t N>
struct ParamterValues {
public:
  explicit ParamterValues(const float t) noexcept {
    t_values[0] = 1;
    one_minus_t_values[N] = 1;
    for (int i = 1; i <= N; i++) {
      t_values[i] = t_values[i-1] * t;
      one_minus_t_values[N-i] = one_minus_t_values[N-i+1] * (1-t);
    }
  }

  std::array<float, N+1> t_values;
  std::array<float, N+1> one_minus_t_values;
};

template <size_t degree, size_t n_points, size_t dimensions>
class BezierCurve {
public:
  explicit BezierCurve(const std::array<Point<dimensions>, degree +1>& weights) noexcept : weights(weights) {
    std::array<int, degree + 1> binomial_coeffs = find_all_binomial_coefficients<degree>();
    for (int i = 0; i < n_points + 1; i++) {
      ParamterValues<degree> p( i * 1./n_points);
      for (int k = 0; k < dimensions; k++) {
        points[i][k] = 0;
        for (int j = 0; j < degree + 1; j++) {
          points[i][k] += binomial_coeffs[j] * p.t_values[j] * p.one_minus_t_values[j] * weights[j][k];
        }
      }
    }
  }

void print() {
  for (int i = 0; i < n_points + 1; i++) {
    for (int j = 0; j < dimensions; j++) {
      std::cout <<  points[i][j] << ", ";
    }
    std::cout << '\n';
  }
}

private:
  std::array<Point<dimensions>, degree + 1> weights;
  std::array<Point<dimensions>, n_points + 1> points;
  std::array<float, n_points + 1> coefficients;
};

#endif  // BEZIER_CURVE_HPP
