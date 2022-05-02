#include <array>
#include <iostream>

#include "utils.hpp"

template<size_t degree, size_t n_control_points, size_t n_points, size_t dimensions,
    std::enable_if_t<std::greater<size_t>()(degree, 0), bool> = true,
    std::enable_if_t<std::less<size_t>()(degree, n_control_points), bool> = true>
class BSplineCurve {
public:
  explicit BSplineCurve (
    const std::array<std::array<double, dimensions>, n_control_points>& control_points,
    std::array<int, n_control_points + degree + 1>& knots,
    std::array<int, n_control_points>& weights
  ) : control_points(control_points), knots(knots), weights(weights) {

    const std::array<size_t, 2> domain = {
      degree,
      knots.size() - 1 - degree
    };

    const auto low  = knots[domain[0]];
    const auto high = knots[domain[1]];

    auto homogeneous_coordinates = std::array<Point<dimensions+1>, n_control_points>();
    for(size_t i = 0; i < n_control_points; i++) {
      for(size_t j = 0; j < dimensions; j++) {
        homogeneous_coordinates[i][j] = control_points[i][j] * weights[i];
      }
      homogeneous_coordinates[i][dimensions] = weights[i];
    }

    double t = 0;
    const double delta = 1./(n_points - 1);
    for(size_t i = 0; i < n_points; i++, t += delta) {
      const double t_normalized = t * (high - low) + low;
      points[i] = interpolate(t_normalized, domain, homogeneous_coordinates);
    }
  }

  void print() const {
    for (size_t i = 0; i < n_points; i++) {
      for (size_t j = 0; j < dimensions; j++) {
        std::cout <<  points[i][j] << ", ";
      }
      std::cout << '\n';
    }
  }

  std::array<std::array<double, dimensions>, n_control_points> get_control_points() const {
    return control_points;
  }

private:

  Point<dimensions> interpolate(const double& t,
    const std::array<size_t, 2>& domain,
    std::array<Point<dimensions+1>, n_control_points> homogeneous_coordinates
  ) {

    size_t s = domain[0];
    for (; s < domain[1]; s++) {
      if(t >= knots[s] && t <= knots[s+1]) {
        break;
      }
    }

    for(size_t l = 1; l <= degree + 1; l++) {
      for(size_t i = s; i > s - degree - 1 + l; i--) {
        const double alpha = (t - knots[i]) / (knots[i + degree + 1 - l] - knots[i]);
        for(size_t j = 0; j < dimensions + 1; j++) {
          homogeneous_coordinates[i][j] = (1 - alpha) * homogeneous_coordinates[i-1][j] + alpha * homogeneous_coordinates[i][j];
        }
      }
    }

    Point<dimensions> result;
    for(size_t i = 0; i < dimensions; i++) {
      result[i] = homogeneous_coordinates[s][i] / homogeneous_coordinates[s][dimensions];
    }

    return result;
  }

  const std::array<std::array<double, dimensions>, n_control_points> control_points;
  const size_t order = degree + 1;
  const size_t n_knots = n_control_points + degree + 1;
  const std::array<int, n_control_points + degree + 1> knots;
  const std::array<int, n_control_points> weights;
  std::array<Point<dimensions>, n_points> points;
};
