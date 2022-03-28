#ifndef CATMULL_ROM_SPLINE_HPP
#define CATMULL_ROM_SPLINE_HPP

#include <array>

// Note: a catmull-rum spline iss basically
// a set of hermite splines connecting sequental points,
// using tao to calculate the slopes at each ot the intermediate points

// n_points here is n_points per segment not total number of points
template<size_t N, size_t n_points, size_t dimensions>
class CatmullRomSpline {
public:
  explicit CatmullRomSpline(const std::array<Point<dimensions>, N>& p, const double tao) : p(p), tao(tao) {
    const double interval = (1.)/n_points;
    std::array<std::array<double, 4>,  n_points+1> coefficients;
    const std::array<std::array<double,4>, 4> tao_matrix {
      std::array<double,4>{0., 1., 0., 0.},
      std::array<double,4>{-tao, 0., tao, 0.},
      std::array<double,4>{2.*tao, tao-3., 3.-2.*tao, -tao},
      std::array<double,4>{-tao, 2.-tao, tao-2., tao}
    };
    double u = 0;
    for (int i = 0; i < n_points + 1; i++) {
      const std::array<double, 4> u_matrix {1, u, u*u, u*u*u};
      for (int j = 0; j < 4; j++) {
        coefficients[i][j] = 0;
        for (int k = 0; k < 4;  k++) {
          coefficients[i][j] += u_matrix[k] * tao_matrix[k][j];
        }
      }
      u += interval;
    }

    for (int i = 0; i <= N-4; i++) {
      for (int j = 0; j < n_points + 1; j++) {
        for (int m = 0; m < dimensions; m++) {
          points[(n_points + 1) * i+j][m] = 0;
          for (int k = 0; k < 4; k++) {
            points[(n_points + 1) * i+j][m] += coefficients[j][k] * p[i+k][m];
          }
        }
      }
    }
  }

  void print() {
    for (const auto& p : points) {
      for (int j = 0; j < dimensions; j++) {
        std::cout << p[j] << ", ";
      }
      std::cout << std::endl;
    }
  }

private:
  std::array<Point<dimensions>, N> p;
  std::array<Point<dimensions>, (N-3)*(n_points+1)> points;
  double tao;
};

#endif  // CATMULL_ROM_SPLINE_HPP
