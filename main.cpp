#include <iostream>

#include "bezier_curve.hpp"
#include "cubic_hermite_spline.hpp"
#include "catmull_rom_spline.hpp"
#include "curve_conversion.hpp"

template<size_t dimensions>
void print(std::array<Point<dimensions>, 4>& a){
  for (const auto& p : a) {
    std::cout << "( ";
    for (const auto ele : p) {
      std::cout << ele << ',';
    }
    std::cout << ")";
  }
  std::cout << '\n';
}

// int main() {
//   std::array<Point<3>, 4> weights = {
//     Point<3>{110,150, 0},
//     Point<3>{25,190,1},
//     Point<3>{210,250,2},
//     Point<3>{210, 30,3}
//   };
//   BezierCurve<3, 10, 3> b (weights);
//   b.print();
//   return 0;
// }

// int main() {
//   std::array<Point<2>, 4> points = {
//     Point<2>{0,0},
//     Point<2>{1,1},
//     Point<2>{-1,0},
//     Point<2>{1,0}
//   };
//
//   CubicHermiteSpline<10, 2> hs (points);
//   hs.print();
//   return 0;
// }

// int main() {
//   constexpr int n_points = 100;
//   constexpr int dimensions = 2;
//   std::array<Point<dimensions>, 4> weights = {
//     Point<dimensions>{110,150},//0},
//     Point<dimensions>{25,190}, //1},
//     Point<dimensions>{210,250}, //2},
//     Point<dimensions>{210, 30} //3}
//   };
//   print(weights);
//
//   BezierCurve<3, n_points, dimensions> b (weights);
//   b.print();
//   std::cout << '\n';
//
//   std::array<Point<dimensions>, 4>  hermite_points = convertBezierPointsToHermitePoints(weights);
//   print(hermite_points);
//   CubicHermiteSpline<n_points, dimensions> hs (hermite_points);
//   hs.print();
//   std::cout << '\n';
//
//   std::array<Point<dimensions>, 4>  weights_returned = convertHermitePointsToBezierPoints(hermite_points);
//   print(weights_returned);
//   BezierCurve<3, n_points, dimensions> b_returned (weights_returned);
//   b_returned.print();
//   std::cout << '\n';
//
//   return 0
// }

int main() {
  constexpr int n_points = 10;
  constexpr int dimensions = 2;
  constexpr int N = 9;
  std::array<Point<dimensions>, N> p = {
      Point<dimensions>{-1., -1.5},
      Point<dimensions>{0, 0.5},
      Point<dimensions>{1, 1.5},
      Point<dimensions>{4, 3.2},
      Point<dimensions>{3.5, -0.1},
      Point<dimensions>{2.5, 0.7},
      Point<dimensions>{-1., -1.5},
      Point<dimensions>{0, 0.5},
      Point<dimensions>{1, 1.5}
  };

  const double tao = 0.5;
  CatmullRomSpline<N, n_points, dimensions> cms(p, tao);
  cms.print();

  return 0;
}
