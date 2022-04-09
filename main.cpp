#include <array>
#include <iostream>

#include "bezier_curve.hpp"
#include "cubic_hermite_spline.hpp"
#include "catmull_rom_spline.hpp"
#include "curve_conversion.hpp"
#include "hermite_splines.hpp"
#include "polynomials.hpp"
#include "utils.hpp"


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
//   constexpr int n_points = 5;
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
//   return 0;
// }

// int main() {
//   constexpr int n_points = 10;
//   constexpr int dimensions = 2;
//   constexpr int N = 7;
//
//   std::array<Point<dimensions>, N> p = {
//       // Point<dimensions>{-4., 1.},
//       Point<dimensions>{-4., 1.},
//       Point<dimensions>{-4., 1.},
//       Point<dimensions>{-2., -2.},
//       Point<dimensions>{0, 0},
//       Point<dimensions>{2, -3},
//       Point<dimensions>{3, 1},
//       Point<dimensions>{3, 1},
//       // Point<dimensions>{3, 1},
//   };
//
//   const double tao = 1;
//   CatmullRomSpline<N, n_points, dimensions> cms(p, tao);
//   cms.print();
//
//   return 0;
// }

int main () {
  constexpr size_t degree = 3;
  constexpr size_t dimensions = 2;
  constexpr size_t n_points = 10;
  constexpr std::array<Point<dimensions>, degree+1> p = {
    Point<dimensions>{0., 0.},
    Point<dimensions>{4., 1.},
    Point<dimensions>{0., 1.},
    Point<dimensions>{0., 1.},
    // Point<dimensions>{0., 1.},
    // Point<dimensions>{0., 1.},
  };

  const auto chs = CubicHermiteSpline<n_points, dimensions>(p);
  std::cout << "Cubic Hermite Spline" << '\n';
  chs.print();
  std::cout << '\n';

  const auto hs = HermiteSplines<degree, n_points, dimensions>(p);
  std::cout << "Hermite Spline" << '\n';
  hs.print();
  std::cout << '\n';
  return 0;
}
