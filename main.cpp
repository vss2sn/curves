#include <array>
#include <iostream>

#include "bezier_curve.hpp"
#include "b_spline_curve.hpp"
#include "catmull_rom_spline.hpp"
#include "cubic_hermite_spline.hpp"
#include "curve_conversion.hpp"
#include "hermite_splines.hpp"
#include "polynomials.hpp"
#include "utils.hpp"

// int main() {
//   constexpr size_t degree = 3;
//   std::array<Point<2>, degree+1> weights = {
//     Point<2>{-1.0,  0.0},
//     Point<2>{-0.5,  0.5},
//     Point<2>{ 0.5, -0.5},
//     Point<2>{ 1.0,  0.0},
//   };
//   BezierCurve<degree, 100, 2> b (weights);
//   b.print();
//   std::cout << '\n';
//
//
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
//
// int main () {
//   constexpr size_t degree = 3;
//   constexpr size_t dimensions = 2;
//   constexpr size_t n_points = 10;
//   constexpr std::array<Point<dimensions>, degree+1> p = {
//     Point<dimensions>{0., 0.},
//     Point<dimensions>{4., 1.},
//     Point<dimensions>{0., 1.},
//     Point<dimensions>{0., 1.},
//     // Point<dimensions>{0., 1.},
//     // Point<dimensions>{0., 1.},
//   };
//
//   const auto chs = CubicHermiteSpline<n_points, dimensions>(p);
//   std::cout << "Cubic Hermite Spline" << '\n';
//   chs.print();
//   std::cout << '\n';
//
//   constexpr auto hs = HermiteSplines<degree, n_points, dimensions>(p);
//   std::cout << "Hermite Spline" << '\n';
//   hs.print();
//   std::cout << '\n';
//   return 0;
// }
//
// int main () {
//   std::array<std::array<double, 4>, 1> m1 = {
//     {1,2,3,4}
//   };
//
//   std::array<std::array<double, 2>, 4> m2 = {
//     std::array<double, 2>{1,2},
//     std::array<double, 2>{1,2},
//     std::array<double, 2>{1,2},
//     std::array<double, 2>{1,2}
//   };
//
//   debug_print(multiply_two_matrices(m1, m2));
// }

// int main() {
//   constexpr size_t n_control_points = 4;
//   constexpr size_t degree = 3;
//   constexpr size_t dimensions = 2;
//   constexpr size_t n_points = 100;
//
//   std::array<std::array<double, dimensions>, n_control_points> control_points = {
//     Point<dimensions>{-1.0,  0.0},
//     Point<dimensions>{-0.5,  0.5},
//     Point<dimensions>{ 0.5, -0.5},
//     Point<dimensions>{ 1.0,  0.0},
//   };
//
//   // B-splines with clamped knot vectors pass through
//   // the two end control points.
//   // A clamped knot vector must have `degree + 1` equal knots
//   // at both its beginning and end.
//   // This is done to ensure that the bspline collapses into a bezier curve
//   // and a comparison and sanity check for the code can be run
//   std::array<int, n_control_points+degree+1> knots = {0,0,0,0,3,3,3,3};
//   std::array<int, n_control_points> weights;
//   std::fill(std::begin(weights), std::end(weights), 1);
//   BSplineCurve<degree, n_control_points, n_points, dimensions> bspline(control_points, knots, weights);
//   bspline.print();
//
//   // assert (degree + 1) == n_control_points;
//   BezierCurve<degree, n_points, dimensions> bezier (control_points);
//   bezier.print();
//
//   std::cout << '\n';
//   return 0;
// }

// constexpr size_t n_control_points = 4;
// constexpr size_t degree = 3;
// constexpr size_t dimensions = 2;
// constexpr size_t n_points = 10;
//
// // To check if compile time evaluation is performed
// constexpr BSplineCurve<degree, n_control_points, n_points, dimensions> f() {
//     if (std::is_constant_evaluated()) {
//         constexpr std::array<std::array<double, dimensions>, n_control_points> control_points = {
//             Point<dimensions>{-1.0,  0.0},
//             Point<dimensions>{-0.5,  0.5},
//             Point<dimensions>{ 0.5, -0.5},
//             Point<dimensions>{ 1.0,  0.0},
//         };
//
//         constexpr std::array<int, n_control_points+degree+1> knots = {0,0,0,0,3,3,3,3};
//         constexpr std::array<int, n_control_points> weights = {1,1,1,1};
//         constexpr BSplineCurve<degree, n_control_points, n_points, dimensions> bspline(control_points, knots, weights);
//         return bspline;
//     } else {
//       std::array<std::array<double, dimensions>, n_control_points> control_points = {
//           Point<dimensions>{-2.0,  0.0},
//           Point<dimensions>{-1.0,  1.0},
//           Point<dimensions>{ 1.0, -1.0},
//           Point<dimensions>{ 2.0,  0.0},
//       };
//       std::array<int, n_control_points+degree+1> knots = {0,0,0,0,3,3,3,3};
//       std::array<int, n_control_points> weights = {1,1,1,1};
//       BSplineCurve<degree, n_control_points, n_points, dimensions> bspline(control_points, knots, weights);
//       return bspline;
//     }
//
// }
//
// int main() {
//     constexpr auto bspline_constexpr = f();
//     bspline_constexpr.print();
//     std::cout << '\n';
//     auto bspline = f();
//     bspline.print();
//     return 0;
// }

// constexpr size_t degree = 3;
// constexpr size_t dimensions = 2;
// constexpr size_t n_points = 10;
//
// // To check if compile time evaluation is performed
// constexpr BezierCurve<degree, n_points, dimensions> f() {
//     if (std::is_constant_evaluated()) {
//         constexpr std::array<std::array<double, dimensions>, degree+1> weights = {
//             Point<dimensions>{-1.0,  0.0},
//             Point<dimensions>{-0.5,  0.5},
//             Point<dimensions>{ 0.5, -0.5},
//             Point<dimensions>{ 1.0,  0.0},
//         };
//
//         constexpr BezierCurve<degree, n_points, dimensions> bezier(weights);
//         return bezier;
//     } else {
//       std::array<std::array<double, dimensions>, degree+1> weights = {
//           Point<dimensions>{-2.0,  0.0},
//           Point<dimensions>{-1.0,  1.0},
//           Point<dimensions>{ 1.0, -1.0},
//           Point<dimensions>{ 2.0,  0.0},
//       };
//       BezierCurve<degree, n_points, dimensions> bezier(weights);
//       return bezier;
//     }
//
// }
//
// int main() {
//     constexpr auto bezier_constexpr = f();
//     bezier_constexpr.print();
//     std::cout << '\n';
//     auto bezier = f();
//     bezier.print();
//     return 0;
// }

// constexpr size_t degree = 3;
// constexpr size_t dimensions = 2;
// constexpr size_t n_points = 10;
//
// // To check if compile time evaluation is performed
// constexpr HermiteSplines<degree, n_points, dimensions> f() {
//     if (std::is_constant_evaluated()) {
//       constexpr std::array<Point<dimensions>, degree+1> p = {
//         Point<dimensions>{0., 0.},
//         Point<dimensions>{4., 1.},
//         Point<dimensions>{0., 1.},
//         Point<dimensions>{0., 1.},
//       };
//
//       constexpr auto hs = HermiteSplines<degree, n_points, dimensions>(p);
//       return hs;
//     } else {
//       constexpr std::array<Point<dimensions>, degree+1> p = {
//         Point<dimensions>{0., 0.},
//         Point<dimensions>{1., 1.},
//         Point<dimensions>{0., 1.},
//         Point<dimensions>{0., 1.},
//       };
//
//       constexpr auto hs = HermiteSplines<degree, n_points, dimensions>(p);
//       return hs;
//     }
//
// }

constexpr int n_points = 10;
constexpr int dimensions = 2;
constexpr int N = 7;
// To check if compile time evaluation is performed
constexpr CatmullRomSpline<N, n_points, dimensions> f() {
    if (std::is_constant_evaluated()) {
      constexpr std::array<Point<dimensions>, N> p = {
          Point<dimensions>{-4., 1.},
          Point<dimensions>{-4., 1.},
          Point<dimensions>{-2., -2.},
          Point<dimensions>{0, 0},
          Point<dimensions>{2, -3},
          Point<dimensions>{3, 1},
          Point<dimensions>{3, 1},
      };

      constexpr double tao = 1;
      constexpr auto cms = CatmullRomSpline<N, n_points, dimensions>(p, tao);
      return cms;
    } else {
      std::array<Point<dimensions>, N> p = {
          Point<dimensions>{-5., 1.},
          Point<dimensions>{-5., 1.},
          Point<dimensions>{-2., -2.},
          Point<dimensions>{0, 0},
          Point<dimensions>{2, -3},
          Point<dimensions>{4, 1},
          Point<dimensions>{4, 1},
      };

      const double tao = 1;
      auto cms = CatmullRomSpline<N, n_points, dimensions>(p, tao);
      return cms;
    }

}


int main() {
  constexpr auto cms_constexpr = f();
  cms_constexpr.print();
  std::cout << '\n';
  auto cms = f();
  cms.print();
  return 0;
}
