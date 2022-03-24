#include <iostream>
#include "bezier_curve.hpp"
#include "cubic_hermite_spline.hpp"

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

int main() {
  std::array<Point<2>, 4> points = {
    Point<2>{0,0},
    Point<2>{1,1},
    Point<2>{-1,0},
    Point<2>{1,0}
  };

  CubicHermiteSpline<10, 2> hs (points);
  hs.print();
  return 0;
}
