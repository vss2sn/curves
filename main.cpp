#include <iostream>
#include "bezier_curve.hpp"

int main() {
  std::array<Point<3>, 4> weights = {
    Point<3>{110,150, 0},
    Point<3>{25,190,1},
    Point<3>{210,250,2},
    Point<3>{210, 30,3}
  };
  BezierCurve<3, 10, 3> b (weights);
  b.print();
  return 0;
}
