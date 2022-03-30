#ifndef UTILS_H
#define UTILS_H

#include <array>
#include <iostream>

template<size_t N>
using Point = std::array<double, N>;

// template<size_t dimensions>
// void print(std::array<Point<dimensions>, 4>& a){
//   for (const auto& p : a) {
//     std::cout << "( ";
//     for (const auto ele : p) {
//       std::cout << ele << ',';
//     }
//     std::cout << ")";
//   }
//   std::cout << '\n';
// }


template<size_t N, size_t M>
void print(std::array<std::array<double, N>, M>& a){
  for (const auto& p : a) {
    std::cout << "( ";
    for (const auto ele : p) {
      std::cout << ele << ',';
    }
    std::cout << ")\n";
  }
  std::cout << '\n';
}

#endif  // UTILS_H
