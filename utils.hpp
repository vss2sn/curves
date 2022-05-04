#ifndef UTILS_H
#define UTILS_H
#include <algorithm>
#include <cmath>
#include <numeric>
#include <array>
#include <iostream>

template<size_t N>
using Point = std::array<double, N>;

template<typename T, size_t N, size_t M>
void debug_print(const std::array<std::array<T, N>, M>& a){
  for (const auto& p : a) {
    std::cout << "( ";
    for (const auto ele : p) {
      std::cout << ele << ',';
    }
    std::cout << ")\n";
  }
  std::cout << '\n';
}

template<typename T, size_t N>
void debug_print(const std::array<T, N>& a){
  for (const auto ele : a) {
    std::cout << ele << ',';
  }
  std::cout << '\n';
}

template<size_t N>
constexpr std::array<std::array<double, N>, N> inverse_using_LU_decomp(std::array<std::array<double, N>, N> mat) {
  std::array<std::array<double, N>, N> p;
  std::array<std::array<double, N>, N> inversed_data;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i == j) p[i][j] = 1;
      else p[i][j] = 0;
    }
  }
  for (int i = 0; i < N - 1; i++){
    auto max_idx = i;
    for (int k  = i+1; k < N; k++) {
      if (std::abs(mat[k][i]) > std::abs(mat[max_idx][i])) {
        max_idx = k;
      }
    }

    std::swap(mat[i], mat[max_idx]);
    std::swap(p[i], p[max_idx]);

    for (int j = i+1; j < N; j++) {
      mat[j][i] /= mat[i][i];
      for (int k = i+1; k < N; k++) {
        mat[j][k] -= mat[j][i] * mat[i][k];
      }
    }
  }

  for (int i_main = 0; i_main < N; i_main++) {
    std::array<double, N> b;
    for (int j = 0; j < N; j++) {
      if (i_main == j) {
        b[j] = 1;
      } else {
        b[j] = 0;
      }
    }

    auto temp = b;
    for (int i = 0; i < N; i++) {
      temp[i] = std::inner_product(std::begin(p[i]), std::end(p[i]), std::begin(b), 0.);
    }
    b = temp;

    for (int i = 0; i < N-1; i++){
      for (int j = i+1; j < N; j++)  {
        b[j] -= mat[j][i] * b[i];
      }
    }

    for (int i = N-1; i >=0; i--) {
      b[i] /= mat[i][i];
      for (int j = 0; j < i; j++) {
        b[j] -= mat[j][i] * b[i];
      }
    }

    for (int k =0; k < N; k++) {
      inversed_data[k][i_main] = b[k];
    }
  }
  return inversed_data;
}

template<typename T, size_t A, size_t B, size_t C>
std::array<std::array<T, C>, A> multiply_two_matrices(
  const std::array<std::array<T, B>, A>& m1,
  const std::array<std::array<T, C>, B>& m2)
{
  std::array<std::array<T, C>, A> ans;
  for(auto& row : ans) {
    for (auto& ele : row) {
      ele = 0;
    }
  }

  for (int i = 0; i < A; i++) {
    for (int j = 0; j < B; j++) {
      for (int k = 0; k <C; k++) {
        ans[i][k] +=  m1[i][j] * m2[j][k];
      }
    }
  }
  return ans;
}

#endif  // UTILS_H
