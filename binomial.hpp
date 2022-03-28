#ifndef BINOMAIAL_HPP
#define BINOMAIAL_HPP

#include <array>
#include <iostream>

// optimize to use binomial property to reduce number of computations
template<int i, int N>
struct Binomial {
  static constexpr int val = static_cast<double>(Binomial<i-1, N-1>::val) * static_cast<double>(N) / static_cast<double>(i);
};

template<int N>
struct Binomial<0, N> {
  static constexpr int val = 1;
};

template<int N>
struct Binomial<N, N> {
  static constexpr int val = 1;
};

template<>
struct Binomial<0, 0> {
  static constexpr int val = 1;
};

template<size_t ... N>
std::array<int, sizeof ... (N) + 1> find_all_binomial_coefficients_impl (std::index_sequence<N ...>) {
    constexpr std::array<int, sizeof ... (N) + 1> a =
      {
        Binomial<N, sizeof ... (N)>::val ...,
        Binomial<sizeof ... (N), sizeof ... (N)>::val
      };
    return a;
}

template<int N>
std::array<int, N+1> find_all_binomial_coefficients() {
  return find_all_binomial_coefficients_impl(std::make_index_sequence<N>());
}

template<int N>
std::array<int, ((N+1)*(N+2))/2> pascals_triangle() {
  std::array<int, (N*(N+1))/2> a = pascals_triangle<N-1>();
  std::array<int, ((N+1)*(N+2))/2> ans;
  for (int i = 0; i < (N*(N+1))/2; i++) {
     ans[i] = a[i];
  }
  const auto b = find_all_binomial_coefficients<N>();
  for (int i = 0; i < N+1; i++) {
     ans[i + N*(N+1)/2] = b[i];
  }
  return ans;
}

template<>
std::array<int, 1> pascals_triangle<0>() {
  return find_all_binomial_coefficients<0>();
}

template<int start, int end>
void compile_time_for() {
  if constexpr (start < end) {
    const auto ans = find_all_binomial_coefficients<start>();
    for (const auto& ele : ans) {
      std::cout << ele << ' ';
    }
    std::cout << '\n';
    compile_time_for<start+1, end>();
  }
}

template<size_t N>
struct BinomialParamterValues {
public:
  explicit BinomialParamterValues(const float t) noexcept {
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

#endif  // BINOMAIAL_HPP

// int main() {
//   constexpr int N = 5;
//   // compile_time_for<0, N+1>();
//   for (const auto& ele : pascals_triangle<N>()) {
//     std::cout << ele << ' ';
//   }
//   std::cout << '\n';
//   return 0;
// }
