#ifndef POLYNOMIALS_HPP
#define POLYNOMIALS_HPP

#include "utils.hpp"

template<size_t degree>
class Polynomial {
public:
  constexpr explicit Polynomial(const std::array<double, degree+1>& coefficients = {}) : coefficients(coefficients) {}

  constexpr double get_value(const double x) const {
    double ans = 0;
    double variable = 1;
    for (auto it= coefficients.rbegin(); it != coefficients.rend(); it++) {
      ans += (*it) * variable;
      variable *= x;
    }
    return ans;
  }

  constexpr std::array<double, degree + 1> get_component_value(const double x) const {
    std::array<double, degree + 1> ans;
    double variable = 1;
    for (int i = degree; i >= 0; i--) {
      ans[i] = coefficients[i] * variable;
      variable *= x;
    }
    return ans;
  }

  constexpr Polynomial<degree - 1> get_derivative () const {
    std::array<double, degree> new_coefficients;
    for (int i = 0; i < degree; i++) {
      new_coefficients[i] = coefficients[i] * (degree - i);
    }
    return Polynomial<degree - 1>(new_coefficients);
  }

  constexpr std::array<double, degree + 1> get_coefficients () const {
    return coefficients;
  }

  constexpr size_t get_degree() const {
    return degree;
  }

  void print() const {
    for (int i = 0; i < degree; i++) {
      std::cout << coefficients[i] << "x" << degree - i << " + ";
    }
    std::cout << coefficients[degree] << '\n';
  }

private:
  std::array<double, degree + 1> coefficients;
};

template<>
class Polynomial<0> {
public:
  constexpr explicit Polynomial(const std::array<double, 1>& coefficients = {0}) : coefficients(coefficients) {}

  constexpr double get_value(const double x) const {
    return coefficients[0];
  }

  constexpr std::array<double, 1> get_component_value(const double x) const {
    return coefficients;
  }

  constexpr Polynomial<0> get_derivative () const {
    return Polynomial<0>({0});
  }

  constexpr std::array<double, 1> get_coefficients () const {
    return coefficients;
  }

  constexpr size_t get_degree() const {
    return 0;
  }

  void print() const {
    std::cout << coefficients[0] << '\n';
  }

private:
  std::array<double, 1> coefficients;
};

template<size_t degree>
constexpr std::array<std::array<double, degree+1>, degree+1> get_coefficients_of_poly_and_all_derivatives(const Polynomial<degree>& poly) {
  std::array<std::array<double, degree+1>, degree+1> ans;
  ans[0] = poly.get_coefficients();
  const auto poly_d = poly.get_derivative();
  const auto lower_order_coeffs = get_coefficients_of_poly_and_all_derivatives(poly_d);
  for (int row = 1; row < degree+1; row++) {
    ans[row][0] = 0;
  }
  for (int row = 1; row < degree+1; row++) {
    for (int col = 1; col < degree+1; col++) {
      ans[row][col] = lower_order_coeffs[row - 1][col-1];
    }
  }
  return ans;
}

template<>
constexpr std::array<std::array<double, 1>, 1> get_coefficients_of_poly_and_all_derivatives<0>(const Polynomial<0>& poly) {
  return {0};
}

#endif  // POLYNOMIALS_HPP
