#ifndef HERMITE_SPLINES_HPP
#define HERMITE_SPLINES_HPP

#include <iostream>

#include "polynomials.hpp"
#include "utils.hpp"

template<size_t degree, size_t n_points, size_t dimensions>
class HermiteSplines {
public:
  explicit HermiteSplines(const std::array<Point<dimensions>, degree + 1> & p) : p(p) {
    calculate_coefficients();
    const double delta_u = 1./(n_points - 1);
    double u = 0;
    std::array<double, degree+1> powers_of_u;
    // std::cout << __LINE__ << '\n';
    for (int i = 0; i <= (n_points - 1); i++) {
      powers_of_u[degree] = 1;
      for (int i = degree-1; i >=0; i--) {
        powers_of_u[i] = powers_of_u[i+1] * u;
      }
      for (int m = 0; m < degree+1; m++) {
        coefficients[i][m] = 0;
        for (int n = 0; n < degree+1; n++) {
          coefficients[i][m] += powers_of_u[n] * coefficients_of_basis_curves[n][m];
        }
      }
      for (int j = 0; j < dimensions; j++) {
        points[i][j] = 0;
        for (int k = 0; k < degree+1; k++) {
          points[i][j] += coefficients[i][k] * p[k][j];
        }
      }
      u += delta_u;
    }
    // std::cout << __LINE__ << '\n';
  }

  void print() const {
    for (const auto& p : points) {
      for (int j = 0; j < dimensions; j++) {
        std::cout << p[j] << ", ";
      }
      std::cout << '\n';
    }
  }

private:
  std::array<Point<dimensions>, degree + 1> p;
  std::array<std::array<double, degree+1>, n_points> coefficients;
  std::array<std::array<double, degree+1>, degree+1> coefficients_of_basis_curves;
  std::array<Point<dimensions>, n_points> points;

  void calculate_coefficients() {
    // Explanation:
    // Given a polynomial Cn * x^(n) + Cn-1 * x^(n-1) + ... + C0
    // Need to find the coeffs Cn, Cn-1, ..., C0
    // The derivate of the polynomial are n * Cn * x^(n-1) + (n-1) * Cn-1 * x^(n-2) + ... + C1
    // The compile time for sustitutes the values of the variable at the start and end points of
    // for these polynomials starting with degree n to degree n/2 and then sthe start point at degree n/2-1 if
    // n is odd. Note that right now the  first 2 rows of the coeffecient_matrix_of_p_to_dnp would be [[Cn, Cn-1, ..., C0], [0, n * Cn, (n-1) * Cn-1, ..., C1], ...]
    // The compile for loop converts it to [[Cn, Cn-1, ..., C0], [n * Cn, (n-1) * Cn-1, ..., C1, 0], ...]
    // Which can then be used in AX=b for to calculate the coefficients. Here X is the column matrix [Cn, C-1, ... C0]
    std::array<double, degree+1> coeffs_for_ploy;
    std::fill(coeffs_for_ploy.begin(), coeffs_for_ploy.end(), 1);
    const Polynomial<degree> poly(coeffs_for_ploy);
    std::array<Polynomial<degree>, degree+1> polys;
    const auto coeffecient_matrix_of_p_to_dnp = get_coefficients_of_poly_and_all_derivatives(poly); // dnp = (d)^n p
    std::array<std::array<double, degree+1>, degree+1> coeffs;
    compile_for<0>(coeffecient_matrix_of_p_to_dnp, coeffs);
    coefficients_of_basis_curves = inverse_using_LU_decomp(coeffs);
  }

  template<size_t I>
  void compile_for(
    const std::array<std::array<double, degree+1>, degree+1>& coeffecient_matrix_of_p_to_dnp,
    std::array<std::array<double, degree+1>, degree+1>& coeffs
  ) {
    if constexpr (degree < I) {
      return;
    } else {
      std::array<double, degree-I/2 + 1> coeffs_at_row_I_by_two;
      std::copy(
        std::begin(coeffecient_matrix_of_p_to_dnp[I/2]) + I/2,
        std::end(coeffecient_matrix_of_p_to_dnp[I/2]),
        std::begin(coeffs_at_row_I_by_two)
      );
      const auto poly_at_i = Polynomial<degree - I/2>(coeffs_at_row_I_by_two);
      const auto components_at_row_I_by_two = poly_at_i.get_component_value(double(I%2));
      std::copy(std::begin(components_at_row_I_by_two), std::end(components_at_row_I_by_two), std::begin(coeffs[I]));
      std::fill(std::begin(coeffs[I]) + components_at_row_I_by_two.size(), std::end(coeffs[I]), 0.);
      compile_for<I+1>(coeffecient_matrix_of_p_to_dnp, coeffs);
    }
  }
};
#endif // HERMITE_SPLINES_HPP
