#include <iostream>


#include "polynomials.hpp"
#include "utils.hpp"

template<size_t degree, size_t n_points, size_t dimensions>
class HermiteSplines {
public:
  explicit HermiteSplines(const std::array<Point<dimensions>, degree + 1> & p) : p(p) {
    calculate_coefficients();
    const double delta_u = 1./n_points;
    double u = 0;
    std::array<double, degree+1> powers_of_u;
    for (int i = 0; i <= n_points; i++) {
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
  }

  void print() {
    for (const auto& p : points) {
      for (int j = 0; j < dimensions; j++) {
        std::cout << p[j] << ", ";
      }
      std::cout << '\n';
    }
  }

private:
  std::array<Point<dimensions>, degree + 1> p;
  std::array<std::array<double, degree+1>, n_points + 1> coefficients;
  std::array<std::array<double, degree+1>, degree+1> coefficients_of_basis_curves;
  std::array<Point<dimensions>, n_points + 1> points;

  void calculate_coefficients() {
    std::array<double, degree+1> coeffs_for_ploy;
    std::fill(coeffs_for_ploy.begin(), coeffs_for_ploy.end(), 1);
    Polynomial<degree> poly(coeffs_for_ploy);
    std::array<Polynomial<degree>, degree+1> polys;
    auto coeffecient_matrix_of_p_to_dnp = get_coefficients_of_poly_and_all_derivatives(poly); // dnp = (d)^n p
    std::array<std::array<double, degree+1>, degree+1> coeffs;
    compile_for<0>(coeffecient_matrix_of_p_to_dnp, coeffs);
    coefficients_of_basis_curves = inverse_using_LU_decomp(coeffs);
  }

  template<size_t I>
  void compile_for(
    std::array<std::array<double, degree+1>, degree+1>& coeffecient_matrix_of_p_to_dnp,
    std::array<std::array<double, degree+1>, degree+1>& coeffs
  ) {
    if constexpr  (degree < I) {
      return;
    } else {
      std::array<double, degree-I/2 + 1> temp;
      std::copy(
        std::begin(coeffecient_matrix_of_p_to_dnp[I/2]) + I/2,
        std::end(coeffecient_matrix_of_p_to_dnp[I/2]),
        std::begin(temp)
      );
      auto poly_at_i = Polynomial<degree - I/2>(temp);
      auto temp2 = poly_at_i.get_component_value(double(I%2));
      for (int i = 0; i < temp2.size(); i++) {
        coeffs[I][i] = temp2[i];
      }
      for (int i = temp2.size(); i < coeffs[I].size(); i++) {
        coeffs[I][i] = 0;
      }
      compile_for<I+1>(coeffecient_matrix_of_p_to_dnp, coeffs);
    }
  }
};
