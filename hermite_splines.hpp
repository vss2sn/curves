#include <iostream>

#include "polynomials.hpp"

template<size_t degree, size_t dimensions>
class HermiteSplines {
public:
  explicit HermiteSplines(const std::array<Point<dimensions>, degree + 1> & p) : p(p) {
    std::array<double, degree+1> coeffs;
    std::fill(coeffs.begin(), coeffs.end(), 1);
    std::array<Polynomial<degree>, degree+1> polys;
    Polynomial<degree> poly(coeffs);
    auto coeffecient_matrix_of_p_to_dnp = get_coefficients_of_poly_and_all_derivatives(poly); // dnp = (d)^n p
    compile_for<0>(coeffecient_matrix_of_p_to_dnp);
    print(coefficients);
    std::cout << "Done" << '\n';
  }

private:
  std::array<Point<dimensions>, degree + 1> p;
  std::array<std::array<double, degree+1>, degree+1> coefficients;

  template<size_t I>
  void compile_for(std::array<std::array<double, degree+1>, degree+1>& coeffecient_matrix_of_p_to_dnp) {
    if constexpr  (degree < I) {
      return;
    } else {
      std::cout << I/2 << ' ' << I%2 << '\n';
      std::array<double, degree-I/2 + 1> temp;
      std::copy(
        std::begin(coeffecient_matrix_of_p_to_dnp[I/2]) + I/2,
        std::end(coeffecient_matrix_of_p_to_dnp[I/2]),
        std::begin(temp)
      );
      auto poly_at_i = Polynomial<degree - I/2>(temp);
      auto temp2 = poly_at_i.get_component_value(double(I%2));
      for (int i = 0; i < temp2.size(); i++) {
        coefficients[I][i] = temp2[i];
      }
      for (int i = temp2.size(); i < coefficients[I].size(); i++) {
        coefficients[I][i] = 0;
      }
      compile_for<I+1>(coeffecient_matrix_of_p_to_dnp);
    }
  }

  // template<>
  // void compile_for<(degree+1)/2>(std::array<std::array<double, degree+1>, degree+1>& coeffecient_matrix_of_p_to_dnp) {}
};
