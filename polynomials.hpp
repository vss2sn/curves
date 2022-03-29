
template<size_t degree>
class Polynomial {
public:
  explicit Polynomial(const std::array<double, degree+1>& coefficients) : coefficients(coefficients) {}

  double get_value(const double x) {
    double ans = 0;
    double variable = 1;
    for (auto it= coefficients.rbegin(); it != coefficients.rend(); it++) {
      ans += (*it) * variable;
      variable *= x;
    }
    return ans;
  }

  Polynomial<degree - 1> get_derivative () {
    std::array<double, degree> new_coefficients;
    for (int i = 0; i < degree; i++) {
      new_coefficients[i] = coefficients[i] * (degree - i);
    }
    return Polynomial<degree - 1>(new_coefficients);
  }

  std::array<double, degree + 1> get_coefficients () {
    return coefficients;
  }

  size_t get_degree() {
    return degree;
  }

  void print() {
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
  explicit Polynomial(const std::array<double, 1>& coefficients) : coefficients(coefficients) {}

  double get_value(const double x) {
    return coefficients[0];
  }

  Polynomial<0> get_derivative () {
    return Polynomial<0>({0});
  }

  std::array<double, 1> get_coefficients () {
    return coefficients;
  }

  size_t get_degree() {
    return 0;
  }

  void print() {
    std::cout << coefficients[0] << '\n';
  }

private:
  std::array<double, 1> coefficients;
};
