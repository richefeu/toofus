#ifndef RIGIDITY_HPP
#define RIGIDITY_HPP

#include "mat9.hpp"

class Rigidity {
public:
  double C1, C2, C3;

  // C1 C2 C2  0  0  0
  // C2 C1 C2  0  0  0
  // C2 C2 C1  0  0  0
  // 0  0  0	 C3 0  0
  // 0  0  0   0  C3 0
  // 0	0  0   0  0  C3

  Rigidity() { set(5e6, 0.3); }
  Rigidity(double Young, double Poisson) { set(Young, Poisson); }

  void set(double Young, double Poisson) {
    double factor = Young / ((1.0 + Poisson) * (1.0 - 2.0 * Poisson));
    C1 = factor * (1.0 - Poisson);
    C2 = factor * Poisson;
    C3 = factor * (1.0 - 2.0 * Poisson) / 2.0;
  }

  mat9r getStress(const mat9r &epsilon) {
    double Sxx = C1 * epsilon.xx + C2 * (epsilon.yy + epsilon.zz);
    double Syy = C1 * epsilon.yy + C2 * (epsilon.xx + epsilon.zz);
    double Szz = C1 * epsilon.zz + C2 * (epsilon.xx + epsilon.yy);
    double Syz = 2.0 * C3 * epsilon.yz;
    double Sxz = 2.0 * C3 * epsilon.xz;
    double Sxy = 2.0 * C3 * epsilon.xy;
    // clang-format off
		return mat9r(
			Sxx, Sxy, Sxz,
		  Sxy, Syy, Syz,
		  Sxz, Syz, Szz
		);
    // clang-format on
  }

  double bigDenum(const mat9r &left, const mat9r &right) {
    mat9r right_part = getStress(right); // C:right
		
		// left:C:right
    // clang-format off		
		return (
			left.xx * right_part.xx +
			left.yy * right_part.yy +
			left.zz * right_part.zz +
			2.0 * left.yz * right_part.yz +
			2.0 * left.xz * right_part.xz +
			2.0 * left.xy * right_part.xy
		);
    // clang-format on
  }
};

#endif /* end of include guard: RIGIDITY_HPP */

#if 0
#include <iostream>
int main(int argc, char const *argv[]) {
  Rigidity C(1e6, 0.3);
  mat9r eps(0.1, 0, 0, 0, 0, 0, 0, 0, 0);
	mat9r Sig = C.getStress(eps);
	mat9r dfds(0.1, 0.2, 0.1, 0.6, 0.3, 0.2, 0.1, 0.4, 0.1);
	
	std::cout << Sig << '\n';
	std::cout << C.bigDenum(dfds, Sig) << '\n';
  return 0;
}

#endif
