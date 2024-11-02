// Copyright (C) <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of TOOFUS (TOols OFten USued)
//
// It can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

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

  /**
   * @brief Default constructor for the Rigidity class.
   *
   * @details Initializes the Rigidity object with default values for Young's modulus and Poisson's ratio.
   * The default Young's modulus is set to 5e6, and the default Poisson's ratio is set to 0.3.
   */
  Rigidity() { set(5e6, 0.3); }

  /**
   * @brief Constructs a Rigidity object with specified material properties.
   *
   * @param Young Young's modulus of the material.
   * @param Poisson Poisson's ratio of the material.
   *
   * @details This constructor initializes the Rigidity object using the provided
   * Young's modulus and Poisson's ratio, setting up internal parameters based on
   * these values.
   */
  Rigidity(double Young, double Poisson) { set(Young, Poisson); }

  /**
   * @brief Sets the material properties for the Rigidity object.
   *
   * @param Young Young's modulus of the material.
   * @param Poisson Poisson's ratio of the material.
   *
   * @details This function sets the material properties for the Rigidity object,
   * given the specified Young's modulus and Poisson's ratio.
   *
   * The internal parameters of the Rigidity object are then set up based on
   * these material properties.
   */
  void set(double Young, double Poisson) {
    double factor = Young / ((1.0 + Poisson) * (1.0 - 2.0 * Poisson));
    C1 = factor * (1.0 - Poisson);
    C2 = factor * Poisson;
    C3 = factor * (1.0 - 2.0 * Poisson) / 2.0;
  }

  /**
   * @brief Computes the stress matrix from the strain matrix.
   *
   * @param epsilon Strain matrix.
   * @return Stress matrix.
   *
   * @details This function calculates the stress matrix based on the provided strain matrix
   * using the material properties defined within the Rigidity class. The stress components
   * are determined by multiplying the strain components with the rigidity coefficients.
   */
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

  /**
   * @brief Computes the big denominator of the compliance matrix.
   *
   * @param left Left matrix part of the compliance matrix.
   * @param right Right matrix part of the compliance matrix.
   * @return The big denominator of the compliance matrix.
   *
   * @details This function calculates the big denominator component of the compliance matrix
   * by computing the dot product of the two provided matrices. The result is the big denominator
   * component of the compliance matrix.
   */
  double bigDenum(const mat9r &left, const mat9r &right) {
    mat9r right_part = getStress(right); // C:right

    // left:C:right
    // clang-format off
    return (left.xx * right_part.xx + left.yy * right_part.yy + left.zz * right_part.zz +
            2.0 * left.yz * right_part.yz + 2.0 * left.xz * right_part.xz + 2.0 * left.xy * right_part.xy);
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
