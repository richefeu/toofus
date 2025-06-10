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

#ifndef COMPLIANCE_HPP
#define COMPLIANCE_HPP

#include "mat9.hpp"

/**
 * @brief Compliance matrix
 *
 * Compliance matrix is the inverse of the Rigidity matrix.
 * It is used to compute the strain from the stress.
 * The compliance matrix is a 3x3 symmetric matrix.
 */
class Compliance {
public:
  double Cinv1, Cinv2, Cinv3;

  // Cinv1 Cinv2 Cinv2  0     0     0
  // Cinv2 Cinv1 Cinv2  0     0     0
  // Cinv2 Cinv2 Cinv1  0     0     0
  // 0     0     0	    Cinv3 0     0
  // 0     0     0      0     Cinv3 0
  // 0	   0     0      0     0     Cinv3

  /**
   * @brief Constructor with default values
   */
  Compliance() { set(5e6, 0.3); }
  /**
   * @brief Constructor with given values
   * @param Young Young modulus
   * @param Poisson Poisson's ratio
   */
  Compliance(double Young, double Poisson) { set(Young, Poisson); }

  /**
   * @brief Set the compliance matrix
   * @param Young Young modulus
   * @param Poisson Poisson's ratio
   */
  void set(double Young, double Poisson) {
    double factor = 1.0 / Young;
    Cinv1 = factor;
    Cinv2 = -factor * Poisson;
    Cinv3 = factor * (1.0 + Poisson); // Il faudrait multiplier par 2 ????????????
  }

  /**
   * @brief Compute the strain from the stress
   * @param stress Stress matrix
   * @return Strain matrix
   */
  mat9r getStrain(const mat9r &stress) {
    double Exx = Cinv1 * stress.xx + Cinv2 * (stress.yy + stress.zz);
    double Eyy = Cinv1 * stress.yy + Cinv2 * (stress.xx + stress.zz);
    double Ezz = Cinv1 * stress.zz + Cinv2 * (stress.xx + stress.yy);
    double Eyz = Cinv3 * stress.yz;
    double Exz = Cinv3 * stress.xz;
    double Exy = Cinv3 * stress.xy;
    // clang-format off
		return mat9r(
			Exx, Exy, Exz,
		  Exy, Eyy, Eyz,
		  Exz, Eyz, Ezz
		);
    // clang-format on
  }
};

#endif /* end of include guard: COMPLIANCE_HPP */

#if 0
#include <iostream>
#include "Rigidity.hpp"
int main(int argc, char const *argv[]) {
  Rigidity C(1e6, 0.3);
  Compliance Cinv(1e6, 0.3);
  mat9r eps(0.1, 0.01, 0, 0.01, 0.12, 0, 0, 0, -0.01);
	mat9r Sig = C.getStress(eps);
	mat9r dfds(0.1, 0.2, 0.1, 0.6, 0.3, 0.2, 0.1, 0.4, 0.1);
	
	std::cout << Sig << '\n';
	std::cout << C.bigDenum(dfds, Sig) << '\n';
	
	std::cout << eps << '\n';
	std::cout << Cinv.getStrain(Sig) << '\n'; 
  return 0;
}

#endif

