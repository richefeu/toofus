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

#ifndef TRANSIT_FUNC_HPP
#define TRANSIT_FUNC_HPP

#include <cmath>

class linearToPlateau {
public:

  linearToPlateau() { set(1.0, 0.05, 0.01); }

  void set(double Plateau, double Transit, double Smoothness) {
    plateau = Plateau;
    transit = Transit;
    smoothness = Smoothness;

    abruptness = 1.0 / smoothness;
    Lfactor = plateau / transit;
    C = fabs(func(0));
    Cfactor = 1.0 / (1.0 + C);
  }

  double func(double x) {
    double F = Lfactor * (x - smoothness * log(1.0 + exp((x - transit) * abruptness)));
    return (F + C) * Cfactor;
  }

  double deriv(double x) {
    double expTransit = exp(transit * abruptness);
    double dfdx = plateau * expTransit / (transit * (exp(x * abruptness) + expTransit));
    return dfdx * Cfactor;
  }

private:
  // the 3 parameters
  double transit{0.0};
  double plateau{0.0};
  double smoothness{0.0};

  // some precomputed values
  double abruptness{0.0};
  double Lfactor{0.0};
  double C{0.0};
  double Cfactor{0.0};
};

class plateauToPlateau {
public:

  plateauToPlateau() { set(0.0, 1.0, 0.5, 0.02); }

  void set(double Plateau_beg, double Plateau_end, double Transit, double Smoothness) {
    plateau_beg = Plateau_beg;
    plateau_end = Plateau_end;
    transit = Transit;
    smoothness = Smoothness;

    abruptness = 1.0 / smoothness;
    Ampl = plateau_beg - plateau_end;
  }

  double func(double x) {
    return Ampl * (exp(-transit * abruptness) + 1.0) / (exp((x - transit) * abruptness) + 1.0) + plateau_end;
  }

  double deriv(double x) {
    double a = exp((x - transit) * abruptness);
    double a1 = a + 1.0;
    return -Ampl * (exp(-transit * abruptness) + 1.0) * a / (smoothness * a1 * a1);
  }

private:
  // the 3 parameters
  double transit{0.0};
  double plateau_beg{0.0};
  double plateau_end{0.0};
  double smoothness{0.0};

  // some precomputed values
  double abruptness{0.0};
  double Ampl{0.0};
};

#endif /* end of include guard: TRANSIT_FUNC_HPP */

#if 0
#include <fstream>
#include <iostream>

int main(int argc, char const *argv[]) {
  std::ofstream file("test.txt");

  linearToPlateau evol;
  evol.set(1.0, 0.05, 0.01); // double Plateau, double Transit, double Smoothness
  for (double x = 0.0; x < 1.0; x += 0.005) {
    file << x << "\t" << evol.func(x) << "\t" << evol.deriv(x) << std::endl;
    std::cout << x << "\t" << evol.func(x) << "\t" << evol.deriv(x) << std::endl;
  }

  /*
plateauToPlateau evol;
evol.set(0., 1., .5, 0.04);
for (double x = 0.0; x < 2.0; x += 0.01) {
file << x << "\t" << evol.func(x) << "\t" << evol.deriv(x) << std::endl;

*/
  return 0;
}

#endif
