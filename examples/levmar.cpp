#include "../Levenberg-Marquardt.hpp"
#include <iostream>

int main() {

  const int m = 15;
  const int n = 3;

  int info;
  int lwa = levmar::get_lwa(m, n);
  int iwa[n];
  double tol;
  double x[n];
  double fvec[m];
  double wa[lwa];

  double y[15];
  for (size_t i = 0; i < m; i++) {
    double x = (double)i;
    y[i] = 0.2307638 * x * x * x * x - 5.7365402 * x * x * x + 0.0027007 +
           0.01 * ((rand() / RAND_MAX) - 0.5);
  }

  x[0] = 1.;
  x[1] = 1.;
  x[2] = 1.;

  auto fcn = [](void *p, int m, int n, const double *x, double *fvec,
                int iflag) {
    const double *y = (double *)p;
    (void)iflag;

    for (size_t i = 0; i < m; ++i) {
      double xv = (double)i;
      fvec[i] = y[i] - (x[0] * xv * xv * xv * xv + x[1] * xv * xv * xv + x[2]);
    }
    return 0;
  };

  info = levmar::lmdif1(fcn, y, m, n, x, fvec, DBL_EPSILON, iwa, wa, lwa);

  std::cout << "exit state: " << levmar::getInfo(info) << std::endl;
  std::cout << "final approximate solution: x0=" << (double)x[0]
            << ", x1=" << (double)x[1] << ", x2=" << (double)x[2] << std::endl;
  return 0;
}
