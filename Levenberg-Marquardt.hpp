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

#pragma once

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <functional>
#include <string>

class levmar {
public:
  static double enorm(int n, const double *x) {
    double ret_val, d1;

    int i;
    double s1, s2, s3, xabs, x1max, x3max, agiant;

    s1 = 0.0;
    s2 = 0.0;
    s3 = 0.0;
    x1max = 0.0;
    x3max = 0.0;
    agiant = DBL_MAX / (double)n;
    for (i = 0; i < n; ++i) {
      xabs = fabs(x[i]);
      if (xabs >= agiant) {
        // sum for large components
        if (xabs > x1max) {
          // Computing 2nd power
          d1 = x1max / xabs;
          s1 = 1.0 + s1 * (d1 * d1);
          x1max = xabs;
        } else {
          // Computing 2nd power
          d1 = xabs / x1max;
          s1 += d1 * d1;
        }
      } else if (xabs <= DBL_MIN) {
        // sum for small components
        if (xabs > x3max) {
          // Computing 2nd power
          d1 = x3max / xabs;
          s3 = 1.0 + s3 * (d1 * d1);
          x3max = xabs;
        } else if (xabs != 0.0) {
          // Computing 2nd power
          d1 = xabs / x3max;
          s3 += d1 * d1;
        }
      } else {
        // sum for intermediate components
        // Computing 2nd power
        s2 += xabs * xabs;
      }
    }

    // calculation of norm

    if (s1 != 0.0) {
      ret_val = x1max * sqrt(s1 + (s2 / x1max) / x1max);
    } else if (s2 != 0.0) {
      if (s2 >= x3max) {
        ret_val = sqrt(s2 * (1.0 + (x3max / s2) * (x3max * s3)));
      } else {
        ret_val = sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
      }
    } else {
      ret_val = x3max * sqrt(s3);
    }
    return ret_val;
  }

  static int fdjac2(std::function<int(void *, int, int, const double *, double *, int)> fcn_mn, void *p, int m, int n,
                    double *x, const double *fvec, double *fjac, int ldfjac, double epsfcn, double *wa) {
    double h;
    int i, j;
    double eps, temp, epsmch;
    int iflag;

    // epsmch is the machine precision
    epsmch = DBL_EPSILON;

    eps = sqrt((std::max(epsfcn, epsmch)));
    for (j = 0; j < n; ++j) {
      temp = x[j];
      h = eps * fabs(temp);
      if (h == 0.0) {
        h = eps;
      }
      x[j] = temp + h;
      // the last parameter of fcn_mn() is set to 2 to differentiate
      // calls made to compute the function from calls made to compute
      // the Jacobian (see fcn() in examples/lmfdrv.c, and how njev
      // is used to compute the number of Jacobian evaluations)
      iflag = fcn_mn(p, m, n, x, wa, 2);
      if (iflag < 0) {
        return iflag;
      }
      x[j] = temp;
      for (i = 0; i < m; ++i) {
        fjac[i + j * ldfjac] = (wa[i] - fvec[i]) / h;
      }
    }

    return 0;
  }

  static void qrfac(int m, int n, double *a, int lda, int pivot, int *ipvt, int lipvt, double *rdiag, double *acnorm,
                    double *wa) {

    double d1;

    int i, j, k, jp1;
    double sum;
    double temp;
    int minmn;
    double epsmch;
    double ajnorm;

    (void)lipvt;

    // epsmch is the machine precision
    epsmch = DBL_EPSILON;

    // compute the initial column norms and initialize several arrays

    for (j = 0; j < n; ++j) {
      acnorm[j] = enorm(m, &a[j * lda + 0]);
      rdiag[j] = acnorm[j];
      wa[j] = rdiag[j];
      if (pivot) {
        ipvt[j] = j + 1;
      }
    }

    // reduce a to r with householder transformations

    minmn = std::min(m, n);
    for (j = 0; j < minmn; ++j) {
      if (pivot) {

        // bring the column of largest norm into the pivot position

        int kmax = j;
        for (k = j; k < n; ++k) {
          if (rdiag[k] > rdiag[kmax]) {
            kmax = k;
          }
        }
        if (kmax != j) {
          for (i = 0; i < m; ++i) {
            temp = a[i + j * lda];
            a[i + j * lda] = a[i + kmax * lda];
            a[i + kmax * lda] = temp;
          }
          rdiag[kmax] = rdiag[j];
          wa[kmax] = wa[j];
          k = ipvt[j];
          ipvt[j] = ipvt[kmax];
          ipvt[kmax] = k;
        }
      }

      // compute the householder transformation to reduce the
      // j-th column of a to a multiple of the j-th unit vector

      ajnorm = enorm(m - (j + 1) + 1, &a[j + j * lda]);
      if (ajnorm != 0.) {
        if (a[j + j * lda] < 0.0) {
          ajnorm = -ajnorm;
        }
        for (i = j; i < m; ++i) {
          a[i + j * lda] /= ajnorm;
        }
        a[j + j * lda] += 1;

        // apply the transformation to the remaining columns
        // and update the norms

        jp1 = j + 1;
        if (n > jp1) {
          for (k = jp1; k < n; ++k) {
            sum = 0.0;
            for (i = j; i < m; ++i) {
              sum += a[i + j * lda] * a[i + k * lda];
            }
            temp = sum / a[j + j * lda];
            for (i = j; i < m; ++i) {
              a[i + k * lda] -= temp * a[i + j * lda];
            }
            if (pivot && rdiag[k] != 0.0) {
              temp = a[j + k * lda] / rdiag[k];
              // Computing MAX
              d1 = 1 - temp * temp;
              rdiag[k] *= sqrt((std::max((double)0.0, d1)));
              // Computing 2nd power
              d1 = rdiag[k] / wa[k];
              if (0.05 * (d1 * d1) <= epsmch) {
                rdiag[k] = enorm(m - (j + 1), &a[jp1 + k * lda]);
                wa[k] = rdiag[k];
              }
            }
          }
        }
      }
      rdiag[j] = -ajnorm;
    }
  }

  static void qrsolv(int n, double *r, int ldr, const int *ipvt, const double *diag, const double *qtb, double *x,
                     double *sdiag, double *wa) {

    int i, j, k, l;
    double Cos, Sin, sum, temp;
    int nsing;
    double qtbpj;

    // copy r and (q transpose)*b to preserve input and initialize s
    // in particular, save the diagonal elements of r in x

    for (j = 0; j < n; ++j) {
      for (i = j; i < n; ++i) {
        r[i + j * ldr] = r[j + i * ldr];
      }
      x[j] = r[j + j * ldr];
      wa[j] = qtb[j];
    }

    // eliminate the diagonal matrix d using a givens rotation

    for (j = 0; j < n; ++j) {

      // prepare the row of d to be eliminated, locating the
      // diagonal element using p from the qr factorization

      l = ipvt[j] - 1;
      if (diag[l] != 0.0) {
        for (k = j; k < n; ++k) {
          sdiag[k] = 0.0;
        }
        sdiag[j] = diag[l];

        // the transformations to eliminate the row of d
        // modify only a single element of (q transpose)*b
        // beyond the first n, which is initially zero

        qtbpj = 0.0;
        for (k = j; k < n; ++k) {

          // determine a givens rotation which eliminates the
          // appropriate element in the current row of d

          if (sdiag[k] != 0.0) {

            if (fabs(r[k + k * ldr]) < fabs(sdiag[k])) {
              double cotan;
              cotan = r[k + k * ldr] / sdiag[k];
              Sin = 0.5 / sqrt(0.25 + 0.25 * (cotan * cotan));
              Cos = Sin * cotan;
            } else {
              double Tan;
              Tan = sdiag[k] / r[k + k * ldr];
              Cos = 0.5 / sqrt(0.25 + 0.25 * (Tan * Tan));
              Sin = Cos * Tan;
            }

            // compute the modified diagonal element of r and
            // the modified element of ((q transpose)*b,0)

            temp = Cos * wa[k] + Sin * qtbpj;
            qtbpj = -Sin * wa[k] + Cos * qtbpj;
            wa[k] = temp;

            // accumulate the tranformation in the row of s

            r[k + k * ldr] = Cos * r[k + k * ldr] + Sin * sdiag[k];
            if (n > k + 1) {
              for (i = k + 1; i < n; ++i) {
                temp = Cos * r[i + k * ldr] + Sin * sdiag[i];
                sdiag[i] = -Sin * r[i + k * ldr] + Cos * sdiag[i];
                r[i + k * ldr] = temp;
              }
            }
          }
        }
      }

      // store the diagonal element of s and restore
      // the corresponding diagonal element of r

      sdiag[j] = r[j + j * ldr];
      r[j + j * ldr] = x[j];
    }

    // solve the triangular system for z. if the system is
    // singular, then obtain a least squares solution

    nsing = n;
    for (j = 0; j < n; ++j) {
      if (sdiag[j] == 0.0 && nsing == n) {
        nsing = j;
      }
      if (nsing < n) {
        wa[j] = 0.0;
      }
    }
    if (nsing >= 1) {
      for (k = 1; k <= nsing; ++k) {
        j = nsing - k;
        sum = 0.0;
        if (nsing > j + 1) {
          for (i = j + 1; i < nsing; ++i) {
            sum += r[i + j * ldr] * wa[i];
          }
        }
        wa[j] = (wa[j] - sum) / sdiag[j];
      }
    }

    // permute the components of z back to components of x

    for (j = 0; j < n; ++j) {
      l = ipvt[j] - 1;
      x[l] = wa[j];
    }
    return;
  }

  static void lmpar(int n, double *r, int ldr, const int *ipvt, const double *diag, const double *qtb, double delta,
                    double *par, double *x, double *sdiag, double *wa1, double *wa2) {

    double d1, d2;

    int j, l;
    double fp;
    double parc, parl;
    int iter;
    double temp, paru, dwarf;
    int nsing;
    double gnorm;
    double dxnorm;

    // dwarf is the smallest positive magnitude
    dwarf = DBL_MIN;

    // compute and store in x the gauss-newton direction. if the
    // jacobian is rank-deficient, obtain a least squares solution

    nsing = n;
    for (j = 0; j < n; ++j) {
      wa1[j] = qtb[j];
      if (r[j + j * ldr] == 0.0 && nsing == n) {
        nsing = j;
      }
      if (nsing < n) {
        wa1[j] = 0.0;
      }
    }

    if (nsing >= 1) {
      int k;
      for (k = 1; k <= nsing; ++k) {
        j = nsing - k;
        wa1[j] /= r[j + j * ldr];
        temp = wa1[j];
        if (j >= 1) {
          int i;
          for (i = 0; i < j; ++i) {
            wa1[i] -= r[i + j * ldr] * temp;
          }
        }
      }
    }

    for (j = 0; j < n; ++j) {
      l = ipvt[j] - 1;
      x[l] = wa1[j];
    }

    // initialize the iteration counter
    // evaluate the function at the origin, and test
    // for acceptance of the gauss-newton direction

    iter = 0;
    for (j = 0; j < n; ++j) {
      wa2[j] = diag[j] * x[j];
    }
    dxnorm = enorm(n, wa2);
    fp = dxnorm - delta;
    if (fp <= 0.1 * delta) {
      goto TERMINATE;
    }

    // if the jacobian is not rank deficient, the newton
    // step provides a lower bound, parl, for the zero of
    // the function. otherwise set this bound to zero

    parl = 0.0;
    if (nsing >= n) {
      for (j = 0; j < n; ++j) {
        l = ipvt[j] - 1;
        wa1[j] = diag[l] * (wa2[l] / dxnorm);
      }

      for (j = 0; j < n; ++j) {
        double sum = 0.0;
        if (j >= 1) {
          int i;
          for (i = 0; i < j; ++i) {
            sum += r[i + j * ldr] * wa1[i];
          }
        }
        wa1[j] = (wa1[j] - sum) / r[j + j * ldr];
      }

      temp = enorm(n, wa1);
      parl = fp / delta / temp / temp;
    }

    // calculate an upper bound, paru, for the zero of the function

    for (j = 0; j < n; ++j) {
      double sum;

      int i;
      sum = 0.0;
      for (i = 0; i <= j; ++i) {
        sum += r[i + j * ldr] * qtb[i];
      }

      l = ipvt[j] - 1;
      wa1[j] = sum / diag[l];
    }
    gnorm = enorm(n, wa1);
    paru = gnorm / delta;
    if (paru == 0.0) {
      paru = dwarf / std::min(delta, (double)0.1);
    }

    // if the input par lies outside of the interval (parl,paru),
    // set par to the closer endpoint

    *par = std::max(*par, parl);
    *par = std::min(*par, paru);
    if (*par == 0.0) {
      *par = gnorm / dxnorm;
    }

    // beginning of an iteration

    for (;;) {
      ++iter;

      // evaluate the function at the current value of par

      if (*par == 0.0) {
        // Computing MAX
        d1 = dwarf, d2 = 0.001 * paru;
        *par = std::max(d1, d2);
      }
      temp = sqrt(*par);
      for (j = 0; j < n; ++j) {
        wa1[j] = temp * diag[j];
      }
      qrsolv(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2);
      for (j = 0; j < n; ++j) {
        wa2[j] = diag[j] * x[j];
      }
      dxnorm = enorm(n, wa2);
      temp = fp;
      fp = dxnorm - delta;

      // if the function is small enough, accept the current value
      // of par. also test for the exceptional cases where parl
      // is zero or the number of iterations has reached 10

      if (fabs(fp) <= 0.1 * delta || (parl == 0.0 && fp <= temp && temp < 0.0) || iter == 10) {
        goto TERMINATE;
      }

      // compute the newton correction

      for (j = 0; j < n; ++j) {
        l = ipvt[j] - 1;
        wa1[j] = diag[l] * (wa2[l] / dxnorm);
      }
      for (j = 0; j < n; ++j) {
        wa1[j] /= sdiag[j];
        temp = wa1[j];
        if (n > j + 1) {
          int i;
          for (i = j + 1; i < n; ++i) {
            wa1[i] -= r[i + j * ldr] * temp;
          }
        }
      }

      temp = enorm(n, wa1);
      parc = fp / delta / temp / temp;

      // depending on the sign of the function, update parl or paru

      if (fp > 0.) {
        parl = std::max(parl, *par);
      }
      if (fp < 0.) {
        paru = std::min(paru, *par);
      }

      // compute an improved estimate for par
      d1 = parl, d2 = *par + parc;
      *par = std::max(d1, d2);

      // end of an iteration
    }

  TERMINATE:

    if (iter == 0) {
      *par = 0.0;
    }
  }

  static int lmdif(std::function<int(void *, int, int, const double *, double *, int)> fcn_mn, void *p, int m, int n,
                   double *x, double *fvec, double ftol, double xtol, double gtol, int maxfev, double epsfcn,
                   double *diag, int mode, double factor, int nprint, int *nfev, double *fjac, int ldfjac, int *ipvt,
                   double *qtf, double *wa1, double *wa2, double *wa3, double *wa4) {

    // System generated locals
    double d1, d2;

    // Local variables
    int i, j, l;
    double par, sum;
    int iter;
    double temp, temp1, temp2;
    int iflag;
    double delta = 0.0;
    double ratio;
    double fnorm, gnorm;
    double pnorm, xnorm = 0.0, fnorm1, actred, dirder, epsmch, prered;
    int info;

    epsmch = DBL_EPSILON;

    info = 0;
    iflag = 0;
    *nfev = 0;

    // check the input parameters for errors

    if (n <= 0 || m < n || ldfjac < m || ftol < 0.0 || xtol < 0.0 || gtol < 0.0 || maxfev <= 0 || factor <= 0.0) {
      goto TERMINATE;
    }
    if (mode == 2) {
      for (j = 0; j < n; ++j) {
        if (diag[j] <= 0.0) {
          goto TERMINATE;
        }
      }
    }

    // evaluate the function at the starting point
    // and calculate its norm

    iflag = fcn_mn(p, m, n, x, fvec, 1);
    *nfev = 1;
    if (iflag < 0) {
      goto TERMINATE;
    }
    fnorm = enorm(m, fvec);

    // initialize levenberg-marquardt parameter and iteration counter

    par = 0.0;
    iter = 1;

    // beginning of the outer loop

    for (;;) {

      // calculate the jacobian matrix

      iflag = fdjac2(fcn_mn, p, m, n, x, fvec, fjac, ldfjac, epsfcn, wa4);
      *nfev += n;
      if (iflag < 0) {
        goto TERMINATE;
      }

      // if requested, call fcn to enable printing of iterates

      if (nprint > 0) {
        iflag = 0;
        if ((iter - 1) % nprint == 0) {
          iflag = fcn_mn(p, m, n, x, fvec, 0);
        }
        if (iflag < 0) {
          goto TERMINATE;
        }
      }

      // compute the qr factorization of the jacobian

      qrfac(m, n, fjac, ldfjac, 1, ipvt, n, wa1, wa2, wa3);

      // on the first iteration and if mode is 1, scale according
      // to the norms of the columns of the initial jacobian

      if (iter == 1) {
        if (mode != 2) {
          for (j = 0; j < n; ++j) {
            diag[j] = wa2[j];
            if (wa2[j] == 0.0) {
              diag[j] = 1.0;
            }
          }
        }

        // on the first iteration, calculate the norm of the scaled x
        // and initialize the step bound delta

        for (j = 0; j < n; ++j) {
          wa3[j] = diag[j] * x[j];
        }
        xnorm = enorm(n, wa3);
        delta = factor * xnorm;
        if (delta == 0.0) {
          delta = factor;
        }
      }

      // form (q transpose)*fvec and store the first n components in qtf

      for (i = 0; i < m; ++i) {
        wa4[i] = fvec[i];
      }
      for (j = 0; j < n; ++j) {
        if (fjac[j + j * ldfjac] != 0.0) {
          sum = 0.0;
          for (i = j; i < m; ++i) {
            sum += fjac[i + j * ldfjac] * wa4[i];
          }
          temp = -sum / fjac[j + j * ldfjac];
          for (i = j; i < m; ++i) {
            wa4[i] += fjac[i + j * ldfjac] * temp;
          }
        }
        fjac[j + j * ldfjac] = wa1[j];
        qtf[j] = wa4[j];
      }

      // compute the norm of the scaled gradient

      gnorm = 0.0;
      if (fnorm != 0.0) {
        for (j = 0; j < n; ++j) {
          l = ipvt[j] - 1;
          if (wa2[l] != 0.0) {
            sum = 0.0;
            for (i = 0; i <= j; ++i) {
              sum += fjac[i + j * ldfjac] * (qtf[i] / fnorm);
            }
            // Computing MAX
            d1 = fabs(sum / wa2[l]);
            gnorm = std::max(gnorm, d1);
          }
        }
      }

      // test for convergence of the gradient norm

      if (gnorm <= gtol) {
        info = 4;
      }
      if (info != 0) {
        goto TERMINATE;
      }

      // rescale if necessary

      if (mode != 2) {
        for (j = 0; j < n; ++j) {
          // Computing MAX
          d1 = diag[j], d2 = wa2[j];
          diag[j] = std::max(d1, d2);
        }
      }

      // beginning of the inner loop

      do {

        // determine the levenberg-marquardt parameter
        lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, &par, wa1, wa2, wa3, wa4);

        // store the direction p and x + p. calculate the norm of p
        for (j = 0; j < n; ++j) {
          wa1[j] = -wa1[j];
          wa2[j] = x[j] + wa1[j];
          wa3[j] = diag[j] * wa1[j];
        }
        pnorm = enorm(n, wa3);

        // on the first iteration, adjust the initial step bound

        if (iter == 1) {
          delta = std::min(delta, pnorm);
        }

        // evaluate the function at x + p and calculate its norm

        iflag = fcn_mn(p, m, n, wa2, wa4, 1);
        ++(*nfev);
        if (iflag < 0) {
          goto TERMINATE;
        }
        fnorm1 = enorm(m, wa4);

        // compute the scaled actual reduction

        actred = -1.0;
        if (0.1 * fnorm1 < fnorm) {
          // Computing 2nd power
          d1 = fnorm1 / fnorm;
          actred = 1 - d1 * d1;
        }

        // compute the scaled predicted reduction and
        // the scaled directional derivative

        for (j = 0; j < n; ++j) {
          wa3[j] = 0.0;
          l = ipvt[j] - 1;
          temp = wa1[l];
          for (i = 0; i <= j; ++i) {
            wa3[i] += fjac[i + j * ldfjac] * temp;
          }
        }
        temp1 = enorm(n, wa3) / fnorm;
        temp2 = (sqrt(par) * pnorm) / fnorm;
        prered = temp1 * temp1 + temp2 * temp2 / 0.5;
        dirder = -(temp1 * temp1 + temp2 * temp2);

        // compute the ratio of the actual to the predicted reduction

        ratio = 0.0;
        if (prered != 0.0) {
          ratio = actred / prered;
        }

        // update the step bound

        if (ratio <= 0.25) {
          if (actred >= 0.0) {
            temp = 0.5;
          } else {
            temp = 0.5 * dirder / (dirder + 0.5 * actred);
          }
          if (0.1 * fnorm1 >= fnorm || temp < 0.1) {
            temp = 0.1;
          }
          // Computing MIN
          d1 = pnorm / 0.1;
          delta = temp * std::min(delta, d1);
          par /= temp;
        } else {
          if (par == 0.0 || ratio >= 0.75) {
            delta = pnorm / 0.5;
            par = 0.5 * par;
          }
        }

        // test for successful iteration

        if (ratio >= 1e-4) {

          // successful iteration. update x, fvec, and their norms

          for (j = 0; j < n; ++j) {
            x[j] = wa2[j];
            wa2[j] = diag[j] * x[j];
          }
          for (i = 0; i < m; ++i) {
            fvec[i] = wa4[i];
          }
          xnorm = enorm(n, wa2);
          fnorm = fnorm1;
          ++iter;
        }

        // tests for convergence

        if (fabs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1.0) {
          info = 1;
        }
        if (delta <= xtol * xnorm) {
          info = 2;
        }
        if (fabs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1.0 && info == 2) {
          info = 3;
        }
        if (info != 0) {
          goto TERMINATE;
        }

        // tests for termination and stringent tolerances

        if (*nfev >= maxfev) {
          info = 5;
        }
        if (fabs(actred) <= epsmch && prered <= epsmch && 0.5 * ratio <= 1.0) {
          info = 6;
        }
        if (delta <= epsmch * xnorm) {
          info = 7;
        }
        if (gnorm <= epsmch) {
          info = 8;
        }
        if (info != 0) {
          goto TERMINATE;
        }

        // end of the inner loop. repeat if iteration unsuccessful

      } while (ratio < 1e-4);

      // end of the outer loop
    }

  TERMINATE:

    // termination, either normal or user imposed

    if (iflag < 0) {
      info = iflag;
    }
    if (nprint > 0) {
      fcn_mn(p, m, n, x, fvec, 0);
    }

    return info;
  }

  static int lmdif1(std::function<int(void *, int, int, const double *, double *, int)> fcn_mn, void *p, int m, int n,
                    double *x, double *fvec, double tol, int *iwa, double *wa, int lwa) {

    const double factor = 100.0;

    int mp5n, mode, nfev;
    double ftol, gtol, xtol;
    double epsfcn;
    int maxfev, nprint;
    int info;

    // check the input parameters for errors
    if (n <= 0 || m < n || tol < 0.0 || lwa < m * n + n * 5 + m) {
      return 0;
    }

    // call lmdif
    maxfev = (n + 1) * 200;
    ftol = tol;
    xtol = tol;
    gtol = 0.0;
    epsfcn = 0.0;
    mode = 1;
    nprint = 0;
    mp5n = m + n * 5;
    info = lmdif(fcn_mn, p, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, wa, mode, factor, nprint, &nfev, &wa[mp5n],
                 m, iwa, &wa[n], &wa[(n << 1)], &wa[n * 3], &wa[(n << 2)], &wa[n * 5]);
    if (info == 8) {
      info = 4;
    }
    return info;
  }

  static int get_lwa(int m, int n) { return m * n + 5 * n + m; }

  static std::string getInfo(int info) {
    std::string txt;

    switch (info) {
    case 0: {
      txt = "(info = 0) improper input parameters";
    } break;
    case 1: {
      txt = "(info = 1) algorithm estimates that the relative error in the sum of squares is at most tol";
    } break;
    case 2: {
      txt = "(info = 2) algorithm estimates that the relative error between x and the solution is at most tol";
    } break;
    case 3: {
      txt = "(info = 3) conditions for info = 1 and info = 2 both hold. algorithm estimates that the relative error "
            "(1) in the sum of squares is at most tol, (2) between x and the solution is at most tol";
    } break;
    case 4: {
      txt = "(info = 4) fvec is orthogonal to the columns of the jacobian to machine precision";
    } break;
    case 5: {
      txt = "(info = 5) number of calls to fcn has reached or exceeded 200*(n+1)";
    } break;
    case 6: {
      txt = "(info = 6) tol is too small. no further reduction in the sum of squares is possible";
    } break;
    case 7: {
      txt = "(info = 7) tol is too small. no further improvement in the approximate solution x is possible";
    } break;
    default: {
      txt = "this value of info is unknown";
    } break;
    }
    return txt;
  }
};
