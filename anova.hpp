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

#ifndef ANOVA_HPP
#define ANOVA_HPP

// =======================================================
// Simple One Way ANOVA (Analysis of Variance) in c++
// =======================================================

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

struct anovaResult {
  double ss_treatment;
  size_t df_treatment;
  double ms_treatment;
  double F;
  double p;

  double ss_residual;
  size_t df_residual;
  double ms_residual;

  double ss_total;
  size_t df_total;
};

class anova {
public:

  /// Compute the mean of all the samples
  /// @param samples 2D vector of samples
  /// @return the mean of all the samples
  static double fullMean(std::vector<std::vector<double>> &samples) {
    double total = 0;
    size_t count = 0;
    for (size_t i = 0; i < samples.size(); i++) {
      for (size_t j = 0; j < samples[i].size(); j++) {
        total += samples[i][j];
      }
      count += samples[i].size();
    }
    return (double)total / (double)count;
  }

  /**
   * @brief Compute the mean of a sample
   * @param sample The sample
   * @return The mean of the sample
   */
  static double mean(std::vector<double> &sample) {
    double total = 0;
    for (size_t i = 0; i < sample.size(); i++) {
      total += sample[i];
    }
    return (double)total / (double)sample.size();
  }

  // Compute the Sum of Squares
  static void SS(std::vector<std::vector<double>> &samples, double &total, double &treatment, double &residual) {
    treatment = 0.0;
    residual = 0.0;
    double fullMean = anova::fullMean(samples);
    for (size_t i = 0; i < samples.size(); i++) {
      double mu = anova::mean(samples[i]);
      for (size_t j = 0; j < samples[i].size(); j++) {
        residual += (samples[i][j] - mu) * (samples[i][j] - mu);
      }
      treatment += (mu - fullMean) * (mu - fullMean) * samples[i].size();
    }
    total = treatment + residual;
  }

  // Compute the Degrees of Freedom
  static void DF(std::vector<std::vector<double>> &samples, size_t &total, size_t &treatment, size_t &residual) {
    size_t tot = 0;
    for (size_t i = 0; i < samples.size(); i++) {
      tot += samples[i].size();
    }
    treatment = samples.size() - 1; // numerator
    residual = tot - samples.size();
    total = treatment + residual;
  }

/**
 * @brief Perform a one-way ANOVA test on the given samples.
 * 
 * This function calculates the sum of squares, degrees of freedom, mean squares,
 * F-statistic, and p-value for the provided samples.
 * 
 * @param samples A 2D vector containing the samples for analysis.
 * @param result An anovaResult struct to store the computed ANOVA results.
 */  
  static void test(std::vector<std::vector<double>> &samples, anovaResult &result) {
    double ss_total;
    double ss_treatment;
    double ss_residual;
    SS(samples, ss_total, ss_treatment, ss_residual);
    size_t df_total;
    size_t df_treatment;
    size_t df_residual;
    DF(samples, df_total, df_treatment, df_residual);

    result.ss_treatment = ss_treatment;
    result.df_treatment = df_treatment;
    result.ms_treatment = ss_treatment / (double)df_treatment;
    result.ss_residual = ss_residual;
    result.df_residual = df_residual;
    result.ms_residual = ss_residual / (double)df_residual;
    result.ss_total = ss_total;
    result.df_total = df_total;
    result.F = result.ms_treatment / result.ms_residual;
    result.p = FDist(result.F, result.df_treatment, result.df_residual);
  }

/**
 * @brief Prints the ANOVA results in a tabular format.
 * 
 * This function outputs the source of variance, degrees of freedom, 
 * sum of squares, mean square, F-statistic, and p-value for both 
 * between-group and within-group variances, as well as the total variance.
 * 
 * @param result An anovaResult struct containing the computed ANOVA results.
 */
  static void print(anovaResult &result) {
    std::printf("Source   DF    Sum of Square  Mean Square  F-stat       p-value\n");
    std::printf("Between  %-5zu %0.4f         %0.4f       %0.4f     %0.8f\n", result.df_treatment, result.ss_treatment,
                result.ms_treatment, result.F, result.p);
    std::printf("Within   %-5zu %0.4f         %0.4f\n", result.df_residual, result.ss_residual, result.ms_residual);
    std::printf("Total    %-5zu %0.4f         %0.4f\n", result.df_total, result.ss_total,
                result.ss_total / result.df_total);
  }

  /**
   * @brief Computes the cumulative distribution function of the F distribution.
   * 
   * @param F The F-statistic.
   * @param m The degrees of freedom in the numerator.
   * @param n The degrees of freedom in the denominator.
   * 
   * @return The cumulative distribution function of the F distribution with degrees
   * of freedom m and n at F.
   */
  static double FDist(double F, double m, double n) {
    double xx = 0.0, p = 0.0;

    if (m <= 0 || n <= 0)
      p = -1;
    else if (F > 0) {
      xx = F / (F + n / m);
      p = anova::betainc(xx, m / 2, n / 2);
    }
    return (1 - p);
  }

  /**
   * @brief Computes the cumulative distribution function of the beta distribution.
   * 
   * @param x A value in the range [0, 1].
   * @param a The first shape parameter of the beta distribution.
   * @param b The second shape parameter of the beta distribution.
   * 
   * @return The cumulative distribution function of the beta distribution evaluated at x.
   */
  static double betainc(double x, double a, double b) {
    double y, BT, AAA;

    if (x == 0 || x == 1)
      BT = 0;
    else {
      AAA = anova::gamma(a + b) - anova::gamma(a) - anova::gamma(b);
      BT = exp(AAA + a * log(x) + b * log(1 - x));
    }
    if (x < (a + 1) / (a + b + 2))
      y = BT * beta_cf(a, b, x) / a;
    else
      y = 1 - BT * beta_cf(b, a, 1 - x) / b;

    return y;
  }

  /**
   * @brief Computes the continued fraction representation of the incomplete beta
   *        function.
   * 
   * @param a The first shape parameter of the beta distribution.
   * @param b The second shape parameter of the beta distribution.
   * @param x A value in the range [0, 1].
   * 
   * @return The continued fraction representation of the incomplete beta
   *         function evaluated at x.
   */
  static double beta_cf(double a, double b, double x) {
    int count, count_max = 100;
    double eps = 0.0000001;
    double AM = 1;
    double BM = 1;
    double AZ = 1;
    double QAB;
    double QAP;
    double QAM;
    double BZ, EM, TEM, D, AP, BP, AAP, BPP, AOLD;

    QAB = a + b;
    QAP = a + 1;
    QAM = a - 1;
    BZ = 1 - QAB * x / QAP;

    for (count = 1; count <= count_max; count++) {
      EM = count;
      TEM = EM + EM;
      D = EM * (b - count) * x / ((QAM + TEM) * (a + TEM));
      AP = AZ + D * AM;
      BP = BZ + D * BM;
      D = -(a + EM) * (QAB + EM) * x / ((a + TEM) * (QAP + TEM));
      AAP = AP + D * AZ;
      BPP = BP + D * BZ;
      AOLD = AZ;
      AM = AP / BPP;
      BM = BP / BPP;
      AZ = AAP / BPP;
      BZ = 1;
      if (fabs(AZ - AOLD) < eps * fabs(AZ))
        return (AZ);
    }
    return AZ;
  }

/**
 * Computes the logarithm of the gamma function for a given input.
 *
 * This function uses the Lanczos approximation to compute the natural logarithm
 * of the gamma function for the input value `xx`. It is particularly useful for
 * large values of `xx`, where the gamma function itself would result in overflow.
 *
 * @param xx The input value for which to compute the logarithm of the gamma function.
 *           It should be a positive real number.
 * @return The natural logarithm of the gamma function evaluated at `xx`.
 */
  static double gamma(double xx) {
    double coef_const[7];
    double step = 2.50662827465;
    double HALF = 0.5;
    double ONE = 1;
    double FPF = 5.5;
    double SER, temp, x, y;
    int j;

    coef_const[1] = 76.18009173;
    coef_const[2] = -86.50532033;
    coef_const[3] = 24.01409822;
    coef_const[4] = -1.231739516;
    coef_const[5] = 0.00120858003;
    coef_const[6] = -0.00000536382;

    x = xx - ONE;
    temp = x + FPF;
    temp = (x + HALF) * log(temp) - temp;
    SER = ONE;
    for (j = 1; j <= 6; j++) {
      x = x + ONE;
      SER = SER + coef_const[j] / x;
    }
    y = temp + log(step * SER);

    return y;
  }
};

#endif /* end of include guard: ANOVA_HPP */

