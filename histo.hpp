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

#ifndef HISTO_HPP
#define HISTO_HPP

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <vector>

#include "cubicSpline.hpp"

/**
@file histo.hpp
EXAMPLE:
@code{.cpp}
int main (int argc, char const *argv[])
{
        std::vector<double> v;
        //for (size_t i = 0 ; i < 1000 ; i++) v.push_back(rand()/(double)RAND_MAX);
        for (size_t i = 0 ; i < 1000 ; i++) v.push_back(i/1000.);

        //histo H = histo::histoNumBins(v, 20, false);
        //histo H = histo::pdfNumBins(v, 4);
        histo H = histo::pdfMaxPerBin(v, 200);
        for (size_t i = 0 ; i < H.data.size() ; i++) std::cout << H.data[i].X << " " << H.data[i].ProbDensity << " " <<
H.data[i].Width << std::endl;

        return 0;
}
@endcode
*/

// +-----+  ^
// |     |  |
// |<-W->|  | P
// |     |  |
// +--+--+  v
//    X
struct bar {
  double X;           // mean position of the bar
  double ProbDensity; // height of the bar
  double Width;       // width of the bar
};

class histo {
public:
  double min;
  double max;
  std::vector<bar> data;

  // Ctors
  histo(size_t nbins) { data.resize(nbins); }
  histo() {}
  void clear() { data.clear(); }

  double entropyShannon() {
    double S = 0.0;
    for (size_t i = 0; i < data.size(); i++) {
      if (data[i].ProbDensity > 0.0)
        S += data[i].ProbDensity * log10(data[i].ProbDensity);
    }
    return -S;
  }

  double entropyFisher() {
    double S = 0.0;
    for (size_t i = 0; i < data.size(); i++) {
      // TODO
    }
    return S;
  }

  void save(const char *name) {
    std::ofstream out(name);
    for (size_t i = 0; i < data.size(); i++) {
      if (isfinite(data[i].ProbDensity)) {
        out << data[i].X << ' ' << data[i].ProbDensity << ' ' << 0.5 * data[i].Width << '\n';
      }
    }
  }

  // histogram where the number of bins have the same width.
  // The histogram can be normalized (number per bin / total number)
  static histo histoNumBins(std::vector<double> &value, int nbins, bool normalize = true) {
    histo H(nbins);
    size_t nb = value.size();
    std::sort(value.begin(), value.end());
    double fact;
    if (normalize == true && nb > 0) {
      fact = 1.0 / nb;
    } else {
      fact = 1.0;
    }
    H.min = value[0];
    H.max = value[nb - 1];
    double binWidth = (H.max - H.min) / (double)nbins;
    double binAbscise;
    double threshold;
    size_t count = 0;
    for (int b = 0; b < nbins; b++) {
      binAbscise = H.min + (b + 0.5) * binWidth;
      size_t amount = 0;
      threshold = binAbscise + 0.5 * binWidth;
      while (count < nb && value[count] <= threshold) {
        ++amount;
        ++count;
      }
      H.data[b].Width = binWidth;
      H.data[b].X = binAbscise;
      H.data[b].ProbDensity = amount * fact;
    }
    return H;
  }

  // The pdf is normalized (integration from -infty to infty = 1)
  // and the width of all bins is (max - min) / nbins
  static histo pdfNumBins(std::vector<double> &value, int nbins) {
    histo H(nbins);
    std::sort(value.begin(), value.end());
    H.min = value[0];
    H.max = value[value.size() - 1];
    double binWidth = (H.max - H.min) / (double)nbins;
    double binAbscise;
    size_t amount;
    double threshold;
    size_t count = 0;
    for (int b = 0; b < nbins; b++) {
      binAbscise = H.min + (b + 0.5) * binWidth;
      amount = 0;
      threshold = binAbscise + 0.5 * binWidth;
      while (count < value.size() && value[count] <= threshold) {
        ++amount;
        ++count;
      }
      H.data[b].Width = binWidth;
      H.data[b].X = binAbscise;
      H.data[b].ProbDensity = amount / binWidth; // a density (not yet normalized)
    }

    // normalization: \int P dx = 1
    double sum = 0.0;
    for (int b = 0; b < nbins; b++) {
      if (H.data[b].Width > 0.0)
        sum += H.data[b].ProbDensity * H.data[b].Width;
    }
    double invSum = 1.0;
    if (sum > 0.0)
      invSum = 1.0 / sum;
    for (int b = 0; b < nbins; b++)
      H.data[b].ProbDensity *= invSum;

    return H;
  }

  // EXPERIMENTAL !!!!
  // remarque : une meilleur solution serait de fixer nbins à une valeur pas trop grande pour limiter les "vagues"
  // et de faire une optimisation des positions des quelques points de la spline,
  // puis on dérive une fois pour obtenir le pdf
  static void pdfNumBinsSpline(std::vector<double> &value, std::vector<double> &xs, std::vector<double> &ys,
                               int nbSlices) {
    std::sort(value.begin(), value.end());
    std::vector<double> x, y;
    // double deltax = (value[value.size() - 1] - value[0]) / (double)nbSlices;
    double deltax = (value.back() - value.front()) / (double)nbSlices;
    double deltap = 1.0 / (double)nbSlices;
    double pval0 = 0.0;
    double xval0 = value[0];
    for (size_t i = 0; i < value.size(); i++) {
      double pval = (double)i / (double)value.size();
      double xval = value[i];
      if ((pval - pval0) >= deltap || (xval - xval0) >= deltax || i == 0 || i == value.size() - 1) {
        x.push_back(xval);
        y.push_back(pval);
        pval0 = pval;
        xval0 = xval;
      }
    }

    std::vector<SplineSet> cs = spline(x, y);
    std::vector<SplineSet> deriv;
    derivSpline(cs, deriv);
    getSlineCurve(deriv, xs, ys, 10);
  }

  /*
  static histo pdfSpline(std::vector<double> &value, int nbPoints = 30) {
    std::sort(value.begin(), value.end());
    std::vector<double> xcdf, ycdf;
    //double deltax = (value.back() - value.front()) / (double)nbSlices;
    for (size_t i = 0; i < value.size(); i++) {
      double pval = i / (double)value.size();
      double xval = value[i];
      xcdf.push_back(xval);
      ycdf.push_back(pval);
    }

    std::vector<SplineSet> cs = spline(xcdf, ycdf);
    std::vector<SplineSet> deriv;
    derivSpline(cs, deriv);
    std::vector<double> xs, ys;
    getSlineCurve(deriv, xs, ys, 10);
  }
  */

  static histo pdfNumBinsRange(std::vector<double> &value, int nbins, double min, double max) {
    histo H(nbins);
    std::sort(value.begin(), value.end());
    H.min = min;
    H.max = max;
    double binWidth = (H.max - H.min) / (double)nbins;
    double binAbscise;
    size_t amount;
    double threshold;
    size_t count = 0;
    for (int b = 0; b < nbins; b++) {
      binAbscise = H.min + (b + 0.5) * binWidth;
      amount = 0;
      threshold = binAbscise + 0.5 * binWidth;
      while (count < value.size() && value[count] <= threshold) {
        ++amount;
        ++count;
      }
      H.data[b].Width = binWidth;
      H.data[b].X = binAbscise;
      H.data[b].ProbDensity = amount / binWidth; // a density (not yet normalized)
    }

    // normalization: \int P dx = 1
    double sum = 0.0;
    for (int b = 0; b < nbins; b++) {
      if (H.data[b].Width > 0.0)
        sum += H.data[b].ProbDensity * H.data[b].Width;
    }

    double invSum = 1.0;
    if (sum > 0.0) {
      invSum = 1.0 / sum;
    }
    for (int b = 0; b < nbins; b++)
      H.data[b].ProbDensity *= invSum;

    return H;
  }

  static histo pdfMaxPerBin(std::vector<double> &value, int maxEltPerBin) {
    histo H;
    std::sort(value.begin(), value.end());
    H.min = value[0];
    H.max = value[value.size() - 1];
    double x0 = H.min;
    double x1 = H.min;
    size_t ivalue = 0;
    size_t prev_ivalue = 0, amount = 0;
    bar B;
    while (ivalue < value.size()) {
      ivalue = prev_ivalue + maxEltPerBin;
      if (ivalue >= value.size()) {
        ivalue = value.size() - 1;
      }
      x1 = value[ivalue];
      amount = ivalue - prev_ivalue;
      if (amount == 0)
        break;
      B.Width = x1 - x0;
      B.X = 0.5 * (x0 + x1);
      if (B.Width > 0.0)
        B.ProbDensity = (double)amount / B.Width; // a density (not yet normalized)
      else
        B.ProbDensity = NAN;
      H.data.push_back(B);
      x0 = x1;
      prev_ivalue = ivalue;
    }

    // normalization: \int P dx = 1
    double sum = 0.0;
    for (size_t b = 0; b < H.data.size(); b++) {
      if (H.data[b].Width > 0.0)
        sum += H.data[b].ProbDensity * H.data[b].Width;
    }
    double invSum = 1.0;
    if (sum > 0.0)
      invSum = 1.0 / sum;
    for (size_t b = 0; b < H.data.size(); b++)
      H.data[b].ProbDensity *= invSum;

    return H;
  }

  static histo pdfMaxPerBin_minWidth(std::vector<double> &value, int maxEltPerBin, double minWidth) {
    histo H;
    std::sort(value.begin(), value.end());
    H.min = value[0];
    H.max = value[value.size() - 1];
    double x0 = H.min;
    double x1 = H.min;
    size_t ivalue = 0;
    size_t prev_ivalue = 0, amount = 0;
    bar B;

    while (ivalue < value.size()) {
      double width = 0.0;
      bool shouldBreak = false;
      while ((amount < maxEltPerBin || width < minWidth)) {
        ++ivalue;
        if (ivalue >= value.size()) {
          ivalue = value.size() - 1;
          shouldBreak = true;
        }
        x1 = value[ivalue];
        amount = ivalue - prev_ivalue;
        width = x1 - x0;
        if (shouldBreak == true)
          break;
      }
      if (amount == 0)
        break;
      B.Width = x1 - x0;
      B.X = 0.5 * (x0 + x1);
      if (B.Width > 0.0)
        B.ProbDensity = (double)amount / B.Width; // a density (not yet normalized)
      else
        B.ProbDensity = 0.0;
      H.data.push_back(B);
      x0 = x1;
      prev_ivalue = ivalue;
    }

    // normalization: \int P dx = 1
    double sum = 0.0;
    for (size_t b = 0; b < H.data.size(); b++) {
      if (H.data[b].Width > 0.0)
        sum += H.data[b].ProbDensity * H.data[b].Width;
    }
    double invSum = 1.0;
    if (sum > 0.0)
      invSum = 1.0 / sum;
    for (size_t b = 0; b < H.data.size(); b++)
      H.data[b].ProbDensity *= invSum;

    return H;
  }

  static histo pdf(std::vector<double> &value, int quality = 100) {
    histo H;
    std::sort(value.begin(), value.end());
    H.min = value[0];
    H.max = value[value.size() - 1];

    int maxEltPerBin = value.size() / quality;
    if (maxEltPerBin < 5)
      maxEltPerBin = 5;
    double minWidth = (H.max - H.min) / (double)quality;
    if (minWidth <= 0.0)
      minWidth = 0.1;

    double x0 = H.min;
    double x1 = H.min;
    size_t ivalue = 0;
    size_t prev_ivalue = 0, amount = 0;
    bar B;

    while (ivalue < value.size()) {
      double width = 0.0;
      bool shouldBreak = false;
      while (amount < maxEltPerBin || width < minWidth) {
        ++ivalue;
        if (ivalue >= value.size()) {
          ivalue = value.size() - 1;
          shouldBreak = true;
        }
        x1 = value[ivalue];
        amount = ivalue - prev_ivalue;
        width = x1 - x0;
        if (shouldBreak == true)
          break;
      }
      if (amount == 0)
        break;
      B.Width = x1 - x0;
      B.X = 0.5 * (x0 + x1);
      if (B.Width > 0.0)
        B.ProbDensity = (double)amount / B.Width; // a density (not yet normalized)
      else
        B.ProbDensity = 0.0;
      H.data.push_back(B);
      x0 = x1;
      prev_ivalue = ivalue;
    }

    // normalization: \int P dx = 1
    double sum = 0.0;
    for (size_t b = 0; b < H.data.size(); b++) {
      if (H.data[b].Width > 0.0)
        sum += H.data[b].ProbDensity * H.data[b].Width;
    }
    double invSum = 1.0;
    if (sum > 0.0)
      invSum = 1.0 / sum;
    for (size_t b = 0; b < H.data.size(); b++)
      H.data[b].ProbDensity *= invSum;

    return H;
  }
};

#endif /* end of include guard: HISTO_HPP */
