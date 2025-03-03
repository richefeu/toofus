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

#ifndef GEOPACK2D_HPP
#define GEOPACK2D_HPP

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

class GeoPack2D {
public:
  struct Disk {
    double x, y;
    double r;
    Disk(double x_, double y_, double r_) {
      x = x_;
      y = y_;
      r = r_;
    }
  };

  double rmin;
  double rmax;
  double gapTol;
  double distNeighbor;
  double distMin;
  int max;
  int k;

  std::vector<Disk> sample;
  std::vector<std::vector<int>> prox;
  std::vector<int> boundaries;
  std::vector<int> active;

  double xmin, xmax;
  double ymin, ymax;

  // Ctor
  GeoPack2D() {
    rmin = rmax = 0.0;
    gapTol = distNeighbor = distMin = 0.0;
    max = 0;
    k = 0;
    xmin = xmax = ymin = ymax = 0.0;
  }

  /**
   * @brief Constructs a GeoPack2D object with specified parameters.
   *
   * @param rmin_ Minimum radius of disks.
   * @param rmax_ Maximum radius of disks.
   * @param k_ Parameter related to the number of disks or iterations.
   * @param xmin_ Minimum x-coordinate of the domain.
   * @param xmax_ Maximum x-coordinate of the domain.
   * @param ymin_ Minimum y-coordinate of the domain.
   * @param ymax_ Maximum y-coordinate of the domain.
   * @param gap_ (Optional) Gap tolerance between disks, default is 0.0.
   * @param max_ (Optional) Maximum number of iterations or disks, default is 0.
   */
  GeoPack2D(double rmin_, double rmax_, int k_, double xmin_, double xmax_, double ymin_, double ymax_,
            double gap_ = 0.0, int max_ = 0) {
    parameters(rmin_, rmax_, k_, xmin_, xmax_, ymin_, ymax_, gap_, max_);
  }

  /**
   * @brief Sets the parameters of the GeoPack2D instance.
   *
   * @param rmin_ Minimum radius of disks.
   * @param rmax_ Maximum radius of disks.
   * @param k_ Parameter related to the number of disks or iterations.
   * @param xmin_ Minimum x-coordinate of the domain.
   * @param xmax_ Maximum x-coordinate of the domain.
   * @param ymin_ Minimum y-coordinate of the domain.
   * @param ymax_ Maximum y-coordinate of the domain.
   * @param gap_ (Optional) Gap tolerance between disks, default is 0.0.
   * @param max_ (Optional) Maximum number of iterations or disks, default is 0.
   */
  void parameters(double rmin_, double rmax_, int k_, double xmin_, double xmax_, double ymin_, double ymax_,
                  double gapTol_ = 0.0, int max_ = 0) {
    rmin = rmin_;
    rmax = rmax_;
    gapTol = gapTol_;
    distNeighbor = 2 * rmax + gapTol;
    distMin = 0.0;
    k = k_;

    // Domain
    xmin = xmin_;
    xmax = xmax_;
    ymin = ymin_;
    ymax = ymax_;

    max = max_;
    if (max == 0) {
      int nx = (int)floor((xmax - xmin) / rmin);
      int ny = (int)floor((ymax - ymin) / rmin);
      max = nx * ny;
    }
  }
	
  
  /**
   * @brief Calculates the solid fraction of the 2D domain.
   *
   * @details This function computes the overall solid fraction by assuming 
   * there is no overlap between disks. It calculates the total area of the 
   * disks and divides it by the area of the domain. If the area of the 
   * domain is non-positive, it returns 0.0.
   *
   * @return The solid fraction as a double value, representing the ratio of 
   * the total area of the disks to the total area of the domain.
   */
  double getSolidFraction() {
    double Vtot = (xmax - xmin) * (ymax - ymin);
    if (Vtot <= 0.0)
      return 0.0;
    double Vs = 0.0;

    for (size_t i = 0; i < sample.size(); i++) {
      Vs += M_PI * sample[i].r * sample[i].r;
    }
    return Vs / Vtot;
  }

  /**
   * @brief Seeds the random number generator with the current time.
   *
   * This function calls srand() with the current time to seed the random
   * number generator. This is useful for generating different sequences of
   * random numbers.
   */
  void seedTime() 
  { srand(time(NULL)); }

  /**
   * @brief Adds a disk at given coordinates and radius.
   *
   * @details This function adds a disk at the specified coordinates and radius
   * to the sample and proximity data structures. The proximity data structure
   * is a vector of vectors, where each inner vector contains the indices of
   * disks that are within distNeighbor + distMin of the disk at the same index
   * in the sample vector. The active vector is also updated to include the
   * index of the newly added disk.
   *
   * @param x_ x-coordinate of the disk.
   * @param y_ y-coordinate of the disk.
   * @param r_ radius of the disk.
   */
  void appendDisk(double x_, double y_, double r_) {
    Disk D(x_, y_, r_);
    sample.push_back(D);
    prox.push_back(std::vector<int>());
    int particleIndex = sample.size() - 1;

    for (int i = 0; i < particleIndex; i++) {
      double dx = fabs(sample[i].x - x_);
      double dy = fabs(sample[i].y - y_);
      double d = sqrt(dx * dx + dy * dy);
      if (d < sample[i].r + r_ + distNeighbor + distMin) {
        prox[i].push_back(particleIndex);
        prox[particleIndex].push_back(i);
      }
    }

    active.push_back(particleIndex);
  }

  /**
   * @brief Reactivates a range of indices in the active list.
   *
   * @details This function appends indices from a specified range [from, to) 
   * into the active vector. If the 'to' parameter is not provided or is zero, 
   * it defaults to the size of the sample vector. This effectively reactivates 
   * all indices from the 'from' position to the end of the sample vector.
   *
   * @param from The starting index of the range to reactivate (inclusive).
   * @param to The ending index of the range to reactivate (exclusive). 
   * Defaults to the size of the sample vector if set to 0.
   */
  void reActivate(int from = 0, int to = 0) {
    if (to == 0)
      to = sample.size();
    for (int i = from; i < to; i++) {
      active.push_back(i);
    }
  }

  
  /**
   * @brief Executes the packing algorithm.
   *
   * @details This function starts the iterative process of adding disks
   * to the sample vector. It stops when the sample vector has reached the
   * maximum size or when no more disks can be added.
   *
   * @sa Packing parameters (e.g. k, max, gapTol, distMin, distNeighbor)
   */
  void exec() {
    // step 1
    if (sample.empty()) {
      double firstr = ran(rmin, rmax);
      Disk P(ran(xmin + firstr, xmax - firstr), ran(ymin + firstr, ymax - firstr), firstr);
      sample.push_back(P);
      prox.push_back(std::vector<int>());
      active.push_back(sample.size() - 1);
    }

    // step 2
    while (active.size() > 0 && sample.size() < (size_t)max) {
      int randIndex = rand() % active.size();
      int currentDisk = active[randIndex];
      double packedx = sample[currentDisk].x;
      double packedy = sample[currentDisk].y;
      double packedr = sample[currentDisk].r;

      bool found = false;

      for (int n = 0; n < k; n++) {
        double testr = ran(rmin, rmax);
        double angle = ran(0, 2.0 * M_PI);
        double m = ran(testr + packedr + distMin, testr + packedr + distMin + gapTol);
        double testx = packedx + m * cos(angle);
        double testy = packedy + m * sin(angle);

        bool ok = true;
        // boundaries
        if (testx < xmin + testr || testx > xmax - testr || testy < ymin + testr || testy > ymax - testr) {
          ok = false;
        }

        // inter-particles
        if (ok == true) {
          for (size_t i = 0; i < prox[currentDisk].size(); i++) {

            int neighborDisk = prox[currentDisk][i];

            double dx = sample[neighborDisk].x - testx;
            double dy = sample[neighborDisk].y - testy;
            double d = sqrt(dx * dx + dy * dy);
            if (d < testr + sample[neighborDisk].r + distMin) {
              ok = false;
              break;
            }
          }
        }

        if (ok == true) {
          found = true;
          Disk P(testx, testy, testr);
          sample.push_back(P);
          prox.push_back(std::vector<int>());
          int particleIndex = sample.size() - 1;

          for (int i = 0; i < particleIndex; i++) {
            // if (i == particleIndex)
            //              continue;
            double dx = fabs(sample[i].x - testx);
            double dy = fabs(sample[i].y - testy);
            double d = sqrt(dx * dx + dy * dy);
            if (d < sample[i].r + testr + distNeighbor + distMin) {
              prox[i].push_back(particleIndex);
              prox[particleIndex].push_back(i);
            }
          }

          active.push_back(particleIndex);
          break;
        }
      } // n-loop

      if (!found) {
        for (int i = randIndex; i < (int)active.size() - 1; i++) {
          active[i] = active[i + 1];
        }
        active.pop_back();
      }
    } // end-while
  }   // end-method-run

  /**
   * @brief Periodic execution of the packing algorithm.
   *
   * This is the version of the algorithm with periodic boundary conditions.
   *
   * @note This method is more efficient than the non-periodic version, but
   * is only suitable for packing disks in a periodic domain (i.e. a torus).
   */
  void execPeriodic() {
    // step 1
    if (sample.empty()) {
      double firstr = ran(rmin, rmax);
      Disk P(ran(xmin + firstr, xmax - firstr), ran(ymin + firstr, ymax - firstr), firstr);
      sample.push_back(P);
      prox.push_back(std::vector<int>());
      active.push_back(sample.size() - 1);
    }

    // step 2
    while (active.size() > 0 && sample.size() < (size_t)max) {
      int randIndex = rand() % active.size();
      int currentDisk = active[randIndex];
      double packedx = sample[currentDisk].x;
      double packedy = sample[currentDisk].y;
      double packedr = sample[currentDisk].r;

      bool found = false;

      for (int n = 0; n < k; n++) {
        double testr = ran(rmin, rmax);
        double angle = ran(0, 2.0 * M_PI);
        double m = ran(testr + packedr + distMin, testr + packedr + distMin + gapTol);
        double testx = packedx + m * cos(angle);
        double testy = packedy + m * sin(angle);

        bool ok = true;
        // boundaries
        if (testx < xmin || testx > xmax || testy < ymin || testy > ymax) {
          ok = false;
        }
        double dv = 2.0 * rmax + distMin + gapTol;
        if (testx < xmin + dv || testx > xmax - dv || testy < ymin + dv || testy > ymax - dv) {
          double lx = xmax - xmin;
          double half_lx = 0.5 * lx;
          double ly = ymax - ymin;
          double half_ly = 0.5 * ly;

          for (size_t i = 0; i < boundaries.size(); i++) {

            int neighborDisk = boundaries[i];

            double dx = sample[neighborDisk].x - testx;
            double dy = sample[neighborDisk].y - testy;

            if (dx > half_lx) {
              dx -= lx;
            } else if (dx < -half_lx) {
              dx += lx;
            }
            if (dy > half_ly) {
              dy -= ly;
            } else if (dy < -half_ly) {
              dy += ly;
            }

            double d = sqrt(dx * dx + dy * dy);
            if (d < testr + sample[neighborDisk].r + distMin) {
              ok = false;
              break;
            }
          }
        }

        // inter-particles
        if (ok == true) {
          for (size_t i = 0; i < prox[currentDisk].size(); i++) {
            int neighborDisk = prox[currentDisk][i];

            double dx = sample[neighborDisk].x - testx;
            double dy = sample[neighborDisk].y - testy;
            double d = sqrt(dx * dx + dy * dy);
            if (d < testr + sample[neighborDisk].r + distMin) {
              ok = false;
              break;
            }
          }
        }

        if (ok == true) {
          found = true;
          Disk P(testx, testy, testr);
          sample.push_back(P);
          prox.push_back(std::vector<int>());

          int particleIndex = sample.size() - 1;
          double dv = 2.0 * rmax + distMin + gapTol; // a verifier
          if (testx < xmin + dv || testx > xmax - dv || testy < ymin + dv || testy > ymax - dv) {
            boundaries.push_back(particleIndex);
          }

          for (int i = 0; i < particleIndex; i++) {
            double dx = fabs(sample[i].x - testx);
            double dy = fabs(sample[i].y - testy);
            double d = sqrt(dx * dx + dy * dy);
            if (d < sample[i].r + testr + distNeighbor + distMin) {
              prox[i].push_back(particleIndex);
              prox[particleIndex].push_back(i);
            }
          }

          active.push_back(particleIndex);
          break;
        }
      } // n-loop

      if (!found) {
        for (int i = randIndex; i < (int)active.size() - 1; i++) {
          active[i] = active[i + 1];
        }
        active.pop_back();
      }
    } // end-while
  }   // end-execPeriodic


  // With gnuplot, the disks can be plotted as follows:
  // set size ratio -1
  // plot "points.txt" u 1:2:4 w circles lc rgb "forest-green" fs transparent solid 0.15  notitle
  void save(const char *name, size_t from = 0) {
    std::ofstream file(name);
    for (size_t i = from; i < sample.size(); i++) {
      file << sample[i].x << ' ' << sample[i].y << " 0 " << sample[i].r << '\n';
    }
  }

protected:
  double ran(double min_, double max_) { return min_ + (rand() / (double)RAND_MAX) * (max_ - min_); }
};

#endif /* end of include guard: GEOPACK2D_HPP */

#if 0

int main(int argc, char const *argv[]) {
  GeoPack2D GP(0.15, 0.15, 5000, 0, 20, 0, 20, 0.0, 0);
  GP.seedTime();
  
  GP.appendDisk(10, 10, 5);
  GP.appendDisk(4, 4, 2);
  GP.appendDisk(16, 4, 2);
  GP.appendDisk(16, 16, 2);
  GP.appendDisk(4, 16, 2);  
  
  GP.reActivate();
  GP.execPeriodic();
  
  GP.parameters(0.05, 0.15, 5000, 0, 20, 0, 20, 0.0, 0);
  GP.reActivate();
  GP.execPeriodic();
  
  GP.save("disks.txt", 5);
  
  return 0;
}

#endif

