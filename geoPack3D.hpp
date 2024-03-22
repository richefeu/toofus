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

#ifndef GEOPACK3D_HPP
#define GEOPACK3D_HPP

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

class GeoPack3D {
public:
  struct Sphere {
    double x, y, z;
    double r;
    size_t nbNeighbors;
    Sphere(double x_, double y_, double z_, double r_) {
      x = x_;
      y = y_;
      z = z_;
      r = r_;
      nbNeighbors = 0;
    }
  };

  bool verbose{true};  // verbosity during packing
  double rmin;         // minimuum radius
  double rmax;         // maximum radius
  double gapTol{0.0};  // a (positive) tolerance on the space between the spheres
  double distNeighbor; // inter-sphere distance for neighbourhood lists
  double distMin{0.0}; // minimum distance for placing spheres
  size_t max;          // maximum number of spheres placed
  size_t k;            // number of placement attempts around a placed sphere before removing it from the asset list

  size_t limitLocalNumberNeighbors;
  size_t localNumberNeighborsMax;
  double distProbingConnectivity;

  size_t limitLocalSolidFraction;
  double distProbingSolidFraction;
  double localSolidFractionMax;

  std::vector<Sphere> sample;            // the packed spheres
  std::vector<std::vector<size_t>> prox; // each sphere has a list of neighbors
  std::vector<size_t> boundaries;        // indices of spheres that are close to the periodic boundaries
  std::vector<size_t> active;            // indices of active spheres

  double xmin{0.0}, xmax{1.0};
  double ymin{0.0}, ymax{1.0};
  double zmin{0.0}, zmax{1.0};

  // Ctor
  GeoPack3D() {
    rmin = rmax = 0.0;
    gapTol = distNeighbor = distMin = 0.0;
    max = 0;
    k = 0;
    xmin = xmax = ymin = ymax = zmin = zmax = 0.0;
    limitLocalNumberNeighbors = 0;
    limitLocalSolidFraction = 0;
    verbose = true;
  }

  GeoPack3D(double rmin_, double rmax_, size_t k_, double xmin_, double xmax_, double ymin_, double ymax_, double zmin_,
            double zmax_, double gap_ = 0.0, size_t max_ = 0) {
    parameters(rmin_, rmax_, k_, xmin_, xmax_, ymin_, ymax_, zmin_, zmax_, gap_, max_);
    verbose = true;
  }

  void parameters(double rmin_, double rmax_, size_t k_, double xmin_, double xmax_, double ymin_, double ymax_,
                  double zmin_, double zmax_, double gapTol_ = 0.0, size_t max_ = 0) {
    rmin = rmin_;
    rmax = rmax_;
    gapTol = gapTol_;
    distNeighbor = 2.0 * rmax + gapTol;
    distMin = 0.0;
    k = k_;
    localNumberNeighborsMax = 12;

    // Domain
    xmin = xmin_;
    xmax = xmax_;
    ymin = ymin_;
    ymax = ymax_;
    zmin = zmin_;
    zmax = zmax_;

    max = max_;
    if (max == 0) {
      size_t nx = static_cast<size_t>(floor((xmax - xmin) / rmin));
      size_t ny = static_cast<size_t>(floor((ymax - ymin) / rmin));
      size_t nz = static_cast<size_t>(floor((zmax - zmin) / rmin));
      max = nx * ny * nz;
    }
  }

  // Overall solid fraction by assuming there is no overlap
  double getSolidFraction() {
    double Vtot = (xmax - xmin) * (ymax - ymin) * (zmax - zmin);
    if (Vtot <= 0.0)
      return 0.0;
    double Vs = 0.0;

    for (size_t i = 0; i < sample.size(); i++) {
      Vs += 4.0 * M_PI * sample[i].r * sample[i].r * sample[i].r / 3.0;
    }
    return Vs / Vtot;
  }

  double sphereLensVolume(double R1, double R2, double dist) {
    if (dist > R1 + R2)
      return 0.0;

    double Rmin = std::min(R1, R2);
    double Rmax = std::max(R1, R2);

    if (dist + Rmin < Rmax) {
      return 4.0 * M_PI * Rmin * Rmin * Rmin / 3.0;
    }
    double a = R1 + R2 - dist;
    return (M_PI * (a * a) *
            (dist * dist + 2.0 * dist * R2 - 3.0 * R2 * R2 + 2.0 * dist * R1 + 6.0 * R2 * R1 - 3.0 * R1 * R1) /
            (12.0 * dist));
  }

  // it returns the solide fraction in a layer of 'width' around the sphere i
  // BEFORE another sphere is added
  double localSolidFraction(size_t i, double width, bool periodic = false) {

    double Rprob = sample[i].r + width;
    double Vprob = 4.0 * M_PI * (Rprob * Rprob * Rprob - sample[i].r * sample[i].r * sample[i].r) / 3.0;
    double Vs = 0.0;
    double Lx = (xmax - xmin), Ly = (ymax - ymin), Lz = (zmax - zmin);

    for (size_t n = 0; n < prox[i].size(); n++) {
      size_t j = prox[i][n];
      double dx = sample[j].x - sample[i].x;
      double dy = sample[j].y - sample[i].y;
      double dz = sample[j].z - sample[i].z;
      if (periodic == true) {
        dx -= floor(dx / Lx + 0.5) * Lx;
        dy -= floor(dy / Ly + 0.5) * Ly;
        dz -= floor(dz / Lz + 0.5) * Lz;
      }
      double dst = sqrt(dx * dx + dy * dy + dz * dz);
      Vs += sphereLensVolume(Rprob, sample[j].r, dst);
    }

    return Vs / Vprob;
  }

  void limit_localNumberNeighbors(double dst, size_t value) {
    limitLocalNumberNeighbors = 1;
    distProbingConnectivity = dst;
    localNumberNeighborsMax = value;
    std::cout << "localNumberNeighborsMax = " << localNumberNeighborsMax << '\n';
  }

  void limit_localSolidFraction(double dst, double value) {
    limitLocalSolidFraction = 1;
    distProbingSolidFraction = dst;
    localSolidFractionMax = value;
  }

  void seedTime() { srand(static_cast<unsigned int>(time(NULL))); }

  void reActivate(size_t from = 0, size_t to = 0) {
    if (to == 0)
      to = sample.size();
    for (size_t i = from; i < to; i++) {
      active.push_back(i);
    }
  }

  void deActivate(size_t index) {
    for (size_t i = index; i < active.size() - 1; i++) {
      active[i] = active[i + 1];
    }
    active.pop_back();
  }

  // execute the Poisson Sampling algorithm
  void exec(bool skipFirstStep = false) {
    // step 1
    if (skipFirstStep == false) {
      double firstr = ran(rmin, rmax);
      Sphere P(ran(xmin + firstr, xmax - firstr), ran(ymin + firstr, ymax - firstr), ran(zmin + firstr, zmax - firstr),
               firstr);
      sample.push_back(P);
      prox.push_back(std::vector<size_t>());
      active.push_back(sample.size() - 1);
    }

    // step 2
    size_t count = 0;
    size_t countMax = (size_t)floor(2.0 * (xmax - xmin) / (rmin + rmax));
    while (active.size() > 0 && sample.size() < (size_t)max) {
      size_t randIndex = static_cast<size_t>(rand()) % active.size();
      if (limitLocalNumberNeighbors == 1 && sample[active[randIndex]].nbNeighbors >= localNumberNeighborsMax) {
        deActivate(randIndex);
        continue;
      }
      if (limitLocalSolidFraction == 1 &&
          localSolidFraction(active[randIndex], distProbingSolidFraction) > localSolidFractionMax) {
        deActivate(randIndex);
        continue;
      }
      size_t currentSphere = active[randIndex];

      bool found = false;

      for (size_t n = 0; n < k; n++) {
        double testr = ran(rmin, rmax);
        double angle1 = ran(0, 2.0 * M_PI);
        double angle2 = ran(-0.5 * M_PI, 0.5 * M_PI);
        double m = ran(testr + sample[currentSphere].r + distMin, testr + sample[currentSphere].r + distMin + gapTol);
        double ux = cos(angle1);
        double uy = sin(angle1);
        double uz = sin(angle2);
        double sc = 1.0f / sqrt(ux * ux + uy * uy + uz * uz);
        double testx = sample[currentSphere].x + m * sc * ux;
        double testy = sample[currentSphere].y + m * sc * uy;
        double testz = sample[currentSphere].z + m * sc * uz;

        bool ok = true;

        // boundaries
        if (ok == true) {
          if (testx < xmin + testr || testx > xmax - testr || testy < ymin + testr || testy > ymax - testr ||
              testz < zmin + testr || testz > zmax - testr) {
            ok = false;
          }
        }

        // inter-particles
        if (ok == true) {
          for (size_t i = 0; i < prox[currentSphere].size(); i++) {
            size_t neighborSphere = prox[currentSphere][i];

            double dx = sample[neighborSphere].x - testx;
            double dy = sample[neighborSphere].y - testy;
            double dz = sample[neighborSphere].z - testz;
            double d = sqrt(dx * dx + dy * dy + dz * dz);
            if (d < testr + sample[neighborSphere].r + distMin) {
              ok = false;
              break;
            }
          }
        }

        // add if ok
        if (ok == true) {
          found = true;
          Sphere P(testx, testy, testz, testr);
          sample.push_back(P);
          prox.push_back(std::vector<size_t>());
          size_t particleIndex = sample.size() - 1;

          for (size_t i = 0; i < sample.size(); i++) {
            if (i == particleIndex)
              continue;
            double dx = fabs(sample[i].x - testx);
            double dy = fabs(sample[i].y - testy);
            double dz = fabs(sample[i].z - testz);
            double d = sqrt(dx * dx + dy * dy + dz * dz);
            if (d <= sample[i].r + testr + distNeighbor + distMin) {
              prox[i].push_back(particleIndex);
              prox[particleIndex].push_back(i);
            }
            if (d <= sample[i].r + testr + distProbingConnectivity) {
              sample[i].nbNeighbors += 1;
              sample[particleIndex].nbNeighbors += 1;
            }
          }

          active.push_back(particleIndex);
          break;
        }
      } // n-loop

      if (!found) {
        deActivate(randIndex);
      }

      count++;
      if (count >= countMax) {
        count = 0;

        if (verbose == true) {
          std::cout << "packed: " << sample.size() << ", ";
          std::cout << "active: " << active.size() << std::endl;
        }
      }
    } // end-while
  }   // end-method-run

  // Periodic version of the Poisson Sampling algorithm
  void execPeriodic(bool skipFirstStep = false) {
    double Lx = (xmax - xmin), Ly = (ymax - ymin), Lz = (zmax - zmin);

    // step 1
    if (skipFirstStep == false) {
      double firstr = ran(rmin, rmax);
      Sphere P(0.5 * Lx, 0.5 * Ly, 0.5 * Lz, firstr);
      sample.push_back(P);
      prox.push_back(std::vector<size_t>());
      active.push_back(sample.size() - 1);
    }

    // step 2
    size_t count = 0;
    size_t countMax = static_cast<size_t>(floor((xmax - xmin) / (0.5 * (rmin + rmax))));
    while (active.size() > 0 && sample.size() < (size_t)max) {

      size_t randIndex = static_cast<size_t>(rand()) % active.size();

      // std::cout << sample[active[randIndex]].nbNeighbors << '\n';
      if (limitLocalNumberNeighbors == 1 && sample[active[randIndex]].nbNeighbors >= localNumberNeighborsMax) {
        // std::cout << "**** deactivate " << active[randIndex] << '\n';
        deActivate(randIndex);
        continue;
      }

      if (limitLocalSolidFraction == 1 &&
          localSolidFraction(active[randIndex], distProbingSolidFraction, true) > localSolidFractionMax) {
        deActivate(randIndex);
        continue;
      }
      size_t currentSphere = active[randIndex];

      bool found = false;

      // try k times to add a sphere arround it
      for (size_t n = 0; n < k; n++) {

        // an attempt of positionning
        double testr = ran(rmin, rmax);
        double angle1 = ran(0, 2.0 * M_PI);
        double angle2 = ran(-0.5 * M_PI, 0.5 * M_PI);
        double m = ran(testr + sample[currentSphere].r + distMin, testr + sample[currentSphere].r + distMin + gapTol);
        double ux = cos(angle1);
        double uy = sin(angle1);
        double uz = sin(angle2);
        double sc = 1.0f / sqrt(ux * ux + uy * uy + uz * uz);
        double testx = sample[currentSphere].x + m * sc * ux;
        double testy = sample[currentSphere].y + m * sc * uy;
        double testz = sample[currentSphere].z + m * sc * uz;

        if (testx < xmin) {
          testx += Lx;
        }
        if (testx > xmax) {
          testx -= Lx;
        }
        if (testy < ymin) {
          testy += Ly;
        }
        if (testy > ymax) {
          testy -= Ly;
        }
        if (testz < zmin) {
          testz += Lz;
        }
        if (testz > zmax) {
          testz -= Lz;
        }

        bool ok = true;

        // boundary limits
        if (testx < xmin || testx > xmax || testy < ymin || testy > ymax || testz < zmin || testz > zmax) {
          ok = false;
        }

        // inter-spheres thoughout the periodicity
        if (ok == true) {
          double dv = 2.0 * rmax + distMin + gapTol;
          if (testx < xmin + dv || testx > xmax - dv || testy < ymin + dv || testy > ymax - dv || testz < zmin + dv ||
              testz > zmax - dv) {

            for (size_t i = 0; i < boundaries.size(); i++) {
              size_t neighborDisk = boundaries[i];

              double dx = sample[neighborDisk].x - testx;
              double dy = sample[neighborDisk].y - testy;
              double dz = sample[neighborDisk].z - testz;

              dx -= floor(dx / Lx + 0.5) * Lx;
              dy -= floor(dy / Ly + 0.5) * Ly;
              dz -= floor(dz / Lz + 0.5) * Lz;

              double d = sqrt(dx * dx + dy * dy + dz * dz);
              if (d < testr + sample[neighborDisk].r + distMin) {
                ok = false;
                break;
              }
            }
          }
        }

        // inter-particles
        if (ok == true) {
          for (size_t i = 0; i < prox[currentSphere].size(); i++) {
            size_t neighborSphere = prox[currentSphere][i];

            double dx = sample[neighborSphere].x - testx;
            double dy = sample[neighborSphere].y - testy;
            double dz = sample[neighborSphere].z - testz;

            dx -= floor(dx / Lx + 0.5) * Lx;
            dy -= floor(dy / Ly + 0.5) * Ly;
            dz -= floor(dz / Lz + 0.5) * Lz;

            double d = sqrt(dx * dx + dy * dy + dz * dz);
            if (d < testr + sample[neighborSphere].r + distMin) {
              ok = false;
              break;
            }
          }
        }

        if (ok == true) {
          found = true;
          Sphere P(testx, testy, testz, testr);
          sample.push_back(P);
          prox.push_back(std::vector<size_t>());
          size_t particleIndex = sample.size() - 1;

          double dv = 2.0 * rmax + distMin + gapTol;
          if (testx < xmin + dv || testx > xmax - dv || testy < ymin + dv || testy > ymax - dv || testz < zmin + dv ||
              testz > zmax - dv) {
            boundaries.push_back(particleIndex);
          }

          for (size_t i = 0; i < particleIndex; i++) {
            double dx = sample[i].x - testx;
            double dy = sample[i].y - testy;
            double dz = sample[i].z - testz;

            dx -= floor(dx / Lx + 0.5) * Lx;
            dy -= floor(dy / Ly + 0.5) * Ly;
            dz -= floor(dz / Lz + 0.5) * Lz;

            double d = sqrt(dx * dx + dy * dy + dz * dz);
            if (d < sample[i].r + testr + distNeighbor + distMin) {
              prox[i].push_back(particleIndex);
              prox[particleIndex].push_back(i);
            }
            if (d <= sample[i].r + testr + distProbingConnectivity) {
              sample[i].nbNeighbors += 1;
              sample[particleIndex].nbNeighbors += 1;
            }
          }

          active.push_back(particleIndex);
          break;
        }
      } // n-loop (k times)

      // if not found, remove the sphere from the active ones
      if (!found) {
        deActivate(randIndex);
      }

      count++;
      if (count >= countMax) {
        count = 0;
        if (verbose == true) {
          std::cout << "packed: " << sample.size() << ", ";
          std::cout << "active: " << active.size() << std::endl;
        }
      }
    } // end-while
  }   // end-execPeriodic

  // Pour visualiser :
  // compiler 'seeSpheres' puis > seeSpheres points.txt
  void save(const char *name) {
    std::ofstream file(name);
    for (size_t i = 0; i < sample.size(); i++) {
      file << sample[i].x << ' ' << sample[i].y << ' ' << sample[i].z << ' ' << sample[i].r << '\n';
    }
  }

  void saveBox(const char *name) {
    std::ofstream file(name);
    file << xmin << ' ' << ymin << ' ' << zmin << '\n';
    file << xmax << ' ' << ymax << ' ' << zmax << '\n';
  }

  void saveVTK(const char *name) {
    std::ofstream fog(name, std::ios::out);
    if (fog) {
      fog.precision(5);
      fog << std::scientific;
      fog << "# vtk DataFile Version 3.0\n";
      fog << "My grains\n";
      fog << "ASCII\n";
      fog << "DATASET UNSTRUCTURED_GRID\n";
      fog << "POINTS " << sample.size() << " float\n";
      for (size_t i = 0; i < sample.size(); i++)
        fog << sample[i].x << " " << sample[i].y << " " << sample[i].z << '\n';
      fog << "POINT_DATA " << sample.size() << '\n';
      fog << "SCALARS Diameter float\n";
      fog << "LOOKUP_TABLE default\n";
      for (size_t i = 0; i < sample.size(); i++)
        fog << 2 * sample[i].r << '\n';
    }
  }

protected:
  double ran(double min_, double max_) { return min_ + (rand() / (double)RAND_MAX) * (max_ - min_); }
};

#endif /* end of include guard: GEOPACK2D_HPP */

#if 0

int main(int argc, char const *argv[]) {
  
  GeoPack3D GP(0.8, 3, 100, 0, 100, 0, 100, 0, 100);
  GP.seedTime();
  GP.exec();
  GP.save("packing.txt");
  
  return 0;
}

#endif
