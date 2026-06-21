// STATUS: [ ] STABLE  [x] EXPERIMENTAL  [ ] DRAFT

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

#include <algorithm>
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
      x           = x_;
      y           = y_;
      z           = z_;
      r           = r_;
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

  size_t limitLocalNumberNeighbors{0};
  size_t localNumberNeighborsMax{12};
  double distProbingConnectivity{0.0};

  size_t limitLocalSolidFraction{0};
  double distProbingSolidFraction{0.0};
  double localSolidFractionMax{0.0};

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
    max                             = 0;
    k                               = 0;
    xmin = xmax = ymin = ymax = zmin = zmax = 0.0;
    limitLocalNumberNeighbors               = 0;
    limitLocalSolidFraction                 = 0;
    verbose                                 = true;
  }

  GeoPack3D(double rmin_, double rmax_, size_t k_, double xmin_, double xmax_, double ymin_, double ymax_, double zmin_,
            double zmax_, double gap_ = 0.0, size_t max_ = 0) {
    parameters(rmin_, rmax_, k_, xmin_, xmax_, ymin_, ymax_, zmin_, zmax_, gap_, max_);
    verbose = true;
  }

  void parameters(double rmin_, double rmax_, size_t k_, double xmin_, double xmax_, double ymin_, double ymax_,
                  double zmin_, double zmax_, double gapTol_ = 0.0, size_t max_ = 0) {
    rmin                    = rmin_;
    rmax                    = rmax_;
    gapTol                  = gapTol_;
    distNeighbor            = 2.0 * rmax + gapTol;
    distMin                 = 0.0;
    k                       = k_;
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
      max       = nx * ny * nz;
    }
  }

  // Overall solid fraction by assuming there is no overlap
  double getSolidFraction() {
    double Vtot = (xmax - xmin) * (ymax - ymin) * (zmax - zmin);
    if (Vtot <= 0.0) return 0.0;
    double Vs = 0.0;

    for (size_t i = 0; i < sample.size(); i++) { Vs += 4.0 * M_PI * sample[i].r * sample[i].r * sample[i].r / 3.0; }
    return Vs / Vtot;
  }

  // --- Diagnostics for static stability ------------------------------------
  // Two spheres are considered in contact when their (periodic) gap is <= contactGap.
  // A small positive contactGap (e.g. a fraction of rmin) tolerates the near-tangencies
  // produced by the geometric packing.

  // Mean coordination number over the whole sample (periodic, rattlers included).
  // For frictional spheres, static rigidity requires z >~ 4 (z = 6 frictionless),
  // so the loose tangent packing (z ~ 2) is far below the isostatic threshold.
  double getCoordinationNumber(double contactGap = 0.0) {
    if (sample.empty()) return 0.0;
    double Lx = (xmax - xmin), Ly = (ymax - ymin), Lz = (zmax - zmin);
    size_t nbContacts = 0;
    for (size_t i = 0; i < sample.size(); i++) {
      for (size_t n = 0; n < prox[i].size(); n++) {
        size_t j = prox[i][n];
        if (j <= i) continue; // count each undirected pair once
        double dx = sample[j].x - sample[i].x;
        double dy = sample[j].y - sample[i].y;
        double dz = sample[j].z - sample[i].z;
        dx -= floor(dx / Lx + 0.5) * Lx;
        dy -= floor(dy / Ly + 0.5) * Ly;
        dz -= floor(dz / Lz + 0.5) * Lz;
        double d = sqrt(dx * dx + dy * dy + dz * dz);
        if (d <= sample[i].r + sample[j].r + contactGap) { nbContacts++; }
      }
    }
    return 2.0 * static_cast<double>(nbContacts) / static_cast<double>(sample.size());
  }

  // Number of connected components of the (periodic) contact graph: 1 means the
  // whole assembly is a single connected cluster, > 1 means there are floating
  // islands / isolated spheres that the DEM cannot hold in place.
  size_t getNumberOfConnectedComponents(double contactGap = 0.0) {
    size_t N = sample.size();
    if (N == 0) return 0;
    std::vector<size_t> parent(N);
    for (size_t i = 0; i < N; i++) { parent[i] = i; }
    double Lx = (xmax - xmin), Ly = (ymax - ymin), Lz = (zmax - zmin);
    for (size_t i = 0; i < N; i++) {
      for (size_t n = 0; n < prox[i].size(); n++) {
        size_t j = prox[i][n];
        if (j <= i) continue;
        double dx = sample[j].x - sample[i].x;
        double dy = sample[j].y - sample[i].y;
        double dz = sample[j].z - sample[i].z;
        dx -= floor(dx / Lx + 0.5) * Lx;
        dy -= floor(dy / Ly + 0.5) * Ly;
        dz -= floor(dz / Lz + 0.5) * Lz;
        double d = sqrt(dx * dx + dy * dy + dz * dz);
        if (d <= sample[i].r + sample[j].r + contactGap) {
          size_t ri = i;
          while (parent[ri] != ri) { ri = parent[ri]; }
          size_t rj = j;
          while (parent[rj] != rj) { rj = parent[rj]; }
          if (ri != rj) { parent[ri] = rj; }
        }
      }
    }
    size_t comps = 0;
    for (size_t i = 0; i < N; i++) {
      if (parent[i] == i) { comps++; }
    }
    return comps;
  }

  // --- Densification helpers ------------------------------------------------
  // Scale all radii by a common factor, centers unchanged. With factor > 1 this
  // turns the near-tangencies into overlaps, building a dense and highly
  // connected initial state to be relaxed by a periodic DEM compression.
  void inflateRadii(double factor) {
    for (size_t i = 0; i < sample.size(); i++) { sample[i].r *= factor; }
  }

  // Inflate radii to reach a target overall solid fraction (no-overlap formula).
  // Returns the applied factor. After inflation, getSolidFraction() reports
  // ~phiTarget but the true (overlap-corrected) fraction is slightly lower; the
  // DEM removes the overlaps at nearly constant volume.
  double inflateToSolidFraction(double phiTarget) {
    double phi = getSolidFraction();
    if (phi <= 0.0) return 1.0;
    double factor = std::cbrt(phiTarget / phi);
    inflateRadii(factor);
    return factor;
  }

  double sphereLensVolume(double R1, double R2, double dist) {
    if (dist > R1 + R2) return 0.0;

    double Rmin = std::min(R1, R2);
    double Rmax = std::max(R1, R2);

    if (dist + Rmin < Rmax) { return 4.0 * M_PI * Rmin * Rmin * Rmin / 3.0; }
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
    double Vs    = 0.0;
    double Lx = (xmax - xmin), Ly = (ymax - ymin), Lz = (zmax - zmin);

    for (size_t n = 0; n < prox[i].size(); n++) {
      size_t j  = prox[i][n];
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
    distProbingConnectivity   = dst;
    localNumberNeighborsMax   = value;
    if (verbose == true) { std::cout << "localNumberNeighborsMax = " << localNumberNeighborsMax << '\n'; }
  }

  // Hard caps evaluated on a candidate sphere (cx, cy, cz, cr) BEFORE it is accepted, by scanning
  // ALL spheres (periodic minimum image if periodic == true). Returns false when accepting the
  // candidate would violate an active hard cap:
  //   - coordination: the candidate would get more than localNumberNeighborsMax contacts, or it
  //     would push an already-saturated neighbor above the limit. Combined with this test applied
  //     to every placement, it guarantees the invariant z_i <= localNumberNeighborsMax for all i;
  //   - local solid fraction: the candidate's shell (width distProbingSolidFraction) would be
  //     filled above localSolidFractionMax by the already-placed neighbors. This is enforced at
  //     the candidate's placement; the shell of a neighbor may still grow as the packing fills.
  bool candidateRespectsHardCaps(double cx, double cy, double cz, double cr, bool periodic) {
    bool checkZ   = (limitLocalNumberNeighbors == 1);
    bool checkPhi = (limitLocalSolidFraction == 1);
    if (checkZ == false && checkPhi == false) { return true; }

    double Lx = (xmax - xmin), Ly = (ymax - ymin), Lz = (zmax - zmin);
    double Rprob      = cr + distProbingSolidFraction;
    double Vprob      = 4.0 * M_PI * (Rprob * Rprob * Rprob - cr * cr * cr) / 3.0;
    double Vs         = 0.0;
    size_t nbContacts = 0;

    for (size_t i = 0; i < sample.size(); i++) {
      double dx = sample[i].x - cx;
      double dy = sample[i].y - cy;
      double dz = sample[i].z - cz;
      if (periodic == true) {
        dx -= floor(dx / Lx + 0.5) * Lx;
        dy -= floor(dy / Ly + 0.5) * Ly;
        dz -= floor(dz / Lz + 0.5) * Lz;
      }
      double d = sqrt(dx * dx + dy * dy + dz * dz);
      if (checkZ == true && d <= sample[i].r + cr + distProbingConnectivity) {
        nbContacts++;
        if (sample[i].nbNeighbors >= localNumberNeighborsMax) { return false; } // neighbor would exceed
      }
      if (checkPhi == true) { Vs += sphereLensVolume(Rprob, sample[i].r, d); }
    }

    if (checkZ == true && nbContacts > localNumberNeighborsMax) { return false; }
    if (checkPhi == true && Vprob > 0.0 && Vs / Vprob > localSolidFractionMax) { return false; }
    return true;
  }

  void limit_localSolidFraction(double dst, double value) {
    limitLocalSolidFraction  = 1;
    distProbingSolidFraction = dst;
    localSolidFractionMax    = value;
    // The prox lists must reach at least as far as the probing shell, otherwise
    // neighbors located between distNeighbor and distProbingSolidFraction are
    // missed and the local solid fraction is underestimated.
    if (distProbingSolidFraction > distNeighbor) { distNeighbor = distProbingSolidFraction; }
  }

  void seedTime() {
    srand(static_cast<unsigned int>(time(NULL)));
  }

  void reActivate(size_t from = 0, size_t to = 0) {
    if (to == 0) to = sample.size();
    for (size_t i = from; i < to; i++) { active.push_back(i); }
  }

  void deActivate(size_t index) {
    for (size_t i = index; i < active.size() - 1; i++) { active[i] = active[i + 1]; }
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
    size_t count    = 0;
    size_t countMax = (size_t)floor(2.0 * (xmax - xmin) / (rmin + rmax));
    while (active.size() > 0 && sample.size() < (size_t)max) {
      size_t randIndex = static_cast<size_t>(rand()) % active.size();
      if (limitLocalNumberNeighbors == 1 && sample[active[randIndex]].nbNeighbors >= localNumberNeighborsMax) {
        deActivate(randIndex);
        continue;
      }
      size_t currentSphere = active[randIndex];

      bool found = false;

      for (size_t n = 0; n < k; n++) {
        double testr  = ran(rmin, rmax);
        double angle1 = ran(0, 2.0 * M_PI);
        double uz     = ran(-1.0, 1.0); // cos(theta) uniform => isotropic directions on the unit sphere
        double sxy    = sqrt(1.0 - uz * uz);
        double m  = ran(testr + sample[currentSphere].r + distMin, testr + sample[currentSphere].r + distMin + gapTol);
        double ux = sxy * cos(angle1);
        double uy = sxy * sin(angle1);
        double sc = 1.0; // (ux, uy, uz) is already a unit vector
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
            double d  = sqrt(dx * dx + dy * dy + dz * dz);
            if (d < testr + sample[neighborSphere].r + distMin) {
              ok = false;
              break;
            }
          }
        }

        // hard caps on the candidate (coordination and/or local solid fraction)
        if (ok == true && candidateRespectsHardCaps(testx, testy, testz, testr, false) == false) { ok = false; }

        // add if ok
        if (ok == true) {
          found = true;
          Sphere P(testx, testy, testz, testr);
          sample.push_back(P);
          prox.push_back(std::vector<size_t>());
          size_t particleIndex = sample.size() - 1;

          for (size_t i = 0; i < sample.size(); i++) {
            if (i == particleIndex) continue;
            double dx = fabs(sample[i].x - testx);
            double dy = fabs(sample[i].y - testy);
            double dz = fabs(sample[i].z - testz);
            double d  = sqrt(dx * dx + dy * dy + dz * dz);
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

      if (!found) { deActivate(randIndex); }

      count++;
      if (count >= countMax) {
        count = 0;

        if (verbose == true) {
          std::cout << "packed: " << sample.size() << ", ";
          std::cout << "active: " << active.size() << std::endl;
        }
      }
    } // end-while
  } // end-method-run

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
    size_t count    = 0;
    size_t countMax = static_cast<size_t>(floor((xmax - xmin) / (0.5 * (rmin + rmax))));
    while (active.size() > 0 && sample.size() < (size_t)max) {

      size_t randIndex = static_cast<size_t>(rand()) % active.size();

      // std::cout << sample[active[randIndex]].nbNeighbors << '\n';
      if (limitLocalNumberNeighbors == 1 && sample[active[randIndex]].nbNeighbors >= localNumberNeighborsMax) {
        // std::cout << "**** deactivate " << active[randIndex] << '\n';
        deActivate(randIndex);
        continue;
      }

      size_t currentSphere = active[randIndex];

      bool found = false;

      // try k times to add a sphere arround it
      for (size_t n = 0; n < k; n++) {

        // an attempt of positionning
        double testr  = ran(rmin, rmax);
        double angle1 = ran(0, 2.0 * M_PI);
        double uz     = ran(-1.0, 1.0); // cos(theta) uniform => isotropic directions on the unit sphere
        double sxy    = sqrt(1.0 - uz * uz);
        double m  = ran(testr + sample[currentSphere].r + distMin, testr + sample[currentSphere].r + distMin + gapTol);
        double ux = sxy * cos(angle1);
        double uy = sxy * sin(angle1);
        double sc = 1.0; // (ux, uy, uz) is already a unit vector
        double testx = sample[currentSphere].x + m * sc * ux;
        double testy = sample[currentSphere].y + m * sc * uy;
        double testz = sample[currentSphere].z + m * sc * uz;

        if (testx < xmin) { testx += Lx; }
        if (testx > xmax) { testx -= Lx; }
        if (testy < ymin) { testy += Ly; }
        if (testy > ymax) { testy -= Ly; }
        if (testz < zmin) { testz += Lz; }
        if (testz > zmax) { testz -= Lz; }

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

        // hard caps on the candidate (coordination and/or local solid fraction)
        if (ok == true && candidateRespectsHardCaps(testx, testy, testz, testr, true) == false) { ok = false; }

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
      if (!found) { deActivate(randIndex); }

      count++;
      if (count >= countMax) {
        count = 0;
        if (verbose == true) {
          std::cout << "packed: " << sample.size() << ", ";
          std::cout << "active: " << active.size() << std::endl;
        }
      }
    } // end-while
  } // end-execPeriodic

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
      for (size_t i = 0; i < sample.size(); i++) fog << sample[i].x << " " << sample[i].y << " " << sample[i].z << '\n';
      fog << "POINT_DATA " << sample.size() << '\n';
      fog << "SCALARS Diameter float\n";
      fog << "LOOKUP_TABLE default\n";
      for (size_t i = 0; i < sample.size(); i++) fog << 2 * sample[i].r << '\n';
    }
  }

  // Options for the cavalier-perspective SVG export (see saveSVG).
  struct SvgOptions {
    double alphaDeg   = 30.0; // angle of the receding z-axis, in degrees (typically 30 or 45)
    double depthScale = 0.5;  // foreshortening factor along z (typically 0.5)
    double ppl        = 60.0; // pixels per unit length

    bool drawCell        = true;  // edges of the periodic cell
    bool drawParticles   = true;  // the spheres, as disks
    bool drawConnections = false; // the contact network (segments between touching spheres)
    bool drawGhosts      = false; // periodic duplicates of the boundary particles (context layer)
    double contactGap    = 0.0;   // gap tolerance defining a contact for the network
  };

  // Cavalier-perspective SVG of the packing for quick visual checks.
  // Projection: Xscreen = x + depthScale*z*cos(alpha), Yscreen = y + depthScale*z*sin(alpha).
  // Particles are drawn back-to-front (painter's algorithm) and shaded by depth (near = light,
  // far = dark). Disk radii are kept true (cavalier preserves lengths in the x-y plane), so
  // Rmin/Rmax are read directly on the figure.
  // Layers are independent (like seeSpheres: 'l' cell, 'k' connections, 'i' ghosts):
  //  - the contact network is anchored on the REAL particle positions: an internal contact is a
  //    solid segment, a contact crossing the period is drawn as two short stubs, one from each
  //    real particle toward its partner's image, clipped at the cell wall (no duplicate is used);
  //  - ghosts are periodic duplicates of the boundary particles, drawn for context only and
  //    never used as connection endpoints.
  void saveSVG(const char *name) {
    saveSVG(name, SvgOptions());
  }

  void saveSVG(const char *name, const SvgOptions &opt) {
    double ca         = cos(opt.alphaDeg * M_PI / 180.0);
    double sa         = sin(opt.alphaDeg * M_PI / 180.0);
    double ppl        = opt.ppl;
    double depthScale = opt.depthScale;

    // cavalier projection of a 3D point onto the drawing plane (world units, origin at min corner)
    auto projX = [&](double x, double z) { return (x - xmin) + depthScale * (z - zmin) * ca; };
    auto projY = [&](double y, double z) { return (y - ymin) + depthScale * (z - zmin) * sa; };

    // drawing bounds from the 8 corners of the cell, with the largest radius as margin
    double cx[2] = {xmin, xmax}, cy[2] = {ymin, ymax}, cz[2] = {zmin, zmax};
    double sxmin = 1e30, sxmax = -1e30, symin = 1e30, symax = -1e30;
    for (int ix = 0; ix < 2; ix++) {
      for (int iy = 0; iy < 2; iy++) {
        for (int iz = 0; iz < 2; iz++) {
          double sx = projX(cx[ix], cz[iz]);
          double sy = projY(cy[iy], cz[iz]);
          sxmin     = std::min(sxmin, sx);
          sxmax     = std::max(sxmax, sx);
          symin     = std::min(symin, sy);
          symax     = std::max(symax, sy);
        }
      }
    }
    sxmin -= rmax;
    symin -= rmax;
    sxmax += rmax;
    symax += rmax;

    double W = (sxmax - sxmin) * ppl;
    double H = (symax - symin) * ppl;

    // map a world-projected point to SVG pixels (y flipped: SVG y grows downward)
    auto px = [&](double sx) { return (sx - sxmin) * ppl; };
    auto py = [&](double sy) { return H - (sy - symin) * ppl; };
    auto PX = [&](double x, double z) { return px(projX(x, z)); }; // 3D point -> pixel x
    auto PY = [&](double y, double z) { return py(projY(y, z)); }; // 3D point -> pixel y

    double Lx = (xmax - xmin), Ly = (ymax - ymin), Lz = (zmax - zmin);

    std::ofstream file(name);
    file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    file << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W << "\" height=\"" << H << "\" viewBox=\"0 0 " << W
         << ' ' << H << "\">\n";
    file << "<rect width=\"" << W << "\" height=\"" << H << "\" fill=\"white\"/>\n";

    // edges of the periodic cell (corner index = ix + 2*iy + 4*iz)
    if (opt.drawCell == true) {
      double X[8], Y[8];
      int c = 0;
      for (int iz = 0; iz < 2; iz++) {
        for (int iy = 0; iy < 2; iy++) {
          for (int ix = 0; ix < 2; ix++) {
            X[c] = PX(cx[ix], cz[iz]);
            Y[c] = PY(cy[iy], cz[iz]);
            c++;
          }
        }
      }
      static const int edges[12][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7},  // along x
                                       {0, 2}, {1, 3}, {4, 6}, {5, 7},  // along y
                                       {0, 4}, {1, 5}, {2, 6}, {3, 7}}; // along z
      for (int e = 0; e < 12; e++) {
        file << "<line x1=\"" << X[edges[e][0]] << "\" y1=\"" << Y[edges[e][0]] << "\" x2=\"" << X[edges[e][1]]
             << "\" y2=\"" << Y[edges[e][1]] << "\" stroke=\"#bbbbbb\" stroke-width=\"1\"/>\n";
      }
    }

    // periodic ghost particles: independent context layer (like seeSpheres key 'i').
    // Boundary particles are duplicated on the adjacent side(s); ghosts are NEVER used as
    // connection endpoints.
    if (opt.drawGhosts == true) {
      for (size_t i = 0; i < sample.size(); i++) {
        for (int sx = -1; sx <= 1; sx++) {
          for (int sy = -1; sy <= 1; sy++) {
            for (int sz = -1; sz <= 1; sz++) {
              if (sx == 0 && sy == 0 && sz == 0) continue;
              double gx = sample[i].x + sx * Lx;
              double gy = sample[i].y + sy * Ly;
              double gz = sample[i].z + sz * Lz;
              // keep only the duplicates adjacent to the cell (within the rmax margin of the view)
              if (gx < xmin - rmax || gx > xmax + rmax || gy < ymin - rmax || gy > ymax + rmax || gz < zmin - rmax ||
                  gz > zmax + rmax) {
                continue;
              }
              drawGhostDisk_SVG(file, PX(gx, gz), PY(gy, gz), sample[i].r * ppl);
            }
          }
        }
      }
    }

    // contact network, always anchored on the REAL particle positions.
    // An internal contact is a solid segment between the two centers. A contact crossing the
    // period is drawn as two short stubs: from each real particle toward its partner's image,
    // each clipped at the cell wall it exits (no duplicated-particle position is ever used).
    if (opt.drawConnections == true) {
      for (size_t i = 0; i < sample.size(); i++) {
        for (size_t n = 0; n < prox[i].size(); n++) {
          size_t j = prox[i][n];
          if (j <= i) continue; // each undirected pair once
          double dx  = sample[j].x - sample[i].x;
          double dy  = sample[j].y - sample[i].y;
          double dz  = sample[j].z - sample[i].z;
          double shx = floor(dx / Lx + 0.5) * Lx; // periodic shift removed to get the minimum image
          double shy = floor(dy / Ly + 0.5) * Ly;
          double shz = floor(dz / Lz + 0.5) * Lz;
          double dxm = dx - shx, dym = dy - shy, dzm = dz - shz;
          double d = sqrt(dxm * dxm + dym * dym + dzm * dzm);
          if (d > sample[i].r + sample[j].r + opt.contactGap) continue;

          bool wrap = (shx != 0.0 || shy != 0.0 || shz != 0.0);
          if (wrap == false) {
            // internal contact: full segment between the two real centers
            file << "<line x1=\"" << PX(sample[i].x, sample[i].z) << "\" y1=\"" << PY(sample[i].y, sample[i].z)
                 << "\" x2=\"" << PX(sample[j].x, sample[j].z) << "\" y2=\"" << PY(sample[j].y, sample[j].z)
                 << "\" stroke=\"#cc0000\" stroke-width=\"1.2\"/>\n";
          } else {
            // stub from i toward the image of j, clipped at the cell wall
            double ex, ey, ez;
            clipToCell_SVG(sample[i].x, sample[i].y, sample[i].z, sample[i].x + dxm, sample[i].y + dym,
                           sample[i].z + dzm, ex, ey, ez);
            file << "<line x1=\"" << PX(sample[i].x, sample[i].z) << "\" y1=\"" << PY(sample[i].y, sample[i].z)
                 << "\" x2=\"" << PX(ex, ez) << "\" y2=\"" << PY(ey, ez)
                 << "\" stroke=\"#cc0000\" stroke-width=\"1.2\"/>\n";
            // stub from j toward the image of i, clipped at the cell wall
            clipToCell_SVG(sample[j].x, sample[j].y, sample[j].z, sample[j].x - dxm, sample[j].y - dym,
                           sample[j].z - dzm, ex, ey, ez);
            file << "<line x1=\"" << PX(sample[j].x, sample[j].z) << "\" y1=\"" << PY(sample[j].y, sample[j].z)
                 << "\" x2=\"" << PX(ex, ez) << "\" y2=\"" << PY(ey, ez)
                 << "\" stroke=\"#cc0000\" stroke-width=\"1.2\"/>\n";
          }
        }
      }
    }

    // particles, back-to-front (painter's algorithm), shaded by depth
    if (opt.drawParticles == true) {
      std::vector<size_t> idx(sample.size());
      for (size_t i = 0; i < sample.size(); i++) { idx[i] = i; }
      std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return sample[a].z > sample[b].z; });

      for (size_t n = 0; n < idx.size(); n++) {
        size_t i = idx[n];
        double t = (Lz > 0.0) ? (sample[i].z - zmin) / Lz : 0.0; // 0 = near, 1 = far
        int g    = static_cast<int>(220.0 - 120.0 * t);
        file << "<circle cx=\"" << PX(sample[i].x, sample[i].z) << "\" cy=\"" << PY(sample[i].y, sample[i].z)
             << "\" r=\"" << sample[i].r * ppl << "\" fill=\"rgb(" << g << ',' << g << ",255)\" stroke=\"#333333\" "
             << "stroke-width=\"0.6\"/>\n";
      }
    }

    file << "</svg>\n";
  }

protected:
  double ran(double min_, double max_) {
    return min_ + (rand() / (double)RAND_MAX) * (max_ - min_);
  }

  // a ghost (periodic image) disk: dashed gray outline, light fill
  void drawGhostDisk_SVG(std::ofstream &file, double cx, double cy, double r) {
    file << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << r
         << "\" fill=\"#f0f0f5\" stroke=\"#999999\" stroke-width=\"0.6\" stroke-dasharray=\"3,2\"/>\n";
  }

  // exit point where the segment (p -> q) crosses the cell AABB, assuming p is inside the cell.
  void clipToCell_SVG(double p0, double p1, double p2, double q0, double q1, double q2, double &e0, double &e1,
                      double &e2) {
    double d0 = q0 - p0, d1 = q1 - p1, d2 = q2 - p2;
    double t = 1.0;
    if (d0 > 1e-12) {
      t = std::min(t, (xmax - p0) / d0);
    } else if (d0 < -1e-12) {
      t = std::min(t, (xmin - p0) / d0);
    }
    if (d1 > 1e-12) {
      t = std::min(t, (ymax - p1) / d1);
    } else if (d1 < -1e-12) {
      t = std::min(t, (ymin - p1) / d1);
    }
    if (d2 > 1e-12) {
      t = std::min(t, (zmax - p2) / d2);
    } else if (d2 < -1e-12) {
      t = std::min(t, (zmin - p2) / d2);
    }
    e0 = p0 + t * d0;
    e1 = p1 + t * d1;
    e2 = p2 + t * d2;
  }
};

#endif /* end of include guard: GEOPACK3D_HPP */

#if 0

int main(int argc, char const *argv[]) {
  
  GeoPack3D GP(0.8, 3, 100, 0, 100, 0, 100, 0, 100);
  GP.seedTime();
  GP.exec();
  GP.save("packing.txt");
  
  return 0;
}

#endif
