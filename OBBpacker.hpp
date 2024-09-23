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

#ifndef OBB_PACKER_HPP
#define OBB_PACKER_HPP

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "OBB.hpp"

class OBBPacker {
public:
  enum ContainerType { CUBOID, CYLINDER };
  enum PackingDirection { X, Y, Z };

private:
  unsigned int seed = 143502;
  double LX = 20e-3;
  double LY = 20e-3;
  double LZ = 20e-3;

  double xOBB = 0.8e-3;
  double yOBB = 0.8e-3;
  double zOBB = 0.8e-3;

  double homothety_min = 1.0;
  double homothety_max = 1.0;

  size_t nb_obbs = 4000;
  size_t nb_trials_max = nb_obbs * 30;
  double SolidFractionTarget = 0.29;

  ContainerType containerType = CUBOID;
  PackingDirection packingDirection = Z;

  std::default_random_engine generator;
  std::uniform_real_distribution<double> random01;

public:
  OBBPacker() : generator(seed), random01(0.0, 1.0) {}

  void setOBBDimensions(double xobb, double yobb, double zobb) {
    xOBB = xobb;
    yOBB = yobb;
    zOBB = zobb;
  }

  void setHomothetyRange(double from, double to) {
    homothety_min = from;
    homothety_max = to;
  }

  void setNbOBBsTarget(size_t nb_obbs) {
    this->nb_obbs = nb_obbs;
    nb_trials_max = nb_obbs * 30;
  }

  void setNbTrialsMax(size_t nb_trials_max) { this->nb_trials_max = nb_trials_max; }

  void setSolidFractionTarget(double solid_fraction_target) { SolidFractionTarget = solid_fraction_target; }

  void setContainerLength(double t_LX, double t_LY, double t_LZ) {
    LX = t_LX;
    LY = t_LY;
    LZ = t_LZ;
  }

  void setContainerType(ContainerType container_type) { containerType = container_type; }
  void setPackingDirection(PackingDirection D) { packingDirection = D; }

  void save_pack(const char *fname, const std::vector<OBB> &obbs) {
    std::ofstream fog(fname, std::ios::out);
    if (fog) {
      fog.precision(5);
      fog << std::scientific;
      for (const auto &obb : obbs) {
        fog << obb.center << ' ' << obb.e[0] << ' ' << obb.e[1] << ' ' << obb.e[2] << ' ' << obb.extent << std::endl;
      }
    }
  }

  void pack(std::vector<OBB> &obbs) {

    double Vbox = 0.0;
    if (containerType == CUBOID) {
      Vbox = LX * LY * LZ;
    } else {
      switch (packingDirection) {
      case X: {
        double Sellipse = M_PI * LY * LZ / 4.0;
        Vbox = Sellipse * LX;
      } break;
      case Y: {
        double Sellipse = M_PI * LX * LZ / 4.0;
        Vbox = Sellipse * LY;
      } break;
      case Z: {
        double Sellipse = M_PI * LX * LY / 4.0;
        Vbox = Sellipse * LZ;
      } break;
      }
    }

    std::vector<vec3r> Ltab;
    double Vs = 0.0;
    for (size_t i = 0; i < nb_obbs; i++) {
      vec3r full_extent;
      double homothety = homothety_min + (homothety_max - homothety_min) * random01(generator);
      full_extent.x = xOBB * homothety;
      full_extent.y = yOBB * homothety;
      full_extent.z = zOBB * homothety;
      Ltab.push_back(full_extent);
      Vs += full_extent.x * full_extent.y * full_extent.z;
      if (Vs / Vbox >= SolidFractionTarget) {
        break;
      }
    }
    std::sort(Ltab.begin(), Ltab.end(),
              [](const vec3r &a, const vec3r &b) { return (a.x * a.y * a.z) > (b.x * b.y * b.z); });

    obbs.clear();

    OBB obb1;
    obb1.extent = 0.5 * Ltab[0];
    quat ori;
    ori.randomize();
    obb1.rotate(ori);
    obb1.center.x = 0.5 * LX;
    obb1.center.y = 0.5 * LY;
    obb1.center.z = 0.5 * LZ;
    obbs.push_back(obb1);
    double sumV = Ltab[0].x * Ltab[0].y * Ltab[0].z;
    size_t cumul_trials = 1;

    std::cerr << "OBBPacker is Running...\n" << std::flush;
    for (size_t i = 1; i < Ltab.size(); i++) {

      OBB obbi;
      double Vi = Ltab[i].x * Ltab[i].y * Ltab[i].z;

      size_t nb_trials = 1;
    retry:

      ori.randomize();
      obbi.rotate(ori);

      double distBound = 0.5 * Ltab[i].length();
      if (containerType == CUBOID) {

        obbi.center.x = distBound + (LX - 2.0 * distBound) * random01(generator);
        obbi.center.y = distBound + (LY - 2.0 * distBound) * random01(generator);
        obbi.center.z = distBound + (LZ - 2.0 * distBound) * random01(generator);

      } else {
        double angle = 2 * M_PI * random01(generator);
        switch (packingDirection) {
        case X: {
          double a = LY / 2.0;
          double b = LZ / 2.0;
          double r = a * b / sqrt(b * b * cos(angle) * cos(angle) + a * a * sin(angle) * sin(angle));
          obbi.center.x = distBound + (LX - 2 * distBound) * random01(generator);
          obbi.center.y = 0.5 * LY + (r - distBound) * cos(angle) * random01(generator);
          obbi.center.z = 0.5 * LZ + (r - distBound) * sin(angle) * random01(generator);
        } break;
        case Y: {
          double a = LX / 2.0;
          double b = LZ / 2.0;
          double r = a * b / sqrt(b * b * cos(angle) * cos(angle) + a * a * sin(angle) * sin(angle));
          obbi.center.x = 0.5 * LX + (r - distBound) * cos(angle) * random01(generator);
          obbi.center.y = distBound + (LY - 2 * distBound) * random01(generator);
          obbi.center.z = 0.5 * LZ + (r - distBound) * sin(angle) * random01(generator);
        } break;
        case Z: {
          double a = LX / 2.0;
          double b = LY / 2.0;
          double r = a * b / sqrt(b * b * cos(angle) * cos(angle) + a * a * sin(angle) * sin(angle));
          obbi.center.x = 0.5 * LX + (r - distBound) * cos(angle) * random01(generator);
          obbi.center.y = 0.5 * LY + (r - distBound) * sin(angle) * random01(generator);
          obbi.center.z = distBound + (LZ - 2 * distBound) * random01(generator);
        } break;
        }
      }
      obbi.extent = 0.5 * Ltab[i];

      for (long k = 0; k < obbs.size(); k++) {
        if (obbi.intersect(obbs[k])) {
          nb_trials++;
          if (nb_trials >= nb_trials_max) {
            std::cerr << "Maximum number of trials reached !" << std::endl;
            std::cerr << "Number of shapes packed = " << i << "\n";
            std::cerr << "Solid fraction = " << sumV / Vbox << "\n";
            return;
          }
          goto retry;
        }
      }

      obbs.push_back(obbi);
      cumul_trials += nb_trials;
      sumV += Vi;
    }

    std::cerr << "Okay! " << obbs.size() << " OBBs have been packed" << std::endl;
    std::cerr << "Solid fraction = " << sumV / Vbox << std::endl;
  }

  void sort(std::vector<OBB> &obbs) {
    std::sort(obbs.begin(), obbs.end(), [&](const OBB &a, const OBB &b) {
      switch (packingDirection) {
      case X:
        return a.center.x < b.center.x;
      case Y:
        return a.center.y < b.center.y;
      case Z:
        return a.center.z < b.center.z;
      }
      return false;
    });
  }
};

#endif // OBB_PACKER_HPP

#if 0

int main() {
  OBBPacker packer;
  packer.setOBBDimensions(1e-3, 1e-3, 1e-3);
  packer.setHomothetyRange(1.0, 1.0);
  packer.setNbOBBsTarget(8000);
  packer.setSolidFractionTarget(0.7);
  packer.setContainerLength(20.0e-3, 20.0e-3, 30.0e-3);
  packer.setContainerType(OBBPacker::CUBOID);
  packer.setPackingDirection(OBBPacker::Z);
  std::vector<OBB> obbs;
  packer.pack(obbs);
  packer.sort(obbs);
  packer.save_pack("obbs.txt", obbs);
  return 0;
}

#endif