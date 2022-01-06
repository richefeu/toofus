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

#ifndef STICK_BEAM_HPP
#define STICK_BEAM_HPP

#include <cmath>
#include <iostream>
#include <vector>

#include "vec3.hpp"

/*
 *   This small struct can be used for analysis bases on beam theory.
 *   The analysis is restricted to nodal forces that MUST be statically balanced.
 *
 *   y
 *   |
 *   |
 *   o-----x
 *   z
 *
 *   5 node  0            1       2         3                 4
 *           +------------+-------+---------+-----------------+
 *            o          e o     e o       e o               e
 *   4 elem        (0)        (1)      (2)           (3)
 *
 */
struct beam {
  double radius;

  // nodes
  std::vector<double> xpos;
  std::vector<vec3r> force;

  // elements (o means origine, e means extremity)
  std::vector<double> Nxo;
  std::vector<double> Nxe;
  std::vector<double> Vyo;
  std::vector<double> Vye;
  std::vector<double> Vzo;
  std::vector<double> Vze;
  std::vector<double> Mzo;
  std::vector<double> Mze;
  std::vector<double> Myo;
  std::vector<double> Mye;
  std::vector<double> sigmaSupo;
  std::vector<double> sigmaSupe;
  std::vector<double> sigmaInfo;
  std::vector<double> sigmaInfe;
  std::vector<double> tauMido;
  std::vector<double> tauMide;

  void clear() {
    xpos.clear();
    force.clear();
    Nxo.clear();
    Nxe.clear();
    Vyo.clear();
    Vye.clear();
    Vzo.clear();
    Vze.clear();
    Mzo.clear();
    Mze.clear();
    Myo.clear();
    Mye.clear();

    sigmaSupo.clear();
    sigmaSupe.clear();
    sigmaInfo.clear();
    sigmaInfe.clear();
    tauMido.clear();
    tauMide.clear();
  }

  void addNode(double x, double fx, double fy, double fz) {
    xpos.push_back(x);
    vec3r f(fx, fy, fz);
    force.push_back(f);

    if (xpos.size() >= 2) {
      size_t ilast = xpos.size() - 1;
      size_t ip = ilast;
      for (size_t p = 0; p < xpos.size(); ++p) {
        if (x <= xpos[p]) {
          ip = p;
          break;
        }
      }

      if (ip != ilast) {
        for (size_t p = xpos.size() - 1; p > ip; --p) {
          xpos[p] = xpos[p - 1];
          force[p] = force[p - 1];
        }
        xpos[ip] = x;
        force[ip] = f;
      }
    }
  }

  void connect() {
    if (xpos.empty()) {
      std::cout << "Cannot connect the nodes!\n";
      return;
    }

    size_t nbNodes = xpos.size();
    size_t nbElems = nbNodes - 1;
    Nxo.resize(nbElems, 0.0);
    Nxe.resize(nbElems, 0.0);
    Vyo.resize(nbElems, 0.0);
    Vye.resize(nbElems, 0.0);
    Vzo.resize(nbElems, 0.0);
    Vze.resize(nbElems, 0.0);
    Mzo.resize(nbElems, 0.0);
    Mze.resize(nbElems, 0.0);
    Myo.resize(nbElems, 0.0);
    Mye.resize(nbElems, 0.0);

    sigmaSupo.resize(nbElems, 0.0);
    sigmaSupe.resize(nbElems, 0.0);
    sigmaInfo.resize(nbElems, 0.0);
    sigmaInfe.resize(nbElems, 0.0);
    tauMido.resize(nbElems, 0.0);
    tauMide.resize(nbElems, 0.0);
  }

  // The static force balance is supposed to be statisfied
  // in this method
  void computeInternalActions() {
    size_t nbNodes = xpos.size();
    size_t nbElems = nbNodes - 1;

    for (size_t e = 0; e < nbElems; e++) {
      vec3r sumForceLeft;
      vec3r sumMomentLefto;
      vec3r sumMomentLefte;
      for (size_t nleft = 0; nleft <= e; nleft++) {
        sumForceLeft += force[nleft];
        sumMomentLefto += (xpos[e] - xpos[nleft]) * force[nleft];
        sumMomentLefte += (xpos[e + 1] - xpos[nleft]) * force[nleft];
      }
      Nxo[e] = Nxe[e] = -sumForceLeft.x;
      Vyo[e] = Vye[e] = -sumForceLeft.y;
      Vzo[e] = Vze[e] = -sumForceLeft.z;
      Mzo[e] = -sumMomentLefto.y;
      Mze[e] = -sumMomentLefte.y;
      Myo[e] = -sumMomentLefto.z;
      Mye[e] = -sumMomentLefte.z;
    }
  }

  // The method computeInternalActions needs to be called
  // before calling this method
  void computeNodeStress() {
    size_t nbNodes = xpos.size();
    size_t nbElems = nbNodes - 1;
    double S = M_PI * radius * radius;
    double d = 2.0 * radius;
    double I = M_PI * d * d * d * d / 64.0;
    for (size_t e = 0; e < nbElems; e++) {

      double Mo = sqrt(Myo[e] * Myo[e] + Mzo[e] * Mzo[e]);
      sigmaInfo[e] = Nxo[e] / S - (Mo / I) * (-radius);
      sigmaSupo[e] = Nxo[e] / S - (Mo / I) * (radius);
      double Vo = sqrt(Vyo[e] * Vyo[e] + Vzo[e] * Vzo[e]);
      tauMido[e] = (4.0 * Vo) / (3.0 * S);

      double Me = sqrt(Mye[e] * Mye[e] + Mze[e] * Mze[e]);
      sigmaInfe[e] = Nxe[e] / S - (Me / I) * (-radius);
      sigmaSupe[e] = Nxe[e] / S - (Me / I) * (radius);
      double Ve = sqrt(Vye[e] * Vye[e] + Vze[e] * Vze[e]);
      tauMide[e] = (4.0 * Ve) / (3.0 * S);
    }
  }

  // sigmaMax = yielding tensile normal stress
  // tauMax = yielding shear stress
  //
  // The cause (tensile or shear yielding) of the breakage is not saved,
  // The exact location within the section is not saved also
  std::vector<bool> getNodeBreakageStatus(double sigmaMax, double tauMax) {
    size_t nbNodes = xpos.size();
    size_t nbElems = nbNodes - 1;

    std::vector<bool> broken(nbNodes);
    for (size_t i = 0; i < nbNodes; i++)
      broken[i] = false;

    for (size_t e = 0; e < nbElems; e++) {
      if (sigmaInfo[e] > sigmaMax || sigmaSupo[e] > sigmaMax)
        broken[e] = true;
      if (sigmaInfe[e] > sigmaMax || sigmaSupe[e] > sigmaMax)
        broken[e + 1] = true;
      if (fabs(tauMido[e]) > tauMax)
        broken[e] = true;
      if (fabs(tauMide[e]) > tauMax)
        broken[e + 1] = true;
    }

    return broken;
  }

  std::vector<double> getBrokenParts(std::vector<bool> &brk) {
    size_t nbNodes = xpos.size();
    size_t nbElems = nbNodes - 1;

    std::vector<double> length;
    double x0 = xpos[0];
    for (size_t i = 1; i < nbNodes - 1; i++) {
      if (brk[i] == true) {
        length.push_back(xpos[i] - x0);
        x0 = xpos[i];
      }
    }
    length.push_back(xpos[nbNodes - 1] - x0);
    return length;
  }

  void print_sumF() {
    vec3r sumf;
    for (size_t i = 0; i < force.size(); i++) {
      sumf += force[i];
    }
    std::cout << "sum f = " << sumf << "\n";
  }

  void print_node_force() {

    for (size_t i = 0; i < xpos.size(); i++) {
      std::cout << "x =  " << xpos[i] << ", force = " << force[i] << "\n";
    }
  }

  void print() {
    size_t nbNodes = xpos.size();
    size_t nbElems = nbNodes - 1;

    for (size_t e = 0; e < nbElems; e++) {
      std::cout << "Elem " << e << "\n";
      std::cout << "     Nxo = " << Nxo[e] << "\t\t\tNxe = " << Nxe[e] << '\n';
      std::cout << "     Vyo = " << Vyo[e] << "\t\t\tVye = " << Vye[e] << '\n';
      std::cout << "     Vzo = " << Vzo[e] << "\t\t\tVze = " << Vze[e] << '\n';
      std::cout << "     Myo = " << Myo[e] << "\t\t\tMye = " << Mye[e] << '\n';
      std::cout << "     Mzo = " << Mzo[e] << "\t\t\tMze = " << Mze[e] << '\n';
      std::cout << "     sigmaInfo = " << sigmaInfo[e] << "\t\t\tsigmaInfe = " << sigmaInfe[e] << '\n';
      std::cout << "     sigmaSupo = " << sigmaSupo[e] << "\t\t\tsigmaSupe = " << sigmaSupe[e] << '\n';
      std::cout << "     tauMido = " << tauMido[e] << "\t\t\ttauMide = " << tauMide[e] << '\n';
    }
  }
};

#if 1
int main(int argc, char const *argv[]) {

  beam B;
  
  B.radius = 0.1;
  B.addNode(7.6, 0, -0.5, 0);
  B.addNode(6.8, 0, -2 + 5.5, 0);
  B.addNode(2.3, 0, -2, 0);
  B.addNode(5.3, 0, -2, 0);
  B.addNode(3.8, 0, -2, 0);
  B.addNode(0.8, 0, -2 + 5.5, 0);
  B.addNode(0, 0, -0.5, 0);
  
  
  /*
  B.radius = 0.1;
  B.addNode(0.0, 0, 5, 0);
  B.addNode(1.0, 0, -10, 0);
  B.addNode(2.0, 0, 5, 0);
  */
  
  B.print_sumF();
  B.print_node_force();

  B.connect();
  B.computeInternalActions();
  B.computeNodeStress();
  B.print();

  /*
  B.addNode(0.0, 0.0, 0.5, 0.0);
  B.addNode(0.5, 0.0, -1.0, 0.0);
  B.addNode(1.0, 0.0, 0.5, 0.0);
  B.connect();
  B.computeInternalActions();
  B.computeNodeStress();
  B.print();

  std::vector<bool> brk = B.getNodeBreakageStatus(320.0, 22.0);
  for (size_t i = 0; i < B.xpos.size(); i++) {
    std::cout << "node " << i << ", broken = " << brk[i] << '\n';
  }

  std::vector<double> len = B.getBrokenParts(brk);
  std::cout << "nb Parts = " << len.size() << '\n';
  for (size_t i = 0; i < len.size(); i++) {
    std::cout << "length " << i << " = " << len[i] << '\n';
  }
*/
  return 0;
}
#endif

#endif /* end of include guard */