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

  /*!
   * \brief Clear all the vectors (nodes and elements)
   *
   * The following vectors are cleared:
   * - xpos
   * - force
   * - Nxo
   * - Nxe
   * - Vyo
   * - Vye
   * - Vzo
   * - Vze
   * - Mzo
   * - Mze
   * - Myo
   * - Mye
   * - sigmaSupo
   * - sigmaSupe
   * - sigmaInfo
   * - sigmaInfe
   * - tauMido
   * - tauMide
   */
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

  /**
   * Adds a node to the stick beam model.
   *
   * @param[in] x the position of the node in the direction of the beam.
   * @param[in] fx the x component of the force at the node.
   * @param[in] fy the y component of the force at the node.
   * @param[in] fz the z component of the force at the node.
   */
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

  /**
   * Connects the nodes by assigning the appropriate values to the
   * member variables Nxo, Nxe, Vyo, Vye, Vzo, Vze, Mzo, Mze, Myo,
   * Mye, sigmaSupo, sigmaSupe, sigmaInfo, sigmaInfe, tauMido, and
   * tauMide.
   *
   * @pre The nodes must be sorted by increasing x position.
   * @pre The nodes must be added using the addNode() member function.
   * @post The nodes are connected.
   */
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

  /**
   * Computes the internal actions of the beam.
   *
   * The method computes the internal actions Nxo, Nxe, Vyo, Vye, Vzo,
   * Vze, Mzo, Mze, Myo, and Mye of the beam elements. The internal
   * actions are computed by summing the external forces and moments
   * on the left of each node and on the right of each node. The
   * results are stored in the member variables Nxo, Nxe, Vyo, Vye,
   * Vzo, Vze, Mzo, Mze, Myo, and Mye.
   *
   * @pre The nodes must be sorted by increasing x position.
   * @pre The nodes must be added using the addNode() member function.
   * @pre The connect() member function must be called before this
   * member function.
   * @post The internal actions are computed.
   */
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

  /**
   * Computes the stress at the nodes of the beam.
   *
   * The method computes the stress at the nodes of the beam. The
   * stress is computed using the internal actions computed by the
   * computeInternalActions() member function. The results are stored
   * in the member variables sigmaInfo, sigmaInfe, sigmaSupo, and
   * sigmaSupe.
   *
   * @pre The nodes must be sorted by increasing x position.
   * @pre The nodes must be added using the addNode() member function.
   * @pre The connect() member function must be called before this
   * member function.
   * @pre The computeInternalActions() member function must be called
   * before this member function.
   * @post The stress at the nodes is computed.
   */
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

  /**
   * Determines the breakage status of nodes based on stress limits.
   *
   * This function evaluates each element in the beam to determine if it
   * exceeds the specified maximum tensile normal stress (`sigmaMax`) or
   * shear stress (`tauMax`). It returns a vector indicating which nodes
   * are considered broken.
   *
   * @param sigmaMax The maximum allowable tensile normal stress.
   * @param tauMax The maximum allowable shear stress.
   * @return A vector of boolean values indicating the breakage status
   *         of each node. `true` indicates the node is broken.
   *
   * @pre The stress at each node must be computed before calling this
   * function.
   * @post The vector returned represents the breakage status of nodes
   *       based on the given stress limits.
   */
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

  /**
   * Computes the lengths of the broken parts of the beam.
   *
   * This function takes the breakage status of nodes as input and returns
   * a vector of the lengths of the broken parts of the beam.
   *
   * @param brk A vector of boolean values indicating the breakage status
   *            of each node. `true` indicates the node is broken.
   * @return A vector of the lengths of the broken parts of the beam.
   *
   * @pre The `getNodeBreakageStatus` function must be called before this
   * function to compute the breakage status of nodes.
   * @post The vector returned represents the lengths of the broken parts
   *       of the beam.
   */
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

  /**
   * Prints the external forces applied to each node of the beam.
   *
   * This function simply prints the external forces applied to each node of
   * the beam. The forces are printed at the same x-coordinate as the node.
   *
   * @post The external forces at each node of the beam are printed.
   */
  void print_node_force() {

    for (size_t i = 0; i < xpos.size(); i++) {
      std::cout << "x =  " << xpos[i] << ", force = " << force[i] << "\n";
    }
  }

  /**
   * Prints the values of the internal forces and moments at the nodes of each
   * element of the beam.
   *
   * This function prints the values of the internal forces and moments at the
   * nodes of each element of the beam. The values are printed at the same
   * x-coordinate as the node.
   *
   * @post The values of the internal forces and moments at the nodes of each
   * element of the beam are printed.
   */
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

#if 0
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