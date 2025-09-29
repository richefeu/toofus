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

#ifndef STLMESH_HPP
#define STLMESH_HPP

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>

#include <map>
#include <set>
#include <vector>

#include <cstdint>
#include <sstream>
#include <string>

#define ZERO 1.0e-12

typedef double real;
typedef unsigned int uint;
typedef float real32;
typedef double real64;

#include "vec3.hpp"

class STLmesh {
public:
  // point
  struct point {
    double x, y, z; // coordinates
    point() : x(0.0), y(0.0), z(0.0) {}
    point(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
  };

  // triangle
  struct triangle {
    int i0, i1, i2; // ID-numbers of points in the mesh
    point normal;
    triangle() : i0(0), i1(0), i2(0), normal() {}
    triangle(int I0, int I1, int I2) : i0(I0), i1(I1), i2(I2), normal() {}
    triangle(int I0, int I1, int I2, double Nx, double Ny, double Nz) : i0(I0), i1(I1), i2(I2), normal(Nx, Ny, Nz) {}
  };

  // edge
  struct edge {
    int i0, i1; // ID-numbers of points in the mesh
    edge() : i0(0), i1(0) {}
    edge(int I0, int I1) : i0(I0), i1(I1) {}
    /**
     * @brief Compares two edges by comparing their node indices.
     *
     * @param rhs Edge to compare to.
     * @return true if the edge is smaller than rhs, false otherwise.
     */
    bool operator<(const edge &rhs) const {
      if (i0 < rhs.i0) return true;
      if ((i0 == rhs.i0) && (i1 < rhs.i1)) return true;
      return false;
    }
  };

  std::vector<point> points;
  std::vector<triangle> triangles;
  std::vector<edge> edges;

  // Ctor
  STLmesh() {}

  /**
   * @brief Clears the mesh.
   *
   * This method clears all the data structure of the mesh, i.e.
   * points, triangles and edges. It is mainly used to reset the
   * mesh to its initial state.
   */
  void clear() {
    points.clear();
    triangles.clear();
    edges.clear();
  }

  /**
   * @brief Removes duplicated edges in the mesh.
   *
   * This method cleans the edge list by removing duplicate edges.
   * An edge is considered duplicated if it is represented more than
   * once in the mesh with the same pair of node indices. The method
   * ensures that each edge is unique by using a set to store edges
   * as unique pairs of node indices. Duplicates are identified and
   * removed, and the total number of duplicated edges removed is
   * printed to the console.
   */
  void cleanEdges() {
    std::set<std::pair<size_t, size_t>> edgeSet;
    for (size_t e = 0; e < edges.size(); e++) {
      std::pair<size_t, size_t> P;
      if (edges[e].i0 < edges[e].i1) {
        P.first  = edges[e].i0;
        P.second = edges[e].i1;
      } else {
        P.first  = edges[e].i1;
        P.second = edges[e].i0;
      }
      edgeSet.insert(P); // That way, duplication is not allowed
    }

    size_t nbePrev = edges.size();
    size_t nbeNew  = edgeSet.size();

    edges.clear();
    for (auto it : edgeSet) { edges.push_back(edge(it.first, it.second)); }

    std::cout << nbePrev - nbeNew << " duplicated edges have been removed\n";
  }

  /**
   * @brief Resets the normals of each triangle in the mesh.
   *
   * This method resets the normal of each triangle in the mesh to
   * the cross product of the two vectors defined by the first and
   * second points of the triangle, and the first and third points
   * of the triangle. The resulting normal is then normalized and
   * stored in the triangle structure.
   */
  void resetNormals() {
    for (size_t t = 0; t < triangles.size(); t++) {
      vec3r a(points[triangles[t].i1].x - points[triangles[t].i0].x,
              points[triangles[t].i1].y - points[triangles[t].i0].y,
              points[triangles[t].i1].z - points[triangles[t].i0].z);
      vec3r b(points[triangles[t].i2].x - points[triangles[t].i0].x,
              points[triangles[t].i2].y - points[triangles[t].i0].y,
              points[triangles[t].i2].z - points[triangles[t].i0].z);
      vec3r n = cross(a, b);
      n.normalize();
      triangles[t].normal.x = n.x;
      triangles[t].normal.y = n.y;
      triangles[t].normal.z = n.z;
    }
  }

  /**
   * @brief Saves the mesh to a binary STL file.
   *
   * This method writes the mesh data, including normals and vertices of each
   * triangle, to a specified file in binary STL format. The STL file
   * begins with an 80-byte header followed by the number of triangles.
   * Each triangle's normal and vertices are written as 32-bit floating-point
   * numbers, followed by an attribute byte count set to zero.
   *
   * @param name The name of the file to save the mesh to.
   */
  void save(const char *name) {
    std::ofstream file(name, std::ofstream::binary | std::ofstream::out);
    if (!file) {
      std::cerr << "Cannot open " << name << std::endl;
      return;
    }

    uint8_t head[80];
    std::strncpy((char *)&head, "STLmesh", sizeof(head) - 1);
    file.write((char *)&head, 80 * sizeof(uint8_t));

    uint32_t nb_triangles = (uint32_t)triangles.size();
    file.write((char *)&nb_triangles, sizeof(uint32_t));

    vec3<real32> vnormal, vertex1, vertex2, vertex3;
    uint16_t attByte = 0;

    for (size_t i = 0; i < triangles.size(); i++) {
      vnormal.x = (real32)(triangles[i].normal.x);
      vnormal.y = (real32)(triangles[i].normal.y);
      vnormal.z = (real32)(triangles[i].normal.z);
      vertex1.x = (real32)(points[triangles[i].i0].x);
      vertex1.y = (real32)(points[triangles[i].i0].y);
      vertex1.z = (real32)(points[triangles[i].i0].z);
      vertex2.x = (real32)(points[triangles[i].i1].x);
      vertex2.y = (real32)(points[triangles[i].i1].y);
      vertex2.z = (real32)(points[triangles[i].i1].z);
      vertex3.x = (real32)(points[triangles[i].i2].x);
      vertex3.y = (real32)(points[triangles[i].i2].y);
      vertex3.z = (real32)(points[triangles[i].i2].z);

      file.write((char *)&vnormal, sizeof(real32) * 3);
      file.write((char *)&vertex1, sizeof(real32) * 3);
      file.write((char *)&vertex2, sizeof(real32) * 3);
      file.write((char *)&vertex3, sizeof(real32) * 3);
      file.write((char *)&attByte, sizeof(uint16_t));
    }
    file.close();
  }

  /**
   * @brief Read a mesh file (stl BINARY format)
   *
   * @param name name of the file to read
   *
   * @return none
   *
   * This function reads a mesh file in stl BINARY format and initialize
   * the STLmesh object with the read data.
   *
   * A STL file is a binary file that contains a header of 80 bytes,
   * a 32 bit integer giving the number of triangles, and for each triangle:
   * - a 32 bit floating point number for the normal (x, y, z)
   * - three 32 bit floating point numbers for each vertex (x, y, z)
   * - a 16 bit integer for the attribute byte count
   *
   * The STLmesh object is initialized with the read data and the edges
   * are computed and stored in the object.
   *
   * @note The function does not check if the file is a valid STL file.
   *       It is the responsability of the caller to ensure that the file
   *       is valid before calling this function.
   */
  void read(const char *name) {

    // UINT8[80]  Header
    // UINT32     Number of triangles
    //
    // foreach triangle
    // REAL32[3]  Normal vector
    // REAL32[3]  Vertex 1
    // REAL32[3]  Vertex 2
    // REAL32[3]  Vertex 3
    // UINT16     Attribute byte count
    // end

    std::ifstream file(name, std::ifstream::binary | std::ifstream::in);
    if (!file) {
      std::cerr << "Cannot read " << name << std::endl;
      return;
    }

    // std::cerr << "Read mesh file (stl BINARY format)... " << std::flush;

    uint8_t head[80];
    file.read((char *)&head, sizeof(uint8_t) * 80);
    std::cout << head << '\n';

    uint32_t nb_triangles;
    file.read((char *)&nb_triangles, sizeof(uint32_t));

    vec3<real32> P1, P2, P3;
    triangle T;
    std::map<vec3<real32>, size_t, std::less<vec3<real32>>> map_points;
    std::map<vec3<real32>, size_t, std::less<vec3<real32>>>::iterator itMap;
    size_t ivertex = 0;
    uint16_t attByte;
    int swp;

    triangles.clear();
    points.clear();

    vec3<real32> vnormal, vertex1, vertex2, vertex3;

    for (uint32_t t = 0; t < nb_triangles; t++) {
      file.read((char *)&vnormal, sizeof(real32) * 3);
      file.read((char *)&vertex1, sizeof(real32) * 3);
      file.read((char *)&vertex2, sizeof(real32) * 3);
      file.read((char *)&vertex3, sizeof(real32) * 3);
      file.read((char *)&attByte, sizeof(uint16_t));

      T.normal.x = (double)(vnormal.x);
      T.normal.y = (double)(vnormal.y);
      T.normal.z = (double)(vnormal.z);

      // std::cout << vnormal.x << ' ' << vnormal.y << ' ' << vnormal.z << " | ";
      //  FIXME: it seems that normal are finally wrong
      //         (don't know why)

      itMap = map_points.find(vertex1);
      if (itMap != map_points.end()) {
        T.i0 = itMap->second;
      } else {
        map_points[vertex1] = ivertex;
        T.i0                = ivertex;
        ivertex++;
      }

      itMap = map_points.find(vertex2);
      if (itMap != map_points.end()) {
        T.i1 = itMap->second;
      } else {
        map_points[vertex2] = ivertex;
        T.i1                = ivertex;
        ivertex++;
      }

      itMap = map_points.find(vertex3);
      if (itMap != map_points.end()) {
        T.i2 = itMap->second;
      } else {
        map_points[vertex3] = ivertex;
        T.i2                = ivertex;
        ivertex++;
      }

      if (T.i0 > T.i1) {
        swp  = T.i0;
        T.i0 = T.i1;
        T.i1 = swp;
      }
      if (T.i1 > T.i2) {
        swp  = T.i1;
        T.i1 = T.i2;
        T.i2 = swp;
      }
      if (T.i0 > T.i1) {
        swp  = T.i0;
        T.i0 = T.i1;
        T.i1 = swp;
      }

      triangles.push_back(T);
    }

    std::map<size_t, vec3<real32>> map_id;
    for (itMap = map_points.begin(); itMap != map_points.end(); ++itMap) { map_id[itMap->second] = itMap->first; }
    map_points.clear();

    point Pt;
    for (std::map<size_t, vec3<real32>>::iterator it = map_id.begin(); it != map_id.end(); ++it) {
      Pt.x = (real)((it->second).x);
      Pt.y = (real)((it->second).y);
      Pt.z = (real)((it->second).z);
      points.push_back(Pt);
    }

    // Build the set of edges
    edges.clear();
    std::set<edge> tmp_edges;
    edge E;
    for (size_t t = 0; t < triangles.size(); t++) {
      E.i0 = triangles[t].i0;
      E.i1 = triangles[t].i1;
      tmp_edges.insert(E);
      E.i0 = triangles[t].i1;
      E.i1 = triangles[t].i2;
      tmp_edges.insert(E);
      E.i0 = triangles[t].i2;
      E.i1 = triangles[t].i0;
      tmp_edges.insert(E);
    }

    for (std::set<edge>::iterator it = tmp_edges.begin(); it != tmp_edges.end(); ++it) { edges.push_back(*it); }

    // std::cout << "done." << std::endl;
    std::cout << "  Number of vertices  " << points.size() << std::endl;
    std::cout << "  Number of edges     " << edges.size() << std::endl;
    std::cout << "  Number of triangles " << triangles.size() << std::endl;
  }
};

#endif /* end of include guard: STLMESH_HPP */

#if 0
#include <iostream>

int main(int argc, char const *argv[]) {
  STLmesh M;

  M.read("examples/cubo5mm.stl");
  M.cleanEdges();
  M.resetNormals();
  M.save("examples/cubo.stl");

  //STLmesh M2;
  //M2.read("examples/cubo.stl");
  //M2.save("examples/cubo2.stl");

  return 0;
}
#endif
