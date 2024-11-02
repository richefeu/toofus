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

#ifndef DATATABLE_HPP
#define DATATABLE_HPP

// ===================================================================
// Classes to store and manage parameters of groups and between groups
// ===================================================================

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

/**
Usage example:
@code
DataTable dt;

dt.set("mu", 0, 0, 0.7);
dt.set("kn", 0, 1, 1000);
dt.set("kn", 1, 2, 200.8);

dt.write(std::cout);
@endcode
*/

/// @brief A class for definition and storage of data for interaction between each group.
/// The class is designed for a rapid access to the data with the method get
class DataTable {
public:
  size_t ngroup; ///< number of groups
  std::vector<std::vector<std::vector<double>>> tables; ///< tables for storing data
  std::map<std::string, size_t> data_id; ///< map for quick access to the data (name to id)
  std::set<std::tuple<size_t, size_t, size_t>> defined; ///< set of defined parameters

public:
  DataTable() { set_ngroup(1); } ///< constructor

  /// clear all data
  void clear() {
    for (size_t t = 0; t < tables.size(); ++t) {
      for (size_t i = 0; i < tables[t].size(); ++i)
        tables[t][i].clear();
      tables[t].clear();
    }
    tables.clear();
    data_id.clear();
    set_ngroup(1);
  }

  /// set the number of groups
  void set_ngroup(size_t n) {
    ngroup = n;
    for (size_t t = 0; t < tables.size(); ++t) {
      tables[t].resize(ngroup);
      for (size_t i = 0; i < tables[t].size(); ++i)
        tables[t][i].resize(ngroup, 0.0);
    }
  }

  /// get the number of groups
  size_t get_ngroup() const { return ngroup; }

  /// check if a parameter exists
  bool exists(const std::string &name) const {
    std::map<std::string, size_t>::const_iterator ip = data_id.find(name);
    if (ip == data_id.end()) {
      return false;
    }
    return true;
  }

  /// add a parameter and return its id
  size_t add(const std::string &name) {
    std::map<std::string, size_t>::const_iterator ip = data_id.find(name);
    if (ip != data_id.end()) {
      return (size_t)ip->second;
    }

    data_id[name] = tables.size();
    std::vector<std::vector<double>> table;
    table.resize(ngroup);
    for (size_t i = 0; i < table.size(); i++)
      table[i].resize(ngroup, 0.0);
    tables.push_back(table);
    return (tables.size() - 1);
  }

  /// get the id of a parameter
  size_t get_id(const std::string &name) const {
    std::map<std::string, size_t>::const_iterator ip = data_id.find(name);
    if (ip != data_id.end()) {
      return (size_t)ip->second;
    } else {
      return 0;
    }
  }

  /// get the value of a parameter
  double get(size_t id, size_t g1, size_t g2) const { return tables[id][g1][g2]; }

  /// check if a parameter is defined
  bool isDefined(size_t id, size_t g1, size_t g2) const {
    return (defined.find(std::tuple<size_t, size_t, size_t>(id, g1, g2)) != defined.end());
  }

  /// set the value of a parameter
  void set(size_t id, size_t g1, size_t g2, double val) {
    if (g1 >= ngroup)
      set_ngroup(g1 + 1);
    if (g2 >= ngroup)
      set_ngroup(g2 + 1);
    if (id < tables.size()) {
      tables[id][g1][g2] = val;
      tables[id][g2][g1] = val;
      defined.insert(std::tuple<size_t, size_t, size_t>(id, g1, g2));
    }
  }

  /// A self-add method to set a parameter
  size_t set(const std::string &name, size_t g1, size_t g2, double val) {
    size_t id = add(name);
    set(id, g1, g2, val);
    return id;
  }
};

#endif /* end of include guard: DATATABLE_HPP */

