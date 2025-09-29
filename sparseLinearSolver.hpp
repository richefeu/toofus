#ifndef SPARSE_LINEAR_SOLVER_H
#define SPARSE_LINEAR_SOLVER_H

#include <cmath>
#include <iostream>
#include <numeric>
#include <unordered_set>
#include <vector>

class SparseLinearSolver {
private:
  std::vector<double> values;   // Non-zero values of the matrix
  std::vector<int> col_indices; // Corresponding column indices
  std::vector<int> row_ptrs;    // Row pointers
  int n;                        // Matrix size (n x n)

  std::unordered_set<int> hashes; // Set to store unique (row, col) hashes

public:
  SparseLinearSolver(int size) : n(size) {
    row_ptrs.resize(n + 1, 0);
  }

  // Method to list all values in the sparse matrix
  void list_values() const {
    for (int row = 0; row < n; row++) {
      for (int j = row_ptrs[row]; j < row_ptrs[row + 1]; j++) {
        std::cout << "(" << row << ", " << col_indices[j] << "): " << values[j] << std::endl;
      }
    }
  }

  // Add a value to the matrix at the specified (row, col)
  void add_value(int row, int col, double value) {
    if (value == 0.0) { return; }

    int hash = row * n + col; // Create a unique hash for each position (row, col)

    // Check if this hash already exists in the set
    if (hashes.contains(hash)) {
      // If the hash exists, we need to find the index where the (row, col) is stored
      for (int i = row_ptrs[row]; i < row_ptrs[row + 1]; ++i) {
        if (col_indices[i] == col) {
          // Found the existing position, add the new value to the current value
          values[i] += value;
          return; // No need to add a new entry
        }
      }
    } else {
      hashes.insert(hash);
      values.push_back(value);
      col_indices.push_back(col);
      row_ptrs[row + 1]++; // Increment the number of elements in this row
    }
  }

  // Finalize the row pointers to point to the correct ends of each row
  void finalize() {
    // Now accumulate the row_ptrs array
    for (int i = 1; i <= n; i++) { row_ptrs[i] += row_ptrs[i - 1]; }
  }

  // Matrix-vector multiplication: result = A * x
  void multiply(const std::vector<double> &x, std::vector<double> &result) const {
    std::fill(result.begin(), result.end(), 0.0);
    for (int row = 0; row < n; row++) {
      for (int j = row_ptrs[row]; j < row_ptrs[row + 1]; j++) { result[row] += values[j] * x[col_indices[j]]; }
    }
  }

  // Conjugate Gradient Solver (CG) for solving Ax = b
  bool solve_CG(const std::vector<double> &b, std::vector<double> &x, double tol = 1e-6, int max_iter = 1000) const {
    std::vector<double> r(n), p(n), Ap(n);
    multiply(x, Ap);
    for (int i = 0; i < n; i++) { r[i] = b[i] - Ap[i]; }
    p             = r;
    double rs_old = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);

    for (int iter = 0; iter < max_iter; iter++) {
      multiply(p, Ap);
      double alpha = rs_old / std::inner_product(p.begin(), p.end(), Ap.begin(), 0.0);
      for (int i = 0; i < n; i++) { x[i] += alpha * p[i]; }
      for (int i = 0; i < n; i++) { r[i] -= alpha * Ap[i]; }

      double rs_new = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
      if (std::sqrt(rs_new) < tol) { return true; }

      double beta = rs_new / rs_old;
      for (int i = 0; i < n; i++) { p[i] = r[i] + beta * p[i]; }

      rs_old = rs_new;
    }
    return false;
  }
};

#endif // SPARSE_LINEAR_SOLVER_H

#if 1

#include <iostream>

int main() {
  int size = 3;
  SparseLinearSolver solver(size);

  solver.add_value(0, 0, 4);
  solver.add_value(0, 1, 1);
  solver.add_value(1, 0, 1);
  solver.add_value(1, 1, 3);
  solver.add_value(1, 2, 1);
  solver.add_value(2, 1, 1);
  solver.add_value(2, 2, 2);

  solver.finalize();
  solver.list_values();

  std::vector<double> b = {1, 2, 3};
  std::vector<double> x(size, 0);

  if (solver.solve_CG(b, x)) {
    std::cout << "Solution trouvée : ";
    for (double xi : x) std::cout << xi << " ";
    std::cout << std::endl;
  } else {
    std::cout << "Échec de convergence.\n";
  }

  return 0;
}

#endif
