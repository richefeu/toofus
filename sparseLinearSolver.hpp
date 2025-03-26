#ifndef SPARSE_LINEAR_SOLVER_H
#define SPARSE_LINEAR_SOLVER_H

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

class SparseLinearSolver {
private:
  std::vector<double> values;   // Valeurs non nulles de la matrice
  std::vector<int> col_indices; // Indices de colonnes correspondants
  std::vector<int> row_ptrs;    // Pointeurs de début de chaque ligne
  int n;                        // Taille de la matrice (n x n)

public:
  SparseLinearSolver(int size) : n(size) { row_ptrs.resize(n + 1, 0); }

  void add_value(int row, int col, double value) {
    if (value != 0.0) {
      values.push_back(value);
      col_indices.push_back(col);
      row_ptrs[row + 1]++;
    }
  }

  void finalize() {
    for (int i = 1; i <= n; i++) {
      row_ptrs[i] += row_ptrs[i - 1];
    }
  }

  void multiply(const std::vector<double> &x, std::vector<double> &result) const {
    std::fill(result.begin(), result.end(), 0.0);
    for (int row = 0; row < n; row++) {
      for (int j = row_ptrs[row]; j < row_ptrs[row + 1]; j++) {
        result[row] += values[j] * x[col_indices[j]];
      }
    }
  }

  bool solve_CG(const std::vector<double> &b, std::vector<double> &x, double tol = 1e-6, int max_iter = 1000) const {
    std::vector<double> r(n), p(n), Ap(n);
    multiply(x, Ap);
    for (int i = 0; i < n; i++) {
      r[i] = b[i] - Ap[i];
    }
    p = r;
    double rs_old = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);

    for (int iter = 0; iter < max_iter; iter++) {
      multiply(p, Ap);
      double alpha = rs_old / std::inner_product(p.begin(), p.end(), Ap.begin(), 0.0);
      for (int i = 0; i < n; i++) {
        x[i] += alpha * p[i];
      }
      for (int i = 0; i < n; i++) {
        r[i] -= alpha * Ap[i];
      }

      double rs_new = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
      if (std::sqrt(rs_new) < tol) {
        return true;
      }

      double beta = rs_new / rs_old;
      for (int i = 0; i < n; i++) {
        p[i] = r[i] + beta * p[i];
      }

      rs_old = rs_new;
    }
    return false;
  }
};

#endif // SPARSE_LINEAR_SOLVER_H

#if 0

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

  std::vector<double> b = {1, 2, 3};
  std::vector<double> x(size, 0);

  if (solver.solve_CG(b, x)) {
    // bonne solution (2/9, 1/9, 13/9)
    std::cout << "Solution trouvée : ";
    for (double xi : x)
      std::cout << xi << " ";
    std::cout << std::endl;
  } else {
    std::cout << "Échec de convergence.\n";
  }

  return 0;
}

#endif
