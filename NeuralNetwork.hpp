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
//
// This is a stand-alone header file that implement a Neural Network
// with a single hidden layer (Non-deep Neural Network)

#ifndef NN_HPP
#define NN_HPP

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

typedef std::vector<double> cvector;
typedef std::vector<cvector> cmatrix;

void cvector_add_inplace(cvector &m, cvector &b) {
  if (m.size() != b.size()) { std::cerr << "@cvector_add_inplace, m and b must have the same size\n"; }
  for (size_t i = 0; i < m.size(); i++) { m[i] += b[i]; }
}

cvector cvector_substract(cvector &a, cvector &b) {
  if (a.size() != b.size()) { std::cerr << "@cvector_substract, a and b must have the same size\n"; }
  cvector result(a.size());
  for (size_t i = 0; i < a.size(); i++) { result[i] = a[i] - b[i]; }
  return result;
}

cmatrix cmatrix_create(size_t rows, size_t cols) {
  cmatrix M;
  M.resize(rows);
  for (size_t i = 0; i < rows; i++) {
    M[i].resize(cols);
    for (size_t j = 0; j < cols; j++) { M[i][j] = 0.0; }
  }
  return M;
}

void cmatrix_getSizes(cmatrix &m, size_t &rows, size_t &cols) {
  rows = m.size();
  cols = m[0].size();
}

void cmatrix_randomize_gaussian_inplace(cmatrix &m, double stddev = 1.0) {
  size_t rows, cols;
  cmatrix_getSizes(m, rows, cols);

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0, stddev);

  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) { m[i][j] = distribution(generator); }
  }
}

void cmatrix_randomize_inplace(cmatrix &m, double min = -1.0, double max = 1.0) {
  size_t rows, cols;
  cmatrix_getSizes(m, rows, cols);

  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) { m[i][j] = min + (rand() / (double)RAND_MAX) * (max - min); }
  }
}

void cvector_randomize_gaussian_inplace(cvector &v, double stddev = 1.0) {
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0, stddev);
  for (size_t i = 0; i < v.size(); i++) { v[i] = distribution(generator); }
}

void cvector_randomize_inplace(cvector &v, double min = -1.0, double max = 1.0) {
  for (size_t i = 0; i < v.size(); i++) { v[i] = min + (rand() / (double)RAND_MAX) * (max - min); }
}

// example with lambda: cvector_map_inplace(m, [](double x) -> double { return 2.0 * x; });
void cvector_map_inplace(cvector &m, std::function<double(double)> func) {
  for (size_t i = 0; i < m.size(); i++) { m[i] = func(m[i]); }
}

cvector cmatrix_mult(cmatrix &a, cvector &b) {
  size_t a_rows, a_cols;
  cmatrix_getSizes(a, a_rows, a_cols);

  size_t b_rows = b.size();

  if (a_cols != b_rows) {
    std::cerr << "@cmatrix_mult(a, b): number of columns of a must equals number of rows of b\n";
  }

  cvector result(a_rows);
  for (size_t i = 0; i < a_rows; i++) {
    double sum = 0.0;
    for (size_t k = 0; k < a_cols; k++) { sum += a[i][k] * b[k]; }
    result[i] = sum;
  }
  return result;
}

void cmatrix_add_inplace(cmatrix &a, cmatrix &b) {
  size_t a_rows, a_cols;
  cmatrix_getSizes(a, a_rows, a_cols);

  size_t b_rows, b_cols;
  cmatrix_getSizes(b, b_rows, b_cols);

  if (a_rows != b_rows || a_cols != b_cols) {
    std::cerr << "@cmatrix_add_inplace(a, b): a and b must have the same sizes\n";
  }

  for (size_t i = 0; i < a_rows; i++) {
    for (size_t j = 0; j < a_cols; j++) { a[i][j] += b[i][j]; }
  }
}

void print(cmatrix &m) {
  size_t rows, cols;
  cmatrix_getSizes(m, rows, cols);
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) { std::cout << m[i][j] << ' '; }
    std::cout << '\n';
  }
  std::cout << '\n';
}

void print(cvector &v) {
  for (size_t i = 0; i < v.size(); i++) { std::cout << v[i] << '\n'; }
  std::cout << '\n';
}

// ==================== ACTIVATORS ============================

// SIGMOID
int SIGMOID_ACTIVATOR_ID = 0;
double sigmoid(double x) {
  return (1.0 / (1.0 + exp(-x)));
}
double dsigmoid(double s) {
  return (s * (1.0 - s));
}
double invSigmoid(double y) {
  return (log(y / (1.0 - y)));
}

// TANH
int TANH_ACTIVATOR_ID = 1;
double Tanh(double x) {
  return tanh(x);
}
double dTanh(double s) {
  double th = tanh(s);
  return (1.0 - th * th);
}
double invTanh(double y) {
  return 0.5 * log((1.0 + y) / (1.0 - y));
}

// BINARY STEP
int BINARY_STEP_ACTIVATOR_ID = 2;
double binary(double x) {
  if (fabs(x) < 1.0) return 0.5 * (x + 1.0);
  if (x < 0.0) return 0.0;
  return 1.0;
}
double dBinary(double x) {
  if (fabs(x) < 1.0) return 0.5;
  return 0.0;
}
double invBinary(double x) {
  return 0.0;
} // TODO

// IDENTITY
int IDENTITY_ACTIVATOR_ID = 3;
double identity(double x) {
  return x;
}
double dIdentity(double s) {
  return 1.0;
}
double invIdentity(double y) {
  return y;
}

// BINARY STEP
int RELU_ACTIVATOR_ID = 4;
double reLU(double x) {
  if (x < 0.0) return 0.0;
  return x;
}
double dReLU(double x) {
  if (x < 0.0) return 0.0;
  return 1.0;
}
double invReLU(double y) {
  if (y < 0.0) return -1.0;
  return y;
}

// ===================================================================

// Stand alone Neural Network
// This version is with a single hidden layer (Non-deep Neural Network)
class NeuralNetwork {

  std::function<double(double)> Activatorfunc;
  std::function<double(double)> dActivatorfunc;
  std::function<double(double)> invActivatorfunc;

public:
  int activatorId;
  size_t nbInputs;
  size_t nbHiddenNodes;
  size_t nbOutputs;
  double learnRate;

  cmatrix W_ih;
  cvector bias_h;
  cmatrix W_ho;
  cvector bias_o;

  NeuralNetwork(size_t nbInputs_, size_t nbHiddenNodes_, size_t nbOutputs_, double learnRate_,
                int activatorId_ = SIGMOID_ACTIVATOR_ID) {
    nbInputs      = nbInputs_;
    nbHiddenNodes = nbHiddenNodes_;
    nbOutputs     = nbOutputs_;
    learnRate     = learnRate_;
    setActivator(activatorId_);

    // initiation (TODO: add possibility to use different strategies)
    W_ih = cmatrix_create(nbHiddenNodes, nbInputs);
    cmatrix_randomize_gaussian_inplace(W_ih, pow((double)nbHiddenNodes, -0.5));
    W_ho = cmatrix_create(nbOutputs, nbHiddenNodes);
    cmatrix_randomize_gaussian_inplace(W_ho, pow((double)nbOutputs, -0.5));
    bias_h.resize(nbHiddenNodes);
    cvector_randomize_gaussian_inplace(bias_h, pow((double)nbHiddenNodes, -0.5));
    bias_o.resize(nbOutputs);
    cvector_randomize_gaussian_inplace(bias_o, pow((double)nbOutputs, -0.5));
  }

  // construct the Neural Network by reading a text file
  NeuralNetwork(const char *filename) {
    std::ifstream file(filename);
    int actId = 0;
    file >> actId;
    setActivator(actId);
    file >> nbInputs >> nbHiddenNodes >> nbOutputs >> learnRate;

    W_ih = cmatrix_create(nbHiddenNodes, nbInputs);
    for (size_t i = 0; i < nbHiddenNodes; i++) {
      for (size_t j = 0; j < nbInputs; j++) { file >> W_ih[i][j]; }
    }

    bias_h.resize(nbHiddenNodes);
    for (size_t i = 0; i < bias_h.size(); i++) { file >> bias_h[i]; }

    W_ho = cmatrix_create(nbOutputs, nbHiddenNodes);
    for (size_t i = 0; i < nbOutputs; i++) {
      for (size_t j = 0; j < nbHiddenNodes; j++) { file >> W_ho[i][j]; }
    }

    bias_o.resize(nbOutputs);
    for (size_t i = 0; i < bias_o.size(); i++) { file >> bias_o[i]; }
  }

  void setActivator(int activatorId_ = 0) {
    activatorId = activatorId_;
    if (activatorId == SIGMOID_ACTIVATOR_ID) {
      Activatorfunc    = sigmoid;
      dActivatorfunc   = dsigmoid;
      invActivatorfunc = invSigmoid;
    } else if (activatorId == BINARY_STEP_ACTIVATOR_ID) {
      Activatorfunc    = binary;
      dActivatorfunc   = dBinary;
      invActivatorfunc = invBinary;
    } else if (activatorId == RELU_ACTIVATOR_ID) {
      Activatorfunc    = reLU;
      dActivatorfunc   = dReLU;
      invActivatorfunc = invReLU;
    } else if (activatorId == TANH_ACTIVATOR_ID) {
      Activatorfunc    = Tanh;
      dActivatorfunc   = dTanh;
      invActivatorfunc = invTanh;
    } else if (activatorId == IDENTITY_ACTIVATOR_ID) {
      Activatorfunc    = identity;
      dActivatorfunc   = dIdentity;
      invActivatorfunc = invIdentity;
    } else { // default activator is sigmoid
      Activatorfunc    = sigmoid;
      dActivatorfunc   = dsigmoid;
      invActivatorfunc = invSigmoid;
    }
  }

  void query(cvector &inputs, cvector &outputs) {
    if (inputs.size() != nbInputs) { std::cerr << "@NeuralNetwork::predict, size of input vector doesn't match\n"; }
    if (outputs.size() != nbOutputs) { std::cerr << "@NeuralNetwork::predict, size of output vector doesn't match\n"; }

    cvector H = cmatrix_mult(W_ih, inputs);
    cvector_add_inplace(H, bias_h);
    cvector_map_inplace(H, Activatorfunc);

    outputs = cmatrix_mult(W_ho, H);
    cvector_add_inplace(outputs, bias_o);
    cvector_map_inplace(outputs, Activatorfunc);
  }

  void train(cvector &inputs, cvector &targets) {
    if (inputs.size() != nbInputs) { std::cerr << "@NeuralNetwork::train, size of input vector doesn't match\n"; }
    if (targets.size() != nbOutputs) { std::cerr << "@NeuralNetwork::train, size of target vector doesn't match\n"; }

    // feedforward
    cvector H = cmatrix_mult(W_ih, inputs);
    cvector_add_inplace(H, bias_h);
    cvector_map_inplace(H, Activatorfunc);

    cvector outputs = cmatrix_mult(W_ho, H);
    cvector_add_inplace(outputs, bias_o);
    cvector_map_inplace(outputs, Activatorfunc);

    // errors in the output layer
    cvector output_errors = cvector_substract(targets, outputs);

    // gradient_ho
    cvector grad_ho(nbOutputs);
    for (size_t i = 0; i < nbOutputs; i++) { grad_ho[i] = learnRate * output_errors[i] * dActivatorfunc(outputs[i]); }
    // delta_ho
    cmatrix delta_ho = cmatrix_create(nbOutputs, nbHiddenNodes);
    for (size_t i = 0; i < nbOutputs; i++) {
      for (size_t j = 0; j < nbHiddenNodes; j++) { delta_ho[i][j] = grad_ho[i] * H[j]; }
    }
    // corrections
    cmatrix_add_inplace(W_ho, delta_ho);
    for (size_t i = 0; i < bias_o.size(); i++) { bias_o[i] += grad_ho[i]; }

    // errors in the hidden layer: W_ho^T * output_errors
    cvector hidden_errors(nbHiddenNodes);
    for (size_t i = 0; i < nbHiddenNodes; i++) {
      double sum = 0.0;
      for (size_t j = 0; j < nbOutputs; j++) { sum += W_ho[j][i] * output_errors[j]; }
      hidden_errors[i] = sum;
    }

    // gradient_ih
    cvector grad_ih(nbHiddenNodes);
    for (size_t i = 0; i < nbHiddenNodes; i++) { grad_ih[i] = learnRate * hidden_errors[i] * dActivatorfunc(H[i]); }
    // delta_ih
    cmatrix delta_ih = cmatrix_create(nbHiddenNodes, nbInputs);
    for (size_t i = 0; i < nbHiddenNodes; i++) {
      for (size_t j = 0; j < nbInputs; j++) { delta_ih[i][j] = grad_ih[i] * inputs[j]; }
    }
    // corrections
    cmatrix_add_inplace(W_ih, delta_ih);
    for (size_t i = 0; i < bias_h.size(); i++) { bias_h[i] += grad_ih[i]; }
  }

  void backQuery(cvector &ouputs, cvector &inputs) {
    // todo
  }

  // save the current Neural Network (may be already trained) into a text file
  void save(const char *filename) {
    std::ofstream file(filename);
    file << activatorId << '\n';
    file << nbInputs << ' ' << nbHiddenNodes << ' ' << nbOutputs << ' ' << learnRate << '\n';

    for (size_t i = 0; i < nbHiddenNodes; i++) {
      for (size_t j = 0; j < nbInputs; j++) { file << W_ih[i][j] << ' '; }
      file << '\n';
    }
    file << '\n';

    for (size_t i = 0; i < bias_h.size(); i++) { file << bias_h[i] << '\n'; }
    file << '\n';

    for (size_t i = 0; i < nbOutputs; i++) {
      for (size_t j = 0; j < nbHiddenNodes; j++) { file << W_ho[i][j] << ' '; }
      file << '\n';
    }
    file << '\n';

    for (size_t i = 0; i < bias_o.size(); i++) { file << bias_o[i] << '\n'; }
    file << '\n';
  }
};

#endif /* end of include guard: NN_HPP */
