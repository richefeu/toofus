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
// This is a stand-alone header file that implement a Deep Neural Network

#ifndef DEEPNEURALNETWORK_HPP
#define DEEPNEURALNETWORK_HPP

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <vector>

#define HIDDEN_LAYERS 0
#define OUTPUT_LAYER 1

struct Connection {
  double weight;
  double deltaWeight;
  Connection() : weight(0.0), deltaWeight(0.0) {}
};

class Neuron;
typedef std::vector<Neuron> Layer;

// ***************************
// ********* NEURON **********
// ***************************
class Neuron {
public:
  Neuron(unsigned numOutputs, unsigned myIndex) {
    for (unsigned c = 0; c < numOutputs; ++c) {
      m_outputWeights.push_back(Connection());
      m_outputWeights.back().weight = randomWeight();
    }

    m_myIndex = myIndex;
  }

  void setOutputVal(double val) { m_outputVal = val; }
  double getOutputVal() const { return m_outputVal; }

  void feedForward(const Layer &prevLayer, bool isLastLayer = false) {
    double sum = 0.0;

    for (unsigned n = 0; n < prevLayer.size(); ++n) {
      sum += prevLayer[n].getOutputVal() * prevLayer[n].m_outputWeights[m_myIndex].weight;
    }

    if (isLastLayer == true) {
      m_outputVal = Neuron::outputLayerActivationFunction(sum);
    } else {
      m_outputVal = Neuron::hiddenLayerActivationFunction(sum);
    }
  }

  void calcOutputGradients(double targetVal) {
    double delta = targetVal - m_outputVal;
    m_gradient = delta * Neuron::outputLayerActivationFunctionDerivative(m_outputVal);
  }

  void calcHiddenGradients(const Layer &nextLayer) {
    double dow = sumDOW(nextLayer);
    m_gradient = dow * Neuron::hiddenLayerActivationFunctionDerivative(m_outputVal);
  }

  void updateInputWeights(Layer &prevLayer) {
    for (unsigned n = 0; n < prevLayer.size(); ++n) {
      Neuron &neuron = prevLayer[n];
      double oldDeltaWeight = neuron.m_outputWeights[m_myIndex].deltaWeight;

      double newDeltaWeight = eta * neuron.getOutputVal() * m_gradient + alpha * oldDeltaWeight;

      neuron.m_outputWeights[m_myIndex].deltaWeight = newDeltaWeight;
      neuron.m_outputWeights[m_myIndex].weight += newDeltaWeight;
    }
  }

  // overall net learning rate, [0.0..1.0]
  static void setLearningRate(double Eta) {
    if (Eta >= 0.0 && Eta <= 1.0) {
      Neuron::eta = Eta;
    } else {
      std::cerr << "Learning rate must be in the range [0.0 .. 1.0]\n";
    }
  }

  // momentum, multiplier of last deltaWeight, [0.0..1.0]
  static void setMomentum(double Mom) {
    if (Mom >= 0.0 && Mom <= 1.0) {
      Neuron::alpha = Mom;
    } else {
      std::cerr << "Momentum must be in the range [0.0 .. 1.0]\n";
    }
  }

  /*
   * Example:

    Neuron::setHiddenLayerActivation(
      [](double x) -> double {return tanh(x);},
      [](double x) -> double {return 1.0 - x*x;}
    );

  */
  static void setHiddenLayerActivation(std::function<double(double)> func,
                                       std::function<double(double)> funcDerivative) {
    hiddenLayerActivationFunction = func;
    hiddenLayerActivationFunctionDerivative = funcDerivative;
  }

  static void setOutputLayerActivation(std::function<double(double)> func,
                                       std::function<double(double)> funcDerivative) {
    outputLayerActivationFunction = func;
    outputLayerActivationFunctionDerivative = funcDerivative;
  }

  /*
   * Example:
     Neuron::setActivation("ReLU", HIDDEN_LAYERS);
  */
  static void setActivation(std::string activationName, int hiddenOrOutput = HIDDEN_LAYERS,
                            Layer *outputLayer = nullptr) {

    std::function<double(double)> Function;
    std::function<double(double)> FunctionDerivative;
    double *param = nullptr;
    if (hiddenOrOutput == HIDDEN_LAYERS) {
      param = &(Neuron::hiddenLayerActivationFunctionParameter);
    } else if (hiddenOrOutput == OUTPUT_LAYER) {
      param = &(Neuron::outputLayerActivationFunctionParameter);
    }

    if (activationName == "ReLU") {
      Function = [](double x) -> double { return std::max(0.0, x); };
      FunctionDerivative = [](double x) -> double {
        if (x < 0.0)
          return 0.0;
        return 1.0;
      };
    } else if (activationName == "LReLU") {
      Function = [](double x) -> double { return std::max(0.1 * x, x); };
      FunctionDerivative = [](double x) -> double {
        if (x < 0.0)
          return 0.1;
        return 1.0;
      };
    } else if (activationName == "ELU") {
      assert(param != nullptr);
      Function = [param](double x) -> double {
        if (x < 0.0)
          return *param * (exp(x) - 1.0);
        return x;
      };
      FunctionDerivative = [param](double x) -> double {
        if (x < 0.0)
          return *param * exp(x);
        return 1.0;
      };
    } else if (activationName == "softmax") {
      assert(outputLayer != nullptr);
      Function = [outputLayer](double x) -> double {
        double sum = 0.0;
        for (unsigned i = 0; i < outputLayer->size() - 1; ++i) {
          sum += exp((*outputLayer)[i].getOutputVal());
        }
        return x / sum;
      };
      FunctionDerivative = [outputLayer](double x) -> double {
        double sum = 0.0;
        for (unsigned i = 0; i < outputLayer->size() - 1; ++i) {
          sum += exp((*outputLayer)[i].getOutputVal());
        }
        return 1.0 / sum;
      };
    } else if (activationName == "tanh") {
      Function = [](double x) -> double { return tanh(x); };
      FunctionDerivative = [](double x) -> double { return 1.0 - x * x; };
    } else if (activationName == "sigmoid") {
      Function = [](double x) -> double { return 1.0 / (1.0 + exp(-x)); };
      FunctionDerivative = [](double x) -> double {
        double sig = 1.0 / (1.0 + exp(-x));
        return sig * (1.0 - sig);
      };
    } else if (activationName == "none") {
      Function = [](double x) -> double { return x; };
      FunctionDerivative = [](double x) -> double { return 1.0; };
    }

    if (hiddenOrOutput == HIDDEN_LAYERS) {
      Neuron::hiddenLayerActivationFunction = Function;
      Neuron::hiddenLayerActivationFunctionDerivative = FunctionDerivative;
    } else if (hiddenOrOutput == OUTPUT_LAYER) {
      Neuron::outputLayerActivationFunction = Function;
      Neuron::outputLayerActivationFunctionDerivative = FunctionDerivative;
    }
  }

  static void setWeightinitialization(double range, double shift) {
    Neuron::weightInitMultiplier = range;
    Neuron::weightInitShift = shift;
  }

  static void setWeightinitialization(std::string model) {
    if (model == "unif[0..1]") {
      Neuron::weightInitMultiplier = 1.0;
      Neuron::weightInitShift = 0.0;
    } else if (model == "unif[-1..1]") {
      Neuron::weightInitMultiplier = 2.0;
      Neuron::weightInitShift = 0.5;
    }
  }

  static void setHiddenLayerActivationFunctionParameter(double value) {
    Neuron::hiddenLayerActivationFunctionParameter = value;
  }

  static void setOutputLayerActivationFunctionParameter(double value) {
    Neuron::outputLayerActivationFunctionParameter = value;
  }

private:
  static double eta;
  static double alpha;

  static double hiddenLayerActivationFunctionParameter;
  static double outputLayerActivationFunctionParameter;
  static std::function<double(double)> hiddenLayerActivationFunction;
  static std::function<double(double)> hiddenLayerActivationFunctionDerivative;
  static std::function<double(double)> outputLayerActivationFunction;
  static std::function<double(double)> outputLayerActivationFunctionDerivative;

  static double weightInitMultiplier;
  static double weightInitShift;
  static double randomWeight(void) { return weightInitMultiplier * (rand() / double(RAND_MAX) - weightInitShift); }

  double sumDOW(const Layer &nextLayer) const {
    double sum = 0.0;

    for (unsigned n = 0; n < nextLayer.size() - 1; ++n) {
      sum += m_outputWeights[n].weight * nextLayer[n].m_gradient;
    }

    return sum;
  }

  double m_outputVal;
  std::vector<Connection> m_outputWeights;
  unsigned m_myIndex;
  double m_gradient;
};

// default activation functions (and derivatives)
// hidden layers = tanh
// output layer = none (linear)
std::function<double(double)> Neuron::hiddenLayerActivationFunction = [](double x) -> double { return tanh(x); };
std::function<double(double)> Neuron::hiddenLayerActivationFunctionDerivative = [](double x) -> double {
  return 1.0 - x * x;
};
std::function<double(double)> Neuron::outputLayerActivationFunction = [](double x) -> double { return x; };
std::function<double(double)> Neuron::outputLayerActivationFunctionDerivative = [](double x) -> double { return 1.0; };

double Neuron::eta = 0.1;
double Neuron::alpha = 0.5;

double Neuron::hiddenLayerActivationFunctionParameter = 0.0;
double Neuron::outputLayerActivationFunctionParameter = 0.0;

double Neuron::weightInitMultiplier = 1.0;
double Neuron::weightInitShift = 0.0;

// ***********************************************************
// ********* DANN = DEEP ARTIFICIAL NEURAL NETWORK ***********
// ***********************************************************
class DANN {
public:
  DANN(const std::vector<unsigned> &topology) {
    unsigned numLayers = topology.size();
    for (unsigned layerNum = 0; layerNum < numLayers; ++layerNum) {
      m_layers.push_back(Layer());
      unsigned numOutputs = (layerNum == topology.size() - 1) ? 0 : topology[layerNum + 1];

      for (unsigned neuronNum = 0; neuronNum <= topology[layerNum]; ++neuronNum) {
        m_layers.back().push_back(Neuron(numOutputs, neuronNum));
      }

      m_layers.back().back().setOutputVal(DANN::m_biasValue);
    }
  }

  void feedForward(const std::vector<double> &inputVals) {
    assert(inputVals.size() == m_layers[0].size() - 1);

    for (unsigned i = 0; i < inputVals.size(); ++i) {
      m_layers[0][i].setOutputVal(inputVals[i]);
    }

    for (unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum) {
      Layer &prevLayer = m_layers[layerNum - 1];
      bool isLastLayer = (layerNum == m_layers.size() - 1);
      for (unsigned n = 0; n < m_layers[layerNum].size() - 1; ++n) {
        m_layers[layerNum][n].feedForward(prevLayer, isLastLayer);
      }
    }
  }

  void backProp(const std::vector<double> &targetVals) {
    Layer &outputLayer = m_layers.back();
    m_error = 0.0;

    for (unsigned n = 0; n < outputLayer.size() - 1; ++n) {
      double delta = targetVals[n] - outputLayer[n].getOutputVal();
      m_error += delta * delta;
    }
    m_error /= (double)(outputLayer.size() - 1); // get average error squared
    m_error = sqrt(m_error);                     // RMS

    m_recentAverageError =
        (m_recentAverageError * m_recentAverageSmoothingFactor + m_error) / (m_recentAverageSmoothingFactor + 1.0);

    for (unsigned n = 0; n < outputLayer.size() - 1; ++n) {
      outputLayer[n].calcOutputGradients(targetVals[n]);
    }

    for (unsigned layerNum = m_layers.size() - 2; layerNum > 0; --layerNum) {
      Layer &hiddenLayer = m_layers[layerNum];
      Layer &nextLayer = m_layers[layerNum + 1];

      for (unsigned n = 0; n < hiddenLayer.size(); ++n) {
        hiddenLayer[n].calcHiddenGradients(nextLayer);
      }
    }

    for (unsigned layerNum = m_layers.size() - 1; layerNum > 0; --layerNum) {
      Layer &layer = m_layers[layerNum];
      Layer &prevLayer = m_layers[layerNum - 1];

      for (unsigned n = 0; n < layer.size() - 1; ++n) {
        layer[n].updateInputWeights(prevLayer);
      }
    }
  }

  void getResults(std::vector<double> &resultVals) const {
    resultVals.clear();

    for (unsigned n = 0; n < m_layers.back().size() - 1; ++n) {
      resultVals.push_back(m_layers.back()[n].getOutputVal());
    }
  }

  double getRecentAverageError() const { return m_recentAverageError; }
  double getCurrentError() const { return m_error; }

  static void setAverageSmoothingFactor(double factor) { DANN::m_recentAverageSmoothingFactor = factor; }
  static void setBias(double value) { DANN::m_biasValue = value; }

private:
  std::vector<Layer> m_layers; // m_layers[layerNum][neuronNum]
  double m_error;
  double m_recentAverageError;
  static double m_recentAverageSmoothingFactor;
  static double m_biasValue;
};

double DANN::m_recentAverageSmoothingFactor = 100.0;
double DANN::m_biasValue = 1.0;

// ****************************************************

/*
class DANNData {
public:

  DANNData() {

  }


};
*/

#endif /* end of include guard: DEEPNEURALNETWORK_HPP */