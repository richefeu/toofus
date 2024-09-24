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

#ifndef NANO_EXPR_PARSE_HPP
#define NANO_EXPR_PARSE_HPP

#include <cctype>
#include <cmath>
#include <map>
#include <string>
#include <vector>

// #include <iostream>

template <typename T> class nanoExprParser {
public:
  nanoExprParser() {
    constants.reserve(5); // so that the address T* in the map 'variables' will not change
    addConstant("pi", 3.14159265358979323846);
    addConstant("e", 2.71828182845904523536);
    addConstant("phi", 1.61803398874989484820);
  }

  bool parse(const std::string &expr, T &result) {
    this->expr = expr;
    skip_whitespace();
    pos = 0;
    result = parse_expression();
    return pos == (expr.size() - 1);
  }

  void addVariable(const std::string &name, T *value) {
    variables[name] = value;
  }

  void addConstant(const std::string &name, T value) {
    if (variables.find(name) != variables.end()) {
      return;
    }
    constants.push_back(value);
    variables[name] = &(constants[constants.size() - 1]);
  }

private:
  std::string expr;
  size_t pos;
  std::map<std::string, T *> variables;
  std::vector<T> constants;

  T parse_expression() {
    T result = parse_term();
    while (true) {
      char c = get_next_char();
      if (c == '+') {
        result += parse_term();
      } else if (c == '-') {
        result -= parse_term();
      } else {
        unget_char();
        break;
      }
    }
    return result;
  }

  T parse_term() {
    T result = parse_factor();
    while (true) {
      char c = get_next_char();
      if (c == '^') {
        T exponent = parse_factor();
        result = std::pow(result, exponent);
      } else if (c == '*') {
        result *= parse_factor();
      } else if (c == '/') {
        T denominator = parse_factor();
        if (denominator == 0) {
          pos = expr.size(); // set position to end of expression to indicate error
          return 0;
        }
        result /= denominator;
      } else {
        unget_char();
        break;
      }
    }
    return result;
  }

  T parse_factor() {
    char c = get_next_char();
    if (c == '(') {
      T result = parse_expression();
      if (get_next_char() != ')') {
        pos = expr.size(); // set position to end of expression to indicate error
        return 0;
      }
      return result;
    } else if (c == '-') {
      return -parse_factor();
    } else if (c == '+') {
      return parse_factor();
    } else if (isalpha(c)) {
      std::string func;
      func += c;
      while (isalpha(c = get_next_char())) {
        func += c;
      }
      if (c == '(') {
        if (func == "sqrt") {
          return std::sqrt(parse_factor());
        } else if (func == "sin") {
          return std::sin(parse_factor());
        } else if (func == "cos") {
          return std::cos(parse_factor());
        } else if (func == "tan") {
          return std::tan(parse_factor());
        } else if (func == "exp") {
          return std::exp(parse_factor());
        } else if (func == "log") {
          return std::log(parse_factor());
        } else if (func == "log10") {
          return std::log10(parse_factor());
        } else if (func == "pow") {
          T base = parse_factor();
          if (get_next_char() != ',') {
            pos = expr.size(); // set position to end of expression to indicate error
            return 0;
          }
          T exponent = parse_factor();
          if (get_next_char() != ')') {
            pos = expr.size(); // set position to end of expression to indicate error
            return 0;
          }
          return std::pow(base, exponent);
        }
      } else {
        unget_char();
        // return parse_number();
        if (variables.find(func) != variables.end()) {
          //std::cout << "'" << func << "' = " << *variables[func] << std::endl;
          return *variables[func];
        } else {
          return parse_number();
        }
      }
    } else {
      unget_char();
      return parse_number();
    }

    // Add this line to handle unknown function names
    pos = expr.size(); // set position to end of expression to indicate error
    return 0;
  }

  T parse_number() {
    T result = 0;
    bool negative = false;
    char c = get_next_char();
    if (c == '-') {
      negative = true;
      c = get_next_char();
    }
    while (c >= '0' && c <= '9') {
      result = result * 10 + (c - '0');
      c = get_next_char();
    }
    if (c == '.') {
      T fraction = 0;
      T denominator = 1;
      c = get_next_char();
      while (c >= '0' && c <= '9') {
        fraction = fraction * 10 + (c - '0');
        denominator *= 10;
        c = get_next_char();
      }
      result += static_cast<T>(fraction) / static_cast<T>(denominator);
    }
    if (c == 'e' || c == 'E') {
      bool exp_negative = false;
      c = get_next_char();
      if (c == '-') {
        exp_negative = true;
        c = get_next_char();
      } else if (c == '+') {
        c = get_next_char();
      }
      int exponent = 0;
      while (c >= '0' && c <= '9') {
        exponent = exponent * 10 + (c - '0');
        c = get_next_char();
      }
      result *= std::pow(10, exp_negative ? -exponent : exponent);
    }
    unget_char();
    return negative ? -result : result;
  }

  char get_next_char() {
    skip_whitespace();
    if (pos >= expr.size()) {
      return '\0';
    }
    return expr[pos++];
  }

  void unget_char() {
    if (pos > 0) {
      pos--;
    }
  }

  void skip_whitespace() {
    while (pos < expr.size() && std::isspace(expr[pos])) {
      pos++;
    }
  }
};

#endif // NANO_EXPR_PARSE_HPP

#if 0

#include <iostream>

int main() {
  nanoExprParser<double> parser;
  double x = 3.0;
  parser.addVariable("x", &x);
  std::string expr = "2.1 + pow(x,2) + cos(pi)";
  //std::string expr = "pi";
  double result;
  if (parser.parse(expr, result)) {
    std::cout << expr << " = " << result << std::endl;
  } else {
    std::cout << "Error parsing expression: " << expr << std::endl;
    std::cout << expr << " = " << result << std::endl;
  }
  return 0;
}

#endif
