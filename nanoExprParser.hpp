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

template <typename T> class nanoExprParser {
public:
  /**
   * Constructs a nanoExprParser object and initializes predefined constants.
   *
   * This constructor reserves space in the constants vector to prevent
   * reallocation, which could invalidate pointers stored in the variables map.
   * It also initializes the parser with three mathematical constants:
   * - pi: 3.14159265358979323846
   * - e: 2.71828182845904523536
   * - phi: 1.61803398874989484820
   */
  nanoExprParser() {
    constants.reserve(5); // so that the address T* in the map 'variables' will not change
    addConstant("pi", 3.14159265358979323846);
    addConstant("e", 2.71828182845904523536);
    addConstant("phi", 1.61803398874989484820);
  }

  /**
   * Parses a mathematical expression and stores the result in a variable.
   *
   * The expression is given as a string. The result is stored in the variable
   * passed by reference. The function returns true if the expression has been
   * parsed successfully, and false otherwise.
   */
  bool parse(const std::string &expr, T &result) {
    this->expr = expr;
    skip_whitespace();
    pos = 0;
    result = parse_expression();
    return pos == (expr.size() - 1);
  }

  /**
   * Adds a variable to the parser.
   *
   * The variable is identified by its name, and its value is passed by pointer.
   * The function stores the address of the variable in the variables map.
   */
  void addVariable(const std::string &name, T *value) { 
    variables[name] = value; }

  /**
   * Adds a constant to the parser.
   *
   * The constant is identified by its name, and its value is passed by value.
   * The function checks if the constant is already defined in the variables
   * map. If it is, the function does nothing. Otherwise, the constant is added
   * to the constants vector and the address of the last element of the vector
   * is stored in the variables map.
   */
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

  /**
   * Parses a mathematical expression and returns the result.
   *
   * This function parses a mathematical expression given in the string
   * 'expr'. The expression is expected to be a sequence of one or more
   * terms separated by '+' or '-' operators. The function returns the
   * result of the expression.
   *
   * The function sets the position 'pos' to the end of the expression if
   * the expression is not valid.
   */
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

  /**
   * Parses a term and returns the result.
   *
   * This function parses a term given in the string 'expr'. The term is
   * expected to be a sequence of one or more factors separated by '*' or '/'
   * operators. The function returns the result of the expression.
   *
   * The function sets the position 'pos' to the end of the expression if
   * the expression is not valid.
   */
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

  /**
   * Parses a factor and returns the result.
   *
   * This function parses a factor given in the string 'expr'. The factor is
   * expected to be a number, a variable, or a function of the form func(number)
   * where func is 'sqrt', 'sin', 'cos', 'tan', 'exp', 'log', 'log10', or 'pow'.
   * The function returns the result of the expression.
   *
   * The function sets the position 'pos' to the end of the expression if
   * the expression is not valid.
   */
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
          // std::cout << "'" << func << "' = " << *variables[func] << std::endl;
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

  /**
   * Parses a number from the expression and returns it.
   *
   * This can be either an integer or a floating point number, and can be
   * optionally followed by an exponent.
   *
   * @returns The parsed number. If the number is not valid, returns 0 and
   *   sets pos to the end of the expression to indicate an error.
   */
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

  /**
   * Retrieves the next character in the expression, skipping over any whitespace.
   *
   * @returns The next character in the expression, or '\0' if the end of the
   *   expression has been reached.
   */
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

  /**
   * Skips over any whitespace characters in the expression, incrementing the
   * parse position as it goes.
   */
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
