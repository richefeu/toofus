#pragma once

#include <cctype>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

template <typename T> class nanoExprParser {
public:
  nanoExprParser() {
    constants.reserve(5);
    addConstant("pi", 3.14159265358979323846);
    addConstant("e", 2.71828182845904523536);
    addConstant("phi", 1.61803398874989484820);
  }

  bool parse(const std::string &expr, T &result) {
    this->expr = expr;
    pos        = 0;
    has_error  = false;

    result = parse_expression();

    skip_whitespace();
    return !has_error && pos == this->expr.size();
  }

  void addVariable(const std::string &name, T *value) {
    variables[name] = value;
  }

  void addConstant(const std::string &name, T value) {
    if (variables.find(name) != variables.end()) return;
    constants.push_back(value);
    variables[name] = &constants.back();
  }

  void getValue(std::istream &is, T &value) {
    char C;
    is >> C;

    if (C == '$') {
      std::string expr_string;

      // read until closing $
      while (is.get(C) && C != '$') { expr_string += C; }

      // parse expression
      if (!parse(expr_string, value)) { value = std::nan(""); }

      if (std::isnan(value)) { std::cout << "Warning: Expression $" << expr_string << "$ resulted in NaN.\n"; }

    } else {
      is.putback(C);
      is >> value;
    }
  }

  void getValue(std::istream &is, int &value) {
    double dvalue;
    getValue(is, dvalue);
    value = static_cast<int>(std::round(dvalue));
  }

  void getValue(std::istream &is, size_t &value) {
    double dvalue;
    getValue(is, dvalue);
    value = static_cast<size_t>(std::round(fabs(dvalue)));
  }

private:
  std::string expr;
  size_t pos     = 0;
  bool has_error = false;

  std::map<std::string, T *> variables;
  std::vector<T> constants;

  /* ================= PARSER ================= */

  T parse_expression() {
    T result = parse_term();

    while (!has_error) {
      skip_whitespace();

      if (match('+')) {
        result += parse_term();
      } else if (match('-')) {
        result -= parse_term();
      } else {
        break;
      }
    }

    return result;
  }

  T parse_term() {
    T result = parse_factor();

    while (!has_error) {
      skip_whitespace();

      if (match('*')) {
        result *= parse_factor();
      } else if (match('/')) {
        T denom = parse_factor();
        if (denom == 0) return error();
        result /= denom;
      } else if (match('^')) {
        T exponent = parse_factor();
        result     = std::pow(result, exponent);
      } else {
        break;
      }
    }

    return result;
  }

  T parse_factor() {
    skip_whitespace();

    if (match('+')) return parse_factor();
    if (match('-')) return -parse_factor();

    if (match('(')) {
      T value = parse_expression();
      if (!match(')')) return error();
      return value;
    }

    if (std::isalpha(peek())) {
      std::string name = parse_identifier();

      if (match('(')) { return parse_function(name); }

      // variable / constant
      auto it = variables.find(name);
      if (it != variables.end()) { return *(it->second); }

      // unknown variable -> return 0 (as required by tests)
      return static_cast<T>(0);
    }

    return parse_number();
  }

  /* ================= FUNCTIONS ================= */

  T parse_function(const std::string &name) {
    if (name == "sqrt") {
      T v = parse_expression();
      if (!match(')')) return error();
      return std::sqrt(v);
    }

    if (name == "sin") {
      T v = parse_expression();
      if (!match(')')) return error();
      return std::sin(v);
    }

    if (name == "cos") {
      T v = parse_expression();
      if (!match(')')) return error();
      return std::cos(v);
    }

    if (name == "tan") {
      T v = parse_expression();
      if (!match(')')) return error();
      return std::tan(v);
    }

    if (name == "exp") {
      T v = parse_expression();
      if (!match(')')) return error();
      return std::exp(v);
    }

    if (name == "log") {
      T v = parse_expression();
      if (!match(')')) return error();
      return std::log(v);
    }

    if (name == "log10") {
      T v = parse_expression();
      if (!match(')')) return error();
      return std::log10(v);
    }

    if (name == "pow") {
      T base = parse_expression();
      if (!match(',')) return error();
      T exp = parse_expression();
      if (!match(')')) return error();
      return std::pow(base, exp);
    }

    return error(); // unknown function
  }

  /* ================= LEXER ================= */

  char peek() const {
    if (pos >= expr.size()) return '\0';
    return expr[pos];
  }

  bool match(char c) {
    skip_whitespace();
    if (peek() == c) {
      pos++;
      return true;
    }
    return false;
  }

  std::string parse_identifier() {
    std::string name;
    while (std::isalnum(peek())) { name += expr[pos++]; }
    return name;
  }

  T parse_number() {
    skip_whitespace();

    size_t start = pos;
    bool has_dot = false;

    if (peek() == '.') {
      has_dot = true;
      pos++;
    }

    while (std::isdigit(peek()) || (!has_dot && peek() == '.')) {
      if (peek() == '.') has_dot = true;
      pos++;
    }

    if (peek() == 'e' || peek() == 'E') {
      pos++;
      if (peek() == '+' || peek() == '-') pos++;
      while (std::isdigit(peek())) pos++;
    }

    if (start == pos) { return error(); }

    return static_cast<T>(std::stod(expr.substr(start, pos - start)));
  }

  void skip_whitespace() {
    while (pos < expr.size() && std::isspace(expr[pos])) pos++;
  }

  T error() {
    has_error = true;
    return 0;
  }
};
