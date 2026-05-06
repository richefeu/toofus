#include "../nanoExprParser2.hpp"
#include <iomanip>
#include <iostream>

int main() {
  nanoExprParser<double> parser;

  double result;

  std::cout << std::fixed << std::setprecision(6);

  // ================= BASIC ARITHMETIC =================
  std::cout << "=== Basic Arithmetic ===\n";

  parser.parse("1+2", result);
  std::cout << "1 + 2 = " << result << "\n";

  parser.parse("10 - 3 * 2", result);
  std::cout << "10 - 3 * 2 = " << result << "\n";

  parser.parse("(10 - 3) * 2", result);
  std::cout << "(10 - 3) * 2 = " << result << "\n\n";

  // ================= POW / EXPONENTS =================
  std::cout << "=== Powers ===\n";

  parser.parse("2^3", result);
  std::cout << "2^3 = " << result << "\n";

  parser.parse("pow(2,3)", result);
  std::cout << "pow(2,3) = " << result << "\n\n";

  // ================= CONSTANTS =================
  std::cout << "=== Constants ===\n";

  parser.parse("pi", result);
  std::cout << "pi = " << result << "\n";

  parser.parse("2*pi", result);
  std::cout << "2*pi = " << result << "\n";

  parser.parse("e", result);
  std::cout << "e = " << result << "\n\n";

  // ================= FUNCTIONS =================
  std::cout << "=== Functions ===\n";

  parser.parse("sqrt(9)", result);
  std::cout << "sqrt(9) = " << result << "\n";

  parser.parse("sin(0)", result);
  std::cout << "sin(0) = " << result << "\n";

  parser.parse("cos(0)", result);
  std::cout << "cos(0) = " << result << "\n";

  parser.parse("log(e)", result);
  std::cout << "log(e) = " << result << "\n";

  parser.parse("log10(100)", result);
  std::cout << "log10(100) = " << result << "\n\n";

  // ================= VARIABLES =================
  std::cout << "=== Variables ===\n";

  double x = 5;
  double y = 3;

  parser.addVariable("x", &x);
  parser.addVariable("y", &y);

  parser.parse("x + y", result);
  std::cout << "x + y = " << result << "\n";

  parser.parse("x * y + 2", result);
  std::cout << "x * y + 2 = " << result << "\n\n";

  // ================= UNKNOWN VARIABLE BEHAVIOR =================
  std::cout << "=== Unknown Variable (returns 0) ===\n";

  parser.parse("unknown + 5", result);
  std::cout << "unknown + 5 = " << result << "\n\n";

  // ================= COMPLEX EXPRESSION =================
  std::cout << "=== Complex Expression ===\n";

  parser.parse("sqrt(16) + pow(2,3) * (3+1)", result);
  std::cout << "sqrt(16) + pow(2,3) * (3+1) = " << result << "\n\n";

  // ================= ERROR CASES =================
  std::cout << "=== Invalid Expressions ===\n";

  bool ok;

  ok = parser.parse("5/0", result);
  std::cout << "5/0 valid? " << (ok ? "yes" : "no") << "\n";

  ok = parser.parse("2+", result);
  std::cout << "2+ valid? " << (ok ? "yes" : "no") << "\n";

  ok = parser.parse("(2+3", result);
  std::cout << "(2+3 valid? " << (ok ? "yes" : "no") << "\n";

  // ================= integration in conf-files =================
  std::stringstream ss1;
  ss1 << "$2-3.1*4$";
  std::stringstream ss2;
  ss2 << "$2-3.1*4$";
  std::stringstream ss3;
  ss3 << "$2-3.1*4$";

  double v_double;
  int v_int;
  size_t v_st;
  parser.getValue(ss1, v_double);
  std::cout << "v_double = " << v_double << "\n";
  parser.getValue(ss2, v_int);
  std::cout << "v_int = " << v_int << "\n";
  parser.getValue(ss3, v_st);
  std::cout << "v_st = " << v_st << "\n";

  return 0;
}