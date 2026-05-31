#include "../exprParser.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

static int g_passed = 0, g_failed = 0;

static bool near(double a, double b, double tol = 1e-9) {
  return std::fabs(a - b) <= tol;
}

static void section(const char *title) {
  std::cout << "\n" << title << "\n" << std::string(std::string(title).size(), '-') << "\n";
}

static void eval(const char *expr, double expected, double tol = 1e-9) {
  exprParser<double> p;
  double result = 0.0;
  bool parsed   = p.parse(expr, result);
  bool pass     = parsed && near(result, expected, tol);

  std::cout << (pass ? "  OK" : "FAIL") << "  " << std::left << std::setw(34) << expr << "  =  ";
  if (!parsed) {
    std::cout << "parse error\n";
    std::ostringstream oss;
    p.printError(oss);
    std::istringstream iss(oss.str());
    std::string line;
    while (std::getline(iss, line)) std::cout << "        " << line << "\n";
  } else {
    std::cout << result;
    if (!pass) std::cout << "   (expected " << expected << ")";
    std::cout << "\n";
  }
  pass ? g_passed++ : g_failed++;
}

static void eval_fail(const char *expr) {
  exprParser<double> p;
  double result = 0.0;
  bool parsed   = p.parse(expr, result);
  bool pass     = !parsed;
  std::string q = "\"" + std::string(expr) + "\"";

  std::cout << (pass ? "  OK" : "FAIL") << "  " << std::left << std::setw(20) << q;
  if (!pass) {
    std::cout << "  should have failed, got " << result << "\n";
  } else {
    std::cout << "\n";
    std::ostringstream oss;
    p.printError(oss);
    std::istringstream iss(oss.str());
    std::string line;
    while (std::getline(iss, line)) std::cout << "        " << line << "\n";
  }
  pass ? g_passed++ : g_failed++;
}

// Runs expr through both interpreter and compiler, checks both give expected.
static void eval_both(const char *expr, double expected, double tol = 1e-9) {
  exprParser<double> pi, pc;
  double ri = 0.0;
  bool ok_i = pi.parse(expr, ri);

  auto ce   = pc.compile(expr);
  bool ok_c = !!ce;
  double rc = ok_c ? ce.eval() : 0.0;

  bool pass = ok_i && ok_c && near(ri, expected, tol) && near(rc, expected, tol);

  std::cout << (pass ? "  OK" : "FAIL") << "  " << std::left << std::setw(34) << expr << "  =  ";
  if (pass) {
    std::cout << ri << "\n";
  } else {
    if (!ok_i) std::cout << "interp:err";
    else std::cout << "interp=" << ri;
    if (!ok_c) std::cout << "  compiled:err";
    else std::cout << "  compiled=" << rc;
    std::cout << "  (expected " << expected << ")\n";
  }
  pass ? g_passed++ : g_failed++;
}

static void check(const std::string &desc, bool cond) {
  std::cout << (cond ? "  OK" : "FAIL") << "  " << desc << "\n";
  cond ? g_passed++ : g_failed++;
}

/* ================================================================ */

int main() {
  const double PI = std::acos(-1.0);

  section("Literals");
  eval("42", 42.0);
  eval("3.14", 3.14);
  eval(".5", 0.5);
  eval("1.5e2", 150.0);
  eval("2.5e-1", 0.25);
  eval("1E3", 1000.0);

  section("Basic arithmetic");
  eval("1 + 2", 3.0);
  eval("5 - 3", 2.0);
  eval("2 * 3", 6.0);
  eval("10 / 4", 2.5);
  eval("10 % 3", 1.0);

  section("Operator precedence  (* > +, ^ > *)");
  eval("2 + 3 * 4", 14.0);
  eval("2 * 3^2", 18.0);
  eval("1 + 2^3", 9.0);

  section("Power right-associativity  (2^3^2 = 2^(3^2) = 512)");
  eval("2^3^2", 512.0);

  section("Unary operators");
  eval("-3", -3.0);
  eval("+3", 3.0);
  eval("--3", 3.0);
  eval("-(2+3)", -5.0);
  eval("2 * -3", -6.0);

  section("Parentheses");
  eval("(2+3)*4", 20.0);
  eval("((2+3))*2", 10.0);

  section("Built-in constants");
  eval("pi", PI, 1e-12);
  eval("e", 2.71828182845904523536, 1e-12);
  eval("phi", 1.61803398874989484820, 1e-12);
  eval("2*pi", 2.0 * PI, 1e-12);

  section("User variables and constants");
  {
    exprParser<double> p;
    double x = 4.0;
    p.addVariable("x", &x);
    double r;
    check("addVariable(x=4):  x*2 = 8", p.parse("x*2", r) && near(r, 8.0));
    x = 9.0;
    check("after x=9:         x   = 9", p.parse("x", r) && near(r, 9.0));
  }
  {
    exprParser<double> p;
    p.addConstant("K", 10.0);
    double x = 3.0;
    p.addVariable("x", &x);
    double r;
    check("addConstant(K=10): K+x = 13", p.parse("K+x", r) && near(r, 13.0));
  }
  {
    exprParser<double> p;
    double x = 5.0;
    p.addVariable("x", &x);
    p.addConstant("x", 99.0);
    double r;
    check("addConstant cannot shadow existing variable", p.parse("x", r) && near(r, 5.0));
  }

  section("removeVariable / clear");
  {
    exprParser<double> p;
    double x = 7.0;
    p.addVariable("x", &x);
    double r;
    check("before remove: x = 7", p.parse("x", r) && near(r, 7.0));
    p.removeVariable("x");
    check("after remove:  x → error", !p.parse("x", r));
  }
  {
    exprParser<double> p;
    double x = 3.0;
    p.addVariable("x", &x);
    p.addConstant("K", 5.0);
    p.clear();
    double r;
    check("after clear: x → error", !p.parse("x", r));
    check("after clear: K → error", !p.parse("K", r));
    check("after clear: pi restored", p.parse("pi", r) && near(r, PI, 1e-12));
  }

  section("Pointer stability: 20 user constants");
  {
    exprParser<double> p;
    for (int i = 0; i < 20; i++) p.addConstant("c" + std::to_string(i), double(i));
    double r;
    check("c0  = 0", p.parse("c0", r) && near(r, 0.0));
    check("c9  = 9", p.parse("c9", r) && near(r, 9.0));
    check("c19 = 19", p.parse("c19", r) && near(r, 19.0));
  }

  section("Identifier with underscore");
  {
    exprParser<double> p;
    double my_var = 7.0;
    p.addVariable("my_var", &my_var);
    double r;
    check("my_var + 1 = 8", p.parse("my_var+1", r) && near(r, 8.0));
  }

  section("Whitespace tolerance");
  eval("  2 + 3  ", 5.0);
  eval("2  *  3", 6.0);

  section("One-argument functions");
  eval("sqrt(4)", 2.0);
  eval("abs(-5)", 5.0);
  eval("sin(0)", 0.0);
  eval("cos(0)", 1.0);
  eval("tan(0)", 0.0);
  eval("asin(1)", PI / 2, 1e-12);
  eval("acos(1)", 0.0, 1e-12);
  eval("atan(1)", PI / 4, 1e-12);
  eval("exp(0)", 1.0);
  eval("log(e)", 1.0, 1e-12);
  eval("log10(100)", 2.0, 1e-12);
  eval("floor(2.9)", 2.0);
  eval("ceil(2.1)", 3.0);
  eval("round(2.5)", 3.0);
  eval("round(2.4)", 2.0);

  section("Two-argument functions");
  eval("pow(2,10)", 1024.0);
  eval("atan2(1,1)", PI / 4, 1e-12);
  eval("min(3,5)", 3.0);
  eval("max(3,5)", 5.0);
  eval("fmod(10,3)", 1.0);

  section("Composed expressions");
  eval("sin(pi/4)^2 + cos(pi/4)^2", 1.0, 1e-12);
  eval("sqrt(3^2 + 4^2)", 5.0, 1e-12);
  eval("max(abs(-3), sqrt(4))", 3.0);

  /* -------- COMPARISON OPERATORS -------- */

  section("Comparison operators  (return 0 or 1)");
  eval("3 > 2", 1.0);
  eval("2 > 3", 0.0);
  eval("2 < 3", 1.0);
  eval("3 < 2", 0.0);
  eval("2 <= 2", 1.0);
  eval("3 <= 2", 0.0);
  eval("2 >= 2", 1.0);
  eval("1 >= 2", 0.0);
  eval("2 == 2", 1.0);
  eval("2 == 3", 0.0);
  eval("2 != 3", 1.0);
  eval("2 != 2", 0.0);
  eval("(3 > 2) + (1 < 2)", 2.0); // comparisons usable in arithmetic

  section("Comparisons with variables");
  {
    exprParser<double> p;
    double x = 5.0;
    p.addVariable("x", &x);
    double r;
    check("x=5:  x > 0 = 1", p.parse("x > 0", r) && near(r, 1.0));
    x = -1.0;
    check("x=-1: x > 0 = 0", p.parse("x > 0", r) && near(r, 0.0));
  }

  section("Comparison inside parentheses");
  eval("(2 > 1) * 5", 5.0); // (true) * 5
  eval("(1 > 2) * 5", 0.0); // (false) * 5

  /* -------- if() FUNCTION -------- */

  section("if(cond, a, b)  — conditional expression");
  eval("if(1, 2, 3)", 2.0);
  eval("if(0, 2, 3)", 3.0);
  eval("if(-1, 2, 3)", 2.0);     // non-zero → true
  eval("if(1>2, 10, 20)", 20.0); // comparison as condition
  eval("if(2>1, 10, 20)", 10.0);

  section("if() with variables");
  {
    exprParser<double> p;
    double x = 5.0;
    p.addVariable("x", &x);
    double r;
    check("x=5:  if(x>0,x,-x) = 5", p.parse("if(x>0,x,-x)", r) && near(r, 5.0));
    x = -3.0;
    check("x=-3: if(x>0,x,-x) = 3", p.parse("if(x>0,x,-x)", r) && near(r, 3.0));
  }
  {
    exprParser<double> p;
    double x = 4.0;
    p.addVariable("x", &x);
    double r;
    check("x=4 >= 0:  if(x>=0,sqrt(x),0) = 2", p.parse("if(x>=0,sqrt(x),0)", r) && near(r, 2.0));
    x = -1.0;
    check("x=-1 < 0:  if(x>=0,sqrt(x),0) = 0", p.parse("if(x>=0,sqrt(x),0)", r) && near(r, 0.0));
  }

  /* -------- EXPECTED ERRORS -------- */

  section("Expected parse errors (with diagnostic)");
  eval_fail("foo");
  eval_fail("blah(1)");
  eval_fail("(2+3");
  eval_fail("1/0");
  eval_fail("5%0");
  eval_fail("");
  eval_fail("2+3 xyz");
  eval_fail("*3");
  eval_fail("pow(2 3)");
  eval_fail(".");
  eval_fail("if(1,2)"); // too few args

  section("Depth limit (max 200 nested levels)");
  {
    // 150 levels should succeed
    std::string e150(150, '(');
    e150 += "1";
    e150 += std::string(150, ')');
    exprParser<double> p;
    double r;
    check("150 levels of '(' → OK", p.parse(e150, r) && near(r, 1.0));
  }
  {
    // 201 levels should error
    std::string e201(201, '(');
    e201 += "1";
    e201 += std::string(201, ')');
    exprParser<double> p;
    double r;
    check("201 levels of '(' → error", !p.parse(e201, r));
  }

  section("printError / errorMsg");
  {
    exprParser<double> p;
    double r;
    p.parse("2 + foo", r);
    check("errorMsg: \"" + p.errorMsg() + "\"", p.errorMsg() == "unknown identifier 'foo'");
    check("errorPos = 4", p.errorPos() == 4);
  }
  {
    exprParser<double> p;
    double r;
    p.parse("2+3", r);
    check("errorMsg clear on success", p.errorMsg().empty());
    check("errorPos npos on success", p.errorPos() == std::string::npos);
  }

  /* -------- COMPILE: CORRECTNESS -------- */

  section("Compile: same results as interpreter");
  eval_both("42", 42.0);
  eval_both("2 + 3 * 4", 14.0);
  eval_both("2 * 3^2", 18.0);
  eval_both("2^3^2", 512.0);
  eval_both("-(2+3)", -5.0);
  eval_both("(2+3)*4", 20.0);
  eval_both("pi", PI, 1e-12);
  eval_both("sin(pi/4)^2 + cos(pi/4)^2", 1.0, 1e-12);
  eval_both("sqrt(3^2 + 4^2)", 5.0, 1e-12);
  eval_both("3 > 2", 1.0);
  eval_both("2 >= 3", 0.0);
  eval_both("if(1>0, 10, 20)", 10.0);
  eval_both("if(0,   10, 20)", 20.0);

  section("Compile: variables are live");
  {
    exprParser<double> p;
    double x = 3.0;
    p.addVariable("x", &x);
    auto ce = p.compile("x^2 + 1");
    check("x=3:  x^2+1 = 10", !!ce && near(ce.eval(), 10.0));
    x = 5.0;
    check("x=5:  x^2+1 = 26", !!ce && near(ce.eval(), 26.0));
    x = 0.0;
    check("x=0:  x^2+1 =  1", !!ce && near(ce.eval(), 1.0));
  }
  {
    exprParser<double> p;
    double x = 4.0;
    p.addVariable("x", &x);
    auto ce = p.compile("if(x>=0, sqrt(x), 0)");
    check("x=4:  if → sqrt(4)=2", !!ce && near(ce.eval(), 2.0));
    x = -1.0;
    check("x=-1: if → 0", !!ce && near(ce.eval(), 0.0));
  }

  section("Compile: constant folding");
  {
    exprParser<double> p;
    {
      auto ce = p.compile("2 * pi");
      check("\"2*pi\"        → 1 instruction (folded)", !!ce && ce.size() == 1);
      check("\"2*pi\"        = 6.28318…", !!ce && near(ce.eval(), 2.0 * PI, 1e-12));
    }
    {
      auto ce = p.compile("sin(pi/6)");
      check("\"sin(pi/6)\"   → 1 instruction (folded)", !!ce && ce.size() == 1);
      check("\"sin(pi/6)\"   = 0.5", !!ce && near(ce.eval(), 0.5, 1e-12));
    }
    {
      // with a variable: not fully foldable
      double x = 1.0;
      p.addVariable("x", &x);
      auto ce = p.compile("x + 2*pi");
      check("\"x + 2*pi\"    → 3 instructions (PUSH_V + PUSH_K(2π) + ADD)", !!ce && ce.size() == 3);
    }
    {
      auto ce = p.compile("if(1>0, sin(pi/6), cos(pi/3))");
      check("\"if(1>0,…)\"   → 1 instruction (all const, folded)", !!ce && ce.size() == 1);
      check("\"if(1>0,…)\"   = 0.5", !!ce && near(ce.eval(), 0.5, 1e-12));
    }
  }

  section("Compile: error handling");
  {
    exprParser<double> p;
    auto ce = p.compile("foo+1");
    check("compile(\"foo+1\") → !ce", !ce);
    check("errorMsg = unknown identifier 'foo'", p.errorMsg() == "unknown identifier 'foo'");
  }
  {
    exprParser<double> p;
    auto ce = p.compile("sin(pi/4");
    std::ostringstream oss;
    p.printError(oss);
    check("compile(\"sin(pi/4\") → !ce", !ce);
    check("printError contains '^'", oss.str().find('^') != std::string::npos);
  }

  section("getValue ($…$ and type overloads)");
  {
    exprParser<double> p;
    double v;
    std::istringstream ss("$2+3$");
    p.getValue(ss, v);
    check("\"$2+3$\"  → double 5", near(v, 5.0));
  }
  {
    exprParser<double> p;
    double v;
    std::istringstream ss("42.5");
    p.getValue(ss, v);
    check("\"42.5\"   → double 42.5", near(v, 42.5));
  }
  {
    exprParser<double> p;
    int v;
    std::istringstream ss("3.7");
    p.getValue(ss, v);
    check("\"3.7\"    → int 4 (round)", v == 4);
  }
  {
    exprParser<double> p;
    size_t v;
    std::istringstream ss("-5.2");
    p.getValue(ss, v);
    check("\"-5.2\"   → size_t 5 (fabs+round)", v == 5);
  }

  /* -------- PERFORMANCE -------- */

  section("Performance: interpreted vs compiled");
  {
    exprParser<double> p;
    double x = 0.0;
    p.addVariable("x", &x);
    const std::string expr = "if(x>=0, sqrt(x), -x) + sin(pi/4)^2";
    const int N            = 500'000;
    auto ce                = p.compile(expr);

    auto t0   = std::chrono::steady_clock::now();
    double s1 = 0.0;
    for (int i = 0; i < N; ++i) {
      x = (i - N / 2.0) * 1e-4;
      double r;
      p.parse(expr, r);
      s1 += r;
    }
    auto t1   = std::chrono::steady_clock::now();
    double s2 = 0.0;
    for (int i = 0; i < N; ++i) {
      x = (i - N / 2.0) * 1e-4;
      s2 += ce.eval();
    }
    auto t2 = std::chrono::steady_clock::now();

    using us  = std::chrono::duration<double, std::micro>;
    double ti = us(t1 - t0).count(), tc = us(t2 - t1).count();
    std::cout << "  " << N << " evals of  \"" << expr << "\"\n";
    std::cout << "  interpreted: " << std::fixed << std::setprecision(1) << ti / 1000 << " ms\n";
    std::cout << "  compiled:    " << tc / 1000 << " ms\n";
    std::cout << "  speedup:     " << std::setprecision(1) << ti / tc << "x\n";
    check("sums match (interp == compiled)", near(s1, s2, 1e-4));
    check("compiled is faster", ti / tc > 2.0);
  }

  std::cout << "\n" << g_passed << " passed, " << g_failed << " failed.\n";
  return g_failed > 0 ? 1 : 0;
}
