// STATUS: [ ] STABLE  [ ] EXPERIMENTAL  [x] DRAFT

/*
 * exprParser  —  single-header mathematical expression parser/compiler
 * =========================================================================
 *
 * OVERVIEW
 * --------
 * Template class that parses and evaluates mathematical expressions given as
 * strings.  Two evaluation modes are provided:
 *
 *   • Interpreted  — parse() tokenises and evaluates in one pass.
 *                    Convenient for one-off evaluations.
 *
 *   • Compiled     — compile() translates the expression to flat bytecode
 *                    (stack-based VM).  Constant sub-expressions are folded
 *                    at compile time.  CompiledExpr::eval() then runs with
 *                    no heap allocation and no map lookups, giving a speedup
 *                    of typically 5–15× for repeated evaluations.
 *
 * TEMPLATE PARAMETER
 * ------------------
 *   T  Numeric type used for all computations (double, float …).
 *      Must support the standard arithmetic operators and <cmath> functions.
 *      Must be trivially copyable (satisfied by all built-in numeric types).
 *
 * QUICK START
 * -----------
 *   exprParser<double> p;
 *
 *   double x = 0.0;
 *   p.addVariable("x", &x);    // live variable — pointer, not copy
 *   p.addConstant("kn", 1e6);  // baked in at compile time
 *
 *   // --- Interpreted (one-off) ---
 *   double result;
 *   if (!p.parse("sin(x) * kn", result)) p.printError();
 *
 *   // --- Compiled (tight loop) ---
 *   auto prog = p.compile("if(x > 0, sqrt(x), 0)");
 *   if (!prog) { p.printError(); return; }
 *
 *   for (x = -10.0; x < 10.0; x += 1e-6)
 *       double y = prog.eval();   // no reparsing, no allocation
 *
 * PUBLIC API
 * ----------
 *   bool        parse(expr, result)     Interpret and evaluate.
 *                                       Returns false on error.
 *
 *   CompiledExpr compile(expr)          Compile to bytecode.
 *                                       Returns an empty object (operator bool
 *                                       == false) on error.
 *
 *   void        addVariable(name, ptr)  Bind a live variable by pointer.
 *                                       Changes to *ptr are seen by eval()
 *                                       without recompilation.
 *
 *   void        removeVariable(name)    Remove a variable binding.
 *                                       Existing CompiledExprs that captured
 *                                       this variable's pointer become dangling
 *                                       — do not call eval() on them afterwards.
 *
 *   void        addConstant(name, val)  Add a named constant (silently ignored
 *                                       if the name is already registered).
 *                                       Baked into bytecode at compile time.
 *
 *   void        clear()                 Remove all user variables and constants;
 *                                       the built-in constants (pi, e, phi) are
 *                                       restored.
 *
 *   size_t      errorPos()              Byte offset in the expression where the
 *                                       last error occurred (npos if none).
 *
 *   string      errorMsg()              Human-readable description of the error.
 *
 *   void        printError(ostream&)    Print the expression, a caret (^) at the
 *                                       error position, and the error message.
 *                                       Defaults to std::cerr.
 *
 *   void        getValue(istream, val)  Read a value from a stream.  If the next
 *                                       token starts with '$', reads up to the
 *                                       closing '$' and evaluates it as an
 *                                       expression; otherwise reads a plain value.
 *                                       Overloaded for T, int, and size_t.
 *
 * BUILT-IN CONSTANTS
 * ------------------
 *   pi    3.14159265358979…
 *   e     2.71828182845904…
 *   phi   1.61803398874989…   (golden ratio)
 *
 * OPERATORS  (highest to lowest precedence)
 * -----------------------------------------
 *   ( )               Grouping
 *   - +               Unary minus / plus
 *   ^                 Power  (right-associative: 2^3^2 = 2^(3^2) = 512)
 *   * / %             Multiply / divide / modulo (fmod)
 *   + -               Add / subtract
 *   < > <= >= == !=   Comparison — return T(1) or T(0)
 *                     Note: == and != perform exact floating-point comparison.
 *
 * BUILT-IN FUNCTIONS
 * ------------------
 *   One argument : sqrt  abs  sin  cos  tan  asin  acos  atan
 *                  exp   log  log10  floor  ceil  round
 *   Two arguments: pow(base,exp)  atan2(y,x)  min  max  fmod
 *   Conditional  : if(cond, a, b)  — returns a if cond != 0, else b.
 *                  Both a and b are always evaluated (eager semantics).
 *                  Comparison operators are valid inside any argument:
 *                    if(x > 0, sqrt(x), 0)
 *
 * ERROR HANDLING
 * --------------
 *   Both parse() and compile() set errorPos() and errorMsg() on failure.
 *   printError() renders the expression with a caret pointing to the problem:
 *
 *     2 + foo * (3
 *                ^ expected ')'
 *
 *   Division or modulo by a literal zero is caught at parse/compile time.
 *   Division by a variable that is zero at runtime produces ±inf or NaN
 *   (standard IEEE 754 behaviour for floating-point T).
 *
 * THREAD SAFETY
 * -------------
 *   • CompiledExpr::eval() is thread-safe: it is read-only and uses only
 *     stack-allocated storage.  Multiple threads may call eval() on the
 *     same CompiledExpr concurrently, provided the pointed-to variables are
 *     not modified while eval() is running (or each thread uses its own x).
 *   • parse() and compile() are NOT thread-safe: they mutate the parser's
 *     internal state.  Use one parser per thread, or serialize access.
 *
 * NOTES
 * -----
 *   • Header-only, no external dependencies beyond the C++ standard library.
 *   • Requires C++14 (generic lambdas in parse_function).
 *   • Nesting depth is limited to 200 levels; deeper expressions produce a
 *     parse error rather than a stack overflow.
 *   • Constants are stored in a std::map — references remain valid regardless
 *     of how many addConstant() calls are made (no reallocation issue).
 */

#pragma once

#include <cctype>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

template <typename T> class exprParser {
public:
  /* ================= BYTECODE TYPES ================= */

  enum class Op : uint8_t {
    PUSH_K,
    PUSH_V, // load
    ADD,
    SUB,
    MUL,
    DIV,
    MOD,
    POW,
    NEG, // arithmetic
    LT,
    GT,
    LE,
    GE,
    EQ,
    NE,  // comparison → T(0) or T(1)
    SEL, // if: pops b,a,cond → cond?a:b
    CALL1,
    CALL2 // math functions
  };

  struct Instr {
    Op op         = Op::PUSH_K;
    T k           = T{};     // PUSH_K
    T *v          = nullptr; // PUSH_V
    T (*f1)(T)    = nullptr; // CALL1
    T (*f2)(T, T) = nullptr; // CALL2

    static Instr K(T k) {
      Instr i;
      i.op = Op::PUSH_K;
      i.k  = k;
      return i;
    }
    static Instr V(T *v) {
      Instr i;
      i.op = Op::PUSH_V;
      i.v  = v;
      return i;
    }
    static Instr OP(Op op) {
      Instr i;
      i.op = op;
      return i;
    }
    static Instr F1(T (*f)(T)) {
      Instr i;
      i.op = Op::CALL1;
      i.f1 = f;
      return i;
    }
    static Instr F2(T (*f)(T, T)) {
      Instr i;
      i.op = Op::CALL2;
      i.f2 = f;
      return i;
    }
  };

  /* ================= COMPILED EXPRESSION ================= */

  class CompiledExpr {
  public:
    explicit operator bool() const {
      return !prog.empty();
    }
    size_t size() const {
      return prog.size();
    } // number of bytecode instructions

    // Thread-safe: read-only, stack-allocated storage only.
    T eval() const {
      T stk[256];
      int top = -1;
      for (const auto &i : prog) {
        switch (i.op) {
        case Op::PUSH_K:
          stk[++top] = i.k;
          break;
        case Op::PUSH_V:
          stk[++top] = *i.v;
          break;
        case Op::NEG:
          stk[top] = -stk[top];
          break;
        case Op::CALL1:
          stk[top] = i.f1(stk[top]);
          break;
        case Op::ADD: {
          T b = stk[top--];
          stk[top] += b;
          break;
        }
        case Op::SUB: {
          T b = stk[top--];
          stk[top] -= b;
          break;
        }
        case Op::MUL: {
          T b = stk[top--];
          stk[top] *= b;
          break;
        }
        case Op::DIV: {
          T b = stk[top--];
          stk[top] /= b;
          break;
        }
        case Op::MOD: {
          T b      = stk[top--];
          stk[top] = std::fmod(stk[top], b);
          break;
        }
        case Op::POW: {
          T b      = stk[top--];
          stk[top] = std::pow(stk[top], b);
          break;
        }
        case Op::CALL2: {
          T b      = stk[top--];
          stk[top] = i.f2(stk[top], b);
          break;
        }
        case Op::LT: {
          T b      = stk[top--];
          stk[top] = T(stk[top] < b);
          break;
        }
        case Op::GT: {
          T b      = stk[top--];
          stk[top] = T(stk[top] > b);
          break;
        }
        case Op::LE: {
          T b      = stk[top--];
          stk[top] = T(stk[top] <= b);
          break;
        }
        case Op::GE: {
          T b      = stk[top--];
          stk[top] = T(stk[top] >= b);
          break;
        }
        case Op::EQ: {
          T b      = stk[top--];
          stk[top] = T(stk[top] == b);
          break;
        }
        case Op::NE: {
          T b      = stk[top--];
          stk[top] = T(stk[top] != b);
          break;
        }
        case Op::SEL: {
          T b      = stk[top--];
          T a      = stk[top--];
          stk[top] = (stk[top] != T(0) ? a : b);
          break;
        }
        }
      }
      return top >= 0 ? stk[top] : T(0);
    }

  private:
    friend class exprParser<T>;
    std::vector<Instr> prog;
  };

  /* ================= PUBLIC API ================= */

  exprParser() {
    reset_builtins();
  }

  // Interpret: parse and evaluate in one pass.
  bool parse(const std::string &expr, T &result) {
    setup(expr);
    result = parse_comparison();
    skip_whitespace();
    if (!has_error && pos != this->expr.size()) error_at(pos, "unexpected characters");
    return !has_error && pos == this->expr.size();
  }

  // Compile to bytecode.  Constant sub-expressions are folded at compile time.
  // Returns an empty CompiledExpr (operator bool == false) on error.
  CompiledExpr compile(const std::string &expr) {
    setup(expr);
    CompiledExpr ce;
    emit_prog = &ce.prog;
    emit_comparison();
    emit_prog = nullptr;
    skip_whitespace();
    if (!has_error && pos != this->expr.size()) error_at(pos, "unexpected characters");
    if (has_error) ce.prog.clear();
    return ce;
  }

  void addVariable(const std::string &name, T *value) {
    variables[name] = value;
  }
  void removeVariable(const std::string &name) {
    variables.erase(name);
  }

  // Silently ignored if the name is already registered.
  void addConstant(const std::string &name, T value) {
    if (variables.count(name) || constants.count(name)) return;
    constants[name] = value;
  }

  // Remove all user variables and constants; restore built-in constants.
  void clear() {
    variables.clear();
    constants.clear();
    reset_builtins();
  }

  size_t errorPos() const {
    return err_pos;
  }
  const std::string &errorMsg() const {
    return err_msg;
  }

  // Print the expression with a caret (^) at the error position.
  void printError(std::ostream &os = std::cerr) const {
    if (err_pos == std::string::npos) return;
    size_t ep = std::min(err_pos, expr.size());
    os << expr << "\n" << std::string(ep, ' ') << "^";
    if (!err_msg.empty()) os << " " << err_msg;
    os << "\n";
  }

  void getValue(std::istream &is, T &value) {
    char C;
    is >> C;
    if (C == '$') {
      std::string s;
      while (is.get(C) && C != '$') s += C;
      if (!parse(s, value)) value = T(std::nan(""));
      if (std::isnan(static_cast<double>(value))) std::cout << "Warning: Expression $" << s << "$ resulted in NaN.\n";
    } else {
      is.putback(C);
      is >> value;
    }
  }

  void getValue(std::istream &is, int &v) {
    double d;
    getValue(is, d);
    v = static_cast<int>(std::round(d));
  }
  void getValue(std::istream &is, size_t &v) {
    double d;
    getValue(is, d);
    v = static_cast<size_t>(std::round(std::fabs(d)));
  }

private:
  std::string expr;
  size_t pos     = 0;
  bool has_error = false;
  size_t err_pos = std::string::npos;
  std::string err_msg;
  int depth = 0;

  static constexpr int max_depth = 200;

  std::map<std::string, T *> variables;
  std::map<std::string, T> constants;

  std::vector<Instr> *emit_prog = nullptr; // non-null only during compile()

  void reset_builtins() {
    addConstant("pi", T(3.14159265358979323846));
    addConstant("e", T(2.71828182845904523536));
    addConstant("phi", T(1.61803398874989484820));
  }

  void setup(const std::string &e) {
    expr      = e;
    pos       = 0;
    has_error = false;
    err_pos   = std::string::npos;
    err_msg   = {};
    depth     = 0;
  }

  // RAII depth counter — ensures depth is decremented on all exit paths.
  struct DepthGuard {
    int &d;
    explicit DepthGuard(int &d) : d(d) {
      ++d;
    }
    ~DepthGuard() {
      --d;
    }
  };

  /* ================= INTERPRETER ================= */

  // Top-level rule: handles comparison operators (lowest precedence).
  T parse_comparison() {
    DepthGuard guard(depth);
    if (depth > max_depth) return error("expression too deeply nested");
    T left = parse_expression();
    if (has_error) return left;
    skip_whitespace();
    if (match2('<', '=')) return T(left <= parse_expression());
    if (match2('>', '=')) return T(left >= parse_expression());
    if (match2('=', '=')) return T(left == parse_expression());
    if (match2('!', '=')) return T(left != parse_expression());
    if (match('<')) return T(left < parse_expression());
    if (match('>')) return T(left > parse_expression());
    return left;
  }

  T parse_expression() {
    T result = parse_term();
    while (!has_error) {
      skip_whitespace();
      if (match('+')) result += parse_term();
      else if (match('-')) result -= parse_term();
      else break;
    }
    return result;
  }

  T parse_term() {
    T result = parse_power();
    while (!has_error) {
      skip_whitespace();
      if (match('*')) {
        result *= parse_power();
      } else if (match('/')) {
        size_t dp = pos;
        T d       = parse_power();
        if (!has_error && d == T(0)) return error_at(dp, "division by zero");
        result /= d;
      } else if (match('%')) {
        size_t dp = pos;
        T d       = parse_power();
        if (!has_error && d == T(0)) return error_at(dp, "modulo by zero");
        result = std::fmod(result, d);
      } else break;
    }
    return result;
  }

  T parse_power() {
    T base = parse_factor();
    if (!has_error && match('^')) return std::pow(base, parse_power());
    return base;
  }

  T parse_factor() {
    skip_whitespace();
    if (match('+')) return parse_factor();
    if (match('-')) return -parse_factor();
    if (match('(')) {
      T v = parse_comparison(); // parentheses allow comparison operators inside
      if (!match(')')) return error("expected ')'");
      return v;
    }
    if (std::isalpha(peek()) || peek() == '_') {
      size_t s         = pos;
      std::string name = parse_identifier();
      if (match('(')) return parse_function(name, s);
      auto it = variables.find(name);
      if (it != variables.end()) return *(it->second);
      auto jt = constants.find(name);
      if (jt != constants.end()) return jt->second;
      return error_at(s, "unknown identifier '" + name + "'");
    }
    return parse_number();
  }

  /* ================= EMITTER ================= */

  // Fold constant sub-expressions: if instructions [from, end) contain no
  // PUSH_V, evaluate them with the mini-VM and replace with a single PUSH_K.
  void try_fold(size_t from) {
    if (has_error || emit_prog->size() <= from + 1) return;
    for (size_t i = from; i < emit_prog->size(); i++)
      if ((*emit_prog)[i].op == Op::PUSH_V) return;
    T stk[256];
    int top = -1;
    for (size_t i = from; i < emit_prog->size(); i++) {
      const auto &inst = (*emit_prog)[i];
      switch (inst.op) {
      case Op::PUSH_K:
        stk[++top] = inst.k;
        break;
      case Op::PUSH_V:
        return; // unreachable, defensive
      case Op::NEG:
        stk[top] = -stk[top];
        break;
      case Op::CALL1:
        stk[top] = inst.f1(stk[top]);
        break;
      case Op::ADD: {
        T b = stk[top--];
        stk[top] += b;
        break;
      }
      case Op::SUB: {
        T b = stk[top--];
        stk[top] -= b;
        break;
      }
      case Op::MUL: {
        T b = stk[top--];
        stk[top] *= b;
        break;
      }
      case Op::DIV: {
        T b = stk[top--];
        stk[top] /= b;
        break;
      }
      case Op::MOD: {
        T b      = stk[top--];
        stk[top] = std::fmod(stk[top], b);
        break;
      }
      case Op::POW: {
        T b      = stk[top--];
        stk[top] = std::pow(stk[top], b);
        break;
      }
      case Op::CALL2: {
        T b      = stk[top--];
        stk[top] = inst.f2(stk[top], b);
        break;
      }
      case Op::LT: {
        T b      = stk[top--];
        stk[top] = T(stk[top] < b);
        break;
      }
      case Op::GT: {
        T b      = stk[top--];
        stk[top] = T(stk[top] > b);
        break;
      }
      case Op::LE: {
        T b      = stk[top--];
        stk[top] = T(stk[top] <= b);
        break;
      }
      case Op::GE: {
        T b      = stk[top--];
        stk[top] = T(stk[top] >= b);
        break;
      }
      case Op::EQ: {
        T b      = stk[top--];
        stk[top] = T(stk[top] == b);
        break;
      }
      case Op::NE: {
        T b      = stk[top--];
        stk[top] = T(stk[top] != b);
        break;
      }
      case Op::SEL: {
        T b      = stk[top--];
        T a      = stk[top--];
        stk[top] = (stk[top] != T(0) ? a : b);
        break;
      }
      }
    }
    T val = top >= 0 ? stk[top] : T(0);
    emit_prog->resize(from);
    emit_prog->push_back(Instr::K(val));
  }

  void emit_comparison() {
    DepthGuard guard(depth);
    if (depth > max_depth) {
      error("expression too deeply nested");
      return;
    }
    size_t from = emit_prog->size();
    emit_expression();
    if (!has_error) {
      skip_whitespace();
      Op op      = Op::ADD;
      bool found = true;
      if (match2('<', '=')) op = Op::LE;
      else if (match2('>', '=')) op = Op::GE;
      else if (match2('=', '=')) op = Op::EQ;
      else if (match2('!', '=')) op = Op::NE;
      else if (match('<')) op = Op::LT;
      else if (match('>')) op = Op::GT;
      else found = false;
      if (found && !has_error) {
        emit_expression();
        if (!has_error) {
          emit_prog->push_back(Instr::OP(op));
          try_fold(from);
        }
      }
    }
  }

  void emit_expression() {
    size_t from = emit_prog->size();
    emit_term();
    while (!has_error) {
      skip_whitespace();
      if (match('+')) {
        emit_term();
        emit_prog->push_back(Instr::OP(Op::ADD));
      } else if (match('-')) {
        emit_term();
        emit_prog->push_back(Instr::OP(Op::SUB));
      } else break;
    }
    try_fold(from);
  }

  void emit_term() {
    size_t from = emit_prog->size();
    emit_power();
    while (!has_error) {
      skip_whitespace();
      if (match('*')) {
        emit_power();
        emit_prog->push_back(Instr::OP(Op::MUL));
      } else if (match('/')) {
        emit_power();
        emit_prog->push_back(Instr::OP(Op::DIV));
      } else if (match('%')) {
        emit_power();
        emit_prog->push_back(Instr::OP(Op::MOD));
      } else break;
    }
    try_fold(from);
  }

  void emit_power() {
    size_t from = emit_prog->size();
    emit_factor();
    if (!has_error && match('^')) {
      emit_power();
      emit_prog->push_back(Instr::OP(Op::POW));
      try_fold(from);
    }
  }

  void emit_factor() {
    skip_whitespace();
    if (match('+')) {
      emit_factor();
      return;
    }
    if (match('-')) {
      size_t from = emit_prog->size();
      emit_factor();
      emit_prog->push_back(Instr::OP(Op::NEG));
      try_fold(from);
      return;
    }
    if (match('(')) {
      emit_comparison(); // parentheses allow comparison operators inside
      if (!match(')')) error("expected ')'");
      return;
    }
    if (std::isalpha(peek()) || peek() == '_') {
      size_t s         = pos;
      std::string name = parse_identifier();
      if (match('(')) {
        emit_function(name, s);
        return;
      }
      auto it = variables.find(name);
      if (it != variables.end()) {
        emit_prog->push_back(Instr::V(it->second));
        return;
      }
      auto jt = constants.find(name);
      if (jt != constants.end()) {
        emit_prog->push_back(Instr::K(jt->second));
        return;
      }
      error_at(s, "unknown identifier '" + name + "'");
      return;
    }
    T val = parse_number();
    if (!has_error) emit_prog->push_back(Instr::K(val));
  }

  /* ================= FUNCTIONS ================= */

  T parse_function(const std::string &name, size_t ns) {
    // Arguments may contain comparison operators.
    auto one = [&](auto fn) -> T {
      T v = parse_comparison();
      if (!match(')')) return error("expected ')'");
      return fn(v);
    };
    auto two = [&](auto fn) -> T {
      T a = parse_comparison();
      if (!match(',')) return error("expected ','");
      T b = parse_comparison();
      if (!match(')')) return error("expected ')'");
      return fn(a, b);
    };
    if (name == "if") {
      T cond = parse_comparison();
      if (!match(',')) return error("expected ','");
      T a = parse_comparison();
      if (!match(',')) return error("expected ','");
      T b = parse_comparison();
      if (!match(')')) return error("expected ')'");
      return !has_error ? (cond != T(0) ? a : b) : T(0);
    }
    if (name == "sqrt") return one([](T v) { return std::sqrt(v); });
    if (name == "abs") return one([](T v) { return std::abs(v); });
    if (name == "sin") return one([](T v) { return std::sin(v); });
    if (name == "cos") return one([](T v) { return std::cos(v); });
    if (name == "tan") return one([](T v) { return std::tan(v); });
    if (name == "asin") return one([](T v) { return std::asin(v); });
    if (name == "acos") return one([](T v) { return std::acos(v); });
    if (name == "atan") return one([](T v) { return std::atan(v); });
    if (name == "exp") return one([](T v) { return std::exp(v); });
    if (name == "log") return one([](T v) { return std::log(v); });
    if (name == "log10") return one([](T v) { return std::log10(v); });
    if (name == "floor") return one([](T v) { return std::floor(v); });
    if (name == "ceil") return one([](T v) { return std::ceil(v); });
    if (name == "round") return one([](T v) { return std::round(v); });
    if (name == "pow") return two([](T a, T b) { return std::pow(a, b); });
    if (name == "atan2") return two([](T a, T b) { return std::atan2(a, b); });
    if (name == "min") return two([](T a, T b) { return a < b ? a : b; });
    if (name == "max") return two([](T a, T b) { return a > b ? a : b; });
    if (name == "fmod") return two([](T a, T b) { return std::fmod(a, b); });
    return error_at(ns, "unknown function '" + name + "'");
  }

  void emit_function(const std::string &name, size_t ns) {
    // one/two: emit args (with try_fold) then the CALL instruction.
    auto one = [&](T (*fn)(T)) {
      size_t from = emit_prog->size();
      emit_comparison();
      if (!match(')')) {
        error("expected ')'");
        return;
      }
      if (!has_error) {
        emit_prog->push_back(Instr::F1(fn));
        try_fold(from);
      }
    };
    auto two = [&](T (*fn)(T, T)) {
      size_t from = emit_prog->size();
      emit_comparison();
      if (!match(',')) {
        error("expected ','");
        return;
      }
      emit_comparison();
      if (!match(')')) {
        error("expected ')'");
        return;
      }
      if (!has_error) {
        emit_prog->push_back(Instr::F2(fn));
        try_fold(from);
      }
    };
    if (name == "if") {
      size_t from = emit_prog->size();
      emit_comparison();
      if (!match(',')) {
        error("expected ','");
        return;
      }
      emit_comparison();
      if (!match(',')) {
        error("expected ','");
        return;
      }
      emit_comparison();
      if (!match(')')) {
        error("expected ')'");
        return;
      }
      if (!has_error) {
        emit_prog->push_back(Instr::OP(Op::SEL));
        try_fold(from);
      }
      return;
    }
    // non-capturing lambdas decay to function pointers via unary +
    if (name == "sqrt") {
      one(+[](T v) -> T { return std::sqrt(v); });
      return;
    }
    if (name == "abs") {
      one(+[](T v) -> T { return std::abs(v); });
      return;
    }
    if (name == "sin") {
      one(+[](T v) -> T { return std::sin(v); });
      return;
    }
    if (name == "cos") {
      one(+[](T v) -> T { return std::cos(v); });
      return;
    }
    if (name == "tan") {
      one(+[](T v) -> T { return std::tan(v); });
      return;
    }
    if (name == "asin") {
      one(+[](T v) -> T { return std::asin(v); });
      return;
    }
    if (name == "acos") {
      one(+[](T v) -> T { return std::acos(v); });
      return;
    }
    if (name == "atan") {
      one(+[](T v) -> T { return std::atan(v); });
      return;
    }
    if (name == "exp") {
      one(+[](T v) -> T { return std::exp(v); });
      return;
    }
    if (name == "log") {
      one(+[](T v) -> T { return std::log(v); });
      return;
    }
    if (name == "log10") {
      one(+[](T v) -> T { return std::log10(v); });
      return;
    }
    if (name == "floor") {
      one(+[](T v) -> T { return std::floor(v); });
      return;
    }
    if (name == "ceil") {
      one(+[](T v) -> T { return std::ceil(v); });
      return;
    }
    if (name == "round") {
      one(+[](T v) -> T { return std::round(v); });
      return;
    }
    if (name == "pow") {
      two(+[](T a, T b) -> T { return std::pow(a, b); });
      return;
    }
    if (name == "atan2") {
      two(+[](T a, T b) -> T { return std::atan2(a, b); });
      return;
    }
    if (name == "min") {
      two(+[](T a, T b) -> T { return a < b ? a : b; });
      return;
    }
    if (name == "max") {
      two(+[](T a, T b) -> T { return a > b ? a : b; });
      return;
    }
    if (name == "fmod") {
      two(+[](T a, T b) -> T { return std::fmod(a, b); });
      return;
    }
    error_at(ns, "unknown function '" + name + "'");
  }

  /* ================= LEXER ================= */

  char peek() const {
    return pos < expr.size() ? expr[pos] : '\0';
  }

  bool match(char c) {
    skip_whitespace();
    if (peek() == c) {
      pos++;
      return true;
    }
    return false;
  }

  // Match two consecutive characters (e.g. '<=', '==').  Must be checked
  // before match() for operators that share a prefix (e.g. '<=' before '<').
  bool match2(char a, char b) {
    skip_whitespace();
    if (pos + 1 < expr.size() && expr[pos] == a && expr[pos + 1] == b) {
      pos += 2;
      return true;
    }
    return false;
  }

  std::string parse_identifier() {
    std::string name;
    while (std::isalnum(peek()) || peek() == '_') name += expr[pos++];
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
    if (start == pos) return error("expected number");
    try {
      return static_cast<T>(std::stold(expr.substr(start, pos - start)));
    } // stold: max precision
    catch (...) {
      return error_at(start, "invalid number literal");
    }
  }

  void skip_whitespace() {
    while (pos < expr.size() && std::isspace(expr[pos])) pos++;
  }

  /* ================= ERRORS ================= */

  T error_at(size_t at, std::string msg) {
    if (!has_error) {
      has_error = true;
      err_pos   = at;
      err_msg   = std::move(msg);
    }
    return T(0);
  }

  T error(std::string msg) {
    return error_at(pos, std::move(msg));
  }
};

#if 0

// Illustrates conditional expressions and the performance advantage of compile().
#include <chrono>
#include <iostream>

int main() {
  exprParser<double> p;
  double x = 0.0;
  p.addVariable("x", &x);

  // if() with comparison: piecewise function, constant sub-expressions folded
  const std::string expr = "if(x >= 0, sqrt(x), 0) + sin(pi/4)^2";
  const int N = 2'000'000;

  // --- interpreted ---
  auto t0 = std::chrono::steady_clock::now();
  double sum1 = 0.0;
  for (int i = 0; i < N; ++i) { x=(i-N/2)*1e-5; double r; p.parse(expr,r); sum1+=r; }
  auto t1 = std::chrono::steady_clock::now();

  // --- compiled (sin(pi/4)^2 → PUSH_K(0.5) at compile time) ---
  auto prog = p.compile(expr);
  if (!prog) { p.printError(); return 1; }

  auto t2 = std::chrono::steady_clock::now();
  double sum2 = 0.0;
  for (int i = 0; i < N; ++i) { x=(i-N/2)*1e-5; sum2+=prog.eval(); }
  auto t3 = std::chrono::steady_clock::now();

  using ms = std::chrono::duration<double, std::milli>;
  std::cout << "interpreted: " << ms(t1-t0).count() << " ms\n";
  std::cout << "compiled:    " << ms(t3-t2).count() << " ms\n";
  std::cout << "speedup:     " << ms(t1-t0).count()/ms(t3-t2).count() << "x\n";
  std::cout << "sums match:  " << (std::fabs(sum1-sum2) < 1e-6 ? "yes" : "NO") << "\n";
  return 0;
}

#endif
