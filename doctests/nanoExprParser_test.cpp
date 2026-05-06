#include "doctest.h"
#include "nanoExprParser2.hpp"

TEST_CASE("nanoExprParser") {

/* ---------------- BASIC OPERATIONS ---------------- */

SUBCASE("Addition") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("1+2", result));
    CHECK(result == doctest::Approx(3));
}

SUBCASE("Subtraction") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("5-3", result));
    CHECK(result == doctest::Approx(2));
}

SUBCASE("Multiplication") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("4*3", result));
    CHECK(result == doctest::Approx(12));
}

SUBCASE("Division") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("10/2", result));
    CHECK(result == doctest::Approx(5));
}

/* ---------------- PRECEDENCE ---------------- */

SUBCASE("Operator precedence") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("2+3*4", result));
    CHECK(result == doctest::Approx(14));
}

SUBCASE("Parentheses override precedence") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("(2+3)*4", result));
    CHECK(result == doctest::Approx(20));
}

/* ---------------- POWER ---------------- */

SUBCASE("Power operator") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("2^3", result));
    CHECK(result == doctest::Approx(8));
}

/* ---------------- CONSTANTS ---------------- */

SUBCASE("Pi constant") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("pi", result));
    CHECK(result == doctest::Approx(3.141592653589793));
}

SUBCASE("E constant") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("e", result));
    CHECK(result == doctest::Approx(2.718281828459045));
}

/* ---------------- VARIABLES ---------------- */

SUBCASE("Variable usage") {
    nanoExprParser<double> parser;
    double result;
    double x = 5;
    parser.addVariable("x", &x);

    CHECK(parser.parse("x*2", result));
    CHECK(result == doctest::Approx(10));
}

SUBCASE("Multiple variables") {
    nanoExprParser<double> parser;
    double result;
    double x = 2, y = 3;
    parser.addVariable("x", &x);
    parser.addVariable("y", &y);

    CHECK(parser.parse("x+y*2", result));
    CHECK(result == doctest::Approx(8));
}

/* ---------------- FUNCTIONS ---------------- */

SUBCASE("Sqrt function") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("sqrt(9)", result));
    CHECK(result == doctest::Approx(3));
}

SUBCASE("Sin function") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("sin(0)", result));
    CHECK(result == doctest::Approx(0));
}

SUBCASE("Cos function") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("cos(0)", result));
    CHECK(result == doctest::Approx(1));
}

SUBCASE("Log function") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("log(e)", result));
    CHECK(result == doctest::Approx(1));
}

SUBCASE("Log10 function") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("log10(100)", result));
    CHECK(result == doctest::Approx(2));
}

SUBCASE("Exp function") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("exp(1)", result));
    CHECK(result == doctest::Approx(std::exp(1)));
}

SUBCASE("Pow function") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("pow(2,3)", result));
    CHECK(result == doctest::Approx(8));
}

/* ---------------- FLOATING POINT ---------------- */

SUBCASE("Floating point number") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("3.14", result));
    CHECK(result == doctest::Approx(3.14));
}

SUBCASE("Scientific notation") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("1e3", result));
    CHECK(result == doctest::Approx(1000));
}

/* ---------------- UNARY OPERATORS ---------------- */

SUBCASE("Unary minus") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("-5+3", result));
    CHECK(result == doctest::Approx(-2));
}

SUBCASE("Unary plus") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("+5", result));
    CHECK(result == doctest::Approx(5));
}

/* ---------------- EDGE CASES ---------------- */

SUBCASE("Division by zero") {
    nanoExprParser<double> parser;
    double result;
    CHECK_FALSE(parser.parse("5/0", result));
}

SUBCASE("Invalid expression") {
    nanoExprParser<double> parser;
    double result;
    CHECK_FALSE(parser.parse("2+", result));
}

SUBCASE("Unknown variable") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("unknown", result));
    CHECK(result == doctest::Approx(0));
}

SUBCASE("Unknown function") {
    nanoExprParser<double> parser;
    double result;
    CHECK_FALSE(parser.parse("foo(2)", result));
}

SUBCASE("Missing parenthesis") {
    nanoExprParser<double> parser;
    double result;
    CHECK_FALSE(parser.parse("(2+3", result));
}

/* ---------------- COMPLEX EXPRESSIONS ---------------- */

SUBCASE("Complex expression") {
    nanoExprParser<double> parser;
    double result;
    CHECK(parser.parse("sqrt(16) + pow(2,3) * (3+1)", result));
    CHECK(result == doctest::Approx(36));
}

}