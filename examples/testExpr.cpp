// on mac:
// brew install googletest
// g++ -std=c++17 testExpr.cpp -I/opt/homebrew/include -L/opt/homebrew/lib -lgtest -lgtest_main -pthread -o testExpr

#include <gtest/gtest.h>
#include <cmath>
#include "../nanoExprParser2.hpp"

class NanoExprParserTest : public ::testing::Test {
protected:
    nanoExprParser<double> parser;
    double result;
};

/* ---------------- BASIC OPERATIONS ---------------- */

TEST_F(NanoExprParserTest, Addition) {
    EXPECT_TRUE(parser.parse("1+2", result));
    EXPECT_DOUBLE_EQ(result, 3);
}

TEST_F(NanoExprParserTest, Subtraction) {
    EXPECT_TRUE(parser.parse("5-3", result));
    EXPECT_DOUBLE_EQ(result, 2);
}

TEST_F(NanoExprParserTest, Multiplication) {
    EXPECT_TRUE(parser.parse("4*3", result));
    EXPECT_DOUBLE_EQ(result, 12);
}

TEST_F(NanoExprParserTest, Division) {
    EXPECT_TRUE(parser.parse("10/2", result));
    EXPECT_DOUBLE_EQ(result, 5);
}

/* ---------------- PRECEDENCE ---------------- */

TEST_F(NanoExprParserTest, OperatorPrecedence) {
    EXPECT_TRUE(parser.parse("2+3*4", result));
    EXPECT_DOUBLE_EQ(result, 14);
}

TEST_F(NanoExprParserTest, ParenthesesOverridePrecedence) {
    EXPECT_TRUE(parser.parse("(2+3)*4", result));
    EXPECT_DOUBLE_EQ(result, 20);
}

/* ---------------- POWER ---------------- */

TEST_F(NanoExprParserTest, PowerOperator) {
    EXPECT_TRUE(parser.parse("2^3", result));
    EXPECT_DOUBLE_EQ(result, 8);
}

/* ---------------- CONSTANTS ---------------- */

TEST_F(NanoExprParserTest, PiConstant) {
    EXPECT_TRUE(parser.parse("pi", result));
    EXPECT_NEAR(result, 3.141592653589793, 1e-9);
}

TEST_F(NanoExprParserTest, EConstant) {
    EXPECT_TRUE(parser.parse("e", result));
    EXPECT_NEAR(result, 2.718281828459045, 1e-9);
}

/* ---------------- VARIABLES ---------------- */

TEST_F(NanoExprParserTest, VariableUsage) {
    double x = 5;
    parser.addVariable("x", &x);

    EXPECT_TRUE(parser.parse("x*2", result));
    EXPECT_DOUBLE_EQ(result, 10);
}

TEST_F(NanoExprParserTest, MultipleVariables) {
    double x = 2, y = 3;
    parser.addVariable("x", &x);
    parser.addVariable("y", &y);

    EXPECT_TRUE(parser.parse("x+y*2", result));
    EXPECT_DOUBLE_EQ(result, 8);
}

/* ---------------- FUNCTIONS ---------------- */

TEST_F(NanoExprParserTest, SqrtFunction) {
    EXPECT_TRUE(parser.parse("sqrt(9)", result));
    EXPECT_DOUBLE_EQ(result, 3);
}

TEST_F(NanoExprParserTest, SinFunction) {
    EXPECT_TRUE(parser.parse("sin(0)", result));
    EXPECT_NEAR(result, 0, 1e-9);
}

TEST_F(NanoExprParserTest, CosFunction) {
    EXPECT_TRUE(parser.parse("cos(0)", result));
    EXPECT_NEAR(result, 1, 1e-9);
}

TEST_F(NanoExprParserTest, LogFunction) {
    EXPECT_TRUE(parser.parse("log(e)", result));
    EXPECT_NEAR(result, 1, 1e-9);
}

TEST_F(NanoExprParserTest, Log10Function) {
    EXPECT_TRUE(parser.parse("log10(100)", result));
    EXPECT_NEAR(result, 2, 1e-9);
}

TEST_F(NanoExprParserTest, ExpFunction) {
    EXPECT_TRUE(parser.parse("exp(1)", result));
    EXPECT_NEAR(result, std::exp(1), 1e-9);
}

TEST_F(NanoExprParserTest, PowFunction) {
    EXPECT_TRUE(parser.parse("pow(2,3)", result));
    EXPECT_DOUBLE_EQ(result, 8);
}

/* ---------------- FLOATING POINT ---------------- */

TEST_F(NanoExprParserTest, FloatingPointNumber) {
    EXPECT_TRUE(parser.parse("3.14", result));
    EXPECT_NEAR(result, 3.14, 1e-9);
}

TEST_F(NanoExprParserTest, ScientificNotation) {
    EXPECT_TRUE(parser.parse("1e3", result));
    EXPECT_DOUBLE_EQ(result, 1000);
}

/* ---------------- UNARY OPERATORS ---------------- */

TEST_F(NanoExprParserTest, UnaryMinus) {
    EXPECT_TRUE(parser.parse("-5+3", result));
    EXPECT_DOUBLE_EQ(result, -2);
}

TEST_F(NanoExprParserTest, UnaryPlus) {
    EXPECT_TRUE(parser.parse("+5", result));
    EXPECT_DOUBLE_EQ(result, 5);
}

/* ---------------- EDGE CASES ---------------- */

TEST_F(NanoExprParserTest, DivisionByZero) {
    EXPECT_FALSE(parser.parse("5/0", result));
}

TEST_F(NanoExprParserTest, InvalidExpression) {
    EXPECT_FALSE(parser.parse("2+", result));
}

TEST_F(NanoExprParserTest, UnknownVariable) {
    EXPECT_TRUE(parser.parse("unknown", result));
    EXPECT_DOUBLE_EQ(result, 0); // fallback to parse_number()
}

TEST_F(NanoExprParserTest, UnknownFunction) {
    EXPECT_FALSE(parser.parse("foo(2)", result));
}

TEST_F(NanoExprParserTest, MissingParenthesis) {
    EXPECT_FALSE(parser.parse("(2+3", result));
}

/* ---------------- COMPLEX EXPRESSIONS ---------------- */

TEST_F(NanoExprParserTest, ComplexExpression) {
    EXPECT_TRUE(parser.parse("sqrt(16) + pow(2,3) * (3+1)", result));
    EXPECT_DOUBLE_EQ(result, 4 + 8 * 4); // 36
}