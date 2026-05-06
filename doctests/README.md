# Doctest Usage Guide

This guide explains how to use **Doctest**, a lightweight C++ testing framework, in your project.

---

## 📦 What is Doctest?

Doctest is a fast, single-header C++ testing framework designed to be:

* Easy to integrate (just one header file)
* Very fast to compile
* Simple and expressive

---

## 🚀 Getting Started

### 1. Add Doctest to your project

Download the header file:

```
doctest.h
```

Place it somewhere in your project (e.g., `third_party/doctest/`), then include it:

```cpp
#include "doctest.h"
```

---

### 2. Create a test executable

In **one source file only**, define:

```cpp
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
```

This provides the `main()` function for running tests.

---

## ✍️ Writing Tests

### Basic test case

```cpp
TEST_CASE("Basic arithmetic") {
    CHECK(1 + 1 == 2);
    CHECK(2 * 3 == 6);
}
```

---

### Assertions

Common assertions:

```cpp
CHECK(expr);          // continues on failure
REQUIRE(expr);        // stops test on failure

CHECK_FALSE(expr);
REQUIRE_FALSE(expr);
```

---

### Floating point comparisons

Use `doctest::Approx`:

```cpp
CHECK(result == doctest::Approx(3.14));
```

You can adjust tolerance:

```cpp
CHECK(result == doctest::Approx(3.14).epsilon(1e-6));
```

---

## 🧩 Test Suites

Test suites group related tests.

### Grouping tests

```cpp
TEST_SUITE("Math") {

TEST_CASE("Addition") {
    CHECK(1 + 2 == 3);
}

TEST_CASE("Multiplication") {
    CHECK(2 * 3 == 6);
}

}
```

---

### Tagging individual tests

```cpp
TEST_CASE("Division" * doctest::test_suite("Math")) {
    CHECK(10 / 2 == 5);
}
```

---

## 🔁 Reusing Setup Code

Unlike some frameworks, Doctest does not require fixtures. You can:

### Option 1: Use helper functions

```cpp
int add(int a, int b) { return a + b; }

TEST_CASE("Add function") {
    CHECK(add(2, 3) == 5);
}
```

---

### Option 2: Use `TEST_CASE_FIXTURE`

```cpp
struct Fixture {
    int x = 5;
};

TEST_CASE_FIXTURE(Fixture, "Using fixture") {
    CHECK(x == 5);
}
```

---

### Option 3: Use `SUBCASE`

```cpp
TEST_CASE("Vector tests") {
    std::vector<int> v;

    SUBCASE("Push back") {
        v.push_back(1);
        CHECK(v.size() == 1);
    }

    SUBCASE("Empty initially") {
        CHECK(v.empty());
    }
}
```

Each `SUBCASE` runs independently.

---

## ▶️ Running Tests

Compile your test executable, then run:

```bash
./tests
```

---

### Filter tests

Run only a specific test:

```bash
./tests --test-case="Addition"
```

Run a specific suite:

```bash
./tests --test-suite="Math"
```

Exclude a suite:

```bash
./tests --test-suite-exclude="Edge Cases"
```

---

### List tests

```bash
./tests --list-test-cases
```

---

## 📊 Output Options

Verbose output:

```bash
./tests --success
```

Minimal output:

```bash
./tests --quiet
```

---

## ⚡ Best Practices

* Keep tests small and focused
* Use descriptive test names
* Group related tests with test suites
* Prefer `CHECK` unless failure should stop execution
* Use `SUBCASE` to avoid duplication

---

## 🧪 Example

```cpp
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

int add(int a, int b) {
    return a + b;
}

TEST_SUITE("Math") {

TEST_CASE("Addition works") {
    CHECK(add(2, 3) == 5);
}

TEST_CASE("Negative numbers") {
    CHECK(add(-1, -2) == -3);
}

}
```

---

## 📚 Summary

* **TEST_CASE** defines a test
* **CHECK / REQUIRE** perform assertions
* **TEST_SUITE** groups tests
* **SUBCASE** helps structure scenarios
* No heavy setup required — keep it simple

---

Doctest is designed to stay out of your way while giving you powerful testing capabilities with minimal overhead.
