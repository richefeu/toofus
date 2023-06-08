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

// Credits
// This class is issued from a tutorial course of Krishna Kumar
// We acknowledge and are grateful to this developer for his contributions.

#ifndef FACTORY_HPP
#define FACTORY_HPP

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

/// The factory - implements singleton pattern!
template <class BaseClass, typename Key = std::string> class Factory {
public:
  /// Get the single instance of the factory
  static Factory *Instance() {
    static Factory factory;
    return &factory;
  }

  /// register a factory function to create an instance of className
  void RegisterFactoryFunction(Key key, std::function<BaseClass *(void)> classFactoryFunction) {
    // register the class factory function
    factoryFunctionRegistry[key] = classFactoryFunction;
  }

  /// create an instance of a registered class
  BaseClass *Create(Key key) {
    // find name in the registry and call factory method.
    auto it = factoryFunctionRegistry.find(key);
    if (it != factoryFunctionRegistry.end())
      return it->second();
    return nullptr;
  }

private:
  /// a private ctor
  Factory() {}

  /// the registry of factory functions
  std::map<Key, std::function<BaseClass *(void)>> factoryFunctionRegistry;
};

/// A helper class to register a factory function
template <class BaseClass, class DerivedClass, typename Key = std::string> class Registrar {
public:
  explicit Registrar(Key key) {
    // register the class factory function
    Factory<BaseClass, Key>::Instance()->RegisterFactoryFunction(
        key, [](void) -> BaseClass * { return new DerivedClass(); });
  }
};

#define REGISTRER_BASE_DERIVED(BASE_CLASS, DERIVED_CLASS)                                                              \
  Factory<BASE_CLASS, std::string>::Instance()->RegisterFactoryFunction(                                               \
      #DERIVED_CLASS, [](void) -> BASE_CLASS * { return new DERIVED_CLASS(); })

#endif /* end of include guard: FACTORY_HPP */

#if 0

#include <iostream>

// Base class
class Product {
public:
  virtual void doSomething() = 0;
};

// Derived class 1
class ConcreteProduct1 : public Product {
public:
  void doSomething() override { std::cout << "ConcreteProduct1 doing something." << std::endl; }
};

// Derived class 2
class ConcreteProduct2 : public Product {
public:
  void doSomething() override { std::cout << "ConcreteProduct2 doing something." << std::endl; }
};

int main() {

  REGISTRER_BASE_DERIVED(Product, ConcreteProduct1);
  REGISTRER_BASE_DERIVED(Product, ConcreteProduct2);

  // Create instances of registered classes using the factory
  Product *product1 = Factory<Product>::Instance()->Create("ConcreteProduct1");
  Product *product2 = Factory<Product>::Instance()->Create("ConcreteProduct2");

  // Call the doSomething() method of the created objects
  if (product1)
    product1->doSomething();

  if (product2)
    product2->doSomething();

  return 0;
}

#endif
