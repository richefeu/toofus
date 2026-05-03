#ifndef NAMEDPARAMETERS_HPP
#define NAMEDPARAMETERS_HPP

#include <any>
#include <iostream>
#include <optional>
#include <string>
#include <typeinfo>
#include <unordered_map>

#define WITH_NP(name)                                                                                                  \
  auto name = NamedParameters {}

class NamedParameters {
private:
  std::unordered_map<std::string, std::any> params;

public:
  // Add a parameter
  template <typename T> NamedParameters &add(const std::string &name, const T &value) {
    params[name] = value;
    return *this;
  }

  // get a parameter
  template <typename T> std::optional<T> get(const std::string &name) const {
    auto it = params.find(name);
    if (it == params.end()) {
      return std::nullopt; //  missing parameter
    }

    // Check the type before casting
    if (it->second.type() != typeid(T)) {
      return std::nullopt; // wrong type
    }

    // if correct, we securly cast
    return std::any_cast<T>(it->second);
  }

  // Chech wether a parameter exists
  bool has(const std::string &name) const {
    return params.find(name) != params.end();
  }
};

#endif /* end of include guard: NAMEDPARAMETERS_HPP */

// TESTS
#if 0

#include <iostream>

// --- Fonction exemple : création d'une fenêtre ---
void createWindow(const NamedParameters &params) {
  int width   = params.get<int>("width").value_or(800);
  auto height = params.get<int>("height");
  auto title  = params.get<std::string>("title");

  if (!height || !title) {
    std::cout << "Erreur : paramètres manquants ou invalides pour createWindow.\n";
    return;
  }

  std::cout << "Fenêtre créée : " << width << "x" << *height << ", titre : \"" << *title << "\"\n";
}

// --- Fonction exemple : affichage d'un point ---
void printPoint(const NamedParameters &params) {
  auto x = params.get<int>("x");
  auto y = params.get<int>("y");

  if (!x || !y) {
    std::cout << "Erreur : paramètres 'x' ou 'y' manquants ou invalides.\n";
    return;
  }

  std::cout << "Point : (" << *x << ", " << *y << ")\n";
}

// --- Fonction principale ---
int main() {
  std::cout << "=== Test 1 : Création d'une fenêtre (OK) ===\n";
  auto windowParams = NamedParameters{} //.add("width", 800)
                          .add("height", 600)
                          .add("title", std::string("Ma Fenêtre"));
  createWindow(windowParams);

  std::cout << "\n=== Test 2 : Point (OK) ===\n";
  auto pointParams = NamedParameters{}.add("x", 10).add("y", 20);
  printPoint(pointParams);

  std::cout << "\n=== Test 3 : Paramètre manquant ===\n";
  auto incompleteParams = NamedParameters{}.add("x", 5);
  printPoint(incompleteParams); // Affiche une erreur

  std::cout << "\n=== Test 4 : Mauvais type (string au lieu de int) ===\n";
  auto badTypeParams = NamedParameters{}
                           .add("width", std::string("800")) // "width" est un string, pas un int
                           .add("height", 600)
                           .add("title", std::string("Test"));
  createWindow(badTypeParams); // Affiche une erreur

  std::cout << "\n=== Test 5 : Vérification de l'existence ===\n";
  auto testParams = NamedParameters{}.add("name", std::string("Vincent")).add("age", 30);
  std::cout << "Le paramètre 'name' existe ? " << (testParams.has("name") ? "Oui" : "Non") << "\n";
  std::cout << "Le paramètre 'address' existe ? " << (testParams.has("address") ? "Oui" : "Non") << "\n";

  std::cout << "\n=== Test 6 : Chaînage des appels ===\n";
  NamedParameters chainedParams;
  chainedParams.add("a", 1).add("b", 2).add("c", 3);
  std::cout << "a = " << *chainedParams.get<int>("a") << "\n";
  std::cout << "b = " << *chainedParams.get<int>("b") << "\n";
  std::cout << "c = " << *chainedParams.get<int>("c") << "\n";

  std::cout << "\n=== Test 7 : Valeur par défaut ===\n";
  auto defaultParams = NamedParameters{}.add("width", 800);
  int height         = defaultParams.get<int>("height").value_or(480);
  std::cout << "Hauteur (par défaut si manquant) : " << height << "\n";

  std::cout << "\n=== Test 8 : Utilisation de WITH_NP  ===\n";

  WITH_NP(macroParams).add("width", 400).add("height", 300).add("title", std::string("Ma Fenêtre WITH_NP"));
  createWindow(macroParams);

  std::cout << "\n=== Test 9 : Utilisation avec bloque pour utiliser le même nom de variable  ===\n";

  {
    WITH_NP(macroParams).add("width", 500).add("height", 400).add("title", std::string("Ma Fenêtre same macroParams"));
    createWindow(macroParams);
  }

  std::cout << "\nTous les tests sont terminés !\n";
  return 0;
}

#endif /* END TEST */
