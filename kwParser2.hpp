// STATUS: [ ] STABLE  [x] EXPERIMENTAL  [ ] DRAFT

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

#pragma once

#include <fstream>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <string>
#include <string_view>
#include <limits>

// ---------------------------------------------------------------------------
// Helper : remplace les macros __DO__ / __GET__ dont les noms étaient réservés
// par le standard C++ (double underscore = UB).
//
//   Avant : parser.kwMap["k"] = __DO__(is) { is >> v; };
//   Après : parser.kwMap["k"] = [&](std::istream& is) { is >> v; };
//
//   Avant : parser.kwMap["k"] = __GET__(is, v);
//   Après : parser.kwMap["k"] = kwGet(v);
// ---------------------------------------------------------------------------
template <typename T>
auto kwGet(T &target) {
  return [&target](std::istream &is) { is >> target; };
}

// ---------------------------------------------------------------------------

class kwParser {
public:
  // -------------------------------------------------------------------------
  // Types
  // -------------------------------------------------------------------------
  using Action = std::function<void(std::istream &)>;

  // -------------------------------------------------------------------------
  // Constructeurs
  // -------------------------------------------------------------------------
  kwParser() = default;
  explicit kwParser(bool warn) : warn_(warn) {}

  // -------------------------------------------------------------------------
  // Enregistrement de mots-clés
  //
  // On accepte std::string_view pour éviter une copie inutile lors de l'appel
  // avec un littéral ou un std::string. La clé est stockée en std::string dans
  // la map (std::map n'accepte pas string_view comme clé hétérogène sans
  // comparateur dédié, qu'on évite ici pour rester simple).
  // -------------------------------------------------------------------------
  void addKw(std::string_view kw, Action func) {
    kwMap_.emplace(std::string(kw), std::move(func));
  }

  // -------------------------------------------------------------------------
  // parse(filename) – ouvre un fichier et délègue à parse(istream&)
  // -------------------------------------------------------------------------
  void parse(std::string_view filename) {
    std::ifstream file{std::string(filename)};
    if (file) {
      parse(file);
    } else {
      std::cerr << "@kwParser::parse, cannot open file: " << filename << '\n';
    }
  }

  // -------------------------------------------------------------------------
  // parseString(chain) – crée un istringstream et délègue à parse(istream&)
  // -------------------------------------------------------------------------
  void parseString(std::string_view chain) {
    std::istringstream ss{std::string(chain)};
    parse(ss);
  }

  // -------------------------------------------------------------------------
  // parse(istream&) – cœur du parseur
  //
  // Lit le flux token par token :
  //   - '/', '#', '!' en début de token  → ignore le reste de la ligne
  //   - breakStr_                         → arrêt anticipé
  //   - mot-clé connu                     → appel de l'action associée
  //   - mot-clé inconnu                   → warning optionnel
  // -------------------------------------------------------------------------
  void parse(std::istream &is) {
    std::string token;
    while (is >> token) {
      // Token vide : ne devrait pas arriver avec operator>>, mais on protège
      // quand même token[0] contre un accès hors-bornes (UB).
      if (token.empty()) {
        continue;
      }

      // Commentaire : ignorer jusqu'à la fin de la ligne
      if (token[0] == '/' || token[0] == '#' || token[0] == '!') {
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        continue;
      }

      // Mot-clé d'arrêt anticipé
      if (token == breakStr_) {
        break;
      }

      // Recherche et exécution de l'action — syntaxe if-init C++17
      if (auto it = kwMap_.find(token); it != kwMap_.end()) {
        it->second(is);
      } else if (warn_) {
        std::cout << "@kwParser::parse, unknown token: " << token << '\n';
      }
    }
  }

  // -------------------------------------------------------------------------
  // Accès direct à la map (rétro-compatibilité avec l'API publique d'origine)
  //
  //   parser.kwMap["key"] = [&](std::istream& is) { is >> value; };
  //
  // On expose une référence sur la map interne pour conserver exactement la
  // même syntaxe d'utilisation qu'avant.
  // -------------------------------------------------------------------------
  //std::map<std::string, Action> &kwMap = kwMap_;
  std::unordered_map<std::string, Action> &kwMap = kwMap_;

  // -------------------------------------------------------------------------
  // Configuration publique
  // -------------------------------------------------------------------------
  std::string breakStr_{"EOF"};
  bool        warn_{true};

private:
  //std::map<std::string, Action> kwMap_;
  std::unordered_map<std::string, Action> kwMap_;
};

// ---------------------------------------------------------------------------
// EXEMPLE D'UTILISATION
// ---------------------------------------------------------------------------
/**
@file

Fichier texte à parser :
@code
# Données à lire
myClass.value  123.456
myClass.value2 654.321
myClass.str    coucou
@endcode
*/

#if 0
struct myClass {
  double      value{};
  double      value2{};
  std::string str;
} mc;

int main() {
  kwParser parser(/*warn=*/false);
  //parser.breakStr_ = "@endcode";

  // Syntaxe lambda explicite (remplace __DO__)
  parser.kwMap["myClass.value"]  = [&](std::istream &is) { is >> mc.value; };

  // Helper kwGet<T> (remplace __GET__)
  parser.kwMap["myClass.value2"] = kwGet(mc.value2);

  // Via addKw
  parser.addKw("myClass.str", [&](std::istream &is) { is >> mc.str; });

  parser.parse("kwParser2.hpp");

  std::cout << "value  = " << mc.value  << '\n';
  std::cout << "value2 = " << mc.value2 << '\n';
  std::cout << "str    = " << mc.str    << '\n';

  parser.parseString("myClass.str toto");
  std::cout << "str    = " << mc.str    << '\n';

  parser.parseString("myClass.value2 6.1 myClass.value 13.6");
  std::cout << "value  = " << mc.value  << '\n';
  std::cout << "value2 = " << mc.value2 << '\n';

  return 0;
}
#endif