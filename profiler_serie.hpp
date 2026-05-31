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

// =============================================================================
// profiler_serie.hpp — minuteur hiérarchique pour code série
// =============================================================================
//
// PRINCIPE
//   Les timers forment un arbre : chaque START_TIMER crée un nœud fils du
//   timer actif. La durée d'un nœud cumule le temps de tous ses appels.
//   À la fin du programme, PRINT_TIMERS affiche le tableau et écrit un
//   fichier texte "perflog_<nom>.txt".
//
// ACTIVATION
//   Compiler avec -DENABLE_PROFILING pour activer les mesures.
//   Sans ce flag, toutes les macros sont des no-ops (zéro overhead).
//
// MACROS
//
//   INIT_TIMERS()
//     À placer une seule fois au début du bloc de mesure (typiquement main()).
//     Crée le nœud racine et démarre le chronomètre global.
//     Déclare la variable locale 'roottimer' utilisée par PRINT_TIMERS.
//
//   START_TIMER("nom")
//     Démarre un timer nommé "nom" comme fils du timer courant.
//     S'arrête automatiquement à la sortie du bloc (RAII).
//     Peut être imbriqué à volonté.
//     Déclare une variable locale : ne pas placer deux START_TIMER
//     au même niveau de portée dans le même bloc sans accolades.
//
//   PRINT_TIMERS("nom_base")
//     Arrête le timer racine, écrit "perflog_nom_base.txt" et affiche
//     le tableau trié par durée décroissante.
//     À placer à la fin du bloc ouvert par INIT_TIMERS().
//
//   WRITE_TIMERS("nom_base")
//     Comme PRINT_TIMERS mais sans arrêter le timer racine.
//     Utile pour des sorties intermédiaires en cours de calcul.
//
// EXEMPLE
//
//   #define ENABLE_PROFILING
//   #include "profiler_serie.hpp"
//
//   void assemble() {
//     START_TIMER("assemble");
//     // ...
//   }
//
//   void solve() {
//     START_TIMER("solve");
//     {
//       START_TIMER("factorize");
//       // ...
//     }
//     {
//       START_TIMER("back_sub");
//       // ...
//     }
//   }
//
//   int main() {
//     INIT_TIMERS();
//     for (int i = 0; i < 100; i++) {
//       assemble();
//       solve();
//     }
//     PRINT_TIMERS("mysim");
//   }
//
// TABLEAU DE SORTIE
//   Colonnes : calls | total (s) | average (s) | stddev (s)
//   - total   : temps cumulé de tous les appels
//   - average : total / calls
//   - stddev  : écart-type inter-appels (algorithme de Welford, stable)
//   Les enfants sont triés par durée décroissante à chaque niveau.
//
// =============================================================================

#ifndef PROFILER_SERIE_HPP
#define PROFILER_SERIE_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace profiler_detail {

static const size_t nColumns             = 4;
static const size_t cWidths[nColumns]    = {9, 14, 14, 14}; // calls, total, average, stddev
static const std::string cName[nColumns] = {"calls", "total", "average", "stddev"};

// Always produces a fixed-width 10-char string: "X.XXXXe+YY"
static std::string fmtTime(double v) {
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%.4e", v);
  return buf;
}

class Node {
public:
  Node()
      : m_name("root"), m_iteration(1), m_level(0), m_mother(nullptr), m_duration(0), m_welford_mean(0.0),
        m_welford_M2(0.0) {}

  Node(const std::string &name, Node *mother)
      : m_name(name), m_iteration(1), m_level(mother->m_level + 1), m_mother(mother), m_duration(0),
        m_welford_mean(0.0), m_welford_M2(0.0) {}

  ~Node() {
    for (auto d : m_daughter) delete d;
  }

  Node *find(const std::string &name) {
    for (auto it : m_daughter) {
      if (it->m_name == name) {
        it->m_iteration++;
        return it;
      }
    }
    Node *child = new Node(name, this);
    m_daughter.push_back(child);
    return child;
  }

  void printReplicate(size_t begin, size_t end, const std::string &motif) {
    for (size_t i = begin; i < end; i++) std::cout << motif;
  }

  // Total visible width of a content row
  size_t totalWidth(size_t shift) const {
    size_t w = shift + 1; // name column + closing │
    for (size_t i = 0; i < nColumns; i++) w += cWidths[i] + 1;
    return w;
  }

  void printBanner(size_t shift) {
    if (m_name != "root") return;
    size_t tw = totalWidth(shift);
    std::cout << " ┌";
    printReplicate(2, tw - 1, "─");
    std::cout << "┐\n";
    // " │    name" = 10 visible columns
    std::cout << " │    name";
    printReplicate(10, shift, " ");
    for (size_t i = 0; i < nColumns; i++) {
      std::cout << "│";
      size_t size = cName[i].size();
      printReplicate(0, cWidths[i] - size - 1, " ");
      std::cout << cName[i] << " ";
    }
    std::cout << "│\n";
    std::cout << " │";
    printReplicate(2, tw - 1, "─");
    std::cout << "│\n";
  }

  void printEnding(size_t shift) {
    if (m_name != "root") return;
    size_t tw = totalWidth(shift);
    std::cout << " └";
    printReplicate(2, tw - 1, "─");
    std::cout << "┘\n";
  }

  void print(size_t shift) {
    size_t realShift    = shift;
    size_t currentShift = 3;
    std::cout << " │ ";
    for (int i = 0; i < int(m_level) - 1; i++) {
      std::cout << "   ";
      currentShift += 3;
    }
    if (m_level > 0) {
      std::cout << "└──";
      currentShift += 3;
    }
    std::cout << "► " << m_name;
    currentShift += m_name.size() + 2;
    printReplicate(currentShift, realShift, " ");

    double n       = double(m_iteration);
    double total   = m_duration.count();
    double mean    = total / n;
    double std_dev = std::sqrt(m_welford_M2 / n); // M2 >= 0 par construction (Welford)

    std::string cValue[nColumns] = {std::to_string(m_iteration), fmtTime(total), fmtTime(mean), fmtTime(std_dev)};
    for (size_t i = 0; i < nColumns; i++) {
      std::cout << "│";
      size_t size = cValue[i].size();
      printReplicate(0, cWidths[i] - size - 1, " ");
      std::cout << cValue[i] << " ";
    }
    std::cout << "│\n";
  }

  std::string m_name;
  std::size_t m_iteration;
  std::size_t m_level;
  std::vector<Node *> m_daughter;
  Node *m_mother;
  std::chrono::duration<double> m_duration;
  double m_welford_mean; // moyenne courante (algorithme de Welford)
  double m_welford_M2;   // somme des écarts quadratiques (Welford M2)
};

enum TimerTag { CURRENT, ROOT };

template <TimerTag T> Node *&getTimer() {
  static Node *ptr = nullptr;
  return ptr;
}

class ScopedTimer {
public:
  ScopedTimer(Node *node) : m_node(node), m_stopped(false) {
    m_start = std::chrono::steady_clock::now();
  }

  void end() {
    if (m_stopped) return;
    m_stopped      = true;
    double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - m_start).count();
    m_node->m_duration += std::chrono::duration<double>(elapsed);
    // Welford's online algorithm: numerically stable mean and M2
    double n     = double(m_node->m_iteration);
    double delta = elapsed - m_node->m_welford_mean;
    m_node->m_welford_mean += delta / n;
    m_node->m_welford_M2 += delta * (elapsed - m_node->m_welford_mean);
    auto &current = getTimer<CURRENT>();
    current       = current->m_mother;
  }

  ~ScopedTimer() {
    end();
  }

private:
  std::chrono::time_point<std::chrono::steady_clock> m_start;
  Node *m_node;
  bool m_stopped;
};

class OutputManager {
public:
  OutputManager(const char *name = "unnamed") : base_name(name) {}

  template <typename Func, typename... Args> void recursiveCall(Func &func, Node *ptr, Args &...arg) {
    func(ptr, arg...);
    for (auto &it : ptr->m_daughter) recursiveCall(func, it, arg...);
  }

  template <typename Func, typename Sort, typename... Args>
  void recursiveSortedCall(Func &func, Sort mySort, Node *ptr, Args &...arg) {
    func(ptr, arg...);
    std::sort(ptr->m_daughter.begin(), ptr->m_daughter.end(), mySort);
    for (auto &it : ptr->m_daughter) recursiveSortedCall(func, mySort, it, arg...);
  }

  void printTimeTable() {
    auto myPrint   = [](Node *a_ptr, size_t a_shift) { a_ptr->print(a_shift); };
    auto sortComp  = [](Node *a_ptr, Node *b_ptr) { return a_ptr->m_duration.count() > b_ptr->m_duration.count(); };
    auto maxLength = [](Node *a_ptr, size_t &a_count, size_t &) {
      size_t length = a_ptr->m_level * 3 + a_ptr->m_name.size() + 2;
      a_count       = std::max(a_count, length);
    };

    Node *root = getTimer<ROOT>();
    size_t count(0), nbElem(0);
    recursiveCall(maxLength, root, count, nbElem);
    count += 6;
    root->printBanner(count);
    recursiveSortedCall(myPrint, sortComp, root, count);
    root->printEnding(count);
  }

  void writeFile() {
    std::string name = "perflog_" + base_name + ".txt";
    std::ofstream myFile(name, std::ofstream::out);
    Node *root   = getTimer<ROOT>();
    auto myWrite = [](Node *a_ptr, std::ofstream &a_file) {
      std::string space;
      for (std::size_t i = 0; i < a_ptr->m_level; i++) space += "   ";
      double std_dev = std::sqrt(a_ptr->m_welford_M2 / double(a_ptr->m_iteration));
      a_file << space << a_ptr->m_name << " " << a_ptr->m_iteration << " " << a_ptr->m_duration.count() << " "
             << std_dev << "\n";
    };
    recursiveCall(myWrite, root, myFile);
  }

private:
  std::string base_name;
};

} // namespace profiler_detail

#ifdef ENABLE_PROFILING

#define INIT_TIMERS()                                                                                                  \
  profiler_detail::Node *&root_timer_ptr = profiler_detail::getTimer<profiler_detail::ROOT>();                         \
  root_timer_ptr                         = new profiler_detail::Node();                                                \
  profiler_detail::Node *&current        = profiler_detail::getTimer<profiler_detail::CURRENT>();                      \
  current                                = root_timer_ptr;                                                             \
  profiler_detail::ScopedTimer roottimer(current);

#define START_TIMER(name)                                                                                              \
  profiler_detail::Node *&current = profiler_detail::getTimer<profiler_detail::CURRENT>();                             \
  current                         = profiler_detail::getTimer<profiler_detail::CURRENT>()->find(name);                 \
  profiler_detail::ScopedTimer tttt(current);

#define WRITE_TIMERS(...)                                                                                              \
  {                                                                                                                    \
    profiler_detail::OutputManager tmp_out(__VA_ARGS__);                                                               \
    tmp_out.writeFile();                                                                                               \
    tmp_out.printTimeTable();                                                                                          \
  }

#define PRINT_TIMERS(...)                                                                                              \
  roottimer.end();                                                                                                     \
  WRITE_TIMERS(__VA_ARGS__)

#else

#define INIT_TIMERS()
#define START_TIMER(name)
#define WRITE_TIMERS(...)
#define PRINT_TIMERS(...)

#endif

#endif /* PROFILER_SERIE_HPP */
