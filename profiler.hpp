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

// Credit
// This tool has been written by Raphael (CEA)
// and then adapted and added into TOOFUS
#ifndef PROFILER_HPP
#define PROFILER_HPP

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// variables
const size_t nColumns = 4;
const size_t cWidth = 20;
const std::string cName[nColumns] = {"number Of Calls", "min(s)", "mean(s)", "max(s)"};

class Profiler {
public:
  Profiler() : m_daughter() { // only used for root
    m_name = "root";
    m_iteration = 1;
    m_level = 0;
    m_mother = nullptr;
  }

  Profiler(std::string name, Profiler *mother) : m_daughter(), m_duration(0) {
    m_name = name;
    m_iteration = 1;
    m_level = mother->m_level + 1;
    m_mother = mother;
  }

  Profiler *find(std::string name) {
    for (auto it = m_daughter.begin(); it < m_daughter.end(); it++) {
      if ((*it)->m_name == name) {
        (*it)->m_iteration++;
        return (*it);
      }
    }
    Profiler *myTmp = new Profiler(name, this);
    m_daughter.push_back(myTmp);
    return myTmp;
  }

  void printReplicate(size_t begin, size_t end, std::string motif) {
    for (size_t i = begin; i < end; i++)
      std::cout << motif;
  }

  void printBanner(size_t shift) {
    if (m_name == "root") {
#ifdef __MPI
      size_t rank = MPI_Comm_rank(&rank, MPI_COMM_WORLD);
      if (rank == 0) {
        std::cout << " MPI feature activated, rank 0:\n";
#else
      std::cout << " MPI feature is disable for timers, if you use MPI please add -D__MPI\n";
#endif
        std::string start_name = " ┌";
        std::cout << start_name;
        printReplicate(start_name.size(), shift + nColumns * (cWidth + 2) - 1, "─");
        std::cout << "┐\n";
        std::cout << " │    name";
        printReplicate(9, shift, " ");
        for (size_t i = 0; i < nColumns; i++) {
          std::cout << "│";
          int size = cName[i].size();
          printReplicate(0, (int(cWidth) - size - 1), " ");
          std::cout << cName[i] << " ";
        }
        std::cout << "│\n";
        shift += nColumns * (cWidth + 2) - 2;
        std::cout << " │";
        printReplicate(2, shift - 1, "─");
        std::cout << "│\n";
#ifdef __MPI
      }
#endif
    }
  }

  void printEnding(size_t shift) {
    if (m_name == "root") {
#ifdef __MPI
      size_t rank = MPI_Comm_rank(&rank, MPI_COMM_WORLD);
      if (rank == 0) {
#endif
        shift += nColumns * (cWidth + 2);
        std::string end_name = " └";
        std::cout << end_name;
        printReplicate(end_name.size(), shift - 1, "─");
        std::cout << "┘\n";
#ifdef __MPI
      }
#endif
    }
  }

  void print(size_t shift) {
#ifdef __MPI
    size_t rank = MPI_Comm_rank(&rank, MPI_COMM_WORLD);
    if (rank == 0) {
#endif
      size_t realShift = shift;
      std::cout << " │ ";
      size_t currentShift = 3;
      for (int i = 0; i < int(m_level) - 1; i++) {
        int spaceSize = 3;
        for (int j = 0; j < spaceSize; j++)
          std::cout << " ";
        currentShift += spaceSize;
      }
      if (m_level > 0) {
        std::cout << "└──";
        currentShift += 3;
      }
      std::cout << "► " << m_name;
      currentShift += m_name.size() + 1;
      printReplicate(currentShift, realShift, " ");
      std::string cValue[nColumns];

      cValue[0] = std::to_string(m_iteration);
#ifdef __MPI
    }
    double local_min = std::to_string(m_duration.count());
    double local_mean = std::to_string(m_duration.count());
    double local_max = std::to_string(m_duration.count());
    double global_min, global_max, global_mean;

    MPI_Reduce(&global_min, &local_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&global_mean, &local_mean, 1, MPI_DOUBLE, MPI_MEAN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&global_max, &local_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {

      cValue[1] = std::to_string(global_min);
      cValue[2] = std::to_string(global_mean);
      cValue[3] = std::to_string(global_max);
#else
    cValue[1] = std::to_string(m_duration.count());
    cValue[2] = std::to_string(m_duration.count());
    cValue[3] = std::to_string(m_duration.count());
#endif
#ifdef __MPI
      size_t rank = MPI_Comm_rank(&rank, MPI_COMM_WORLD);
#endif
      for (size_t i = 0; i < nColumns; i++) {
        std::cout << "│";
        int size = cValue[i].size();
        printReplicate(0, (int(cWidth) - size - 1), " ");
        std::cout << cValue[i] << " ";
      }
      std::cout << "│\n";
#ifdef __MPI
    }
#endif
  }

  std::string m_name;
  std::size_t m_iteration;
  std::size_t m_level;
  std::vector<Profiler *> m_daughter;
  Profiler *m_mother;
  std::chrono::duration<double> m_duration;
};

enum enumTimer { CURRENT, ROOT };

namespace PROFILER {
template <enumTimer T> Profiler *&getTimer() {
  static Profiler *__current;
  return __current;
}
}; // namespace PROFILER

class chronoTime {
public:
  chronoTime(std::chrono::duration<double> *acc) {
    m_duration = acc;
    m_start = std::chrono::steady_clock::now();
  }
  void end() {
    m_stop = std::chrono::steady_clock::now();
    *m_duration += m_stop - m_start;
  }
  ~chronoTime() {
    end();
    auto &current_timer = PROFILER::getTimer<CURRENT>();
    current_timer = current_timer->m_mother;
  }
  std::chrono::time_point<std::chrono::steady_clock> m_start;
  std::chrono::time_point<std::chrono::steady_clock> m_stop;
  std::chrono::duration<double> *m_duration;
};

// speedup in openmp
// This feature uses outputs of previous simulations
// Do not change Timers in the code, otherwise the results will be wrong

// output in txt
class outputManager {
private:
  std::string base_name;

public:
  outputManager() { base_name = "unnamed"; }

  outputManager(const char *name) { base_name = name; }

  std::string buildName() {
    std::size_t nthreads = 1;

#ifdef _OPENMP
#pragma omp parallel
    { nthreads = omp_get_num_threads(); }
#endif

    std::string file_name = "perflog_" + base_name + "_" + std::to_string(nthreads) + ".txt";
    return file_name;
  }

  template <typename Func, typename... Args> void recursiveCall(Func &func, Profiler *ptr, Args &...arg) {
    func(ptr, arg...);
    for (auto &it : ptr->m_daughter)
      recursiveCall(func, it, arg...);
  }

  template <typename Func, typename Sort, typename... Args>
  void recursiveSortedCall(Func &func, Sort mySort, Profiler *ptr, Args &...arg) {
    func(ptr, arg...);
    std::sort(ptr->m_daughter.begin(), ptr->m_daughter.end(), mySort);
    for (auto &it : ptr->m_daughter)
      recursiveSortedCall(func, mySort, it, arg...);
  }

  void printTimeTable() {
    auto myPrint = [](Profiler *a_ptr, size_t a_shift) { a_ptr->print(a_shift); };

    auto sortComp = [](Profiler *a_ptr, Profiler *b_ptr) {
      return a_ptr->m_duration.count() > b_ptr->m_duration.count();
    };

    auto maxLength = [](Profiler *a_ptr, size_t &a_count, size_t &a_nbElem) {
      size_t length = a_ptr->m_level * 3 + a_ptr->m_name.size();
      a_count = std::max(a_count, length);
      a_nbElem++;
    };
    Profiler *root_timer = PROFILER::getTimer<ROOT>();
    size_t count(0), nbElem(0);

    recursiveCall(maxLength, root_timer, count, nbElem);
    count += 6;
    root_timer->printBanner(count);
    recursiveSortedCall(myPrint, sortComp, root_timer, count);
    root_timer->printEnding(count);
  }

  void writeFile() {
    std::string name = buildName();
    std::ofstream myFile(name, std::ofstream::out);
    Profiler *root_timer = PROFILER::getTimer<ROOT>();
    auto myWrite = [](Profiler *a_ptr, std::ofstream &a_file) {
      std::string space;
      std::string motif = "   ";

      for (std::size_t i = 0; i < a_ptr->m_level; i++)
        space += motif;

      a_file << space << a_ptr->m_name << " " << a_ptr->m_iteration << " " << a_ptr->m_duration.count() << std::endl;
    };

    recursiveCall(myWrite, root_timer, myFile);
  }
};

#ifdef ENABLE_PROFILING
#define INIT_TIMERS()                                                                                                  \
  Profiler *&root_timer_ptr = PROFILER::getTimer<ROOT>();                                                              \
  root_timer_ptr = new Profiler();                                                                                     \
  Profiler *&current = PROFILER::getTimer<CURRENT>();                                                                  \
  current = root_timer_ptr;                                                                                            \
  chronoTime roottimer(&(current->m_duration));

#define START_TIMER(name)                                                                                              \
  Profiler *&current = PROFILER::getTimer<CURRENT>();                                                                  \
  current = PROFILER::getTimer<CURRENT>()->find(name);                                                                 \
  chronoTime tttt(&(current->m_duration));

#define WRITE_TIMERS(...)                                                                                              \
  outputManager tmp_out(__VA_ARGS__);                                                                                  \
  tmp_out.writeFile();                                                                                                 \
  tmp_out.printTimeTable();

#define PRINT_TIMERS(...)                                                                                              \
  roottimer.end();                                                                                                     \
  WRITE_TIMERS(__VA_ARGS__);

#else

#define INIT_TIMERS()
#define START_TIMER(name)
#define WRITE_TIMERS(...)
#define PRINT_TIMERS(...)

#endif

#endif /* end of include guard: PROFILER_HPP */