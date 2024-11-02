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
// This tool has been initially written by Raphael Prat (CEA)
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
  /**
   * @brief Constructor for Profiler objects at the root level.
   *
   * This constructor is used only for the root Profiler object.
   * It sets the Profiler name to "root", its iteration number to 1,
   * its level to 0 and its mother to nullptr.
   *
   * @param none
   *
   * @return none
   */
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

  /**
   * @brief Find a Profiler object.
   *
   * This function finds a Profiler object in the m_daughter vector
   * by its name. If the object is found, it increments the iteration
   * number of the object and returns a pointer to it. If the object
   * is not found, a new Profiler object is created with the given name
   * and the current object as its mother. The new object is then added
   * to the m_daughter vector and a pointer to it is returned.
   *
   * @param name string giving the name of the Profiler object to find.
   *
   * @return A pointer to the Profiler object.
   */
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

  /**
   * @brief Prints a motif repeatedly within a specified range.
   *
   * This function outputs the given motif to the standard output
   * starting from the 'begin' index up to, but not including, the 'end' index.
   *
   * @param begin The starting index for the repetition.
   * @param end The ending index, non-inclusive, for the repetition.
   * @param motif The string motif to be printed repeatedly.
   */
  void printReplicate(size_t begin, size_t end, std::string motif) {
    for (size_t i = begin; i < end; i++)
      std::cout << motif;
  }

  /**
   * @brief Prints a banner in the console.
   *
   * This function prints a banner in the console which is used to
   * delimit the output of the print() function. The banner is composed
   * of a dashed line, a line with the column names and another dashed
   * line.
   *
   * @param shift The number of spaces to shift the banner to the right.
   *
   * @return none
   */
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
          size_t size = cName[i].size();
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

  /**
   * @brief Prints the ending banner of the print() function.
   *
   * @details The ending banner is composed of a dashed line and a line with
   * the column names.
   *
   * @param shift The number of spaces to shift the banner to the right.
   *
   * @return none
   */
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

  /**
   * @brief Prints the profiling information.
   *
   * @details This function prints a table containing the profiling
   * information of the profiler. The table is printed with the first column
   * (the name of the profiler) shifted to the right by the given amount of
   * spaces. The table is composed of a header line, a dashed line and a line
   * for each profiler with the same level as the current one.
   *
   * @param shift The number of spaces to shift the table to the right.
   *
   * @return none
   */
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
        size_t size = cValue[i].size();
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
  /**
   * @brief Constructs a chronoTime object and starts the timer.
   *
   * This constructor initializes the chronoTime object with a pointer to a
   * duration variable where the elapsed time will be accumulated. It also
   * records the current time as the start time.
   *
   * @param acc Pointer to a std::chrono::duration<double> where the elapsed
   * time will be accumulated.
   */
  chronoTime(std::chrono::duration<double> *acc) {
    m_duration = acc;
    m_start = std::chrono::steady_clock::now();
  }

  /**
   * @brief Stops the timer and updates the duration.
   *
   * This function records the current time as the stop time, computes the
   * elapsed time since the last start, and adds it to the accumulated duration.
   */
  void end() {
    m_stop = std::chrono::steady_clock::now();
    *m_duration += m_stop - m_start;
  }

  /**
   * @brief Destructor for chronoTime that stops the timer and goes up one level.
   *
   * This destructor calls the end() function to stop the timer and update the
   * duration. It also sets the current timer to the current timer's mother,
   * effectively going up one level in the profiler's tree structure.
   */
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

  /**
   * @brief Builds a file name for the profiler's output.
   *
   * This function builds a file name based on the base name and the number of
   * threads used in the simulation. The file name is of the form
   * "perflog_<base_name>_<nthreads>.txt".
   *
   * @return The constructed file name.
   */
  std::string buildName() {
    std::size_t nthreads = 1;

#ifdef _OPENMP
#pragma omp parallel
    {
      nthreads = omp_get_num_threads();
    }
#endif

    std::string file_name = "perflog_" + base_name + "_" + std::to_string(nthreads) + ".txt";
    return file_name;
  }

  /**
   * @brief Recursively calls a function on all nodes of a tree.
   *
   * This function takes a function `func` and its arguments `arg...` and calls
   * `func` on the current node and all of its daughter nodes recursively. This
   * is a convenient way to traverse the profiler's tree structure.
   */
  template <typename Func, typename... Args> void recursiveCall(Func &func, Profiler *ptr, Args &...arg) {
    func(ptr, arg...);
    for (auto &it : ptr->m_daughter)
      recursiveCall(func, it, arg...);
  }

  /**
   * @brief Recursively calls a function on all nodes of a tree, sorting children nodes.
   *
   * This function takes a function `func`, a sorting criterion `mySort`, and its arguments `arg...`.
   * It calls `func` on the current node, sorts its daughter nodes using `mySort`, and recursively
   * applies the same process to each daughter node. This allows for a depth-first traversal of the
   * profiler's tree structure, with nodes processed in a sorted order.
   *
   * @param func A callable object that operates on a Profiler node.
   * @param mySort A sorting criterion used to order the daughter nodes.
   * @param ptr A pointer to the current Profiler node.
   * @param arg Additional arguments passed to `func`.
   */
  template <typename Func, typename Sort, typename... Args>
  void recursiveSortedCall(Func &func, Sort mySort, Profiler *ptr, Args &...arg) {
    func(ptr, arg...);
    std::sort(ptr->m_daughter.begin(), ptr->m_daughter.end(), mySort);
    for (auto &it : ptr->m_daughter)
      recursiveSortedCall(func, mySort, it, arg...);
  }

  /**
   * @brief Prints a sorted table of profiling results.
   *
   * This function prints a table of profiling results, sorted in descending order of time.
   * The table is indented to show the hierarchy of timers, and shows the name of each timer,
   * its duration, and the number of calls.
   */
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

  /**
   * @brief Writes profiling results to a file.
   *
   * This function writes the profiling results to a file, with each line containing
   * the name of the timer, the number of iterations, and the total duration. The
   * file is named according to the build name, and is written in the current working
   * directory. The file is overwritten if it already exists.
   */
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