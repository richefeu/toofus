#include "../profiler.hpp"

#include <cstdlib>
#include <iostream>
#include <vector>

int main(int argc, char const *argv[]) {
  INIT_TIMERS();

  {
    START_TIMER("MAIN");

    std::vector<int> v;

    for (size_t i = 0; i < 1000000; i++) {
      {
        START_TIMER("FIRST");
        v.push_back(rand());
      }
    }

    for (size_t i = 1; i < 1000000; i++) {
      {
        START_TIMER("SECOND");
        v[i] -= v[i - 1];
      }
    }
  }

  PRINT_TIMERS_NAME("coucou");

  return 0;
}