#define ENABLE_PROFILING
#include "../profiler_serie.hpp"
#include <cmath>

void heavy_work(int n) {
  START_TIMER("heavy_work");
  volatile double x = 0.0;
  for (int i = 0; i < n; i++) x += std::sin(i) * std::cos(i);
}

void light_work(int n) {
  START_TIMER("light_work_bilou_bilou long_name");
  volatile double x = 0.0;
  for (int i = 0; i < n; i++) x += i * 0.001;
}

void step(int iter) {
  START_TIMER("step");
  heavy_work(iter % 2 == 0 ? 500000 : 200000);
  light_work(100000);
}

int main() {
  INIT_TIMERS();

  {
    START_TIMER("init");
    volatile double x = 0.0;
    for (int i = 0; i < 100000; i++) x += i;
  }

  for (int i = 0; i < 10; i++) { step(i); }

  PRINT_TIMERS("example");
  return 0;
}
