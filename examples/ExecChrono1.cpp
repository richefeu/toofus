#include "../ExecChrono.hpp"
#include <iostream>
#include <vector>

void bubbleSort(std::vector<double> &arr) {
  int i, j;
  for (i = 0; i < (int)arr.size() - 1; i++) {
    // Last i elements are already in place
    for (j = 0; j < (int)arr.size() - i - 1; j++) {
      if (arr[j] > arr[j + 1]) {
        double tmp = arr[j];
        arr[j] = arr[j + 1];
        arr[j + 1] = tmp;
      }
    }
  }
}

int main(int argc, char const *argv[]) {

  int n = 25000;
  std::vector<double> tab(n);
  for (int i = 0; i < n; i++) {
    tab[i] = 10.0 * (double)rand() / (double)RAND_MAX;
  }

  for (int i = 0; i < 10; i++) {
    std::cout << tab[i] << ' ';
  }
  std::cout << '\n';

  ExecChrono MyChrono; // the chrono will start automatically (but we can use 'start' if we want)
  bubbleSort(tab);
  MyChrono.stop();

  for (int i = 0; i < 10; i++) {
    std::cout << tab[i] << ' ';
  }
  std::cout << '\n';

  return 0;
}