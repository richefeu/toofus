#include "../consoleProgressBar.hpp"

int main(int argc, char const *argv[]) {
  ConsoleProgressBar cpb(100);
  cpb.setTitle("Title: ");
  cpb.setWidth(50);
  cpb.setProgressChar(124);
  cpb.setVoidChar('-');
  cpb.update(10, std::cout);
  return 0;
}
