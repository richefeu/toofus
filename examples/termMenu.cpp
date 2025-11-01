#include "../termMenu.hpp"
#include "../termPrompt.hpp"
#include <functional>
#include <vector>

int main() {

  TermMenu menu;
  menu.setPrompt("Select a color");
  menu.addChoice("Blue");
  menu.addChoice("Red");
  menu.addChoice("Green");
  menu.addChoice("Orange");
  menu.addChoice("None");

  int selection = menu.run();
  std::cout << "You selected: " << selection << std::endl;

  TermMenu menuF("Select a fruit", {"Apple", "Cherry", "Orange"});
  selection = menuF.run();
  std::cout << "You selected: " << selection << std::endl;

  TermMenu menuAction("Select an action", {"print coucou", "print toto", "print tata"});
  std::vector<std::function<void()>> actions = {[]() { std::cout << "coucou" << std::endl; },
                                                []() { std::cout << "toto" << std::endl; },
                                                []() { std::cout << "tata" << std::endl; }};
  selection                                  = menuAction.run();
  if (selection >= 0 && selection < actions.size()) { actions[selection](); }

  Prompt<int> myPrompt("Enter a number");
  int N = myPrompt.run();
  std::cout << "You entered: " << N << std::endl;
  
  Prompt<size_t> myPrompt2("Enter a positive number");
  size_t N2 = myPrompt2.run();
  std::cout << "You entered: " << N2 << std::endl;
  
  Prompt<double> myPrompt3("Enter a float number");
  double D = myPrompt3.run();
  std::cout << "You entered: " << D << std::endl;
  
  Prompt<std::string> myPrompt4("Enter a string");
  std::string S = myPrompt4.run();
  std::cout << "You entered: '" << S << "'" <<  std::endl;

  return 0;
}