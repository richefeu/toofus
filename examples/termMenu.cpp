#include "../termMenu.hpp"


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
    return 0;
}