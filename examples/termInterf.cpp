#include "../termInterf.hpp"

int main() {
    TermInterf app;

    // Menus
    TermMenu mainMenu("select an action", {"select fruit","select color", "quit"});
    TermMenu fruitMenu("Select a fruit", {"Apple", "Cherry", "Orange"});
    TermMenu colorMenu("Pick a color", {"Red", "Green", "Blue"});

    app.addMenu("main", mainMenu);
    app.addMenu("fruit", fruitMenu);
    app.addMenu("color", colorMenu);

    // Actions
    app.addAction("apple_action", []() -> std::string {
        Prompt<int> p("How many?");
        int number = p.run();
        std::cout << "You want: " << number << std::endl;
        return "main"; // return to previous menu
    });

    app.addAction("exit", []() -> std::string {
        std::cout << "Goodbye!\n";
        exit(0);
        return ""; // stops program (no previous menu)
    });

    // Transitions
    app.link("main", 0, "fruit");
    app.link("main", 1, "color");
    app.link("main", 2, "exit");

    app.link("fruit", 0, "apple_action");
    app.link("fruit", 1, "apple_action");
    app.link("fruit", 2, "apple_action");

    app.link("color", 0, "main");
    app.link("color", 1, "main");
    app.link("color", 2, "main");

    app.setStart("main");
    app.run();

    return 0;
}
