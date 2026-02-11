// Example usage of CmdLine-single.h
// This demonstrates that the single-header version works correctly

#include "CmdLine.h"
#include <iostream>

int main(int argc, char* argv[]) {
    try {
        // Create a command line parser
        TCLAP::CmdLine cmd("Example program using CmdLine-single.h", ' ', "1.0");
        
        // Add some arguments
        TCLAP::ValueArg<std::string> nameArg("n", "name", "Your name", true, "default", "string");
        TCLAP::SwitchArg verboseArg("v", "verbose", "Verbose output", false);
        TCLAP::MultiArg<int> numbersArg("i", "numbers", "Numbers to process", false, "int");
        
        // Add arguments to the command line
        cmd.add(nameArg);
        cmd.add(verboseArg);
        cmd.add(numbersArg);
        
        // Parse the command line
        cmd.parse(argc, argv);
        
        // Get the values
        std::string name = nameArg.getValue();
        bool verbose = verboseArg.getValue();
        std::vector<int> numbers = numbersArg.getValue();
        
        // Output the results
        std::cout << "Name: " << name << std::endl;
        std::cout << "Verbose: " << (verbose ? "true" : "false") << std::endl;
        std::cout << "Numbers: ";
        for (size_t i = 0; i < numbers.size(); i++) {
            std::cout << numbers[i] << " ";
        }
        std::cout << std::endl;
        
    } catch (TCLAP::ArgException &e) {
        std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
        return 1;
    }
    
    return 0;
}