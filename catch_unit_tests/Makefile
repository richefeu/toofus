# Declaration of variables
CC = g++-12
CC_FLAGS = -O3 -std=c++11 -I.. -Wall -Wextra -Wshadow -Wconversion -pedantic


# File names
EXEC = tests
SOURCES = main.cpp $(wildcard test_*.cpp)

OBJECTS = $(SOURCES:.cpp=.o)

# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(CC_FLAGS) $(OBJECTS) -o $(EXEC)

# To obtain object files
%.o: %.cpp
	$(CC) -c $(CC_FLAGS) $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)
