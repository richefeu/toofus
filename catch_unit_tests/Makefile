# Declaration of variables
CC = clang++
CC_FLAGS = -O3 -std=c++11 -I.. -Wall -Wextra -Wshadow -Wconversion -Wpedantic

# File names
EXEC = tests
SOURCES = main.cpp $(wildcard test_*.cpp)

OBJECTS = $(SOURCES:.cpp=.o)

# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXEC)

# To obtain object files
%.o: %.cpp
	$(CC) -c $(CC_FLAGS) $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)
