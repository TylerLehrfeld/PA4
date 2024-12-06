# Compiler
CXX = g++

#paths for dependencies
PATHS = ./dependencies/eigen-3.4.0

# Flags
CXXFLAGS = -g -I $(PATHS)


# Program files
SOURCES = Matrix.cpp Transform.cpp 
HEADERS = data_reader.h helperFunctions.h Matrix.h PointCloudTransform.h transformTransform.h main-helpers.h triangle-functions.h
MAIN_OBJECTS = main.o $(SOURCES:.cpp=.o)
TEST_OBJECTS = test.o $(SOURCES:.cpp=.o)
PERFORMANCE_OBJECTS = performance.o $(SOURCES:.cpp=.o)

# Main compilation
main: $(MAIN_OBJECTS)
	$(CXX) $(CXXFLAGS) $(MAIN_OBJECTS) -o main

performance: $(PERFORMANCE_OBJECTS)
	$(CXX) $(CXXFLAGS) $(PERFORMANCE_OBJECTS) -o performance

# Test compilation
test: $(TEST_OBJECTS)
	$(CXX) $(CXXFLAGS) $(TEST_OBJECTS) -o test

# Rule for compiling .cpp files into .o files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up object files and executables
clean:
	rm -f $(MAIN_OBJECTS) $(TEST_OBJECTS) main test Matrix-test.o

clean-nonexecutables:
	rm -f $(MAIN_OBJECTS) $(TEST_OBJECTS)  ${PERFORMANCE_OBJECTS}
