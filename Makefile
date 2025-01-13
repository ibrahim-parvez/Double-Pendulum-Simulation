# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -O2 -std=c++11
LIBS = -lcpgplot -lpgplot -lX11 -lm
# Executables
EXEC_PENDULUM = pendulum
EXEC_POSITION = position
EXEC_GRAPHS = graphs
EXEC_DEMO = demo
EXEC_ENERGY = energy
# Source files
SRC_PENDULUM = pendulum.cpp
SRC_POSITION = position.cpp
SRC_GRAPHS = graphs.cpp
SRC_DEMO = demo.cpp
SRC_ENERGY = energy.cpp
# Default target: build and run pendulum, and build position
all: $(EXEC_PENDULUM) $(EXEC_POSITION) $(EXEC_GRAPHS) $(EXEC_DEMO) $(EXEC_ENERGY)
./$(EXEC_PENDULUM)
./$(EXEC_POSITION)
./$(EXEC_GRAPHS)
./$(EXEC_DEMO)
./$(EXEC_ENERGY)
# Compile the pendulum executable
$(EXEC_PENDULUM): $(SRC_PENDULUM)
$(CXX) $(CXXFLAGS) $(SRC_PENDULUM) -o $(EXEC_PENDULUM) $(LIBS)
# Compile the position executable
$(EXEC_POSITION): $(SRC_POSITION)
$(CXX) $(CXXFLAGS) $(SRC_POSITION) -o $(EXEC_POSITION) $(LIBS)
# Compile the graphs executable
$(EXEC_GRAPHS): $(SRC_GRAPHS)
$(CXX) $(CXXFLAGS) $(SRC_GRAPHS) -o $(EXEC_GRAPHS) $(LIBS)
# Compile the graphs executable
$(EXEC_DEMO): $(SRC_DEMO)
$(CXX) $(CXXFLAGS) $(SRC_DEMO) -o $(EXEC_DEMO) $(LIBS)
# Compile the energy graph executable
$(EXEC_ENERGY): $(SRC_ENERGY)
$(CXX) $(CXXFLAGS) $(SRC_ENERGY) -o $(EXEC_ENERGY) $(LIBS)
# Clean up the build files
clean:
rm -f $(EXEC_PENDULUM) $(EXEC_POSITION) $(EXEC_GRAPHS) $(EXEC_DEMO)
$(EXEC_ENERGY)
