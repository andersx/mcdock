CXX = g++

OBDIR = /home/andersx/opt/openbabel2.4

INCLUDE = -I$(OBDIR)/include/openbabel-2.0/
LIBS = -L$(OBDIR)/lib

CXX_FLAGS = -std=c++11 -O3 -march=native -Wall
LINKER_FLAGS = -lopenbabel

all: mcdock mcdock_pm6 trajectory trajectory_pm6 conformers

mcdock: src/mcdock.cpp src/utils.hpp
	$(CXX) $(INCLUDE) $(LIBS) $(CXX_FLAGS) src/mcdock.cpp -o mcdock $(LINKER_FLAGS)

mcdock_pm6: src/mcdock_pm6.cpp src/utils.hpp
	$(CXX) $(INCLUDE) $(LIBS) $(CXX_FLAGS) src/mcdock_pm6.cpp -o mcdock_pm6 $(LINKER_FLAGS)

trajectory: src/trajectory.cpp
	$(CXX) $(INCLUDE) $(LIBS) $(CXX_FLAGS) src/trajectory.cpp -o trajectory $(LINKER_FLAGS)

trajectory_pm6: src/trajectory_pm6.cpp
	$(CXX) $(INCLUDE) $(LIBS) $(CXX_FLAGS) src/trajectory_pm6.cpp -o trajectory_pm6 $(LINKER_FLAGS)

conformers: src/conformers.cpp
	$(CXX) $(INCLUDE) $(LIBS) $(CXX_FLAGS) src/conformers.cpp -o conformers $(LINKER_FLAGS)


clean:
	rm -f mcdock
	rm -f mcdock_pm6
	rm -f trajectory 
	rm -f trajectory_pm6
	rm -f conformers


