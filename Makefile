CXX = g++

OBDIR = /home/andersx/opt/openbabel2.4

INCLUDE = -I$(OBDIR)/include/openbabel-2.0/
LIBS = -L$(OBDIR)/lib

CXX_FLAGS = -std=c++11 -O3 -march=native -Wall
LINKER_FLAGS = -lopenbabel

all: mcdock main conformers

mcdock: src/mcdock.cpp src/utils.hpp
	$(CXX) $(INCLUDE) $(LIBS) $(CXX_FLAGS) src/mcdock.cpp -o mcdock $(LINKER_FLAGS)

main: src/main_ob.cpp
	$(CXX) $(INCLUDE) $(LIBS) $(CXX_FLAGS) src/main_ob.cpp -o main $(LINKER_FLAGS)

conformers: src/conformers.cpp
	$(CXX) $(INCLUDE) $(LIBS) $(CXX_FLAGS) src/conformers.cpp -o conformers $(LINKER_FLAGS)


clean:
	rm -f main
	rm -f mcdock
	rm -f conformers


