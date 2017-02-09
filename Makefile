CXX = g++

OBDIR = /home/andersx/opt/openbabel2.4

INCLUDE = -I$(OBDIR)/include/openbabel-2.0/
LIBS = -L$(OBDIR)/lib

CXX_FLAGS = -std=c++11 -O3 -march=native -Wall
LINKER_FLAGS = -lopenbabel

all: mcdock

mcdock: src/mcdock.cpp src/utils.hpp
	$(CXX) $(INCLUDE) $(LIBS) $(CXX_FLAGS) src/mcdock.cpp -o mcdock $(LINKER_FLAGS)

clean:
	rm -f mcdock
	rm -f temp.mop
	rm -f temp.out
	rm -f temp.arc
	rm -f min.xyz
	rm -f out.xyz
	rm -f conformers.xyz


