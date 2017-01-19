CXX = g++

OBDIR = /home/andersx/opt/openbabel2.4

INCLUDE = -I$(OBDIR)/include/openbabel-2.0/
LIBS = -L$(OBDIR)/lib

all: main

main: src/main_ob.cpp
	$(CXX) $(INCLUDE) $(LIBS) -std=c++11 -O3 -march=native src/main_ob.cpp -o main -lopenbabel

clean:
	rm -f main


