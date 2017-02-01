# mcdock
McDock: Simple Monte Carlo docking algorithm in C++

## How to install:
Note: you need Open Babel installed (and working) to do the following.

    git clone https://github.com/andersx/mcdock.git

If you have Open Babel installed in a non-standard path, add this to the `OBDIR` variable in the `Makefile`.

Now compile McDock:

    cd mcdock
    make

## Running McDock:
Running McDock for 2000 steps at temperature T = 0.1 (in units of kB T):

    ./mcdock xyz/glucosepane.xyz xyz/water.xyz 2000 0.1

First McDock will minimize the two individual molecules, and then run a constant-temperature (i.e. Metropolis-Hastings) Monte Carlo simulation. The last step is further minimized. Every accepted MC step and the final geometry is dumped in a file named `out.xyz`.

The energy is printed in units of kJ/mol using the MMFF94 force field.

Cheers!
