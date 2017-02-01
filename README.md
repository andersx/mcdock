# mcdock
McDock: Simple Monte Carlo docking algorithm in C++

## How to install:
Note: you need Open Babel installed (and working) to do the following.

    git clone https://github.com/andersx/mcdock.git

If you have Open Babel installed in a non-standard path, add this to the `OBDIR` variable in the `Makefile`.

Now compile McDock:

    cd mcdock
    make

## How it works:

1. The host molecule is minimized
2. A rotamer search is performed for the ligand
   * Each rotamer is energy minimized
   * A number of Monte Carlo Metropolis-Hastings (i.e. constant temperature) docking simulations is performed for each rotamer - the number of trajectories and steps for each trajectory is given in the input line, as well as the MC temperature.
   * The final docking conformation in minimized
   * The binding energy is calculated

The strongest binding motif is saved in the file `min.xyz`. 

The binding energy is calculated as:

>  Ebind = E - [Ehost + Eligand\_min]

Where `E` is the energy of a specific binding conformation, `Ehost` is the energy of the minimized host molecule, and `Eligand_min` is the energy of the ligand in the lowest-energy conformation from the minimzed rotamer search.

The MC move is combination of random translations and rotations.

## Running McDock:
Running McDock for 2000 steps, with 10 trajectories at temperature T = 0.1 (in units of kB T):

    ./mcdock xyz/glucosepane.xyz xyz/water.xyz 2000 0.1 10

The energy is printed in units of **kcal/mol** using the MMFF94 force field. The MC temperature unit is **kB T**.

## Future developments:

A MOPAC interface is being implemented to screen at the semiempirical level.


Cheers!
