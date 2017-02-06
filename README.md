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
Running McDock for 2000 steps, with 10 trajectories at temperature T = 1.0 (in units of kB T):

    ./mcdock xyz/glucosepane.xyz xyz/water.xyz 2000 1.0 10

The energy is printed in units of **kcal/mol** using the MMFF94 force field. The MC temperature unit is **kB T**.

This will generate something like the following output:

    Running McDock 0.2 alpha
    Pocket (minimized) E =   -33.4960 kcal/mol    file = xyz/glucosepane.xyz
    Performing rotor search for ligand molecule
    Found   0 rotatable bonds
    Rotamer    0     E =     0.0000 kcal/mol
    Lowest energy conformation  E =     0.0000 kcal/mol
    Rotamer:   1 /   1   Trajectory:   1 /  10    acceptance =  17.44 %   E_bind =    -9.3618 kcal/mol      <---- New best
    Rotamer:   1 /   1   Trajectory:   2 /  10    acceptance =  17.69 %   E_bind =    -6.8568 kcal/mol
    Rotamer:   1 /   1   Trajectory:   3 /  10    acceptance =  28.39 %   E_bind =    -4.7040 kcal/mol
    Rotamer:   1 /   1   Trajectory:   4 /  10    acceptance =  16.34 %   E_bind =    -8.6960 kcal/mol
    Rotamer:   1 /   1   Trajectory:   5 /  10    acceptance =  15.94 %   E_bind =    -9.0156 kcal/mol
    Rotamer:   1 /   1   Trajectory:   6 /  10    acceptance =  20.94 %   E_bind =    -9.8126 kcal/mol      <---- New best
    Rotamer:   1 /   1   Trajectory:   7 /  10    acceptance =  19.64 %   E_bind =    -9.3822 kcal/mol
    Rotamer:   1 /   1   Trajectory:   8 /  10    acceptance =  25.84 %   E_bind =    -5.4531 kcal/mol
    Rotamer:   1 /   1   Trajectory:   9 /  10    acceptance =  15.09 %   E_bind =    -9.3082 kcal/mol
    Rotamer:   1 /   1   Trajectory:  10 /  10    acceptance =  18.14 %   E_bind =    -8.7319 kcal/mol
    Optimized E_bind =    -9.8126 kcal/mol    Elapsed time =   3.01 seconds

A "good" temperature will usually yield an acceptance rate around 20 %.
If you wish to use PM6 instead of MMFF94 through MOPAC (requires MOPAC installed), you can use the `mcdock_pm6` executable instead.

## Other McDock tools:

You can run a single MC trajectory with the `trajectory` or `trajectory_pm6` tools.

    ./trajectory xyz/glucosepane.xyz xyz/water.xyz 2000 1.0

Another tool in the `conformers` which runs a conformational search through Open babel and produces a file named `conformers.xyz`

    ./conformers xyz/conformer_test.xyz

