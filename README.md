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


The binding energy is calculated as:

>  Ebind = E - [Ehost + Eligand\_min]

Where `E` is the energy of a specific binding conformation, `Ehost` is the energy of the minimized host molecule, and `Eligand_min` is the energy of the ligand in the lowest-energy conformation from the minimzed rotamer search.

The MC move is combination of random translations and rotations.

* The conformers found from the rotamer search is saved in the file `conformers.xyz`. 
* The last motif from each trajectory is saved in the file `out.xyz`. 
* The strongest binding motif is saved in the file `min.xyz`. 


## Program options:

The possible program options are:

    Running McDock 0.2 alpha
    Usage: ./mcdock --target file1.xyz --ligand file2.xyz [--args]
    
    Optional arguments:
    --energy [string]       Potential energy function "MMFF94" (default),
                            "UFF", "PM6-D3H+" (requires MOPAC). 
    --temperature [float]   Temperature in units of [R T] (default = 1.0).
    --trajectories [int]    Number of independent trajectories (default = 10).
    --steps [int]           Number of Monte Carlo steps in each trajectories (default = 1000).
    --no-rotor-search       Disable rotor search.  (default = perform rotor search).

The possible potential energy functions are the `MMFF` and `UFF` force fields as implemented in Open Babel.
If MOPAC is installed in `/opt/mopac/` it is possible to use the semi-empirical method `PM6-D3H+` as well.

## Running McDock:

Running McDock for 2000 steps, with 5 trajectories at temperature T = 0.3 (in units of R T), using the MMFF94 force field:

    ./mcdock --target xyz/glucosepane.xyz --ligand xyz/small_conformer.xyz \
        --trajectories 5 --temperature 0.3 --mc-steps 2000 --energy MMFF94

The energy is printed in units of **kcal/mol**. The MC temperature unit is **R T** where **R** is the ideal gas constant.

This will generate something like the following output:

    Running McDock 0.2 alpha
    Target (minimized) E =   -33.4960 kcal/mol     file: xyz/glucosepane.xyz 
    Performing rotor search for ligand molecule    file: xyz/small_conformer.xyz
    ..tot conformations = 3
    ..tot confs tested = 3
    ..below energy threshold = 3
    Found   1 rotatable bonds
    Rotamer    0     E =    20.2736 kcal/mol
    Rotamer    1     E =    19.9836 kcal/mol
    Rotamer    2     E =    20.2721 kcal/mol
    Rotamer    3     E =    24.2927 kcal/mol
    Lowest energy conformation  E =    19.9836 kcal/mol
    Running   5 trajectories for       2000 steps.
    MC temperature (tau) =     0.3000
    
    Conformation:       Trajectory:         Acceptance rate:    Final Ebind:
    ---------------------------------------------------------------------------
       1 /   4            1 /   5             7.40 %           -9.0204 kcal/mol <---- New lowest
       1 /   4            2 /   5             9.50 %          -11.8772 kcal/mol <---- New lowest
       1 /   4            3 /   5            19.74 %           -2.4114 kcal/mol
       1 /   4            4 /   5             9.15 %          -10.1515 kcal/mol
       1 /   4            5 /   5             6.25 %          -10.1287 kcal/mol
       2 /   4            1 /   5             8.80 %           -8.3443 kcal/mol
       2 /   4            2 /   5            12.99 %          -10.4063 kcal/mol
       2 /   4            3 /   5             9.40 %          -10.6525 kcal/mol
       2 /   4            4 /   5            15.54 %           -8.1072 kcal/mol
       2 /   4            5 /   5             7.95 %           -8.9660 kcal/mol
       3 /   4            1 /   5             5.85 %          -11.0246 kcal/mol
       3 /   4            2 /   5            12.14 %          -11.3170 kcal/mol
       3 /   4            3 /   5             6.95 %          -10.3502 kcal/mol
       3 /   4            4 /   5            11.84 %           -8.1763 kcal/mol
       3 /   4            5 /   5             8.10 %           -6.1784 kcal/mol
       4 /   4            1 /   5             8.40 %           -4.8038 kcal/mol
       4 /   4            2 /   5             7.55 %           -3.6463 kcal/mol
       4 /   4            3 /   5             7.20 %           -3.3909 kcal/mol
       4 /   4            4 /   5            24.74 %           -1.4923 kcal/mol
       4 /   4            5 /   5             9.90 %           -1.9934 kcal/mol
    
    Optimized E_bind =   -11.8772 kcal/mol    Elapsed time =   8.47 seconds

A "good" temperature will usually yield an acceptance rate around 20 %.


## License:

McDock is licensed under the MIT open source license.
