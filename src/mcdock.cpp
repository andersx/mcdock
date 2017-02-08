// MIT License
// 
// Copyright (c) 2017 Anders Steen Christensen
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Conformer search requires HAVE_EIGEN defined
// due to bug/quirk in Open Babel
#define HAVE_EIGEN
#include <stdio.h>

#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <random>
#include <getopt.h>

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/obutil.h>

#include "utils.hpp"

// Prints help
void print_help() {

        printf("Usage: ./mcdock --target file1.xyz --ligand file2.xyz [--args]\n\n");

        printf("Optional arguments:\n");

        printf("--energy [string]       Potential energy function \"MMFF94\" (default),\n");
        printf("                        \"UFF\", \"PM6-D3H+\" (requires MOPAC). \n");
        printf("--temperature [float]   Temperature in units of [R T] (default = 1.0).\n");
        printf("--trajectories [int]    Number of independent trajectories (default = 10).\n");
        printf("--steps [int]           Number of Monte Carlo steps in each trajectories (default = 1000).\n");
        printf("--no-rotor-search       Disable rotor search.  (default = perform rotor search).\n");

        printf("\n");
}


// Container for command-line options
struct Option {

    std::string ff = "MMFF94";
    std::string target;
    std::string ligand;
    bool use_mopac = false;
    unsigned int trajectories = 1;
    int steps = 1000;
    double temperature = 1.0;
    double verbose_flag = false;
    double rotor_flag = false;
    double help_flag= false;

};

// Option parser
Option get_options (int argc, char **argv) {

    if (argc < 2) {
        print_help();
        exit(0);
    }

    int c;

    Option opts;

    int verbose_flag = 0;
    int rotor_flag = 0;
    int help_flag = 0;

    while (1) {
        static struct option long_options[] = {
            {"verbose",           no_argument,    &verbose_flag,   1},
            {"help",              no_argument,    &help_flag,      1},
            {"no-rotor-search",   no_argument,    &rotor_flag,     1},
            {"target",        required_argument, 0, 'a'},
            {"ligand",        required_argument, 0, 'b'},
            {"trajectories",  required_argument, 0, 'c'},
            {"energy",        required_argument, 0, 'd'},
            {"temperature",   required_argument, 0, 'e'},
            {"mc-steps",      required_argument, 0, 'f'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "a:b:c:d:f:", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

            case 0:
                if (long_options[option_index].flag != 0) break;
            break;

            case 'a':
              opts.target = optarg;
              break;

            case 'b':
              opts.ligand = optarg;
              break;

            case 'e':
              opts.temperature = std::stod(optarg);
              break;

            case 'c':
              opts.trajectories = atoi(optarg);
              break;

            case 'f':
              opts.steps = atoi(optarg);
              break;

            case 'd': {
                std::string ff = optarg;

                if (ff.compare("MMFF94") == 0) {
                    opts.ff = optarg;
                    opts.use_mopac = false;
                } else if (ff.compare("MMFF94s") == 0) {
                    opts.ff = optarg;
                    opts.use_mopac = false;
                } else if (ff.compare("GAFF") == 0) {
                    opts.ff = optarg;
                    opts.use_mopac = false;
                } else if (ff.compare("UFF") == 0) {
                    opts.ff = optarg;
                    opts.use_mopac = false;
                } else if (ff.compare("PM6-D3H4") == 0) {
                    opts.ff = optarg;
                    opts.use_mopac = true;

                }  else {
                    printf ("ERROR: Unsupported force field `%s'\n", optarg);
                    exit(0);
                }
                break;
            }
            case '?':
              /* getopt_long already printed an error message. */
              break;

            default:
              abort ();
            }
    }

  if (verbose_flag) opts.verbose_flag = true;
  if (rotor_flag) opts.rotor_flag = true;
  if (help_flag) {
        print_help();
        exit(0);
    }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
    {
      printf ("non-option ARGV-elements: ");
      while (optind < argc)
        printf ("%s ", argv[optind++]);
      putchar ('\n');
    }

  if (opts.temperature < 0) { 
      std::cout << "ERROR: Positive temperature required" << std::endl;
      exit(0);
  }

    if (opts.ligand.compare("") == 0) {
        printf("ERROR: No ligand xyz-file specified.\n");

        print_help();
        exit(0);
    }

    if (opts.target.compare("") == 0) {
        printf("ERROR: No target xyz-file specified.\n");

        print_help();
        exit(0);
    }

    return opts;

}



int main(int argc, char *argv[]) {

    // Start timer
    OpenBabel::OBStopwatch timer;
    timer.Start();

    printf("Running McDock 0.2 alpha\n");
    Option opts = get_options(argc, argv);

    // Get OBMols
    OpenBabel::OBMol mol = readfile(opts.target);
    OpenBabel::OBMol ligand = readfile(opts.ligand);
    OpenBabel::vector3 com; 

    // Center ligand molecule
    com = get_com(ligand);
    move_molecule(ligand, -com);

    // Gent energies of target molecule
    double ea;
    if (opts.use_mopac == false) {
        ea = minimize_molecule(mol, ff);
    } else {
        ea = mopac_optimize(mol);
    }

    printf("Target (minimized) E = %10.4f kcal/mol     file: %s \n", ea, opts.target.c_str());
    
    // Center target molecule
    com = get_com(mol);
    move_molecule(mol, -com);

    // Do rotor search for ligand molecule
    if (opts.rotor_flag == 0) {
        printf("Performing rotor search for ligand molecule    file: %s\n", opts.ligand.c_str());
        set_conformations(ligand);
        printf("Found %3i rotatable bonds\n", ligand.NumRotors());
    } else {
        printf("No rotors search performed for ligand molecule file: %s\n", opts.ligand.c_str());
    }

    // Initialize random number generators
    std::default_random_engine generator;
    std::uniform_real_distribution<double> random_length(0.0, 0.25);
    std::uniform_real_distribution<double> random_angle(0.0, 90.0/180.0*M_PI);
    std::uniform_real_distribution<double> uniform1(0.0, 1.0);

    // Initialize variables
    double e_low = std::numeric_limits<double>::infinity();
    double eb_min = std::numeric_limits<double>::infinity();

    // Initialize force field
    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);

    // Initialize variables
    double theta;
    OpenBabel::vector3 temp;
    OpenBabel::vector3 dir;
    OpenBabel::vector3 rot;
    OpenBabel::vector3 move;

    double eb;

    // Initialize file output 
    OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("xyz", "xyz");
    std::remove("conformers.xyz");
    std::ofstream ofs_conf("conformers.xyz");

    // Minimize all conformers from rotamer search and dump to file
    for (int c = 0; c < ligand.NumConformers(); ++c) {

        ligand.SetConformer(c);

        if (opts.use_mopac == false) {
            eb = minimize_molecule(ligand, ff);
        } else {
            eb = mopac_optimize(ligand);
        }

        if (eb < eb_min) eb_min = eb;

        printf("Rotamer %4i     E = %10.4f kcal/mol\n", c, eb);
        conv.Write(&ligand, &ofs_conf);
        }

    ofs_conf.close();
    printf("Lowest energy conformation  E = %10.4f kcal/mol\n", eb_min);
    printf("Running %3i trajectories for %10i steps.\n", opts.trajectories, opts.steps);
    printf("MC temperature (tau) = %10.4f\n", opts.temperature);

    // Initialize output for trajectory
    std::remove("out.xyz");
    std::ofstream ofs("out.xyz");

    // Print header
    printf("\nConformation:       Trajectory:         Acceptance rate:    Final Ebind:\n");
    printf("---------------------------------------------------------------------------\n");
    
    for (int c = 0; c < ligand.NumConformers(); ++c) {

        ligand.SetConformer(c);

        // Minimize conformer
        if (opts.use_mopac == false) {
            eb = minimize_molecule(ligand, ff);
        } else {
            eb = mopac_optimize(ligand);
        }

        for (unsigned int n = 0; n < opts.trajectories; n++) {

            // Translate molecule randomly
            dir.randomUnitVector();
            com = get_com(ligand);
            temp = dir * 4.0 - com;
            move_molecule(ligand, temp);

            // Rotate molecule randomly
            rot.randomUnitVector();
            theta = random_angle(generator); 
            rotate_molecule(ligand, rot, theta);

            // Make object to roll-back rejected MC moves 
            OpenBabel::OBMol mol_ligand = mol;
            mol_ligand += ligand;
            OpenBabel::OBMol mol_old = mol_ligand;
            mol_old.SetCoordinates(mol_ligand.GetCoordinates());

            // Initialize MC energy
            pFF->Setup(mol_ligand);
            double e;

            // Calculate energy of complex
            if (opts.use_mopac == false) {
                e = pFF->Energy();

            } else {
                e = mopac_energy(mol_ligand);
            }

            // Initialize variables
            double energy_old = e;
            double delta_e = 0.0;
            int accept = 0;
            double acceptance_ratio;

            // Start and end ID of ligand in complex
            int startid = mol_ligand.NumAtoms() - ligand.NumAtoms() + 1;
            int endid = mol_ligand.NumAtoms() + 1;

            // Begin MC simulation
            for (int step = 0; step < opts.steps; step++) {

                // Translation move
                move.randomUnitVector();
                move *= random_length(generator);
                move_molecule(mol_ligand, move, startid=startid, endid=endid);

                // Rotation move
                rot.randomUnitVector();
                theta = random_angle(generator); 
                rotate_molecule(mol_ligand, rot, theta, startid=startid, endid=endid);

                // Evaluate total energy
                pFF->SetCoordinates(mol_ligand);
                if (opts.use_mopac == false) {
                    e = pFF->Energy();
                } else {
                    e = mopac_energy(mol_ligand);
                }

                delta_e =  e - energy_old;

                // Metropolis-Hastings MC criterion, accept ...
                if (std::exp( - delta_e / opts.temperature) >= uniform1(generator)) {
                    mol_old.SetCoordinates(mol_ligand.GetCoordinates());
                    energy_old = e;
                    // printf(" Etot = %10.4f    Ea = %10.4f Eb = %10.4f   E_bind = %10.4f\n", e, ea, eb_min, e - (ea + eb_min));
                    // printf("Step: %6i       Etotal = %10.4f kcal/mol\n", step, e);
                    accept += 1;
                // ... or reject.
                } else {
                    mol_ligand.SetCoordinates(mol_old.GetCoordinates());
                    e = energy_old;
                }

                // acceptance_ratio = accept * 100.0 / (step + 1);
                // printf("Step: %6i   acceptance = %6.2f %%   Etotal = %10.4f kcal/mol\n", step + 1, acceptance_ratio, e);
            }

            double ec;
            if (opts.use_mopac == false) {
                ec = minimize_molecule(mol_ligand, ff);
            } else {
                ec = mopac_optimize(mol_ligand);
            }

            // double ec = mopac_optimize(mol_ligand);
            double e_bind = ec - (ea + eb_min);
            acceptance_ratio = accept * 100.0 / (opts.steps + 1);

            // printf("Rotamer: %3i / %3i   Trajectory: %3i / %3i    acceptance = %6.2f %%   E_bind = %10.4f kcal/mol", c + 1, ligand.NumConformers(), n + 1, opts.trajectories, acceptance_ratio, e_bind);
            printf(" %3i / %3i          %3i / %3i           %6.2f %%        %10.4f kcal/mol", c + 1, ligand.NumConformers(), n + 1, opts.trajectories, acceptance_ratio, e_bind);
            conv.Write(&mol_ligand, &ofs);

            // Check if we found the lowest energy
            if (e_bind < e_low) {

                printf(" <---- New lowest\n");
                OpenBabel::OBConversion conv2;
                conv2.SetInAndOutFormats("xyz", "xyz");
                e_low = e_bind;
                std::remove("min.xyz");
                std::ofstream ofs_min("min.xyz");
                std::string title = std::to_string(e_bind);
                mol_ligand.SetTitle(title);
                conv2.Write(&mol_ligand, &ofs_min);
                ofs_min.close();

            } else {
                printf("\n");
            }
        }
    }

    ofs.close();

    double time_elapsed = timer.Elapsed();
    printf("\nOptimized E_bind = %10.4f kcal/mol    Elapsed time = %6.2f seconds\n", 
            e_low, time_elapsed);


    return 0;

}
