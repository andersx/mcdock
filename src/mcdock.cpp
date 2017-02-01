#define HAVE_EIGEN
#include <stdio.h>

#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <random>

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/obutil.h>

#include "utils.hpp"

int main(int argc, char *argv[]) {

    OpenBabel::OBStopwatch timer;
    timer.Start();

    printf("Running McDock 0.2 alpha\n");

    if (argc < 5) {
        printf("\nUsage: ./mcdock file1.xyz file2.xyz STEPS TAU NMC\n\n");
        printf("STEPS  = Number of MC steps\n");
        printf("TAU    = Temperature in units of [kB T]\n");
        printf("NMC    = Number of MC trajectories for each conformation \n\n");
        return 1;
    }

    const std::string base_file = argv[1];
    const std::string dock_file = argv[2];

    const int nsteps = atoi(argv[3]);
    const double tau = atoi(argv[4]);
    
    const unsigned int starting_points = atoi(argv[5]);

    // Get OBMols
    OpenBabel::OBMol mol;
    OpenBabel::OBMol ligand;

    OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("xyz", "xyz");
    std::ifstream ifs;

    ifs.open(base_file.c_str());
    conv.Read(&mol, &ifs);
    ifs.close();

    ifs.open(dock_file.c_str());
    conv.Read(&ligand, &ifs);
    ifs.close();

    // Gent energies of pocket molecule
    double ea = minimize_molecule(mol, ff);
    printf("Pocket (minimized) E = %10.4f kcal/mol    file = ", ea / 4.184);
    std::cout << base_file << std::endl;

    OpenBabel::vector3 com = get_com(mol);
    OpenBabel::vector3 temp;
    OpenBabel::vector3 dir;
    OpenBabel::vector3 rot;
    OpenBabel::vector3 move;
    OpenBabel::OBAtom *atom;

    // Center pocket
    for (unsigned int i = 1; i < mol.NumAtoms() + 1; i++) {
       atom = mol.GetAtom(i);
       temp = atom->GetVector();
       temp -= com;
       atom->SetVector(temp);
    }
    
    // Do rotor search for ligand molecule
    printf("Performing rotor search for ligand molecule\n");
    set_conformations(ligand);
    printf("Found %3i rotatable bonds\n", ligand.NumRotors());


    // Make sure output file doesn't already exist.
    std::remove("out.xyz");
    std::ofstream ofs("out.xyz");
  
    // Initialize random number generators
    std::default_random_engine generator;
    std::uniform_real_distribution<double> random_length(0.0, 0.25);
    std::uniform_real_distribution<double> random_angle(0.0, 90.0/180.0*M_PI);
    std::uniform_real_distribution<double> uniform1(0.0, 1.0);

    double theta;
    double eb;  

    double e_low = std::numeric_limits<double>::infinity();
    double eb_min = std::numeric_limits<double>::infinity();

    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);

    for (int c = 0; c < ligand.NumConformers(); ++c) {

        ligand.SetConformer(c);
        eb = minimize_molecule(ligand, ff);
        if (eb < eb_min) eb_min = eb;
        printf("Rotamer %4i     E = %10.4f kcal/mol\n", c, eb / 4.184);
    }
    
    printf("Lowest energy conformation  E = %10.4f kcal/mol\n", eb_min / 4.184);

    for (int c = 0; c < ligand.NumConformers(); ++c) {

        ligand.SetConformer(c);
        eb = minimize_molecule(ligand, ff);
        // printf("Rotamer %4i     E = %10.4f kJ/mol\n", c, eb);

        for (unsigned int n = 0; n < starting_points; n++) {

            dir.randomUnitVector();
            rot.randomUnitVector();
            theta = random_angle(generator); 
            com = get_com(ligand);

            // Center ligand, and give random rotation
            for (unsigned int i = 1; i < ligand.NumAtoms() + 1; i++) {

                atom = ligand.GetAtom(i);
                temp = atom->GetVector();
                temp -= com;
                temp = rotate(temp, rot, theta);
                temp += dir * 4.0;
                atom->SetVector(temp);

            }
           
            // Make object to roll-back rejected MC moves 
            OpenBabel::OBMol mol_ligand = mol;
            mol_ligand += ligand;
            OpenBabel::OBMol mol_old = mol_ligand;
            mol_old.SetCoordinates(mol_ligand.GetCoordinates());

            // Initialize MC energy
            pFF->Setup(mol_ligand);
            double e = pFF->Energy();
            double energy_old = e;
            double delta_e = 0.0;

            // Begin MC simulation
            for (int step = 0; step < nsteps; step++) {

                // Translation move
                move.randomUnitVector();
                move *= random_length(generator);

                for (unsigned int i = mol_ligand.NumAtoms() - ligand.NumAtoms() + 1;
                    i < mol_ligand.NumAtoms() + 1; i++) {
                    atom = mol_ligand.GetAtom(i);
                    atom->SetVector(atom->GetVector() + move);
                }

                // Rotation move
                rot.randomUnitVector();
                theta = random_angle(generator); 

                // Get center of mass
                com.Set(0.0, 0.0, 0.0);

                for (unsigned int i = mol_ligand.NumAtoms() - ligand.NumAtoms() + 1;
                    i < mol_ligand.NumAtoms() + 1; i++) {
                    atom = mol_ligand.GetAtom(i);
                    com += atom->GetVector();
                }

                com /= mol_ligand.NumAtoms() - ligand.NumAtoms();

                // Rotate around center of mass
                for (unsigned int i = mol_ligand.NumAtoms() - ligand.NumAtoms() + 1;
                    i < mol_ligand.NumAtoms() + 1; i++) {

                    atom = mol_ligand.GetAtom(i);
                    temp = atom->GetVector();

                    temp -= com;
                    temp = rotate(temp, rot, theta);
                    temp += com;

                    atom->SetVector(temp);
                }

                // Evaluate total energy
                pFF->SetCoordinates(mol_ligand);
                e = pFF->Energy();
                delta_e =  e - energy_old;

                // Metropolis-Hastings MC criterion, accept ...
                if (std::exp( - delta_e / tau) >= uniform1(generator)) {
                    mol_old.SetCoordinates(mol_ligand.GetCoordinates());
                    energy_old = e;
                    // printf(" Etot = %10.4f    Ea = %10.4f Eb = %10.4f   E_bind = %10.4f\n", e, ea, eb_min, e - (ea + eb_min));
                // ... or reject.
                } else {
                    mol_ligand.SetCoordinates(mol_old.GetCoordinates());
                    e = energy_old;
                }
            }

            double ec = minimize_molecule(mol_ligand, ff);
            double e_bind = ec - (ea + eb_min);

            printf("Rotamer: %3i / %3i   Trajectory: %3i / %3i   E_bind = %10.4f kcal/mol", c + 1, ligand.NumConformers(), n + 1, starting_points, e_bind / 4.184);
            conv.Write(&mol_ligand, &ofs);

            if (e_bind < e_low) {

                OpenBabel::OBConversion conv2;
                conv2.SetInAndOutFormats("xyz", "xyz");
                printf("      <---- New best\n");
                e_low = e_bind;
                std::remove("min.xyz");
                std::ofstream ofs_min("min.xyz");
                conv2.Write(&mol_ligand, &ofs_min);
                ofs_min.close();

            } else {
                printf("\n");
            }
        }
    }
    double time_elapsed = timer.Elapsed();
    printf("Optimized E_bind = %10.4f kcal/mol    Elapsed time = %4.2f seconds\n", 
            e_low / 4.184, time_elapsed);


    ofs.close();

    return 0;

}
