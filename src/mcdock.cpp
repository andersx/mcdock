#define HAVE_EIGEN
#include <stdio.h>

#include <iostream>
#include <vector>
#include <math.h>
// #include <cmath>
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

std::string ff = "MMFF94";

void set_conformations(OpenBabel::OBMol &mol){

    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);
    pFF->Setup(mol);

    double rmsd_cutoff = 0.5;
    double energy_cutoff = 50.0;
    unsigned int conf_cutoff = 1000000; // 1 Million
    bool verbose = false;

    pFF->DiverseConfGen(rmsd_cutoff, conf_cutoff, energy_cutoff, verbose);
    pFF->GetConformers(mol);


}


OpenBabel::vector3 get_com(OpenBabel::OBMol mol) {

    OpenBabel::vector3 com;
    com.Set(0.0, 0.0, 0.0);
    OpenBabel::OBAtom *atom;

    OpenBabel::vector3 temp;
    for (unsigned int i = 1; i < mol.NumAtoms() + 1; i++) {

        atom = mol.GetAtom(i);
        temp = atom->GetVector();
        com += temp;
    }
    
    com /= mol.NumAtoms();

    return com;

}


int main(int argc, char *argv[]) {

    // OpenBabel::OBPlugin::List("forcefields", "verbose");
    // std::string ff = "UFF";

    printf("Running McDock 0.2 alpha\n");

    if (argc < 5) {
        printf("\nUsage: ./mcdock file1.xyz file2.xyz STEPS TAU\n\n");
        printf("STEPS  = Number of MC steps\n");
        printf("TAU    = Temperature in units of [kB T]\n\n");
        return 1;
    }

    const std::string base_file = argv[1];
    const std::string dock_file = argv[2];

    const int nsteps = atoi(argv[3]);
    const double tau = atoi(argv[4]);
    
    const unsigned int starting_points = 20;

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

    // Gent energies of monomers Mol1 and Mol2
    double ea = minimize_molecule(mol, ff);
    printf("Pocket (minimized) E = %10.4f kJ/mol    file = ", ea);
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

    // double eb = minimize_molecule(ligand, ff);
    // printf("Ligand (minimized) E = %10.4f kJ/mol    file = ", eb);
    // std::cout << dock_file << std::endl;

    printf("Performing rotor search for ligand molecule\n");
    set_conformations(ligand);
    printf("Found %3i rotatable bonds\n", ligand.NumRotors());


    std::remove("out.xyz");
    std::ofstream ofs("out.xyz");
  
    std::default_random_engine generator;
    std::uniform_real_distribution<double> random_length(0.0, 0.25);
    std::uniform_real_distribution<double> random_angle(0.0, 90.0/180.0*M_PI);

    std::uniform_real_distribution<double> uniform1(0.0, 1.0);


    double theta;
    double eb;  

    double e_low = std::numeric_limits<double>::infinity();

    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);

    for (int c = 0; c < ligand.NumConformers(); ++c) {

        ligand.SetConformer(c);
        eb = minimize_molecule(ligand, ff);
        printf("Rotamer %4i     E = %10.4f kJ/mol\n", c, eb);

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
            
            OpenBabel::OBMol mol_ligand = mol;
            mol_ligand += ligand;
            OpenBabel::OBMol mol_old = mol_ligand;

            mol_old.SetCoordinates(mol_ligand.GetCoordinates());

            pFF->Setup(mol_ligand);
            double e = pFF->Energy();
            // printf("Starting conformation %4i Initial binding energy of complex: %10.3f kJ/mol\n", n, e - (ea + eb));

            double energy_old = e;
            double delta_e = 0.0;

            // conv.Write(&mol_ligand, &ofs);

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

                pFF->SetCoordinates(mol_ligand);
                e = pFF->Energy();

                delta_e =  e - energy_old;

                if (std::exp( - delta_e / tau) >= uniform1(generator)) {

                    mol_old.SetCoordinates(mol_ligand.GetCoordinates());
                    energy_old = e;

                    // printf("Step: %10u  E_bind = %10.4f kJ/mol\n", step, e - (ea + eb));

                    // conv.Write(&mol_ligand, &ofs);

                } else {
                    mol_ligand.SetCoordinates(mol_old.GetCoordinates());
                    e = energy_old;
                }
            }

            double ec = minimize_molecule(mol_ligand, ff);
            double e_bind = ec - (ea + eb);
            printf("Final E_bind = %10.4f kJ/mol", e_bind);
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

        // conv.Write(&ligand, &ofs);
    }


    ofs.close();

    return 0;

}
