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


void minimize_molecule(OpenBabel::OBMol &mol, const std::string &ff) {

    OpenBabel::OBStopwatch timer;
    timer.Start();
    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);
    pFF->Setup(mol);

    double e = pFF->Energy();

    printf("Initial energy of molecule: %10.3f kJ/mol\n", e);
    printf("Minimizing ....");

    const int steps = 50;
    const double crit = 5.0e-4;

    // pFF->SteepestDescentInitialize(steps, crit);
    pFF->ConjugateGradientsInitialize(steps, crit);

    bool done = true;

    while (done) {

        // done = pFF->SteepestDescentTakeNSteps(1);
        done = pFF->ConjugateGradientsTakeNSteps(1);

        if (pFF->DetectExplosion()) {

            std::cerr << "explosion has occured!" << std::endl;
            exit(1);

        } else {

            pFF->GetCoordinates(mol);
        }

    }
    e = pFF->Energy();
    double time_elapsed = timer.Elapsed();
    printf(" %5.2f seconds\n", time_elapsed);
    printf("Minimzed energy of molecule: %10.3f kJ/mol\n", e);

}


int main(int argc, char *argv[]) {

    std::string ff = "MMFF94";

    printf("Running McDock 0.1 alpha Conformer Search\n");

    if (argc < 2) {
        printf("Not enough CLI arguments!\n");
        return 1;
    }

    std::string base_file = argv[1];

    OpenBabel::OBMol mol;

    OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("xyz", "xyz");
    std::ifstream ifs;

    ifs.open(base_file.c_str());
    conv.Read(&mol, &ifs);
    ifs.close();

    minimize_molecule(mol, ff);

    std::remove("conformers.xyz");
    std::ofstream ofs("conformers.xyz");

    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);
    pFF->Setup(mol);

    std::cout << "..title = " << mol.GetTitle() << std::endl;
    std::cout << "..number of rotatable bonds = " << mol.NumRotors() << std::endl;

    double rmsd_cutoff;
    double energy_cutoff;
    unsigned int conf_cutoff;
    bool verbose;
    bool include_original;

    rmsd_cutoff = 0.5;
    energy_cutoff = 50.0;
    conf_cutoff = 1000; // 1 Million
    verbose = false;
    include_original = false;

    pFF->DiverseConfGen(rmsd_cutoff, conf_cutoff, energy_cutoff, true);

    pFF->GetConformers(mol);
    int nconfs = include_original ? mol.NumConformers() : mol.NumConformers() - 1;

    std::cout << "..generated " << nconfs << " conformers" << std::endl;

    unsigned int c = include_original ? 0 : 1;

    for (; c < mol.NumConformers(); ++c) {
      mol.SetConformer(c);
        minimize_molecule(mol, ff);
       double e = pFF->Energy();
        std::cout << c << " id " << e << " energy" << std::endl;

      // if((conv.GetOutFormat())->WriteMolecule(mol, conv))
             conv.Write(&mol, &ofs);
      //  break;
    }

    // conv.Write(&mol, &ofs);

    // double e = pFF->Energy();
    // printf("Initial energy of complex: %10.3f kJ/mol\n", e);

    // double energy_old = e;
    // double delta_e = 0.0;

    // OpenBabel::OBMol mol_old = mol;
    // mol_old.SetCoordinates(mol.GetCoordinates());

    // std::default_random_engine generator;
    // std::uniform_real_distribution<double> random_length(0.0, 0.25);
    // std::uniform_real_distribution<double> random_angle(0.0, 90.0/180.0*M_PI);

    // std::uniform_real_distribution<double> uniform1(0.0, 1.0);

    // OpenBabel::OBAtom *atom;
    // OpenBabel::vector3 move;
    // OpenBabel::vector3 com;
    // OpenBabel::vector3 dir;
    // OpenBabel::vector3 temp;
    // double t;

    // double tau = 0.1;
    // for (unsigned int step = 0; step < 5000; step++) {

    //     // Translation move
    //     move.randomUnitVector();
    //     move *= random_length(generator);

    //     for (unsigned int i = mol.NumAtoms() - mol2.NumAtoms() + 1;
    //         i < mol.NumAtoms() + 1; i++) {
    //         atom = mol.GetAtom(i);
    //         atom->SetVector(atom->GetVector() + move);
    //     }

    //     // Rotation move
    //     dir.randomUnitVector();
    //     t = random_angle(generator); 

    //     // Get center of mass
    //     com.Set(0.0, 0.0, 0.0);

    //     for (unsigned int i = mol.NumAtoms() - mol2.NumAtoms() + 1;
    //         i < mol.NumAtoms() + 1; i++) {
    //         atom = mol.GetAtom(i);
    //         com += atom->GetVector();
    //     }

    //     com /= mol.NumAtoms() - mol2.NumAtoms();

    //     // Rotate around center of mass
    //     for (unsigned int i = mol.NumAtoms() - mol2.NumAtoms() + 1;
    //         i < mol.NumAtoms() + 1; i++) {

    //         atom = mol.GetAtom(i);
    //         temp = atom->GetVector();

    //         temp -= com;
    //         temp = rotate(temp, dir, t);
    //         temp += com;

    //         atom->SetVector(temp);
    //     }

    //     pFF->SetCoordinates(mol);
    //     e = pFF->Energy();

    //     delta_e =  e - energy_old;

    //     if (std::exp( - delta_e / tau) >= uniform1(generator)) {

    //         mol_old.SetCoordinates(mol.GetCoordinates());
    //         energy_old = e;
    //         std::cout << step << "  " << e << "    accept    " <<  std::exp( - delta_e / tau) << "    " << delta_e << std::endl;
    //         conv.Write(&mol, &ofs);

    //     } else {
    //         mol.SetCoordinates(mol_old.GetCoordinates());
    //         e = energy_old;
    //     }
    // }
   
    // conv.Write(&mol, &ofs);
    // minimize_molecule(mol, ff);
    // ofs.close();

    return 0;

}
