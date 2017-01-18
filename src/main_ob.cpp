#include <stdio.h>

#include <iostream>
#include <vector>
// #include <math.h>
// #include <cmath>
#include <string.h>
#include <fstream>
#include<sstream>

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/obutil.h>

int main(int argc, char *argv[]) {

    // OpenBabel::OBPlugin::List("forcefields", "verbose");
    // std::string ff = "MMFF94";
    std::string ff = "UFF";
    const int steps = 200;
    const double crit = 5.0e-4;

    printf("Running McDock 0.1 alpha\n");

    if (argc < 3) {
        printf("Not enough CLI arguments!\n");
        return 1;
    }

    std::string base_file = argv[1];
    std::string dock_file = argv[2];

    OpenBabel::OBMol mol;
    OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("xyz", "xyz");

    std::ifstream ifs;
    ifs.open(base_file.c_str());
    conv.Read(&mol, &ifs);
    ifs.close();


    OpenBabel::OBMol mol2;

    ifs.open(dock_file.c_str());
    conv.Read(&mol2, &ifs);
    ifs.close();

    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);
    if (!pFF) {
        std::cerr << ": could not find forcefield '" << ff << "'." << std::endl;
    }

    pFF->Setup(mol);
    pFF->SteepestDescentInitialize(steps, crit);

    double e = pFF->Energy();
    // std::cout << e << std::endl;
    unsigned int total_steps = 0;
    bool done = true;
    // bool done = false;
    // conv.Write(&mol, &std::cout);
    while (done) {
        done = pFF->SteepestDescentTakeNSteps(1);
        total_steps++;

      if (pFF->DetectExplosion()) {
        std::cerr << "explosion has occured!" << std::endl;
        // conv.Write(&mol, &std::cout);
        return(1);
      } else
        pFF->GetCoordinates(mol);
    }
    // std::cout << total_steps << std::endl;
    e = pFF->Energy();
    // std::cout << e << std::endl;
    // conv.Write(&mol, &std::cout);
   
    e = pFF->Energy();
    // std::cout << e << std::endl;

    OpenBabel::OBAtom *atom;
    std::ofstream ofs("out.xyz");

    mol += mol2;
    pFF->Setup(mol);
    conv.Write(&mol, &std::cout);
    conv.Write(&mol, &ofs);
    e = pFF->Energy();
    std::cout << e << std::endl;

    double energy_new = e;
    double energy_old = energy_new;

    OpenBabel::OBMol mol_old = mol;
    mol_old.SetCoordinates(mol.GetCoordinates());

    for (unsigned int step = 0; step < 1000; step++) {

        OpenBabel::vector3 move;
        move.randomUnitVector();

        for (unsigned int i = mol.NumAtoms() - mol2.NumAtoms() + 1;
            i < mol.NumAtoms() + 1; i++) {
            atom = mol.GetAtom(i);
            atom->SetVector(atom->GetVector() + move);
        }

        pFF->SetCoordinates(mol);
        e = pFF->Energy();

        conv.Write(&mol, &ofs);

        if (e < energy_old) {

            mol_old.SetCoordinates(mol.GetCoordinates());
            energy_old = e;
            std::cout << e << "    accept" << std::endl;
        } else {
            mol.SetCoordinates(mol_old.GetCoordinates());
            e = energy_old;
            std::cout << e << "    reject" << std::endl;
        }
    }
   
   
   
   
    ofs.close();
    // e = pFF->Energy();
    // std::cout << e << std::endl;
    // pFF->SetCoordinates(mol);
    // e = pFF->Energy();
    // std::cout << e << std::endl;
    // Initialize energy function

    // Minimize monomer1
    // Minimize monomer2

    // Start initial conformation loop
    //   Put in same coordinate frame

    //   Start MC loop
    //     Make MC step
    //     Calculate single-point energy
    //     MC-criterion
    //     Save Markov chain

    //   Output lowest for each initial conformation

    return 0;

}
