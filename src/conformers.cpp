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

    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);
    pFF->Setup(mol);

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

    printf("Found %3i rotatable bonds\n", mol.NumRotors());

    double rmsd_cutoff = 0.5;
    double energy_cutoff = 50.0;
    unsigned int conf_cutoff = 1000000; // 1 Million
    bool verbose = false;

    pFF->DiverseConfGen(rmsd_cutoff, conf_cutoff, energy_cutoff, verbose);
    pFF->GetConformers(mol);

    double e;
    for (int c = 0; c < mol.NumConformers(); ++c) {
        mol.SetConformer(c);
        minimize_molecule(mol, ff);
        e = pFF->Energy();
        printf("Rotamer %4i     E = %10.4f kJ/mol\n", c, e);
        conv.Write(&mol, &ofs);
    }

    return 0;

}
