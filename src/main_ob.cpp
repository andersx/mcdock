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
    std::string ff = "MMFF94";
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

    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);
    if (!pFF) {
        std::cerr << ": could not find forcefield '" << ff << "'." << std::endl;
    }

    pFF->Setup(mol);
    pFF->SteepestDescentInitialize(steps, crit);

    double e = pFF->Energy();
    std::cout << e << std::endl;
    unsigned int total_steps = 0;
    bool done = true;
    // bool done = false;
    while (done) {
        done = pFF->SteepestDescentTakeNSteps(1);
        total_steps++;

      if (pFF->DetectExplosion()) {
        std::cerr << "explosion has occured!" << std::endl;
        conv.Write(&mol, &std::cout);
        return(1);
      } else
        pFF->GetCoordinates(mol);
    }
    std::cout << total_steps << std::endl;
    e = pFF->Energy();
    std::cout << e << std::endl;
    conv.Write(&mol, &std::cout);
   
    OpenBabel::OBAtom *atom;

    for (OpenBabel::OBMolAtomIter a(mol); a; ++a ) {
      atom = mol.GetAtom(a->GetIdx());
      std::cout << atom->GetVector() << std::endl;
    }

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
