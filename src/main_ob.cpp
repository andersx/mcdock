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

OpenBabel::vector3 rotate(const OpenBabel::vector3 &V, const OpenBabel::vector3 &J, const double T) {

    double x = V.x();
    double y = V.y();
    double z = V.z();

    double u = J.x();
    double v = J.y();
    double w = J.z();

    double a = (u*(u*x + v*y + w*z) + (x * (v*v + w*w) - u *(v*y + w*z)) * std::cos(T) + std::sqrt(u*u + v*v + w*w) * (-w*y + v*z) * std::sin(T)) / (u*u + v*v + w*w);
    double b = (v*(u*x + v*y + w*z) + (y * (u*u + w*w) - v *(u*x + w*z)) * std::cos(T) + std::sqrt(u*u + v*v + w*w) * ( w*x - u*z) * std::sin(T)) / (u*u + v*v + w*w);
    double c = (w*(u*x + v*y + w*z) + (z * (u*u + v*v) - w *(u*x + v*y)) * std::cos(T) + std::sqrt(u*u + v*v + w*w) * (-v*x + u*y) * std::sin(T)) / (u*u + v*v + w*w);

    OpenBabel::vector3 rotated;
    rotated.Set(a, b, c);

    return rotated;
}

int main(int argc, char *argv[]) {

    // OpenBabel::OBPlugin::List("forcefields", "verbose");
    std::string ff = "MMFF94";
    // std::string ff = "UFF";
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

    double energy_old = e;
    double delta_e = 0.0;

    OpenBabel::OBMol mol_old = mol;
    mol_old.SetCoordinates(mol.GetCoordinates());


    std::default_random_engine generator;
    std::uniform_real_distribution<double> random_length(0.0, 0.25);
    std::uniform_real_distribution<double> random_angle(0.0, 10.0);

    std::uniform_real_distribution<double> uniform1(0.0, 1.0);

    OpenBabel::vector3 move;
    OpenBabel::vector3 com;
    OpenBabel::vector3 dir;
    OpenBabel::vector3 temp;
    double t;

    double tau = 0.25;
    for (unsigned int step = 0; step < 100000; step++) {

        // Translation move
        move.randomUnitVector();
        move *= random_length(generator);

        for (unsigned int i = mol.NumAtoms() - mol2.NumAtoms() + 1;
            i < mol.NumAtoms() + 1; i++) {
            atom = mol.GetAtom(i);
            atom->SetVector(atom->GetVector() + move);
        }

        // Rotation move

        dir.randomUnitVector();
        t = random_angle(generator); 

        // Get center of mass
        com.Set(0.0, 0.0, 0.0);

        for (unsigned int i = mol.NumAtoms() - mol2.NumAtoms() + 1;
            i < mol.NumAtoms() + 1; i++) {
            atom = mol.GetAtom(i);
            com += atom->GetVector();
        }

        com /= mol.NumAtoms() - mol2.NumAtoms();

        // Rotate around center of mass
        for (unsigned int i = mol.NumAtoms() - mol2.NumAtoms() + 1;
            i < mol.NumAtoms() + 1; i++) {

            atom = mol.GetAtom(i);
            temp = atom->GetVector();

            temp -= com;
            temp = rotate(temp, dir, t);
            temp += com;

            atom->SetVector(temp);
        }

        pFF->SetCoordinates(mol);
        e = pFF->Energy();

        delta_e =  e - energy_old;

            if (std::exp( - delta_e / tau) >= uniform1(generator)) {

            mol_old.SetCoordinates(mol.GetCoordinates());
            energy_old = e;
            std::cout << step << "  " << e << "    accept    " <<  std::exp( - delta_e / tau) << "    " << delta_e << std::endl;
            conv.Write(&mol, &ofs);

        } else {
            mol.SetCoordinates(mol_old.GetCoordinates());
            e = energy_old;
            // std::cout << step << "  " << e << "    reject" << std::endl;
        }
    }
   
    ofs.close();

    return 0;

}
