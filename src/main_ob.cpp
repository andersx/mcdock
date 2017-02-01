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


// Returns a rotated vector, rotated by T
OpenBabel::vector3 rotate(const OpenBabel::vector3 &V, const OpenBabel::vector3 &J, const double T) {

    double x = V.x();
    double y = V.y();
    double z = V.z();

    double u = J.x();
    double v = J.y();
    double w = J.z();

    double norm = std::sqrt(u*u + v*v + w*w);
    double inv_norm_sqrt = 1.0 / (norm * norm);
    double sint = std::sin(T);
    double cost = std::cos(T);

    double a = (u * (u*x + v*y + w*z) + (x * (v*v + w*w) - u * (v*y + w*z)) * cost + norm * (-w*y + v*z) * sint) * inv_norm_sqrt;
    double b = (v * (u*x + v*y + w*z) + (y * (u*u + w*w) - v * (u*x + w*z)) * cost + norm * ( w*x - u*z) * sint) * inv_norm_sqrt;
    double c = (w * (u*x + v*y + w*z) + (z * (u*u + v*v) - w * (u*x + v*y)) * cost + norm * (-v*x + u*y) * sint) * inv_norm_sqrt;

    OpenBabel::vector3 rotated;
    rotated.Set(a, b, c);

    return rotated;
}

double minimize_molecule(OpenBabel::OBMol &mol, const std::string &ff) {

    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);
    pFF->Setup(mol);

    double e = pFF->Energy();

    const int steps = 200;
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

    return e;

}


int main(int argc, char *argv[]) {

    // OpenBabel::OBPlugin::List("forcefields", "verbose");
    std::string ff = "MMFF94";
    // std::string ff = "UFF";

    printf("Running McDock 0.1 alpha\n");

    if (argc < 5) {
        printf("\nUsage: ./mcdock file1.xyz file2.xyz STEPS TAU\n\n");
        printf("STEPS  = Number of MC steps\n");
        printf("TAU    = Temperature in units of [kB T]\n\n");
        return 1;
    }

    std::string base_file = argv[1];
    std::string dock_file = argv[2];

    const int nsteps = atoi(argv[3]);
    const double tau = atoi(argv[4]);

    OpenBabel::OBMol mol;
    OpenBabel::OBMol mol2;

    OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("xyz", "xyz");
    std::ifstream ifs;

    ifs.open(base_file.c_str());
    conv.Read(&mol, &ifs);
    ifs.close();

    ifs.open(dock_file.c_str());
    conv.Read(&mol2, &ifs);
    ifs.close();

    double ea = minimize_molecule(mol, ff);
    printf("Molecule A (minimized) E = %10.4f kcal/mol    file = ", ea);
    std::cout << base_file << std::endl;

    double eb = minimize_molecule(mol2, ff);
    printf("Molecule B (minimized) E = %10.4f kcal/mol    file = ", eb);
    std::cout << dock_file << std::endl;

    std::remove("out.xyz");
    std::ofstream ofs("out.xyz");

    mol += mol2;

    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);
    pFF->Setup(mol);

    conv.Write(&mol, &ofs);

    double e = pFF->Energy();
    printf("Initial binding energy of complex: %10.3f kJ/mol\n", e - (ea + eb));

    double energy_old = e;
    double delta_e = 0.0;

    OpenBabel::OBMol mol_old = mol;
    mol_old.SetCoordinates(mol.GetCoordinates());

    std::default_random_engine generator;
    std::uniform_real_distribution<double> random_length(0.0, 0.25);
    std::uniform_real_distribution<double> random_angle(0.0, 90.0/180.0*M_PI);

    std::uniform_real_distribution<double> uniform1(0.0, 1.0);

    OpenBabel::OBAtom *atom;
    OpenBabel::vector3 move;
    OpenBabel::vector3 com;
    OpenBabel::vector3 dir;
    OpenBabel::vector3 temp;
    double t;

    // double tau = 0.1;

    for (int step = 0; step < nsteps; step++) {

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

            printf("Step: %10u  E_bind = %10.4f kcal/mol\n", step, (e - (ea + eb)));

            conv.Write(&mol, &ofs);

        } else {
            mol.SetCoordinates(mol_old.GetCoordinates());
            e = energy_old;
        }
    }
   

    printf("Minimizing final conformation ... ");
    conv.Write(&mol, &ofs);

    double ec = minimize_molecule(mol, ff);
    ofs.close();
    printf("done!\n");
    printf("Wrote out.xyz\n");
    printf("Final E_bind = %10.4f kcal/mol\n", (ec - (ea + eb)));

    return 0;

}
