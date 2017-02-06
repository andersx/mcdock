#ifndef SRC_UTILS_HPP
#define SRC_UTILS_HPP

std::string ff = "MMFF94";

double mopac_energy(OpenBabel::OBMol &mol) {

    std::remove("temp.mop");
    std::remove("temp.out");
    std::remove("temp.arc");

    OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("moo", "mop");

    std::ofstream ofs("temp.mop");
    conv.Write(&mol, &ofs);
    ofs.close();

    int success = system("./run_mopac_1scf");
    if (success > 1) printf("Error running MOPAC!");

    OpenBabel::OBMol molout;

    std::ifstream ifs;
    ifs.open("temp.out");
    conv.Read(&molout, &ifs);
    ifs.close();

    std::remove("temp.mop");
    std::remove("temp.out");
    std::remove("temp.arc");

    return  molout.GetEnergy();

}
                    
double mopac_optimize(OpenBabel::OBMol &mol) {

    std::remove("temp.mop");
    std::remove("temp.out");
    std::remove("temp.arc");

    OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("moo", "mop");

    std::ofstream ofs("temp.mop");
    conv.Write(&mol, &ofs);
    ofs.close();

    int success = system("./run_mopac_opt");
    if (success > 1) printf("Error running MOPAC!");

    OpenBabel::OBMol molout;

    std::ifstream ifs;
    ifs.open("temp.out");
    conv.Read(&molout, &ifs);
    ifs.close();

    // Segfaults for unknown random reasons??
    // mol.SetCoordinates(molout.GetCoordinates());

    // Workaround for above problem
    for (unsigned int i = 1; i < mol.NumAtoms() + 1; i++) {
        OpenBabel::OBAtom *atom = mol.GetAtom(i);
        atom->SetVector((molout.GetAtom(i))->GetVector());
    }

    mol.SetEnergy(molout.GetEnergy());

    std::remove("temp.mop");
    std::remove("temp.out");
    std::remove("temp.arc");

    return  molout.GetEnergy();

}

// Returns a rotated vector, rotated by T
OpenBabel::vector3 rotate(const OpenBabel::vector3 &V, 
        const OpenBabel::vector3 &J, const double T) {

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

    double a = (u * (u*x + v*y + w*z) + (x * (v*v + w*w) - u * (v*y + w*z)) 
            * cost + norm * (-w*y + v*z) * sint) * inv_norm_sqrt;
    double b = (v * (u*x + v*y + w*z) + (y * (u*u + w*w) - v * (u*x + w*z)) 
            * cost + norm * ( w*x - u*z) * sint) * inv_norm_sqrt;
    double c = (w * (u*x + v*y + w*z) + (z * (u*u + v*v) - w * (u*x + v*y)) 
            * cost + norm * (-v*x + u*y) * sint) * inv_norm_sqrt;

    OpenBabel::vector3 rotated;
    rotated.Set(a, b, c);

    return rotated;
}
double minimize_molecule(OpenBabel::OBMol &mol, const std::string &ff) {

    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);
    pFF->Setup(mol);

    double e = pFF->Energy();
    // printf("E_before = %10.4f kJ/mol   %i\n", e, mol.NumAtoms());

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

    // printf("E_after = %10.4f kJ/mol %i\n", e, mol.NumAtoms());
    return e;

}

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

#endif
