#ifndef SRC_UTILS_HPP
#define SRC_UTILS_HPP
#define HAVE_EIGEN

namespace McDock {

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

void move_molecule(OpenBabel::OBMol &mol, OpenBabel::vector3 move,
        int startid = 1, int endid = -1) {

    if (endid == -1) endid = mol.NumAtoms() + 1;

    OpenBabel::vector3 temp;
    OpenBabel::OBAtom *atom;

    for (int i = startid; i < endid; i++) {
       atom = mol.GetAtom(i);
       temp = atom->GetVector();
       temp += move;
       atom->SetVector(temp);
    }

}

void rotate_molecule(OpenBabel::OBMol &mol, OpenBabel::vector3 direction, 
        double theta, int startid = 1, int endid = -1) {

    if (endid == -1) endid = mol.NumAtoms() + 1;
    
    OpenBabel::vector3 com;
    com.Set(0.0, 0.0, 0.0);
    OpenBabel::OBAtom *atom;

    for (int i = startid; i < endid; i++) {

        atom = mol.GetAtom(i);
        com += atom->GetVector();
    }
    
    com /= (double)(endid - startid);


    OpenBabel::vector3 temp;
    // Center ligand, and give random rotation
    for (int i = startid; i < endid; i++) {

        atom = mol.GetAtom(i);
        temp = atom->GetVector();
        temp -= com;
        temp = rotate(temp, direction, theta);
        temp += com;
        atom->SetVector(temp);

    }


}


OpenBabel::OBMol readfile(std::string filename) {

    OpenBabel::OBMol mol;

    OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("xyz", "xyz");
    std::ifstream ifs;

    ifs.open(filename.c_str());
    conv.Read(&mol, &ifs);
    ifs.close();

    return mol;

}

double mopac_energy(OpenBabel::OBMol &mol) {

    std::remove("temp.mop");
    std::remove("temp.out");
    std::remove("temp.arc");

    OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("moo", "mop");

    std::ofstream ofs("temp.mop");
    conv.Write(&mol, &ofs);
    ofs.close();
    
    // int success = system("./run_mopac_1scf");
    int success = 0;
    success += system("sed -i \"s/PUT KEYWORDS HERE/PM6-D3H4 1SCF/\" temp.mop");
    success += system("LD_LIBRARY_PATH=/opt/mopac:$LD_LIBRARY_PATH OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/mopac/MOPAC2016.exe temp.mop 2> /dev/null ");
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

    int success = 0;
    success += system("sed -i \"s/PUT KEYWORDS HERE/PM6-D3H4 EF/\" temp.mop");
    success += system("LD_LIBRARY_PATH=/opt/mopac:$LD_LIBRARY_PATH OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /opt/mopac/MOPAC2016.exe temp.mop 2> /dev/null ");
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

void set_conformations(OpenBabel::OBMol &mol, std::string ff){

    // if (ff.compare("PM6-D3H4") == 0) ff = "MMFF94";
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


void save_xyz(OpenBabel::OBMol mol, std::string filename) {

    OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("xyz", "xyz");

    std::ofstream ofs(filename);
    
    //std::string title = std::to_string(e_bind);
    //mol_ligand.SetTitle(title);
    
    conv.Write(&mol, &ofs);
    ofs.close();


}

// Prints help
void print_help() {

        printf("Usage: ./mcdock --target file1.xyz --ligand file2.xyz [--args]\n\n");

        printf("Optional arguments:\n");

        printf("--energy [string]       Potential energy function \"MMFF94\" (default),\n");
        printf("                        \"UFF\", \"PM6-D3H4\" (requires MOPAC). \n");
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
                    opts.ff = "MMFF94";
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



} // Namespace McDock

#endif
