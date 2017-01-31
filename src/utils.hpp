#ifndef SRC_UTILS_HPP
#define SRC_UTILS_HPP


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

    OpenBabel::OBStopwatch timer;
    timer.Start();
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

#endif
