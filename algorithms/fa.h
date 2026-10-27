/*
 * fa.h  --  Firefly Algorithm (declaration only)
 */
#ifndef FA_H
#define FA_H

#include "core.h"

class FA : public RNA_Core {
public:
    // FA hyper-parameters ("parameter" supplies gamma, the light absorption
    // coefficient -- the one free constant of canonical FA)
    double BETA0            = 1.0;    // attractiveness at zero distance
    double ALPHA0           = 0.20;   // initial randomization strength
    double ALPHA_DAMP       = 0.97;   // geometric damping applied to alpha per iteration
    double SELECT_THRESHOLD = 0.50;   // decoder threshold (kept equal to PSO's configured value)

    FA(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("FA", seqID, seq, known, popSize, iterations, param) {}

    Individual decode(const vector<double> &position);
    void run() override;
};

#endif // FA_H
