/*
 * Greywolf.h  --  Grey Wolf Optimizer (declaration only)
 *
 * Formulated for RNA secondary structure prediction with explicit Wolf class.
 * Inherits from RNA_Core and follows OPTIMIZER_CONTRACT.md.
 */
#ifndef GREYWOLF_H
#define GREYWOLF_H

#include "core.h"

struct Wolf {
    Individual position;
    double fitness = 1e20;
    string structure;
};

class Greywolf : public RNA_Core {
public:
    vector<Wolf> wolves;
    Individual alphaPos, betaPos, deltaPos;
    double alphaFE = 1e20, betaFE = 1e20, deltaFE = 1e20;

    Greywolf(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("Greywolf", seqID, seq, known, popSize, iterations, param) {}

    void run() override;
};

#endif // GREYWOLF_H
