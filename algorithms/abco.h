/*
 * abco.h  --  Artificial Bee Colony Optimization (declaration only)
 */
#ifndef ABCO_H
#define ABCO_H

#include "core.h"

class ABCO : public RNA_Core {
public:
    // ABCO hyper-parameter
    int ABANDON_LIMIT = 20;   // scout threshold: trials before a source resets

    ABCO(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("ABCO", seqID, seq, known, popSize, iterations, param) {}

    double     computeFit(double f);
    Individual generateNeighbour(int i);
    void run() override;
};

#endif // ABCO_H
