/*
 * sa.h  --  Simulated Annealing (declaration only)
 */
#ifndef SA_H
#define SA_H

#include "core.h"

class SA : public RNA_Core {
public:
    // SA hyper-parameters
    double T0_MIN = 1.0;     // lower bound on the runtime-derived start temperature
    double T_MIN  = 1e-8;    // temperature floor (keeps exp() well defined)

    SA(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("SA", seqID, seq, known, popSize, iterations, param) {}

    Individual neighbour(const Individual &ind);
    void run() override;
};

#endif // SA_H
