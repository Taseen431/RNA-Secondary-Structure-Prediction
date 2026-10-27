/*
 * boa.h  --  Butterfly Optimization Algorithm (declaration only)
 */
#ifndef BOA_H
#define BOA_H

#include "core.h"

class BOA : public RNA_Core {
public:
    BOA(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("BOA", seqID, seq, known, popSize, iterations, param) {}

    // BOA-specific search operators
    pair<Individual, Individual> separateGlobalSearch(const Individual &B);
    Individual reverseGlobalSearch(const Individual &B, const Individual &Bp);
    void exchangeLocalSearch(const Individual &B1, const Individual &B2,
                             Individual &N1, Individual &N2, Individual &N3, Individual &N4);
    pair<Individual, Individual> mergeLocalSearch(const Individual &B1, const Individual &B2);

    void run() override;
};

#endif // BOA_H
