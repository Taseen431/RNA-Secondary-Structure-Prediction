/*
 * ga.h  --  Genetic Algorithm (declaration only)
 */
#ifndef GA_H
#define GA_H

#include "core.h"

class GA : public RNA_Core {
public:
    // GA hyper-parameters
    double CROSSOVER_RATE = 0.80;
    double MUTATION_RATE = 0.20;
    int    TOURNAMENT_SIZE = 3;

    GA(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("GA", seqID, seq, known, popSize, iterations, param) {}

    int tournamentSelect(const vector<double> &fitVals);
    Individual crossover(const Individual &a, const Individual &b);
    Individual mutate(const Individual &ind);
    void run() override;
};

#endif // GA_H
