/*
 * GeneticAlgorithm.h  --  Genetic Algorithm (declaration only)
 *
 * Formulated for RNA secondary structure prediction with explicit Chromosome class.
 * Inherits from RNA_Core and follows OPTIMIZER_CONTRACT.md.
 */
#ifndef GENETICALGORITHM_H
#define GENETICALGORITHM_H

#include "core.h"

struct Chromosome {
    Individual genes;
    double fitness = 1e20;
    string structure;
};

class GeneticAlgorithm : public RNA_Core {
public:
    double crossoverRate = 0.80;
    double mutationRate = 0.20;
    int tournamentSize = 3;
    vector<Chromosome> population_pool;

    GeneticAlgorithm(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("GeneticAlgorithm", seqID, seq, known, popSize, iterations, param) {}

    int tournamentSelect();
    Individual crossover(const Individual &p1, const Individual &p2);
    Individual mutate(const Individual &ind);

    void run() override;
};

#endif // GENETICALGORITHM_H
