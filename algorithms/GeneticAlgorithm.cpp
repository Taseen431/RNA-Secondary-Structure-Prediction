/*
 * GeneticAlgorithm.cpp  --  Genetic Algorithm (implementation)
 *
 * Formulated for RNA secondary structure prediction with explicit Chromosome class.
 * Inherits from RNA_Core and follows OPTIMIZER_CONTRACT.md.
 */

#include "GeneticAlgorithm.h"

int GeneticAlgorithm::tournamentSelect() {
    int n = (int)population_pool.size();
    int bestIdx = randInt(0, n - 1);
    double bestVal = population_pool[bestIdx].fitness;

    for (int t = 1; t < tournamentSize; t++) {
        int idx = randInt(0, n - 1);
        if (population_pool[idx].fitness < bestVal) {
            bestVal = population_pool[idx].fitness;
            bestIdx = idx;
        }
    }
    return bestIdx;
}

Individual GeneticAlgorithm::crossover(const Individual &p1, const Individual &p2) {
    if (p1.empty()) return p2;
    if (p2.empty()) return p1;

    int c1 = randInt(0, (int)p1.size());
    int c2 = randInt(0, (int)p2.size());

    Individual child(p1.begin(), p1.begin() + c1);
    for (size_t i = c2; i < p2.size(); i++) {
        child.push_back(p2[i]);
    }

    if (child.empty()) child.push_back(randInt(1, (int)stemList.size()));
    return optimizeButterfly(child);
}

Individual GeneticAlgorithm::mutate(const Individual &ind) {
    Individual c = ind;
    int m = (int)stemList.size();
    if (m == 0) return c;

    int op = randInt(0, 2);
    if (op == 0 && !c.empty()) {
        c[randInt(0, (int)c.size() - 1)] = randInt(1, m);
    } else if (op == 1) {
        if ((int)c.size() < MAX_LEN) c.push_back(randInt(1, m));
    } else {
        if (c.size() > 1) c.erase(c.begin() + randInt(0, (int)c.size() - 1));
    }
    return optimizeButterfly(c);
}

void GeneticAlgorithm::run() {
    extractStems();
    int m = (int)stemList.size();
    if (m == 0) return;

    int popSize = populationSize;
    if (popSize < 4) popSize = 4;

    population_pool.clear();
    population_pool.resize(popSize);
    population.clear();

    for (int i = 0; i < popSize; i++) {
        population_pool[i].genes = randomPermutation();
        population_pool[i].fitness = calculateFitness(population_pool[i].genes, population_pool[i].structure);
        population.push_back(population_pool[i].genes);
    }

    int bestIdx = 0;
    for (int i = 1; i < popSize; i++) {
        if (population_pool[i].fitness < population_pool[bestIdx].fitness) bestIdx = i;
    }

    globalBest = population_pool[bestIdx].genes;
    globalBestFE = population_pool[bestIdx].fitness;
    globalBestStructure = population_pool[bestIdx].structure;

    for (int gen = 0; gen < maxIterations; gen++) {
        vector<Chromosome> newPool;
        newPool.reserve(popSize);

        // Elitism: Preserve current best chromosome
        int curBest = 0;
        for (int i = 1; i < popSize; i++) {
            if (population_pool[i].fitness < population_pool[curBest].fitness) curBest = i;
        }
        newPool.push_back(population_pool[curBest]);

        while ((int)newPool.size() < popSize) {
            int p1Idx = tournamentSelect();
            int p2Idx = tournamentSelect();

            Chromosome child;
            if (randDouble(0.0, 1.0) < crossoverRate) {
                child.genes = crossover(population_pool[p1Idx].genes, population_pool[p2Idx].genes);
            } else {
                child.genes = population_pool[p1Idx].genes;
            }

            if (randDouble(0.0, 1.0) < mutationRate) {
                child.genes = mutate(child.genes);
            }

            child.genes = optimizeButterfly(child.genes);
            child.fitness = calculateFitness(child.genes, child.structure);
            newPool.push_back(child);

            if (child.fitness < globalBestFE) {
                globalBestFE = child.fitness;
                globalBest = child.genes;
                globalBestStructure = child.structure;
            }
        }

        population_pool = newPool;
        population.clear();
        for (int i = 0; i < popSize; i++) {
            population.push_back(population_pool[i].genes);
        }

        int iterBest = 0;
        for (int i = 1; i < popSize; i++) {
            if (population_pool[i].fitness < population_pool[iterBest].fitness) iterBest = i;
        }

        logIteration(gen + 1, population_pool[iterBest].fitness, population_pool[iterBest].structure);
    }

    hillClimbBest();
}