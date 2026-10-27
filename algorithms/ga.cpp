/*
 * ga.cpp  --  Genetic Algorithm (implementation)
 *
 * Tournament selection, one-point crossover, mutation and elitism are
 * IDENTICAL to the original; only per-iteration logging is added.
 * "parameter" is parsed but unused by GA.
 */

#include "ga.h"

int GA::tournamentSelect(const vector<double> &fitVals) {
    int bestIdx = randInt(0, (int)fitVals.size() - 1);
    double bestVal = fitVals[bestIdx];
    for (int t = 1; t < TOURNAMENT_SIZE; t++) {
        int idx = randInt(0, (int)fitVals.size() - 1);
        if (fitVals[idx] < bestVal) { bestVal = fitVals[idx]; bestIdx = idx; }
    }
    return bestIdx;
}

Individual GA::crossover(const Individual &a, const Individual &b) {
    if (a.empty()) return b;
    if (b.empty()) return a;
    int c1 = randInt(0, (int)a.size());
    int c2 = randInt(0, (int)b.size());
    Individual child(a.begin(), a.begin() + c1);
    for (int i = c2; i < (int)b.size(); i++) child.push_back(b[i]);
    if (child.empty()) child.push_back(randInt(1, (int)stemList.size()));
    return optimizeButterfly(child);
}

Individual GA::mutate(const Individual &ind) {
    Individual c = ind;
    int m = (int)stemList.size();
    int op = randInt(0, 2);
    int cap = MAX_LEN;
    if (op == 0 && !c.empty()) c[randInt(0, (int)c.size() - 1)] = randInt(1, m);
    else if (op == 1) { if ((int)c.size() < cap) c.push_back(randInt(1, m)); }
    else { if (c.size() > 1) c.erase(c.begin() + randInt(0, (int)c.size() - 1)); }
    return optimizeButterfly(c);
}

void GA::run() {
    extractStems();
    int popSize = populationSize;

    population.clear();
    for (int i = 0; i < popSize; i++) population.push_back(randomPermutation());

    vector<double> fitVals(popSize);
    vector<string> structs(popSize);
    for (int i = 0; i < popSize; i++) fitVals[i] = calculateFitness(population[i], structs[i]);

    int bestIdx = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
    globalBest = population[bestIdx];
    globalBestFE = fitVals[bestIdx];
    globalBestStructure = structs[bestIdx];

    for (int gen = 0; gen < maxIterations; gen++) {
        vector<Individual> newPop;
        newPop.reserve(popSize);

        int curBest = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
        newPop.push_back(population[curBest]); // elitism

        while ((int)newPop.size() < popSize) {
            const Individual &p1 = population[tournamentSelect(fitVals)];
            const Individual &p2 = population[tournamentSelect(fitVals)];
            Individual child;
            if (randDouble(0.0, 1.0) < CROSSOVER_RATE) child = crossover(p1, p2);
            else child = p1;
            if (randDouble(0.0, 1.0) < MUTATION_RATE) child = mutate(child);
            child = optimizeButterfly(child);
            newPop.push_back(child);
        }

        population = newPop;

        for (int i = 0; i < popSize; i++) {
            fitVals[i] = calculateFitness(population[i], structs[i]);
            if (fitVals[i] < globalBestFE) { globalBestFE = fitVals[i]; globalBest = population[i]; globalBestStructure = structs[i]; }
        }
        // log the best of the CURRENT population at this generation
        int iterBest = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
        logIteration(gen + 1, fitVals[iterBest], structs[iterBest]);
    }
    hillClimbBest();
}
