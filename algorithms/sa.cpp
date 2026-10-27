/*
 * sa.cpp  --  Simulated Annealing (implementation)
 *
 * populationSize independent Metropolis chains share one geometric cooling
 * schedule; one trial move per chain per iteration, so the cost is pop
 * evaluations per iteration (same as GA and PSO).  Single-chain SA would give
 * 1 evaluation per iteration and need ~45000 iterations to match that budget.
 * The start temperature is derived from the spread of the initial energies.
 * "parameter" is the cooling rate alpha.
 */

#include "sa.h"

// Three-way neighbour move (replace / append / delete), capped at MAX_LEN.
Individual SA::neighbour(const Individual &ind) {
    Individual c = ind;
    int m = (int)stemList.size();
    int op = randInt(0, 2);
    int cap = MAX_LEN;
    if (op == 0 && !c.empty()) c[randInt(0, (int)c.size() - 1)] = randInt(1, m);
    else if (op == 1) { if ((int)c.size() < cap) c.push_back(randInt(1, m)); }
    else { if (c.size() > 1) c.erase(c.begin() + randInt(0, (int)c.size() - 1)); }
    return optimizeButterfly(c);
}

void SA::run() {
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

    // Start temperature from the energy spread of the initial chains: a fixed T0
    // would accept everything or nothing, since energies span roughly -50..+100.
    double eMax = *max_element(fitVals.begin(), fitVals.end());
    double T0 = eMax - fitVals[bestIdx];
    if (T0 < T0_MIN) T0 = T0_MIN;
    double T = T0;

    for (int iter = 0; iter < maxIterations; iter++) {
        T *= parameter;                 // geometric cooling
        if (T < T_MIN) T = T_MIN;

        for (int i = 0; i < popSize; i++) {
            Individual cand = neighbour(population[i]);
            string cst; double cf = calculateFitness(cand, cst);
            double delta = cf - fitVals[i];
            // Metropolis: downhill always, uphill with probability exp(-delta/T)
            bool accept = (delta <= 0.0) || (randDouble(0.0, 1.0) < exp(-delta / T));
            if (accept) { population[i] = cand; fitVals[i] = cf; structs[i] = cst; }
            if (cf < globalBestFE) { globalBestFE = cf; globalBest = cand; globalBestStructure = cst; }
        }

        // log the best of the CURRENT chain states at this iteration
        int curBest = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
        logIteration(iter + 1, fitVals[curBest], structs[curBest]);
    }
    hillClimbBest();
}
