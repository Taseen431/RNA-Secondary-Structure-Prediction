/*
 * de.cpp  --  Differential Evolution (implementation)
 *
 * DE/rand/1/bin over m-dimensional continuous position vectors, decoded into a
 * stem-index list by the same threshold + priority-ordering scheme PSO uses, so
 * DE and PSO differ only in how positions move.  One trial vector per target per
 * iteration, so the cost is pop evaluations per iteration (same as GA, PSO, SA
 * and HS).  Greedy selection never accepts a worse member, so the logged curve
 * is monotone.
 * "parameter" is F, the differential weight.
 */

#include "de.h"

// Continuous position -> stem-index list: keep the components above the
// threshold, order them by value (the ordering is the conflict priority used by
// buildStructure), and truncate to MAX_LEN.
Individual DE::decode(const vector<double> &position) {
    int m = (int)stemList.size();
    vector<pair<int, double>> selected;
    for (int i = 0; i < m; i++)
        if (position[i] > SELECT_THRESHOLD) selected.push_back({i + 1, position[i]});
    sort(selected.begin(), selected.end(),
         [](const pair<int, double> &a, const pair<int, double> &b) { return a.second > b.second; });
    if (selected.empty()) {
        int bestI = 0;
        for (int i = 1; i < m; i++) if (position[i] > position[bestI]) bestI = i;
        selected.push_back({bestI + 1, position[bestI]});
    }
    Individual ind;
    int cap = MAX_LEN;
    for (int i = 0; i < (int)selected.size() && i < cap; i++) ind.push_back(selected[i].first);
    return optimizeButterfly(ind);
}

void DE::run() {
    extractStems();
    int m = (int)stemList.size();
    // DE/rand/1 draws three donors distinct from the target, so the population
    // cannot be smaller than four.
    int popSize = populationSize < 4 ? 4 : populationSize;

    vector<vector<double>> positions(popSize, vector<double>(m));
    for (int i = 0; i < popSize; i++)
        for (int d = 0; d < m; d++) positions[i][d] = randDouble(0.0, 1.0);

    population.clear();
    vector<double> fitVals(popSize);
    vector<string> structs(popSize);
    for (int i = 0; i < popSize; i++) {
        Individual b = decode(positions[i]);
        population.push_back(b);
        fitVals[i] = calculateFitness(b, structs[i]);
    }

    int bestIdx = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
    globalBest = population[bestIdx];
    globalBestFE = fitVals[bestIdx];
    globalBestStructure = structs[bestIdx];

    vector<double> trial(m);
    for (int iter = 0; iter < maxIterations; iter++) {
        for (int i = 0; i < popSize; i++) {
            int r1, r2, r3;
            do { r1 = randInt(0, popSize - 1); } while (r1 == i);
            do { r2 = randInt(0, popSize - 1); } while (r2 == i || r2 == r1);
            do { r3 = randInt(0, popSize - 1); } while (r3 == i || r3 == r1 || r3 == r2);

            // mutation v = x_r1 + F * (x_r2 - x_r3), then binomial crossover;
            // jrand forces at least one component to come from the mutant
            int jrand = randInt(0, m - 1);
            for (int d = 0; d < m; d++) {
                if (d == jrand || randDouble(0.0, 1.0) < CROSSOVER_CR) {
                    double v = positions[r1][d] + parameter * (positions[r2][d] - positions[r3][d]);
                    if (v < 0.0) v = 0.0; else if (v > 1.0) v = 1.0;
                    trial[d] = v;
                } else {
                    trial[d] = positions[i][d];
                }
            }

            Individual b = decode(trial);
            string st; double fe = calculateFitness(b, st);
            // greedy selection: the trial replaces its target when not worse
            if (fe <= fitVals[i]) { positions[i] = trial; population[i] = b; fitVals[i] = fe; structs[i] = st; }
            if (fe < globalBestFE) { globalBestFE = fe; globalBest = b; globalBestStructure = st; }
        }
        // log the best of the CURRENT population at this iteration
        int curBest = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
        logIteration(iter + 1, fitVals[curBest], structs[curBest]);
    }
    hillClimbBest();
}
