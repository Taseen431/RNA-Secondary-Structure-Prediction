/*
 * abco.cpp  --  Artificial Bee Colony Optimization (implementation)
 *
 * Employed / onlooker / scout bee phases are IDENTICAL to the original;
 * only per-iteration logging is added.
 */

#include "abco.h"

// Objective (free energy, lower better) -> ABC fitness (higher better).
double ABCO::computeFit(double f) {
    if (f >= 0.0) return 1.0 / (1.0 + f);
    return 1.0 + fabs(f);
}

// Neighbour of food source i, influenced by a random partner k.
Individual ABCO::generateNeighbour(int i) {
    Individual cand = population[i];
    int m = (int)stemList.size();
    int popSize = (int)population.size();

    int k = i;
    if (popSize > 1) { do { k = randInt(0, popSize - 1); } while (k == i); }
    const Individual &partner = population[k];

    int cap = MAX_LEN;
    int op = randInt(0, 2);
    if (op == 0 && !cand.empty()) {
        int val = partner.empty() ? randInt(1, m) : partner[randInt(0, (int)partner.size() - 1)];
        cand[randInt(0, (int)cand.size() - 1)] = val;
    } else if (op == 1) {
        if ((int)cand.size() < cap) {
            int val = partner.empty() ? randInt(1, m) : partner[randInt(0, (int)partner.size() - 1)];
            cand.push_back(val);
        }
    } else {
        if (cand.size() > 1) cand.erase(cand.begin() + randInt(0, (int)cand.size() - 1));
    }
    return optimizeButterfly(cand);
}

void ABCO::run() {
    extractStems();
    int popSize = populationSize;

    population.clear();
    for (int i = 0; i < popSize; i++) population.push_back(randomPermutation());

    vector<double> fVals(popSize);
    vector<string> structs(popSize);
    vector<int> trials(popSize, 0);
    for (int i = 0; i < popSize; i++) fVals[i] = calculateFitness(population[i], structs[i]);

    int bestIdx = (int)(min_element(fVals.begin(), fVals.end()) - fVals.begin());
    globalBest = population[bestIdx];
    globalBestFE = fVals[bestIdx];
    globalBestStructure = structs[bestIdx];

    for (int cycle = 0; cycle < maxIterations; cycle++) {
        // ----- Employed bees -----
        for (int i = 0; i < popSize; i++) {
            Individual cand = generateNeighbour(i);
            string cst; double cf = calculateFitness(cand, cst);
            if (cf < fVals[i]) { population[i] = cand; fVals[i] = cf; structs[i] = cst; trials[i] = 0; }
            else trials[i]++;
            if (cf < globalBestFE) { globalBestFE = cf; globalBest = cand; globalBestStructure = cst; }
        }

        // ----- Probabilities -----
        vector<double> fit(popSize);
        double sumFit = 0.0;
        for (int i = 0; i < popSize; i++) { fit[i] = computeFit(fVals[i]); sumFit += fit[i]; }
        vector<double> prob(popSize);
        for (int i = 0; i < popSize; i++) prob[i] = (sumFit > 0.0) ? fit[i] / sumFit : 1.0 / popSize;

        // ----- Onlooker bees -----
        int idx = 0, count = 0;
        while (count < popSize) {
            if (randDouble(0.0, 1.0) < prob[idx]) {
                count++;
                Individual cand = generateNeighbour(idx);
                string cst; double cf = calculateFitness(cand, cst);
                if (cf < fVals[idx]) { population[idx] = cand; fVals[idx] = cf; structs[idx] = cst; trials[idx] = 0; }
                else trials[idx]++;
                if (cf < globalBestFE) { globalBestFE = cf; globalBest = cand; globalBestStructure = cst; }
            }
            idx = (idx + 1) % popSize;
        }

        // ----- Scout bees -----
        int maxTrialIdx = 0;
        for (int i = 1; i < popSize; i++) if (trials[i] > trials[maxTrialIdx]) maxTrialIdx = i;
        int limit = (int)(parameter * popSize);
        if (limit < 5) limit = 5; // safety lower bound
        if (trials[maxTrialIdx] > limit) {
            population[maxTrialIdx] = randomPermutation();
            fVals[maxTrialIdx] = calculateFitness(population[maxTrialIdx], structs[maxTrialIdx]);
            trials[maxTrialIdx] = 0;
            if (fVals[maxTrialIdx] < globalBestFE) { globalBestFE = fVals[maxTrialIdx]; globalBest = population[maxTrialIdx]; globalBestStructure = structs[maxTrialIdx]; }
        }

        // log the best of the CURRENT population at this cycle
        int curBest = (int)(min_element(fVals.begin(), fVals.end()) - fVals.begin());
        logIteration(cycle + 1, fVals[curBest], structs[curBest]);
    }
    hillClimbBest();
}
