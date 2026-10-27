/*
 * hs.cpp  --  Harmony Search (implementation)
 *
 * populationSize harmonies are improvised per iteration (the documented
 * parallel-HS variant), so the cost is pop evaluations per iteration, the same
 * as GA, PSO and SA.  Canonical HS improvises ONE harmony per iteration, which
 * would reproduce BOA's ~1 evaluation/iteration budget problem.
 * Each new harmony replaces the worst member of the harmony memory if better,
 * so the logged curve is monotone.
 * "parameter" is HMCR, the harmony memory considering rate.
 */

#include "hs.h"

// Improvise one harmony gene by gene: memory consideration (HMCR) with pitch
// adjustment (PAR), otherwise a random stem index.  Length is borrowed from a
// random memory member and capped at MAX_LEN.
Individual HS::improvise() {
    int m = (int)stemList.size();
    int hms = (int)population.size();
    int cap = MAX_LEN;

    int len = (int)population[randInt(0, hms - 1)].size();
    if (len > cap) len = cap;
    if (len < 1) len = 1;

    Individual nb(len);
    for (int p = 0; p < len; p++) {
        if (randDouble(0.0, 1.0) < parameter) {
            const Individual &src = population[randInt(0, hms - 1)];
            int val;
            if (src.empty()) val = randInt(1, m);
            else val = (p < (int)src.size()) ? src[p] : src[randInt(0, (int)src.size() - 1)];
            // pitch adjustment: stemList is end-index-major, so a +/-1 shift of
            // the index is a genuine local move in the search space
            if (randDouble(0.0, 1.0) < PAR) {
                val += (randInt(0, 1) == 0) ? -PITCH_BW : PITCH_BW;
                if (val < 1) val = 1;
                if (val > m) val = m;
            }
            nb[p] = val;
        } else {
            nb[p] = randInt(1, m);
        }
    }
    return optimizeButterfly(nb);
}

void HS::run() {
    extractStems();
    int hms = populationSize;

    population.clear();
    for (int i = 0; i < hms; i++) population.push_back(randomPermutation());

    vector<double> fitVals(hms);
    vector<string> structs(hms);
    for (int i = 0; i < hms; i++) fitVals[i] = calculateFitness(population[i], structs[i]);

    int bestIdx = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
    globalBest = population[bestIdx];
    globalBestFE = fitVals[bestIdx];
    globalBestStructure = structs[bestIdx];

    for (int iter = 0; iter < maxIterations; iter++) {
        for (int i = 0; i < hms; i++) {
            Individual nb = improvise();
            string st; double fe = calculateFitness(nb, st);
            // memory update: a new harmony displaces the worst one if better
            int worstIdx = 0; double worstVal = fitVals[0];
            for (int k = 1; k < hms; k++) if (fitVals[k] > worstVal) { worstVal = fitVals[k]; worstIdx = k; }
            if (fe < worstVal) { population[worstIdx] = nb; fitVals[worstIdx] = fe; structs[worstIdx] = st; }
            if (fe < globalBestFE) { globalBestFE = fe; globalBest = nb; globalBestStructure = st; }
        }
        // log the best of the CURRENT harmony memory at this iteration
        int curBest = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
        logIteration(iter + 1, fitVals[curBest], structs[curBest]);
    }
    hillClimbBest();
}
