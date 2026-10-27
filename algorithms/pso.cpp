/*
 * pso.cpp  --  Particle Swarm Optimization (implementation)
 *
 * A particle's continuous position vector is decoded into a stem-index list
 * (threshold + priority ordering), which the shared framework then scores.
 * PSO update logic is IDENTICAL to the original; only per-iteration logging
 * is added.  "parameter" is parsed but unused by PSO.
 */

#include "pso.h"

Individual PSO::decode(const vector<double> &position) {
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

void PSO::run() {
    extractStems();
    int m = (int)stemList.size();
    int popSize = populationSize;

    vector<vector<double>> positions(popSize, vector<double>(m));
    vector<vector<double>> velocities(popSize, vector<double>(m));
    for (int i = 0; i < popSize; i++)
        for (int d = 0; d < m; d++) { positions[i][d] = randDouble(0.0, 1.0); velocities[i][d] = randDouble(-VMAX, VMAX); }

    vector<vector<double>> pbestPos = positions;
    vector<double> pbestFE(popSize);
    vector<string> pbestStruct(popSize);

    population.clear();
    globalBestFE = numeric_limits<double>::max();
    vector<double> gbestPos(m, 0.0);
    bool gbestSet = false;

    for (int i = 0; i < popSize; i++) {
        Individual b = decode(positions[i]);
        population.push_back(b);
        string st; double fe = calculateFitness(b, st);
        pbestFE[i] = fe; pbestStruct[i] = st;
        if (fe < globalBestFE) { globalBestFE = fe; globalBest = b; globalBestStructure = st; gbestPos = positions[i]; gbestSet = true; }
    }
    if (!gbestSet) gbestPos = positions[0];

    for (int iter = 0; iter < maxIterations; iter++) {
        double iterBestFE = numeric_limits<double>::max();
        string iterBestStruct;
        for (int i = 0; i < popSize; i++) {
            for (int d = 0; d < m; d++) {
                double r1 = randDouble(0.0, 1.0), r2 = randDouble(0.0, 1.0);
                double vel = INERTIA_W * velocities[i][d]
                           + COGNITIVE_C1 * r1 * (pbestPos[i][d] - positions[i][d])
                           + SOCIAL_C2 * r2 * (gbestPos[d] - positions[i][d]);
                if (vel > VMAX) vel = VMAX; else if (vel < -VMAX) vel = -VMAX;
                velocities[i][d] = vel;
                double np = positions[i][d] + vel;
                if (np < 0.0) np = 0.0; else if (np > 1.0) np = 1.0;
                positions[i][d] = np;
            }
            Individual b = decode(positions[i]);
            string st; double fe = calculateFitness(b, st);
            if (fe < pbestFE[i]) { pbestFE[i] = fe; pbestPos[i] = positions[i]; pbestStruct[i] = st; }
            if (fe < globalBestFE) { globalBestFE = fe; globalBest = b; globalBestStructure = st; gbestPos = positions[i]; }
            if (fe < iterBestFE) { iterBestFE = fe; iterBestStruct = st; }  // current-iteration best
        }
        // log the best of the CURRENT swarm at this iteration
        logIteration(iter + 1, iterBestFE, iterBestStruct);
    }
    hillClimbBest();
}
