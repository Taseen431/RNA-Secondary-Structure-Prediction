/*
 * woa.cpp  --  Whale Optimization Algorithm (implementation)
 *
 * Canonical WOA over m-dimensional continuous position vectors, decoded into a
 * stem-index list by the same threshold + priority-ordering scheme PSO, DE and
 * GWO use, so the four differ only in how positions move.  Every whale is moved
 * and re-evaluated once per iteration, so the cost is pop evaluations per
 * iteration (same as GA, PSO, SA, HS, DE and GWO).  Whales move unconditionally
 * -- there is no greedy accept -- so the logged curve is non-monotone, as PSO's
 * and GWO's are.  The leader is a persistent best-so-far position, refreshed
 * whale by whale inside the pod loop exactly as GWO refreshes alpha.
 * A, C, l and the branch draw p are drawn once per whale and held across the
 * dimension loop, following the reference implementation.
 * "parameter" is p, the encircling / spiral branch probability.
 */

#include "woa.h"

// Continuous position -> stem-index list: keep the components above the
// threshold, order them by value (the ordering is the conflict priority used by
// buildStructure), and truncate to MAX_LEN.
Individual WOA::decode(const vector<double> &position) {
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

void WOA::run() {
    extractStems();
    int m = (int)stemList.size();
    int popSize = populationSize;
    const double twoPi = 2.0 * acos(-1.0);

    vector<vector<double>> positions(popSize, vector<double>(m));
    for (int i = 0; i < popSize; i++)
        for (int d = 0; d < m; d++) positions[i][d] = randDouble(0.0, 1.0);

    // the leader: the best position seen so far, the prey the pod encircles
    vector<double> leaderPos(m, 0.0);
    double leaderFE = numeric_limits<double>::max();

    population.clear();
    vector<double> fitVals(popSize);
    vector<string> structs(popSize);
    for (int i = 0; i < popSize; i++) {
        Individual b = decode(positions[i]);
        population.push_back(b);
        fitVals[i] = calculateFitness(b, structs[i]);
        if (fitVals[i] < leaderFE) { leaderFE = fitVals[i]; leaderPos = positions[i]; }
    }

    int bestIdx = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
    globalBest = population[bestIdx];
    globalBestFE = fitVals[bestIdx];
    globalBestStructure = structs[bestIdx];

    for (int iter = 0; iter < maxIterations; iter++) {
        // encircling coefficient, decayed linearly from 2 down to 0
        double a = 2.0 - 2.0 * (double)iter / (double)maxIterations;

        for (int i = 0; i < popSize; i++) {
            double r1 = randDouble(0.0, 1.0), r2 = randDouble(0.0, 1.0);
            double A = 2.0 * a * r1 - a, C = 2.0 * r2;
            double l = randDouble(-1.0, 1.0);   // spiral parameter
            double p = randDouble(0.0, 1.0);    // encircling / spiral branch draw
            // exploration branch picks any whale of the pod, self included
            int rw = randInt(0, popSize - 1);

            for (int d = 0; d < m; d++) {
                double np;
                if (p < parameter) {
                    if (fabs(A) < 1.0) {
                        // shrinking encircling: close in on the leader
                        np = leaderPos[d] - A * fabs(C * leaderPos[d] - positions[i][d]);
                    } else {
                        // search for prey: move relative to a random whale
                        np = positions[rw][d] - A * fabs(C * positions[rw][d] - positions[i][d]);
                    }
                } else {
                    // spiral: bubble-net attack around the leader
                    double dLeader = fabs(leaderPos[d] - positions[i][d]);
                    np = dLeader * exp(SPIRAL_B * l) * cos(twoPi * l) + leaderPos[d];
                }
                if (np < 0.0) np = 0.0; else if (np > 1.0) np = 1.0;
                positions[i][d] = np;
            }

            // the whale always moves: there is no greedy accept in WOA
            Individual b = decode(positions[i]);
            string st; double fe = calculateFitness(b, st);
            population[i] = b; fitVals[i] = fe; structs[i] = st;

            if (fe < leaderFE) { leaderFE = fe; leaderPos = positions[i]; }
            if (fe < globalBestFE) { globalBestFE = fe; globalBest = b; globalBestStructure = st; }
        }
        // log the best of the CURRENT pod at this iteration
        int curBest = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
        logIteration(iter + 1, fitVals[curBest], structs[curBest]);
    }
    hillClimbBest();
}
