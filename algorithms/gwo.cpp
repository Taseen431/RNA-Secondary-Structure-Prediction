/*
 * gwo.cpp  --  Grey Wolf Optimizer (implementation)
 *
 * Canonical GWO over m-dimensional continuous position vectors, decoded into a
 * stem-index list by the same threshold + priority-ordering scheme PSO and DE
 * use, so the three differ only in how positions move.  Every wolf is moved and
 * re-evaluated once per iteration, so the cost is pop evaluations per iteration
 * (same as GA, PSO, SA, HS and DE).  Wolves move unconditionally -- there is no
 * greedy accept -- so the logged curve is non-monotone, as PSO's is.
 * The alpha / beta / delta leaders are a persistent best-so-far top-3 archive,
 * refreshed wolf by wolf inside the pack loop exactly as PSO refreshes gbest.
 * "parameter" is parsed but unused by GWO.
 */

#include "gwo.h"

// Continuous position -> stem-index list: keep the components above the
// threshold, order them by value (the ordering is the conflict priority used by
// buildStructure), and truncate to MAX_LEN.
Individual GWO::decode(const vector<double> &position) {
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

void GWO::run() {
    extractStems();
    int m = (int)stemList.size();
    int popSize = populationSize;

    vector<vector<double>> positions(popSize, vector<double>(m));
    for (int i = 0; i < popSize; i++)
        for (int d = 0; d < m; d++) positions[i][d] = randDouble(0.0, 1.0);

    // social hierarchy: the three best positions seen so far
    vector<double> alphaPos(m, 0.0), betaPos(m, 0.0), deltaPos(m, 0.0);
    double alphaFE = numeric_limits<double>::max();
    double betaFE  = numeric_limits<double>::max();
    double deltaFE = numeric_limits<double>::max();

    population.clear();
    vector<double> fitVals(popSize);
    vector<string> structs(popSize);
    for (int i = 0; i < popSize; i++) {
        Individual b = decode(positions[i]);
        population.push_back(b);
        fitVals[i] = calculateFitness(b, structs[i]);
        // rank into the hierarchy: alpha displaces beta, beta displaces delta
        if (fitVals[i] < alphaFE) {
            deltaFE = betaFE;  deltaPos = betaPos;
            betaFE  = alphaFE; betaPos  = alphaPos;
            alphaFE = fitVals[i]; alphaPos = positions[i];
        } else if (fitVals[i] < betaFE) {
            deltaFE = betaFE; deltaPos = betaPos;
            betaFE  = fitVals[i]; betaPos = positions[i];
        } else if (fitVals[i] < deltaFE) {
            deltaFE = fitVals[i]; deltaPos = positions[i];
        }
    }

    int bestIdx = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
    globalBest = population[bestIdx];
    globalBestFE = fitVals[bestIdx];
    globalBestStructure = structs[bestIdx];

    for (int iter = 0; iter < maxIterations; iter++) {
        // exploration coefficient, decayed linearly from 2 down to 0
        double a = 2.0 - 2.0 * (double)iter / (double)maxIterations;

        for (int i = 0; i < popSize; i++) {
            for (int d = 0; d < m; d++) {
                double r1 = randDouble(0.0, 1.0), r2 = randDouble(0.0, 1.0);
                double A1 = 2.0 * a * r1 - a, C1 = 2.0 * r2;
                double X1 = alphaPos[d] - A1 * fabs(C1 * alphaPos[d] - positions[i][d]);

                r1 = randDouble(0.0, 1.0); r2 = randDouble(0.0, 1.0);
                double A2 = 2.0 * a * r1 - a, C2 = 2.0 * r2;
                double X2 = betaPos[d] - A2 * fabs(C2 * betaPos[d] - positions[i][d]);

                r1 = randDouble(0.0, 1.0); r2 = randDouble(0.0, 1.0);
                double A3 = 2.0 * a * r1 - a, C3 = 2.0 * r2;
                double X3 = deltaPos[d] - A3 * fabs(C3 * deltaPos[d] - positions[i][d]);

                double np = (X1 + X2 + X3) / 3.0;
                if (np < 0.0) np = 0.0; else if (np > 1.0) np = 1.0;
                positions[i][d] = np;
            }

            // the wolf always moves: there is no greedy accept in GWO
            Individual b = decode(positions[i]);
            string st; double fe = calculateFitness(b, st);
            population[i] = b; fitVals[i] = fe; structs[i] = st;

            if (fe < alphaFE) {
                deltaFE = betaFE;  deltaPos = betaPos;
                betaFE  = alphaFE; betaPos  = alphaPos;
                alphaFE = fe; alphaPos = positions[i];
            } else if (fe < betaFE) {
                deltaFE = betaFE; deltaPos = betaPos;
                betaFE  = fe; betaPos = positions[i];
            } else if (fe < deltaFE) {
                deltaFE = fe; deltaPos = positions[i];
            }

            if (fe < globalBestFE) { globalBestFE = fe; globalBest = b; globalBestStructure = st; }
        }
        // log the best of the CURRENT pack at this iteration
        int curBest = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
        logIteration(iter + 1, fitVals[curBest], structs[curBest]);
    }
    hillClimbBest();
}
