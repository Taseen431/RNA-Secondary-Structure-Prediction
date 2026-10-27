/*
 * fa.cpp  --  Firefly Algorithm (implementation)
 *
 * Canonical FA over m-dimensional continuous position vectors, decoded into a
 * stem-index list by the same threshold + priority-ordering scheme PSO, DE, GWO
 * and WOA use, so the five differ only in how positions move.  Attraction is
 * pairwise: every firefly is drawn towards each brighter firefly in turn and is
 * re-evaluated after each such move (Yang's formulation -- "move firefly i
 * towards j, evaluate new solutions and update light intensity" sits inside the
 * inner loop), so the cost is quadratic, ~pop*(pop-1)/2 evaluations per
 * iteration against pop for every other optimizer here.  FA is therefore run at
 * a reduced population to hold the evaluation budget equal; see
 * BENCHMARK_PROTOCOL.md 7.4.  Fireflies move unconditionally -- there is no
 * greedy accept -- so the logged curve is non-monotone, as PSO's, GWO's and
 * WOA's are.  Brightness needs no transform: j is brighter than i exactly when
 * its free energy is lower.
 * "parameter" is gamma, the light absorption coefficient.
 *
 * One scaling convention: the squared separation r^2 is the MEAN squared
 * component difference, not the raw sum.  m ranges from 68 to 5151 stems across
 * the twenty datasets, and a raw sum would drive exp(-gamma*r^2) to zero on every
 * long sequence -- attraction would vanish and FA would degenerate into a random
 * walk.  Normalising keeps r^2 in [0,1], so a single gamma has the same meaning
 * on every sequence.
 */

#include "fa.h"

// Continuous position -> stem-index list: keep the components above the
// threshold, order them by value (the ordering is the conflict priority used by
// buildStructure), and truncate to MAX_LEN.
Individual FA::decode(const vector<double> &position) {
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

void FA::run() {
    extractStems();
    int m = (int)stemList.size();
    int popSize = populationSize;

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

    for (int iter = 0; iter < maxIterations; iter++) {
        // randomization strength, damped geometrically as the swarm settles
        double alpha = ALPHA0 * pow(ALPHA_DAMP, (double)iter);

        for (int i = 0; i < popSize; i++) {
            bool moved = false;
            for (int j = 0; j < popSize; j++) {
                // brightness is the free energy itself: lower is brighter
                if (j == i || fitVals[j] >= fitVals[i]) continue;

                double sum = 0.0;
                for (int d = 0; d < m; d++) {
                    double diff = positions[i][d] - positions[j][d];
                    sum += diff * diff;
                }
                double r2 = sum / (double)m;              // mean squared separation
                double beta = BETA0 * exp(-parameter * r2);

                for (int d = 0; d < m; d++) {
                    double np = positions[i][d] + beta * (positions[j][d] - positions[i][d])
                                + alpha * (randDouble(0.0, 1.0) - 0.5);
                    if (np < 0.0) np = 0.0; else if (np > 1.0) np = 1.0;
                    positions[i][d] = np;
                }
                moved = true;

                // the firefly always moves: there is no greedy accept in FA, and
                // its light intensity is refreshed before the next partner is tried
                Individual b = decode(positions[i]);
                string st; double fe = calculateFitness(b, st);
                population[i] = b; fitVals[i] = fe; structs[i] = st;

                if (fe < globalBestFE) { globalBestFE = fe; globalBest = b; globalBestStructure = st; }
            }

            if (!moved) {
                // no brighter partner exists: the brightest firefly moves at random
                for (int d = 0; d < m; d++) {
                    double np = positions[i][d] + alpha * (randDouble(0.0, 1.0) - 0.5);
                    if (np < 0.0) np = 0.0; else if (np > 1.0) np = 1.0;
                    positions[i][d] = np;
                }
                Individual b = decode(positions[i]);
                string st; double fe = calculateFitness(b, st);
                population[i] = b; fitVals[i] = fe; structs[i] = st;

                if (fe < globalBestFE) { globalBestFE = fe; globalBest = b; globalBestStructure = st; }
            }
        }
        // log the best of the CURRENT swarm at this iteration
        int curBest = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
        logIteration(iter + 1, fitVals[curBest], structs[curBest]);
    }
    hillClimbBest();
}
