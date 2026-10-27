/*
 * boa.cpp  --  Butterfly Optimization Algorithm (implementation)
 *
 * Search operators and the BOA loop are IDENTICAL to the original program;
 * the only addition is one logIteration() call per iteration (logging only).
 * "parameter" is the BOA switch probability (5th input field).
 */

#include "boa.h"

pair<Individual, Individual> BOA::separateGlobalSearch(const Individual &B) {
    int n_ = (int)B.size();
    if (n_ == 0) return {B, B};
    Individual New_B1(n_), New_B2(n_);
    for (int p = 0; p < n_; p++) {
        int x = randInt(1, max(1, (int)stemList.size()));
        if ((p + 1) % 2 == 1) { New_B1[p] = B[p]; New_B2[p] = x; }
        else { New_B1[p] = x; New_B2[p] = B[p]; }
    }
    return {optimizeButterfly(New_B1), optimizeButterfly(New_B2)};
}

Individual BOA::reverseGlobalSearch(const Individual &B, const Individual &Bp) {
    int n_ = (int)B.size();
    if (n_ == 0) return B;
    Individual New_B(n_);
    for (int p = 0; p < n_; p++) {
        double x = randDouble(0.0, 1.0);
        if (x < 0.5) New_B[p] = B[p];
        else {
            int mirrorIdx = n_ - 1 - p;
            New_B[p] = (mirrorIdx >= 0 && mirrorIdx < (int)Bp.size()) ? Bp[mirrorIdx] : B[p];
        }
    }
    return optimizeButterfly(New_B);
}

void BOA::exchangeLocalSearch(const Individual &B1, const Individual &B2,
                              Individual &New_B1, Individual &New_B2,
                              Individual &New_B3, Individual &New_B4) {
    int len = min(B1.size(), B2.size());
    if (len == 0) { New_B1 = B1; New_B2 = B2; New_B3 = B1; New_B4 = B2; return; }
    int mid = len / 2;
    New_B1.assign(len, 0); New_B2.assign(len, 0); New_B3.assign(len, 0); New_B4.assign(len, 0);
    for (int p = 0; p < mid; p++) {
        New_B1[p] = B1[p]; New_B2[p] = B2[p]; New_B3[p] = B1[p];
        New_B4[p] = (mid + p) < (int)B1.size() ? B1[mid + p] : B1[p];
    }
    for (int p = mid; p < len; p++) {
        New_B1[p] = B2[p]; New_B2[p] = B1[p];
        New_B3[p] = (p - mid) < (int)B2.size() ? B2[p - mid] : B2[p];
        New_B4[p] = B2[p];
    }
    New_B1 = optimizeButterfly(New_B1); New_B2 = optimizeButterfly(New_B2);
    New_B3 = optimizeButterfly(New_B3); New_B4 = optimizeButterfly(New_B4);
}

pair<Individual, Individual> BOA::mergeLocalSearch(const Individual &B1, const Individual &B2) {
    int len = min(B1.size(), B2.size());
    if (len == 0) return {B1, B2};
    Individual New_B1(len), New_B2(len);
    for (int p = 0; p < len; p++) {
        if ((p + 1) % 2 == 1) { New_B1[p] = B1[p]; New_B2[p] = B2[p]; }
        else { New_B1[p] = B2[p]; New_B2[p] = B1[p]; }
    }
    return {optimizeButterfly(New_B1), optimizeButterfly(New_B2)};
}

void BOA::run() {
    extractStems();

    population.clear();
    for (int i = 0; i < populationSize; i++) population.push_back(randomPermutation());

    vector<double> fitVals(population.size());
    vector<string> structs(population.size());
    for (size_t i = 0; i < population.size(); i++) fitVals[i] = calculateFitness(population[i], structs[i]);

    int bestIdx = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
    globalBest = population[bestIdx];
    globalBestFE = fitVals[bestIdx];
    globalBestStructure = structs[bestIdx];

    for (int iter = 0; iter < maxIterations; iter++) {
        vector<Individual> newInds;
        double r = randDouble(0.0, 1.0);
        if (r < parameter) {
            double r1 = randDouble(0.0, 1.0);
            if (r1 < 0.5) { auto pr = separateGlobalSearch(globalBest); newInds.push_back(pr.first); newInds.push_back(pr.second); }
            else { int idx2 = randInt(0, (int)population.size() - 1); newInds.push_back(reverseGlobalSearch(globalBest, population[idx2])); }
        } else {
            int i1 = randInt(0, (int)population.size() - 1);
            int i2 = randInt(0, (int)population.size() - 1);
            double r1 = randDouble(0.0, 1.0);
            if (r1 < 0.5) {
                Individual n1, n2, n3, n4;
                exchangeLocalSearch(population[i1], population[i2], n1, n2, n3, n4);
                newInds.push_back(n1); newInds.push_back(n2); newInds.push_back(n3); newInds.push_back(n4);
            } else { auto pr = mergeLocalSearch(population[i1], population[i2]); newInds.push_back(pr.first); newInds.push_back(pr.second); }
        }
        for (auto &nb : newInds) {
            if (randDouble(0.0, 1.0) < 0.2) {
                int m_ = (int)stemList.size();
                int op = randInt(0, 2);
                int cap = MAX_LEN;
                if (op == 0 && !nb.empty()) nb[randInt(0, (int)nb.size() - 1)] = randInt(1, m_);
                else if (op == 1) { if ((int)nb.size() < cap) nb.push_back(randInt(1, m_)); }
                else { if (nb.size() > 1) nb.erase(nb.begin() + randInt(0, (int)nb.size() - 1)); }
                nb = optimizeButterfly(nb);
            }
            string st; double fe = calculateFitness(nb, st);
            int worstIdx = 0; double worstVal = fitVals[0];
            for (size_t k = 1; k < fitVals.size(); k++) if (fitVals[k] > worstVal) { worstVal = fitVals[k]; worstIdx = (int)k; }
            if (fe < worstVal) { population[worstIdx] = nb; fitVals[worstIdx] = fe; structs[worstIdx] = st; }
            if (fe < globalBestFE) { globalBestFE = fe; globalBest = nb; globalBestStructure = st; }
        }
        // log the best of the CURRENT population at this iteration
        int curBest = (int)(min_element(fitVals.begin(), fitVals.end()) - fitVals.begin());
        logIteration(iter + 1, fitVals[curBest], structs[curBest]);
    }
    hillClimbBest();
}
