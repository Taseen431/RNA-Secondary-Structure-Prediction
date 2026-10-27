/*
 * Greywolf.cpp  --  Grey Wolf Optimizer (implementation)
 *
 * Formulated for RNA secondary structure prediction with explicit Wolf class.
 * Inherits from RNA_Core and follows OPTIMIZER_CONTRACT.md.
 */

#include "Greywolf.h"

static Individual updateWolfPosition(const Individual &current, const Individual &alpha, const Individual &beta, const Individual &delta, int m, double aFactor, RNA_Core *core) {
    Individual cand = current;
    int nStems = (int)cand.size();

    // Guided by Alpha wolf
    if (randDouble(0.0, 1.0) < (aFactor / 2.0) && !alpha.empty()) {
        if (nStems > 0) cand[randInt(0, nStems - 1)] = alpha[randInt(0, (int)alpha.size() - 1)];
        else cand.push_back(alpha[0]);
    }

    // Guided by Beta wolf
    if (randDouble(0.0, 1.0) < (aFactor / 2.0) && !beta.empty()) {
        if ((int)cand.size() < 15) cand.push_back(beta[randInt(0, (int)beta.size() - 1)]);
        else if (!cand.empty()) cand[randInt(0, (int)cand.size() - 1)] = beta[randInt(0, (int)beta.size() - 1)];
    }

    // Guided by Delta wolf
    if (randDouble(0.0, 1.0) < (aFactor / 2.0) && !delta.empty()) {
        if ((int)cand.size() < 15) cand.push_back(delta[randInt(0, (int)delta.size() - 1)]);
    }

    if (randDouble(0.0, 1.0) < 0.10 && (int)cand.size() >= 2) {
        swap(cand[randInt(0, (int)cand.size() - 1)], cand[randInt(0, (int)cand.size() - 1)]);
    }

    return core->optimizeButterfly(cand);
}

void Greywolf::run() {
    extractStems();
    int m = (int)stemList.size();
    if (m == 0) return;

    int nWolves = populationSize;
    if (nWolves < 4) nWolves = 4;

    wolves.clear();
    wolves.resize(nWolves);
    population.clear();

    for (int i = 0; i < nWolves; i++) {
        wolves[i].position = randomPermutation();
        wolves[i].fitness = calculateFitness(wolves[i].position, wolves[i].structure);
        population.push_back(wolves[i].position);
    }

    vector<int> idxs(nWolves);
    for (int i = 0; i < nWolves; i++) idxs[i] = i;
    sort(idxs.begin(), idxs.end(), [this](int a, int b) {
        return wolves[a].fitness < wolves[b].fitness;
    });

    alphaPos = wolves[idxs[0]].position; alphaFE = wolves[idxs[0]].fitness;
    betaPos = wolves[idxs[1]].position;  betaFE = wolves[idxs[1]].fitness;
    deltaPos = wolves[idxs[2]].position; deltaFE = wolves[idxs[2]].fitness;

    globalBest = alphaPos;
    globalBestFE = alphaFE;
    globalBestStructure = wolves[idxs[0]].structure;

    for (int iter = 0; iter < maxIterations; iter++) {
        double aFactor = 2.0 - 2.0 * ((double)iter / maxIterations);

        for (int i = 0; i < nWolves; i++) {
            Individual cand = updateWolfPosition(wolves[i].position, alphaPos, betaPos, deltaPos, m, aFactor, this);
            string cst;
            double cf = calculateFitness(cand, cst);

            if (cf <= wolves[i].fitness) {
                wolves[i].position = cand;
                wolves[i].fitness = cf;
                wolves[i].structure = cst;
            }
        }

        for (int i = 0; i < nWolves; i++) idxs[i] = i;
        sort(idxs.begin(), idxs.end(), [this](int a, int b) {
            return wolves[a].fitness < wolves[b].fitness;
        });

        alphaPos = wolves[idxs[0]].position; alphaFE = wolves[idxs[0]].fitness;
        betaPos = wolves[idxs[1]].position;  betaFE = wolves[idxs[1]].fitness;
        deltaPos = wolves[idxs[2]].position; deltaFE = wolves[idxs[2]].fitness;

        if (alphaFE < globalBestFE) {
            globalBestFE = alphaFE;
            globalBest = alphaPos;
            globalBestStructure = wolves[idxs[0]].structure;
        }

        population.clear();
        for (int i = 0; i < nWolves; i++) {
            population.push_back(wolves[i].position);
        }

        logIteration(iter + 1, alphaFE, wolves[idxs[0]].structure);
    }

    hillClimbBest();
}
