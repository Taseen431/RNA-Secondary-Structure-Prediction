/*
 * ABCSS.cpp  --  Artificial Bee Colony Selective Search (implementation)
 *
 * Formulated for RNA secondary structure prediction with explicit Solution class.
 * Inherits from RNA_Core and follows OPTIMIZER_CONTRACT.md.
 */

#include "ABCSS.h"

double ABCSS::computeFit(double f) {
    if (f >= 0.0) return 1.0 / (1.0 + f);
    return 1.0 + fabs(f);
}

Individual ABCSS::applyRuleUpdate(const Individual &ind, const Individual &partner, int ruleType) {
    Individual cand = ind;
    int m = (int)stemList.size();
    if (m == 0) return cand;

    int nStems = (int)cand.size();

    switch (ruleType % 5) {
    case 0:
        if (nStems >= 2) {
            int i1 = randInt(0, nStems - 1);
            int i2 = randInt(0, nStems - 1);
            swap(cand[i1], cand[i2]);
        }
        break;
    case 1:
        if (nStems > 0 && !partner.empty()) {
            int i1 = randInt(0, nStems - 1);
            int pVal = partner[randInt(0, (int)partner.size() - 1)];
            cand[i1] = pVal;
        }
        break;
    case 2:
        if (!partner.empty()) {
            int cut = randInt(1, max(1, nStems));
            Individual newCand(cand.begin(), cand.begin() + min(cut, nStems));
            for (int val : partner) {
                if (find(newCand.begin(), newCand.end(), val) == newCand.end()) {
                    newCand.push_back(val);
                }
                if ((int)newCand.size() >= MAX_LEN) break;
            }
            cand = newCand;
        }
        break;
    case 3:
        if (nStems < MAX_LEN && !partner.empty()) {
            int pVal = partner[randInt(0, (int)partner.size() - 1)];
            cand.push_back(pVal);
        }
        break;
    case 4:
        if (nStems >= 3) {
            int l = randInt(0, nStems - 2);
            int r = randInt(l + 1, nStems - 1);
            reverse(cand.begin() + l, cand.begin() + r + 1);
        }
        break;
    }

    return optimizeButterfly(cand);
}

Individual ABCSS::kOptSearch(const Individual &ind, int maxTrial) {
    Individual bestCand = ind;
    string bestStr;
    double bestFE = calculateFitness(bestCand, bestStr);

    for (int t = 0; t < maxTrial; t++) {
        Individual trial = bestCand;
        if ((int)trial.size() >= 2) {
            int i1 = randInt(0, (int)trial.size() - 1);
            int i2 = randInt(0, (int)trial.size() - 1);
            swap(trial[i1], trial[i2]);
        } else {
            trial.push_back(randInt(1, (int)stemList.size()));
        }
        trial = optimizeButterfly(trial);
        string st;
        double fe = calculateFitness(trial, st);
        if (fe < bestFE) {
            bestFE = fe;
            bestCand = trial;
        }
    }
    return bestCand;
}

void ABCSS::run() {
    extractStems();
    int popSize = populationSize;
    if (popSize < 4) popSize = 4;

    int limit = (int)(parameter * popSize);
    if (limit < 5) limit = 5;
    updateTryLimit = limit;

    particle.clear();
    particle.resize(popSize);

    population.clear();

    for (int i = 0; i < popSize; i++) {
        particle[i].Xid = randomPermutation();
        particle[i].Pid = particle[i].Xid;
        particle[i].fitness = calculateFitness(particle[i].Xid, particle[i].structure);
        particle[i].Pid_val = particle[i].fitness;
        particle[i].trial = 0;
        population.push_back(particle[i].Xid);
    }

    int bestIdx = 0;
    for (int i = 1; i < popSize; i++) {
        if (particle[i].fitness < particle[bestIdx].fitness) bestIdx = i;
    }

    globalBest = particle[bestIdx].Xid;
    globalBestFE = particle[bestIdx].fitness;
    globalBestStructure = particle[bestIdx].structure;

    for (int cycle = 0; cycle < maxIterations; cycle++) {
        // ----- 1. Employed Bees Phase -----
        for (int i = 0; i < popSize; i++) {
            int k = i;
            do { k = randInt(0, popSize - 1); } while (k == i && popSize > 1);

            int rule = randInt(0, 4);
            Individual cand = applyRuleUpdate(particle[i].Xid, particle[k].Xid, rule);
            string cst;
            double cf = calculateFitness(cand, cst);

            if (cf <= particle[i].fitness) {
                particle[i].Xid = cand;
                particle[i].fitness = cf;
                particle[i].structure = cst;
                particle[i].trial = 0;
                if (cf < particle[i].Pid_val) {
                    particle[i].Pid = cand;
                    particle[i].Pid_val = cf;
                }
            } else {
                particle[i].trial++;
            }

            if (cf < globalBestFE) {
                globalBestFE = cf;
                globalBest = cand;
                globalBestStructure = cst;
            }
        }

        // ----- 2. Probabilities for Onlookers -----
        vector<double> fit(popSize);
        double sumFit = 0.0;
        for (int i = 0; i < popSize; i++) {
            fit[i] = computeFit(particle[i].fitness);
            sumFit += fit[i];
        }
        vector<double> prob(popSize);
        for (int i = 0; i < popSize; i++) {
            prob[i] = (sumFit > 0.0) ? fit[i] / sumFit : 1.0 / popSize;
        }

        // ----- 3. Onlooker Bees Phase -----
        int idx = 0, count = 0;
        while (count < popSize) {
            if (randDouble(0.0, 1.0) < prob[idx]) {
                count++;
                int k = idx;
                do { k = randInt(0, popSize - 1); } while (k == idx && popSize > 1);

                int rule = randInt(0, 4);
                Individual cand = applyRuleUpdate(particle[idx].Xid, particle[k].Xid, rule);
                string cst;
                double cf = calculateFitness(cand, cst);

                if (cf <= particle[idx].fitness) {
                    particle[idx].Xid = cand;
                    particle[idx].fitness = cf;
                    particle[idx].structure = cst;
                    particle[idx].trial = 0;
                    if (cf < particle[idx].Pid_val) {
                        particle[idx].Pid = cand;
                        particle[idx].Pid_val = cf;
                    }
                } else {
                    particle[idx].trial++;
                }

                if (cf < globalBestFE) {
                    globalBestFE = cf;
                    globalBest = cand;
                    globalBestStructure = cst;
                }
            }
            idx = (idx + 1) % popSize;
        }

        // ----- 4. Scout Bees Phase -----
        int maxTrialIdx = 0;
        for (int i = 1; i < popSize; i++) {
            if (particle[i].trial > particle[maxTrialIdx].trial) maxTrialIdx = i;
        }

        if (particle[maxTrialIdx].trial > updateTryLimit) {
            Individual scoutCand = kOptSearch(globalBest, KOptMaxTrial);
            string sst;
            double sf = calculateFitness(scoutCand, sst);

            if (sf >= particle[maxTrialIdx].fitness) {
                scoutCand = randomPermutation();
                sf = calculateFitness(scoutCand, sst);
            }

            particle[maxTrialIdx].Xid = scoutCand;
            particle[maxTrialIdx].fitness = sf;
            particle[maxTrialIdx].structure = sst;
            particle[maxTrialIdx].trial = 0;

            if (sf < globalBestFE) {
                globalBestFE = sf;
                globalBest = scoutCand;
                globalBestStructure = sst;
            }
        }

        // Update population vector for RNA_Core
        population.clear();
        for (int i = 0; i < popSize; i++) {
            population.push_back(particle[i].Xid);
        }

        int curBest = 0;
        for (int i = 1; i < popSize; i++) {
            if (particle[i].fitness < particle[curBest].fitness) curBest = i;
        }
        logIteration(cycle + 1, particle[curBest].fitness, particle[curBest].structure);
    }

    hillClimbBest();
}