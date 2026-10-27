/*
 * SMO.cpp  --  Spider Monkey Optimization (implementation)
 *
 * Formulated for RNA secondary structure prediction with explicit SpiderMonkey class.
 * Inherits from RNA_Core and follows OPTIMIZER_CONTRACT.md.
 */

#include "SMO.h"

// Helper function to update monkey position guided by a leader and partner
static Individual updateMonkeyPosition(const Individual &current, const Individual &leader, const Individual &partner, int m, double pr, RNA_Core *core) {
    Individual cand = current;
    int nStems = (int)cand.size();

    if (randDouble(0.0, 1.0) < pr && !leader.empty()) {
        if (nStems > 0) {
            int idx = randInt(0, nStems - 1);
            cand[idx] = leader[randInt(0, (int)leader.size() - 1)];
        } else {
            cand.push_back(leader[0]);
        }
    }

    if (randDouble(0.0, 1.0) < pr && !partner.empty()) {
        if ((int)cand.size() < 15) {
            cand.push_back(partner[randInt(0, (int)partner.size() - 1)]);
        } else if (!cand.empty()) {
            cand[randInt(0, (int)cand.size() - 1)] = partner[randInt(0, (int)partner.size() - 1)];
        }
    }

    if (randDouble(0.0, 1.0) < 0.15 && (int)cand.size() >= 2) {
        swap(cand[randInt(0, (int)cand.size() - 1)], cand[randInt(0, (int)cand.size() - 1)]);
    }

    return core->optimizeButterfly(cand);
}

void SMO::run() {
    extractStems();
    int m = (int)stemList.size();
    if (m == 0) return;

    int popSize = populationSize;
    if (popSize < 6) popSize = 6;

    monkeys.clear();
    monkeys.resize(popSize);
    population.clear();

    for (int i = 0; i < popSize; i++) {
        monkeys[i].position = randomPermutation();
        monkeys[i].fitness = calculateFitness(monkeys[i].position, monkeys[i].structure);
        monkeys[i].groupID = 0;
        population.push_back(monkeys[i].position);
    }

    int bestIdx = 0;
    for (int i = 1; i < popSize; i++) {
        if (monkeys[i].fitness < monkeys[bestIdx].fitness) bestIdx = i;
    }

    globalBest = monkeys[bestIdx].position;
    globalBestFE = monkeys[bestIdx].fitness;
    globalBestStructure = monkeys[bestIdx].structure;

    int numGroups = 1;
    vector<vector<int>> groups(1);
    for (int i = 0; i < popSize; i++) groups[0].push_back(i);

    vector<Individual> localLeaders(1, globalBest);
    vector<double> localLeadersFE(1, globalBestFE);
    vector<int> localStagnation(1, 0);
    int globalStagnation = 0;

    for (int iter = 0; iter < maxIterations; iter++) {
        // ----- 1. Local Leader Phase (LLP) -----
        for (int g = 0; g < numGroups; g++) {
            const vector<int> &groupMembers = groups[g];
            int gSize = (int)groupMembers.size();

            for (int k = 0; k < gSize; k++) {
                int idx = groupMembers[k];
                int partnerIdx = groupMembers[randInt(0, gSize - 1)];

                Individual cand = updateMonkeyPosition(monkeys[idx].position, localLeaders[g], monkeys[partnerIdx].position, m, pr, this);
                string cst;
                double cf = calculateFitness(cand, cst);

                if (cf <= monkeys[idx].fitness) {
                    monkeys[idx].position = cand;
                    monkeys[idx].fitness = cf;
                    monkeys[idx].structure = cst;
                }
            }
        }

        // ----- 2. Global Leader Phase (GLP) -----
        vector<double> fit(popSize);
        double sumFit = 0.0;
        for (int i = 0; i < popSize; i++) {
            fit[i] = (monkeys[i].fitness >= 0.0) ? 1.0 / (1.0 + monkeys[i].fitness) : 1.0 + fabs(monkeys[i].fitness);
            sumFit += fit[i];
        }

        for (int i = 0; i < popSize; i++) {
            double prob = (sumFit > 0.0) ? (0.9 * fit[i] / sumFit + 0.1) : 0.2;
            if (randDouble(0.0, 1.0) < prob) {
                int partnerIdx = randInt(0, popSize - 1);
                Individual cand = updateMonkeyPosition(monkeys[i].position, globalBest, monkeys[partnerIdx].position, m, pr, this);
                string cst;
                double cf = calculateFitness(cand, cst);

                if (cf <= monkeys[i].fitness) {
                    monkeys[i].position = cand;
                    monkeys[i].fitness = cf;
                    monkeys[i].structure = cst;
                }
            }
        }

        // ----- 3. Local Leader Learning (LLL) & Global Leader Learning (GLL) -----
        bool globalImproved = false;
        for (int g = 0; g < numGroups; g++) {
            bool localImproved = false;
            for (int idx : groups[g]) {
                if (monkeys[idx].fitness < localLeadersFE[g]) {
                    localLeaders[g] = monkeys[idx].position;
                    localLeadersFE[g] = monkeys[idx].fitness;
                    localImproved = true;
                }
                if (monkeys[idx].fitness < globalBestFE) {
                    globalBest = monkeys[idx].position;
                    globalBestFE = monkeys[idx].fitness;
                    globalBestStructure = monkeys[idx].structure;
                    globalImproved = true;
                }
            }
            if (localImproved) localStagnation[g] = 0;
            else localStagnation[g]++;
        }

        if (globalImproved) globalStagnation = 0;
        else globalStagnation++;

        // ----- 4. Local Leader Decision (LLD) -----
        for (int g = 0; g < numGroups; g++) {
            if (localStagnation[g] > localLeaderLimit) {
                localStagnation[g] = 0;
                for (int idx : groups[g]) {
                    if (monkeys[idx].position != localLeaders[g]) {
                        monkeys[idx].position = updateMonkeyPosition(monkeys[idx].position, globalBest, localLeaders[g], m, 0.5, this);
                        monkeys[idx].fitness = calculateFitness(monkeys[idx].position, monkeys[idx].structure);
                    }
                }
            }
        }

        // ----- 5. Global Leader Decision (GLD) -----
        if (globalStagnation > globalLeaderLimit) {
            globalStagnation = 0;
            if (numGroups < maxGroups) {
                numGroups++;
                groups.resize(numGroups);
                localLeaders.resize(numGroups);
                localLeadersFE.resize(numGroups);
                localStagnation.resize(numGroups, 0);

                for (int g = 0; g < numGroups; g++) groups[g].clear();
                for (int i = 0; i < popSize; i++) {
                    monkeys[i].groupID = i % numGroups;
                    groups[i % numGroups].push_back(i);
                }
                for (int g = 0; g < numGroups; g++) {
                    int bIdx = groups[g][0];
                    for (int idx : groups[g]) {
                        if (monkeys[idx].fitness < monkeys[bIdx].fitness) bIdx = idx;
                    }
                    localLeaders[g] = monkeys[bIdx].position;
                    localLeadersFE[g] = monkeys[bIdx].fitness;
                }
            } else {
                numGroups = 1;
                groups.resize(1);
                groups[0].clear();
                for (int i = 0; i < popSize; i++) {
                    monkeys[i].groupID = 0;
                    groups[0].push_back(i);
                }
                localLeaders.resize(1); localLeadersFE.resize(1); localStagnation.resize(1, 0);
                localLeaders[0] = globalBest;
                localLeadersFE[0] = globalBestFE;
            }
        }

        population.clear();
        for (int i = 0; i < popSize; i++) {
            population.push_back(monkeys[i].position);
        }

        int curBest = 0;
        for (int i = 1; i < popSize; i++) {
            if (monkeys[i].fitness < monkeys[curBest].fitness) curBest = i;
        }

        logIteration(iter + 1, monkeys[curBest].fitness, monkeys[curBest].structure);
    }

    hillClimbBest();
}