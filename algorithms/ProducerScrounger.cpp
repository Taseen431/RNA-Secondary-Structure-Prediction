/*
 * ProducerScrounger.cpp  --  Group Search Optimizer (GSO) (implementation)
 *
 * Formulated for RNA secondary structure prediction with explicit GSOMember class.
 * Inherits from RNA_Core and follows OPTIMIZER_CONTRACT.md.
 */

#include "ProducerScrounger.h"

void PrdScrg::run() {
    extractStems();
    int m = (int)stemList.size();
    if (m == 0) return;

    int popSize = populationSize;
    if (popSize < 5) popSize = 5;

    members.clear();
    members.resize(popSize);
    population.clear();

    for (int i = 0; i < popSize; i++) {
        members[i].position = randomPermutation();
        members[i].fitness = calculateFitness(members[i].position, members[i].structure);
        population.push_back(members[i].position);
    }

    int bestIdx = 0;
    for (int i = 1; i < popSize; i++) {
        if (members[i].fitness < members[bestIdx].fitness) bestIdx = i;
    }

    globalBest = members[bestIdx].position;
    globalBestFE = members[bestIdx].fitness;
    globalBestStructure = members[bestIdx].structure;

    int nProducers = max(1, (int)(producerRatio * popSize));
    int nRangers = max(1, (int)(rangerRatio * popSize));

    for (int iter = 0; iter < maxIterations; iter++) {
        // Sort members by fitness
        vector<int> sortedIndices(popSize);
        for (int i = 0; i < popSize; i++) sortedIndices[i] = i;
        sort(sortedIndices.begin(), sortedIndices.end(), [this](int a, int b) {
            return members[a].fitness < members[b].fitness;
        });

        // Assign roles: Top nProducers -> Producers, Last nRangers -> Rangers, Middle -> Scroungers
        for (int i = 0; i < popSize; i++) {
            int idx = sortedIndices[i];
            if (i < nProducers) members[idx].type = 0; // Producer
            else if (i >= popSize - nRangers) members[idx].type = 2; // Ranger
            else members[idx].type = 1; // Scrounger
        }

        // ----- 1. Producer Phase (Greedy local scanning around producer) -----
        for (int i = 0; i < nProducers; i++) {
            int pIdx = sortedIndices[i];
            Individual cand = members[pIdx].position;

            // Perform local scanning / stem swap or insertion
            if (randDouble(0.0, 1.0) < 0.50 && (int)cand.size() >= 2) {
                swap(cand[randInt(0, (int)cand.size() - 1)], cand[randInt(0, (int)cand.size() - 1)]);
            } else if ((int)cand.size() < MAX_LEN) {
                cand.push_back(randInt(1, m));
            }

            cand = optimizeButterfly(cand);
            string cst;
            double cf = calculateFitness(cand, cst);

            if (cf <= members[pIdx].fitness) {
                members[pIdx].position = cand;
                members[pIdx].fitness = cf;
                members[pIdx].structure = cst;
            }

            if (cf < globalBestFE) {
                globalBestFE = cf;
                globalBest = cand;
                globalBestStructure = cst;
            }
        }

        // ----- 2. Scrounger Phase (Follow chosen producer) -----
        for (int i = nProducers; i < popSize - nRangers; i++) {
            int sIdx = sortedIndices[i];
            // Pick a producer to follow
            int chosenProducerIdx = sortedIndices[randInt(0, nProducers - 1)];

            Individual cand = members[sIdx].position;
            const Individual &prodPos = members[chosenProducerIdx].position;

            // Move towards producer (Crossover / Stem Replacement)
            if (!prodPos.empty()) {
                if (randDouble(0.0, 1.0) < 0.50 && !cand.empty()) {
                    cand[randInt(0, (int)cand.size() - 1)] = prodPos[randInt(0, (int)prodPos.size() - 1)];
                } else if ((int)cand.size() < MAX_LEN) {
                    cand.push_back(prodPos[randInt(0, (int)prodPos.size() - 1)]);
                }
            }

            cand = optimizeButterfly(cand);
            string cst;
            double cf = calculateFitness(cand, cst);

            if (cf <= members[sIdx].fitness) {
                members[sIdx].position = cand;
                members[sIdx].fitness = cf;
                members[sIdx].structure = cst;
            }

            if (cf < globalBestFE) {
                globalBestFE = cf;
                globalBest = cand;
                globalBestStructure = cst;
            }
        }

        // ----- 3. Ranger Phase (Random walk for global exploration) -----
        for (int i = popSize - nRangers; i < popSize; i++) {
            int rIdx = sortedIndices[i];
            Individual cand;

            if (randDouble(0.0, 1.0) < 0.50) {
                cand = randomPermutation();
            } else {
                cand = members[rIdx].position;
                if ((int)cand.size() >= 2) {
                    swap(cand[randInt(0, (int)cand.size() - 1)], cand[randInt(0, (int)cand.size() - 1)]);
                }
                if (randDouble(0.0, 1.0) < 0.30 && !cand.empty()) {
                    cand[randInt(0, (int)cand.size() - 1)] = randInt(1, m);
                }
                cand = optimizeButterfly(cand);
            }

            string cst;
            double cf = calculateFitness(cand, cst);

            if (cf <= members[rIdx].fitness) {
                members[rIdx].position = cand;
                members[rIdx].fitness = cf;
                members[rIdx].structure = cst;
            }

            if (cf < globalBestFE) {
                globalBestFE = cf;
                globalBest = cand;
                globalBestStructure = cst;
            }
        }

        population.clear();
        for (int i = 0; i < popSize; i++) {
            population.push_back(members[i].position);
        }

        int curBest = sortedIndices[0];
        logIteration(iter + 1, members[curBest].fitness, members[curBest].structure);
    }

    hillClimbBest();
}
