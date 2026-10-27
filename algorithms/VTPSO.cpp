/*
 * VTPSO.cpp  --  Variable Threshold Particle Swarm Optimization (implementation)
 *
 * Formulated for RNA secondary structure prediction with explicit VTPSOParticle class.
 * Inherits from RNA_Core and follows OPTIMIZER_CONTRACT.md.
 */

#include "VTPSO.h"

void VTPSO::run() {
    extractStems();
    int m = (int)stemList.size();
    if (m == 0) return;

    int popSize = populationSize;
    if (popSize < 4) popSize = 4;

    particle.clear();
    particle.resize(popSize);
    population.clear();

    for (int i = 0; i < popSize; i++) {
        particle[i].Xid = randomPermutation();
        particle[i].Pid = particle[i].Xid;
        particle[i].fitness = calculateFitness(particle[i].Xid, particle[i].structure);
        particle[i].Pid_val = particle[i].fitness;
        population.push_back(particle[i].Xid);
    }

    int bestIdx = 0;
    for (int i = 1; i < popSize; i++) {
        if (particle[i].fitness < particle[bestIdx].fitness) bestIdx = i;
    }

    globalBest = particle[bestIdx].Xid;
    globalBestFE = particle[bestIdx].fitness;
    globalBestStructure = particle[bestIdx].structure;

    for (int iter = 0; iter < maxIterations; iter++) {
        double thresholdProb = 0.80 - (0.50 * (double)iter / maxIterations);

        for (int i = 0; i < popSize; i++) {
            Individual cand = particle[i].Xid;
            int nStems = (int)cand.size();

            // Cognitive Component (Move towards personal best Pid)
            if (randDouble(0.0, 1.0) < c1 / 2.0 && !particle[i].Pid.empty()) {
                if (nStems > 0) {
                    int idx = randInt(0, nStems - 1);
                    cand[idx] = particle[i].Pid[randInt(0, (int)particle[i].Pid.size() - 1)];
                }
            }

            // Social Component (Move towards global best globalBest)
            if (randDouble(0.0, 1.0) < c2 / 2.0 && !globalBest.empty()) {
                if ((int)cand.size() < MAX_LEN) {
                    cand.push_back(globalBest[randInt(0, (int)globalBest.size() - 1)]);
                } else if (!cand.empty()) {
                    cand[randInt(0, (int)cand.size() - 1)] = globalBest[randInt(0, (int)globalBest.size() - 1)];
                }
            }

            // Variable Threshold Block Adjustment (Sub-sequence reversal / 2-opt)
            if (randDouble(0.0, 1.0) < thresholdProb && (int)cand.size() >= 3) {
                int l = randInt(0, (int)cand.size() - 2);
                int r = randInt(l + 1, (int)cand.size() - 1);
                reverse(cand.begin() + l, cand.begin() + r + 1);
            }

            particle[i].Xid = optimizeButterfly(cand);
            particle[i].fitness = calculateFitness(particle[i].Xid, particle[i].structure);

            if (particle[i].fitness < particle[i].Pid_val) {
                particle[i].Pid = particle[i].Xid;
                particle[i].Pid_val = particle[i].fitness;
            }

            if (particle[i].fitness < globalBestFE) {
                globalBestFE = particle[i].fitness;
                globalBest = particle[i].Xid;
                globalBestStructure = particle[i].structure;
            }
        }

        population.clear();
        for (int i = 0; i < popSize; i++) {
            population.push_back(particle[i].Xid);
        }

        int curBest = 0;
        for (int i = 1; i < popSize; i++) {
            if (particle[i].fitness < particle[curBest].fitness) curBest = i;
        }

        logIteration(iter + 1, particle[curBest].fitness, particle[curBest].structure);
    }

    hillClimbBest();
}