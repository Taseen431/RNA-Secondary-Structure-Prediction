/*
 * pso.h  --  Particle Swarm Optimization (declaration only)
 */
#ifndef PSO_H
#define PSO_H

#include "core.h"

class PSO : public RNA_Core {
public:
    // PSO hyper-parameters (instance members: no linkage issues)
    double INERTIA_W = 0.7;
    double COGNITIVE_C1 = 1.5;
    double SOCIAL_C2 = 1.5;
    double VMAX = 1.0;
    double SELECT_THRESHOLD = 0.5;

    PSO(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("PSO", seqID, seq, known, popSize, iterations, param) {
        SELECT_THRESHOLD = param;
    }

    Individual decode(const vector<double> &position);
    void run() override;
};

#endif // PSO_H
