/*
 * ParticleSwarm.h  --  Particle Swarm Optimization (declaration only)
 *
 * Formulated for RNA secondary structure prediction with explicit PSOParticle class.
 * Inherits from RNA_Core and follows OPTIMIZER_CONTRACT.md.
 */
#ifndef PARTICLESWARM_H
#define PARTICLESWARM_H

#include "core.h"

struct PSOParticle {
    Individual Xid; // Current position
    Individual Pid; // Personal best position
    double fitness = 1e20;
    double Pid_val = 1e20;
    string structure;
};

class ParticleSwarm : public RNA_Core {
public:
    double w = 0.729;    // Inertia weight
    double c1 = 1.49445;  // Cognitive factor
    double c2 = 1.49445;  // Social factor
    vector<PSOParticle> particle;

    ParticleSwarm(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("ParticleSwarm", seqID, seq, known, popSize, iterations, param) {}

    void run() override;
};

#endif // PARTICLESWARM_H
