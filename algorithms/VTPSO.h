/*
 * VTPSO.h  --  Variable Threshold Particle Swarm Optimization (declaration only)
 *
 * Formulated for RNA secondary structure prediction with explicit VTPSOParticle class.
 * Inherits from RNA_Core and follows OPTIMIZER_CONTRACT.md.
 */
#ifndef VTPSO_H
#define VTPSO_H

#include "core.h"

struct VTPSOParticle {
    Individual Xid; // Current stem permutation
    Individual Pid; // Personal best stem permutation
    double fitness = 1e20;
    double Pid_val = 1e20;
    string structure;
};

class VTPSO : public RNA_Core {
public:
    double w = 0.729;    // Inertia weight
    double c1 = 1.49445;  // Cognitive coefficient
    double c2 = 1.49445;  // Social coefficient
    vector<VTPSOParticle> particle;

    VTPSO(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("VTPSO", seqID, seq, known, popSize, iterations, param) {}

    void run() override;
};

#endif // VTPSO_H
