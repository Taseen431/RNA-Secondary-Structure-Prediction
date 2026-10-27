/*
 * ABCSS.h  --  Artificial Bee Colony Selective Search (declaration only)
 *
 * Formulated for RNA secondary structure prediction with explicit Solution class.
 * Inherits from RNA_Core and follows OPTIMIZER_CONTRACT.md.
 */
#ifndef ABCSS_H
#define ABCSS_H

#include "core.h"

struct ABCSSParticle {
    Individual Xid;     // Current stem permutation (position)
    Individual Pid;     // Personal best stem permutation
    double fitness = 1e20;
    double Pid_val = 1e20;
    string structure;
    int trial = 0;
};

class ABCSS : public RNA_Core {
public:
    int updateTryLimit = 5;
    int KOptMaxTrial = 10;
    vector<ABCSSParticle> particle;

    ABCSS(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("ABCSS", seqID, seq, known, popSize, iterations, param) {}

    double computeFit(double f);
    Individual applyRuleUpdate(const Individual &ind, const Individual &partner, int ruleType);
    Individual kOptSearch(const Individual &ind, int maxTrial);

    void run() override;
};

#endif // ABCSS_H
