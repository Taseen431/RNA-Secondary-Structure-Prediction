/*
 * SMO.h  --  Spider Monkey Optimization (declaration only)
 *
 * Formulated for RNA secondary structure prediction with explicit SpiderMonkey class.
 * Inherits from RNA_Core and follows OPTIMIZER_CONTRACT.md.
 */
#ifndef SMO_H
#define SMO_H

#include "core.h"

struct SpiderMonkey {
    Individual position;
    double fitness = 1e20;
    string structure;
    int groupID = 0;
};

class SMO : public RNA_Core {
public:
    int maxGroups = 5;
    int localLeaderLimit = 20;
    int globalLeaderLimit = 10;
    double pr = 0.2; // Perturbation probability
    vector<SpiderMonkey> monkeys;

    SMO(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("SMO", seqID, seq, known, popSize, iterations, param) {}

    void run() override;
};

#endif // SMO_H