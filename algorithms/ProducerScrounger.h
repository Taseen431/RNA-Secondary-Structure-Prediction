/*
 * ProducerScrounger.h  --  Group Search Optimizer (GSO) (declaration only)
 *
 * Formulated for RNA secondary structure prediction with explicit GSOMember class.
 * Inherits from RNA_Core and follows OPTIMIZER_CONTRACT.md.
 */
#ifndef PRODUCERSCROUNGER_H
#define PRODUCERSCROUNGER_H

#include "core.h"

struct GSOMember {
    Individual position;
    double fitness = 1e20;
    string structure;
    int type = 1; // 0: Producer, 1: Scrounger, 2: Ranger
};

class PrdScrg : public RNA_Core {
public:
    double producerRatio = 0.20; // 20% producers
    double rangerRatio = 0.20;   // 20% rangers
    vector<GSOMember> members;

    PrdScrg(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("GSO", seqID, seq, known, popSize, iterations, param) {}

    void run() override;
};

#endif // PRODUCERSCROUNGER_H
