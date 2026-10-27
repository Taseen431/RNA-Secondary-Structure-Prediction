/*
 * woa.h  --  Whale Optimization Algorithm (declaration only)
 */
#ifndef WOA_H
#define WOA_H

#include "core.h"

class WOA : public RNA_Core {
public:
    // WOA hyper-parameters ("parameter" supplies p, the encircling / spiral
    // branch probability -- the same role parameter plays for BOA)
    double SPIRAL_B         = 1.0;    // shape constant of the logarithmic spiral
    double SELECT_THRESHOLD = 0.50;   // decoder threshold (kept equal to PSO's configured value)

    WOA(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("WOA", seqID, seq, known, popSize, iterations, param) {}

    Individual decode(const vector<double> &position);
    void run() override;
};

#endif // WOA_H
