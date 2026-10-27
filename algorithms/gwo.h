/*
 * gwo.h  --  Grey Wolf Optimizer (declaration only)
 */
#ifndef GWO_H
#define GWO_H

#include "core.h"

class GWO : public RNA_Core {
public:
    // GWO hyper-parameters ("parameter" is parsed but unused: canonical GWO has
    // no free parameter beyond the linear decay of a)
    double SELECT_THRESHOLD = 0.50;   // decoder threshold (kept equal to PSO's configured value)

    GWO(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("GWO", seqID, seq, known, popSize, iterations, param) {}

    Individual decode(const vector<double> &position);
    void run() override;
};

#endif // GWO_H
