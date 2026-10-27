/*
 * de.h  --  Differential Evolution (declaration only)
 */
#ifndef DE_H
#define DE_H

#include "core.h"

class DE : public RNA_Core {
public:
    // DE hyper-parameters ("parameter" supplies F, the differential weight)
    double CROSSOVER_CR     = 0.90;   // binomial crossover rate
    double SELECT_THRESHOLD = 0.50;   // decoder threshold (kept equal to PSO's configured value)

    DE(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("DE", seqID, seq, known, popSize, iterations, param) {}

    Individual decode(const vector<double> &position);
    void run() override;
};

#endif // DE_H
