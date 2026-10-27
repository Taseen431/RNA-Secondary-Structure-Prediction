/*
 * hs.h  --  Harmony Search (declaration only)
 */
#ifndef HS_H
#define HS_H

#include "core.h"

class HS : public RNA_Core {
public:
    // HS hyper-parameters
    double PAR      = 0.30;   // pitch adjustment rate
    int    PITCH_BW = 1;      // pitch bandwidth: +/- this many stem indices

    HS(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("HS", seqID, seq, known, popSize, iterations, param) {}

    Individual improvise();
    void run() override;
};

#endif // HS_H
