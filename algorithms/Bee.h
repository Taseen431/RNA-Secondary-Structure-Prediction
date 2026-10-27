/*
 * Bee.h  --  A single artificial bee for Bee Colony Optimization (declaration only)
 *
 * One bee owns the solution it is building this cycle (tour), the solution it has
 * chosen to follow while foraging (preferredPath), the stems still available to it
 * during construction (allowedStems), and its own best-so-far.  A "tour" here is an
 * RNA candidate: an ordered list of 1-based stem indices, at most MAX_LEN long,
 * which the shared framework decodes and scores.
 */

#ifndef BEE_H
#define BEE_H

#include "core.h"

class Bee {
public:
    Individual tour;            // stem indices chosen this cycle (the candidate)
    Individual preferredPath;   // the tour this bee follows while foraging
    Individual allowedStems;    // stem indices still selectable during construction
    vector<double> p;           // per-stem transition probabilities (scratch)
    int    danceDuration;       // cycles this bee's tour stays in the dance pool
    double fitness;             // free energy of tour  (lower is better)
    double pf;                  // advertising strength (higher is better)
    Individual bestPath;        // best tour this bee has held
    double bestPathLength;      // free energy of bestPath
    string structure;           // dot-bracket decoding of tour, kept for logging
};

#endif // BEE_H
