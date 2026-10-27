/*
 * BeeColony.h  --  Bee Colony Optimization (declaration only)
 *
 * Lucic & Teodorovic BCO: every cycle the bees advertise good solutions in a
 * waggle-dance pool, decide whose solution to follow, then rebuild their own
 * solution step by step under a transition rule that mixes loyalty to the followed
 * path with a problem heuristic.  This is NOT abco.*, which is Karaboga's
 * employed / onlooker / scout ABC -- different algorithm, different name.
 *
 * RNA formulation: a "city" is a candidate stem, a "tour" is an ordered list of at
 * most MAX_LEN stem indices, and tour length is the free energy returned by the
 * shared calculateFitness().  The TSP distance matrix is replaced by a per-stem
 * heuristic (stem length), since stems have no pairwise distance and an m x m
 * matrix would be 212 MB at m = 5151.
 * "parameter" is lambda, the loyalty probability.
 */

#ifndef BEECOLONY_H
#define BEECOLONY_H

#include "Bee.h"

class BeeColony : public RNA_Core {
public:
    // BCO hyper-parameters
    double alpha = 1.0;           // exponent on the transition probability
    double beta1 = 2.0;           // exponent on the stem-length heuristic
    int    k = 5000;              // dance-duration scale
    int    CANDIDATE_LIST = 64;   // stems weighed per placement (see the .cpp)

    BeeColony(string seqID, string seq, string known, int popSize, int iterations, double param)
        : RNA_Core("BCO", seqID, seq, known, popSize, iterations, param) {}

    vector<Bee> allBees;
    double lambda;                // loyalty to the preferred path, from parameter
    double pfColony;              // mean advertising strength over the dance pool
    vector<pair<int, int> > dance; // (bee index, cycles remaining in the pool)
    int n;                        // how many bees are currently dancing
    double best;
    vector<double> pathCost;
    double maxCost;
    int nBee;
    int stemCount;
    vector<double> stemHeuristic; // replaces AdjMat: eta(s) = stem length, 1-based
    vector<double> stemEta;       // eta(s)^beta1, precomputed once per run
    set<int> allStems;
    vector<int> posInPreferred;   // stem index -> its position in preferredPath
    vector<int> posInAllowed;     // stem index -> its slot in a bee's allowedStems
    vector<char> inAllowed;       // stem index -> still selectable this construction
    vector<double> pho;           // transition weights, indexed by stem
    Individual besttour;
    vector<Individual> Tour;

    double sumCost;
    double worstCost;
    double SD;
    double a;
    double b;

    void initialize(int, int);
    void waggle_Dance();
    void observe_Dance();
    void forage_byTranserule();
    void Fitness();
    void maintenance_ofDancePool();
    void operation();
    void allClear();

    void run() override;
};

#endif // BEECOLONY_H
