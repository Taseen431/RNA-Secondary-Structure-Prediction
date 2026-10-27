/*
 * core.h  --  Shared RNA framework (declarations only)
 *
 * Contains: includes, structs, the RNA_Core base class declaration,
 * the extern random generator, and result-I/O function declarations.
 * The four algorithms (BOA / PSO / GA / ABCO) inherit RNA_Core and only
 * implement their own run() + operators; every RNA-model routine lives here.
 */

#ifndef CORE_H
#define CORE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <limits>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <cctype>

using namespace std;

// ----------------------------- RNG (defined in core.cpp) ------------------
extern mt19937 rng;
double randDouble(double lo = 0.0, double hi = 1.0);
int    randInt(int lo, int hi);   // inclusive

// ----------------------------- Data structures ----------------------------
struct Stem {
    int start;   // start index (0-based)
    int end;     // end index (0-based)
    int length;  // number of consecutive base pairs
};

// An "individual" (butterfly / particle / chromosome / food source) is a
// sequence of 1-based stem indices.
typedef vector<int> Individual;

struct Metrics {
    int TP = 0, FP = 0, FN = 0;
    double sensitivity = 0, specificity = 0, fmeasure = 0, inf = 0;
};

// One logged snapshot: the best-so-far result at a given iteration.
struct IterationRecord {
    int    iteration;
    double energy;
    string structure;
    Metrics metrics;
};

// ----------------------------- Free helper --------------------------------
bool isValidPair(char a, char b);   // canonical base pair check

// ==========================================================================
//  RNA_Core : shared RNA model + state.  Algorithms derive from this.
// ==========================================================================
class RNA_Core {
public:
    static const int MAX_LEN = 15;   // max individual length (shared cap)

    // ---- run configuration / identity ----
    string algorithm;
    string sequenceID;
    string sequence;
    string knownStructure;
    int    n;
    int    populationSize;
    int    maxIterations;
    double parameter;                // 5th input field (BOA switch prob.)

    // ---- shared state ----
    vector<vector<int>> matrix;
    vector<Stem>        stemList;
    vector<Individual>  population;

    Individual globalBest;
    double     globalBestFE;
    string     globalBestStructure;

    vector<IterationRecord> iterationLog;   // one entry per iteration

    RNA_Core(string algo, string seqID, string seq, string known, int popSize, int iterations, double param);
    virtual ~RNA_Core() {}

    // ---- shared RNA framework (implemented in core.cpp) ----
    void        extractStems();
    int         getPairType(char a, char b);
    double      getStackingEnergy(char a, char b, char c, char d);
    double      calculateFreeEnergy(const string &structure);
    string      buildStructure(const Individual &B);
    Individual  optimizeButterfly(const Individual &B);
    double      calculateFitness(const Individual &B, string &structOut);
    Individual  randomPermutation();
    void        hillClimbBest();
    vector<pair<int,int>> extractPairs(const string &structure);
    Metrics     calculateMetrics(const string &predicted, const string &known);
    Metrics     evaluateResult();                 // metrics of globalBest vs known
    void        logIteration(int iter, double energy, const string &structure);

    // ---- console diagnostics ----
    void printMatrix();
    void printStems();
    void printResult(double elapsedSec);

    // ---- algorithm hook ----
    virtual void run() = 0;
};

// ----------------------------- Result I/O ---------------------------------
string currentDateTime();
void initCSV(const string &path);
void appendToCSV(const string &path, const string &runID,const string &dateTime, RNA_Core *solver);
void saveResultToText(const string &path, const string &runID,const string &dateTime, RNA_Core *solver,double elapsedSec);

#endif // CORE_H
