/*
 * BeeColony.cpp  --  Bee Colony Optimization (implementation)
 *
 * The five phases and their order are unchanged from the TSP original:
 *   waggle_Dance -> observe_Dance -> forage_byTranserule -> Fitness -> maintenance.
 * What changed is the problem underneath.  A city is now a candidate stem, a tour
 * is an ordered list of at most MAX_LEN stem indices, and tour length is the free
 * energy returned by the shared calculateFitness().  The distance matrix is gone:
 * stems have no pairwise distance, and m x m would be 212 MB at m = 5151, so the
 * transition rule uses a per-stem heuristic instead (eta = stem length, since a
 * longer stem contributes more base pairs).
 * "parameter" is lambda, the loyalty probability.
 */

#include "BeeColony.h"

// Free energy (lower is better) -> advertising strength (higher is better).
// Energies go negative, so the original's 1.0/fitness is undefined here; this is
// the ABCO::computeFit shape (abco.cpp:11-14), copied because the core is frozen.
static double advertisingStrength(double f) {
    if (f >= 0.0) return 1.0 / (1.0 + f);
    return 1.0 + fabs(f);
}

void BeeColony::allClear() {
    allBees.clear();
    dance.clear();
    pathCost.clear();
    stemHeuristic.clear();
    allStems.clear();
    posInPreferred.clear();
    posInAllowed.clear();
    inAllowed.clear();
    pho.clear();
    besttour.clear();
    Tour.clear();
}

void BeeColony::initialize(int cit, int pop) {
    maxCost = 0.0;
    sumCost = 0.0;
    worstCost = 0.0;
    SD = 0.0;
    a = 0.0;
    b = 0.0;
    best = numeric_limits<double>::max();
    n = 0;
    pfColony = 0.0;

    lambda = parameter;                       // loyalty comes from the 5th input field
    if (lambda <= 0.0 || lambda >= 1.0) lambda = 0.99;

    // TSP used one bee per city.  Here a "city" is a candidate stem, so that rule
    // would ask for 5151 bees on seq20 and put BCO on a different evaluation
    // budget from every other optimizer.  One bee per population slot instead.
    nBee = pop < 1 ? 1 : pop;
    stemCount = cit;

    // eta(s) = stem length: a longer stem contributes more base pairs and more
    // stacking energy, so it plays the role a short edge plays in TSP.
    stemHeuristic.assign(stemCount + 1, 0.0);
    stemEta.assign(stemCount + 1, 0.0);
    for (int i = 0; i < stemCount; i++) {
        allStems.insert(i + 1);
        stemHeuristic[i + 1] = (double)stemList[i].length;
        stemEta[i + 1] = pow(stemHeuristic[i + 1], beta1);   // constant per run
    }
    posInPreferred.assign(stemCount + 1, -1);
    posInAllowed.assign(stemCount + 1, -1);
    inAllowed.assign(stemCount + 1, 1);
    pho.assign(stemCount + 1, 0.0);

    Bee tempBee;
    for (int i = 0; i < nBee; i++) {
        allBees.push_back(tempBee);
        Tour.push_back(Individual());
    }

    for (int i = 0; i < nBee; i++) {          // give a tour to each bee randomly
        Individual temp = randomPermutation();  // conflict-free and MAX_LEN-capped
        allBees[i].tour = temp;
        allBees[i].bestPath = temp;
        allBees[i].bestPathLength = numeric_limits<double>::max();
        allBees[i].danceDuration = 0;
        allBees[i].pf = 0.0;
        allBees[i].fitness = numeric_limits<double>::max();
        Tour[i] = temp;
    }
    Fitness();
}

void BeeColony::Fitness() {
    pathCost.clear();
    sumCost = 0.0;
    worstCost = 0.0;
    a = 0.0;
    b = 0.0;

    population.clear();
    for (int i = 0; i < nBee; i++) {
        string st;
        double temp = calculateFitness(allBees[i].tour, st);
        allBees[i].fitness = temp;
        allBees[i].structure = st;
        if (temp < allBees[i].bestPathLength) {
            allBees[i].bestPathLength = temp;
            allBees[i].bestPath = allBees[i].tour;
        }
        population.push_back(allBees[i].tour);
        pathCost.push_back(temp);
        if (temp > maxCost) maxCost = temp;
        if (temp > worstCost) worstCost = temp;
        sumCost += temp;
        a += temp * temp;
        b += temp;
    }
    double mean = b / nBee;
    double var = (a / nBee) - (mean * mean);
    SD = (var > 0.0) ? sqrt(var) : 0.0;
}

void BeeColony::waggle_Dance() {
    // A bee advertises when its solution is at least as good as the colony mean.
    // The original tested fitness < bestPathLength, which Fitness() has just made
    // false for every bee, leaving n == 0 and pfColony /= n dividing by zero.
    double mean = 0.0;
    for (int i = 0; i < nBee; i++) mean += allBees[i].fitness;
    mean /= nBee;

    for (int i = 0; i < nBee; i++) allBees[i].pf = advertisingStrength(allBees[i].fitness);

    // bees already advertising keep their slot until their duration runs out;
    // only new qualifiers join, so the pool cannot fill with duplicates
    vector<char> inPool(nBee, 0);
    for (size_t d = 0; d < dance.size(); d++) inPool[dance[d].first] = 1;

    Individual joining;
    for (int i = 0; i < nBee; i++)
        if (!inPool[i] && allBees[i].fitness <= mean) joining.push_back(i);

    pfColony = 0.0;
    for (size_t d = 0; d < dance.size(); d++) pfColony += allBees[dance[d].first].pf;
    for (size_t j = 0; j < joining.size(); j++) pfColony += allBees[joining[j]].pf;

    n = (int)dance.size() + (int)joining.size();
    if (n == 0) { pfColony = 0.0; return; }
    pfColony /= n;

    for (size_t j = 0; j < joining.size(); j++) {
        int who = joining[j];
        int dur = (int)floor(k * allBees[who].pf / pfColony);
        if (dur < 1) dur = 1;
        allBees[who].danceDuration = dur;
        dance.push_back(make_pair(who, dur));
    }
    n = (int)dance.size();
}

void BeeColony::observe_Dance() {
    for (int i = 0; i < nBee; i++) {
        if (n > 0) {
            // A bee advertising weakly relative to the colony is the one that
            // should follow somebody else.  The original compared a tour length
            // against a mean of transformed fitness -- different units, opposite
            // directions -- so the comparison is made on pf here.
            double rel = (pfColony > 0.0) ? allBees[i].pf / pfColony : 1.0;
            double follow = -1.0;
            if (rel < 0.50)      follow = 0.60;
            else if (rel < 0.65) follow = 0.20;
            else if (rel < 0.85) follow = 0.02;

            if (follow < 0.0) {
                allBees[i].preferredPath = allBees[i].tour;      // strong: keeps its own
            } else if (randDouble(0.0, 1.0) <= follow) {         // follow a dance
                int whom = randInt(0, n - 1);
                allBees[i].preferredPath = allBees[dance[whom].first].tour;
            } else {                                             // follow the previous bee
                allBees[i].preferredPath = (i == 0) ? allBees[nBee - 1].tour : allBees[i - 1].tour;
            }
        } else {
            // nothing in the dance pool, so follow the previous bee in circular order
            allBees[i].preferredPath = (i == 0) ? allBees[nBee - 1].tour : allBees[i - 1].tour;
        }
    }
}

void BeeColony::forage_byTranserule() {
    int cap = MAX_LEN;

    // The allowed list is built once per cycle, not once per bee: a bee removes at
    // most cap stems, so restoring it afterwards costs O(cap) rather than O(m).
    Individual allowed;
    allowed.reserve(stemCount);
    for (int s = 1; s <= stemCount; s++) {
        allowed.push_back(s);
        posInAllowed[s] = s - 1;
        inAllowed[s] = 1;
    }

    // swap two slots of the allowed list, keeping posInAllowed consistent
    auto swapSlots = [&](int x, int y) {
        if (x == y) return;
        int sx = allowed[x], sy = allowed[y];
        allowed[x] = sy; allowed[y] = sx;
        posInAllowed[sy] = x; posInAllowed[sx] = y;
    };

    Individual removed;
    for (int bIdx = 0; bIdx < nBee; bIdx++) {
        Bee &bee = allBees[bIdx];
        bee.tour.clear();
        removed.clear();

        const Individual &pref = bee.preferredPath;
        int prefLen = (int)pref.size();
        for (int j = 0; j < prefLen; j++)
            if (pref[j] >= 1 && pref[j] <= stemCount) posInPreferred[pref[j]] = j;

        int F = (prefLen > 0) ? pref[0] : 0;     // the preferred next stem
        int A = stemCount;                       // live size of the allowed list

        for (int step = 0; step < cap && A > 0; step++) {
            bool fAllowed = (F >= 1 && F <= stemCount && inAllowed[F] != 0);

            // Weighing every allowed stem is O(m) per placement -- 5151 on seq20,
            // billions of weights per run.  Draw a bounded candidate list instead
            // (standard ACO practice).  The sample is uniform, so every stem stays
            // reachable and the search space is unchanged; only the per-step
            // breadth is bounded.  When A <= CANDIDATE_LIST this is exact.
            int L = (A < CANDIDATE_LIST) ? A : CANDIDATE_LIST;
            if (L < A) {                         // sampling only bites when A is large
                for (int t = 0; t < L; t++) swapSlots(t, randInt(t, A - 1));
                if (fAllowed && posInAllowed[F] >= L) swapSlots(L - 1, posInAllowed[F]);
            }

            // the transition probability takes only two values per placement -- the
            // loyal one for F and a shared one for everything else -- so both are
            // raised to alpha once here rather than once per candidate, and the
            // heuristic term is the stemEta lookup precomputed in initialize()
            double wF     = pow(lambda, alpha);
            double wOther = fAllowed
                          ? pow((A > 1) ? (1.0 - lambda) / (A - 1) : (1.0 - lambda), alpha)
                          : pow(1.0 / A, alpha);

            double sum = 0.0;
            for (int t = 0; t < L; t++) {
                int s = allowed[t];
                double w = (s == F ? wF : wOther) * stemEta[s];
                pho[s] = w;
                sum += w;
            }

            // roulette wheel over the candidate list
            int chosenSlot = -1;
            if (sum > 0.0) {
                double see = randDouble(0.0, 1.0) * sum;
                for (int t = 0; t < L; t++) {
                    see -= pho[allowed[t]];
                    if (see <= 0.0) { chosenSlot = t; break; }
                }
                if (chosenSlot < 0) chosenSlot = L - 1;
            } else {
                chosenSlot = randInt(0, L - 1);
            }

            int chosen = allowed[chosenSlot];
            bee.tour.push_back(chosen);
            removed.push_back(chosen);
            inAllowed[chosen] = 0;
            swapSlots(chosenSlot, A - 1);        // swap-remove from the allowed list
            A--;

            int t2 = posInPreferred[chosen];     // advance along the preferred path
            F = (t2 >= 0 && prefLen > 0) ? pref[(t2 + 1) % prefLen] : 0;
        }

        for (int j = 0; j < prefLen; j++)
            if (pref[j] >= 1 && pref[j] <= stemCount) posInPreferred[pref[j]] = -1;

        bee.allowedStems.assign(allowed.begin(), allowed.begin() + A);
        for (size_t r = 0; r < removed.size(); r++) inAllowed[removed[r]] = 1;

        bee.tour = optimizeButterfly(bee.tour);
    }
}

void BeeColony::maintenance_ofDancePool() {
    for (vector<pair<int, int> >::iterator it = dance.begin(); it != dance.end(); ) {
        it->second--;
        if (it->second <= 0) it = dance.erase(it);   // erasing through a live
        else ++it;                                   // iterator is undefined
    }
    n = (int)dance.size();
}

void BeeColony::operation() {
    for (int iter = 0; iter < maxIterations; iter++) {
        waggle_Dance();
        observe_Dance();
        forage_byTranserule();
        Fitness();
        maintenance_ofDancePool();

        int curBest = 0;
        best = allBees[0].fitness;
        for (int i = 1; i < nBee; i++)
            if (allBees[i].fitness < best) { best = allBees[i].fitness; curBest = i; }

        if (best < globalBestFE) {
            globalBestFE = best;
            globalBest = allBees[curBest].tour;
            globalBestStructure = allBees[curBest].structure;
        }
        // log the best of the CURRENT colony at this iteration
        logIteration(iter + 1, allBees[curBest].fitness, allBees[curBest].structure);
    }
    besttour = globalBest;
}

void BeeColony::run() {
    extractStems();
    allClear();
    initialize((int)stemList.size(), populationSize);

    int bestIdx = 0;
    for (int i = 1; i < nBee; i++)
        if (allBees[i].fitness < allBees[bestIdx].fitness) bestIdx = i;
    globalBest = allBees[bestIdx].tour;
    globalBestFE = allBees[bestIdx].fitness;
    globalBestStructure = allBees[bestIdx].structure;

    operation();
    hillClimbBest();
}
