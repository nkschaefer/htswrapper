#ifndef _HTSWRAPPER_FUZZY_H
#define _HTSWRAPPER_FUZZY_H
#include <utility>
#include <cstdlib>
#include <set>
#include <string>
#include <vector>
#include "edlib/edlib.h"

/**
 * This file provides a class that can find the best fuzzy match
 * to a sequence among a set of options, using the edlib sequence
 * alignment library, which reports Levenshtein distance between
 * two sequences.
 */

class seq_fuzzy_match{
    private:
        EdlibEqualityPair eqs[4];
        EdlibAlignConfig config;
        std::vector<std::string> seqs;
        std::string seqs_concat;
        std::vector<int> starts;
        std::vector<int> ends;
        bool global;
    public:
        int mismatches;
        seq_fuzzy_match(std::vector<std::string>& seqs, int max_edit_dist=-1, bool penalize_n=false, bool global=true);
        int match(const char* seq);
};

#endif
