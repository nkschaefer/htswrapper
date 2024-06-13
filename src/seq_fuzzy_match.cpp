#include <utility>
#include <cstdlib>
#include <string>
#include <cstring>
#include <vector>
#include <math.h>
#include "edlib/edlib.h"
#include "seq_fuzzy_match.h"

using namespace std;

/**
 * Configure to match trimmed portions of reads (where we expect to see search seqs)
 * against a concatenated, single reference sequence of all search seqs.
 * Determines the correct sequence based on location.
 *
 * Fastest if many search seqs.
 */
void seq_fuzzy_match::set_global(){
    
    // Tell edlib (optionally) to count N - any base as a match 
    eqs[0] = {'N', 'A'};
    eqs[1] = {'N', 'C'};
    eqs[2] = {'N', 'G'};
    eqs[3] = {'N', 'T'};
    this->global = true;
    this->reverse = false;

    string padding = "";
    // How long should spacers between seqs be?
    // Assume these are barcodes (relatively short), and that 
    // seqs aligning to them will be comparable length.
    // Use padding equal to half the length of the longest sequence
    int n_padding = 0;
    for (int i = 0; i < seqs.size(); ++i){
        if (seqs[i].length() > n_padding){
            n_padding = seqs[i].length();
        }
    }
    n_padding = (int)round((double)n_padding/2.0);
    for (int i = 0; i < n_padding; ++i){
        // Choose a character no bases can match
        padding += "?";
    }
    
    starts.clear();
    ends.clear();

    starts.push_back(0);
    seqs_concat = seqs[0];
    int coord = seqs[0].length();
    ends.push_back(coord-1);
    for (int i = 1; i < seqs.size(); ++i){
        starts.push_back(coord + padding.length());
        seqs_concat += padding + seqs[i];
        coord += padding.length() + seqs[i].length();
        ends.push_back(coord-1);
    }
    if (penalize_n){
        this->config = edlibNewAlignConfig(allowed_edit_dist, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0);
    }
    else{
        this->config = edlibNewAlignConfig(allowed_edit_dist, EDLIB_MODE_HW, EDLIB_TASK_LOC, &eqs[0], 4);
    }
}

/**
 * Configure to search a trimmed read, where we expect to see a search seq, 
 * against each possible sequence in a list (one by one). Slower unless
 * few search seqs.
 */
void seq_fuzzy_match::set_single(){
    // Tell edlib (optionally) to count N - any base as a match 
    eqs[0] = {'N', 'A'};
    eqs[1] = {'N', 'C'};
    eqs[2] = {'N', 'G'};
    eqs[3] = {'N', 'T'};
    this->global = false;
    this->reverse = false;
    if (penalize_n){
        this->config = edlibNewAlignConfig(allowed_edit_dist, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
    }
    else{
        this->config = edlibNewAlignConfig(allowed_edit_dist, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, &eqs[0], 4);
    }
}

/**
 * In cases where the position of the search seq in reads is unknown, searches every
 * seq in a list against an entire read. Will set the match_pos variable to the 
 * index of the match within the read, where it is found.
 */
void seq_fuzzy_match::set_reverse(){
    // Tell edlib (optionally) to count N - any base as a match 
    eqs[0] = {'N', 'A'};
    eqs[1] = {'N', 'C'};
    eqs[2] = {'N', 'G'};
    eqs[3] = {'N', 'T'};
    this->global = false;
    this->reverse = true;
    if (penalize_n){
        this->config = edlibNewAlignConfig(allowed_edit_dist, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0);
    }
    else{
        this->config = edlibNewAlignConfig(allowed_edit_dist, EDLIB_MODE_HW, EDLIB_TASK_LOC, &eqs[0], 4);
    }
}

/**
 * Params:
 * seqs: vector of all sequences to match against
 * edit_dist: maximum allowable edit dist between query and seq
 * penalize_n: true = N/ACGT is mismatch. false = N/ACGT is match
 * global: search by making a single reference sequence containing all
 *   sequences to match against, with padding characters between. For each alignment, 
 *   find the best position in this ref seq of a query seq (faster for many seqs)
 *   Default: test query seq against every ref seq until a match is found
 * reverse: use when query seqs are untrimmed  and you don't know where in a read 
 *   matches to search seqs should occur. Each alignment will be every (short) search
 *   sequence against the (long) query sequence, not penalizing gaps on ends.
 */
seq_fuzzy_match::seq_fuzzy_match(vector<string>& seqs, 
    int edit_dist, 
    bool penalize_n, 
    bool global, 
    bool reverse){
    
    this->seqs = seqs;
    this->mismatches = -1;    
    this->allowed_edit_dist = edit_dist;
    this->match_pos = -1;
    this->penalize_n = penalize_n;

    if (global){
        set_global();    
    }
    else if (reverse){
        set_reverse();
    }
    else{
        set_single();
    }
}

int seq_fuzzy_match::match(const char* seq){
    this->match_pos = -1;
    if (global){
        int sl = strlen(seq);
        EdlibAlignResult result = edlibAlign(seq, 
            sl, 
            seqs_concat.c_str(), 
            seqs_concat.length(),
            config);
        mismatches = -1;
        int mindiff = -1;
        int mindiff_idx = -1;
        bool tie = false;
        if (result.status == EDLIB_STATUS_OK){
            for (int j = 0; j < result.numLocations; ++j){
                for (int i = 0; i < starts.size(); ++i){
                    if (starts[i] <= result.startLocations[j] && ends[i] >= result.endLocations[j]){
                        // Valid mapping
                        int diff = abs(sl-(result.endLocations[j]-result.startLocations[j]+1));
                        if (mindiff_idx == -1 || diff < mindiff){
                            mindiff_idx = i;
                            mindiff = diff;
                            tie = false;
                        }
                        else if (mindiff_idx != -1 && diff == mindiff){
                            tie = true;
                        }
                        if (diff == 0){
                            // Be greedy
                            break;
                        }
                    }
                }
            }
        }
        if (mindiff_idx != -1 && !tie){
            mismatches = result.editDistance;
        }
        else{
            mismatches = -1;
            mindiff_idx = -1;
        }
        edlibFreeAlignResult(result);
        return mindiff_idx;
    }
    else if (reverse){
        int sl = strlen(seq);
        int min_ed = -1;
        int min_idx = -1;
        bool tie = false;
        int min_matchpos = -1;
        for (int i = 0; i < seqs.size(); ++i){
            EdlibAlignResult result = edlibAlign(seqs[i].c_str(), 
                seqs[i].length(), 
                seq, 
                sl,
                config);   
            if (result.status == EDLIB_STATUS_OK){
                int mindiff = -1;
                int mindiff_idx = -1;
                int min_startpos = -1;

                for (int j = 0; j < result.numLocations; ++j){   
                    int startpos = result.startLocations[j];              
                    int diff = abs((int)seqs[i].length()-(result.endLocations[j]-result.startLocations[j]+1));
                    if (mindiff == -1 || diff < mindiff){
                        mindiff = diff;
                        mindiff_idx = i;
                        min_startpos = startpos;
                    }
                    if (diff == 0){
                        break;
                    }
                }
                if (mindiff != -1){
                    int ed = result.editDistance + mindiff;
                    if (ed == 0){
                        // Be greedy about perfect matches
                        min_ed = ed;
                        min_idx = i;
                        tie = false;
                        min_matchpos = min_startpos;
                        break;
                    }
                    else if (ed > 0){
                        if (min_ed == -1 || ed < min_ed){
                            min_ed = ed;
                            min_idx = i;
                            tie = false;
                            min_matchpos = min_startpos;
                        }
                        else if (min_ed != -1 && ed == min_ed){
                            tie = true;
                        }
                    }
                }
            }
        }
        if (min_idx != -1 && !tie){
            this->mismatches = min_ed;
            this->match_pos = min_matchpos;
            return min_idx;
        }
        this->mismatches = -1;
        this->match_pos = -1;
        return -1;
    }
    else{
        int sl = strlen(seq);
        int min_ed = -1;
        int min_idx = -1;
        bool tie = false;
        for (int i = 0; i < seqs.size(); ++i){
            EdlibAlignResult result = edlibAlign(seq,
                sl,
                seqs[i].c_str(),
                seqs[i].length(),
                config);
            if (result.status == EDLIB_STATUS_OK){
                int ed = result.editDistance;
                edlibFreeAlignResult(result);
                if (ed == 0){
                    // Be greedy about perfect matches
                    min_ed = ed;
                    min_idx = i;
                    tie = false;
                    break;
                }
                else if (ed > 0){
                    if (min_ed == -1 || ed < min_ed){
                        min_ed = ed;
                        min_idx = i;
                        tie = false;
                    }
                    else if (min_ed != -1 && ed == min_ed){
                        tie = true;
                    }
                }
            }
            else{
                edlibFreeAlignResult(result);
            }
        }
        if (min_ed != -1 && min_idx != -1 && !tie){
            this->mismatches = min_ed;
            return min_idx;
        }
        this->mismatches = -1;
        return -1;
    }
}

