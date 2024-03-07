#include <utility>
#include <cstdlib>
#include <string>
#include <cstring>
#include <vector>
#include <math.h>
#include "edlib/edlib.h"
#include "seq_fuzzy_match.h"

using namespace std;

seq_fuzzy_match::seq_fuzzy_match(vector<string>& seqs, int edit_dist, bool penalize_n, bool global){
    
    // Tell edlib (optionally) to count N - any base as a match 
    eqs[0] = {'N', 'A'};
    eqs[1] = {'N', 'C'};
    eqs[2] = {'N', 'G'};
    eqs[3] = {'N', 'T'};

    this->global = global;
    this->mismatches = -1;

    if (global){
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
            this->config = edlibNewAlignConfig(edit_dist, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0);
        }
        else{
            this->config = edlibNewAlignConfig(edit_dist, EDLIB_MODE_HW, EDLIB_TASK_LOC, &eqs[0], 4);
        }
    }
    else{
        this->seqs = seqs;
        if (penalize_n){
            this->config = edlibNewAlignConfig(edit_dist, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
        }
        else{
            this->config = edlibNewAlignConfig(edit_dist, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, &eqs[0], 4);
        }
    }
}

int seq_fuzzy_match::match(const char* seq){
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

