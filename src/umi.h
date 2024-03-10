#ifndef _HTSWRAPPER_UMI_H
#define _HTSWRAPPER_UMI_H
#include <utility>
#include <cstdlib>
#include <set>
#include <unordered_set>
#include <bitset>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>
#include "bc.h"
#include "umi.h"
#include "robin_hood/robin_hood.h"

/**
 * This file contains a class designed to help collapse and count UMIs
 */

// Allow UMIs to be up to the length of a barcode.

// UMI masks will mark which sites are visible/non-N, and they only
// need to be half as many bits for that reason.
typedef std::bitset<BC_LENX2/2> umi_mask;

class umi{  
    private:
        void str2umi(const char* str, int len);
    public:
        bc bits;
        umi_mask mask;
        bool missing_sites;
        bool missing_multiple;
        umi(const char* str, int len);
        umi(const umi& other);
};

class umi_set{
   
    private:
        int len;
        
        // Store hash of every UMI with no missing sites 
        robin_hood::unordered_map<unsigned long, unsigned long> nomissing;
        
        // Store k-mers for fuzzy matching
        kmer_lookup kmers;

        // Store all UMIs
        std::vector<umi> umis;
        // Store group index per UMI
        std::vector<int> umi_groups;
        
        // How many UMI groups?
        int ngrp;
        
        bool match(const umi& u1, const umi& u2);
        
        // To use for scanning for k-mers
        std::bitset<BC_LENX2> mask;
        std::bitset<BC_LENX2/2> mask_mask;
        
        // Compute optimal k-mer length given UMI length
        int get_k();
        
        // Do matches have to be exact to be collapsed? (faster)
        bool exact;

    public:
        
        umi_set();

        umi_set(int len);
        
        umi_set(const umi_set& other);

        void exact_matches_only(bool exact);

        void set_len(int len);

        // Put an additional UMI in the set
        // Return whether or not it matched a pre-existing one.
        bool add(const umi& umi); 
        
        // Count all unique UMIs 
        int count();
        
        void print(const umi& umi);
};


#endif
