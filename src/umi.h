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
        std::set<unsigned long> nomissing;

        // Store all UMIs
        std::vector<umi> umis;
        // Store group index per UMI
        std::vector<int> umi_groups;
        
        // How many UMI groups?
        int ngrp;

        bool match(const umi& u1, const umi& u2);

    public:
        
        umi_set();

        umi_set(int len);
        
        void set_len(int len);

        // Put an additional UMI in the set
        void add(const umi& umi); 
        
        // Count all unique UMIs 
        int count();
};


#endif
