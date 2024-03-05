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

using namespace std;

void umi::str2umi(const char* str, int len){
    bits.reset();
    mask.reset();
    int bcbit = 0;
    for (int i = 0; i < len; ++i){
        switch(str[i]){
            case 'A':
            case 'a':
                mask.set(i);
                bcbit += 2;
                break;
            case 'C':
            case 'c':
                mask.set(i);
                bits.set(bcbit);
                bcbit += 2;
                break;
            case 'G':
            case 'g':
                mask.set(i);
                bits.set(bcbit+1);
                bcbit += 2;
                break;
            case 'T':
            case 't':
                mask.set(i);
                bits.set(bcbit);
                bits.set(bcbit+1);
                bcbit += 2;
                break;
            case 'n':
            case 'N':
                // Do not set mask or bits
                bcbit += 2;
                break;
        }
    }
    this->missing_sites = mask.count() < len;
    if (this->missing_sites){
        this->missing_multiple = mask.count() < len-1;
    }
    else{
        this->missing_multiple = false;
    }
}

umi::umi(const char* str, int len){
    this->str2umi(str, len);    
}   

umi::umi(const umi& other){
    this->bits = other.bits;
    this->mask = other.mask;
    this->missing_sites = other.missing_sites;
    this->missing_multiple = other.missing_multiple;
}

/**
 * Returns true if the two UMIs are within edit distance 1 of
 * each other and false otherwise.
 */
bool umi_set::match(const umi& u1, const umi& u2){
    if (u1.missing_multiple || u2.missing_multiple){
        // Can't be edit distance <= 1 
        return false;
    }
    if (!u1.missing_sites && !u2.missing_sites && u1.bits.to_ulong() == u2.bits.to_ulong()){
        return true;
    }
    else{
        int ed = 0;
        for (int i = 0; i < len; ++i){
            if (u1.mask.test(i) && u2.mask.test(i)){
                if (u1.bits.test(i*2) == u2.bits.test(i*2) && 
                    u1.bits.test(i*2+1) == u2.bits.test(i*2+1)){
                    // Pass
                }
                else{
                    ++ed;
                    if (ed > 1){
                        break;
                    }
                }
            }
            else{
                // Count N as mismatch.
                ++ed;
                if (ed > 1){
                    break;
                }
            }
        }
        if (ed <= 1){
            return true;
        }
        else{
            return false;
        }
    } 
}

umi_set::umi_set(int len){
    this->len = len;
    this->ngrp = 0;
}

umi_set::umi_set(){
    this->len = -1;
    this->ngrp = 0;
}

void umi_set::set_len(int len){
    this->len = len;
}

void umi_set::add(const umi& to_add){
    if (this->len < 0){
        fprintf(stderr, "ERROR: initialize length first\n");
        exit(1);
    }
    if (to_add.missing_multiple){
        // Impossible to be within edit dist 1 of another UMI
        return;
    }
    else if (!to_add.missing_sites && nomissing.find(to_add.bits.to_ulong()) != nomissing.end()){
        // Has exact match; no need to add
        return;
    }
    // Look for a match with edit dist <= 1
    int grp_idx = -1;
    for (int i = 0; i < this->umis.size(); ++i){
        if (match(umis[i], to_add)){
            grp_idx = umi_groups[i];
            break;
        }
    }
    if (grp_idx == -1){
        // Need to add a new group.
        grp_idx = ngrp;
        ++ngrp;
    }
    umis.push_back(to_add);
    umi_groups.push_back(grp_idx);
    if (!to_add.missing_sites){
        this->nomissing.insert(to_add.bits.to_ulong());
    }
}

int umi_set::count(){
    return ngrp;
}

