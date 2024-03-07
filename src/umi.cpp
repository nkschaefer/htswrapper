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
    this->ngrp = 0;
    this->exact = false;
    set_len(len);
}

umi_set::umi_set(){
    this->len = -1;
    this->exact = false;
    this->ngrp = 0;
}

umi_set::umi_set(const umi_set& other){
    if (other.len > 0){
        len = other.len;
        kmers = other.kmers;
        nomissing = other.nomissing;
        umis = other.umis;
        umi_groups = other.umi_groups;
        ngrp = other.ngrp;
        mask = other.mask;
        mask_mask = other.mask_mask;
    }
    exact = other.exact;
}

int umi_set::get_k(){
    if (this->len <= 0){
        return -1;
    }
    // Determine fuzzy matching k value
    // k must be < (bc_length + 1)/2
    int k;
    if (len % 2 == 0){
        // Even barcode length
        k = len/2;
    }
    else{
        k = (len-1)/2;
    }
    return k;
}

void umi_set::set_len(int len){
    this->len = len;
    if (!exact){
        int k = get_k();
        kmers.init(k);
        mask.reset();
        for (int i = 0; i < k*2; ++i){
            mask.set(i);
        }
        mask_mask.reset();
        for (int i = 0; i < k; ++i){
            mask_mask.set(i);
        }
    }
}

void umi_set::print(const umi& u){
    string umi_str = bc2str(u.bits, len);
    for (int i = 0; i < len; ++i){
        if (!u.mask[i]){
            umi_str[i] = 'N';
        }
    }    
    fprintf(stderr, "%s\n", umi_str.c_str());
}

void umi_set::exact_matches_only(bool exact){
    this->exact = exact;
}

void umi_set::add(const umi& to_add){
    if (this->len < 0){
        fprintf(stderr, "ERROR: initialize length first\n");
        exit(1);
    }
    unsigned long ul = to_add.bits.to_ulong();
    if (to_add.missing_multiple){
        // Impossible to be within edit dist 1 of another UMI
        // Count as unique but do not store
        ++ngrp;
        return;
    }
    else if (!to_add.missing_sites && nomissing.count(ul) > 0 &&
        nomissing[ul] == to_add.mask.to_ulong()){
        // Has exact match; no need to add
        return;
    }
    
    bool has_match = false;
    bool add_new = true;
    int group_num = -1;
    
    if (!exact){
        // Look for a match with edit dist <= 1
        vector<int> matches(umis.size(), 0);
        set<unsigned long> inds_check;
        
        int k = get_k();

        for (int start = 0; start <= len-k; ++start){
            // Test for number of valid (non-N) bases in this k-mer   
            if ((mask_mask & (to_add.mask >> start)).count() == k){
                
                // No skipped bases in this k-mer.
                // Get bitset representation of current k-mer
                unsigned long kmer = (mask & (to_add.bits >> 2*start)).to_ulong();
                
                // In bc.cpp, kmer_lookup is used to store unsigned long representation of 
                // cell barcodes. Here, the unsigned longs will (somewhat wastefully)
                // represent indices into the list of UMIs. 
                unsigned long list_idx;
                while(kmers.lookup(kmer, list_idx)){
                    matches[list_idx]++;
                    if (matches[list_idx] >= len - 2*k + 1){
                        // Pass
                        inds_check.insert(list_idx);
                    } 
                } 
            }
        }
        
        
        
        for (set<unsigned long>::iterator i = inds_check.begin(); i != inds_check.end(); ++i){
            unsigned long idx = *i;
            // Compare the two.
            // Greedy: just go with the first match.
            if (match(umis[idx], to_add)){
                has_match = true;
                group_num = umi_groups[idx];
                if (umis[idx].missing_sites && !to_add.missing_sites){
                    // Use the new k-mer to fill in the holes in the previous one.
                    for (int i = 0; i < len; ++i){
                        if (!umis[idx].mask.test(i) && to_add.mask.test(i)){
                            umis[idx].bits[2*i] = to_add.bits[2*i];
                            umis[idx].bits[2*i+1] = to_add.bits[2*i+1];
                            umis[idx].mask.set(i);

                            // Add missing k-mer lookup stuff
                            int first_idx = i-k+1;
                            if (first_idx < 0){
                                first_idx = 0;
                            }
                            int last_idx = i;
                            if (last_idx > len-k){
                                last_idx = len-k;
                            }
                            for (int start = first_idx; start <= last_idx; ++start){
                                unsigned long kmer = (mask & (to_add.bits >> 2*start)).to_ulong();
                                kmers.insert(kmer, idx);
                            }
                            break;
                        }
                    }
                    add_new = false;
                }
                else if (to_add.missing_sites && !umis[idx].missing_sites){
                    // The old one is more informative.
                    add_new = false;
                }
                // Stop looking for matches here
                break;
            }
        }
    }

    if (add_new){
        if (group_num == -1){
            // Need to add a new group.
            group_num = ngrp;
            ++ngrp;
        }
        unsigned long new_umi_idx = umis.size();
        umis.push_back(to_add);
        umi_groups.push_back(group_num);
        if (exact || !to_add.missing_sites){
            this->nomissing.emplace(ul, to_add.mask.to_ulong());
        }
        
        if (!exact){
            int k = get_k();
            // Add all k-mers to lookup table.
            for (int start = 0; start <= len-k; ++start){
                // Test for number of valid (non-N) bases in this k-mer   
                if ((mask_mask & (to_add.mask >> start)).count() == k){
                    
                    // No skipped bases in this k-mer.
                    // Get bitset representation of current k-mer
                    unsigned long kmer = (mask & (to_add.bits >> 2*start)).to_ulong();
                    kmers.insert(kmer, new_umi_idx);
                }
            } 
        }
    }
}

int umi_set::count(){
    return ngrp;
}

