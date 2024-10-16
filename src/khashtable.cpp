#include <zlib.h>
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <set>
#include <deque>
#include <math.h>
#include <cstdlib>
#include "bc.h"
#include "robin_hood/robin_hood.h"
#include "khashtable.h"

using namespace std;

bool fill_chunk(std::bitset<KHASH_CHUNK_SIZE>& this_chunk, 
    const char* seq, int start, int len, bool rc){
    if (rc){
        int bcbit = 0;
        for (int base_idx = start + len - 1; base_idx >= start; --base_idx){
            switch(seq[base_idx]){
                case 'T':
                    this_chunk.reset(bcbit);
                    this_chunk.reset(bcbit+1);
                    break;
                case 'G':
                    this_chunk.set(bcbit);
                    this_chunk.reset(bcbit+1);
                    break;
                case 'C':
                    this_chunk.reset(bcbit);
                    this_chunk.set(bcbit+1);
                    break;
                case 'A':
                    this_chunk.set(bcbit);
                    this_chunk.set(bcbit+1);
                    break;
                default:
                    return false;
                    break;
            }
            bcbit += 2;
        }
    }
    else{
        int bcbit = 0;
        for (int base_idx = start; base_idx < start+len; ++base_idx){
            switch(seq[base_idx]){
                case 'A':
                    this_chunk.reset(bcbit);
                    this_chunk.reset(bcbit+1);
                    break;
                case 'C':
                    this_chunk.set(bcbit);
                    this_chunk.reset(bcbit+1);
                    break;
                case 'G':
                    this_chunk.reset(bcbit);
                    this_chunk.set(bcbit+1);
                    break;
                case 'T':
                    this_chunk.set(bcbit);
                    this_chunk.set(bcbit+1);
                    break;
                default:
                    return false;
                    break;
            }
            bcbit += 2;
        }
    }
    return true;
}

void print_chunk(std::bitset<KHASH_CHUNK_SIZE>& chunk, bool rc){
    for (int idx = 0; idx < KHASH_CHUNK_SIZE/2; ++idx){
        int i = idx;
        if (rc){
            i = KHASH_CHUNK_SIZE/2 - 1 - idx;
        }
        if (!chunk.test(2*i) && !chunk.test(2*i+1)){
            if (rc){
                fprintf(stderr, "T");
            }
            else{
                fprintf(stderr, "A");
            }
        }
        else if (chunk.test(2*i) && !chunk.test(2*i+1)){
            if (rc){
                fprintf(stderr, "G");
            }
            else{
                fprintf(stderr, "C");
            }

        }
        else if (!chunk.test(2*i) && chunk.test(2*i+1)){
            if (rc){
                fprintf(stderr, "C");
            }
            else{
                fprintf(stderr, "G");
            }
        }
        else{
            if (rc){
                fprintf(stderr, "A");
            }
            else{
                fprintf(stderr, "T");
            }
        }
    }
    fprintf(stderr, "\n");
}

void khashkey::init(int k){
    chunksize = KHASH_CHUNK_SIZE/2;
    n_chunks = (int)ceil((double)k/(double)chunksize);
   
    last_chunksize = k - (n_chunks-1)*chunksize;
   
    mask = ((1UL << last_chunksize*2) - 1);
    
    this->k = k;
    this->n_pos = -1;
    for (int i = 0; i < n_chunks; ++i){
        std::bitset<KHASH_CHUNK_SIZE> chunk;
        chunks.push_back(chunk);
        chunks_rev.push_back(chunk);
    }
    need_rebuild = true;
}

khashkey::khashkey(int k){
    init(k);
}

khashkey::khashkey(const char* seq, int k){
    init(k);
    need_rebuild = !rebuild(seq);    
}

void khashkey::reset(){
    need_rebuild = true;
    n_pos = -1;
}

bool khashkey::rebuild(const char* seq){
    bc this_bc;
    for (int chunk = 0; chunk < n_chunks; ++chunk){
        int start = chunk*chunksize;
        int len = chunksize;
        if (chunk == n_chunks -1){
            len = last_chunksize;
        }
        
        int start_rev = k-((chunk+1)*chunksize);
        int len_rev = chunksize;
        if (start_rev < 0){
            len_rev += start_rev;
            start_rev = 0;
        }
        
        if (fill_chunk(chunks_rev[chunk], seq, start_rev, len_rev, true)){
            // all good
        }
        else{
            // Find the N. We'll see the furthest-right N this way, which is what we want.
            for (int i = start_rev + len_rev-1; i >= start_rev; --i){
                if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T'){
                    n_pos = i;
                    need_rebuild = true;
                    return false;
                }
            }
        }
        if (!fill_chunk(chunks[chunk], seq, start, len, false)){
            // This shouldn't happen (we should catch it above)
            need_rebuild = true;
            return false;
        }
        if (chunk < n_chunks-1){
            //chunks[chunk] &= mask;
            //chunks_rev[chunk] &= mask;
        }
        else{
            chunks[chunk] &= mask;
            chunks_rev[chunk] &= mask;
        }
    }

    n_pos = -1;
    return true;
}

void khashkey::print(bool rc){
    for (int i = 0; i < chunks.size(); ++i){
        int size = chunksize;
        if (i == chunks.size()-1){
            size = last_chunksize;
        }
        for (int x = 0; x < i; ++x){
            fprintf(stderr, " ");
        }
        for (int j = 0; j < size; ++j){
            if (rc){
                if (chunks_rev[i].test(j*2)){
                    if (chunks_rev[i].test(j*2+1)){
                        fprintf(stderr, "T");
                    }
                    else{
                        fprintf(stderr, "C");
                    }
                }
                else{
                    if(chunks_rev[i].test(j*2+1)){
                        fprintf(stderr, "G");
                    }
                    else{
                        fprintf(stderr, "A");
                    }
                }
            }
            else{
                if (chunks[i].test(j*2)){
                    if (chunks[i].test(j*2 + 1)){
                        fprintf(stderr, "T");
                    }
                    else{
                        fprintf(stderr, "C");
                    }
                }
                else{
                    if (chunks[i].test(j*2+1)){
                        fprintf(stderr, "G");
                    }
                    else{
                        fprintf(stderr, "A");
                    }
                }
            }
        }
        fprintf(stderr, "\n");
    }
}

bool khashkey::advance(char c){
    // A = 00
    // C = 10
    // G = 01
    // T = 11
    
    if (c == 'N' || c == 'n'){
        // Assume char c is at the current right-most end of seq
        n_pos = k-1;
        need_rebuild = true;
        return false;
    }
    else{
        for (int ci = 0; ci < n_chunks; ++ci){
            int ci_rev = n_chunks-1-ci;
            
            bool not_full = (ci == n_chunks-1);
            bool not_full_rev = (ci_rev == n_chunks-1);

            if (ci_rev == 0){
                chunks_rev[ci_rev] = (chunks_rev[ci_rev] << 2);
                if (not_full_rev){
                    chunks_rev[ci_rev] &= mask;
                }
                switch(c){
                    case 'A':
                        chunks_rev[ci_rev].set(0);
                        chunks_rev[ci_rev].set(1);
                        break;
                    case 'C':
                        chunks_rev[ci_rev].reset(0);
                        chunks_rev[ci_rev].set(1);
                        break;
                    case 'G':
                        chunks_rev[ci_rev].set(0);
                        chunks_rev[ci_rev].reset(1);
                        break;
                    case 'T':
                        chunks_rev[ci_rev].reset(0);
                        chunks_rev[ci_rev].reset(1);
                        break;
                    default:
                        n_pos = k-1;
                        need_rebuild = true;
                        return false;
                        break;
                }
            }
            else{
                chunks_rev[ci_rev] <<= 2;
                if (not_full_rev){
                    chunks_rev[ci_rev] &= mask;
                }

                if (chunks_rev[ci_rev-1].test((chunksize-1)*2)){
                    chunks_rev[ci_rev].set(0);
                }
                else{
                    chunks_rev[ci_rev].reset(0);
                }
                if (chunks_rev[ci_rev-1].test((chunksize-1)*2+1)){
                    chunks_rev[ci_rev].set(1);
                }
                else{
                    chunks_rev[ci_rev].reset(1);
                }
            }
            
            if (ci == n_chunks-1){
                chunks[ci] >>= 2;
                if (not_full){
                    chunks[ci] &= mask;
                }
                switch(c){
                    case 'A':
                        chunks[ci].reset((last_chunksize-1)*2);
                        chunks[ci].reset((last_chunksize-1)*2+1);
                        break;
                    case 'C':
                        chunks[ci].set((last_chunksize-1)*2);
                        chunks[ci].reset((last_chunksize-1)*2+1);
                        break;
                    case 'G':
                        chunks[ci].reset((last_chunksize-1)*2);
                        chunks[ci].set((last_chunksize-1)*2+1);
                        break;
                    case 'T':
                        chunks[ci].set((last_chunksize-1)*2);
                        chunks[ci].set((last_chunksize-1)*2+1);
                        break;
                    default:
                        n_pos = k-1;
                        need_rebuild = true;
                        return false;
                }    
            }
            else{
                chunks[ci] >>= 2;
                if (not_full){
                    chunks[ci] &= mask;
                }
                if (chunks[ci+1].test(0)){
                    chunks[ci].set((chunksize-1)*2);
                }
                else{
                    chunks[ci].reset((chunksize-1)*2);
                }
                if (chunks[ci+1].test(1)){
                    chunks[ci].set((chunksize-1)*2+1);
                }
                else{
                    chunks[ci].reset((chunksize-1)*2+1);
                }
                
            }
        }
    }
    return true;
}

bool khashkey::scan_kmers(const char* seq, int len, int& pos){
    while (pos < len-k + 1){
        if (n_pos > 0){
            n_pos--;
        }
        else{
            if (need_rebuild){
                if (rebuild(seq + pos)){
                    // Saw a valid k-mer.
                    need_rebuild = false;
                    ++pos;
                    return true;
                }
            }
            else{
                // Try to advance to next base.
                if (advance(seq[pos + k - 1])){
                    ++pos;
                    return true;
                }
                else{
                    need_rebuild = true;
                }
            }
        }
        ++pos;
    } 
    return false;
}

bool khashkey::scan_kmers(string& seq, int& pos){
    return scan_kmers(seq.c_str(), seq.length(), pos);
}
