#ifndef _HTSWRAPPER_KMHASH_H
#define _HTSWRAPPER_KMHASH_H
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
#include <math.h>
#include <unordered_map>
#include <set>
#include <deque>
#include <cstdlib>
#include <bitset>
#include "emhash/hash_table8.hpp"
#include "robin_hood/robin_hood.h"
#define KHASH_CHUNK_SIZE (64)

class khashkey{
    private:
        std::bitset<KHASH_CHUNK_SIZE> mask;
         
        //bc mask;
        //bc mask_end;
        int n_chunks;
        int chunksize;
        int last_chunksize;
        int k;
        bool need_rebuild;
        int n_pos;
        void init(int k);
        bool rebuild(const char* seq);
        bool advance(const char c);
    public:
        std::vector<std::bitset<KHASH_CHUNK_SIZE> > chunks;
        std::vector<std::bitset<KHASH_CHUNK_SIZE> > chunks_rev;

        khashkey(const char* seq, int k);
        khashkey(int k);
        void reset();
        void print(bool rc);
        bool scan_kmers(const char* seq, int len, int& pos);
        bool scan_kmers(std::string& seq, int& pos);
};

template<typename datatype> 
struct khashnode;

template<typename datatype>
using knode_map = robin_hood::unordered_flat_map<unsigned long long, khashnode<datatype>* >;

template<typename datatype> struct khashnode{
    datatype data;
    bool end;
    knode_map<datatype>* map;

    khashnode(bool e=false){
        this->end = e;
        if (e){
            this->map = NULL;
        }
        else{
            this->map = new knode_map<datatype>;
        }
    };
    
    ~khashnode(){
        if (this->end){
            // Don't do anything weird
        }
        else if (this->map != NULL){
            for (typename knode_map<datatype>::iterator m = this->map->begin(); 
                m != this->map->end(); ++m){        
                delete m->second;
            }
            delete this->map;
        }   
    };
    void set(datatype d){
        this->data = d;   
    };
};



bool fill_chunk(std::bitset<KHASH_CHUNK_SIZE>& this_chunk, const char* seq, int start, int len, bool rc);
void print_chunk(std::bitset<KHASH_CHUNK_SIZE>& this_chunk, bool rc=false);

template<typename datatype>
class khashtable{
    private:
        int n_chunks;
        int chunk_size;
        int k;
        int last_chunksize;
        unsigned long array_size;
        knode_map<datatype> elts;
    public:
        
        khashtable(int k){
            
            chunk_size = KHASH_CHUNK_SIZE/2;
            n_chunks = (int)ceil((double)k/(double)chunk_size);
            
            // Error check
            if (n_chunks >= k){
                fprintf(stderr, "ERROR: num chunks must be < k\n");
                exit(1);
            }
            if (chunk_size > 32){
                fprintf(stderr, "ERROR: maximum chunk size is 16 bp\n");
                fprintf(stderr, "please re-compile with -DKHASH_CHUNK_SIZE=n\n");
                fprintf(stderr, "where n <= 16\n");
                exit(1);
            }
            last_chunksize = k - (n_chunks-1)*chunk_size;
            this->k = k;
        };

        ~khashtable(){
            for (typename knode_map<datatype>::iterator x = this->elts.begin(); 
                x != this->elts.end(); ++x){
                delete x->second;
            }
            this->elts.clear();   
            
        };
        
        void clear(){
            for (typename knode_map<datatype>::iterator x = this->elts.begin(); 
                x != this->elts.end(); ++x){
                delete x->second;
            }
            this->elts.clear();
            
        }

        bool is_rc_first(const char* kmer){
            for (int i = 0; i < k; ++i){
                char char_f = kmer[i];
                char char_rc = 'N';
                switch(kmer[k-1-i]){
                    case 'A':
                        char_rc = 'T';
                        break;
                    case 'C':
                        char_rc = 'G';
                        break;
                    case 'G':
                        char_rc = 'C';
                        break;
                    case 'T':
                        char_rc = 'A';
                        break;
                    default:
                        char_rc = 'N';
                        break;
                }
                if (char_f == 'N' || char_rc == 'N'){
                    // No matter; neither will be in the lookup table.
                    return false;
                }
                if (char_rc < char_f){
                    return true;
                }
                else if (char_f < char_rc){
                    return false;
                }   
            }
            return false;
        };

        bool add(const char* seq, datatype& data){
            unsigned long ul;
            std::deque<unsigned long long> keys_path;
            
            bool rc = is_rc_first(seq); 
            for (int ci = 0; ci < n_chunks; ++ci){
                int start;
                int len;
                if (rc){
                    start = k-((ci+1)*chunk_size);
                    len = chunk_size;
                    
                    if (start < 0){
                        len += start;
                        start = 0;
                    }
                }
                else{
                    start = ci*chunk_size;
                    len = chunk_size;
                    if (ci == n_chunks-1){
                        len = last_chunksize;
                    }
               }

                std::bitset<KHASH_CHUNK_SIZE> this_chunk;
                
                bool success = fill_chunk(this_chunk, seq, start, len, rc);
                
                if (success){
                    unsigned long long ul = this_chunk.to_ullong();
                    keys_path.push_back(ul);
                }
                else{
                    // Do not create.
                    return false;
                }
            }
            
            khashnode<datatype>* prev = NULL;
            for (int i = 0; i < keys_path.size(); ++i){
                khashnode<datatype>* n;
                unsigned long long ul = keys_path[i];
                if (prev == NULL){
                    if (this->elts.count(ul) > 0){
                        n = this->elts[ul];
                    }
                    else{
                        if (i == keys_path.size()-1){
                            n = new khashnode<datatype>(true);
                            n->set(data);
                        }
                        else{
                            n = new khashnode<datatype>;
                        }
                        this->elts.emplace(ul, n);
                    }
                    
                }
                else if (i == keys_path.size()-1){
                    // Final node - should not be in map yet
                    if (prev->map->count(ul) > 0){
                        // Do nothing - something went wrong (duplicate keys)
                    }
                    else{
                        n = new khashnode<datatype>(true);
                        n->set(data);
                        prev->map->emplace(ul, n);
                        //prev->map[ul] = n;
                    }
                }
                else{
                    
                    if (prev->map->count(ul) > 0){
                        n = (*prev->map)[ul];
                    }
                    else{
                        n = new khashnode<datatype>;
                        prev->map->emplace(ul, n);
                    }
                    
                }
                prev = n;
            }
            return true;   
        };
        
        bool lookup(const khashkey& key, datatype& val){
            khashnode<datatype>* n_f = NULL;
            khashnode<datatype>* n_r = NULL;
            
            if (elts.count(key.chunks[0].to_ullong()) > 0){
                n_f = elts[key.chunks[0].to_ullong()];
            }
            if (elts.count(key.chunks_rev[0].to_ullong()) > 0){
                n_r = elts[key.chunks_rev[0].to_ullong()];
            }
            
            if (n_f == NULL && n_r == NULL){
                return false;
            }
            for (int i = 1; i < key.chunks.size(); ++i){
                if (n_f != NULL){
                    
                    if (n_f->map->count(key.chunks[i].to_ullong()) > 0){
                        n_f = (*n_f->map)[key.chunks[i].to_ullong()];
                    }
                    else{
                        n_f = NULL;
                    }
                    
                }
                if (n_r != NULL){
                    
                    if (n_r->map->count(key.chunks_rev[i].to_ullong()) > 0){
                        n_r = (*n_r->map)[key.chunks_rev[i].to_ullong()];
                    }
                    else{
                        n_r = NULL;
                    }
                    
                }
                if (n_f == NULL && n_r == NULL){
                    return false;
                }
            }
            if (n_f != NULL && n_f->end){
                val = n_f->data;
                return true;
            }
            if (n_r != NULL && n_r->end){
                val = n_r->data;
                return true;
            }
            return false;
        };

        bool lookup(const char* key, datatype& val){
            khashkey hkey(key, k);
            return lookup(hkey, val);
        };
        
        double lf(){
            return this->elts.load_factor();
        }     
};

#endif

