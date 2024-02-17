#include <cstdlib>
#include <string.h>
#include <utility>
#include <bitset>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <zlib.h>
#include <limits.h>
#include "gzreader.h"
#include "bc.h"

using namespace std;

/**
 * Convert a DNA sequence from string to bitset
 * representation
 */
bool str2bc(const char* str, bc& this_bc){
    this_bc.reset();
    int bcbit = 0;
    for (int i = 0; i < BC_LENX2/2; ++i){
        switch(str[i]){
            case 'A':
            case 'a':
                bcbit += 2;
                break;
            case 'C':
            case 'c':
                this_bc.set(bcbit);
                bcbit += 2;
                break;
            case 'G':
            case 'g':
                this_bc.set(bcbit+1);
                bcbit += 2;
                break;
            case 'T':
            case 't':
                this_bc.set(bcbit);
                this_bc.set(bcbit+1);
                bcbit += 2;
                break;
            default:
                return false;
                break;
        }
    }
    return true;
}

bool str2kmers(const char* str, int& start, kmer& cur_kmer){
    if (start > BC_LENX2/2 - KX2/2){
        return false;
    }
    cur_kmer.reset();
    int kmerbit = 0;
    for (int i = start; i < start + KX2/2; ++i){
        switch(str[i]){
            case 'A':
            case 'a':
                kmerbit += 2;
                break;
            case 'C':
            case 'c':
                cur_kmer.set(kmerbit);
                kmerbit += 2;
                break;
            case 'G':
            case 'g':
                cur_kmer.set(kmerbit+1);
                kmerbit += 2;
                break;
            case 'T':
            case 't':
                cur_kmer.set(kmerbit);
                cur_kmer.set(kmerbit+1);
                kmerbit += 2;
                break;
            default:
                start += KX2/2;
                return false;
                break;
        }
    }
    start++;
    return true;
}

/**
 * Convert a DNA sequence from string to bitset format,
 * in reverse complement orientation
 */
bool str2bc_rc(const char* str, bc& this_bc){
    this_bc.reset();
    int bcbit = 0;
    for (int i = BC_LENX2/2-1; i >= 0; --i){
        switch(str[i]){
            case 'T':
            case 't':
                bcbit += 2;
                break;
            case 'G':
            case 'g':
                this_bc.set(bcbit);
                bcbit += 2;
                break;
            case 'C':
            case 'c':
                this_bc.set(bcbit+1);
                bcbit += 2;
                break;
            case 'A':
            case 'a':
                this_bc.set(bcbit);
                this_bc.set(bcbit+1);
                bcbit += 2;
                break;
            default:
                return false;
                break;
        }
    }
    return true;
}

/**
 * Convert a DNA sequence from string to bitset format,
 * in reverse complement orientation
 */
bool str2kmers_rc(const char* str, int& start, kmer& cur_kmer){
    cur_kmer.reset();
    int kmerbit = 0;
    if (start == -1){
        // initialize start position
        start = BC_LENX2/2-1;
    }
    else if (start < KX2/2-1){
        return false;
    }
    for (int i = start; i >= start-KX2/2+1; --i){
        switch(str[i]){
            case 'T':
            case 't':
                kmerbit += 2;
                break;
            case 'G':
            case 'g':
                cur_kmer.set(kmerbit);
                kmerbit += 2;
                break;
            case 'C':
            case 'c':
                cur_kmer.set(kmerbit+1);
                kmerbit += 2;
                break;
            case 'A':
            case 'a':
                cur_kmer.set(kmerbit);
                cur_kmer.set(kmerbit+1);
                kmerbit += 2;
                break;
            default:
                return false;
                break;
        }
    }
    start--;
    return true;
}

/**
 * Convert a bitset representation of a DNA sequence
 * to a string representation
 */
string bc2str(bc& this_bc){
    char strbuf[BC_LENX2/2+1];
    strbuf[BC_LENX2/2] = '\0';
    for (int i = 0; i < BC_LENX2/2; ++i){
        if (this_bc.test(i*2)){
            if (this_bc.test(i*2+1)){
                strbuf[i] = 'T';
            }
            else{
                strbuf[i] = 'C';
            }
        }
        else{
            if (this_bc.test(i*2 + 1)){
                strbuf[i] = 'G';    
            }
            else{
                strbuf[i] = 'A';
            }
        }
    }
    return string(strbuf);
}

string bc2str(unsigned long ul){
    bc as_bc(ul);
    return bc2str(as_bc);
}

/**
 * Convert a bitset representation of a DNA sequence
 * to a string representation
 */
string kmer2str(kmer& this_bc){
    char strbuf[KX2/2+1];
    strbuf[KX2/2] = '\0';
    for (int i = 0; i < KX2/2; ++i){
        if (this_bc.test(i*2)){
            if (this_bc.test(i*2+1)){
                strbuf[i] = 'T';
            }
            else{
                strbuf[i] = 'C';
            }
        }
        else{
            if (this_bc.test(i*2 + 1)){
                strbuf[i] = 'G';    
            }
            else{
                strbuf[i] = 'A';
            }
        }
    }
    return string(strbuf);
}

/**
 * Convert a bitset representation of a DNA sequence
 * to a string representation, in reverse complement
 * orientation.
 */
string bc2str_rc(bc& this_bc){
    char strbuf[BC_LENX2/2+1];
    strbuf[BC_LENX2/2] = '\0';
    int buf_idx = 0;
    for (int i = BC_LENX2/2-1; i >= 0; --i){
        if (this_bc.test(i*2)){
            if (this_bc.test(i*2+1)){
                strbuf[buf_idx] = 'A';
            }
            else{
                strbuf[buf_idx] = 'G';
            }
        }
        else{
            if (this_bc.test(i*2 + 1)){
                strbuf[buf_idx] = 'C';    
            }
            else{
                strbuf[buf_idx] = 'T';
            }
        }
        buf_idx += 1;
    }
    return string(strbuf);
}

string bc2str_rc(unsigned long ul){
    bc as_bits(ul);
    return bc2str_rc(as_bits);
}

/**
 * Convert a bitset representation of a DNA sequence
 * to a string representation, in reverse complement
 * orientation.
 */
string kmer2str_rc(kmer& this_bc){
    char strbuf[KX2/2+1];
    strbuf[KX2/2] = '\0';
    int buf_idx = 0;
    for (int i = KX2/2-1; i >= 0; --i){
        if (this_bc.test(i*2)){
            if (this_bc.test(i*2+1)){
                strbuf[buf_idx] = 'A';
            }
            else{
                strbuf[buf_idx] = 'G';
            }
        }
        else{
            if (this_bc.test(i*2 + 1)){
                strbuf[buf_idx] = 'C';    
            }
            else{
                strbuf[buf_idx] = 'T';
            }
        }
        buf_idx += 1;
    }
    return string(strbuf);
}

/**
 * Helper function for below: trims extraneous stuff off the 
 * end of a cell barcode string and returns an unsigned long
 * representation.
 *
 */
unsigned long bc_ul(string& barcode){
    // Look for any extra stuff (i.e. library IDs)
    // common convention is to affix these to the beginning or end, with
    // - or _ as a separator.
    int sep_pos = -1;
    // Look backward, since 10X puts them on the end.
    for (int i = barcode.length()-1; i >= 0; --i){
        if (barcode[i] == '-' || barcode[i] == '_'){
            sep_pos = i;
            break;
        }
    }
    if (sep_pos != -1){
        if (barcode.length()-sep_pos < sep_pos){
            // End of read
            barcode = barcode.substr(0, sep_pos);
        }
        else{
            // Beginning of read.
            barcode = barcode.substr(sep_pos+1, barcode.length()-sep_pos-1);
        }
    }
    bc as_bitset;
    str2bc(barcode.c_str(), as_bitset);
    return as_bitset.to_ulong();
}

/**
 * Same as above, but if working with a character buffer,
 * avoid copying into a std::string for speed
 */
unsigned long bc_ul(char* barcode){
    // Look for any extra stuff (i.e. library IDs)
    // common convention is to affix these to the beginning or end, with
    // - or _ as a separator.
    int sep_pos = -1;
    // Look backward, since 10X puts them on the end.
    for (int i = strlen(barcode)-1; i >= 0; --i){
        if (barcode[i] == '-' || barcode[i] == '_'){
            sep_pos = i;
            break;
        }
    }
    if (sep_pos != -1){
        if (strlen(barcode)-sep_pos < sep_pos){
            // End of read
            barcode[sep_pos] = '\0';
        }
        else{
            // Beginning of read.
            int len = strlen(barcode)-sep_pos-1;
            memmove(&barcode[0], &barcode[sep_pos+1], len);
            barcode[len] = '\0';
        }
    }
    bc as_bitset;
    str2bc(barcode, as_bitset);
    return as_bitset.to_ulong();
}

/**
 * Load a filtered list of cell barcodes, gzipped or not.
 * 
 * Trims any trailing stuff (i.e. -1) off the end of the actual
 * barcode sequence, compresses the barcode, and interprets it 
 * as an unsigned long.
 *
 */
void parse_barcode_file(string& filename, set<unsigned long>& cell_barcodes){
   
    gzreader reader(filename);
    while (reader.next()){
        cell_barcodes.insert(bc_ul(reader.line));
    } 
    
    fprintf(stderr, "Read %ld barcodes from file\n", cell_barcodes.size());
}

kmer_lookup::kmer_lookup(){
    // Determine max # k-mers
    kmer test;
    for (int i = 0; i < KX2; ++i){
        test.set(i);
    }
    unsigned long nk = test.to_ulong();
    n_kmers = nk;
    table = (kmer_lookup_node**)malloc(nk*sizeof(kmer_lookup_node*));
    table_last = (kmer_lookup_node**)malloc(nk*sizeof(kmer_lookup_node*));
    for (int i = 0; i < nk; ++i){
        table[i] = NULL;
        table_last[i] = NULL;
    }
    curnode = NULL;
}

void kmer_lookup::free_members(kmer_lookup_node* n){
    if (n->next != NULL){
        free_members(n->next);
        n->next = NULL;
    }
    delete n;
}

kmer_lookup::~kmer_lookup(){
    // Free all nodes
    for (size_t i = 0; i < n_kmers; ++i){
        if (table[i] != NULL){
            table_last[i] = NULL;
            free_members(table[i]);
        }
    }
    // Free tables
    free(table);
    free(table_last);
}

void kmer_lookup::insert(kmer& k, unsigned long barcode){
    kmer_lookup_node* n = new kmer_lookup_node(barcode);
    unsigned long idx = k.to_ulong();
    if (table[idx] == NULL){
        table[idx] = n;
        table_last[idx] = n;
    }
    else{
        kmer_lookup_node* bucket = table_last[idx];
        bucket->next = n;
        table_last[idx] = n;
    }
}

bool kmer_lookup::lookup(kmer& k, unsigned long& barcode){
    if (curnode != NULL){
        // Assume we're still iterating
        if (curnode->next == NULL){
            // Done.
            curnode = NULL;
            return false;
        }
        else{
            // Get next in list
            curnode = curnode->next;
            barcode = curnode->barcode;
            return true;
        }
    }
    else{
        // Start a new lookup.
        unsigned long ul = k.to_ulong();
        if (table[ul] == NULL){
            curnode = NULL;
            return false;
        }
        else{
            curnode = table[ul];
            barcode = curnode->barcode;
            return true;
        }
    }
}

void bc_whitelist::parse_whitelist(string& name){ 
    fprintf(stderr, "Loading allowed barcode list...\n");
    int bc_idx = 0;
    bc cur_bc;
    kmer cur_kmer;
    gzreader reader(name);
    while(reader.next()){
        if (!str2bc(reader.line, cur_bc)){
            fprintf(stderr, "ERROR: invalid barcode %s in allowed barcode list\n", reader.line);
            exit(1);
        }
        unsigned long bc_ul = cur_bc.to_ulong();
        wl.insert(bc_ul);
        int start = 0;
        while(str2kmers(reader.line, start, cur_kmer)){
            //kmer2bc.insert(make_pair(cur_kmer.to_ulong(), bc_ul));
            kmer2bc.insert(cur_kmer, bc_ul);
        }
        bc_idx++;
    }
    fprintf(stderr, "done\n");
}

/**
 * Second should be secondary whitelist (i.e. ATAC for multiome)
 * Lookups in second whitelist will return value from first whitelist
 */
void bc_whitelist::parse_whitelist_pair(string& name1, string& name2){
    fprintf(stderr, "Loading allowed barcode lists...\n");
    vector<unsigned long> firstwl;
    bc cur_bc;
    kmer cur_kmer;
    gzreader reader(name1);
    while (reader.next()){
        if (!str2bc(reader.line, cur_bc)){
            fprintf(stderr, "ERROR: invalid barcode %s in allowed barcode list\n", reader.line);
            exit(1);
        }
        unsigned long bc_ul = cur_bc.to_ulong();
        wl.insert(bc_ul);
        // We'll need to look this up by index later.
        firstwl.push_back(bc_ul);
        int start = 0;
        while(str2kmers(reader.line, start, cur_kmer)){
            kmer2bc.insert(cur_kmer, bc_ul);
        }
    } 
    gzreader reader2(name2);
    int bc_idx = 0;
    while(reader2.next()){
        if (!str2bc(reader2.line, cur_bc)){
            fprintf(stderr, "ERROR: invalid barcode %s in whitelist\n", reader2.line);
            exit(1);
        }
        unsigned long bc_ul = cur_bc.to_ulong();
        unsigned long bc_ul_first = firstwl[bc_idx];
        wl2.emplace(bc_ul, bc_ul_first);
        int start = 0;
        while(str2kmers(reader2.line, start, cur_kmer)){
            kmer2bc2.insert(cur_kmer, bc_ul);
        }
        ++bc_idx;
    }
    fprintf(stderr, "done\n");
}

void bc_whitelist::check_lengths(){
    if (KX2/2 <= 1){
        fprintf(stderr, "ERROR: compiled with insufficiently large k (%d)\n", KX2/2);
        fprintf(stderr, "Please re-compile and set KX2 between 4 and < %d\n", (BC_LENX2/2+1)/2);
        exit(1);
    }
    // Make sure we compiled with a reasonable k-mer to barcode length ratio
    if (KX2/2 >= (BC_LENX2/2+1)/2){
        fprintf(stderr, "ERROR: you compiled for %d base barcodes and fuzzy matching k-mer k = %d\n",
            BC_LENX2/2, KX2/2);
        fprintf(stderr, "Maximum k-mer length for this barcode length is < %d\n", (BC_LENX2/2+1)/2);
        fprintf(stderr, "Please re-compile and set KX2 to a number between 4 and < %d\n", 
            (BC_LENX2/2+1)/2);
        exit(1);
    }
    if (BC_LENX2 > sizeof(unsigned long) * CHAR_BIT){
        fprintf(stderr, "ERROR: barcode length of %d will overflow numeric representation\n", BC_LENX2);
        fprintf(stderr, "Please recompile with BC_LENX2 set to a smaller value\n");
        exit(1);
    }
}

/**
 * Constructor for single whitelist (i.e. cellranger RNA)
 */
bc_whitelist::bc_whitelist(string name){
    check_lengths(); 
    
    kmer k;
    k.reset();
    for (int i = 0; i < KX2; ++i){
        k.set(i);
    }
    n_kmer_buckets = k.to_ulong();

    // Uncomment if using std::unordered_set
    //wl.max_load_factor(0.8);
    this->parse_whitelist(name);
    this->multiome = false;
}

/**
 * Constructor for two whitelists, where lines correspond to each
 * other in the two files and we want to use the barcode from the
 * first one in the BAM, etc.
 * (i.e. cellranger multiome: wl1 = RNA, wl2 = ATAC)
 */
bc_whitelist::bc_whitelist(string name1, string name2){
    check_lengths();

    kmer k;
    k.reset();
    for (int i = 0; i < KX2; ++i){
        k.set(i);
    }
    n_kmer_buckets = k.to_ulong();

    // Uncomment if using std::unordered_set
    //wl.max_load_factor(0.8);
    // Uncomment if using std::unordered_map
    //wl2.max_load_factor(0.8);
    this->parse_whitelist_pair(name1, name2);
    this->multiome = true;
}

/**
 * If a barcode contains N, try mutating it to every possible base.
 * If it contains more than one N, give up.
 */
bool bc_whitelist::mutate(const char* str, bool rc, vector<unsigned long>& alts){
    char strcpy[strlen(str)+1];
    strcpy[strlen(str)] = '\0';

    int npos = -1;
    for (int i = 0; i < strlen(str); ++i){
        strcpy[i] = str[i];
        if (str[i] == 'N' || str[i] == 'n'){
            if (npos != -1){
                // Multiple N
                return false;
            }
            npos = i;
        }
    }

    strcpy[npos] = 'A';
    if (rc){
        str2bc_rc(strcpy, cur_bc);
    }
    else{
        str2bc(strcpy, cur_bc);
    }
    alts.push_back(cur_bc.to_ulong());
    strcpy[npos] = 'C';
    if (rc){
        str2bc_rc(strcpy, cur_bc);
    }
    else{
        str2bc(strcpy, cur_bc);
    }
    alts.push_back(cur_bc.to_ulong());
    strcpy[npos] = 'G';
    if (rc){
        str2bc_rc(strcpy, cur_bc);
    }
    else{
        str2bc(strcpy, cur_bc);
    }
    alts.push_back(cur_bc.to_ulong());
    strcpy[npos] = 'T';
    if (rc){
        str2bc_rc(strcpy, cur_bc);
    }
    else{
        str2bc(strcpy, cur_bc);
    }
    alts.push_back(cur_bc.to_ulong());
    return true;
}

/**
 * Helper function for fuzzy barcode matching. Given a single
 * barcode that matches a k-mer starting at position i, stores
 * counts. Also tracks whether we haven't seen enough of any
 * particular match to match anything (so we can give up).
 *
 * Returns whether to bail out
 */
bool bc_whitelist::fuzzy_count_barcode(unsigned long barcode,
     robin_hood::unordered_node_map<unsigned long, matchinfo>& matches,
     int i,
     int& maxmatches){
    
    // Don't bother with it if it can't reach threshold.
    if (i > KX2/2 && (matches.count(barcode) == 0 ||
        matches[barcode].count < i-KX2/2)){
        // skip
    }
    else{
        if (matches.count(barcode) == 0){
            matches.emplace(barcode, matchinfo(i));
        }
        matches[barcode].count++;
        if (matches[barcode].last < i-1){
            matches[barcode].gap++;
        }
        matches[barcode].last = i;

        if (matches[barcode].count > maxmatches){
            maxmatches = matches[barcode].count;
            if (maxmatches < i-KX2/2){
                return true;
            }
        }
    }
    return false;
}

/**
 * If a string didn't match, try to match up k-mers in it to known
 * barcodes. 
 *
 * There can be a match to a known barcode with exactly one 
 * mismatch if the number of mismatching k-mers is <= k.
 *
 * There are l-k+1 possible k-mers, if l is the barcode length.
 *
 * This means that the number of matching k-mers should be 
 * >= BC_LEN-K+1 - K 
 * >= BC_LEN-(2K)+1
 *
 * Lower matching rates mean errors occurred near the beginning
 * of the sequence -- if error index is i (0-based), then an error near
 * the beginning (i < k - 1) will wipe out i + 1 k-mers
 *
 * An error in the middle (k - 1 <= i <= l-k-1) will wipe out
 * k k-mers
 *
 * An error near the end (i > l-k) will wipe out i - (l-k) + 1 k-mers
 *
 * If there is only one barcode with the number of matching k-mers 
 * greater than or equal to BC_LEN-K+1, then we accept it.
 *
 * To speed things up, we keep track of the maximum number of matches
 * encountered so far. We know that we can have up to k mismatches
 * for a sequence to pass. And if we're checking k-mer index i,
 * then the number of mismatches a sequences has is i - matches.
 * If we track the maximum matches for any barcode and that number
 * falls below i-k, we can give up altogether.
 */
unsigned long bc_whitelist::fuzzy_match(const char* str, bool rc, bool is_wl2, bool& success){
    
    robin_hood::unordered_node_map<unsigned long, matchinfo> matches;
    int maxmatches = 0;
    unsigned long mm_bc;

    unsigned long cur_bc;

    int i = 0;
    if (rc){
        int start = -1;
        while(str2kmers_rc(str, start, cur_kmer)){
            if (is_wl2){
                while(kmer2bc2.lookup(cur_kmer, cur_bc)){
                    if (fuzzy_count_barcode(cur_bc, matches, i, maxmatches)){
                        success = false;
                        return 0;
                    }
                }
            }
            else{
                while(kmer2bc.lookup(cur_kmer, cur_bc)){
                    if (fuzzy_count_barcode(cur_bc, matches, i, maxmatches)){
                        success = false;
                        return 0;
                    }
                }
            }
            ++i;
        }
    }
    else{
        int start = 0;
        while (str2kmers(str, start, cur_kmer)){
            if (is_wl2){
                while(kmer2bc2.lookup(cur_kmer, cur_bc)){
                    if (fuzzy_count_barcode(cur_bc, matches, i, maxmatches)){
                        success = false;
                        return 0;
                    }
                }
            }
            else{
                while(kmer2bc.lookup(cur_kmer, cur_bc)){
                    if (fuzzy_count_barcode(cur_bc, matches, i, maxmatches)){
                        success = false;
                        return 0;
                    }
                }
            }
            ++i;
        }
    }
    
    // If the best match has the minimum level of matching, then we can be confident it's
    // edit dist 1. If there's only one, take it.
    //
    // If the best match has more than the minimum level of matching, there could be 
    unsigned long match_bc;
    int nmatch = 0;
    
    bc bc_this;
    if (rc){
        str2bc_rc(str, bc_this);
    }
    else{
        str2bc(str, bc_this);
    }

    for (robin_hood::unordered_node_map<unsigned long, matchinfo>::iterator m = 
        matches.begin(); m != matches.end(); ++m){
        if (m->second.count >= BC_LENX2/2-KX2+1 && !(m->second.first > 0 && m->second.last < 
            BC_LENX2/2-KX2/2) && !(m->second.gap > 0 && (m->second.first > 0 || m->second.last < 
            BC_LENX2/2-KX2/2))){
           
            // Do string comparison.
            bc other(m->first);
            
            int edit_dist = 0;
            
            // Use match info to determine which end to start looking from
            //if ((!rc && m->second.first > 0) || (rc && m->second.last < BC_LENX2/2-KX2/2+1)){
             if (m->second.first > 0){
                for (int i = 0; i < BC_LENX2/2; ++i){
                    if (other.test(i*2) == bc_this.test(i*2) && 
                        other.test(i*2+1) == bc_this.test(i*2+1)){
                        // pass
                    }        
                    else{
                        ++edit_dist;
                        if (edit_dist > 1){
                            break;
                        }
                    }
                }
            }
            else{
                for (int i = BC_LENX2/2-1; i >= 0; i--){
                    if (other.test(i*2) == bc_this.test(i*2) && 
                        other.test(i*2+1) == bc_this.test(i*2+1)){
                        // pass
                    }
                    else{
                        ++edit_dist;
                        if (edit_dist > 1){
                            break;
                        }
                    }
                }
            }
            if (edit_dist == 1){
                ++nmatch;
                if (nmatch > 1){
                    // Can't say anything here
                    success = false;
                    return 0.0;
                }
                match_bc = m->first;
            }
        }
    }
    if (nmatch == 1){
        success = true;
        if (is_wl2){
            return wl2[match_bc];
        }
        return match_bc;
    }
    else{
        // Didn't manage to find anything.
        success = false;
        return 0.0;
    }
}

/**
 * Look up a string for matches in a single barcode whitelist.
 */
unsigned long bc_whitelist::lookup_aux(const char* str, bool rc, bool is_wl2, bool& success){
    
    bool try_mut = false;
    bool try_kmers = false;

    if (rc){
        if (str2bc_rc(str, cur_bc)){ 
            unsigned long ul = cur_bc.to_ulong();
            if (is_wl2){
                if (wl2.count(ul) > 0){
                    success = true;
                    return wl2[ul];
                }
                else{
                    try_kmers = true;
                }
            }
            else{
                if (wl.find(ul) != wl.end()){
                    success = true;
                    return ul;
                }
                else{
                    // Try k-mers
                    try_kmers = true;            
                }
            }
        }
        else{
            // Contains N
            try_mut = true;
        }
    }
    else{
        if (str2bc(str, cur_bc)){
            unsigned long ul = cur_bc.to_ulong();
            if (is_wl2){
                if (wl2.count(ul) > 0){
                    success = true;
                    return wl2[ul];
                }
                else{
                    try_kmers = true;
                }
            }
            else{
                if (wl.find(ul) != wl.end()){
                    success = true;
                    return ul;
                }
                else{
                    // Try k-mers
                    try_kmers = true;
                }
            }
        }
        else{
            // Contains N
            try_mut = true;
        }
    }
    if (try_mut){
        vector<unsigned long> alts;
        if (mutate(str, rc, alts)){
            // Accept if exactly one passes.
            int npass = 0;
            int pass_idx = -1;
            for (int i = 0; i < alts.size(); ++i){
                if (is_wl2){
                    if (wl2.count(alts[i]) > 0){
                        pass_idx = i;
                        if (npass > 0){
                            // Give up
                            success = false;
                            return 0.0;
                        }
                        ++npass;
                    }
                }
                else{
                    if (wl.find(alts[i]) != wl.end()){
                        pass_idx = i;
                        if (npass > 0){
                            // Give up
                            success = false;
                            return 0.0;
                        }
                        ++npass;
                    }
                }
            }
            if (pass_idx != -1){
                if (is_wl2){
                    success = true;
                    return wl2[alts[pass_idx]];
                }
                else{
                    success = true;
                    return alts[pass_idx];
                }
            }
        }
    }
    else if (try_kmers){
        return fuzzy_match(str, rc, is_wl2, success);
    }
    return false;
}

/** 
 * Look up in first (or only) whitelist, allowing revcomp or not
 */
bool bc_whitelist::lookup(const char* str, bool rc, unsigned long& ul){
    bool success;
    ul = lookup_aux(str, rc, false, success);
    return success;
}

/**
 * Look up in first (or only) whitelist, assuming forward orientation
 */
bool bc_whitelist::lookup(const char* str, unsigned long& ul){
    return lookup(str, false, ul);
}

/**
 * Look up in second whitelist, allowing revcomp or not
 */
bool bc_whitelist::lookup2(const char* str, bool rc, unsigned long& ul){
    if (!multiome){
        // Can't look up in second whitelist
        return false;
    }
    bool success;
    ul = lookup_aux(str, rc, true, success);
    return success;
}

/**
 * Look up in second whitelist, assuming forward orientation
 */
bool bc_whitelist::lookup2(const char* str, unsigned long& ul){
    return lookup2(str, false, ul);
}