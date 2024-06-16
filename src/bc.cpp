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
bool str2bc(const char* str, bc& this_bc, int len){
    this_bc.reset();
    int bcbit = 0;
    for (int i = 0; i < len; ++i){
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

bool str2kmers(const char* str, int& start, kmer& cur_kmer, int k, int bc_len){
    if (start > bc_len - k){
        return false;
    }
    cur_kmer.reset();
    int kmerbit = 0;
    for (int i = start; i < start + k; ++i){
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
bool str2bc_rc(const char* str, bc& this_bc, int len){
    this_bc.reset();
    int bcbit = 0;
    for (int i = len-1; i >= 0; --i){
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
bool str2kmers_rc(const char* str, int& start, kmer& cur_kmer, int k, int bc_len){
    cur_kmer.reset();
    int kmerbit = 0;
    if (start == -1){
        // initialize start position
        start = bc_len-1;
    }
    else if (start < k-1){
        return false;
    }
    for (int i = start; i >= start-k+1; --i){
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
string bc2str(const bc& this_bc, int len){
    char strbuf[len+1];
    strbuf[len] = '\0';
    for (int i = 0; i < len; ++i){
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

string bc2str(const unsigned long ul, int len){
    bc as_bc(ul);
    return bc2str(as_bc, len);
}

/**
 * Convert a bitset representation of a DNA sequence
 * to a string representation
 */
string kmer2str(const kmer& this_bc){
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
string bc2str_rc(const bc& this_bc){
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

string bc2str_rc(const unsigned long ul){
    bc as_bits(ul);
    return bc2str_rc(as_bits);
}

/**
 * Convert a bitset representation of a DNA sequence
 * to a string representation, in reverse complement
 * orientation.
 */
string kmer2str_rc(const kmer& this_bc){
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
void parse_barcode_file(const string& filename, set<unsigned long>& cell_barcodes){
   
    gzreader reader(filename);
    while (reader.next()){
        cell_barcodes.insert(bc_ul(reader.line));
    } 
    
    fprintf(stderr, "Read %ld barcodes from file\n", cell_barcodes.size());
}

/**
 * Same as above, but populates a vector instead of a set.
 */
void parse_barcode_file(const string& filename, vector<unsigned long>& cell_barcodes){
   
    gzreader reader(filename);
    while (reader.next()){
        cell_barcodes.push_back(bc_ul(reader.line));
    } 
}

void kmer_lookup::init(int k){
    // Determine max # k-mers
    // NOTE: from cppreference:
    // https://en.cppreference.com/w/cpp/utility/bitset
    // the sequence is thought of as having its lowest indexed elements at the 
    // right, as in the binary representation of integers.
    
    // This means bitset[0] is the smallest numeric bit (2^0)
    // This means that the highest index in the bitset is the largest number
    // This means that the maximum unsigned long representation equals all bits
    // up to k*2 turned on -- higher bits are all set to 0
    
    kmer test;
    for (int i = 0; i < k*2; ++i){
        test.set(i);
    }
    
    // ulong representation of test is max possible value
    // need one more bucket than this
    unsigned long nk = test.to_ulong()+1;

    n_kmers = nk;
    table = new kmer_lookup_node*[n_kmers];
    table_last = new kmer_lookup_node*[n_kmers];
    //table = (kmer_lookup_node**)malloc(nk*sizeof(kmer_lookup_node*));
    //table_last = (kmer_lookup_node**)malloc(nk*sizeof(kmer_lookup_node*));
    for (int i = 0; i < nk; ++i){
        table[i] = NULL;
        table_last[i] = NULL;
    }
    curnode = NULL;
    initialized = true;
}

kmer_lookup::kmer_lookup(){
    initialized = false;
}

kmer_lookup::kmer_lookup(int k){
    init(k);
}

kmer_lookup::kmer_lookup(const kmer_lookup& other){
    if (other.initialized){
        n_kmers = other.n_kmers;
        table = new kmer_lookup_node*[n_kmers];
        table_last = new kmer_lookup_node*[n_kmers];
        //table = (kmer_lookup_node**)malloc(n_kmers*sizeof(kmer_lookup_node*));
        //table_last = (kmer_lookup_node**)malloc(n_kmers*sizeof(kmer_lookup_node*));
        for (int i = 0; i < n_kmers; ++i){
            if (other.table[i] != NULL){
                kmer_lookup_node* n = other.table[i];
                while (n != NULL){
                    this->insert(i, n->barcode);
                    n = n->next;
                }
            }
            else{
                table[i] = NULL;
                table_last[i] = NULL;
            }
        }
        curnode = NULL;
    }
    else{
        initialized = false;
    }
}

void kmer_lookup::free_members(kmer_lookup_node* n){
    if (!initialized){
        return;
    }
    if (n->next != NULL){
        free_members(n->next);
        n->next = NULL;
    }
    delete n;
}

kmer_lookup::~kmer_lookup(){
    if (initialized){
        // Free all nodes
        for (size_t i = 0; i < n_kmers; ++i){
            if (table[i] != NULL){
                table_last[i] = NULL;
                free_members(table[i]);
                table[i] = NULL;
            }
        }
        // Free tables
        delete[] this->table;
        delete[] this->table_last;
        //free(table);
        //free(table_last);
    }
}

void kmer_lookup::insert(unsigned long idx, unsigned long barcode){
    if (!initialized){
        fprintf(stderr, "ERROR: kmer_lookup not initialized\n");
        exit(1);
    }
    kmer_lookup_node* n = new kmer_lookup_node(barcode);
    if (table_last[idx] == NULL){
        table[idx] = n;
        table_last[idx] = n;
    }
    else{
        kmer_lookup_node* bucket = table_last[idx];
        bucket->next = n;
        table_last[idx] = n;
    }
}

void kmer_lookup::insert(kmer& k, unsigned long barcode){
    insert(k.to_ulong(), barcode);
}

bool kmer_lookup::lookup(kmer& k, unsigned long& barcode){
    return this->lookup(k.to_ulong(), barcode);
}

bool kmer_lookup::lookup(unsigned long ul, unsigned long& barcode){
    if (!initialized){
        fprintf(stderr, "ERROR: kmer_lookup not initialized\n");
        exit(1);
    }
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

void bc_whitelist::parse_whitelist_line(const char* line){
    if (!str2bc(line, cur_bc, bc_len)){
        fprintf(stderr, "ERROR: invalid barcode %s in allowed barcode list\n", line);
        exit(1);
    }
    unsigned long bc_ul = cur_bc.to_ulong();
    wl.insert(bc_ul);
    int start = 0;
    while(str2kmers(line, start, cur_kmer, k, bc_len)){
        kmer2bc.insert(cur_kmer, bc_ul);
    }
}

void bc_whitelist::parse_whitelist_line(unsigned long bc_ul){
    string bc_str = bc2str(bc_ul, bc_len);
    wl.insert(bc_ul);
    int start = 0;
    while(str2kmers(bc_str.c_str(), start, cur_kmer, k, bc_len)){
        kmer2bc.insert(cur_kmer, bc_ul);
    }
}

void bc_whitelist::parse_whitelist(string& name){ 
    fprintf(stderr, "Loading allowed barcode list...\n");
    bc cur_bc;
    kmer cur_kmer;
    gzreader reader(name);
    while(reader.next()){
        parse_whitelist_line(reader.line);
    }
    fprintf(stderr, "done\n");
}

void bc_whitelist::parse_whitelist(vector<string>& bcs){
    for (int i = 0; i < bcs.size(); ++i){
        parse_whitelist_line(bcs[i].c_str());
    }
}

void bc_whitelist::parse_whitelist(vector<unsigned long>& bcs){
    for (int i = 0; i < bcs.size(); ++i){
        parse_whitelist_line(bcs[i]);
    }
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
        if (!str2bc(reader.line, cur_bc, bc_len)){
            fprintf(stderr, "ERROR: invalid barcode %s in allowed barcode list\n", reader.line);
            exit(1);
        }
        unsigned long bc_ul = cur_bc.to_ulong();
        wl.insert(bc_ul);
        // We'll need to look this up by index later.
        firstwl.push_back(bc_ul);
        int start = 0;
        while(str2kmers(reader.line, start, cur_kmer, k, bc_len)){
            kmer2bc.insert(cur_kmer, bc_ul);
        }
    } 
    gzreader reader2(name2);
    int bc_idx = 0;
    while(reader2.next()){
        if (!str2bc(reader2.line, cur_bc, bc_len)){
            fprintf(stderr, "ERROR: invalid barcode %s in whitelist\n", reader2.line);
            exit(1);
        }
        unsigned long bc_ul = cur_bc.to_ulong();
        unsigned long bc_ul_first = firstwl[bc_idx];
        wl2.emplace(bc_ul, bc_ul_first);
        wl2_rev.emplace(bc_ul_first, bc_ul);
        int start = 0;
        while(str2kmers(reader2.line, start, cur_kmer, k, bc_len)){
            kmer2bc2.insert(cur_kmer, bc_ul);
        }
        ++bc_idx;
    }
    fprintf(stderr, "done\n");
}

void bc_whitelist::parse_whitelist_pair(vector<string>& list1, vector<string>& list2){
    if (list1.size() != list2.size()){
        fprintf(stderr, "ERROR: different whitelist lengths: %ld vs %ld\n", list1.size(), list2.size());
        return;
    }
    kmer cur_kmer;
    bc cur_bc1;
    bc cur_bc2;
    for (int i = 0; i < list1.size(); ++i){
        if (!str2bc(list1[i].c_str(), cur_bc1, bc_len)){
            fprintf(stderr, "ERROR: invalid barcode %s in whitelist1\n", list1[i].c_str());
            return;
        }
        if (!str2bc(list2[i].c_str(), cur_bc2, bc_len)){
            fprintf(stderr, "ERROR: invalid barcode %s in whitelist2\n", list2[i].c_str());
            return;
        }
        unsigned long ul1 = cur_bc1.to_ulong();
        unsigned long ul2 = cur_bc2.to_ulong();
        wl.insert(ul1);
        wl2.emplace(ul2, ul1);
        wl2_rev.emplace(ul1, ul2);
        int start = 0;
        while(str2kmers(list1[i].c_str(), start, cur_kmer, k, bc_len)){
            kmer2bc.insert(cur_kmer, ul1);
        }
        start = 0;
        while(str2kmers(list2[i].c_str(), start, cur_kmer, k, bc_len)){
            kmer2bc2.insert(cur_kmer, ul2);
        }
    }
}

void bc_whitelist::parse_whitelist_pair(vector<unsigned long>& list1, vector<unsigned long>& list2){
    if (list1.size() != list2.size()){
        fprintf(stderr, "ERROR: different whitelist lengths: %ld vs %ld\n", list1.size(), list2.size());
        return;
    }
    kmer cur_kmer;
    for (int i = 0; i < list1.size(); ++i){
        wl.insert(list1[i]);
        wl2.emplace(list2[i], list1[i]);
        wl2_rev.emplace(list1[i], list2[i]);
        int start = 0;
        string bc_str = bc2str(list1[i], bc_len);
        while(str2kmers(bc_str.c_str(), start, cur_kmer, k, bc_len)){
            kmer2bc.insert(cur_kmer, list1[i]);
        }
        bc_str = bc2str(list2[i], bc_len);
        start = 0;
        while(str2kmers(bc_str.c_str(), start, cur_kmer, k, bc_len)){
            kmer2bc2.insert(cur_kmer, list2[i]);
        }
    }
}

unsigned long bc_whitelist::wl2towl1(unsigned long barcode){
    if (!multiome){
        return 0;
    }
    if (this->wl2.count(barcode) > 0){
        return this->wl2[barcode];
    }
    return 0;
}

unsigned long bc_whitelist::wl1towl2(unsigned long barcode){
    if (!multiome){
        return 0;
    }
    if (this->wl2_rev.count(barcode) > 0){
        return this->wl2_rev[barcode];
    }
    return 0;
}

void bc_whitelist::check_lengths(){
    if (KX2/2 <= 1){
        fprintf(stderr, "ERROR: compiled with insufficiently large k (%d)\n", KX2/2);
        fprintf(stderr, "Please re-compile and set KX2 between 4 and %d\n", (BC_LENX2/2+1)/2);
        exit(1);
    }
    // Make sure we compiled with a reasonable k-mer to barcode length ratio
    bool k_size_ok = true;
    bool can_equal = false;
    if ((BC_LENX2/2 % 2 == 0)){
        // even barcode length.
        // k must be less than (bc length + 1)/2, which is a fraction.
        // k can therefore equal the integer division result.
        if (KX2/2 > (BC_LENX2/2+1)/2){
            k_size_ok = false;
            can_equal = true;
        }
    }
    else{
        // odd barcode length.
        // k must be less than (bc length + 1)/2
        if (KX2/2 >= (BC_LENX2/2+1)/2){
            k_size_ok = false;
            can_equal = false;
        }
    }
    if (!k_size_ok){
        fprintf(stderr, "ERROR: you compiled for %d base barcodes and fuzzy matching k-mer k = %d\n",
            BC_LENX2/2, KX2/2);
        string eqstr;
        if (can_equal){
            eqstr = "<=";
        }
        else{
            eqstr = "<";
        }
        fprintf(stderr, "Maximum k-mer length for this barcode length is %s %d\n", eqstr.c_str(),
            (BC_LENX2/2+1)/2);
        fprintf(stderr, "Please re-compile with KX2= a suitable number\n");
        exit(1);
    }
    if (BC_LENX2 > sizeof(unsigned long) * CHAR_BIT){
        fprintf(stderr, "ERROR: barcode length of %d will overflow numeric representation\n", BC_LENX2);
        fprintf(stderr, "Please recompile with BC_LENX2 set to a smaller value\n");
        exit(1);
    }
}

void bc_whitelist::init_aux(int bc_len, int k, bool has_two){
    // Ensure everything was compiled correctly
    check_lengths();

    // Store lengths
    this->bc_len = bc_len;
    this->k = k;
    
    this->kmer2bc.init(k);
    if (has_two){
        this->kmer2bc2.init(k);
    }

    // Default to allow fuzzy matching
    this->exact_only = false;

    // Uncomment if using std::unordered_set
    //wl.max_load_factor(0.8);
    // Uncomment if using std::unordered_map
    //wl2.max_load_factor(0.8);
}

void bc_whitelist::init(string name, int bc_len, int k){
    init_aux(bc_len, k, false);
    this->parse_whitelist(name);
    this->multiome = false;
    initialized = true;
}

void bc_whitelist::init(string name1, string name2, int bc_len, int k){
    init_aux(bc_len, k, true);
    this->parse_whitelist_pair(name1, name2);
    this->multiome = true;
    initialized = true;
}

void bc_whitelist::init(vector<string>& bcs, int bc_len, int k){
    init_aux(bc_len, k, false);
    this->parse_whitelist(bcs);
    this->multiome = false;
    initialized = true;
}

void bc_whitelist::init(vector<string>& bcs1, vector<string>& bcs2, int bc_len, int k){
    init_aux(bc_len, k, true);
    this->parse_whitelist_pair(bcs1, bcs2);
    this->multiome = true;
    initialized = true;
}

void bc_whitelist::init(vector<unsigned long>& bcs, int bc_len, int k){
    init_aux(bc_len, k, false);
    this->parse_whitelist(bcs);
    this->multiome = false;
    this->initialized = true;
}

void bc_whitelist::init(vector<unsigned long>& bcs1, vector<unsigned long>& bcs2, int bc_len, int k){
    init_aux(bc_len, k, true);
    this->parse_whitelist_pair(bcs1, bcs2);
    this->multiome = true;
    this->initialized = true;
}

/**
 * Constructor for single whitelist (i.e. cellranger RNA)
 */
bc_whitelist::bc_whitelist(string name, int bc_len, int k){
   init(name, bc_len, k); 
}

/**
 * Constructor for two whitelists, where lines correspond to each
 * other in the two files and we want to use the barcode from the
 * first one in the BAM, etc.
 * (i.e. cellranger multiome: wl1 = RNA, wl2 = ATAC)
 */
bc_whitelist::bc_whitelist(string name1, string name2, int bc_len, int k){
    init(name1, name2, bc_len, k);
}

bc_whitelist::bc_whitelist(vector<string>& bcs, int bc_len, int k){
    init(bcs, bc_len, k);
}

bc_whitelist::bc_whitelist(vector<string>& bcs1, vector<string>& bcs2, int bc_len, int k){
    init(bcs1, bcs2, bc_len, k);
}

bc_whitelist::bc_whitelist(vector<unsigned long>& bcs, int bc_len, int k){
    init(bcs, bc_len, k);
}

bc_whitelist::bc_whitelist(vector<unsigned long>& bcs1, vector<unsigned long>& bcs2,
    int bc_len, int k){
    init(bcs1, bcs2, bc_len, k);
}

bc_whitelist::bc_whitelist(){
    initialized = false;
}

void bc_whitelist::exact_matches_only(){
    exact_only = true;
}

void bc_whitelist::allow_mismatch(){
    exact_only = false;
}

/**
 * If a barcode contains N, try mutating it to every possible base.
 * If it contains more than one N, give up.
 */
bool bc_whitelist::mutate(const char* str, bool rc, vector<unsigned long>& alts){
    if (this->exact_only){
        return true;
    }

    char strcpy[bc_len+1];
    strcpy[bc_len] = '\0';

    int npos = -1;

    if (rc){
        int i2 = 0;
        for (int i = strlen(str)-bc_len; i < strlen(str); ++i){
            strcpy[i2] = str[i];
            if (str[i] == 'N' || str[i] == 'n'){
                if (npos != -1){
                    // Multiple N
                    return false;
                }
                npos = i2;
            }
            ++i2;
        }
    }
    else{
        for (int i = 0; i < bc_len; ++i){
            strcpy[i] = str[i];
            if (str[i] == 'N' || str[i] == 'n'){
                if (npos != -1){
                    // Multiple N
                    return false;
                }
                npos = i;
            }
        }
    }

    strcpy[npos] = 'A';
    if (rc){
        str2bc_rc(strcpy, cur_bc, bc_len);
    }
    else{
        str2bc(strcpy, cur_bc, bc_len);
    }
    alts.push_back(cur_bc.to_ulong());
    strcpy[npos] = 'C';
    if (rc){
        str2bc_rc(strcpy, cur_bc, bc_len);
    }
    else{
        str2bc(strcpy, cur_bc, bc_len);
    }
    alts.push_back(cur_bc.to_ulong());
    strcpy[npos] = 'G';
    if (rc){
        str2bc_rc(strcpy, cur_bc, bc_len);
    }
    else{
        str2bc(strcpy, cur_bc, bc_len);
    }
    alts.push_back(cur_bc.to_ulong());
    strcpy[npos] = 'T';
    if (rc){
        str2bc_rc(strcpy, cur_bc, bc_len);
    }
    else{
        str2bc(strcpy, cur_bc, bc_len);
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
    if (i > k && (matches.count(barcode) == 0 ||
        matches[barcode].count < i-k)){
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
            if (maxmatches < i-k){
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
    if (this->exact_only){
        success = false;
        return 0;
    }
    robin_hood::unordered_node_map<unsigned long, matchinfo> matches;
    int maxmatches = 0;
    unsigned long mm_bc;

    unsigned long cur_bc;
    kmer cur_kmer;

    int i = 0;
    if (rc){
        int start = -1;
        while(str2kmers_rc(str, start, cur_kmer, k, bc_len)){
            if (is_wl2){
                lock_guard<mutex> lock(kmer2mutex);
                while(kmer2bc2.lookup(cur_kmer, cur_bc)){
                    if (fuzzy_count_barcode(cur_bc, matches, i, maxmatches)){
                        success = false;
                        return 0;
                    }
                }
            }
            else{
                lock_guard<mutex> lock(kmer1mutex);
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
        while (str2kmers(str, start, cur_kmer, k, bc_len)){
            if (is_wl2){
                lock_guard<mutex> lock(kmer2mutex);
                while(kmer2bc2.lookup(cur_kmer, cur_bc)){
                    if (fuzzy_count_barcode(cur_bc, matches, i, maxmatches)){
                        success = false;
                        return 0;
                    }
                }
            }
            else{
                lock_guard<mutex> lock(kmer1mutex);
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
        str2bc_rc(str, bc_this, bc_len);
    }
    else{
        str2bc(str, bc_this, bc_len);
    }

    for (robin_hood::unordered_node_map<unsigned long, matchinfo>::iterator m = 
        matches.begin(); m != matches.end(); ++m){
        if (m->second.count >= bc_len - 2*k + 1 && !(m->second.first > 0 && m->second.last < 
            bc_len-k) && !(m->second.gap > 0 && (m->second.first > 0 || m->second.last < 
            bc_len-k))){
           
            // Do string comparison.
            bc other(m->first);
            
            int edit_dist = 0;
            
            // Use match info to determine which end to start looking from
            //if ((!rc && m->second.first > 0) || (rc && m->second.last < BC_LENX2/2-KX2/2+1)){
             if (m->second.first > 0){
                for (int i = 0; i < bc_len; ++i){
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
                for (int i = bc_len-1; i >= 0; i--){
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
    if (!initialized){
        fprintf(stderr, "ERROR: whitelist not initialized\n");
        return false;
    }
    bool try_mut = false;
    bool try_kmers = false;

    if (rc){
        if (str2bc_rc(str, cur_bc, bc_len)){ 
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
        bc cur_bc;
        if (str2bc(str, cur_bc, bc_len)){
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

/**
 * Convenience functions to help remember how to look up different types of barcodes.
 */

// Remove a barcode from the beginning of a read.
void bc_whitelist::trim_begin(char* str){
    int nremain = strlen(str)-bc_len;
    memmove(&str[0], &str[bc_len], nremain);
    str[nremain] = '\0';
}

// Remove a barcode from the end of a read.
void bc_whitelist::trim_end(char* str){
    str[strlen(str)-bc_len] = '\0';
}

// whitelist 1, beginning of read, forward
bool bc_whitelist::lookup1_bf(char* str, unsigned long& ul, bool trim){
    if (lookup(str, false, ul)){
        if (trim){
            trim_begin(str);
        }
        return true;
    }
    return false;
}

// whitelist 1, beginning of read, reverse complement
bool bc_whitelist::lookup1_br(char* str, unsigned long& ul, bool trim){
    if (lookup(str, true, ul)){
        if (trim){
            trim_begin(str);
        }
        return true;
    }
    return false;
}

// whitelist 1, end of read, forward
bool bc_whitelist::lookup1_ef(char* str, unsigned long& ul, bool trim){
    if (lookup(str + strlen(str) - bc_len, false, ul)){
        if (trim){
            trim_end(str);
        }
        return true;
    }
    return false;
}

// whitelist 1, end of read, reverse complement
bool bc_whitelist::lookup1_er(char* str, unsigned long& ul, bool trim){
    if (lookup(str + strlen(str) - bc_len, true, ul)){
        if (trim){
            trim_end(str);
        }
        return true;
    }
    return false;
}

// whitelist 2, beginning of read, forward
bool bc_whitelist::lookup2_bf(char* str, unsigned long& ul, bool trim){
    if (lookup2(str, false, ul)){
        if (trim){
            trim_begin(str);
        }
        return true;
    }
    return false;
}

// whitelist 2, beginning of read, reverse complement
bool bc_whitelist::lookup2_br(char* str, unsigned long& ul, bool trim){
    if (lookup2(str, true, ul)){
        if (trim){
            trim_begin(str);
        }
        return true;
    }
    return false;
}

// whitelist 2, end of read, forward
bool bc_whitelist::lookup2_ef(char* str, unsigned long& ul, bool trim){
    if (lookup2(str + strlen(str) - bc_len, false, ul)){
        if (trim){
            trim_end(str);
        }
        return true;
    }
    return false;
}

// whitelist 2, end of read, reverse complement
bool bc_whitelist::lookup2_er(char* str, unsigned long& ul, bool trim){
    if (lookup2(str + strlen(str) - bc_len, true, ul)){
        if (trim){
            trim_end(str);
        }
        return true;
    }
    return false;
}

