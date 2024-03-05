#ifndef _HTSWRAPPER_BC_H
#define _HTSWRAPPER_BC_H
#include <utility>
#include <cstdlib>
#include <set>
#include <unordered_set>
#include <bitset>
#include <string>
#include <unordered_map>
#include <vector>
#include "robin_hood/robin_hood.h"

// Represent cell barcodes
typedef std::bitset<BC_LENX2> bc;

// Represent k-mers for fuzzy matching cell barcodes
// NOTE: KX2 must be less than BC_LENX2
// Both must be divisible by 2
typedef std::bitset<KX2> kmer;

// A hash function that could be used to speed up insertion/retrieval of hashed
// barcodes to unordered maps

struct ulong_hash_func {
    // https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
    unsigned long operator()(const unsigned long& x) const{
        unsigned long y = ((x >> 16) ^ x) * 0x45d9f3b;
        y = ((y >> 16) ^ y) * 0x45d9f3b;
        y  = (y >> 16) ^ y;
        return y;
    }
};

struct ulong_equal_func {
    bool operator()( const unsigned long& u1, const unsigned long& u2) const{
        return u1 == u2;
    }
};

bool str2bc(const char* string, bc& this_bc, int len=BC_LENX2/2);
bool str2bc_rc(const char* string, bc& this_bc, int len=BC_LENX2/2);
bool str2kmers(const char* string, int& start, kmer& cur_kmer, int k=KX2/2, int bc_len=BC_LENX2/2);
bool str2kmers_rc(const char* string, int& start, kmer& cur_kmer, int k=KX2/2, int bc_len=BC_LENX2/2);

std::string bc2str(const bc& this_bc, int len=BC_LENX2/2);
std::string bc2str(const unsigned long ul, int len=BC_LENX2/2);
std::string kmer2str(const kmer& this_kmer);
std::string bc2str_rc(const bc& this_bc);
std::string bc2str_rc(const unsigned long ul);
std::string kmer2str_rc(const kmer& this_kmer);

unsigned long bc_ul(std::string& barcode);
unsigned long bc_ul(char* barcode);
void parse_barcode_file(std::string& filename, std::set<unsigned long>& cell_barcodes);

// Class to represent barcode whitelists

// For fuzzy matching
struct matchinfo{
    int first;
    int last;
    int count;
    int gap;
    matchinfo(){ first = -1; last = -1; count = 0; gap = 0; };
    matchinfo(int f){ first = f; last = f; count = 0; gap = 0;};
    matchinfo(const matchinfo& m){ first = m.first; last = m.last; count = m.count; gap = m.gap; };
};

// For storing k-mers for fuzzy matching
// Implement as a hash table: key = unsigned long interpretation of kmer
//  this is index into array of nodes
//  nodes are linked lists

struct kmer_lookup_node{
    kmer_lookup_node* next;
    unsigned long barcode;
    kmer_lookup_node(unsigned long ul){ barcode = ul; next = NULL; };
};

class kmer_lookup{
    private:
        size_t n_kmers;
        kmer_lookup_node** table;
        kmer_lookup_node** table_last;
        void free_members(kmer_lookup_node* n);
            kmer_lookup_node* curnode;
    public:
        kmer_lookup();
        ~kmer_lookup();
        void insert(kmer& k, unsigned long barcode);
        bool lookup(kmer& k, unsigned long& barcode);
};

class bc_whitelist{
    private:
        // How long are cell barcodes? default = BC_LENX2/2        
        int bc_len;
        // How long are fuzzy matching k-mers? default = KX2/2
        int k;
        
        // Set true to skip fuzzy matching 
        bool exact_only;

        bool initialized;

        // One whitelist or two?
        bool multiome;
        
        // For normal whitelists: a list of barcodes
        robin_hood::unordered_flat_set<unsigned long> wl;
        
        // For multiome ATAC: look up the ATAC barcode and get an RNA-seq barcode
        robin_hood::unordered_flat_map<unsigned long, unsigned long> wl2;
        
        // For fuzzy matching: how many indices into kmer lookup arrays do we need?
        int n_kmer_buckets;
        
        // Look up barcode by shorter k-mer for fuzzy matching        
        kmer_lookup kmer2bc;
        kmer_lookup kmer2bc2;

        // Parse a single whitelist
        void parse_whitelist(std::string& filename);
        // Parse two whitelists (multiome)
        void parse_whitelist_pair(std::string& filename1, std::string& filename2);
        // Parse a whitelist already loaded into a string vector
        void parse_whitelist(std::vector<std::string>& bcs);        
        // Helper function for single whitelist
        void parse_whitelist_line(const char* str);

        bc cur_bc;
        kmer cur_kmer;
        
        bool mutate(const char* str, bool rc, std::vector<unsigned long>& alts);
        bool fuzzy_count_barcode(unsigned long barcode, 
            robin_hood::unordered_node_map<unsigned long, matchinfo>& matches,
            int i, int& maxmatches);
        unsigned long fuzzy_match(const char* str, bool rc, bool is_wl2, bool& success);
        
        unsigned long lookup_aux(const char* str, bool rc, bool is_wl2, bool& success);
        
        // Make sure variables were set correctly at compile time
        void check_lengths();
        
        // Helper function to trim barcode off beginning/end of string
        void trim_begin(char* str);
        void trim_end(char* str);
        
        void init_aux(int bc_len, int k);

    public:
        
        bc_whitelist();

        void init(std::string filename, int bclen=BC_LENX2/2, int k=KX2/2);
        void init(std::string filename, std::string filename2, int bclen=BC_LENX2/2, int k=KX2/2);
        void init(std::vector<std::string>& bcs, int bclen=BC_LENX2/2, int k=KX2/2);

        bc_whitelist(std::string filename, int bclen=BC_LENX2/2, int k=KX2/2);
        bc_whitelist(std::string filename, std::string filename2, int bclen=BC_LENX2/2, int k=KX2/2);
        bc_whitelist(std::vector<std::string>& bcs, int bclen=BC_LENX2/2, int k=KX2/2);
        
        // Toggle setting that dictates whether or not we allow up to one N or
        // base mismatch per barcode
        void exact_matches_only();
        void allow_mismatch();

        // Look up a barcode in the first whitelist. Return success/failure
        // and set the resulting key (on success) to the last parameter 
        bool lookup(const char* str, bool rc, unsigned long& bc_ul);
        
        // Version that assumes not reverse complement
        bool lookup(const char* str, unsigned long& bc_ul);

        // Look up a barcode in the second whitelist. Return success/failure
        // and set the resulting key (on success) to the last parameter.
        // The resulting key is the key in the first whitelist corresponding
        // to the barcode in the second whitelist. This is designed the way
        // 10X Genomics does it for multiome data (first whitelist = RNA, 
        // second whitelist = ATAC).
        bool lookup2(const char* str, bool rc, unsigned long& bc_ul);
        
        // Version that assumes not reverse complement
        bool lookup2(const char* str, unsigned long& bc_ul);

        bool lookup1_bf(char* str, unsigned long& bc_ul, bool trim);
        bool lookup1_ef(char* str, unsigned long& bc_ul, bool trim);
        bool lookup1_br(char* str, unsigned long& bc_ul, bool trim);
        bool lookup1_er(char* str, unsigned long& bc_ul, bool trim);
        bool lookup2_bf(char* str, unsigned long& bc_ul, bool trim);
        bool lookup2_ef(char* str, unsigned long& bc_ul, bool trim);
        bool lookup2_br(char* str, unsigned long& bc_ul, bool trim);
        bool lookup2_er(char* str, unsigned long& bc_ul, bool trim);

};


#endif
