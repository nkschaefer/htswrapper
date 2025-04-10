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
#include <mutex>
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
void parse_barcode_file(const std::string& filename, std::set<unsigned long>& cell_barcodes);
void parse_barcode_file(const std::string& filename, std::vector<unsigned long>& cell_barcodes);

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
    kmer_lookup_node(){ barcode = 0; next = NULL; };
};

class kmer_lookup{
    private:
        bool initialized;
        size_t n_kmers;
        kmer_lookup_node** table;
        kmer_lookup_node** table_last;
        void free_members(kmer_lookup_node* n);
            kmer_lookup_node* curnode;
    public:
        kmer_lookup();
        kmer_lookup(int k);
        kmer_lookup(const kmer_lookup& other);

        void init(int k=KX2/2);
        ~kmer_lookup();
        
        void insert(unsigned long ul, unsigned long barcode);
        void insert(kmer& k, unsigned long barcode);
        bool lookup(unsigned long ul, unsigned long& barcode);
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
        
        bc cur_bc;
        kmer cur_kmer;

        // For normal whitelists: a list of barcodes
        robin_hood::unordered_flat_set<unsigned long> wl;
        
        // For multiome ATAC: look up the ATAC barcode and get an RNA-seq barcode
        robin_hood::unordered_flat_map<unsigned long, unsigned long> wl2;
        // Flip the above around (usually not used)
        robin_hood::unordered_flat_map<unsigned long, unsigned long> wl2_rev;
        
        std::vector<std::bitset<BC_LENX2> > fuzzy_masks;
        std::vector<robin_hood::unordered_flat_map<unsigned long, std::vector<unsigned long> > > fuzzy_matches;
        std::vector<robin_hood::unordered_flat_map<unsigned long, std::vector<unsigned long> > > fuzzy_matches2;

        // Look up barcode by shorter k-mer for fuzzy matching        
        kmer_lookup kmer2bc;
        kmer_lookup kmer2bc2;
        
        std::mutex kmer1mutex;
        std::mutex kmer2mutex;

        // Parse a single whitelist
        void parse_whitelist(std::string& filename);
        // Parse two whitelists (multiome)
        void parse_whitelist_pair(std::string& filename1, std::string& filename2);
        
        // Parse a whitelist already loaded into a string vector
        void parse_whitelist(std::vector<std::string>& bcs);        
        void parse_whitelist(std::vector<unsigned long>& bcs);
        void parse_whitelist_pair(std::vector<std::string>& bcs1, std::vector<std::string>& bcs2);
        void parse_whitelist_pair(std::vector<unsigned long>& bcs1, std::vector<unsigned long>& bcs2);

        // Helper function for single whitelist
        void parse_whitelist_line(const char* str);
        void parse_whitelist_line(unsigned long ul);

        bool mutate(const char* str, int len, bool rc, std::vector<unsigned long>& alts);
        bool fuzzy_count_barcode(unsigned long barcode, 
            robin_hood::unordered_node_map<unsigned long, matchinfo>& matches,
            int i, int& maxmatches);
        unsigned long fuzzy_match(const char* str, bool rc, bool is_wl2, bool& success);
        
        unsigned long lookup_aux(const char* str, int len, bool rc, bool is_wl2, bool& success, bool& exact_match);
        
        // Make sure variables were set correctly at compile time
        void check_lengths();
        
        // Helper function to trim barcode off beginning/end of string
        void trim_begin(char* str, int len);
        void trim_end(char* str, int len);
        
        // Helper function to replace barcode sequence at beginning/end of string
        void replace_begin(char* str, bool rc, unsigned long ul);
        void replace_end(char* str, bool rc, unsigned long ul, int len);

        void init_aux(int bc_len, int k, bool has_second_wl);

    public:
        
        bc_whitelist();
        ~bc_whitelist();

        void init(std::string filename, int bclen=BC_LENX2/2, int k=KX2/2);
        void init(std::string filename, std::string filename2, int bclen=BC_LENX2/2, int k=KX2/2);
        void init(std::vector<std::string>& bcs, int bclen=BC_LENX2/2, int k=KX2/2);
        void init(std::vector<std::string>& bcs1, std::vector<std::string>& bcs2, int bclen=BC_LENX2/2, 
            int k=KX2/2);
        void init(std::vector<unsigned long>& bcs, int bclen=BC_LENX2/2, int k=KX2/2);
        void init(std::vector<unsigned long>& bcs1, std::vector<unsigned long>& bcs2, int bclen=BC_LENX2/2,
            int k=KX2/2);
        
        bc_whitelist(std::string filename, int bclen=BC_LENX2/2, int k=KX2/2);
        bc_whitelist(std::string filename, std::string filename2, int bclen=BC_LENX2/2, int k=KX2/2);
        bc_whitelist(std::vector<std::string>& bcs, int bclen=BC_LENX2/2, int k=KX2/2);
        bc_whitelist(std::vector<std::string>& bcs1, std::vector<std::string>& bcs2, int bclen=BC_LENX2/2, 
            int k=KX2/2);
        bc_whitelist(std::vector<unsigned long>& bcs1, int bclen=BC_LENX2/2, int k=KX2/2);
        bc_whitelist(std::vector<unsigned long>& bcs1, std::vector<unsigned long>& bcs2, int bclen=BC_LENX2/2,
            int k=KX2/2);

        // Toggle setting that dictates whether or not we allow up to one N or
        // base mismatch per barcode
        void exact_matches_only();
        void allow_mismatch();
        
        // Convert a (whitelisted) ATAC-seq barcode to an RNA-seq barcode
        unsigned long wl2towl1(unsigned long barcode);
        // Convert a (whitelisted) RNA-seq barcode to an ATAC-seq barcode
        unsigned long wl1towl2(unsigned long barcode);

        // Look up a barcode in the first whitelist. Return success/failure
        // and set the resulting key (on success) to the last parameter 
        bool lookup(const char* str, bool rc, unsigned long& bc_ul, bool& exact_match, int len=-1);
        
        // Version that assumes not reverse complement
        bool lookup(const char* str, unsigned long& bc_ul, bool& exact_match, int len=-1);

        // Look up a barcode in the second whitelist. Return success/failure
        // and set the resulting key (on success) to the last parameter.
        // The resulting key is the key in the first whitelist corresponding
        // to the barcode in the second whitelist. This is designed the way
        // 10X Genomics does it for multiome data (first whitelist = RNA, 
        // second whitelist = ATAC).
        bool lookup2(const char* str, bool rc, unsigned long& bc_ul, bool& exact_match, int len=-1);
        
        // Version that assumes not reverse complement
        bool lookup2(const char* str, unsigned long& bc_ul, bool& exact_match, int len=-1);

        bool lookup1_bf(char* str, int len, unsigned long& bc_ul, bool trim = false, bool replace = false);
        bool lookup1_ef(char* str, int len, unsigned long& bc_ul, bool trim = false, bool replace = false);
        bool lookup1_br(char* str, int len, unsigned long& bc_ul, bool trim = false, bool replace = false);
        bool lookup1_er(char* str, int len, unsigned long& bc_ul, bool trim = false, bool replace = false);
        bool lookup2_bf(char* str, int len, unsigned long& bc_ul, bool trim = false, bool replace = false);
        bool lookup2_ef(char* str, int len, unsigned long& bc_ul, bool trim = false, bool replace = false);
        bool lookup2_br(char* str, int len, unsigned long& bc_ul, bool trim = false, bool replace = false);
        bool lookup2_er(char* str, int len, unsigned long& bc_ul, bool trim = false, bool replace = false);

        // Does this contain a second whitelist?
        bool two_lists();
        int len_bc();
};

// Modifies a string representation of a cell barcode to include unique
// library information before printing.
void mod_bc_libname(std::string& bc, const std::string& libname, bool cellranger,
   bool seurat, bool underscore);


#endif
