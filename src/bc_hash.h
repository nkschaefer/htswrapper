#ifndef BCHASH_H
#define BCHASH_H
#include <utility>
#include <cstdlib>
#include <set>
#include <bitset>
#include <string>

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

// Designed to hold 16-base cell barcodes
typedef std::bitset<32> bc;

bool str2bc(const char* string, bc& this_bc, int len);
bool str2bc_rc(const char* string, bc& this_bc, int len);
std::string bc2str(bc& this_bc, int len);
std::string bc2str_rc(bc& this_bc, int len);
void parse_barcode_file(std::string& filename, std::set<unsigned long>& cell_barcodes);

#endif
