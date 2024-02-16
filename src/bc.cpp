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

/**
 * Helper function for below: trims extraneous stuff off the 
 * end of a cell barcode string and returns an unsigned long
 * representation.
 *
 */
unsigned long hash_bc(string& barcode){
    // Strip anything that's not ACGT off the end of the sequence.
    int ntrim = 0;
    for (int i = barcode.length()-1; i >= 0; --i){
        if (barcode[i] != 'A' && barcode[i] != 'C' &&
            barcode[i] != 'G' && barcode[i] != 'T'){
            // trim it.
            ++ntrim;
        }
        else{
            break;
        }
    }
    if (ntrim > 0){
        barcode = barcode.substr(0, barcode.length()-ntrim);
    }
    bc as_bitset;
    str2bc(barcode.c_str(), as_bitset);
    return as_bitset.to_ulong();
}

/**
 * Same as above, but if working with a character buffer,
 * avoid copying into a std::string for speed
 */
unsigned long hash_bc(char* barcode){
    int ntrim = 0;
    for (int i = strlen(barcode)-1; i >= 0; --i){
        if (barcode[i] != 'A' && barcode[i] != 'C' &&
            barcode[i] != 'G' && barcode[i] != 'T'){
            // trim it.
            ++ntrim;
        }
        else{
            break;
        }
    }
    if (ntrim > 0){
        barcode[strlen(barcode)-ntrim] = '\0';
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
        cell_barcodes.insert(hash_bc(reader.line));
    } 
    
    fprintf(stderr, "Read %ld barcodes from file\n", cell_barcodes.size());
}

/**
 * Read data about whitelisted 10X ATAC / RNA barcodes. Required to map ATAC to RNA barcodes.
 */
void parse_whitelists(string& whitelist_atac_filename,
    string& whitelist_rna_filename, 
    bcset& atac_bc2idx,
    vector<unsigned long>& whitelist_rna,
    bcset& rna_bc2idx){

    if (whitelist_atac_filename != ""){
        ifstream wlfile_atac(whitelist_atac_filename);
        int bc_idx = 0;
        string line;
        bc atac_bc;
        
        gzreader reader(whitelist_atac_filename);
        while (reader.next()){
            if (!str2bc(reader.line, atac_bc)){
                fprintf(stderr, "ERROR: invalid barcode %s in ATAC whitelist\n", reader.line);
                exit(1);
            }
            atac_bc2idx.emplace(atac_bc.to_ulong(), bc_idx);
            bc_idx++;
        }
        
        // Now map these to RNA barcodes
        gzreader reader2(whitelist_rna_filename);
        bc_idx = 0;
        while(reader2.next()){
            if (!str2bc(reader2.line, atac_bc)){
                fprintf(stderr, "ERROR: invalid barcode %s in RNA-seq whitelist\n", reader2.line);
                exit(1);
            }
            whitelist_rna.push_back(atac_bc.to_ulong());
            rna_bc2idx.emplace(atac_bc.to_ulong(), bc_idx);
            bc_idx++;
        }
    }
    else if (whitelist_rna_filename != ""){
        int bc_idx = 0;
        string line;
        bc rna_bc;
        gzreader reader(whitelist_rna_filename);
        while(reader.next()){
            if (!str2bc(reader.line, rna_bc)){
                fprintf(stderr, "ERROR: invalid barcode %s in RNA whitelist\n", reader.line);
                exit(1);
            }
            rna_bc2idx.emplace(rna_bc.to_ulong(), bc_idx);
            whitelist_rna.push_back(rna_bc.to_ulong());
            bc_idx++;
        }
    }    
}

/**
 *  Given a candidate barcode, checks against the set of all legitimate
 *  barcodes. If found, returns true and sets result_bc_str to a string
 *  representation of the correct barcode.
 */
int match_bc(const char* cur_bc,   
    bool rc, 
    bcset& bc2idx,
    multimap<unsigned long, unsigned long>& kmer2bc){
    
    static bc bc_binary;
    static char bc_str[BC_LENX2/2 + 1]; // For replacing N characters

    bool success = false;
    if (rc){
        if (str2bc_rc(cur_bc, bc_binary)){ 
            unsigned long ul = bc_binary.to_ulong();
            if (bc2idx.count(ul) > 0){
                return bc2idx[ul];
            }
            else{
                return -1;
            }
        }
        else{
            // Contains N.
            // Only allow up to one.
            int npos = -1;
            for (int i = 0; i < BC_LENX2/2; ++i){
                if (cur_bc[i] == 'N'){
                    if (npos != -1){
                        // Multiple Ns; bail.
                        return -1;
                    }
                    else{
                        npos = i;
                    }
                }
            }
            
            // Try to mutate to every possible letter.
            bc bc_A;
            bc bc_C;
            bc bc_G;
            bc bc_T;

            bool pass_A = false;
            bool pass_C = false;
            bool pass_G = false;
            bool pass_T = false;
            
            strncpy(&bc_str[0], cur_bc, BC_LENX2/2);
            bc_str[BC_LENX2/2] = '\0';

            bc_str[npos] = 'A';
            if (str2bc_rc(bc_str, bc_A) && bc2idx.count(bc_A.to_ulong()) > 0){
                pass_A = true;
            } 
            bc_str[npos] = 'C';
            if (str2bc_rc(bc_str, bc_C) && bc2idx.count(bc_C.to_ulong()) > 0){
                pass_C = true;
            }
            bc_str[npos] = 'G';
            if (str2bc_rc(bc_str, bc_G) && bc2idx.count(bc_G.to_ulong()) > 0){
                pass_G = true;
            }
            bc_str[npos] = 'T';
            if (str2bc_rc(bc_str, bc_T) && bc2idx.count(bc_T.to_ulong()) > 0){
                pass_T = true;
            }
            if (pass_A && !pass_C && !pass_G && !pass_T){
                return bc2idx[bc_A.to_ulong()];
            }
            else if (pass_C && !pass_A && !pass_G && !pass_T){
                return bc2idx[bc_C.to_ulong()];
            }
            else if (pass_G && !pass_A && !pass_C && !pass_T){
                return bc2idx[bc_G.to_ulong()];
            }
            else if (pass_T && !pass_A && !pass_C && !pass_G){
                return bc2idx[bc_T.to_ulong()];
            }
            else{
                return -1;
            }
        }
    }
    else{
        if (str2bc(cur_bc, bc_binary)){
            if (bc2idx.count(bc_binary.to_ulong()) > 0){
                return bc2idx[bc_binary.to_ulong()];
            }
            else{
                return -1;
            }
        }
        else{
            // Contains N.
            // Only allow up to one.
            int npos = -1;
            for (int i = 0; i < BC_LENX2/2; ++i){
                if (cur_bc[i] == 'N'){
                    if (npos != -1){
                        // Multiple Ns; bail.
                        return -1;
                    }
                    else{
                        npos = i;
                    }
                }
            }
            
            strncpy(&bc_str[0], cur_bc, BC_LENX2/2);
            bc_str[BC_LENX2/2] = '\0';

            // Try to mutate to every possible letter.
            bc bc_A;
            bc bc_C;
            bc bc_G;
            bc bc_T;

            bool pass_A = false;
            bool pass_C = false;
            bool pass_G = false;
            bool pass_T = false;

            bc_str[npos] = 'A';
            if (str2bc(bc_str, bc_A) && bc2idx.count(bc_A.to_ulong()) > 0){
                pass_A = true;
            } 
            bc_str[npos] = 'C';
            if (str2bc(bc_str, bc_C) && bc2idx.count(bc_C.to_ulong()) > 0){
                pass_C = true;
            }
            bc_str[npos] = 'G';
            if (str2bc(bc_str, bc_G) && bc2idx.count(bc_G.to_ulong()) > 0){
                pass_G = true;
            }
            bc_str[npos] = 'T';
            if (str2bc(bc_str, bc_T) && bc2idx.count(bc_T.to_ulong()) > 0){
                pass_T = true;
            }
            if (pass_A && !pass_C && !pass_G && !pass_T){
                return bc2idx[bc_A.to_ulong()];
            }
            else if (pass_C && !pass_A && !pass_G && !pass_T){
                return bc2idx[bc_C.to_ulong()];
            }
            else if (pass_G && !pass_A && !pass_C && !pass_T){
                return bc2idx[bc_G.to_ulong()];
            }
            else if (pass_T && !pass_A && !pass_C && !pass_G){
                return bc2idx[bc_T.to_ulong()];
            }
            else{
                return -1;
            }
        }
    }
    // If we've made it here, there are no Ns and the barcode does not fully match
    // an existing barcode.
    // Check for k-mer matches.

    return -1;
}
