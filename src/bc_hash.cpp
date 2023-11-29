#include <cstdlib>
#include <string.h>
#include <utility>
#include <bitset>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <zlib.h>
#include "bc_hash.h"

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
 * Convert a bitset representation of a DNA sequence
 * to a string representation
 */
string bc2str(bc& this_bc, int len){
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

/**
 * Convert a bitset representation of a DNA sequence
 * to a string representation, in reverse complement
 * orientation.
 */
string bc2str_rc(bc& this_bc, int len){
    char strbuf[len+1];
    strbuf[len] = '\0';
    int buf_idx = 0;
    for (int i = len-1; i >= 0; --i){
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
 * Assumes barcodes are length 16.
 *
 */
unsigned long convert_from_barcode_list(string& barcode){
    int ntrim = 0;
    for (int i = barcode.length()-1; i >= 0; ++i){
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
    str2bc(barcode.c_str(), as_bitset, 16);
    return as_bitset.to_ulong();
}

/**
 * Load a filtered list of cell barcodes, gzipped or not.
 * 
 * Trims any trailing stuff (i.e. -1) off the end of the actual
 * barcode sequence, compresses the barcode, and interprets it 
 * as an unsigned long.
 *
 * File names ending with ".gz" are assumed to be gzipped;
 * File names not ending with ".gz" are assumed to be uncompressed.
 */
void parse_barcode_file(string& filename, set<unsigned long>& cell_barcodes){
    
    if (filename.length() > 3 && filename.substr(filename.length()-3, 3) == ".gz"){
        fprintf(stderr, "gzipped\n");
        // Process gzipped file
        int bufsize = 1024;
        char buf[bufsize];
        char strbuf[100];
        gzFile inf = gzopen(filename.c_str(), "r");
        int eof = 0;
        int idx_start = 0;
        while (!eof){
            int nread = gzread(inf, &buf[idx_start], bufsize-idx_start);
            eof = gzeof(inf);
            if (nread < bufsize-idx_start){
                eof = 1;
            }
            // Parse lines.
            int line_start = 0;
            for (int i = 0; i < nread; ++i){
                if (buf[i] == '\n'){
                    strncpy(&strbuf[0], &buf[line_start], i-line_start-1);
                    string bc_str = strbuf;
                    fprintf(stderr, "bc_str %s\n", bc_str.c_str());
                    cell_barcodes.insert(convert_from_barcode_list(bc_str));
                    line_start = i + 1;
                }
            }
            if (eof && line_start < nread){
                // Get last bit
                strncpy(&strbuf[0], &buf[line_start], nread-line_start);
                string bc_str = strbuf;
                fprintf(stderr, "last bc_str %s\n", bc_str.c_str());
                cell_barcodes.insert(convert_from_barcode_list(bc_str));
            }
            else if (line_start < bufsize){
                // Need to copy what remains in buffer to beginning.
                memmove(&buf[0], &buf[line_start], bufsize-line_start);
                idx_start = bufsize-line_start;
            }
            else{
                idx_start = 0;
            }
        }
        gzclose(inf);
    }
    else{ 
        // Read file the convenient way
        ifstream infile(filename.c_str());
        string bc_str;
        while (infile >> bc_str){
            cell_barcodes.insert(convert_from_barcode_list(bc_str));
        }    
    }
    fprintf(stderr, "Read %ld barcodes from file\n", cell_barcodes.size());
}

