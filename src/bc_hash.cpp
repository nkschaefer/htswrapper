#include <cstdlib>
#include <string.h>
#include <utility>
#include <bitset>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sys/time.h>
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

