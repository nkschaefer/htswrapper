#include <stdio.h>   
#include <stdlib.h> 
#include <ctype.h>
#include <string.h>
#include <map>
#include <string>
#include <algorithm>
#include "bc.h"
#include "khashtable.h"

using namespace std;

string revcomp(const string& seq){
    char seqrc[seq.length() + 1];
    seqrc[seq.length()] = '\0';
    int ind_rc = 0;
    for (int i = seq.length()-1; i >= 0; --i){
        switch(seq[i]){
            case 'A':
                seqrc[ind_rc] = 'T';
            break;
            case 'C':
                seqrc[ind_rc] = 'G';
            break;
            case 'G':
                seqrc[ind_rc] = 'C';
            break;
            case 'T':
                seqrc[ind_rc] = 'A';
            break;
            default:
                seqrc[ind_rc] = 'N';
            break;
        }
        ++ind_rc;
    }
    return string(seqrc);
}

int main(int argc, char *argv[]) {   
    
    // Decide value of k to use
    int k = 44; 
    
    khashtable<int> tab(k);

    khashtable<int> tab2(23);
    
    fprintf(stderr, "===== Testing insertion and lookup: =====\n\n");
    string seq = "ACGTGGGGTAGAGTAGGACCGCGATTTGAGATGCGACCTACTCTACCCCACGTGGTCGCATCTCAAATCGCGG";
    // Create null-terminated buffer to store k-mers as we read from the sequence
    char kmer[k + 1];
    kmer[k] = '\0';
    
    for (int i = 0; i < seq.length() - k + 1; ++i){
        strncpy(kmer, &seq[i], k);
        tab.add(kmer, i);
    }
    for (int i = 0; i < seq.length()-23 + 1; ++i){
        tab2.add(seq.c_str() + i, i);
    }
    for (int i = 0; i < seq.length() - k +1; ++i){
        strncpy(kmer, &seq[i], k);
        khashkey hkey(kmer, k);
        int val;
        if (tab.lookup(hkey, val)){
            fprintf(stderr, "success %s %d\n", kmer, val);
        } 
        else{
            fprintf(stderr, "failure %s %d\n", kmer, val);
        }
    }
    
    fprintf(stderr, "load factor = %f\n", tab.lf());

    fprintf(stderr, "\n===== Testing khashkey::scan_kmers() =====\n");
    
    khashkey hkey(k);
    int pos = 0;
    while(hkey.scan_kmers(seq, pos)){
        hkey.print(false);
        hkey.print(true);
        int val = 0;
        if (tab.lookup(hkey, val)){
            fprintf(stderr, "success2 %d\n", val);
        }
        else{
            fprintf(stderr, "failure2 %d\n", val);
        }
    }

    fprintf(stderr, "\n===== Testing khashkey::scan_kmers() with Ns =====\n");
    string seq2 = "ACGTGGGGTAGAGTAGNACCGCGATTTGAGATGCGACCTACTCTNNCCCACGTGGTCGCATCTCAAATCGCGG";
    
    khashkey hkey2(23);
    pos = 0;
    while(hkey2.scan_kmers(seq2, pos)){
        int val = 0;
        if (tab2.lookup(hkey2, val)){
            fprintf(stderr, "success3 %d\n", val);
        }
        else{
            fprintf(stderr, "failure3 %d\n", val);
            fprintf(stderr, "pos %d\n", pos);
            hkey2.print(false);
        }
    }
}
