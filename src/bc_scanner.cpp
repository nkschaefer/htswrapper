#include <utility>
#include <cstdlib>
#include <set>
#include <unordered_set>
#include <bitset>
#include <string>
#include <unordered_map>
#include <vector>
#include <htslib/kseq.h>
#include <zlib.h>
#include "bc.h"
#include "bc_scanner.h"

using namespace std;

void bc_scanner::set_defaults(){
    kseq1 = NULL;
    kseq2 = NULL;
    kseq3 = NULL;
    read_f = NULL;
    read_r = NULL;
    barcode_read = NULL;
    read_f_qual = NULL;
    read_r_qual = NULL;
    barcode_read_qual = NULL;
    read_f_len = 0;
    read_r_len = 0;
    barcode_read_len = 0;
    has_read_f = false;
    has_read_r = false;
    barcode = 0;
    use_wl2 = false;
    trim_bc = false;
    file_idx_bc = -1;
    file_idx_umi = -1;
    umi_len = 0;
    umi_start = -1;
    has_umi = false;
    umi = NULL;
    initialized = false;
}

void bc_scanner::add_reads(string seqfile){
    if (!initialized){
        fprintf(stderr, "ERROR: bc_scanner not initialized; cannot add reads\n");
        exit(1);
    }
    if (file_idx_bc > 0){
        fprintf(stderr, "ERROR: illegal file index %d with only 1 read\n", file_idx_bc);
        exit(1);
    }
    if (file_idx_umi > 0){
        fprintf(stderr, "ERROR: illegal UMI file index %d with only 1 read\n", file_idx_umi);
        exit(1);
    }
    close_seqs();
    seq1_fp = gzopen(seqfile.c_str(), "r");
    if (!seq1_fp){
        fprintf(stderr, "ERROR opening %s\n", seqfile.c_str());
        exit(1);
    } 
    kseq1 = kseq_init(seq1_fp);
    has_r2 = false;
    has_r3 = false;
    kseq2 = NULL;
    kseq3 = NULL;
    set_seq_pointers();
}

void bc_scanner::add_reads(string seq1file, string seq2file){
    if (!initialized){
        fprintf(stderr, "ERROR: bc_scanner not initialized; cannot add reads\n");
        exit(1);
    }
    if (file_idx_bc > 1){
        fprintf(stderr, "ERROR: illegal file index %d with only 2 reads\n", file_idx_bc);
        exit(1);
    }
    if (file_idx_umi > 1){
        fprintf(stderr, "ERROR: illegal file index %d with only 2 reads\n", file_idx_umi);
    }
    close_seqs();
    seq1_fp = gzopen(seq1file.c_str(), "r");
    if (!seq1_fp){
        fprintf(stderr, "ERROR opening %s\n", seq1file.c_str());
        exit(1);
    }
    seq2_fp = gzopen(seq2file.c_str(), "r");
    if (!seq2_fp){
        fprintf(stderr, "ERROR opening %s\n", seq2file.c_str());
        exit(1);
    }
    kseq1 = kseq_init(seq1_fp);
    kseq2 = kseq_init(seq2_fp);
    has_r2 = true;
    has_r3 = false;
    kseq3 = NULL; 
    set_seq_pointers();
}

void bc_scanner::add_reads(string seq1file, string seq2file, string seq3file){
    if (!initialized){
        fprintf(stderr, "ERROR: bc_scanner not initialized; cannot add reads\n");
        exit(1);
    }
    close_seqs();
    seq1_fp = gzopen(seq1file.c_str(), "r");
    if (!seq1_fp){
        fprintf(stderr, "ERROR opening %s\n", seq1file.c_str());
        exit(1);
    }
    seq2_fp = gzopen(seq2file.c_str(), "r");
    if (!seq2_fp){
        fprintf(stderr, "ERROR opening %s\n", seq2file.c_str());
        exit(1);
    }
    seq3_fp = gzopen(seq3file.c_str(), "r");
    if (!seq3_fp){
        fprintf(stderr, "ERROR opening %s\n", seq3file.c_str());
        exit(1);
    }
    kseq1 = kseq_init(seq1_fp);
    kseq2 = kseq_init(seq2_fp);
    kseq3 = kseq_init(seq3_fp);
    has_r2 = true;
    has_r3 = true;
    set_seq_pointers();
}

bc_scanner::bc_scanner(){
    set_defaults();
}
bc_scanner::bc_scanner(string seqfile){
    add_reads(seqfile);
    set_defaults();
}
bc_scanner::bc_scanner(string seq1file, string seq2file){
    add_reads(seq1file, seq2file);
    set_defaults();
}

bc_scanner::bc_scanner(string seq1file, string seq2file, string seq3file){
    add_reads(seq1file, seq2file, seq3file);
    set_defaults();
}

void bc_scanner::close_seqs(){
    if (kseq1 != NULL){
        kseq_destroy(kseq1);
        kseq1 = NULL;
    }
    if (kseq2 != NULL){
        kseq_destroy(kseq2);
        kseq2 = NULL;
    }
    if (kseq3 != NULL){
        kseq_destroy(kseq3);
        kseq3 = NULL;
    }
}

bc_scanner::~bc_scanner(){
    close_seqs();
    if (umi != NULL){
        free(umi);
    }
}

void bc_scanner::set_seq_pointers(){
    
    this->seq_id = kseq1->name.s;
    this->seq_id_len = kseq1->name.l;

    if (file_idx_bc == 0){
        barcode_read = &kseq1->seq.s[0];
        barcode_read_qual = &kseq1->qual.s[0];
        barcode_read_len = kseq1->seq.l;
        has_read_f = false;
        has_read_r = false;
        if (has_r2){
            has_read_f = true;
            read_f = &kseq2->seq.s[0];
            read_f_qual = &kseq2->qual.s[0];
            read_f_len = kseq2->seq.l;
        }
        if (has_r3){
            has_read_r = true;
            read_r = &kseq3->seq.s[0];
            read_r_qual = &kseq3->qual.s[0];
            read_r_len = kseq3->seq.l;
        }
    }
    else if (file_idx_bc == 1){
        barcode_read = &kseq2->seq.s[0];
        barcode_read_qual = &kseq2->qual.s[0];
        barcode_read_len = kseq2->seq.l;
        read_f = &kseq1->seq.s[0];
        read_f_qual = &kseq1->qual.s[0];
        has_read_f = true;
        if (has_r3){
            has_read_r = true;
            read_r = &kseq3->seq.s[0];
            read_r_qual = &kseq3->qual.s[0];
            read_r_len = kseq3->seq.l;
        }
    }
    else{
        barcode_read = &kseq3->seq.s[0];
        barcode_read_qual = &kseq3->qual.s[0];
        barcode_read_len = kseq3->seq.l;
        read_f = &kseq1->seq.s[0];
        read_f_qual = &kseq1->qual.s[0];
        read_f_len = kseq1->seq.l;
        read_r = &kseq2->seq.s[0];
        read_r_qual = &kseq2->qual.s[0];
        read_r_len = kseq3->seq.l;
        has_read_f = true;
        has_read_r = true;
    }
}

void bc_scanner::init(string whitelistfile, 
    string whitelistfile2, 
    int file_idx_bc, 
    bool at_end, 
    bool rc, 
    bool wl2, 
    int file_idx_umi, 
    int umi_start, 
    int umi_len,
    int bc_len){
    
    // Determine optimal value of k given barcode length.
    // k must be < (bc_length + 1)/2
    // Maximal possible k = best performance (fewer expensive k-mer lookups)
    int k;
    if (bc_len/2 % 2 == 0){
        // Even barcode length.
        k = bc_len/2;
    }
    else{
        k = (bc_len-1)/2;
    }
    
    if (initialized){
        // Can only initialize once.
        fprintf(stderr, "ERROR: bc_scanner already initialized\n");
        exit(1);
    }
    if (file_idx_bc > 2 || file_idx_bc < 0){
        fprintf(stderr, "ERROR: illegal file index %d\n", file_idx_bc);
        exit(1);
    }
    if (file_idx_umi > 2){
        fprintf(stderr, "ERROR: illegal UMI file index %d\n", file_idx_umi);
        exit(1);
    }
    if (wl2){
        wl.init(whitelistfile, whitelistfile2, bc_len, k);
    }
    else{
        wl.init(whitelistfile, bc_len, k);
    }
    
    this->file_idx_bc = file_idx_bc;
    this->at_end = at_end;
    this->rc = rc;
    this->use_wl2 = wl2; 
    
    this->file_idx_umi = file_idx_umi;
    if (file_idx_umi >= 0){
        // Create a buffer to store UMI sequence
        umi = (char*)malloc((umi_len+1)*sizeof(char));
        this->umi_start = umi_start;
        this->umi_len = umi_len;
    }
    else{
        umi = NULL;
        this->umi_start = -1;
        this->umi_len = 0;
    }

    initialized = true;
}

// Create pre-sets
// Each needs to know how to read whitelists (one or two)
// Cell barcode: which file, start or end of read, orientation
// UMI: present/not, which file, position in read, length

// Useful reference (and the only place I could find this information compiled 
// somewhere): kallisto/bustools command kb --list

// v2 has 10 bp UMIs instead of 12 bp
void bc_scanner::init_10x_RNA_v2(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 10);
}

void bc_scanner::init_10x_RNA_v3(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

// Assumes V3
void bc_scanner::init_10x_RNA(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

// Assumes V3
// Ignore second whitelist -- first is for RNA
void bc_scanner::init_10x_multiome_RNA(string wlfile, string wlfile2){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

void bc_scanner::init_10x_featureBarcode_v2(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 10);
}

void bc_scanner::init_10x_featureBarcode_v3(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

// Assumes V3
void bc_scanner::init_10x_featureBarcode(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

void bc_scanner::init_10x_ATAC(string wlfile){
    init(wlfile, "", 1, true, true, false, -1, -1, 0);
}

void bc_scanner::init_10x_multiome_ATAC(string wlfile, string wlfile2){
    init(wlfile, wlfile2, 1, true, true, true, -1, -1, 0);
}

void bc_scanner::init_multiseq_v2(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 10);
}

void bc_scanner::init_multiseq_v3(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

// Assumes v3
void bc_scanner::init_multiseq(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

bool bc_scanner::next(){
    bool has_next = true;
    bool bc_found = false;

    while (has_next && !bc_found){
        int progress = kseq_read(kseq1);
        if (progress < 0){
            has_next = false;
        }
        if (has_r2){
            int prog2 = kseq_read(kseq2);
            if (strcmp(kseq2->name.s, kseq1->name.s) != 0){
                fprintf(stderr, "ERROR: ID mismatch: %s vs %s\n", kseq1->name.s,
                    kseq2->name.s);
                exit(1);
            }
            if (prog2 < 0){
                has_next = false;
            }
        }
        if (has_r3){
            int prog3 = kseq_read(kseq3);
            if (strcmp(kseq3->name.s, kseq1->name.s) != 0){
                fprintf(stderr, "ERROR: ID mismatch: %s vs %s\n", kseq1->name.s,
                    kseq3->name.s);
                exit(1);
            }
            if (prog3 < 0){
                has_next = false;
            }
        }
        // Check for a barcode.
        char* lookupseq;
        if (file_idx_bc == 0){
            lookupseq = kseq1->seq.s;
        }
        else if (file_idx_bc == 1){
            lookupseq = kseq2->seq.s;
        }
        else{
            lookupseq = kseq3->seq.s;
        }
        if (!at_end && !rc && !use_wl2){
            bc_found = wl.lookup1_bf(lookupseq, barcode, trim_bc); 
        }
        else if (!at_end && !rc && use_wl2){
            bc_found = wl.lookup2_bf(lookupseq, barcode, trim_bc);
        }
        else if (!at_end && rc && !use_wl2){
            bc_found = wl.lookup1_br(lookupseq, barcode, trim_bc);
        }
        else if (!at_end && rc && use_wl2){
            bc_found = wl.lookup2_br(lookupseq, barcode, trim_bc);
        }
        else if (at_end && !rc && !use_wl2){
            bc_found = wl.lookup1_ef(lookupseq, barcode, trim_bc);
        }
        else if (at_end && !rc && use_wl2){
            bc_found = wl.lookup2_ef(lookupseq, barcode, trim_bc);
        }
        else if (at_end && rc && !use_wl2){
            bc_found = wl.lookup1_er(lookupseq, barcode, trim_bc);
        }
        else if (at_end && rc && use_wl2){
            bc_found = wl.lookup2_er(lookupseq, barcode, trim_bc);
        }
        if (bc_found){
            // Set all sequence pointers
            set_seq_pointers();
            // Copy UMI from appropriate seq
            if (file_idx_umi == 0){
                if (kseq1->seq.l >= umi_start + umi_len){
                    strncpy(&umi[0], kseq1->seq.s + umi_start, umi_len); 
                    umi[umi_len] = '\0';
                    has_umi = true;
                }
                else{
                    umi[0] = '\0';
                    has_umi = false;
                }
            }
            else if (file_idx_umi == 1){
                if (kseq2->seq.l >= umi_start + umi_len){
                    strncpy(&umi[0], kseq2->seq.s + umi_start, umi_len);
                    umi[umi_len] = '\0';
                    has_umi = true;
                }
                else{
                    umi[0] = '\0';
                    has_umi = false;
                }
            }
            else if (file_idx_umi == 2){
                if (kseq3->seq.l >= umi_start + umi_len){
                    strncpy(&umi[0], kseq3->seq.s + umi_start, umi_len);
                    umi[umi_len] = '\0';
                    has_umi = true;
                }
                else{
                    umi[0] = '\0';
                    has_umi = false;
                }
            }
            else{
                has_umi = false;
            }
            /*
            seq_id = kseq1->name.s;
            seq_id_len = kseq1->name.l;            
            if (file_idx == 0){
                barcode_read = kseq1->seq.s;
                barcode_read_qual = kseq1->qual.s;
                barcode_read_len = kseq1->seq.l;
                if (has_r2){
                    read_f = kseq2->seq.s;
                    read_f_qual = kseq2->qual.s;
                    read_f_len = kseq2->seq.l;
                }
                if (has_r3){
                    read_r = kseq3->seq.s;
                    read_r_qual = kseq3->qual.s;
                    read_r_len = kseq3->seq.l;
                }
            }
            else if (file_idx == 1){
                barcode_read = kseq2->seq.s;
                barcode_read_qual = kseq2->qual.s;
                barcode_read_len = kseq2->seq.l;
                read_f = kseq1->seq.s;
                read_f_qual = kseq1->qual.s;
                read_f_len = kseq1->seq.l;
                if (has_r3){
                    read_r = kseq3->seq.s;
                    read_r_qual = kseq3->qual.s;
                    read_r_len = kseq3->seq.l;
                }
            }
            else{
                barcode_read = kseq3->seq.s;
                barcode_read_qual = kseq3->qual.s;
                barcode_read_len = kseq3->seq.l;
                read_f = kseq1->seq.s;
                read_f_qual = kseq1->qual.s;
                read_f_len = kseq1->seq.l;
                read_r = kseq2->seq.s;
                read_r_qual = kseq2->qual.s;
                read_r_len = kseq2->seq.l;
            }
            */
        }
    }
    return has_next;
}
    
