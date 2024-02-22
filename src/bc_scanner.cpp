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

//KSEQ_INIT(gzFile, gzread);

bc_scanner::bc_scanner(string seqfile){
    set_defaults();
    seq1_fp = gzopen(seqfile.c_str(), "r");
    if (!seq1_fp){
        fprintf(stderr, "ERROR opening %s\n", seqfile.c_str());
        exit(1);
    } 
    kseq1 = kseq_init(seq1_fp);
    has_r2 = false;
    has_r3 = false;
}

bc_scanner::bc_scanner(string seq1file, string seq2file){
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
    set_defaults();
}

bc_scanner::bc_scanner(string seq1file, string seq2file, string seq3file){
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
    set_defaults();
}

void bc_scanner::set_defaults(){
    read1 = NULL;
    read2 = NULL;
    barcode_read = NULL;
    read1_qual = NULL;
    read2_qual = NULL;
    barcode_read_qual = NULL;
    read1_len = 0;
    read2_len = 0;
    barcode_read_len = 0;
    has_read1 = false;
    has_read2 = false;
    barcode = 0;
    use_wl2 = false;
    trim_bc = false;
}

bc_scanner::~bc_scanner(){
    kseq_destroy(kseq1);
    if (has_r2){
        kseq_destroy(kseq2);
    }
    if (has_r3){
        kseq_destroy(kseq3);
    }
}

void bc_scanner::init(string whitelistfile, string whitelistfile2, 
    int file_idx, bool at_end, bool rc, bool wl2){
    
    if (file_idx > 2 || file_idx < 0){
        fprintf(stderr, "ERROR: illegal file index %d\n", file_idx);
        exit(1);
    }
    if (!has_r2 && file_idx > 0){
        fprintf(stderr, "ERROR: illegal file index %d with only 1 read\n", file_idx);
        exit(1);
    }
    if (!has_r3 && file_idx > 1){
        fprintf(stderr, "ERROR: illegal file index %d with only 2 reads\n", file_idx);
        exit(1);
    }
    if (wl2){
        wl.init(whitelistfile, whitelistfile2);
    }
    else{
        wl.init(whitelistfile);
    }

    this->file_idx = file_idx;
    this->at_end = at_end;
    this->rc = rc;
    this->use_wl2 = wl2; 
    
    this->seq_id = kseq1->name.s;
    this->seq_id_len = 0;

    if (file_idx == 0){
        barcode_read = &kseq1->seq.s[0];
        barcode_read_qual = &kseq1->qual.s[0];
        has_read1 = false;
        has_read2 = false;
        if (has_r2){
            has_read1 = true;
            read1 = &kseq2->seq.s[0];
            read1_qual = &kseq2->qual.s[0];
        }
        if (has_r3){
            has_read2 = true;
            read2 = &kseq3->seq.s[0];
            read2_qual = &kseq3->qual.s[0];
        }
    }
    else if (file_idx == 1){
        barcode_read = &kseq2->seq.s[0];
        barcode_read_qual = &kseq2->qual.s[0];
        read1 = &kseq1->seq.s[0];
        read1_qual = &kseq1->qual.s[0];
        has_read1 = true;
        if (has_r3){
            has_read2 = true;
            read2 = &kseq3->seq.s[0];
            read2_qual = &kseq3->qual.s[0];
        }
    }
    else{
        barcode_read = &kseq3->seq.s[0];
        barcode_read_qual = &kseq3->qual.s[0];
        read1 = &kseq1->seq.s[0];
        read1_qual = &kseq1->qual.s[0];
        read2 = &kseq2->seq.s[0];
        read2_qual = &kseq2->qual.s[0];
        has_read1 = true;
        has_read2 = true;
    }

    initialized = true;
}

void bc_scanner::init_10x_RNA(string wlfile){
    init(wlfile, "", 0, false, false, false);    
}

void bc_scanner::init_10x_multiome_RNA(string wlfile, string wlfile2){
    init(wlfile, wlfile2, 0, false, false, false);
}

void bc_scanner::init_10x_featureBarcode(string wlfile){
    init(wlfile, "", 0, false, false, false);
}

void bc_scanner::init_10x_ATAC(string wlfile){
    init(wlfile, "", 1, true, true, false);
}

void bc_scanner::init_10x_multiome_ATAC(string wlfile, string wlfile2){
    init(wlfile, wlfile2, 1, true, true, true);
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
        if (file_idx == 0){
            lookupseq = kseq1->seq.s;
        }
        else if (file_idx == 1){
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
            seq_id = kseq1->name.s;
            seq_id_len = kseq1->name.l;            
            if (file_idx == 0){
                barcode_read = kseq1->seq.s;
                barcode_read_qual = kseq1->qual.s;
                barcode_read_len = kseq1->seq.l;
                if (has_r2){
                    read1 = kseq2->seq.s;
                    read1_qual = kseq2->qual.s;
                    read1_len = kseq2->seq.l;
                }
                if (has_r3){
                    read2 = kseq3->seq.s;
                    read2_qual = kseq3->qual.s;
                    read2_len = kseq3->seq.l;
                }
            }
            else if (file_idx == 1){
                barcode_read = kseq2->seq.s;
                barcode_read_qual = kseq2->qual.s;
                barcode_read_len = kseq2->seq.l;
                read1 = kseq1->seq.s;
                read1_qual = kseq1->qual.s;
                read1_len = kseq1->seq.l;
                if (has_r3){
                    read2 = kseq3->seq.s;
                    read2_qual = kseq3->qual.s;
                    read2_len = kseq3->seq.l;
                }
            }
            else{
                barcode_read = kseq3->seq.s;
                barcode_read_qual = kseq3->qual.s;
                barcode_read_len = kseq3->seq.l;
                read1 = kseq1->seq.s;
                read1_qual = kseq1->qual.s;
                read1_len = kseq1->seq.l;
                read2 = kseq2->seq.s;
                read2_qual = kseq2->qual.s;
                read2_len = kseq2->seq.l;
            }
            
        }
    }
    return has_next;
}
    
