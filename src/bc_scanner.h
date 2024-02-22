#ifndef _HTSWRAPPER_BC_SCANNER_H
#define _HTSWRAPPER_BC_SCANNER_H
#include <utility>
#include <cstdlib>
#include <set>
#include <unordered_set>
#include <bitset>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>
#include <htslib/kseq.h>
#include "bc.h"

KSEQ_INIT(gzFile, gzread);

/**
 * This file contains a class designed to scan FASTQ files for reads
 * matching a barcode whitelist.
 */

class bc_scanner{
    private:
        
        // Store all FASTQ files
        gzFile seq1_fp;
        gzFile seq2_fp;
        gzFile seq3_fp;

        kseq_t* kseq1;
        kseq_t* kseq2;
        kseq_t* kseq3;
        
        // Are barcodes reverse-complemented?
        bool rc;
        
        // Are barcodes at the end of the read?
        bool at_end;

        bool has_r2;
        bool has_r3;
        
        // Which file idx are the barcodes in?
        int file_idx;

        bc_whitelist wl;

        // Is this a 2-part (i.e. 10X multiome) whitelist, where we look up
        // in the second and return the value in the first?
        bool use_wl2;
        
        // Should we trim barcodes off reads?
        bool trim_bc;

        bool initialized;
        
        void set_defaults();

    public:
        
        // Allow users to access sequence data
        char* seq_id;
        int seq_id_len;

        char* barcode_read;
        char* barcode_read_qual;
        int barcode_read_len;
        
        char* read1;
        char* read1_qual;
        int read1_len;

        char* read2;
        char* read2_qual;
        int read2_len;

        bool has_read1;
        bool has_read2;
        
        unsigned long barcode;

        // Constructors - tell it about FASTQ files
        bc_scanner(std::string filename);
        bc_scanner(std::string filename1, std::string filename2);
        bc_scanner(std::string filename1, std::string filename2, std::string filename3);
        
        // Destructor
        ~bc_scanner();

        // tell it where to look for barcodes
        void init(std::string wlfile1, std::string wlfile2,
            int file_idx, bool at_end, bool rc, bool wl2);
        
        // Pre-sets
        void init_10x_RNA(std::string wlfile);
        void init_10x_ATAC(std::string wlfile);
        void init_10x_multiome_RNA(std::string wlfile, std::string wlfile2);
        void init_10x_multiome_ATAC(std::string wlfile, std::string wlfile2);
        void init_10x_featureBarcode(std::string wlfile);

        // Tell whether or not to remove barcodes.
        void trim_barcodes(bool trim);

        // Iterate and return false if done, or true if a next seq is available.
        // Seqs are set to public variables.
        bool next();
};


#endif
