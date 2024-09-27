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
#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include "bc.h"

KSEQ_INIT(gzFile, gzread);

/**
 * This file contains a class designed to scan FASTQ files for reads
 * matching a barcode whitelist.
 */

// For multithreading
struct seq_info{
    unsigned long barcode;
    char* n;
    int nl;
    char* seq1;
    int len1;
    char* qual1;
    char* seq2;
    int len2;
    char* qual2;
    char* seq3;
    int len3;
    char* qual3;
    seq_info();
    seq_info(char* name, int name_l,
            char* s1, int l1, char* q1,
            char* s2, int l2, char* q2,
            char* s3, int l3, char* q3);
    seq_info(seq_info& other);
};

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
        int file_idx_bc;
        int file_idx_umi;
        int umi_start;

        bc_whitelist* wl;
        bool wl_copied;

        // Is this a 2-part (i.e. 10X multiome) whitelist, where we look up
        // in the second and return the value in the first?
        bool use_wl2;
        
        // Should we trim barcodes off reads?
        bool trim_bc;
        
        // Should we replace barcodes with error-corrected barcodes in reads?
        bool corr_bc;

        bool initialized;
        bool reads_init;

        void set_defaults();
        void close_seqs();
        
        // Multithreading related stuff
        int nthreads;
        bool threads_init;
        bool terminate_threads;

        std::deque<seq_info> jobs;
        std::mutex queue_mutex;
        std::deque<seq_info> output;
        std::mutex output_mutex;
        std::condition_variable has_jobs;
        //std::condition_variable has_output;
        std::mutex jobcount_mutex;

        std::deque<std::thread> threads;
        //std::vector<std::thread> threads;

        void launch_threads();
        void close_pool();
        void worker();
        void add_job(char* n, int nl,
            char* s1, int l1, char* q1,
            char* s2, int l2, char* q2,
            char* s3, int l3, char* q3); 
        void pop_output_queue();
        
        long int jobcount;
        long int count_missing_bc;
        long int count_returned;

        bool eof; 
        
        bool check_bc_aux(char* lookupseq, int lookuplen, unsigned long& barcode_retrieve);

        // For non-multithreaded (take from kseq_ts)
        bool check_bc();
        void set_seq_pointers();
        void copy_umi();
        
        // For multithreaded (take from passed parameters)
        unsigned long check_bc(seq_info si);
        void set_seq_pointers(seq_info si);
        void copy_umi(seq_info si);

        void check_uneven_eof(bool term1, bool term2, bool term3);
        
        void init_aux(int file_idx_bc,
            bool at_end,
            bool rc,
            int file_idx_umi = -1,
            int umi_start = -1,
            int umi_len = 0);

    public:
        
        // Allow users to access sequence data
        char* seq_id;
        int seq_id_len;

        // Read containing barcode information
        char* barcode_read;
        char* barcode_read_qual;
        int barcode_read_len;
        
        // Forward (primary non-barcode) read
        char* read_f;
        char* read_f_qual;
        int read_f_len;
    
        // Reverse (secondary non-barcode) read, if it exists
        char* read_r;
        char* read_r_qual;
        int read_r_len;
        
        // Unique molecular identifier (UMI) sequence, if it exists
        char* umi;
        int umi_len;

        bool has_read_f;
        bool has_read_r;
        bool has_umi;

        unsigned long barcode;
        
        // Point it to input FASTQ files 
        void add_reads(std::string filename);
        void add_reads(std::string filename1, std::string filename2);
        void add_reads(std::string filename1, std::string filename2, std::string filename3);
        
        // Default constructor
        bc_scanner();

        // Constructors - tell it about FASTQ files
        bc_scanner(std::string filename);
        bc_scanner(std::string filename1, std::string filename2);
        bc_scanner(std::string filename1, std::string filename2, std::string filename3);
        
        // Destructor
        ~bc_scanner();

        // tell it where to look for barcodes
        void init(bc_whitelist& wl,
            int file_idx_bc,
            bool at_end,
            bool rc,
            int file_idx_umi = -1,
            int umi_start = -1,
            int umi_len = 0);

        void init(std::string wlfile1, 
            std::string wlfile2,
            int file_idx_bc, 
            bool at_end, 
            bool rc, 
            bool wl2, 
            int file_idx_umi = -1, 
            int umi_start = -1,
            int umi_len = 0,
            int bc_len=BC_LENX2/2);
        
        // Pre-sets
        void init_10x_RNA_v2(std::string wlfile);
        void init_10x_RNA_v2(bc_whitelist& w);
        void init_10x_RNA_v3(std::string wlfile);
        void init_10x_RNA_v3(bc_whitelist& w);
        void init_10x_RNA(std::string wlfile);
        void init_10x_RNA(bc_whitelist& w);
        void init_10x_ATAC(std::string wlfile);
        void init_10x_ATAC(bc_whitelist& w);
        void init_10x_multiome_RNA(std::string wlfile, std::string wlfile2);
        void init_10x_multiome_RNA(bc_whitelist& w);
        void init_10x_multiome_ATAC(std::string wlfile, std::string wlfile2);
        void init_10x_multiome_ATAC(bc_whitelist& w);
        void init_10x_featureBarcode_v2(std::string wlfile);
        void init_10x_featureBarcode_v2(bc_whitelist& w);
        void init_10x_featureBarcode_v3(std::string wlfile);
        void init_10x_featureBarcode_v3(bc_whitelist& w);
        void init_10x_featureBarcode(std::string wlfile);
        void init_10x_featureBarcode(bc_whitelist& w);
        void init_multiseq_v2(std::string wlfile);
        void init_multiseq_v2(bc_whitelist& w);
        void init_multiseq_v3(std::string wlfile);
        void init_multiseq_v3(bc_whitelist& w);
        void init_multiseq(std::string wlfile);
        void init_multiseq(bc_whitelist& w); 
        
        // Set to only check for exact matches
        void exact_matches_only();

        // Tell whether or not to remove barcodes.
        void trim_barcodes(bool trim);
        
        // Tell whether or not to correct barcodes.
        void corr_barcodes(bool corr);

        // Iterate and return false if done, or true if a next seq is available.
        // Seqs are set to public variables.
        bool next();
    
        void set_threads(int nt);
};


#endif
