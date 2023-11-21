#ifndef _BAMWRAPPER
#include <string>
#include <deque>
#include <vector>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hts.h>

class bam_reader{
    private:
        std::string filename;
        samFile* fp;
        bam_hdr_t* header;
        std::string bc_tag;
        int prevtid;
        bool cb; // Should we look for cell barcode tags for every entry?
        bool cb_raw; // Should we look for raw cell barcode / quality for every entry?
        void set_defaults();
        int32_t get_query_start();
        int32_t get_query_end();
        bool genes; // Should we look for 10x GX and GN tags to tell what gene we're in? 
        hts_itr_t* itr;
        hts_idx_t* idx;
        bool idx_init;
        bool region_set;
        bool initialized;

        void fill_data();

    public:
        bam1_t* reader;
        std::string bc;
        
        int n_header_seqs;
        char** header_seqs;

        long int reference_start;
        long int reference_end;
        long int reference_length;
        long int query_start;
        long int query_end;
        long int query_length;
    
        int mapq;        

        // Parse gene tags from 10x
        std::vector<std::string> gene_ids;
        std::vector<std::string> gene_names;

        // Define stuff to store all potential 10X barcode data
        char* bc_z;
        char* bx_z;
        char* st_z;
        char* rx_z;
        char* qx_z;
        char* tr_z;
        char* tq_z;
        char* cb_z;
        char* cr_z;
        char* cy_z;
        char* qt_z;

        bool has_bc_z;
        bool has_bx_z;
        bool has_st_z;
        bool has_rx_z;
        bool has_qx_z;
        bool has_tr_z;
        bool has_tq_z;
        bool has_cb_z;
        bool has_cr_z;
        bool has_cy_z;
        bool has_qt_z;

        long int isize;
        bool has_bc_tag;
        
    // Constructor/destructor
    bam_reader();
    bam_reader(std::string&);
    ~bam_reader();
    
    void set_file(std::string& filename);
    void clear_read_groups_hdr();    
    void add_read_group_hdr(const std::string&, const std::string&, 
        const std::string&, const std::string&, const std::string&);
    void add_read_group_read(const std::string&);
    void print_contigs_vcf();
    void get_contigs_vcf(std::vector<std::string>&);
    void set_bc_tag(std::string&);
    void set_cb();
    void set_cb_raw();
    void check_genes();

    bool set_query_site(int, long int);    
    bool set_query_region(const char*, long int, long int);
    void unset_query_region();

    bool next();
    
    std::map<std::string, int> get_seq2tid();

    // Flag functions
    bool paired();
    bool proper_pair();
    bool unmapped();
    bool mate_unmapped();
    bool reverse();
    bool mate_reverse();
    bool read1();
    bool read2();
    bool secondary();
    bool qcfail();
    bool dup();
    bool supplementary();
    
    char* ref_id();
    int tid();
    bool chimeric();
    char* read_id();
    int map_start();
    int map_end();
    void get_seq(char*);
    void get_qual(char*);
    void write_header(BGZF*);
    void write_record(BGZF*);
    char get_base_at(long int);
    int get_pos_in_read(long int); 
};

#define _BAMWRAPPER
#endif
