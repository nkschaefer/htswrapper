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

seq_info::seq_info(){
    n = NULL;
    nl = 0;
    seq1 = NULL;
    len1 = 0;
    qual1 = NULL;
    seq2 = NULL;
    len2 = 0;
    qual2 = NULL;
    seq3 = NULL;
    len3 = 0;
    qual3 = NULL;
    barcode = 0;
}
seq_info::seq_info(char* name, int name_l,
        char* s1, int l1, char* q1,
        char* s2, int l2, char* q2,
        char* s3, int l3, char* q3){
    n = name;
    nl = name_l;
    seq1 = s1;
    len1 = l1;
    qual1 = q1;
    seq2 = s2;
    len2 = l2;
    qual2 = q2;
    seq3 = s3;
    len3 = l3;
    qual3 = q3;
    barcode = 0;
}

seq_info::seq_info(seq_info& other){
    this->n = other.n;
    this->nl = other.nl;
    this->seq1 = other.seq1;
    this->len1 = other.len1;
    this->qual1 = other.qual1;
    this->seq2 = other.seq2;
    this->len2 = other.len2;
    this->qual2 = other.qual2;
    this->seq3 = other.seq3;
    this->len3 = other.len3;
    this->qual3 = other.qual3;
    this->barcode = other.barcode;
}


void bc_scanner::set_defaults(){
    kseq1 = NULL;
    kseq2 = NULL;
    kseq3 = NULL;
    
    seq_id = NULL;
    seq_id_len = 0;
     
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
    
    umi = NULL;
    umi_len = 0;

    use_wl2 = false;
    trim_bc = false;
    file_idx_bc = -1;
    file_idx_umi = -1;
    umi_start = -1;
    has_umi = false;
    initialized = false;
    reads_init = false;
    
    wl_copied = false;

    eof = false; 
    nthreads = 0;
    threads_init = false;
    terminate_threads = false;
    jobcount = 0;
    count_missing_bc = 0;
    count_returned = 0;
}

// ----- Thread-related -----
void bc_scanner::set_threads(int nt){
    if (nt < 2){
        nt = 0;
    }
    this->nthreads = nt;
    //this->threads.reserve(nt);
}

void bc_scanner::launch_threads(){
    terminate_threads = false;
    for (int i = 0; i < nthreads; ++i){
        threads.emplace_back(&bc_scanner::worker, this); 
        //thread* t = new thread(&bc_scanner::worker, this);
        //threads.push_back(t);
    }
    
    threads_init = true;
}

void bc_scanner::close_pool(){
    {
        unique_lock<mutex> lock(queue_mutex);
        terminate_threads = true;
    }
    has_jobs.notify_all();
    while (threads.size() > 0){
        threads[0].join();
        threads.pop_front();
    }
    /*
    for (int i = 0; i < threads.size(); i){
        //threads[i]->join();
        //delete threads[i];
        threads[i].join();
    }
    threads.clear();
    */
    threads_init = false;
}

void bc_scanner::add_job(char* n, int nl,
    char* seq1, int seq1len, char* q1, 
    char* seq2, int seq2len, char* q2, 
    char* seq3, int seq3len, char* q3){
    {
        unique_lock<mutex> lock(queue_mutex);
        jobs.emplace_back(n, nl, 
                seq1, seq1len, q1, 
                seq2, seq2len, q2,
                seq3, seq3len, q3);
    }
    has_jobs.notify_one();
}

void bc_scanner::worker(){
    while(true){
        //if (jobs.size() == 0 && terminate_threads){
        //    return;
        //}
        seq_info si;
        {
            unique_lock<mutex> lock(queue_mutex);
            has_jobs.wait(lock, [this]{ return jobs.size() > 0 || terminate_threads; });
            if (jobs.size() == 0 && terminate_threads){
                return;
            }
            si = jobs[0];
            jobs.pop_front();
        }
        unsigned long bc_result = check_bc(si);
        if (bc_result != 0){
            si.barcode = bc_result;
            unique_lock<mutex> lock(output_mutex);
            output.emplace_back(si);
        }
        else{
            if (si.n != NULL){
                free(si.n);
            }
            if (si.seq1 != NULL){
                free(si.seq1);
            }
            if (si.qual1 != NULL){
                free(si.qual1);
            }
            if (si.seq2 != NULL){
                free(si.seq2);
            }
            if (si.qual2 != NULL){
                free(si.qual2);
            }
            if (si.seq3 != NULL){
                free(si.seq3);
            }
            if (si.qual3 != NULL){
                free(si.qual3);
            }
        }
    }
}

void bc_scanner::trim_barcodes(bool trim){
    this->trim_bc = trim;
    if (trim){
        this->corr_bc = false;
    }
}

void bc_scanner::corr_barcodes(bool corr){
    this->corr_bc = corr;
    if (corr){
        this->trim_bc = false;
    }
}

void bc_scanner::add_reads(string seqfile){
    if (initialized){
        if (file_idx_bc > 0){
            fprintf(stderr, "ERROR: illegal file index %d with only 1 read\n", file_idx_bc);
            exit(1);
        }
        if (file_idx_umi > 0){
            fprintf(stderr, "ERROR: illegal UMI file index %d with only 1 read\n", file_idx_umi);
            exit(1);
        }
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
    reads_init = true;
}

void bc_scanner::add_reads(string seq1file, string seq2file){
    if (initialized){
        if (file_idx_bc > 1){
            fprintf(stderr, "ERROR: illegal file index %d with only 2 reads\n", file_idx_bc);
            exit(1);
        }
        if (file_idx_umi > 1){
            fprintf(stderr, "ERROR: illegal file index %d with only 2 reads\n", file_idx_umi);
        }
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
    reads_init = true;
}

void bc_scanner::add_reads(string seq1file, string seq2file, string seq3file){
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
    reads_init = true;
}

bc_scanner::bc_scanner(){
    set_defaults();
}
bc_scanner::bc_scanner(string seqfile){
    set_defaults();
    add_reads(seqfile);
}
bc_scanner::bc_scanner(string seq1file, string seq2file){
    set_defaults();
    add_reads(seq1file, seq2file);
}

bc_scanner::bc_scanner(string seq1file, string seq2file, string seq3file){
    set_defaults();
    add_reads(seq1file, seq2file, seq3file);
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
    reads_init = false;
}

/**
 * Destructor
 */
bc_scanner::~bc_scanner(){
    close_seqs();
    if (umi != NULL){
        free(umi);
    }
    if (!wl_copied){
        delete wl;
    }
    if (nthreads > 1){
        if (seq_id != NULL){
            free(seq_id);
        }
        if (barcode_read != NULL){
            free(barcode_read);
        }
        if (barcode_read_qual != NULL){
            free(barcode_read_qual);
        }
        if (read_f != NULL){
            free(read_f);
        }
        if (read_r != NULL){
            free(read_r);
        }
        for (int i = 0; i < threads.size(); ++i){
            threads[i].join();
        }
        threads.clear();
        jobs.clear();
        output.clear();
    }
}

void bc_scanner::exact_matches_only(){
    this->wl->exact_matches_only();
}

void bc_scanner::init_aux(int file_idx_bc, 
    bool at_end, 
    bool rc, 
    int file_idx_umi, 
    int umi_start, 
    int umi_len){
    
    if (file_idx_umi != -1 && umi_len <= 0){
        fprintf(stderr, "ERROR: invalid UMI len %d\n", umi_len);
        exit(1);
    }
    if (file_idx_umi != -1 && umi_start < 0){
        fprintf(stderr, "ERROR: invalid UMI start idx %d\n", umi_start);
        exit(1);
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
    if (reads_init){
        if (!has_r2 && file_idx_bc > 0){
            fprintf(stderr, "ERROR: illegal BC file index with 1 read file\n");
            exit(1);
        }
        else if (!has_r3 && file_idx_bc > 1){
            fprintf(stderr, "ERROR: illegal BC file index with 2 read files\n");
            exit(1);
        }
        else if (!has_r2 && file_idx_umi > 0){
            fprintf(stderr, "ERROR: illegal UMI file index with 1 read file\n");
            exit(1);
        }
        else if (!has_r3 && file_idx_umi > 1){
            fprintf(stderr, "ERROR: illegal UMI file index with 2 read files\n");
            exit(1);
        }
    }
    
    this->file_idx_bc = file_idx_bc;
    this->at_end = at_end;
    this->rc = rc;
    
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

void bc_scanner::init(bc_whitelist& wl, 
    int file_idx_bc, 
    bool at_end, 
    bool rc, 
    int file_idx_umi, 
    int umi_start, 
    int umi_len){
    
    this->wl = &wl;
    this->use_wl2 = wl.two_lists();
    wl_copied = true;
    init_aux(file_idx_bc, at_end, rc, file_idx_umi, umi_start, umi_len);

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
    
    if (bc_len <= 0){
        fprintf(stderr, "ERROR: invalid barcode len %d\n", bc_len);
        exit(1);
    }
    
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
    if (wl2){
        this->wl = new bc_whitelist(whitelistfile, whitelistfile2, bc_len, k); 
    }
    else{
        this->wl = new bc_whitelist(whitelistfile, bc_len, k);
    }
    wl_copied = false;

    this->use_wl2 = wl2;
    init_aux(file_idx_bc, at_end, rc, file_idx_umi, umi_start, umi_len);
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
void bc_scanner::init_10x_RNA_v2(bc_whitelist& w){
    if (w.two_lists()){
        fprintf(stderr, "ERROR: two allowed barcode lists detected, but requires 1.\n");
        exit(1);
    }
    init(w, 0, false, false, 0, 16, 10); 
}

void bc_scanner::init_10x_RNA_v3(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

void bc_scanner::init_10x_RNA_v3(bc_whitelist& w){
    if (w.two_lists()){
        fprintf(stderr, "ERROR: two allowed barcode lists detected, but requires 1.\n");
        exit(1);
    }
    init(w, 0, false, false, 0, 16, 12);
}

// Assumes V3
void bc_scanner::init_10x_RNA(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

void bc_scanner::init_10x_RNA(bc_whitelist& w){
    if (w.two_lists()){
        fprintf(stderr, "ERROR: two allowed barcode lists detected, but requires 1.\n");
        exit(1);
    }
    init(w, 0, false, false, 0, 16, 12);
}

// Assumes V3
// Ignore second whitelist -- first is for RNA
void bc_scanner::init_10x_multiome_RNA(string wlfile, string wlfile2){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

void bc_scanner::init_10x_multiome_RNA(bc_whitelist& w){
    init(w, 0, false, false, 0, 16, 12);
}

void bc_scanner::init_10x_featureBarcode_v2(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 10);
}

void bc_scanner::init_10x_featureBarcode_v2(bc_whitelist& w){
    if (w.two_lists()){
        fprintf(stderr, "ERROR: two allowed barcode lists detected, but requires 1.\n");
        exit(1);
    }
    init(w, 0, false, false, 0, 16, 10);
}

void bc_scanner::init_10x_featureBarcode_v3(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

void bc_scanner::init_10x_featureBarcode_v3(bc_whitelist& w){
    if (w.two_lists()){
        fprintf(stderr, "ERROR: two allowed barcode lists detected, but requires 1.\n");
        exit(1);
    }
    init(w, 0, false, false, 0, 16, 12);
}

// Assumes V3
void bc_scanner::init_10x_featureBarcode(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

void bc_scanner::init_10x_featureBarcode(bc_whitelist& w){
    if (w.two_lists()){
        fprintf(stderr, "ERROR: two allowed barcode lists detected, but requires 1.\n");
        exit(1);
    }
    init(w, 0, false, false, 0, 16, 12);
}

void bc_scanner::init_10x_ATAC(string wlfile){
    init(wlfile, "", 1, true, true, false, -1, -1, 0);
}

void bc_scanner::init_10x_ATAC(bc_whitelist& w){
    if (w.two_lists()){
        fprintf(stderr, "ERROR: two allowed barcode lists detected, but requires 1.\n");
        exit(1);
    }
    init(w, 1, true, true, -1, -1, 0);
}

void bc_scanner::init_10x_multiome_ATAC(string wlfile, string wlfile2){
    init(wlfile, wlfile2, 1, true, true, true, -1, -1, 0);
}

void bc_scanner::init_10x_multiome_ATAC(bc_whitelist& w){
    if (!w.two_lists()){
        fprintf(stderr, "ERROR: one allowed barcode list detected, but requires 2.\n");
        exit(1);
    }
    init(w, 1, true, true, -1, -1, 0);
    use_wl2 = true;
}

void bc_scanner::init_multiseq_v2(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 10);
}

void bc_scanner::init_multiseq_v2(bc_whitelist& w){
    init(w, 0, false, false, 0, 16, 10);
    use_wl2 = false;
}

void bc_scanner::init_multiseq_v3(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

void bc_scanner::init_multiseq_v3(bc_whitelist& w){
    init(w, 0, false, false, 0, 16, 12);
    use_wl2 = false;
}

// Assumes v3
void bc_scanner::init_multiseq(string wlfile){
    init(wlfile, "", 0, false, false, false, 0, 16, 12);
}

void bc_scanner::init_multiseq(bc_whitelist& w){
    init(w, 0, false, false, 0, 16, 12);
    use_wl2 = false;
}

bool bc_scanner::check_bc(){
    char* lookupseq;
    int lookuplen;
    if (file_idx_bc == 0){
        lookupseq = kseq1->seq.s;
        lookuplen = kseq1->seq.l;
    }
    else if (file_idx_bc == 1){
        lookupseq = kseq2->seq.s;
        lookuplen = kseq2->seq.l;
    }
    else{
        lookupseq = kseq3->seq.s;
        lookuplen = kseq3->seq.l;
    }
    return check_bc_aux(lookupseq, lookuplen, barcode);
}

unsigned long bc_scanner::check_bc(seq_info si){
    
    // Check for a barcode.
    char* lookupseq;
    int lookuplen;
    if (file_idx_bc == 0){
        lookupseq = si.seq1;
        lookuplen = si.len1;
    }
    else if (file_idx_bc == 1){
        lookupseq = si.seq2;
        lookuplen = si.len2;
    }
    else{
        lookupseq = si.seq3;
        lookuplen = si.len3;
    }
    unsigned long barcode_check = 0;
    if (!check_bc_aux(lookupseq, lookuplen, barcode_check)){
       barcode_check = 0;
    }
    return barcode_check;
}

bool bc_scanner::check_bc_aux(char* lookupseq, int lookuplen, unsigned long& barcode_retrieve){
    bool bc_found = false;
    if (!at_end && !rc && !use_wl2){
        bc_found = wl->lookup1_bf(lookupseq, lookuplen, barcode_retrieve, trim_bc, corr_bc); 
    }
    else if (!at_end && !rc && use_wl2){
        bc_found = wl->lookup2_bf(lookupseq, lookuplen, barcode_retrieve, trim_bc, corr_bc);
    }
    else if (!at_end && rc && !use_wl2){
        bc_found = wl->lookup1_br(lookupseq, lookuplen, barcode_retrieve, trim_bc, corr_bc);
    }
    else if (!at_end && rc && use_wl2){
        bc_found = wl->lookup2_br(lookupseq, lookuplen, barcode_retrieve, trim_bc, corr_bc);
    }
    else if (at_end && !rc && !use_wl2){
        bc_found = wl->lookup1_ef(lookupseq, lookuplen, barcode_retrieve, trim_bc, corr_bc);
    }
    else if (at_end && !rc && use_wl2){
        bc_found = wl->lookup2_ef(lookupseq, lookuplen, barcode_retrieve, trim_bc, corr_bc);
    }
    else if (at_end && rc && !use_wl2){
        bc_found = wl->lookup1_er(lookupseq, lookuplen, barcode_retrieve, trim_bc, corr_bc);
    }
    else if (at_end && rc && use_wl2){
        bc_found = wl->lookup2_er(lookupseq, lookuplen, barcode_retrieve, trim_bc, corr_bc);
    }
    return bc_found;
}

void bc_scanner::copy_umi(){
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
}
void bc_scanner::copy_umi(seq_info si){
    if (file_idx_umi == 0){
        if (si.seq1 != NULL && si.len1 >= umi_start + umi_len){
            strncpy(&umi[0], si.seq1 + umi_start, umi_len); 
            umi[umi_len] = '\0';
            has_umi = true;
        }
        else{
            umi[0] = '\0';
            has_umi = false;
        }
    }
    else if (file_idx_umi == 1){
        if (si.seq2 != NULL && si.len2 >= umi_start + umi_len){
            strncpy(&umi[0], si.seq2 + umi_start, umi_len);
            umi[umi_len] = '\0';
            has_umi = true;
        }
        else{
            umi[0] = '\0';
            has_umi = false;
        }
    }
    else if (file_idx_umi == 2){
        if (si.seq3 != NULL && si.len3 >= umi_start + umi_len){
            strncpy(&umi[0], si.seq3 + umi_start, umi_len);
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
        read_f_len = kseq1->seq.l;
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
        read_r_len = kseq2->seq.l;
        has_read_f = true;
        has_read_r = true;
    }
}

void bc_scanner::set_seq_pointers(seq_info si){
    barcode = si.barcode;
    if (seq_id != NULL){
        free(seq_id);
        seq_id = NULL;
    }
    if (barcode_read != NULL){
        free(barcode_read);
        barcode_read = NULL;
    }
    if (barcode_read_qual != NULL){
        free(barcode_read_qual);
        barcode_read_qual = NULL;
    }
    if (read_f != NULL){
        free(read_f);
        read_f = NULL;
    }
    if (read_f_qual != NULL){
        free(read_f_qual);
        read_f_qual = NULL;
    }
    if (read_r != NULL){
        free(read_r);
        read_r = NULL;
    }
    if (read_r_qual != NULL){
        free(read_r_qual);
        read_r_qual = NULL;
    }

    this->seq_id = si.n;
    this->seq_id_len = si.nl;

    if (file_idx_bc == 0){
        barcode_read = si.seq1;
        barcode_read_qual = si.qual1;
        barcode_read_len = si.len1;
        has_read_f = false;
        has_read_r = false;
        if (has_r2){
            has_read_f = true;
            read_f = si.seq2;
            read_f_qual = si.qual2;
            read_f_len = si.len2;
        }
        if (has_r3){
            has_read_r = true;
            read_r = si.seq3;
            read_r_qual = si.qual3;
            read_r_len = si.len3;
        }
    }
    else if (file_idx_bc == 1){
        barcode_read = si.seq2;
        barcode_read_qual = si.qual2;
        barcode_read_len = si.len2;
        read_f = si.seq1;
        read_f_qual = si.qual1;
        read_f_len = si.len1;
        has_read_f = true;
        if (has_r3){
            has_read_r = true;
            read_r = si.seq3;
            read_r_qual = si.qual3;
            read_r_len = si.len3;
        }
    }
    else{
        barcode_read = si.seq3;
        barcode_read_qual = si.qual3;
        barcode_read_len = si.len3;
        read_f = si.seq1;
        read_f_qual = si.qual1;
        read_f_len = si.len1;
        read_r = si.seq2;
        read_r_qual = si.qual2;
        read_r_len = si.len2;
        has_read_f = true;
        has_read_r = true;
    }
}

void bc_scanner::pop_output_queue(){
    // This function also frees whatever data was pointed to at these
    // pointers before
    seq_info si = output[0];
    set_seq_pointers(si);
    copy_umi(si);
    output.pop_front();
}

void bc_scanner::check_uneven_eof(bool term1, bool term2, bool term3){
    // Check to see if any file(s) ended before the others.
    if (term1 || term2 || term3){
        if (has_r2){
            if (has_r3){
                if (term1 && !term2 && !term3){
                    fprintf(stderr, "WARNING: reached end of R1 file, but R2 and R3 files \
still contain reads. Your R1 file is likely truncated or corrupted.\n");
                }
                else if (term1 && term2 && !term3){
                    fprintf(stderr, "WARNING: reached end of R1 and R2 files, but R3 file \
still contains reads. Your R1 & R2 files may be truncated or corrupted.\n");
                }
                else if (term1 && !term2 && term3){
                    fprintf(stderr, "WARNING: reached end of R1 and R3 files, but R2 file \
still contains reads. Your R1 & R3 files may be truncated or corrupted.\n");
                }
                else if (!term1 && term2 && !term3){
                    fprintf(stderr, "WARNING: reached end of R2 file, but R1 and R3 files \
still contain reads. Your R2 file is likely truncated or corrupted.\n");
                }
                else if (!term1 && term2 && term3){
                    fprintf(stderr, "WARNING: reached end of R2 and R3 files, but R1 file \
still contains reads. Your R2 and R3 files may be truncated or corrupted.\n");
                }
                else if (!term1 && !term2 && term3){
                    fprintf(stderr, "WARNING: reached end of R3 file, but R1 and R2 files \
still contain reads. Your R3 file is likely truncated or corrupted.\n");
                }
            }
            else{
                if (term1 && !term2){
                    fprintf(stderr, "WARNING: reached end of R1 file, but R2 file still \
contains reads. Your R2 file is likely truncated or corrupted.\n");
                }
                else if (term2 && !term1){
                    fprintf(stderr, "WARNING: reached end of R2 file, but R1 file still \
contains reads. Your R1 file is likely truncated or corrupted.\n");
                }
            }
        }
    }
}

bool bc_scanner::next(){
    if (nthreads > 1 && !threads_init && !terminate_threads){
        launch_threads();
    }

    bool has_next = true;
    bool bc_found = false;
    
    bool term1 = false;
    bool term2 = false;
    bool term3 = false;
    

    if (eof && nthreads > 1){
        // This will only be the case when multithreading. Should have already closed the pool.
        unique_lock<mutex> lock(output_mutex);
        if (output.size() > 0){
            pop_output_queue();
            return output.size() > 0;
        }
        else{
            return false;
        }
    }
    else{
        while (has_next && !bc_found){
            // Read the next set of sequences from input files.

            int progress = kseq_read(kseq1);
            if (progress < 0){
                has_next = false;
                term1 = true;
            }
            if (has_r2){
                int prog2 = kseq_read(kseq2);
                if (strcmp(kseq2->name.s, kseq1->name.s) != 0){
                    fprintf(stderr, "ERROR: ID mismatch.\n  R1 has read ID: %s\n  R2 has read ID: %s\n", 
                        kseq1->name.s, kseq2->name.s);
                    exit(1);
                }
                if (prog2 < 0){
                    has_next = false;
                    term2 = true;
                }
            }
            if (has_r3){
                int prog3 = kseq_read(kseq3);
                if (strcmp(kseq3->name.s, kseq1->name.s) != 0){
                    fprintf(stderr, "ERROR: ID mismatch.\n  R1/R2 have read ID:  %s\n  R3 has read ID: %s\n", 
                        kseq1->name.s, kseq3->name.s);
                    exit(1);
                }
                if (prog3 < 0){
                    has_next = false;
                    term3 = true;
                }
            }
            
            if (term1 || term2 || term3){
                check_uneven_eof(term1, term2, term3);
            }

            if (nthreads > 1){
                // Queue a job.
                // Create copies of sequences
                char* n = (char*)malloc((kseq1->name.l+1)*sizeof(char));
                strncpy(n, kseq1->name.s, kseq1->name.l);
                int nl = kseq1->name.l;
                char* s1 = (char*)malloc((kseq1->seq.l+1)*sizeof(char));
                strncpy(s1, kseq1->seq.s, kseq1->seq.l);
                s1[kseq1->seq.l] = '\0';
                char* q1 = (char*)malloc((kseq1->seq.l+1)*sizeof(char));
                strncpy(q1, kseq1->qual.s, kseq1->seq.l);
                int l1 = kseq1->seq.l;
                char* s2 = NULL;
                int l2 = 0;
                char* q2 = NULL;
                if (has_r2){
                    s2 = (char*)malloc((kseq2->seq.l+1)*sizeof(char));
                    strncpy(s2, kseq2->seq.s, kseq2->seq.l);
                    s2[kseq2->seq.l] = '\0';
                    l2 = kseq2->seq.l;
                    q2 = (char*)malloc((kseq2->seq.l+1)*sizeof(char));
                    strncpy(q2, kseq2->qual.s, kseq2->seq.l);
                    q2[kseq2->seq.l] = '\0';
                }
                char* s3 = NULL;
                int l3 = 0;
                char* q3 = NULL;
                if (has_r3){
                    s3 = (char*)malloc((kseq3->seq.l+1)*sizeof(char));
                    strncpy(s3, kseq3->seq.s, kseq3->seq.l);
                    s3[kseq3->seq.l] = '\0';
                    l3 = kseq3->seq.l;
                    q3 = (char*)malloc((kseq3->seq.l+1)*sizeof(char));
                    strncpy(q3, kseq3->qual.s, kseq3->seq.l);
                    q3[kseq3->seq.l] = '\0';
                }
                ++jobcount;
                add_job(n, nl, s1, l1, q1, s2, l2, q2, s3, l3, q3);
                
                { 
                    unique_lock<mutex> lock(output_mutex);
                    if (output.size() > 0){
                        pop_output_queue();
                        bc_found = true;             
                        ++count_returned;
                    }
                }

                if (!has_next){
                    eof = true;
                    has_next = true;
                    
                    // Ask all jobs to finish.
                    close_pool();
                    
                    {
                        unique_lock<mutex> lock(output_mutex);
                        if (output.size() > 0){
                            pop_output_queue();
                            bc_found = true;
                            ++count_returned;
                        }
                        else{
                            bc_found = false;
                            has_next = false;
                        }
                    }
                    check_uneven_eof(term1, term2, term3);
                }
            }
            else{
                // Check for a valid barcode.
                bc_found = check_bc();

                if (bc_found){
                    // Set all sequence pointers
                    set_seq_pointers();
                    copy_umi();
                }

                if (!has_next){
                    check_uneven_eof(term1, term2, term3);
                }
            }
        }
        return has_next;
    }
}

