#ifndef SERIALIZE_H
#define SERIALIZE_H
#define _FILE_OFFSET_BITS 64

#include <string>
#include <algorithm>
#include <vector>
#include <set>
#include <iterator>
#include <map>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <utility>
#include <bitset>
#include <mutex>

// Stores data for reading from a file written in serial binary format.
struct instream_info{
    //std::ifstream input;
    FILE* input;
    int buffer_index;
    int bufsize;
    int buf_stored;
    int last_bytes_read;
    char* bytes;
    bool eof;
    bool initialized;
    std::recursive_mutex m;
    
    int64_t prevoffset;
    int64_t prevoffset_rec;
    int prev_start_idx;
    int prev_buf_idx;
    
    bool finished(){ return ((feof(this->input) || this->eof) && 
        (this->buffer_index == this->buf_stored || this->buf_stored == -1)); };
    void print_bufdata(){ fprintf(stderr, "==========\n"); 
        fprintf(stderr, "bufsize %d\n", this->bufsize); 
        fprintf(stderr, "bytes stored %d\n", this->buf_stored); 
        fprintf(stderr, "buffer index %d\n", this->buffer_index); 
        fprintf(stderr, "==========\n");};
    
    // NOTE: the below function should only be called AFTER unserializing DNA.
    // it will report the information just BEFORE the DNA was unserialized, corresponding
    // to the beginning of the record.
    std::pair<int64_t, int> tell(){ return std::make_pair(this->prevoffset, 
        this->prev_buf_idx-this->prev_start_idx); };
    instream_info(){ this->initialized = false; };
    ~instream_info(){ if (this->initialized){ fclose(this->input); free(this->bytes); } };
};

// Stores data for writing to a file in serial binary format. Allows us to lock
// output files when writing to them.
struct outstream_info{
    std::ofstream output;
    std::recursive_mutex m;
    bool initialized;
    outstream_info(){ this->initialized = false; };
    std::streampos tell(){ return this->output.tellp(); };
    void close(){ if (this->initialized){ this->output.close(); this->initialized = false;} };
    ~outstream_info(){ if (this->initialized){ this->output.close(); } };
};

void instream_init(instream_info&, std::string filename, int);
void outstream_init(outstream_info&, std::string filename);
void outstream_close(outstream_info&);
void check_bytes_left(instream_info&, int);
void seek(instream_info&, int64_t, int);
void serialize_int(outstream_info&, int);
void unserialize_int(instream_info&, int&);
void serialize_ul(outstream_info&, unsigned long);
void unserialize_ul(instream_info&, unsigned long&);

#endif
