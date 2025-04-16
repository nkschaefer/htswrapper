#include <cstdlib>
#include <string.h>
#include <utility>
#include <string>
#include <iostream>
#include <sstream>
#include <zlib.h>
#include <vector>
#include "gzreader.h"

// Class to read through gzipped (or not gzipped) files

using namespace std;

/**
 * Constructor
 */
gzreader::gzreader(string fn){

    this->bufsize = 1048576;
    this->strbufsize = 1024;
    this->buf = (char*)malloc((bufsize)*sizeof(char));    
    this->line = (char*)malloc((strbufsize)*sizeof(char));
    this->split = false;
    this->token = '\t';

    // Check if file is gzipped.
    FILE* ftest = fopen(fn.c_str(), "r");
    if (ftest == NULL){
        fprintf(stderr, "ERROR: could not open file %s\n", fn.c_str());
        exit(1);
    }
    size_t testbytes = fread(&this->buf[0], sizeof(char), 2, ftest);
    if (testbytes < 2){
        fprintf(stderr, "ERROR: could not read file %s\n", fn.c_str());
        exit(1);
    }
    fclose(ftest);
    
    is_gzipped = false;
    if ((unsigned char)buf[0] == 0x1f && (unsigned char)buf[1] == 0x8b){
        is_gzipped = true;
    }
    
    this->eof = false;
    this->idx_start = 0;
    this->nread = 0;
    this->line_start = 0;
    this->init_read_empty = false;
    
    if (is_gzipped){
        this->inf = NULL;
        this->inf_gz = gzopen(fn.c_str(), "r");
        // Do initial read.
        nread = gzread(this->inf_gz, &(this->buf[idx_start]), this->bufsize-idx_start);
        eof = gzeof(this->inf_gz);
    }
    else{
        this->inf = fopen(fn.c_str(), "r");
        // Do initial read.
        nread = fread(&(this->buf[idx_start]), sizeof(char), this->bufsize-idx_start, this->inf);
        //nread = fread(this->inf, &(this->buf[idx_start]), sizeof(char), this->bufsize-idx_start);
        eof = feof(this->inf);
    }
    
    if (nread < bufsize-idx_start){
        eof = true;
    }
    init_read_empty = nread == 0;
    line_start = 0;
}

/**
 * Destructor
 */
gzreader::~gzreader(){
    if (this->is_gzipped){
        gzclose(inf_gz);
    }
    else{
        fclose(inf);
    }
    free(this->buf);
    free(this->line);
}

/**
 * Populates this->line, a char* contanining the current line
 * returns true/false == has next line
 */
bool gzreader::next(){
    if (init_read_empty){
        return false;
    }
    if (eof && line_start >= nread){
        return false;
    }

    while(true){
        
        // Parse lines.
        for (int i = line_start; i < nread + idx_start; ++i){
            if (buf[i] == '\n'){
                // Ensure string buffer is big enough to hold the line that is
                // going to be copied into it
                if (i-line_start > this->strbufsize){
                    this->strbufsize = i-line_start+1;
                    this->line = (char*)realloc(line, this->strbufsize);
                }
                strncpy(&line[0], &buf[line_start], i-line_start);
                line[i-line_start] = '\0';
                line_start = i + 1;
                //return !(eof && line_start >= nread);
                split_fields();
                return true;
            }
        }
        if (eof && line_start < nread + idx_start){
            // Get last bit
            strncpy(&line[0], &buf[line_start], nread+idx_start-line_start);
            line[nread+idx_start-line_start] = '\0';
            line_start = nread;
            split_fields();
            return true;
        }

        if (line_start < bufsize){
            // Need to copy what remains in buffer to beginning.
            memmove(&buf[0], &buf[line_start], bufsize-line_start);
            idx_start = bufsize-line_start;
        }
        else{
            idx_start = 0;
        }

        if (!eof){
            // Read the next chunk from the gzFile
            if (this->is_gzipped){
                nread = gzread(this->inf_gz, &(this->buf[idx_start]), this->bufsize-idx_start);
                eof = gzeof(this->inf_gz);
            }
            else{
                nread = fread(&this->buf[idx_start], sizeof(char), this->bufsize-idx_start, this->inf);
                //nread = fread(this->inf, &(this->buf[idx_start]), sizeof(char), this->bufsize-idx_start);
                eof = feof(this->inf);
            }
            if (nread < bufsize-idx_start){
                eof = true;
            }
            line_start = 0;
        }
        else{
            // Needed to read and can't read anything else.
            return false;
        }
       
    }
    return false;
}

void gzreader::delimited(){
    this->split = true;
}

void gzreader::delimited(bool delim){
    this->split = delim;
}

void gzreader::delimited(char tok){
    this->split = true;
    this->token = tok;
}

void gzreader::split_fields(){
    if (split){
        fields.clear();
        istringstream splitter(line);
        while (getline(splitter, field, token)){
            fields.push_back(field);
        } 
    }
}

