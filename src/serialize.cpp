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
#include <climits>
#include "serialize.h"

using std::cout;
using std::endl;
using namespace std;

int int_bytes = sizeof(int);
int longint_bytes = sizeof(long int);
int int64_bytes = sizeof(int64_t);
int float_bytes = sizeof(float);

/**
 * Initialize file reading stuff
 */
void instream_init(instream_info& is, string filename, int bufsize){
    if (bufsize % CHAR_BIT != 0){
        // Make bufsize a round number of bytes.
        bufsize += (bufsize % CHAR_BIT);
    }
    
    // Open file for reading.
    
    is.input = fopen(filename.c_str(), "rb");
    is.bufsize = bufsize;
    is.bytes = (char*)malloc(bufsize+1);
    is.buffer_index = 0;
    is.eof = false;
    is.initialized = true;
    
    is.prevoffset = 0;
    is.prevoffset_rec = 0;
    is.prev_buf_idx = 0;
    is.prev_start_idx = 0;
    
    // Do a first read.
    is.prevoffset = ftello(is.input);
    int bytes_read = fread(&is.bytes[0], sizeof(char), is.bufsize, is.input);
    if (bytes_read < bufsize){
        is.eof = true;
    }
    is.buf_stored = bytes_read;
}

void outstream_init(outstream_info& os, string filename){
    os.output.open(filename, std::ofstream::binary);
    os.initialized = true;
}
void outstream_close(outstream_info& os){
    os.output.close();
    os.initialized = false;
}

/**
 * Needed by many other functions. If we've reached too close to the end of the buffer
 * to take in the data we need, we have to copy whatever's left at the end of the buffer
 * to the beginning and fill the rest of the buffer from the input stream. If there's a
 * problem doing this, we'll have to bail out.
 */
void check_bytes_left(int bytes_needed, instream_info& is){
    
    if ((bytes_needed + is.buffer_index) > is.buf_stored){
        if (feof(is.input)){
            // We've hit the end of the input stream and the buffer is not full.
            // Unfortunately, asking for the amount of data we need will throw us off the
            // end of not-garbage memory. We need to bail out here.
            is.eof = true;
            is.buf_stored = -1;
            return;
        }
        
        is.prev_start_idx = 0;
        char* read_ptr = &is.bytes[0];
        int readsize = is.bufsize;

        // Copy whatever is left of what we need to the beginning of the input buffer.
        
        if (is.buffer_index < (is.buf_stored)){
            int bytestocopy = is.buf_stored-is.buffer_index;
            memmove(&is.bytes[0], &is.bytes[is.buffer_index], bytestocopy);
            // Get ready for the next file read
            read_ptr += bytestocopy;
            is.prev_start_idx += bytestocopy;
            readsize -= bytestocopy;
            is.buf_stored = bytestocopy;
        }
        else{
            is.buf_stored = 0;
        }
        // No matter what we just did, we want to look at the first byte in the array
        // first
        is.buffer_index = 0;
        
        // Read in more from the input buffer.
        
        is.prevoffset = ftello(is.input);

        int bytes_read = fread(read_ptr, sizeof(char), readsize, is.input);
        
        if (bytes_read < readsize){
            is.eof = true;
            is.last_bytes_read = bytes_read;
        }
        is.buf_stored += bytes_read;
    }    
}

void seek (instream_info& is, int64_t pos, int bufind){ 
    is.buf_stored = 0; 
    is.buffer_index = 0; 
    fflush(is.input); 
    fseeko(is.input, pos, SEEK_SET);  
    int bytes_read = fread(is.bytes, sizeof(char), is.bufsize, is.input); 
    if (bytes_read < is.bufsize){ 
        is.eof = true; 
        is.last_bytes_read = bytes_read; 
    } 
    is.buf_stored += bytes_read; 
    is.buffer_index = bufind; 
}

/**
 * Function to write out an integer in binary format.
 */
void serialize_int(outstream_info& os, int num){
    const char* buf = reinterpret_cast<const char *>(&num);
    
    char bytes[int_bytes + 1];
    memset(&bytes, 0x00, int_bytes);
    bytes[int_bytes] = '\0';
    int byte_index = 0;
    
    for (int i = 0; i < int_bytes; i++){
        bytes[byte_index] = buf[i];
        byte_index++;
    }
    
    os.output.write(&bytes[0], int_bytes);
}

/**
 * Function to read in an integer in binary format.
 */
void unserialize_int(instream_info& is, int& num){
    check_bytes_left(int_bytes, is);
    
    // Make sure there's not extra garbage there in case the system's representation
    // of int is bigger than ours
    memset(&num, 0x00, sizeof(int));
    
    memcpy(&num, &is.bytes[is.buffer_index], int_bytes);
    is.buffer_index += int_bytes;
    
}

/**
 * Function to write out an unsigned long integer in binary format.
 */
void serialize_ul(outstream_info& os, unsigned long num){
    const char* buf = reinterpret_cast<const char *>(&num);
    
    char bytes[sizeof(unsigned long) + 1];
    memset(&bytes, 0x00, sizeof(unsigned long));
    bytes[sizeof(unsigned long)] = '\0';
    int byte_index = 0;
    
    for (int i = 0; i < sizeof(unsigned long); i++){
        bytes[byte_index] = buf[i];
        byte_index++;
    }
    
    os.output.write(&bytes[0], sizeof(unsigned long));
}

/**
 * Function to read in an unsigned long integer in binary format.
 */
void unserialize_ul(instream_info& is, unsigned long& num){
   
    check_bytes_left(sizeof(unsigned long), is);
    
    // Make sure there's not extra garbage there in case the system's representation
    // of int is bigger than ours
    memset(&num, 0x00, sizeof(unsigned long));
    
    memcpy(&num, &is.bytes[is.buffer_index], sizeof(unsigned long));
    is.buffer_index += sizeof(unsigned long);
    
}

