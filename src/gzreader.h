#ifndef _HTSWRAPPER_GZREADER_H
#define _HTSWRAPPER_GZREADER_H
#include <cstdlib>
#include <string.h>
#include <utility>
#include <bitset>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <zlib.h>

class gzreader{
    private:
        bool is_gzipped;
        int bufsize;
        int strbufsize;
        char* buf;
        gzFile inf_gz;
        FILE* inf;
        bool eof;
        int idx_start;
        int nread;
        int line_start;
        bool init_read_empty;
        bool split;
        char token;
        std::string field;
        void split_fields();
    public:
        gzreader(std::string filename);
        ~gzreader();
        bool next();
        char* line;
        std::vector<std::string> fields;
        void delimited();
        void delimited(bool d);
        void delimited(char t);
};

#endif
