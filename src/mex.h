#ifndef _HTSWRAPPER_MEX_H
#define _HTSWRAPPER_MEX_H
#include <cstdlib>
#include <string.h>
#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <zlib.h>
#include <map>
#include "robin_hood/robin_hood.h"
#include <vector>

// Load a market exchange format (MEX) file
bool parse_mex(const std::string& barcodesfile,
    const std::string& featuresfile,
    const std::string& matrixfile, 
    robin_hood::unordered_map<unsigned long, std::map<int, long int> >& counts,
    std::vector<std::string>& labels,
    const std::string& featuretype = "");

// Write a sparse matrix to disk in MEX format
void write_mex(std::string& out_dir,
    robin_hood::unordered_map<unsigned long, std::map<int, double> >& mtx,
    std::vector<std::string>& features,
    bool round_counts = false,
    std::string barcode_group = "",
    bool cellranger = false,
    bool seurat = false,
    bool underscore = false);

#endif
