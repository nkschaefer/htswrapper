#include <cstdlib>
#include <string.h>
#include <utility>
#include <string>
#include <iostream>
#include <zlib.h>
#include <cmath>
#include <sstream>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include "mex.h"
#include "gzreader.h"
#include "bc.h"

// Functions for reading and writing Market exchange (MEX) format
// representations of sparse matrices

using namespace std;

/**
 * Parses input files in Market exchange format (MEX),
 * used to represent sparse matrices. Takes as input the three MEX files
 * (gzipped or not), along with two data structures to store results.
 * counts will store (hashed) cell barcodes mapped to counts of label indices,
 * and labels will store the names of the labels (i.e. genes).
 */
bool parse_mex(const string& barcodesfile,
    const string& featuresfile,
    const string& matrixfile, 
    robin_hood::unordered_map<unsigned long, map<int, long int> >& counts,
    vector<string>& labels,
    const string& featuretype){
     
    // First, parse barcodes
    vector<unsigned long> barcodes;
    parse_barcode_file(barcodesfile, barcodes);
    
    set<string> unique_featuretype; 
    
    // Map feature index in whole file to feature index in subset 
    map<int, int> feature_inds;
    
    // Then, parse genes
    // Feature index in file
    int feature_idx = 0;
    // Feature index in subset
    int feature_idx_subset = 0;

    gzreader labelsreader(featuresfile);
    int nfeatures_match = 0;
    while(labelsreader.next()){
        istringstream splitter(labelsreader.line);
        string field;
        int idx = 0;
        string id = "";
        string name = "";
        string type = "";
        while(getline(splitter, field, '\t')){
            if (idx == 0){
                id = field;
            }
            else if (idx == 1){
                name = field;
            }
            else if (idx == 2){
                type = field;
            }
            idx++;
        }
        if (featuretype == "" || type == "" || featuretype == type){
            if (name != ""){
                labels.push_back(name);
            }
            else{
                labels.push_back(id);
            }
            feature_inds.insert(make_pair(feature_idx, feature_idx_subset));
            ++feature_idx_subset;
            ++nfeatures_match;
        }
        if (featuretype == "" && type != ""){
            unique_featuretype.insert(type);
        }
        ++feature_idx;
    }
    
    if (nfeatures_match == 0){
        fprintf(stderr, "ERROR: 0 features match provided type %s\n", featuretype.c_str());
        fprintf(stderr, "Allowed feature types:\n");
        for (set<string>::iterator ut = unique_featuretype.begin(); ut != unique_featuretype.end();
            ++ut){
            fprintf(stderr, "%s\n", ut->c_str());
        }
        exit(1);
    }
    
    set<int> bciuniq;

    bool barcodes_first = false;
    
    // Then, populate data structure
    gzreader mtxreader(matrixfile);
    int mexline = 0;
    while(mtxreader.next()){
        if (mtxreader.line[0] == '%'){
            continue;
        }
        else{
            if (mexline == 0){
                // Read header.
                istringstream splitter(mtxreader.line);
                long int n1;
                long int n2;
                int idx = 0;
                string token;
                while(getline(splitter, token, ' ' )){
                    if (idx == 0){
                        n1 = atol(token.c_str());
                    }
                    else if (idx == 1){
                        n2 = atol(token.c_str());
                    }
                    ++idx;
                }
                if (n1 == barcodes.size()){
                    barcodes_first = true;
                }
                else if (n2 == barcodes.size()){
                    barcodes_first = false;
                }
                else{
                    fprintf(stderr, "ERROR: %ld barcodes; does not match MTX file: %ld or %ld\n",
                        barcodes.size(), n1, n2);
                    exit(1);
                }
            }
            else if (mexline > 0){
                istringstream splitter(mtxreader.line);
                int idx = 0;
                string token;
                int bc_idx;
                int feature_idx;
                long int count;
                while (getline(splitter, token, ' ')){
                    if ((!barcodes_first && idx == 0) || (barcodes_first && idx == 1)){
                        // feature index
                        feature_idx = atoi(token.c_str());
                        feature_idx--; // Make 0-based
                    }
                    else if ((!barcodes_first && idx == 1) || (barcodes_first && idx == 0)){
                        // barcode index
                        bc_idx = atoi(token.c_str());
                        bc_idx--; // Make 0-based
                    }
                    else if (idx == 2){
                        // UMI count
                        double count_d = atof(token.c_str());
                        count = (long int)round(count_d);
                        
                        if (feature_inds.count(feature_idx) > 0){ 
                            // This feature is in the list
                            bciuniq.insert(bc_idx);
                            unsigned long barcode = barcodes[bc_idx];
                            if (counts.count(barcode) == 0){
                                map<int, long int> m;
                                counts.emplace(barcode, m);
                            }
                            counts[barcode].emplace(feature_inds[feature_idx], count);
                        }
                    }
                    ++idx;
                }
            }
            ++mexline;
        }
    }
    if (featuretype == "" && unique_featuretype.size() > 1){
        fprintf(stderr, "ERROR: no feature type filter provided, but %ld feature types encountered in \
MEX data\n", unique_featuretype.size());
        return false;
    }
    fprintf(stderr, "Loaded %ld barcodes and %d features\n", barcodes.size(), nfeatures_match);
    return true;
}

void write_mex(string& out_dir,
    robin_hood::unordered_map<unsigned long, map<int, double> >& mtx,
    vector<string>& features, 
    bool round_counts,
    string barcode_group, 
    bool cellranger, 
    bool seurat, 
    bool underscore){
    
    if (out_dir[out_dir.length()-1] == '/'){
        out_dir = out_dir.substr(0, out_dir.length()-1);
    }

    if (!mkdir(out_dir.c_str(), 0775)){
        // Assume directory already exists
    }

    string out_barcodes = out_dir + "/barcodes.tsv.gz";
    string out_features = out_dir + "/features.tsv.gz";
    string out_mtx = out_dir + "/matrix.mtx.gz";

    // Write features
    gzFile out_features_f = gzopen(out_features.c_str(), "w");
    for (int i = 0; i < features.size(); ++i){
        gzwrite(out_features_f, features[i].c_str(), features[i].length());
        gzwrite(out_features_f, "\n", 1); 
    }
    gzclose(out_features_f);
    
    gzFile out_barcodes_f = gzopen(out_barcodes.c_str(), "w");
    
    // Indices must be 1-based
    int bc_idx = 1;

    long int n_entries = 0;

    vector<string> mtxlines;
    
    // This is wasteful but it's also 2024 and people have RAM
    char mtxbuf[1024];

    for (robin_hood::unordered_map<unsigned long, map<int, double> >::iterator x = mtx.begin();
        x != mtx.end(); ++x){    
        
        string bc_str = bc2str(x->first);
        mod_bc_libname(bc_str, barcode_group, cellranger, seurat, underscore);
        
        gzwrite(out_barcodes_f, bc_str.c_str(), bc_str.length());
        gzwrite(out_barcodes_f, "\n", 1);
        
        // CellRanger format: features first, then barcodes, then entry
        for (map<int, double>::iterator y = x->second.begin(); y != x->second.end(); ++y){
            if (round_counts){
                sprintf(&mtxbuf[0], "%d %d %d", y->first+1, bc_idx, (int)round(y->second));
            }
            else{
                sprintf(&mtxbuf[0], "%d %d %.4f", y->first+1, bc_idx, y->second);
            }
            mtxlines.push_back(mtxbuf);
            n_entries++;
        }
        ++bc_idx;
    }
    
    gzclose(out_barcodes_f);

    gzFile out_mtx_f = gzopen(out_mtx.c_str(), "w");
    string hdrline = "%%MatrixMarket matrix coordinate real general\n";
    if (round_counts){
        hdrline = "%%MatrixMarket matrix coordinate integer general\n";
    }
    gzwrite(out_mtx_f, hdrline.c_str(), hdrline.length());
    sprintf(&mtxbuf[0], "%ld %d %ld\n", features.size(), bc_idx-1, mtxlines.size());
    string mtxline = mtxbuf;
    gzwrite(out_mtx_f, mtxline.c_str(), mtxline.length());
    for (int i = 0; i < mtxlines.size(); ++i){
        gzwrite(out_mtx_f, mtxlines[i].c_str(), mtxlines[i].length());
        gzwrite(out_mtx_f, "\n", 1);
    } 
    gzclose(out_mtx_f);

}

