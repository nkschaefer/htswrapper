#ifndef _HTSWRAPPER_ANNDATA_H
#define _HTSWRAPPER_ANNDATA_H
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <cstdlib>
#include <utility>
#include <highfive/H5File.hpp>
#include "h5_reader.h"

// Class for reading/manipulating h5ad (h5 files for AnnData) objects.

namespace sch5{
    
    class anndata : public h5_reader{

        public:
            
            anndata(const std::string& fn);
            ~anndata();
            
            void load_meta_col(const std::string& col_name, 
                std::vector<int>& col_vals,
                bool fix_nan = true) override;
            
            void load_meta_col(const std::string& col_name, 
                std::vector<double>& col_vals,
                bool fix_nan = true) override;

            void load_meta_col(const std::string& col_name, 
                std::vector<std::string>& col_vals,
                bool fix_nan = true) override;
            
            // Load expression data in CSR format
            void load_expr(std::vector<double>& X_data,
                std::vector<int32_t>& X_indices,
                std::vector<int64_t>& X_indptr) override;
           
            // Retrieve the names of metadata columns of each type 
            void get_meta_colnames(std::vector<std::string>& str_cols,
                std::vector<std::string>& int_cols,
                std::vector<std::string>& float_cols) override;
            
            // Check for the existence of a metadata column
            bool has_meta_col(const std::string& name) override;
            
            // Return the type of a metadata column
            short meta_col_type(const std::string& name) override;
            
            // Get a list of layer names
            void list_layers(std::vector<std::string>& layernames) override;
            
            bool has_layer(const std::string& layername);

            void load_layer(const std::string& layername,
                std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr);

            void write_layer(const std::string& layername,
                std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr);

    };
};
#endif
