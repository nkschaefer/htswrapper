#ifndef _HTSWRAPPER_LOOM_READER_H
#define _HTSWRAPPER_LOOM_READER_H
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

// For reading h5 loom files (for RNA velocity data)

namespace sch5{
    
    class loom : public h5_reader{
        
        private:

            std::string parse_loom_bc(const std::string& bc);

        public:
            
            loom(const std::string& fn);
            ~loom();
            
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
            
            void load_expr(std::map<int32_t, std::map<int32_t, double> >& mapmtx);

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
            
            void load_layer(const std::string& layername,
                std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr);
        
            void load_layer(const std::string& layername,
                std::map<int32_t, std::map<int32_t, double> >& mtx);
            
            void load_spliced(std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr);

            void load_unspliced(std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr);

            void load_spliced(std::map<int32_t, std::map<int32_t, double> >& mtx);
            
            void load_unspliced(std::map<int32_t, std::map<int32_t, double> >& mtx);

            void write_layer(const std::string& layername,
                std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr);
            
            void write_layer(const std::string& layername,
                std::map<int32_t, std::map<int32_t, double> >& mapmtx);
    };
};
#endif
