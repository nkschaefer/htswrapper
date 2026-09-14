#ifndef _HTSWRAPPER_SEURAT_H
#define _HTSWRAPPER_SEURAT_H
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

// Class for reading/manipulating h5Seurat (h5 files from SeuratDisk) objects.

namespace sch5{

    class seurat : public h5_reader{

        private:

            std::string active_assay;

            // Data slot names within an assay that are not metadata
            static bool is_data_slot(const std::string& name);

            // Load a sparse or dense matrix, handling h5Seurat's convention
            // where encoding-type may be absent (dgCMatrix from SeuratDisk).
            void load_seurat_mtx(const std::string& name,
                std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr);

        public:

            seurat(const std::string& fn);
            ~seurat();

            void set_active_assay(const std::string& name);

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

            // Get a list of layer names (data slots in the active assay)
            void list_layers(std::vector<std::string>& layernames) override;

            bool has_layer(const std::string& layername);

            void load_layer(const std::string& layername,
                std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr) override;

            void write_layer(const std::string& layername,
                std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr) override;

    };
};
#endif
