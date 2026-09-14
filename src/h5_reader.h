#ifndef _HTSWRAPPER_H5_READER_H
#define _HTSWRAPPER_H5_READER_H
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
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <cstdlib>
#include <utility>
#include <highfive/H5File.hpp>

// Parent class (to be extended) for reading/writing single-cell data sets
// to/from hdf5 (h5) format.

namespace sch5{
    
    const short h5_type_unknown = -1;
    const short h5_type_str = 0;
    const short h5_type_int = 1;
    const short h5_type_float = 2;

    class h5_reader{
        
        private:
            
            template <typename T>
            struct int_needs_range_check: std::integral_constant<bool,
                    (sizeof(T) > sizeof(int)) ||
                    (sizeof(T) == sizeof(int) && std::is_unsigned<T>::value)> {};

            // Range-checked narrowing (selected for wide / unsigned-int-width types).
            template <typename T>
            int narrow_to_int(T v, std::true_type) {
                // long long holds every value of any T reaching this overload EXCEPT
                // the top of uint64_t's range — handled separately below.
                if (static_cast<long long>(v) < static_cast<long long>(std::numeric_limits<int>::min()) ||
                    static_cast<long long>(v) > static_cast<long long>(std::numeric_limits<int>::max())) {
                    throw std::runtime_error("column has a value outside int range");
                }
                return static_cast<int>(v);
            }

            // No-check narrowing (selected for int8/int16/uint8 — always fit in int).
            template <typename T>
            int narrow_to_int(T v, std::false_type) {
                return static_cast<int>(v);
            }
            
            inline int narrow_u64_to_int(uint64_t v){
                // Unsigned: no lower bound to check (can't be < INT_MIN).
                // Compare in uint64_t space so there's no signed/unsigned pitfall
                // and no cast that could overflow.
                if (v > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
                    throw std::runtime_error("column has a value outside int range");
                }
                return static_cast<int>(v);
            }

            template <typename T>
            void read_int_flex(const HighFive::DataSet& ds, std::vector<int>& out) {
                static_assert(std::is_integral<T>::value, "read_int_as requires an integral T");
                std::vector<T> buf;
                ds.read(buf);
                out.clear();
                out.reserve(buf.size());
                for (size_t i = 0; i < buf.size(); ++i) {
                    out.push_back(narrow_to_int<T>(buf[i], int_needs_range_check<T>{}));
                }
            }

            template <typename T>
            void read_float_flex(const HighFive::DataSet& ds,
                               std::vector<double>& out,
                               bool fix_nan) {
                static_assert(std::is_floating_point<T>::value, "read_float_flex requires a floating-point T");
                std::vector<T> buf;
                ds.read(buf);
                out.clear();
                out.reserve(buf.size());
                for (size_t i = 0; i < buf.size(); ++i) {
                    double v = static_cast<double>(buf[i]);   // float→double is exact
                    if (fix_nan && !std::isfinite(v)) {
                        v = 0.0;
                    }
                    out.push_back(v);
                }
            }

        protected:
            
            // Name of input h5 file
            std::string filename;
            
            // Underlying HighFive file parser
            HighFive::File file;

            // Name of main count matrix
            std::string countsname;
            
            // Cell barcodes to skip
            std::unordered_set<std::string> excl_bc;
            
            // Look up encoding of h5 group
            bool check_encoding(const HighFive::Group& g);
            
            void load_int_col_aux(HighFive::DataSet& ds,
                std::vector<int>& col_data);

            void load_float_col_aux(HighFive::DataSet& ds,
                std::vector<double>& col_data,
                bool fix_nan = true);
            
            void load_str_col_aux(HighFive::DataSet& ds,
                std::vector<std::string>& col_data);
            
            short type_from_ds(HighFive::DataSet& ds);
        
            // Load a count matrix (in CSR format)
            void load_mtx(const std::string& name,
                std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr);
            
            // Load a dense matrix into CSR format
            void load_mtx_dense(const std::string& name,
                std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr);

            // Load a dense matrix into a map
            void load_mtx_dense(const std::string& name,
                std::map<int32_t, std::map<int32_t, double> >& mtxmap);
            
            bool write_dense_aux(const std::string& path,
                std::vector<double>& data,
                bool rows_are_genes,
                bool force=false);

            bool write_mtx_dense(const std::string& path,
                std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr,
                bool rows_are_genes,
                bool force=false);

            bool write_mtx_dense(const std::string& path,
                std::map<int32_t, std::map<int32_t, double> >& mtxmap,
                bool rows_are_genes,
                bool force=false);

        public:
            
            // How many cells?
            long int n_cells;

            // How many genes?
            long int n_genes;
            
            // Cell barcodes
            std::vector<std::string> cell_names;

            // Gene names
            std::vector<std::string> gene_names;
            
            h5_reader(const std::string& fn);
            virtual ~h5_reader();
            
            void open(const std::string& fn);

            // Tell it the name of the (raw, un-scaled) counts layer
            void set_countsname(const std::string& n);
            
            virtual void load_meta_col(const std::string& col_name, 
                std::vector<int>& col_vals,
                bool fix_nan = true) = 0;
            
            virtual void load_meta_col(const std::string& col_name, 
                std::vector<double>& col_vals,
                bool fix_nan = true) = 0;

            virtual void load_meta_col(const std::string& col_name, 
                std::vector<std::string>& col_vals,
                bool fix_nan = true) = 0;
            
            // Convert a CSC-format sparse matrix to CSR format
            void csc_to_csr(std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr);
            
            // Load expression data in CSR format
            virtual void load_expr(std::vector<double>& X_data,
                std::vector<int32_t>& X_indices,
                std::vector<int64_t>& X_indptr) = 0;
           
            // Retrieve the names of metadata columns of each type 
            virtual void get_meta_colnames(std::vector<std::string>& str_cols,
                std::vector<std::string>& int_cols,
                std::vector<std::string>& float_cols) = 0;
            
            // Check for the existence of a metadata column
            virtual bool has_meta_col(const std::string& name) = 0;
            
            // Return the type of a metadata column
            virtual short meta_col_type(const std::string& name) = 0;
            
            virtual void list_layers(std::vector<std::string>& layers) = 0;

            virtual void load_layer(const std::string& name,
                std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr) = 0;

            void load_layer(const std::string& name,
                std::map<int32_t, std::map<int32_t, double> >& mtx);

            void mtx2map(std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr,
                std::map<int32_t, std::map<int32_t, double> >& mtxmap);
            
            void map2mtx(std::map<int32_t, std::map<int32_t, double> >& mtxmap,
                std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr);
            
            // Write a new layer (same dims as main data matrix) to the file
            virtual void write_layer(const std::string& layername,
                std::vector<double>& data,
                std::vector<int32_t>& indices,
                std::vector<int64_t>& indptr) = 0;
            
            void write_layer(const std::string& layername,
                std::map<int32_t, std::map<int32_t, double> >& mtxmap);

    };
};
#endif
