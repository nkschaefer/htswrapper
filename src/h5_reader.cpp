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
#include <set>
#include <unordered_set>
#include <cstdlib>
#include <utility>
#include <highfive/H5File.hpp>
#include "h5_reader.h"

using std::cout;
using std::endl;
using namespace std;

namespace sch5{

    /**
     * Constructor
     */
    h5_reader::h5_reader(const string& fn) try 
        : file(fn, HighFive::File::ReadOnly){
        n_cells = 0;
        n_genes = 0;
        filename = fn;
    } catch(const HighFive::Exception& e){
        throw runtime_error("Error opening h5 file " + fn); 
    }
    
    /**
     * Destructor
     */
    h5_reader::~h5_reader(){

    }

    /**
     * Tell loader where to look for the expression data
     */
    void h5_reader::set_countsname(const string& name){
        this->countsname = name;
    }
    
    /**
     * Load any int data type column as standard int
     */ 
    void h5_reader::load_int_col_aux(HighFive::DataSet& dataset,
        vector<int>& col_vals){
        
        auto dt = dataset.getDataType();
        if (dt == HighFive::AtomicType<int8_t>()){
            read_int_flex<int8_t>(dataset, col_vals);
        }
        else if (dt == HighFive::AtomicType<int16_t>()){
            read_int_flex<int16_t>(dataset, col_vals);
        }
        else if (dt == HighFive::AtomicType<int32_t>()){
            read_int_flex<int32_t>(dataset, col_vals);
        }
        else if (dt == HighFive::AtomicType<int64_t>()){
            read_int_flex<int64_t>(dataset, col_vals);
        }
        else if (dt == HighFive::AtomicType<uint8_t>()){
            read_int_flex<uint8_t>(dataset, col_vals);
        }
        else if (dt == HighFive::AtomicType<uint16_t>()){
            read_int_flex<uint16_t>(dataset, col_vals);
        }
        else if (dt == HighFive::AtomicType<uint32_t>()){
            read_int_flex<uint32_t>(dataset, col_vals);
        }
        else if (dt == HighFive::AtomicType<uint64_t>()){
            std::vector<uint64_t> buf;
            dataset.read(buf);
            col_vals.clear();
            col_vals.reserve(buf.size());
            for (size_t i = 0; i < buf.size(); ++i){
                col_vals.push_back(narrow_u64_to_int(buf[i]));
            }
            //read_int_flex<uint64_t>(ds, col_vals);
        }
        else{
            throw runtime_error("Encountered unsupported int column type");
        }
    }
    
    /**
     * Load any int data type column as standard int
     */ 
    void h5_reader::load_float_col_aux(HighFive::DataSet& dataset,
        vector<double>& col_vals,
        bool fix_nan){
        
        auto dt = dataset.getDataType();
        if (dt == HighFive::AtomicType<float>()){
            read_float_flex<float>(dataset, col_vals, fix_nan);
        }
        else if (dt == HighFive::AtomicType<double>()){
            read_float_flex<double>(dataset, col_vals, fix_nan);
        }
        else if (dt == HighFive::AtomicType<long double>()){
            read_float_flex<long double>(dataset, col_vals, fix_nan);
        }
        else{
            throw runtime_error("encountered unsupported float column type");
        }
    }

    /**
     * Load fixed & variable-length strings as standard strings
     */
    void h5_reader::load_str_col_aux(HighFive::DataSet& dataset,
        vector<string>& col_vals){
        
        auto dtype = dataset.getDataType();
        if (dtype.getClass() != HighFive::DataTypeClass::String){
            throw runtime_error("encountered unsupported string column type");
        }
        const size_t n = dataset.getSpace().getElementCount();
        if (n == 0){
            return;
        }
        col_vals.clear();
        if (dtype.isVariableStr()){
            dataset.read(col_vals);
        }
        else{
            // Fixed-length strings
            const size_t width = dtype.getSize();
            vector<char> buf(n * width);
            dataset.read_raw(buf.data(), dtype);
            col_vals.reserve(n);
            for (size_t i = 0; i < n; ++i){
                const char* p = buf.data() + i * width;
                // Length = up to first null
                size_t len = 0;
                while (len < width && p[len] != '\0'){
                    ++len;
                }
                while (len > 0 && p[len - 1] == ' '){
                    --len;
                }
                col_vals.emplace_back(p, len);
            }
        }
    }
    
    short h5_reader::type_from_ds(HighFive::DataSet& ds){
        auto dtype = ds.getDataType();
        auto type_class = dtype.getClass(); 
        if (type_class == HighFive::DataTypeClass::String){
            return h5_type_str;
        }
        else if (type_class == HighFive::DataTypeClass::Integer){
            return h5_type_int;
        }
        else if (type_class == HighFive::DataTypeClass::Float){
            return h5_type_float;
        }
        return h5_type_unknown;

        /*
        if (dtype == HighFive::AtomicType<int8_t>()){
            return h5_type_int;    
        } 
        else if (dtype == HighFive::AtomicType<uint8_t>()){
            return h5_type_int;
        }
        else if (dtype == HighFive::AtomicType<int16_t>()){
            return h5_type_int;    
        } 
        else if (dtype == HighFive::AtomicType<uint16_t>()){
            return h5_type_int;
        }
        else if (dtype == HighFive::AtomicType<int32_t>()){
            return h5_type_int;
        } 
        else if (dtype == HighFive::AtomicType<uint32_t>()){
            return h5_type_int;
        }
        else if (dtype == HighFive::AtomicType<int64_t>()){
            return h5_type_int; 
        }
        else if (dtype == HighFive::AtomicType<uint64_t>()){
            return h5_type_int;
        }
        else if (dtype == HighFive::AtomicType<float>()){
            return h5_type_float;
        } 
        else if (dtype == HighFive::AtomicType<double>()){
            return h5_type_float;
        } 
        else if (dtype == HighFive::AtomicType<long double>()){
            return h5_type_float;
        }
        else if (dtype.getClass() == HighFive::DataTypeClass::String){
            return h5_type_str;
        }
        else{
            return h5_type_unknown;
        }
        */
    }

    /**
     * true: csr
     * false: csc
     * error: neither
     */
    bool h5_reader::check_encoding(const HighFive::Group& g){
        string encoding;
        g.getAttribute("encoding-type").read(encoding);
        if (encoding == "csr_matrix"){
            return true;
        }
        else if (encoding == "csc_matrix"){
            return false;
        }
        else{
            string name = g.getPath();
            throw runtime_error("layer " + name + " has unexpected encoding type: " + encoding); 
        }
        return false;
    }

    /**
     * Convert sparse data in CSC format to CSR format
     */
    void h5_reader::csc_to_csr(vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr){

        vector<double> data_csr;
        vector<int32_t> indices_csr;
        vector<int64_t> indptr_csr;
        
        const int64_t num = indptr.back();
        indptr_csr.assign(n_cells + 1, 0);
        for (int64_t k = 0; k < num; ++k){
            indptr_csr[indices[k] + 1]++;
        }
        for (long int i = 0; i < n_cells; ++i){
            indptr_csr[i + 1] += indptr_csr[i];
        }
        
        data_csr.resize(num);
        indices_csr.resize(num);
        vector<int64_t> next(indptr_csr.begin(), indptr_csr.end() - 1);
        
        for (int32_t col = 0; col < (int32_t)n_genes; ++col){
            for (int64_t jj = indptr[col]; jj < indptr[col + 1]; ++jj) {
                int32_t row = indices[jj];
                int64_t dest = next[row]++;
                indices_csr[dest] = col;
                data_csr[dest] = data[jj];
            }
        }
        
        data.clear();
        indices.clear();
        indptr.clear();

        data = data_csr;
        indices = indices_csr;
        indptr = indptr_csr;
    }

    /**
     * Loads count data from /layers, ensures data are in CSR format
     */
    void h5_reader::load_mtx(const string& name,
        vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr){
        
        auto obj_type = file.getObjectType(name);

        if (obj_type == HighFive::ObjectType::Group) {
            // sparse
            auto grp = file.getGroup(name);
            grp.getDataSet("data").read(data);
            grp.getDataSet("indices").read(indices);
            grp.getDataSet("indptr").read(indptr);
            if (!check_encoding(grp)){
                // CSC format
                csc_to_csr(data, indices, indptr);
            }
        }
        else if (obj_type == HighFive::ObjectType::Dataset){
            // dense
            vector<double> mtx(n_cells * n_genes);
            auto ds = file.getDataSet(name);
            ds.read_raw<double>(mtx.data());
            int64_t num = 0; 
            for (long int i = 0; i < n_cells; ++i){
                indptr.push_back(num);
                const double* row = mtx.data() + i * n_genes;
                for (long int j = 0; j < n_genes; ++j) {
                    double count = row[j];
                    if (count > 0.0){
                        indices.push_back((int32_t)j);
                        data.push_back(count);
                        ++num;
                    }
                }
            }
            indptr.push_back(num);
        }
        else{
            throw runtime_error("h5 path " + name + " has unknown type; cannot load");
        }
    }
    
    /**
     * Load a matrix in dense format and convert to CSR.
     */
    void h5_reader::load_mtx_dense(const string& name,
        vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr){
        
        HighFive::DataSet dset = file.getDataSet(name);
        vector<size_t> dims = dset.getDimensions();
    
        vector<double> dense_mtx(n_genes * n_cells);
        dset.read_raw<double>(dense_mtx.data());
         
        if (dims.size() != 2){
            throw runtime_error("Dimensions missing from dense matrix " + name);
        }
        // dense data is always row major
        bool col_major = false;
        if (dims[0] == n_genes && dims[1] == n_cells){
            //col_major = true;
        }
        else if (dims[0] == n_cells && dims[1] == n_genes){
            //col_major = false;
        }
        else{
            throw runtime_error("Dense matrix dimensions don't match metadata dimensions");
        }
        data.clear();
        indices.clear();
        indptr.clear();
        indptr.reserve((size_t)n_cells + 1);
        
        int64_t idx_global = 0;
        for (int32_t c = 0; c < n_cells; ++c){
            indptr.push_back(idx_global);
            for (int32_t g = 0; g < n_genes; ++g){
                double count = 0.0;
                if (!col_major){
                    count = dense_mtx[g * n_cells + c];
                }
                else{
                    count = dense_mtx[c * n_genes + g];
                }
                if (count != 0.0){
                    indices.push_back(g);
                    data.push_back(count);
                    ++idx_global;
                }
            }
        }
        indptr.push_back(idx_global);
    }
    
    /**
     * Load a matrix in dense format and convert to a map.
     */
    void h5_reader::load_mtx_dense(const string& name,
        map<int32_t, map<int32_t, double> >& mtxmap){
        
        HighFive::DataSet dset = file.getDataSet(name);
        vector<size_t> dims = dset.getDimensions();
    
        vector<double> dense_mtx(n_genes * n_cells);
        dset.read_raw<double>(dense_mtx.data());
         
        if (dims.size() != 2){
            throw runtime_error("Dimensions missing from dense matrix " + name);
        }
        // dense data is always row-major
        bool col_major = false;
        if (dims[0] == n_genes && dims[1] == n_cells){
            //col_major = true;
        }
        else if (dims[0] == n_cells && dims[1] == n_genes){
            //col_major = false;
        }
        else{
            throw runtime_error("Dense matrix dimensions don't match metadata dimensions");
        }
        for (int32_t c = 0; c < n_cells; ++c){
            for (int32_t g = 0; g < n_genes; ++g){
                double count = 0.0;
                if (!col_major){
                    count = dense_mtx[g * n_cells + c];
                }
                else{
                    count = dense_mtx[c * n_genes + g];
                }
                if (count > 0.0){
                    if (mtxmap.count(c) == 0){
                        map<int32_t, double> m;
                        mtxmap.insert(make_pair(c, m));
                    }
                    mtxmap[c].insert(make_pair(g, count));
                }
            }
        }
    }

    /**
     * Helper method for writing dense matrices.
     */
    bool h5_reader::write_dense_aux(const string& path,
        vector<double>& data,
        bool rows_are_genes,
        bool force){
        
        if (data.size() != n_cells * n_genes){
            throw runtime_error("Dense matrix dimensions do not agree with metadata");
        }
        if (file.exist(path)){
            if (!force){
                return false;
            }
            H5Ldelete(file.getId(), path.c_str(), H5P_DEFAULT);
        }
        
        vector<size_t> dims{(size_t)n_cells, (size_t)n_genes};
        if (rows_are_genes){
            dims[0] = (size_t)n_genes;
            dims[1] = (size_t)n_cells;
        }
        
        auto dset = file.createDataSet<double>(path, HighFive::DataSpace(dims));
        dset.write_raw(data.data());
        dset.createAttribute<std::string>("encoding-type",
            HighFive::DataSpace::From(string("array"))).write(string("array"));
        dset.createAttribute<std::string>("encoding-version",
            HighFive::DataSpace::From(string("0.2.0"))).write(string("0.2.0"));
        return true;
    }

    /**
     * Write data to file in dense format.
     */
    bool h5_reader::write_mtx_dense(const string& path,
        vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr,
        bool rows_are_genes,
        bool force){
        
        // Create dense data
        vector<double> dense(n_cells * n_genes, 0.0);
        for (int i = 0; i < n_cells; ++i){
            int start = indptr[i];
            int end = indptr[i + 1];
            for (int j = start; j < end; ++j){
                int gene_idx = indices[j];
                double count = data[j];
                int ix;
                if (rows_are_genes){
                    ix = gene_idx * n_cells + i;
                }
                else{
                    ix = i * n_genes + gene_idx;
                }
                dense[ix] = count;
            }
        }
        return write_dense_aux(path, dense, rows_are_genes, force);
    }
    
    /**
     * Write data to file in dense format.
     */
    bool h5_reader::write_mtx_dense(const string& path,
        map<int32_t, map<int32_t, double> >& mtx,
        bool rows_are_genes,
        bool force){
        
        // Create dense data
        vector<double> dense(n_cells * n_genes, 0.0);
        
        for (map<int32_t, map<int32_t, double> >::iterator m = mtx.begin(); m != mtx.end(); ++m){
            for (map<int32_t, double>::iterator m2 = m->second.begin(); m2 != m->second.end(); 
                ++m2){
                int gene_idx = m2->first;
                int obs_idx = m->first;
                // dense mtx must always be row major
                int32_t ix;
                if (rows_are_genes){
                    ix = gene_idx * n_cells + obs_idx;
                }
                else{
                    ix = obs_idx * n_genes + gene_idx;
                }
                dense[ix] = m2->second;
            }
        }

        return write_dense_aux(path, dense, rows_are_genes, force);
    }

    /**
     * Convert a map representation of CSC sparse matrix format to
     * vector representation (used by h5 files)
     */
    void h5_reader::map2mtx(map<int32_t, map<int32_t, double> >& mtx,
        vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr){
        
        int32_t i = 0;
        for (int32_t o = 0; o < n_cells; ++o){
            indptr.push_back(i);
            if (mtx.count(o) > 0){
                for (map<int32_t, double>::iterator y = mtx[o].begin(); 
                    y != mtx[o].end(); ++y){
                    indices.push_back(y->first);
                    data.push_back(y->second);
                    i++;
                }
            }
        }
        indptr.push_back(i);
    }
    
    /**
     * Convert a vector representation of CSC sparse matrix format 
     * (used by h5 files) into a map representation
     */
    void h5_reader::mtx2map(vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr,
        map<int32_t, map<int32_t, double> >& mtx){
        
        for (int32_t i = 0; i < n_cells; ++i){
            int64_t start = indptr[i];
            int64_t end = indptr[i + 1];
            for (int j = start; j < end; ++j){
                if (j == start){
                    map<int32_t, double> m;
                    mtx.insert(make_pair(i, m));
                } 
                int32_t gene_idx = indices[j];
                double count = data[j];
                mtx[i].insert(make_pair(gene_idx, count));
            }
        }
    }
    
    /**
     * Load a layer in map format
     */ 
    void h5_reader::load_layer(const string& layername,
        map<int32_t, map<int32_t, double> >& mtx){
        
        vector<double> data;
        vector<int32_t> indices;
        vector<int64_t> indptr;
        
        load_layer(layername, data, indices, indptr);
        mtx2map(data, indices, indptr, mtx);
    }

    /**
     * Write a new layer (same dimensions as main expression matrix) to the
     * file.
     */
    void h5_reader::write_layer(const std::string& layername,
        map<int32_t, map<int32_t, double> >& mtxmap){
        
        vector<double> data;
        vector<int32_t> indices;
        vector<int64_t> indptr;
        map2mtx(mtxmap, data, indices, indptr);
        write_layer(layername, data, indices, indptr);

    }

}

