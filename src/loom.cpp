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
#include <limits>
#include <unordered_set>
#include <cstdlib>
#include <utility>
#include <regex>
#include <highfive/H5File.hpp>
#include "h5_reader.h"
#include "loom.h"

using std::cout;
using std::endl;
using namespace std;

namespace sch5{
    
    /**
     * Loom files to weird stuff to cell barcode strings
     */
    string loom::parse_loom_bc(const string& loom_bc) {
        static const regex re(R"(^(?:[^:]*:)?([ACGTN]+)x?$)");
        smatch m;
        if (regex_search(loom_bc, m, re)){
            return m[1].str();
        }
        else{
            return loom_bc;
        }
    }

    /**
     * Constructor
     */
    loom::loom(const string& fn) : h5_reader(fn){
        // Read the obs and var names
        try{
            // Load gene and barcode names
            file.getDataSet("/row_attrs/Gene").read(gene_names);
            file.getDataSet("/col_attrs/CellID").read(cell_names);
            
            n_cells = cell_names.size();
            n_genes = gene_names.size();

            // Clean up weird cell barcodes
            for (int i = 0; i < cell_names.size(); ++i){
                cell_names[i] = parse_loom_bc(cell_names[i]);
            }
        }
        catch(exception& e){
            throw runtime_error("ERROR loading cell/gene metadata from loom file " + filename);
        }
    }
    
    /**
     * Destructor
     */
    loom::~loom(){

    }

    /**
     * Load an integer metadata column (to be implemented by child classes)
     */
    void loom::load_meta_col(const string& col_name,
        vector<int>& col_vals,
        bool fix_nan){
        string path = "/col_attrs/" + col_name;
        auto dataset = file.getDataSet(path);
        // loom can use enum to store boolean; we will treat this as integer
        auto datatype = dataset.getDataType();
        if (datatype.getClass() == HighFive::DataTypeClass::Enum){
            vector<long long> codes;
            dataset.read(codes);
            col_vals.clear();
            col_vals.reserve(codes.size());
            for (long long c : codes) {
                col_vals.push_back(static_cast<int>(c));
            }
        }
        else{
            load_int_col_aux(dataset, col_vals);
        }
    }

    /**
     * Load a float metadata column (to be implemented by child classes)
     */
    void loom::load_meta_col(const string& col_name,
        vector<double>& col_vals,
        bool fix_nan){
        
        string path = "/col_attrs/" + col_name;
        auto dataset = file.getDataSet(path);
        load_float_col_aux(dataset, col_vals, fix_nan);
    }

    /**
     * Load a string or categorical column (as string data) from an h5 file.
     * To be implemented by child classes.
     */
    void loom::load_meta_col(const string& col_name,
        vector<string>& col_vals,
        bool fix_nan){
        
        string path = "/col_attrs/" + col_name;
        auto dataset = file.getDataSet(path);
        load_str_col_aux(dataset, col_vals);

    }
    
    /**
     * Retrieve layer names.
     */
    void loom::list_layers(vector<string>& names){
        if (file.exist("layers")){
            auto layers = file.getGroup("layers");
            names = layers.listObjectNames();
        }
        else{
            // Not an error for layers not to exist
        }
    }

    /**
     * Retrieve metadata columns by type.
     */
    void loom::get_meta_colnames(vector<string>& str_cols,
        vector<string>& int_cols,
        vector<string>& float_cols){
        
        // Load obs
        auto obs = file.getGroup("col_attrs");
        
        // Read existing column order
        vector<string> col_order;
        obs.getAttribute("column-order").read(col_order);
        
        for (vector<string>::iterator col = col_order.begin(); col != col_order.end(); ++col){
            bool int_type = false;
            bool float_type = false;
            auto obj_type = obs.getObjectType(*col);
            if (obj_type == HighFive::ObjectType::Group){
                throw runtime_error("Error: column " + *col + " is a group (illegal in loom spec)");
            }
            else if (obj_type == HighFive::ObjectType::Dataset){
                auto dataset = obs.getDataSet(*col);
                auto dtype = dataset.getDataType();
                auto type_class = dtype.getClass();
                if (type_class == HighFive::DataTypeClass::String){
                    str_cols.push_back(*col);
                }
                else if (type_class == HighFive::DataTypeClass::Integer || type_class == HighFive::DataTypeClass::Enum){
                    int_cols.push_back(*col);
                }
                else if (type_class == HighFive::DataTypeClass::Float){
                    float_cols.push_back(*col);
                }
            }
        }
    }
    
    /**
     * Return whether a given metadata column exists.
     */
    bool loom::has_meta_col(const string& colname){
        // Load obs
        auto obs = file.getGroup("col_attrs");
        return obs.exist(colname);
    }

    /**
     * Return the type of a given metadata column.
     */
    short loom::meta_col_type(const string& col_name){
        auto obs = file.getGroup("col_attrs");
        if (!obs.exist(col_name)){
            throw runtime_error("Error: column " + col_name + " does not exist in metadata"); 
        }
        HighFive::ObjectType ot = obs.getObjectType(col_name);
        if (ot == HighFive::ObjectType::Dataset){
            HighFive::DataSet ds = obs.getDataSet(col_name);
            if (ds.hasAttribute("categories")){
                // categorical column, unsupported type
                return h5_type_unknown;
            }
            else{
                auto dtype = ds.getDataType();
                auto type_class = dtype.getClass();
                if (type_class == HighFive::DataTypeClass::Enum){
                    return h5_type_int;
                }
                return type_from_ds(ds);
            }
        }
        return h5_type_unknown;
    }

    void loom::load_expr(vector<double>& X_data,
        vector<int32_t>& X_indices,
        vector<int64_t>& X_indptr){
        load_mtx_dense("matrix", X_data, X_indices, X_indptr);
    }
    
    void loom::load_expr(map<int32_t, map<int32_t, double> >& mtx){
        load_mtx_dense("matrix", mtx);
    }

    void loom::load_layer(const string& layername,
        vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr){
        
        try{
            auto layers = file.getGroup("layers");
        }
        catch(HighFive::Exception& e){
            throw runtime_error("No layers in file");
        }
        string path = "/layers/" + layername;
        load_mtx_dense(path, data, indices, indptr);
    }

    void loom::load_layer(const string& layername,
        map<int32_t, map<int32_t, double> >& mtx){
        try{
            auto layers = file.getGroup("layers");
        }
        catch(HighFive::Exception& e){
            throw runtime_error("No layers in file");
        }
        string path = "/layers/" + layername;
        load_mtx_dense(path, mtx);
    }
    
    void loom::load_spliced(vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr){
        string path = "/spliced";
        load_mtx_dense(path, data, indices, indptr);
    }

    void loom::load_spliced(map<int32_t, map<int32_t, double> >& mtx){
        string path = "/spliced";
        load_mtx_dense(path, mtx);
    }
    
    void loom::load_unspliced(vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr){
        string path = "/unspliced";
        load_mtx_dense(path, data, indices, indptr);
    }

    void loom::load_unspliced(map<int32_t, map<int32_t, double> >& mtx){
        string path = "/unspliced";
        load_mtx_dense(path, mtx);
    }

    void loom::write_layer(const string& layername,
        vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr){

        fprintf(stderr, "Writing layer \"%s\" to h5ad...\n", layername.c_str());

        HighFive::Group layers;
        if (!file.exist("layers")){
           layers = file.getGroup("layers");
        }
        else{
           layers = file.createGroup("layers");
        }
        if (!layers.hasAttribute("encoding-type")) {
            layers.createAttribute<string>("encoding-type", string("dict"));
            layers.createAttribute<string>("encoding-version", string("0.1.0"));
        }
        string path = "layers/" + layername;
        write_mtx_dense(path, data, indices, indptr, true, false);
    }

    void loom::write_layer(const std::string& layername,
        map<int32_t, map<int32_t, double> >& mapmtx){
        
        fprintf(stderr, "Writing layer \"%s\" to h5ad...\n", layername.c_str());

        HighFive::Group layers;
        if (!file.exist("layers")){
            layers = file.getGroup("layers");
        }
        else{
            layers = file.createGroup("layers");
        }
        if (!layers.hasAttribute("encoding-type")) {
            layers.createAttribute<string>("encoding-type", string("dict"));
            layers.createAttribute<string>("encoding-version", string("0.1.0"));
        }
        string path = "layers/" + layername;
        write_mtx_dense(path, mapmtx, true, false);
    }
}

