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
#include <highfive/H5File.hpp>
#include "h5_reader.h"
#include "anndata.h"

using std::cout;
using std::endl;
using namespace std;

namespace sch5{

    /**
     * Constructor
     */
    anndata::anndata(const string& fn) : h5_reader(fn){
        // Read the obs and var names
        try{
            // Load obs
            auto obs = file.getGroup("obs");
            
            // Load index names    
            string obs_idx_name;
            obs.getAttribute("_index").read(obs_idx_name);
            if (obs.getObjectType(obs_idx_name) == HighFive::ObjectType::Dataset){
                // string column
                obs.getDataSet(obs_idx_name).read(cell_names);
            }
            else{
                // categorical column
                auto col_group = obs.getGroup(obs_idx_name);
                col_group.getDataSet("values").read(cell_names);
            }
            n_cells = cell_names.size();

            // Load var
            auto var = file.getGroup("var");
            
            // Load index names
            string var_idx_name;
            var.getAttribute("_index").read(var_idx_name);
            if (var.getObjectType(var_idx_name) == HighFive::ObjectType::Dataset){
                var.getDataSet(var_idx_name).read(gene_names);
            }
            else{
                auto col_group = var.getGroup(var_idx_name);
                col_group.getDataSet("values").read(gene_names);
            }
            n_genes = gene_names.size();
        }
        catch(exception& e){
            throw runtime_error("ERROR loading obs and var metadata from " + filename + ". \
    If this file is old, try opening in a recent version of scanpy and rewriting.");
        }
    }
    
    /**
     * Destructor
     */
    anndata::~anndata(){

    }

    /**
     * Load an integer metadata column (to be implemented by child classes)
     */
    void anndata::load_meta_col(const string& col_name,
        vector<int>& col_vals,
        bool fix_nan){
        string path = "/obs/" + col_name;
        auto dataset = file.getDataSet(path);
        auto dtype = dataset.getDataType();
        
        if (dataset.hasAttribute("categories")){
            throw runtime_error("Error: metadata column has attribute \"categories.\" \
This is not a valid integer column.");
        }
        load_int_col_aux(dataset, col_vals);
    }

    /**
     * Load a float metadata column (to be implemented by child classes)
     */
    void anndata::load_meta_col(const string& col_name,
        vector<double>& col_vals,
        bool fix_nan){
        string path = "/obs/" + col_name;
        auto dataset = file.getDataSet(path);
        auto dtype = dataset.getDataType();
        load_float_col_aux(dataset, col_vals, fix_nan);
    }

    /**
     * Load a string or categorical column (as string data) from an h5 file.
     * To be implemented by child classes.
     */
    void anndata::load_meta_col(const string& col_name,
        vector<string>& col_vals,
        bool fix_nan){
        string path = "/obs/" + col_name;
        HighFive::ObjectType ot = file.getObjectType(path);
        if (ot == HighFive::ObjectType::Dataset){
            HighFive::DataSet ds = file.getDataSet(path);
            if (ds.hasAttribute("categories")){
                // categorical column.
                throw runtime_error("Error: metadata column " + col_name + " has attribute \"categories.\" \
Please open and re-save this file using a recent version of scanpy.");
            }
            else{
                auto dtype = ds.getDataType();
                if (dtype.getClass() == HighFive::DataTypeClass::String){
                    // string column
                    load_str_col_aux(ds, col_vals);
                }
                else{
                    throw runtime_error("Error: metadata column " + col_name + " has unknown type.");
                }
            }
        }
        else if (ot == HighFive::ObjectType::Group){
            // categorical column
            auto col_group = file.getGroup(path);
            auto names =  col_group.listObjectNames();
            bool has_values = find(names.begin(), names.end(), "values") != names.end();
            bool has_mask = find(names.begin(), names.end(), "mask") != names.end();
            bool has_cats = find(names.begin(), names.end(), "categories") != names.end();
            bool has_codes = find(names.begin(), names.end(), "codes") != names.end();
            if (has_values && has_mask){
                // Masked strings
                HighFive::DataSet ds = col_group.getDataSet("values");
                load_str_col_aux(ds, col_vals);
                vector<bool> mask;
                col_group.getDataSet("mask").read(mask);
                for (int i = 0; i < mask.size(); ++i){
                    // Set missing values to empty string
                    if (!mask[i]){
                        col_vals[i] = "";
                    }
                }
            }
            else if (has_cats && has_codes){
                // Factor variable
                vector<string> cats;
                HighFive::DataSet ds = col_group.getDataSet("categories");
                load_str_col_aux(ds, cats);
                auto codes_ds = col_group.getDataSet("codes");
                auto codes_dtype = codes_ds.getDataType();
                std::vector<int> codes;
                if (codes_dtype == HighFive::AtomicType<int8_t>()) {
                    std::vector<int8_t> buf;
                    codes_ds.read(buf);
                    codes.assign(buf.begin(), buf.end());
                } else if (codes_dtype == HighFive::AtomicType<int16_t>()) {
                    std::vector<int16_t> buf;
                    codes_ds.read(buf);
                    codes.assign(buf.begin(), buf.end());
                } else if (codes_dtype == HighFive::AtomicType<int32_t>()) {
                    codes_ds.read(codes);
                }
                col_vals.reserve(codes.size());
                for (int i = 0; i < codes.size(); ++i){
                    if (codes[i] < 0){
                        col_vals.push_back("");
                    }
                    else{
                        col_vals.push_back(cats[codes[i]]);
                    }
                }
            }
            else{
                throw runtime_error("Error: unable to interpret group data for column " + col_name);
            }
            
        }
        else{
            throw runtime_error("Error: unknown object type for column " + col_name);
        }
    }
    
    /**
     * Retrieve layer names.
     */
    void anndata::list_layers(vector<string>& names){
        if (file.exist("layers")){
            auto layers = file.getGroup("layers");
            names = layers.listObjectNames();
        }
        else{
            // Okay for no layers to exist
        }
    }

    /**
     * Retrieve metadata columns by type.
     */
    void anndata::get_meta_colnames(vector<string>& str_cols,
        vector<string>& int_cols,
        vector<string>& float_cols){
        
        // Load obs
        auto obs = file.getGroup("obs");
        
        // Read existing column order
        vector<string> col_order;
        obs.getAttribute("column-order").read(col_order);
        
        for (vector<string>::iterator col = col_order.begin(); col != col_order.end(); ++col){
            bool int_type = false;
            bool float_type = false;
            auto obj_type = obs.getObjectType(*col);
            if (obj_type == HighFive::ObjectType::Group){
                // Categorical
                str_cols.push_back(*col);
            }
            else if (obj_type == HighFive::ObjectType::Dataset){
                auto dataset = obs.getDataSet(*col);
                auto dtype = dataset.getDataType();
                auto type_class = dtype.getClass();
                if (type_class == HighFive::DataTypeClass::String){
                    str_cols.push_back(*col);
                }
                else if (type_class == HighFive::DataTypeClass::Integer){
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
    bool anndata::has_meta_col(const string& colname){
        // Load obs
        auto obs = file.getGroup("obs");
        return obs.exist(colname);
    }

    /**
     * Return the type of a given metadata column.
     */
    short anndata::meta_col_type(const string& col_name){
        auto obs = file.getGroup("obs");
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
                return type_from_ds(ds);
            }
        }
        else if (ot == HighFive::ObjectType::Group){
            // categorical column
            auto col_group = obs.getGroup(col_name);
            auto names =  col_group.listObjectNames();
            bool has_values = find(names.begin(), names.end(), "values") != names.end();
            bool has_mask = find(names.begin(), names.end(), "mask") != names.end();
            bool has_cats = find(names.begin(), names.end(), "categories") != names.end();
            bool has_codes = find(names.begin(), names.end(), "codes") != names.end();
            if (has_values && has_mask){
                return sch5::h5_type_str;
            }
            else if (has_cats && has_codes){
                return sch5::h5_type_str;
            }
            else{
                return sch5::h5_type_unknown;
            }
        }
        return sch5::h5_type_unknown;
    }

    void anndata::load_expr(vector<double>& X_data,
        vector<int32_t>& X_indices,
        vector<int64_t>& X_indptr){
        
        // Strategy: 
        // countsname set? 
        //  yes: load it or fail
        //  no: look for a layer called counts
        //      found: load it
        //      no: look for adata.raw
        //        found: load it
        //        no: load adata.X

        if (file.exist("layers")){
            auto layers = file.getGroup("layers");
            vector<string> layer_names = layers.listObjectNames();
            bool found = false;
            for (int i = 0; i < layer_names.size(); ++i){
                if (layer_names[i] == countsname || 
                    (countsname == "" && layer_names[i] == "counts")){
                    fprintf(stderr, "Loading expression from layer %s\n", 
                        layer_names[i].c_str());
                    // Found
                    string path = "/layers/" + layer_names[i];
                    load_mtx(path, X_data, X_indices, X_indptr);            
                    return;
                }
            }
        }

        if (countsname != ""){
            // Try to retrieve it outside of the "layers" context.
            if (file.exist(countsname) && file.getObjectType(countsname) == HighFive::ObjectType::Group){
                string path = "/" + countsname;
                load_mtx(path, X_data, X_indices, X_indptr);
                return;
            }
            else{
                fprintf(stderr, "ERROR: layer/data set %s not found in anndata.\n", countsname.c_str());
                exit(1);
            }
        }
        else{
            // Look for raw.
            if (file.exist("raw/X") && 
                (file.getObjectType("raw/X") == HighFive::ObjectType::Group || 
                 file.getObjectType("raw/X") == HighFive::ObjectType::Dataset)){
                fprintf(stderr, "Loading expression from /raw/X\n");
                string path = "/raw/X";
                load_mtx(path, X_data, X_indices, X_indptr);
                return;
            }
            else if (file.exist("raw/counts") && 
                (file.getObjectType("raw/counts") == HighFive::ObjectType::Group || 
                 file.getObjectType("raw/counts") == HighFive::ObjectType::Dataset)){
                fprintf(stderr, "Loading expression data from /raw/counts\n");
                string path = "/raw/counts";
                load_mtx(path, X_data, X_indices, X_indptr);    
                return;
            }
            else if (file.exist("raw") && 
                (file.getObjectType("raw") == HighFive::ObjectType::Group || 
                 file.getObjectType("raw") == HighFive::ObjectType::Dataset)){
                fprintf(stderr, "Loading expression data from /raw\n");
                string path = "/raw";
                load_mtx(path, X_data, X_indices, X_indptr);
                return;
            }
            else{
                if (file.exist("X") && 
                    (file.getObjectType("X") == HighFive::ObjectType::Group ||
                     file.getObjectType("X") == HighFive::ObjectType::Dataset)) {
                    fprintf(stderr, "Loading expression from /X\n");
                    string path = "/X";
                    load_mtx(path, X_data, X_indices, X_indptr);
                    return;
                }
                else{
                    fprintf(stderr, "ERROR: could not find expression data. Please specify where \
    raw expression counts are stored in %s\n", filename.c_str());
                    exit(1);
                }
            }
        }
    }
    
    void anndata::load_layer(const string& layername,
        vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr){
        try{
            auto layers = file.getGroup("layers");
        }
        catch (HighFive::Exception& e){
            throw runtime_error("Error: no layers present in file");
        }
        string path = "/layers/" + layername;
        load_mtx(path, data, indices, indptr); 
    }
    
    bool anndata::has_layer(const string& layername){
        string path = "/layers/" + layername;
        return file.exist(path); 
    }

    void anndata::write_layer(const string& layername,
        vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr){

        fprintf(stderr, "Writing layer \"%s\" to h5ad...\n", layername.c_str());

        HighFive::Group layers;
        if (file.exist("layers")){
           layers = file.getGroup("layers");
        }
        else{
           layers = file.createGroup("layers");
        }
        if (!layers.hasAttribute("encoding-type")) {
            layers.createAttribute<string>("encoding-type", string("dict"));
            layers.createAttribute<string>("encoding-version", string("0.1.0"));
        }
        
        HighFive::Group g = layers.createGroup(layername);
        g.createAttribute<string>("encoding-type", string("csr_matrix"));
        g.createAttribute<string>("encoding-version", string("0.1.0"));
        vector<int64_t> shape = {(int64_t)n_cells, (int64_t)n_genes};
        g.createAttribute<vector<int64_t>>("shape", shape);

        // Write it to the new slot.
        g.createDataSet("data", data);
        g.createDataSet("indices", indices);
        g.createDataSet("indptr", indptr);
    }
}

