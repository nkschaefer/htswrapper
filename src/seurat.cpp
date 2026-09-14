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
#include "seurat.h"

using std::cout;
using std::endl;
using namespace std;

namespace sch5{

    bool seurat::is_data_slot(const string& name){
        return (name == "counts" || name == "data" ||
                name == "scale.data");
    }

    /**
     * Constructor
     */
    seurat::seurat(const string& fn) : h5_reader(fn){
        try{
            // Read active assay name
            if (file.exist("active.assay")){
                auto ot = file.getObjectType("active.assay");
                if (ot == HighFive::ObjectType::Dataset){
                    file.getDataSet("active.assay").read(active_assay);
                }
            }
            if (active_assay.empty()){
                active_assay = "RNA";
            }

            // Read cell names
            if (file.exist("cell.names")){
                file.getDataSet("cell.names").read(cell_names);
            }
            else{
                throw runtime_error("cell.names not found");
            }
            n_cells = cell_names.size();

            // Read gene names from active assay
            string feat_path = "assays/" + active_assay + "/features";
            if (file.exist(feat_path)){
                file.getDataSet(feat_path).read(gene_names);
            }
            else{
                throw runtime_error("features not found at " + feat_path);
            }
            n_genes = gene_names.size();
        }
        catch(exception& e){
            throw runtime_error("ERROR loading metadata from h5Seurat file " +
                filename + ": " + e.what());
        }
    }

    /**
     * Destructor
     */
    seurat::~seurat(){

    }

    void seurat::set_active_assay(const string& name){
        active_assay = name;
        // Reload gene names for the new assay
        string feat_path = "assays/" + active_assay + "/features";
        if (file.exist(feat_path)){
            gene_names.clear();
            file.getDataSet(feat_path).read(gene_names);
            n_genes = gene_names.size();
        }
        else{
            throw runtime_error("features not found for assay " + active_assay);
        }
    }

    /**
     * Load an integer metadata column.
     */
    void seurat::load_meta_col(const string& col_name,
        vector<int>& col_vals,
        bool fix_nan){

        string path = "meta.data/" + col_name;
        if (!file.exist(path)){
            throw runtime_error("metadata column " + col_name + " not found");
        }
        auto dataset = file.getDataSet(path);

        if (dataset.hasAttribute("levels")){
            throw runtime_error("metadata column " + col_name +
                " is a factor (categorical); load as string instead");
        }
        load_int_col_aux(dataset, col_vals);
    }

    /**
     * Load a float metadata column.
     */
    void seurat::load_meta_col(const string& col_name,
        vector<double>& col_vals,
        bool fix_nan){

        string path = "meta.data/" + col_name;
        if (!file.exist(path)){
            throw runtime_error("metadata column " + col_name + " not found");
        }
        auto dataset = file.getDataSet(path);
        load_float_col_aux(dataset, col_vals, fix_nan);
    }

    /**
     * Load a string or factor column as string data.
     */
    void seurat::load_meta_col(const string& col_name,
        vector<string>& col_vals,
        bool fix_nan){

        string path = "meta.data/" + col_name;
        if (!file.exist(path)){
            throw runtime_error("metadata column " + col_name + " not found");
        }
        auto ds = file.getDataSet(path);
        auto dtype = ds.getDataType();

        if (dtype.getClass() == HighFive::DataTypeClass::String){
            load_str_col_aux(ds, col_vals);
        }
        else if (ds.hasAttribute("levels")){
            // Factor column: integer codes + levels attribute
            vector<string> levels;
            ds.getAttribute("levels").read(levels);

            vector<int> codes;
            load_int_col_aux(ds, codes);

            col_vals.clear();
            col_vals.reserve(codes.size());
            for (size_t i = 0; i < codes.size(); ++i){
                if (codes[i] < 0 || codes[i] >= (int)levels.size()){
                    col_vals.push_back("");
                }
                else{
                    col_vals.push_back(levels[codes[i]]);
                }
            }
        }
        else{
            throw runtime_error("metadata column " + col_name +
                " is not a string or factor column");
        }
    }

    /**
     * Retrieve metadata column names by type.
     */
    void seurat::get_meta_colnames(vector<string>& str_cols,
        vector<string>& int_cols,
        vector<string>& float_cols){

        if (!file.exist("meta.data")){
            return;
        }
        auto meta = file.getGroup("meta.data");
        vector<string> names = meta.listObjectNames();

        for (size_t i = 0; i < names.size(); ++i){
            auto obj_type = meta.getObjectType(names[i]);
            if (obj_type != HighFive::ObjectType::Dataset){
                continue;
            }
            auto ds = meta.getDataSet(names[i]);
            auto dtype = ds.getDataType();
            auto type_class = dtype.getClass();

            if (type_class == HighFive::DataTypeClass::String){
                str_cols.push_back(names[i]);
            }
            else if (ds.hasAttribute("levels")){
                // Factor → treat as string
                str_cols.push_back(names[i]);
            }
            else if (type_class == HighFive::DataTypeClass::Integer){
                int_cols.push_back(names[i]);
            }
            else if (type_class == HighFive::DataTypeClass::Float){
                float_cols.push_back(names[i]);
            }
        }
    }

    /**
     * Return whether a given metadata column exists.
     */
    bool seurat::has_meta_col(const string& colname){
        string path = "meta.data/" + colname;
        return file.exist(path);
    }

    /**
     * Return the type of a given metadata column.
     */
    short seurat::meta_col_type(const string& col_name){
        string path = "meta.data/" + col_name;
        if (!file.exist(path)){
            throw runtime_error("metadata column " + col_name +
                " does not exist");
        }
        auto ds = file.getDataSet(path);

        if (ds.hasAttribute("levels")){
            return h5_type_str;
        }
        return type_from_ds(ds);
    }

    /**
     * Load a sparse or dense matrix, handling h5Seurat files that lack
     * the anndata-style encoding-type attribute.
     *
     * h5Seurat stores sparse matrices as dgCMatrix (CSC with features as
     * rows, cells as columns).  CSC of (features x cells) is equivalent
     * to CSR of (cells x features), so the data/indices/indptr can be
     * used directly — no conversion needed.
     *
     * If encoding-type IS present (e.g. written by this library), we
     * honour it so a CSC-tagged group still gets converted correctly.
     */
    void seurat::load_seurat_mtx(const string& name,
        vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr){

        auto obj_type = file.getObjectType(name);

        if (obj_type == HighFive::ObjectType::Group){
            auto grp = file.getGroup(name);
            grp.getDataSet("data").read(data);
            grp.getDataSet("indices").read(indices);
            grp.getDataSet("indptr").read(indptr);

            if (grp.hasAttribute("encoding-type")){
                // Attribute present — use standard logic
                if (!check_encoding(grp)){
                    csc_to_csr(data, indices, indptr);
                }
            }
            // No encoding-type: SeuratDisk dgCMatrix.
            // CSC of (features x cells) = CSR of (cells x features).
            // indptr already has n_cells+1 entries and indices are gene
            // indices, so the data is directly usable as CSR.
        }
        else if (obj_type == HighFive::ObjectType::Dataset){
            // Dense — cells x genes, row-major
            vector<double> mtx(n_cells * n_genes);
            auto ds = file.getDataSet(name);
            ds.read_raw<double>(mtx.data());
            int64_t num = 0;
            for (long int i = 0; i < n_cells; ++i){
                indptr.push_back(num);
                const double* row = mtx.data() + i * n_genes;
                for (long int j = 0; j < n_genes; ++j){
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
            throw runtime_error("h5 path " + name +
                " has unknown type; cannot load");
        }
    }

    /**
     * Load expression data.
     */
    void seurat::load_expr(vector<double>& X_data,
        vector<int32_t>& X_indices,
        vector<int64_t>& X_indptr){

        string assay_root = "assays/" + active_assay;

        // If countsname is set, try it first
        if (countsname != ""){
            string path = assay_root + "/" + countsname;
            if (file.exist(path)){
                fprintf(stderr, "Loading expression from %s\n", path.c_str());
                load_seurat_mtx(path, X_data, X_indices, X_indptr);
                return;
            }
            if (file.exist(countsname)){
                fprintf(stderr, "Loading expression from %s\n",
                    countsname.c_str());
                load_seurat_mtx(countsname, X_data, X_indices, X_indptr);
                return;
            }
            throw runtime_error("ERROR: data slot " + countsname +
                " not found in h5Seurat file");
        }

        // Try counts first, then data
        string counts_path = assay_root + "/counts";
        if (file.exist(counts_path)){
            auto ot = file.getObjectType(counts_path);
            if (ot == HighFive::ObjectType::Group ||
                ot == HighFive::ObjectType::Dataset){
                fprintf(stderr, "Loading expression from %s\n",
                    counts_path.c_str());
                load_seurat_mtx(counts_path, X_data, X_indices, X_indptr);
                return;
            }
        }

        string data_path = assay_root + "/data";
        if (file.exist(data_path)){
            auto ot = file.getObjectType(data_path);
            if (ot == HighFive::ObjectType::Group ||
                ot == HighFive::ObjectType::Dataset){
                fprintf(stderr, "Loading expression from %s\n",
                    data_path.c_str());
                load_seurat_mtx(data_path, X_data, X_indices, X_indptr);
                return;
            }
        }

        throw runtime_error("ERROR: could not find expression data in "
            "h5Seurat file " + filename);
    }

    /**
     * List available data layers in the active assay.
     */
    void seurat::list_layers(vector<string>& layernames){
        string assay_root = "assays/" + active_assay;
        if (!file.exist(assay_root)){
            return;
        }
        auto assay_group = file.getGroup(assay_root);
        vector<string> names = assay_group.listObjectNames();
        for (size_t i = 0; i < names.size(); ++i){
            if (is_data_slot(names[i])){
                layernames.push_back(names[i]);
            }
        }
    }

    bool seurat::has_layer(const string& layername){
        string path = "assays/" + active_assay + "/" + layername;
        return file.exist(path);
    }

    void seurat::load_layer(const string& layername,
        vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr){

        string path = "assays/" + active_assay + "/" + layername;
        if (!file.exist(path)){
            throw runtime_error("layer " + layername + " not found in assay " +
                active_assay);
        }
        load_seurat_mtx(path, data, indices, indptr);
    }

    void seurat::write_layer(const string& layername,
        vector<double>& data,
        vector<int32_t>& indices,
        vector<int64_t>& indptr){

        fprintf(stderr, "Writing layer \"%s\" to h5Seurat...\n",
            layername.c_str());

        string assay_root = "assays/" + active_assay;
        if (!file.exist("assays")){
            file.createGroup("assays");
        }
        if (!file.exist(assay_root)){
            file.createGroup(assay_root);
        }

        string path = assay_root + "/" + layername;
        HighFive::Group g = file.createGroup(path);

        // h5Seurat uses CSC (features x cells), but CSR (cells x features)
        // has the same indptr/indices/data layout, so write directly.
        g.createAttribute<string>("encoding-type", string("csr_matrix"));
        g.createAttribute<string>("encoding-version", string("0.1.0"));
        vector<int64_t> shape = {(int64_t)n_cells, (int64_t)n_genes};
        g.createAttribute<vector<int64_t>>("shape", shape);

        g.createDataSet("data", data);
        g.createDataSet("indices", indices);
        g.createDataSet("indptr", indptr);
    }

}
