#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <stdexcept>
#include <highfive/H5File.hpp>
#include "h5_reader.h"
#include "anndata.h"
#include "loom.h"
#include "seurat.h"

using namespace std;

static int tests_run = 0;
static int tests_passed = 0;

#define ASSERT_EQ(a, b, msg) do { \
    if ((a) != (b)) { \
        cerr << "FAIL: " << msg << " (expected " << (b) << ", got " << (a) << ")" << endl; \
        return false; \
    } \
} while(0)

#define ASSERT_NEAR(a, b, tol, msg) do { \
    if (fabs((a) - (b)) > (tol)) { \
        cerr << "FAIL: " << msg << " (expected ~" << (b) << ", got " << (a) << ")" << endl; \
        return false; \
    } \
} while(0)

#define ASSERT_TRUE(cond, msg) do { \
    if (!(cond)) { \
        cerr << "FAIL: " << msg << endl; \
        return false; \
    } \
} while(0)

#define ASSERT_THROW(expr, msg) do { \
    bool threw = false; \
    try { expr; } catch (...) { threw = true; } \
    if (!threw) { \
        cerr << "FAIL: expected exception: " << msg << endl; \
        return false; \
    } \
} while(0)

#define RUN_TEST(fn) do { \
    tests_run++; \
    if (fn()) { \
        tests_passed++; \
        cout << "  PASS: " << #fn << endl; \
    } \
    else { \
        cout << "  FAIL: " << #fn << endl; \
    } \
} while(0)

// =====================================================================
// Helper: create a minimal anndata h5ad file
// =====================================================================
static string create_test_anndata(const string& path,
    int n_cells, int n_genes,
    bool sparse, bool csc = false){

    {
        HighFive::File f(path, HighFive::File::Truncate);

        // --- obs ---
        auto obs = f.createGroup("obs");
        vector<string> barcodes;
        for (int i = 0; i < n_cells; ++i){
            barcodes.push_back("CELL_" + to_string(i));
        }
        obs.createDataSet("_index_col", barcodes);
        obs.createAttribute<string>("_index", string("_index_col"));

        vector<string> col_order;

        // int metadata column
        vector<int32_t> int_col;
        for (int i = 0; i < n_cells; ++i) int_col.push_back(i * 10);
        obs.createDataSet("n_counts", int_col);
        col_order.push_back("n_counts");

        // float metadata column
        vector<double> float_col;
        for (int i = 0; i < n_cells; ++i) float_col.push_back(i * 0.5);
        obs.createDataSet("pct_mito", float_col);
        col_order.push_back("pct_mito");

        // string metadata column
        vector<string> str_col;
        for (int i = 0; i < n_cells; ++i){
            str_col.push_back(i % 2 == 0 ? "typeA" : "typeB");
        }
        obs.createDataSet("cell_type", str_col);
        col_order.push_back("cell_type");

        // categorical column (codes + categories, stored as group)
        auto cat_group = obs.createGroup("batch");
        vector<string> cats = {"batch0", "batch1", "batch2"};
        cat_group.createDataSet("categories", cats);
        vector<int8_t> codes;
        for (int i = 0; i < n_cells; ++i) codes.push_back(i % 3);
        cat_group.createDataSet("codes", codes);
        col_order.push_back("batch");

        obs.createAttribute("column-order", col_order);

        // --- var ---
        auto var = f.createGroup("var");
        vector<string> genes;
        for (int i = 0; i < n_genes; ++i){
            genes.push_back("Gene_" + to_string(i));
        }
        var.createDataSet("_index_col", genes);
        var.createAttribute<string>("_index", string("_index_col"));

        // --- X (expression matrix) ---
        // Diagonal-ish sparse: cell i has a nonzero at gene (i % n_genes)
        if (sparse){
            auto xgrp = f.createGroup("X");
            vector<double> data;
            vector<int32_t> indices;
            vector<int64_t> indptr;

            if (!csc){
                // CSR: rows = cells
                for (int c = 0; c < n_cells; ++c){
                    indptr.push_back((int64_t)data.size());
                    int g = c % n_genes;
                    data.push_back((double)(c + 1));
                    indices.push_back((int32_t)g);
                }
                indptr.push_back((int64_t)data.size());
                xgrp.createAttribute<string>("encoding-type", string("csr_matrix"));
            }
            else{
                // CSC: columns = genes
                for (int g = 0; g < n_genes; ++g){
                    indptr.push_back((int64_t)data.size());
                    for (int c = 0; c < n_cells; ++c){
                        if (c % n_genes == g){
                            data.push_back((double)(c + 1));
                            indices.push_back((int32_t)c);
                        }
                    }
                }
                indptr.push_back((int64_t)data.size());
                xgrp.createAttribute<string>("encoding-type", string("csc_matrix"));
            }
            xgrp.createAttribute<string>("encoding-version", string("0.1.0"));
            vector<int64_t> shape = {(int64_t)n_cells, (int64_t)n_genes};
            xgrp.createAttribute("shape", shape);
            xgrp.createDataSet("data", data);
            xgrp.createDataSet("indices", indices);
            xgrp.createDataSet("indptr", indptr);
        }
        else{
            // Dense: cells x genes (row-major)
            vector<double> dense(n_cells * n_genes, 0.0);
            for (int c = 0; c < n_cells; ++c){
                int g = c % n_genes;
                dense[c * n_genes + g] = (double)(c + 1);
            }
            vector<size_t> dims = {(size_t)n_cells, (size_t)n_genes};
            auto ds = f.createDataSet<double>("X", HighFive::DataSpace(dims));
            ds.write_raw(dense.data());
            ds.createAttribute<string>("encoding-type", string("array"));
            ds.createAttribute<string>("encoding-version", string("0.2.0"));
        }
    }
    return path;
}

// =====================================================================
// Helper: create a minimal loom file
// =====================================================================
static string create_test_loom(const string& path,
    int n_cells, int n_genes){

    {
        HighFive::File f(path, HighFive::File::Truncate);

        // Gene names (row attrs)
        auto row_attrs = f.createGroup("row_attrs");
        vector<string> genes;
        for (int i = 0; i < n_genes; ++i){
            genes.push_back("Gene_" + to_string(i));
        }
        row_attrs.createDataSet("Gene", genes);

        // Cell barcodes (col attrs)
        auto col_attrs = f.createGroup("col_attrs");
        vector<string> barcodes;
        string bases = "ACGT";
        for (int i = 0; i < n_cells; ++i){
            string bc = "sample:";
            for (int j = 0; j < 8; ++j){
                bc += bases[(i + j) % 4];
            }
            bc += "x";
            barcodes.push_back(bc);
        }
        col_attrs.createDataSet("CellID", barcodes);

        // int metadata
        vector<int32_t> clusters;
        for (int i = 0; i < n_cells; ++i) clusters.push_back(i % 3);
        col_attrs.createDataSet("ClusterID", clusters);

        // float metadata
        vector<double> scores;
        for (int i = 0; i < n_cells; ++i) scores.push_back(i * 1.5);
        col_attrs.createDataSet("Score", scores);

        vector<string> col_order = {"ClusterID", "Score"};
        col_attrs.createAttribute("column-order", col_order);

        // Main matrix: genes x cells (loom convention, row-major)
        vector<double> dense(n_genes * n_cells, 0.0);
        for (int c = 0; c < n_cells; ++c){
            int g = c % n_genes;
            dense[g * n_cells + c] = (double)(c + 1);
        }
        vector<size_t> dims = {(size_t)n_genes, (size_t)n_cells};
        auto ds = f.createDataSet<double>("matrix", HighFive::DataSpace(dims));
        ds.write_raw(dense.data());

        // A dense layer (spliced)
        auto layers = f.createGroup("layers");
        vector<double> spliced(n_genes * n_cells, 0.0);
        for (int c = 0; c < n_cells; ++c){
            int g = (c + 1) % n_genes;
            spliced[g * n_cells + c] = (double)(c + 10);
        }
        auto sds = f.createDataSet<double>("layers/spliced",
            HighFive::DataSpace(dims));
        sds.write_raw(spliced.data());
    }
    return path;
}

// =====================================================================
// Tests: anndata
// =====================================================================

bool test_anndata_cell_gene_names(){
    string path = "/tmp/test_anndata_basic.h5ad";
    create_test_anndata(path, 5, 3, true);
    sch5::anndata ad(path);

    ASSERT_EQ(ad.n_cells, 5, "n_cells");
    ASSERT_EQ(ad.n_genes, 3, "n_genes");
    ASSERT_EQ(ad.cell_names[0], "CELL_0", "cell_names[0]");
    ASSERT_EQ(ad.cell_names[4], "CELL_4", "cell_names[4]");
    ASSERT_EQ(ad.gene_names[0], "Gene_0", "gene_names[0]");
    ASSERT_EQ(ad.gene_names[2], "Gene_2", "gene_names[2]");
    remove(path.c_str());
    return true;
}

bool test_anndata_load_expr_sparse_csr(){
    string path = "/tmp/test_anndata_csr.h5ad";
    create_test_anndata(path, 5, 3, true, false);
    sch5::anndata ad(path);

    vector<double> data;
    vector<int32_t> indices;
    vector<int64_t> indptr;
    ad.load_expr(data, indices, indptr);

    ASSERT_EQ((int)indptr.size(), 6, "indptr size (n_cells+1)");
    ASSERT_EQ((int)data.size(), 5, "data size (one nonzero per cell)");

    // Cell 0 → gene 0, value 1
    ASSERT_EQ(indices[0], 0, "cell 0 gene index");
    ASSERT_NEAR(data[0], 1.0, 1e-9, "cell 0 value");

    // Cell 3 → gene 0 (3%3==0), value 4
    ASSERT_EQ(indices[3], 0, "cell 3 gene index");
    ASSERT_NEAR(data[3], 4.0, 1e-9, "cell 3 value");

    remove(path.c_str());
    return true;
}

bool test_anndata_load_expr_sparse_csc(){
    string path = "/tmp/test_anndata_csc.h5ad";
    create_test_anndata(path, 4, 3, true, true);
    sch5::anndata ad(path);

    vector<double> data;
    vector<int32_t> indices;
    vector<int64_t> indptr;
    ad.load_expr(data, indices, indptr);

    // After CSC→CSR conversion, should have same logical entries
    ASSERT_EQ((int)indptr.size(), 5, "indptr size after csc_to_csr");
    ASSERT_EQ((int)data.size(), 4, "data size");

    // Cell 0 → gene 0, value 1 (CSR: first row)
    ASSERT_EQ(indptr[0], 0, "indptr[0]");
    ASSERT_TRUE(indptr[1] - indptr[0] >= 1, "cell 0 has at least one entry");

    remove(path.c_str());
    return true;
}

bool test_anndata_load_expr_dense(){
    string path = "/tmp/test_anndata_dense.h5ad";
    create_test_anndata(path, 4, 3, false);
    sch5::anndata ad(path);

    vector<double> data;
    vector<int32_t> indices;
    vector<int64_t> indptr;
    ad.load_expr(data, indices, indptr);

    ASSERT_EQ((int)indptr.size(), 5, "indptr size");
    ASSERT_EQ((int)data.size(), 4, "data size");

    ASSERT_EQ(indices[0], 0, "cell 0 gene index");
    ASSERT_NEAR(data[0], 1.0, 1e-9, "cell 0 value");

    remove(path.c_str());
    return true;
}

bool test_anndata_load_meta_int(){
    string path = "/tmp/test_anndata_meta_int.h5ad";
    create_test_anndata(path, 4, 2, true);
    sch5::anndata ad(path);

    vector<int> vals;
    ad.load_meta_col("n_counts", vals);
    ASSERT_EQ((int)vals.size(), 4, "int col size");
    ASSERT_EQ(vals[0], 0, "n_counts[0]");
    ASSERT_EQ(vals[2], 20, "n_counts[2]");

    remove(path.c_str());
    return true;
}

bool test_anndata_load_meta_float(){
    string path = "/tmp/test_anndata_meta_float.h5ad";
    create_test_anndata(path, 4, 2, true);
    sch5::anndata ad(path);

    vector<double> vals;
    ad.load_meta_col("pct_mito", vals);
    ASSERT_EQ((int)vals.size(), 4, "float col size");
    ASSERT_NEAR(vals[0], 0.0, 1e-9, "pct_mito[0]");
    ASSERT_NEAR(vals[3], 1.5, 1e-9, "pct_mito[3]");

    remove(path.c_str());
    return true;
}

bool test_anndata_load_meta_string(){
    string path = "/tmp/test_anndata_meta_str.h5ad";
    create_test_anndata(path, 4, 2, true);
    sch5::anndata ad(path);

    vector<string> vals;
    ad.load_meta_col("cell_type", vals);
    ASSERT_EQ((int)vals.size(), 4, "str col size");
    ASSERT_EQ(vals[0], "typeA", "cell_type[0]");
    ASSERT_EQ(vals[1], "typeB", "cell_type[1]");

    remove(path.c_str());
    return true;
}

bool test_anndata_load_meta_categorical(){
    string path = "/tmp/test_anndata_meta_cat.h5ad";
    create_test_anndata(path, 6, 2, true);
    sch5::anndata ad(path);

    vector<string> vals;
    ad.load_meta_col("batch", vals);
    ASSERT_EQ((int)vals.size(), 6, "categorical col size");
    ASSERT_EQ(vals[0], "batch0", "batch[0]");
    ASSERT_EQ(vals[1], "batch1", "batch[1]");
    ASSERT_EQ(vals[2], "batch2", "batch[2]");
    ASSERT_EQ(vals[3], "batch0", "batch[3]");

    remove(path.c_str());
    return true;
}

bool test_anndata_has_meta_col(){
    string path = "/tmp/test_anndata_has_col.h5ad";
    create_test_anndata(path, 3, 2, true);
    sch5::anndata ad(path);

    ASSERT_TRUE(ad.has_meta_col("n_counts"), "has n_counts");
    ASSERT_TRUE(ad.has_meta_col("pct_mito"), "has pct_mito");
    ASSERT_TRUE(ad.has_meta_col("cell_type"), "has cell_type");
    ASSERT_TRUE(ad.has_meta_col("batch"), "has batch");
    ASSERT_TRUE(!ad.has_meta_col("nonexistent"), "!has nonexistent");

    remove(path.c_str());
    return true;
}

bool test_anndata_meta_col_type(){
    string path = "/tmp/test_anndata_coltype.h5ad";
    create_test_anndata(path, 3, 2, true);
    sch5::anndata ad(path);

    ASSERT_EQ(ad.meta_col_type("n_counts"), sch5::h5_type_int, "n_counts type");
    ASSERT_EQ(ad.meta_col_type("pct_mito"), sch5::h5_type_float, "pct_mito type");
    ASSERT_EQ(ad.meta_col_type("cell_type"), sch5::h5_type_str, "cell_type type");
    ASSERT_EQ(ad.meta_col_type("batch"), sch5::h5_type_str, "batch type (categorical)");

    remove(path.c_str());
    return true;
}

bool test_anndata_get_meta_colnames(){
    string path = "/tmp/test_anndata_colnames.h5ad";
    create_test_anndata(path, 3, 2, true);
    sch5::anndata ad(path);

    vector<string> str_cols, int_cols, float_cols;
    ad.get_meta_colnames(str_cols, int_cols, float_cols);

    ASSERT_EQ((int)int_cols.size(), 1, "int cols count");
    ASSERT_EQ(int_cols[0], "n_counts", "int col name");
    ASSERT_EQ((int)float_cols.size(), 1, "float cols count");
    ASSERT_EQ(float_cols[0], "pct_mito", "float col name");
    // str_cols should have cell_type + batch (categorical)
    ASSERT_TRUE((int)str_cols.size() >= 2, "str cols count >= 2");

    remove(path.c_str());
    return true;
}

bool test_anndata_list_layers_empty(){
    string path = "/tmp/test_anndata_nolayers.h5ad";
    create_test_anndata(path, 3, 2, true);
    sch5::anndata ad(path);

    vector<string> layers;
    ad.list_layers(layers);
    ASSERT_EQ((int)layers.size(), 0, "no layers");

    remove(path.c_str());
    return true;
}

bool test_anndata_load_expr_not_found(){
    string path = "/tmp/test_anndata_notfound.h5ad";
    create_test_anndata(path, 3, 2, true);
    sch5::anndata ad(path);
    ad.set_countsname("nonexistent_layer");

    vector<double> data;
    vector<int32_t> indices;
    vector<int64_t> indptr;

    ASSERT_THROW(ad.load_expr(data, indices, indptr),
        "load_expr with bad countsname should throw");

    remove(path.c_str());
    return true;
}

// =====================================================================
// Tests: csc_to_csr conversion
// =====================================================================

bool test_csc_to_csr(){
    // Build a small 3x3 identity matrix in CSC format
    // CSC: column pointers, row indices
    string path = "/tmp/test_csc.h5ad";
    create_test_anndata(path, 3, 3, true);
    sch5::anndata ad(path);

    // Identity matrix in CSC:
    // col 0: row 0, col 1: row 1, col 2: row 2
    vector<double> data = {1.0, 2.0, 3.0};
    vector<int32_t> indices = {0, 1, 2};   // row indices
    vector<int64_t> indptr = {0, 1, 2, 3}; // col pointers

    ad.csc_to_csr(data, indices, indptr);

    // CSR result should be the same for identity matrix
    ASSERT_EQ((int)indptr.size(), 4, "csr indptr size");
    ASSERT_EQ((int)data.size(), 3, "csr data size");
    // Row 0: col 0, Row 1: col 1, Row 2: col 2
    ASSERT_EQ(indices[0], 0, "csr indices[0]");
    ASSERT_EQ(indices[1], 1, "csr indices[1]");
    ASSERT_EQ(indices[2], 2, "csr indices[2]");
    ASSERT_NEAR(data[0], 1.0, 1e-9, "csr data[0]");
    ASSERT_NEAR(data[1], 2.0, 1e-9, "csr data[1]");
    ASSERT_NEAR(data[2], 3.0, 1e-9, "csr data[2]");

    remove(path.c_str());
    return true;
}

bool test_csc_to_csr_nonsquare(){
    // 4 cells x 3 genes, non-trivial CSC
    string path = "/tmp/test_csc2.h5ad";
    create_test_anndata(path, 4, 3, true);
    sch5::anndata ad(path);

    // CSC format (3 gene columns):
    // Gene 0: cells 0,3 with values 1.0,4.0
    // Gene 1: cell 1 with value 2.0
    // Gene 2: cell 2 with value 3.0
    vector<double> data = {1.0, 4.0, 2.0, 3.0};
    vector<int32_t> indices = {0, 3, 1, 2};
    vector<int64_t> indptr = {0, 2, 3, 4};

    ad.csc_to_csr(data, indices, indptr);

    // CSR: 4 cell rows
    ASSERT_EQ((int)indptr.size(), 5, "csr indptr size");
    // Cell 0: gene 0 (val 1.0)
    ASSERT_EQ(indptr[0], 0, "cell 0 start");
    ASSERT_EQ(indptr[1], 1, "cell 0 end");
    ASSERT_EQ(indices[0], 0, "cell 0 gene");
    ASSERT_NEAR(data[0], 1.0, 1e-9, "cell 0 value");
    // Cell 1: gene 1 (val 2.0)
    ASSERT_EQ(indices[1], 1, "cell 1 gene");
    ASSERT_NEAR(data[1], 2.0, 1e-9, "cell 1 value");
    // Cell 2: gene 2 (val 3.0)
    ASSERT_EQ(indices[2], 2, "cell 2 gene");
    ASSERT_NEAR(data[2], 3.0, 1e-9, "cell 2 value");
    // Cell 3: gene 0 (val 4.0)
    ASSERT_EQ(indices[3], 0, "cell 3 gene");
    ASSERT_NEAR(data[3], 4.0, 1e-9, "cell 3 value");

    remove(path.c_str());
    return true;
}

// =====================================================================
// Tests: mtx2map / map2mtx round-trip
// =====================================================================

bool test_mtx2map_map2mtx(){
    string path = "/tmp/test_mtx2map.h5ad";
    create_test_anndata(path, 3, 3, true);
    sch5::anndata ad(path);

    vector<double> data = {1.0, 2.0, 3.0};
    vector<int32_t> indices = {0, 1, 2};
    vector<int64_t> indptr = {0, 1, 2, 3};

    map<int32_t, map<int32_t, double>> mtxmap;
    ad.mtx2map(data, indices, indptr, mtxmap);

    ASSERT_EQ((int)mtxmap.size(), 3, "map has 3 cells");
    ASSERT_NEAR(mtxmap[0][0], 1.0, 1e-9, "map[0][0]");
    ASSERT_NEAR(mtxmap[1][1], 2.0, 1e-9, "map[1][1]");
    ASSERT_NEAR(mtxmap[2][2], 3.0, 1e-9, "map[2][2]");

    // Round-trip back
    vector<double> data2;
    vector<int32_t> indices2;
    vector<int64_t> indptr2;
    ad.map2mtx(mtxmap, data2, indices2, indptr2);

    ASSERT_EQ((int)data2.size(), 3, "round-trip data size");
    ASSERT_EQ((int)indptr2.size(), 4, "round-trip indptr size");
    ASSERT_NEAR(data2[0], 1.0, 1e-9, "round-trip data[0]");
    ASSERT_NEAR(data2[1], 2.0, 1e-9, "round-trip data[1]");
    ASSERT_NEAR(data2[2], 3.0, 1e-9, "round-trip data[2]");

    remove(path.c_str());
    return true;
}

// =====================================================================
// Tests: loom
// =====================================================================

bool test_loom_cell_gene_names(){
    string path = "/tmp/test_loom_basic.loom";
    create_test_loom(path, 4, 3);
    sch5::loom lm(path);

    ASSERT_EQ(lm.n_cells, 4, "n_cells");
    ASSERT_EQ(lm.n_genes, 3, "n_genes");
    // parse_loom_bc strips "sample:" prefix and trailing "x"
    ASSERT_EQ(lm.cell_names[0], "ACGTACGT", "cell 0 barcode parsed");
    ASSERT_EQ(lm.gene_names[0], "Gene_0", "gene 0");

    remove(path.c_str());
    return true;
}

bool test_loom_load_expr(){
    string path = "/tmp/test_loom_expr.loom";
    create_test_loom(path, 4, 3);
    sch5::loom lm(path);

    vector<double> data;
    vector<int32_t> indices;
    vector<int64_t> indptr;
    lm.load_expr(data, indices, indptr);

    ASSERT_EQ((int)indptr.size(), 5, "indptr size");
    ASSERT_EQ((int)data.size(), 4, "data size");

    // Cell 0 → gene 0, value 1
    ASSERT_EQ(indices[0], 0, "cell 0 gene");
    ASSERT_NEAR(data[0], 1.0, 1e-9, "cell 0 value");

    // Cell 2 → gene 2, value 3
    ASSERT_EQ(indices[2], 2, "cell 2 gene");
    ASSERT_NEAR(data[2], 3.0, 1e-9, "cell 2 value");

    remove(path.c_str());
    return true;
}

bool test_loom_load_expr_map(){
    string path = "/tmp/test_loom_expr_map.loom";
    create_test_loom(path, 4, 3);
    sch5::loom lm(path);

    map<int32_t, map<int32_t, double>> mtx;
    lm.load_expr(mtx);

    ASSERT_EQ((int)mtx.size(), 4, "map has 4 cells");
    ASSERT_NEAR(mtx[0][0], 1.0, 1e-9, "cell 0 gene 0");
    ASSERT_NEAR(mtx[1][1], 2.0, 1e-9, "cell 1 gene 1");

    remove(path.c_str());
    return true;
}

bool test_loom_load_meta_int(){
    string path = "/tmp/test_loom_meta_int.loom";
    create_test_loom(path, 4, 3);
    sch5::loom lm(path);

    vector<int> vals;
    lm.load_meta_col("ClusterID", vals);
    ASSERT_EQ((int)vals.size(), 4, "int col size");
    ASSERT_EQ(vals[0], 0, "ClusterID[0]");
    ASSERT_EQ(vals[1], 1, "ClusterID[1]");
    ASSERT_EQ(vals[2], 2, "ClusterID[2]");

    remove(path.c_str());
    return true;
}

bool test_loom_load_meta_float(){
    string path = "/tmp/test_loom_meta_float.loom";
    create_test_loom(path, 4, 3);
    sch5::loom lm(path);

    vector<double> vals;
    lm.load_meta_col("Score", vals);
    ASSERT_EQ((int)vals.size(), 4, "float col size");
    ASSERT_NEAR(vals[0], 0.0, 1e-9, "Score[0]");
    ASSERT_NEAR(vals[2], 3.0, 1e-9, "Score[2]");

    remove(path.c_str());
    return true;
}

bool test_loom_has_meta_col(){
    string path = "/tmp/test_loom_has_col.loom";
    create_test_loom(path, 3, 2);
    sch5::loom lm(path);

    ASSERT_TRUE(lm.has_meta_col("ClusterID"), "has ClusterID");
    ASSERT_TRUE(lm.has_meta_col("Score"), "has Score");
    ASSERT_TRUE(!lm.has_meta_col("Fake"), "!has Fake");

    remove(path.c_str());
    return true;
}

bool test_loom_meta_col_type(){
    string path = "/tmp/test_loom_coltype.loom";
    create_test_loom(path, 3, 2);
    sch5::loom lm(path);

    ASSERT_EQ(lm.meta_col_type("ClusterID"), sch5::h5_type_int, "ClusterID type");
    ASSERT_EQ(lm.meta_col_type("Score"), sch5::h5_type_float, "Score type");

    remove(path.c_str());
    return true;
}

bool test_loom_list_layers(){
    string path = "/tmp/test_loom_layers.loom";
    create_test_loom(path, 3, 2);
    sch5::loom lm(path);

    vector<string> layers;
    lm.list_layers(layers);
    ASSERT_EQ((int)layers.size(), 1, "one layer (spliced)");
    ASSERT_EQ(layers[0], "spliced", "layer name");

    remove(path.c_str());
    return true;
}

bool test_loom_load_layer(){
    string path = "/tmp/test_loom_loadlayer.loom";
    create_test_loom(path, 4, 3);
    sch5::loom lm(path);

    vector<double> data;
    vector<int32_t> indices;
    vector<int64_t> indptr;
    lm.load_layer("spliced", data, indices, indptr);

    ASSERT_EQ((int)indptr.size(), 5, "layer indptr size");
    ASSERT_EQ((int)data.size(), 4, "layer data size");

    // Cell 0 → gene 1 ((0+1)%3=1), value 10
    ASSERT_EQ(indices[0], 1, "layer cell 0 gene");
    ASSERT_NEAR(data[0], 10.0, 1e-9, "layer cell 0 value");

    remove(path.c_str());
    return true;
}

bool test_loom_get_meta_colnames(){
    string path = "/tmp/test_loom_colnames.loom";
    create_test_loom(path, 3, 2);
    sch5::loom lm(path);

    vector<string> str_cols, int_cols, float_cols;
    lm.get_meta_colnames(str_cols, int_cols, float_cols);

    ASSERT_EQ((int)int_cols.size(), 1, "int cols count");
    ASSERT_EQ(int_cols[0], "ClusterID", "int col name");
    ASSERT_EQ((int)float_cols.size(), 1, "float cols count");
    ASSERT_EQ(float_cols[0], "Score", "float col name");

    remove(path.c_str());
    return true;
}

// =====================================================================
// Tests: anndata write_layer + load_layer round-trip
// =====================================================================

bool test_anndata_write_load_layer(){
    string path = "/tmp/test_anndata_wl.h5ad";
    create_test_anndata(path, 4, 3, true);

    // Reopen with write access
    sch5::anndata ad(path);
    ad.open(path);

    // Write a layer
    vector<double> data = {10.0, 20.0, 30.0, 40.0};
    vector<int32_t> indices = {0, 1, 2, 0};
    vector<int64_t> indptr = {0, 1, 2, 3, 4};
    ad.write_layer("test_layer", data, indices, indptr);

    // Read it back
    vector<double> data2;
    vector<int32_t> indices2;
    vector<int64_t> indptr2;
    ad.load_layer("test_layer", data2, indices2, indptr2);

    ASSERT_EQ((int)data2.size(), 4, "layer data size");
    ASSERT_NEAR(data2[0], 10.0, 1e-9, "layer data[0]");
    ASSERT_NEAR(data2[3], 40.0, 1e-9, "layer data[3]");

    // list_layers should now include it
    vector<string> layers;
    ad.list_layers(layers);
    ASSERT_TRUE(find(layers.begin(), layers.end(), "test_layer") != layers.end(),
        "test_layer in layer list");

    remove(path.c_str());
    return true;
}

bool test_anndata_has_layer(){
    string path = "/tmp/test_anndata_haslayer.h5ad";
    create_test_anndata(path, 3, 2, true);
    sch5::anndata ad(path);

    ASSERT_TRUE(!ad.has_layer("nope"), "no layer exists yet");

    remove(path.c_str());
    return true;
}

// =====================================================================
// Helper: create a minimal h5Seurat file
// =====================================================================
static string create_test_seurat(const string& path,
    int n_cells, int n_genes, bool sparse = true){

    {
        HighFive::File f(path, HighFive::File::Truncate);

        // Active assay
        f.createDataSet("active.assay", string("RNA"));

        // Cell names
        vector<string> barcodes;
        for (int i = 0; i < n_cells; ++i){
            barcodes.push_back("CELL_" + to_string(i));
        }
        f.createDataSet("cell.names", barcodes);

        // --- meta.data ---
        auto meta = f.createGroup("meta.data");

        // int column
        vector<int32_t> ncount;
        for (int i = 0; i < n_cells; ++i) ncount.push_back(i * 100);
        meta.createDataSet("nCount_RNA", ncount);

        // float column
        vector<double> pct;
        for (int i = 0; i < n_cells; ++i) pct.push_back(i * 0.25);
        meta.createDataSet("percent.mt", pct);

        // string column
        vector<string> orig;
        for (int i = 0; i < n_cells; ++i){
            orig.push_back(i % 2 == 0 ? "sampleA" : "sampleB");
        }
        meta.createDataSet("orig.ident", orig);

        // factor column (integer codes + "levels" attribute)
        vector<string> levels = {"G1", "G2M", "S"};
        vector<int32_t> codes;
        for (int i = 0; i < n_cells; ++i) codes.push_back(i % 3);
        auto phase_ds = meta.createDataSet("Phase", codes);
        phase_ds.createAttribute("levels", levels);

        // --- assays/RNA ---
        auto assays = f.createGroup("assays");
        auto rna = assays.createGroup("RNA");

        // Gene names
        vector<string> genes;
        for (int i = 0; i < n_genes; ++i){
            genes.push_back("Gene_" + to_string(i));
        }
        rna.createDataSet("features", genes);
        rna.createDataSet("key", string("rna_"));

        if (sparse){
            // Sparse counts (CSR / equivalent to h5Seurat CSC of genes x cells).
            // Each cell c has one nonzero at gene (c % n_genes) with value c+1.
            auto counts_grp = rna.createGroup("counts");
            vector<double> data;
            vector<int32_t> indices;
            vector<int64_t> indptr;

            for (int c = 0; c < n_cells; ++c){
                indptr.push_back((int64_t)data.size());
                int g = c % n_genes;
                data.push_back((double)(c + 1));
                indices.push_back((int32_t)g);
            }
            indptr.push_back((int64_t)data.size());

            counts_grp.createDataSet("data", data);
            counts_grp.createDataSet("indices", indices);
            counts_grp.createDataSet("indptr", indptr);
            counts_grp.createAttribute<string>("encoding-type",
                string("csr_matrix"));
            counts_grp.createAttribute<string>("encoding-version",
                string("0.1.0"));
            vector<int64_t> shape = {(int64_t)n_cells, (int64_t)n_genes};
            counts_grp.createAttribute("shape", shape);

            // A "data" slot (normalized) with doubled values
            auto data_grp = rna.createGroup("data");
            vector<double> norm_data;
            vector<int32_t> norm_indices;
            vector<int64_t> norm_indptr;
            for (int c = 0; c < n_cells; ++c){
                norm_indptr.push_back((int64_t)norm_data.size());
                int g = c % n_genes;
                norm_data.push_back((double)(c + 1) * 2.0);
                norm_indices.push_back((int32_t)g);
            }
            norm_indptr.push_back((int64_t)norm_data.size());

            data_grp.createDataSet("data", norm_data);
            data_grp.createDataSet("indices", norm_indices);
            data_grp.createDataSet("indptr", norm_indptr);
            data_grp.createAttribute<string>("encoding-type",
                string("csr_matrix"));
            data_grp.createAttribute<string>("encoding-version",
                string("0.1.0"));
            data_grp.createAttribute("shape", shape);
        }
        else{
            // Dense counts: cells x genes, row-major
            vector<double> dense(n_cells * n_genes, 0.0);
            for (int c = 0; c < n_cells; ++c){
                int g = c % n_genes;
                dense[c * n_genes + g] = (double)(c + 1);
            }
            vector<size_t> dims = {(size_t)n_cells, (size_t)n_genes};
            auto ds = rna.createDataSet<double>("counts",
                HighFive::DataSpace(dims));
            ds.write_raw(dense.data());
            ds.createAttribute<string>("encoding-type", string("array"));
            ds.createAttribute<string>("encoding-version", string("0.2.0"));
        }
    }
    return path;
}

// =====================================================================
// Tests: seurat
// =====================================================================

bool test_seurat_cell_gene_names(){
    string path = "/tmp/test_seurat_basic.h5seurat";
    create_test_seurat(path, 5, 3);
    sch5::seurat sr(path);

    ASSERT_EQ(sr.n_cells, 5, "n_cells");
    ASSERT_EQ(sr.n_genes, 3, "n_genes");
    ASSERT_EQ(sr.cell_names[0], "CELL_0", "cell_names[0]");
    ASSERT_EQ(sr.cell_names[4], "CELL_4", "cell_names[4]");
    ASSERT_EQ(sr.gene_names[0], "Gene_0", "gene_names[0]");
    ASSERT_EQ(sr.gene_names[2], "Gene_2", "gene_names[2]");
    remove(path.c_str());
    return true;
}

bool test_seurat_load_expr_sparse(){
    string path = "/tmp/test_seurat_expr.h5seurat";
    create_test_seurat(path, 5, 3, true);
    sch5::seurat sr(path);

    vector<double> data;
    vector<int32_t> indices;
    vector<int64_t> indptr;
    sr.load_expr(data, indices, indptr);

    ASSERT_EQ((int)indptr.size(), 6, "indptr size (n_cells+1)");
    ASSERT_EQ((int)data.size(), 5, "data size");

    // Cell 0 → gene 0, value 1
    ASSERT_EQ(indices[0], 0, "cell 0 gene index");
    ASSERT_NEAR(data[0], 1.0, 1e-9, "cell 0 value");

    // Cell 3 → gene 0 (3%3==0), value 4
    ASSERT_EQ(indices[3], 0, "cell 3 gene index");
    ASSERT_NEAR(data[3], 4.0, 1e-9, "cell 3 value");

    remove(path.c_str());
    return true;
}

bool test_seurat_load_expr_dense(){
    string path = "/tmp/test_seurat_dense.h5seurat";
    create_test_seurat(path, 4, 3, false);
    sch5::seurat sr(path);

    vector<double> data;
    vector<int32_t> indices;
    vector<int64_t> indptr;
    sr.load_expr(data, indices, indptr);

    ASSERT_EQ((int)indptr.size(), 5, "indptr size");
    ASSERT_EQ((int)data.size(), 4, "data size");
    ASSERT_EQ(indices[0], 0, "cell 0 gene");
    ASSERT_NEAR(data[0], 1.0, 1e-9, "cell 0 value");

    remove(path.c_str());
    return true;
}

bool test_seurat_load_meta_int(){
    string path = "/tmp/test_seurat_meta_int.h5seurat";
    create_test_seurat(path, 4, 2);
    sch5::seurat sr(path);

    vector<int> vals;
    sr.load_meta_col("nCount_RNA", vals);
    ASSERT_EQ((int)vals.size(), 4, "int col size");
    ASSERT_EQ(vals[0], 0, "nCount_RNA[0]");
    ASSERT_EQ(vals[2], 200, "nCount_RNA[2]");

    remove(path.c_str());
    return true;
}

bool test_seurat_load_meta_float(){
    string path = "/tmp/test_seurat_meta_float.h5seurat";
    create_test_seurat(path, 4, 2);
    sch5::seurat sr(path);

    vector<double> vals;
    sr.load_meta_col("percent.mt", vals);
    ASSERT_EQ((int)vals.size(), 4, "float col size");
    ASSERT_NEAR(vals[0], 0.0, 1e-9, "percent.mt[0]");
    ASSERT_NEAR(vals[3], 0.75, 1e-9, "percent.mt[3]");

    remove(path.c_str());
    return true;
}

bool test_seurat_load_meta_string(){
    string path = "/tmp/test_seurat_meta_str.h5seurat";
    create_test_seurat(path, 4, 2);
    sch5::seurat sr(path);

    vector<string> vals;
    sr.load_meta_col("orig.ident", vals);
    ASSERT_EQ((int)vals.size(), 4, "str col size");
    ASSERT_EQ(vals[0], "sampleA", "orig.ident[0]");
    ASSERT_EQ(vals[1], "sampleB", "orig.ident[1]");

    remove(path.c_str());
    return true;
}

bool test_seurat_load_meta_factor(){
    string path = "/tmp/test_seurat_meta_fac.h5seurat";
    create_test_seurat(path, 6, 2);
    sch5::seurat sr(path);

    vector<string> vals;
    sr.load_meta_col("Phase", vals);
    ASSERT_EQ((int)vals.size(), 6, "factor col size");
    ASSERT_EQ(vals[0], "G1", "Phase[0]");
    ASSERT_EQ(vals[1], "G2M", "Phase[1]");
    ASSERT_EQ(vals[2], "S", "Phase[2]");
    ASSERT_EQ(vals[3], "G1", "Phase[3]");

    remove(path.c_str());
    return true;
}

bool test_seurat_has_meta_col(){
    string path = "/tmp/test_seurat_has_col.h5seurat";
    create_test_seurat(path, 3, 2);
    sch5::seurat sr(path);

    ASSERT_TRUE(sr.has_meta_col("nCount_RNA"), "has nCount_RNA");
    ASSERT_TRUE(sr.has_meta_col("percent.mt"), "has percent.mt");
    ASSERT_TRUE(sr.has_meta_col("orig.ident"), "has orig.ident");
    ASSERT_TRUE(sr.has_meta_col("Phase"), "has Phase");
    ASSERT_TRUE(!sr.has_meta_col("nonexistent"), "!has nonexistent");

    remove(path.c_str());
    return true;
}

bool test_seurat_meta_col_type(){
    string path = "/tmp/test_seurat_coltype.h5seurat";
    create_test_seurat(path, 3, 2);
    sch5::seurat sr(path);

    ASSERT_EQ(sr.meta_col_type("nCount_RNA"), sch5::h5_type_int, "nCount type");
    ASSERT_EQ(sr.meta_col_type("percent.mt"), sch5::h5_type_float, "pct type");
    ASSERT_EQ(sr.meta_col_type("orig.ident"), sch5::h5_type_str, "str type");
    ASSERT_EQ(sr.meta_col_type("Phase"), sch5::h5_type_str, "factor type");

    remove(path.c_str());
    return true;
}

bool test_seurat_get_meta_colnames(){
    string path = "/tmp/test_seurat_colnames.h5seurat";
    create_test_seurat(path, 3, 2);
    sch5::seurat sr(path);

    vector<string> str_cols, int_cols, float_cols;
    sr.get_meta_colnames(str_cols, int_cols, float_cols);

    ASSERT_EQ((int)int_cols.size(), 1, "int cols count");
    ASSERT_EQ(int_cols[0], "nCount_RNA", "int col name");
    ASSERT_EQ((int)float_cols.size(), 1, "float cols count");
    ASSERT_EQ(float_cols[0], "percent.mt", "float col name");
    // str_cols should have orig.ident + Phase (factor)
    ASSERT_TRUE((int)str_cols.size() >= 2, "str cols >= 2");

    remove(path.c_str());
    return true;
}

bool test_seurat_list_layers(){
    string path = "/tmp/test_seurat_layers.h5seurat";
    create_test_seurat(path, 3, 2, true);
    sch5::seurat sr(path);

    vector<string> layers;
    sr.list_layers(layers);
    ASSERT_TRUE((int)layers.size() >= 1, "at least counts layer");
    ASSERT_TRUE(find(layers.begin(), layers.end(), "counts") != layers.end(),
        "has counts");
    ASSERT_TRUE(find(layers.begin(), layers.end(), "data") != layers.end(),
        "has data");

    remove(path.c_str());
    return true;
}

bool test_seurat_has_layer(){
    string path = "/tmp/test_seurat_haslayer.h5seurat";
    create_test_seurat(path, 3, 2, true);
    sch5::seurat sr(path);

    ASSERT_TRUE(sr.has_layer("counts"), "has counts");
    ASSERT_TRUE(sr.has_layer("data"), "has data");
    ASSERT_TRUE(!sr.has_layer("scale.data"), "!has scale.data");

    remove(path.c_str());
    return true;
}

bool test_seurat_load_layer(){
    string path = "/tmp/test_seurat_loadlayer.h5seurat";
    create_test_seurat(path, 4, 3, true);
    sch5::seurat sr(path);

    vector<double> data;
    vector<int32_t> indices;
    vector<int64_t> indptr;
    sr.load_layer("data", data, indices, indptr);

    ASSERT_EQ((int)indptr.size(), 5, "layer indptr size");
    ASSERT_EQ((int)data.size(), 4, "layer data size");

    // Normalized "data" slot has doubled values: cell 0 → value 2.0
    ASSERT_NEAR(data[0], 2.0, 1e-9, "layer cell 0 value");
    ASSERT_NEAR(data[3], 8.0, 1e-9, "layer cell 3 value");

    remove(path.c_str());
    return true;
}

bool test_seurat_load_expr_no_encoding(){
    // Simulate a real SeuratDisk file: sparse group without encoding-type
    string path = "/tmp/test_seurat_noenc.h5seurat";
    {
        HighFive::File f(path, HighFive::File::Truncate);
        f.createDataSet("active.assay", string("RNA"));

        vector<string> barcodes = {"C0", "C1", "C2"};
        f.createDataSet("cell.names", barcodes);

        auto assays = f.createGroup("assays");
        auto rna = assays.createGroup("RNA");
        vector<string> genes = {"G0", "G1"};
        rna.createDataSet("features", genes);

        // Sparse counts WITHOUT encoding-type (dgCMatrix style).
        // CSC of (2 genes x 3 cells) = CSR of (3 cells x 2 genes).
        // Cell 0→gene 0 (val 1), Cell 1→gene 1 (val 2), Cell 2→gene 0 (val 3)
        auto grp = rna.createGroup("counts");
        vector<double> data = {1.0, 2.0, 3.0};
        vector<int32_t> indices = {0, 1, 0};
        vector<int64_t> indptr = {0, 1, 2, 3};
        grp.createDataSet("data", data);
        grp.createDataSet("indices", indices);
        grp.createDataSet("indptr", indptr);
        // No encoding-type attribute — just like SeuratDisk
    }

    sch5::seurat sr(path);
    vector<double> data;
    vector<int32_t> indices;
    vector<int64_t> indptr;
    sr.load_expr(data, indices, indptr);

    ASSERT_EQ((int)indptr.size(), 4, "indptr size");
    ASSERT_EQ((int)data.size(), 3, "data size");
    ASSERT_EQ(indices[0], 0, "cell 0 gene");
    ASSERT_NEAR(data[0], 1.0, 1e-9, "cell 0 value");
    ASSERT_EQ(indices[1], 1, "cell 1 gene");
    ASSERT_NEAR(data[1], 2.0, 1e-9, "cell 1 value");
    ASSERT_EQ(indices[2], 0, "cell 2 gene");
    ASSERT_NEAR(data[2], 3.0, 1e-9, "cell 2 value");

    remove(path.c_str());
    return true;
}

bool test_seurat_load_expr_not_found(){
    string path = "/tmp/test_seurat_notfound.h5seurat";
    create_test_seurat(path, 3, 2);
    sch5::seurat sr(path);
    sr.set_countsname("nonexistent_slot");

    vector<double> data;
    vector<int32_t> indices;
    vector<int64_t> indptr;

    ASSERT_THROW(sr.load_expr(data, indices, indptr),
        "load_expr with bad countsname should throw");

    remove(path.c_str());
    return true;
}

bool test_seurat_write_load_layer(){
    string path = "/tmp/test_seurat_wl.h5seurat";
    create_test_seurat(path, 4, 3);

    sch5::seurat sr(path);
    sr.open(path);

    vector<double> data = {10.0, 20.0, 30.0, 40.0};
    vector<int32_t> indices = {0, 1, 2, 0};
    vector<int64_t> indptr = {0, 1, 2, 3, 4};
    sr.write_layer("scale.data", data, indices, indptr);

    vector<double> data2;
    vector<int32_t> indices2;
    vector<int64_t> indptr2;
    sr.load_layer("scale.data", data2, indices2, indptr2);

    ASSERT_EQ((int)data2.size(), 4, "written layer data size");
    ASSERT_NEAR(data2[0], 10.0, 1e-9, "written layer data[0]");
    ASSERT_NEAR(data2[3], 40.0, 1e-9, "written layer data[3]");

    remove(path.c_str());
    return true;
}

// =====================================================================
// main
// =====================================================================

int main(){
    cout << "=== h5_reader / anndata / loom unit tests ===" << endl;

    cout << endl << "--- anndata ---" << endl;
    RUN_TEST(test_anndata_cell_gene_names);
    RUN_TEST(test_anndata_load_expr_sparse_csr);
    RUN_TEST(test_anndata_load_expr_sparse_csc);
    RUN_TEST(test_anndata_load_expr_dense);
    RUN_TEST(test_anndata_load_meta_int);
    RUN_TEST(test_anndata_load_meta_float);
    RUN_TEST(test_anndata_load_meta_string);
    RUN_TEST(test_anndata_load_meta_categorical);
    RUN_TEST(test_anndata_has_meta_col);
    RUN_TEST(test_anndata_meta_col_type);
    RUN_TEST(test_anndata_get_meta_colnames);
    RUN_TEST(test_anndata_list_layers_empty);
    RUN_TEST(test_anndata_load_expr_not_found);
    RUN_TEST(test_anndata_write_load_layer);
    RUN_TEST(test_anndata_has_layer);

    cout << endl << "--- csc_to_csr ---" << endl;
    RUN_TEST(test_csc_to_csr);
    RUN_TEST(test_csc_to_csr_nonsquare);

    cout << endl << "--- mtx2map / map2mtx ---" << endl;
    RUN_TEST(test_mtx2map_map2mtx);

    cout << endl << "--- loom ---" << endl;
    RUN_TEST(test_loom_cell_gene_names);
    RUN_TEST(test_loom_load_expr);
    RUN_TEST(test_loom_load_expr_map);
    RUN_TEST(test_loom_load_meta_int);
    RUN_TEST(test_loom_load_meta_float);
    RUN_TEST(test_loom_has_meta_col);
    RUN_TEST(test_loom_meta_col_type);
    RUN_TEST(test_loom_list_layers);
    RUN_TEST(test_loom_load_layer);
    RUN_TEST(test_loom_get_meta_colnames);

    cout << endl << "--- seurat ---" << endl;
    RUN_TEST(test_seurat_cell_gene_names);
    RUN_TEST(test_seurat_load_expr_sparse);
    RUN_TEST(test_seurat_load_expr_dense);
    RUN_TEST(test_seurat_load_meta_int);
    RUN_TEST(test_seurat_load_meta_float);
    RUN_TEST(test_seurat_load_meta_string);
    RUN_TEST(test_seurat_load_meta_factor);
    RUN_TEST(test_seurat_has_meta_col);
    RUN_TEST(test_seurat_meta_col_type);
    RUN_TEST(test_seurat_get_meta_colnames);
    RUN_TEST(test_seurat_list_layers);
    RUN_TEST(test_seurat_has_layer);
    RUN_TEST(test_seurat_load_layer);
    RUN_TEST(test_seurat_load_expr_no_encoding);
    RUN_TEST(test_seurat_load_expr_not_found);
    RUN_TEST(test_seurat_write_load_layer);

    cout << endl << "=== Results: " << tests_passed << "/" << tests_run << " passed ===" << endl;
    return (tests_passed == tests_run) ? 0 : 1;
}
