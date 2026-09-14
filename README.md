# htswrapper
Collection of handy things for code dealing with high throughput sequencing data

## Details/installation

### Dependencies
Requires [HTSLib](https://github.com/samtools/htslib) >= 1.10.2 and [zlib](https://www.zlib.net/). Also uses [robin_hood](https://github.com/martinus/robin-hood-hashing) hashing, which is included in this repository.

For hdf5 file reading support, requires libhdf5. If this is not present, the relevant programs will not be compiled.

### Barcode length

Because barcodes are represented as bitsets, their length (and the bitset width) has to be set at compile time. This is done through two variables present in the Makefile: `BC_LENX2`, which is two times the barcode length, and `KX2`, which is two times the k-mer length for fuzzy barcode matching.

Default barcode length is set to the current 10X Genomics standard (and maximum possible value for interpreting a bitset an integer and thus not breaking all of this code): 16. For barcodes of length $L$, K-mers for fuzzy matching should be set at their maximum possible length for best performance, which is $\lfloor\frac{L+1}{2}\rfloor$. In the default case of 16-base barcodes, this is 8. If you want to change these values, you need to tell the compiler to make the corresponding bitsets twice as wide. 

To do this, specify `BC_LENX2=[your value]` and `KX2=[your value]` when running make. If you want to use 14-base barcodes instead of 16-base barcodes, for example, run

`make BC_LENX2=28 KX2=14`

## Features
### Reading BAM files
#### bam.cpp / bam.h
A class that wraps HTSLib's BAM reader to make it easier to use and remember how to access stuff

### Hashing cell barcodes
#### bc.cpp / bc.h
Functionality to make it easier to deal with cell barcodes, bit-packed and interpreted as `unsigned long`. This makes it fast to store, look up, and compare cell barcodes. 
*  Contains functions to interpret strings as barcodes and vice versa, in both forward and reverse complement orientation
*  Contains a class designed to represent allowed barcode lists, or even paired lists (i.e. for 10X Genomics multiome data, where there is one list for RNA-seq and another for ATAC-seq, and each ATAC-seq barcode corresponds to an RNA-seq barcode, which is what ends up in the BAM file
    * When a single `N` is encountered in barcode matching, attempts to mutate it to `A`, `C`, `G`, and `T` and checks if any of these correspond to valid barcodes. If only one does, that barcode is returned.
    * Also stores shorter (length $\lfloor\frac{L+1}{2}\rfloor$) k-mers that make up cell barcodes. If a sequence does not exactly match any known barcode and contains no Ns, then these k-mers are used in a fuzzy matching routine that will return a valid barcode that is edit distance 1 from the sequence, if only one such barcode exists.
    * Barcode lookup is fast: relies on [robin_hood](https://github.com/martinus/robin-hood-hashing) hash maps and hash sets instead of `std::unordered_map` and `std::unordered_set`
    * k-mer lookup is fast: k-mers are directly interpreted as indices into an array of linked lists of cell barcodes

#### bc_scanner.cpp
Builds on `bc.cpp` to provide a class that iterates through a set of FASTQ files and finds (fuzzy) matches to barcodes in an allowed list. Provides access to each matched read's ID, sequence, and barcode via the `next()` method, and provides the option to trim barcodes from sequences.

### Reading gzipped files
#### gzreader.cpp / gzreader.h
A convenience class that makes it easy to read through a file line by line, whether gzipped or not. Automatically detects whether the file is gzipped and handles buffering. Can optionally split each line into fields using a fixed delimiter.

### Counting / hashing k-mers
#### khashtable.cpp / khashtable.h
A data structure for fast lookup of arbitrary-length k-mers. This uses the idea of representing k-mers as bit strings (A = 00, C = 01, G = 10, T = 11), similarly to bc.cpp and umi.cpp. Bit strings then allow quick lookup: they can be interpreted as an integer. This is limited by the size (in bytes) of the integer variable. This data structure gets around this by using a chain of uint64_t. This lets us pack 32 ACGT characters into each uint64_t, so kmers up to length 32 are stored in a single integer, and longer k-mers add a second integer for every additional set of bases up to 32. This currently seems to be the most efficient data structure for storing things mapped to k-mers.

#### kmsuftree.c / kmsuftree.h
This is another k-mer lookup data structure based on the same principle. This data structure is implemented as a suffix tree in C, where the first set of bases are used to build an integer that provides the starting point via an array lookup (the integer is the array index). Each subsequent base then points to the next tree node. This is a modification of a file/data structure originally written by Ed Green. 

### Fuzzy matching for strings
#### seq_fuzzy_match.cpp / seq_fuzzy_match.h
A class that takes a list of reference sequences (i.e. barcodes of some type) and provides a function to match any sequence to one of the references, returning the index of the match and making the edit distance accessible. Can be set to return best match overall or best match within a set edit distance, and can either count Ns as mismatches or as a match to any base. Reports no match in the case of a tie. Uses the [edlib](https://github.com/Martinsos/edlib) library to do fuzzy matching.

### Serializing data in binary format
#### serialize.cpp / serialize.h
A class that makes it quick/easy to read and write binary data. Not much is yet implemented here.

### UMI deduplication
#### umi.cpp
A class that makes it easy to de-duplicate UMIs. Stores bit-packed representations of UMI sequences. Can quickly find exact matches and uses k-mers to find matches with an edit distance of 1 if this fails. Provides a function to count unique UMIs after collapsing.

### Reading/modifying hdf5 files for single-cell data
hdf5 file support requires libhdf5 (without it, these programs will not be compiled). hdf5 reading/writing is done via the included [HighFive](https://github.com/highfive-devs/highfive/) library.

#### h5_reader.cpp / h5_reader.h
Base class providing function templates

#### anndata.cpp / anndata.h
Class for manipulating AnnData h5 files

#### loom.cpp / loom.h
Class for manipulating loom files (hdf5 for storing spliced/unspliced read counts)

#### seurat.cpp / seurat.h
Class for manipulating Seurat h5 files 
