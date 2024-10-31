# htswrapper
Collection of handy things for code dealing with high throughput sequencing data

## Details/installation

### Dependencies
Requires [HTSLib](https://github.com/samtools/htslib) >= 1.10.2 and [zlib](https://www.zlib.net/). Also uses [robin_hood](https://github.com/martinus/robin-hood-hashing) hashing, which is included in this repository.

### Barcode length

Because barcodes are represented as bitsets, their length (and the bitset width) has to be set at compile time. This is done through two variables present in the Makefile: `BC_LENX2`, which is two times the barcode length, and `KX2`, which is two times the k-mer length for fuzzy barcode matching.

Default barcode length is set to the current 10X Genomics standard (and maximum possible value for interpreting a bitset an integer and thus not breaking all of this code): 16. For barcodes of length $L$, K-mers for fuzzy matching should be set at their maximum possible length for best performance, which is $\lfloor\frac{L+1}{2}\rfloor$. In the default case of 16-base barcodes, this is 8. If you want to change these values, you need to tell the compiler to make the corresponding bitsets twice as wide. 

To do this, specify `BC_LENX2=[your value]` and `KX2=[your value]` when running make. If you want to use 14-base barcodes instead of 16-base barcodes, for example, run

`make BC_LENX2=28 KX2=14`

### bam.cpp
A class that wraps HTSLib's BAM reader to make it easier to use and remember how to access stuff

### bc.cpp 
Functionality to make it easier to deal with cell barcodes, bit-packed and interpreted as `unsigned long`. This makes it fast to store, look up, and compare cell barcodes. 
*  Contains functions to interpret strings as barcodes and vice versa, in both forward and reverse complement orientation
*  Contains a class designed to represent allowed barcode lists, or even paired lists (i.e. for 10X Genomics multiome data, where there is one list for RNA-seq and another for ATAC-seq, and each ATAC-seq barcode corresponds to an RNA-seq barcode, which is what ends up in the BAM file
    * When a single `N` is encountered in barcode matching, attempts to mutate it to `A`, `C`, `G`, and `T` and checks if any of these correspond to valid barcodes. If only one does, that barcode is returned.
    * Also stores shorter (length $\lfloor\frac{L+1}{2}\rfloor$) k-mers that make up cell barcodes. If a sequence does not exactly match any known barcode and contains no Ns, then these k-mers are used in a fuzzy matching routine that will return a valid barcode that is edit distance 1 from the sequence, if only one such barcode exists.
    * Barcode lookup is fast: relies on [robin_hood](https://github.com/martinus/robin-hood-hashing) hash maps and hash sets instead of `std::unordered_map` and `std::unordered_set`
    * k-mer lookup is fast: k-mers are directly interpreted as indices into an array of linked lists of cell barcodes

### bc_scanner.cpp
Builds on `bc.cpp` to provide a class that iterates through a set of FASTQ files and finds (fuzzy) matches to barcodes in an allowed list. Provides access to each matched read's ID, sequence, and barcode via the `next()` method, and provides the option to trim barcodes from sequences.

### gzreader.cpp
A convenience class that makes it easy to read through a file line by line, whether gzipped or not. Automatically detects whether the file is gzipped and handles buffering.

### seq_fuzzy_match.cpp
A class that takes a list of reference sequences (i.e. barcodes of some type) and provides a function to match any sequence to one of the references, returning the index of the match and making the edit distance accessible. Can be set to return best match overall or best match within a set edit distance, and can either count Ns as mismatches or as a match to any base. Reports no match in the case of a tie. Uses the [edlib](https://github.com/Martinsos/edlib) library to do fuzzy matching.

### serialize.cpp
A class that makes it quick/easy to read and write binary data.

### umi.cpp
A class that makes it easy to de-duplicate UMIs. Stores bit-packed representations of UMI sequences. Can quickly find exact matches and uses k-mers to find matches with an edit distance of 1 if this fails. Provides a function to count unique UMIs after collapsing.
