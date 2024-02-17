# htswrapper
Collection of handy things for code dealing with high throughput sequencing data

### bam.cpp
A class that wraps HTSLib's BAM reader to make it easier to use and remember how to access stuff

### bc.cpp 
Functionality to make it easier to deal with cell barcodes, bit-packed and interpreted as `unsigned long`. This makes it fast to store, look up, and compare cell barcodes. 
*  Contains functions to interpret strings as barcodes and vice versa, in both forward and reverse complement orientation
*  Contains a class designed to represent allowed barcode lists, or even paired lists (i.e. for 10X Genomics multiome data, where there is one list for RNA-seq and another for ATAC-seq, and each ATAC-seq barcode corresponds to an RNA-seq barcode, which is what ends up in the BAM file
    * When a single `N` is encountered in barcode matching, attempts to mutate it to `A`, `C`, `G`, and `T` and checks if any of these correspond to valid barcodes. If only one does, that barcode is returned.
    * Also stores shorter (length $\lfloor\frac{L+1}{2}\rfloor$) k-mers that make up cell barcodes. If a sequence does not exactly match any known barcode and contains no Ns, then these k-mers are used in a fuzzy matching routine that will return a valid barcode that is edit distance 1 from the sequence, if only one such barcode exists.
    * Barcode lookup is fast: relies on [robin_hood hash maps](https://github.com/martinus/robin-hood-hashing) instead of `std::unordered_map` and `std::unordered_set`
    * k-mer lookup is fast: k-mers are directly interpreted as indices into an array of linked lists of cell barcodes

### gzreader.cpp
A convenience class that makes it quick/easy to read through a file line by line, whether gzipped or not. Automatically detects whether the file is gzipped and handles buffering.

### serialize.cpp
A convenience class that makes it quick/easy to read and write binary data.
