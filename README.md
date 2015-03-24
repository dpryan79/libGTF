libGTF a library for parsing/searching GTF and BED files
========================================================

Note that this library is not yet appropriate for real-world use (the documentation hasn't even been written yet!).

In essense, this library construct an interval tree representation of GTF or BED files and allows overlap searches, similar to GenomicRanges or bedtools. The primary difference here is the convenient C interface that allows incorporation into other programs easily.

Note that murmur3.c and murmur3.h are C implementations of MurmurHash. The C implementation is from [Peter Scott](https://github.com/PeterScott/murmur3) and MurmurHash itself is by [Austin Appleby](https://code.google.com/p/smhasher/wiki/MurmurHash3). Both of these are in the public domain.
