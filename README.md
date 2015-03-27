libGTF a library for parsing/searching GTF and BED files
========================================================

Note that this library is not yet appropriate for real-world use (the documentation hasn't even been written yet!).

In essense, this library construct an interval tree representation of GTF or BED files and allows overlap searches, similar to GenomicRanges or bedtools. The primary difference here is the convenient C interface that allows incorporation into other programs easily.

Note that murmur3.c and murmur3.h are C implementations of MurmurHash. The C implementation is from [Peter Scott](https://github.com/PeterScott/murmur3) and MurmurHash itself is by [Austin Appleby](https://code.google.com/p/smhasher/wiki/MurmurHash3). Both of these are in the public domain.

Examples
========

There are currently a few example programs in the `tests/` directory.

  * `testBED` demonstrates how to parse a BED file and display its tree representation in dot format.
  * `testGTF` is the equivalent program for GTF files.
  * `testFindOverlaps` demonstrates how to find overlaps of alignments in SAM/BAM/CRAM format with a GTF file. This is largely similar to featureCounts and htseq-count, though this program doesn't handle paired-end alignments intelligently (it is, afterall, just a demonstration). The code demonstrates processing sets of overlaps and merging them for further processing. This program also demonstrates how to use different match and strand types. Note that in a real program, it'd be simpler to use the `findOverlapsBAM()` function.

To Do
=====

 - [X] cntGeneIDs()
 - [ ] findOverlaps, but only if feature is XXX or not YYY (allow lists)
 - [ ] uniqueGeneIDs()
 - [ ] cntTranscriptIDs()
 - [ ] uniqueTranscriptIDs()
 - [ ] cntAttributes()
 - [ ] uniqueAttributes()
 - [ ] test various overlap options
   - [X] findOverlaps
   - [ ] countOverlaps
   - [ ] overlapsAny
 - [ ] test strandedness with above
 - [ ] test match type options with above
 - [ ] compare with featureCounts/bedtools and ensure the results are the same
 - [ ] add a findClosest() function
 - [ ] parseGFF?
 - [ ] test for GFF?
 - [ ] intersect/union/diff on overlapSets
