#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include "../gtf.h"
#include "../htslib/htslib/sam.h"

int onlyPseudo(GTFtree *t, GTFentry *e) {
    char *biotype = getAttribute(t, e, "gene_biotype");
    if(!biotype) return 0;
    if(strcmp(biotype, "processed_pseudogene") == 0) return 1;
    return 0;
}
int noPseudo(GTFtree *t, GTFentry *e) {
    char *biotype = getAttribute(t, e, "gene_biotype");
    if(!biotype) return 1;
    if(strcmp(biotype, "processed_pseudogene") == 0) return 0;
    return 1;
}

void countOverlapsBAM2(GTFtree *t, htsFile *fp, bam_hdr_t *hdr, int matchType, int strandType, int minMapq, FILTER_ENTRY_FUNC ffunc) {
    bam1_t *b = bam_init1();
    int32_t cnt = 0;
    kstring_t *ks = calloc(1, sizeof(kstring_t));
    int32_t i, start, end;
    uint32_t *CIGAR, op;
    assert(b);
    assert(ks);

    while(sam_read1(fp, hdr, b) >= 0) {
        if(b->core.tid < 0) continue;
        if(b->core.flag & BAM_FUNMAP) continue;
        if(b->core.qual < minMapq) continue;

        //Iterate over mapped segments, finding the count of each mapped segment
        start = b->core.pos;
        end = b->core.pos-1;
        CIGAR = bam_get_cigar(b);
        for(i=0; i<b->core.n_cigar; i++) {
            op = bam_cigar_op(CIGAR[i]);
            if(bam_cigar_type(op) == 3) { //M, = or X
                end += bam_cigar_oplen(CIGAR[i]);
            } else if(bam_cigar_type(op) == 2) { //D or N
                if(end >= start) {
                    cnt += countOverlaps(t,
                        hdr->target_name[b->core.tid],
                        start,
                        end+1,
                        (b->core.flag&16)?1:0,
                        matchType,
                        strandType,
                        ffunc);
                }
                start = end + bam_cigar_oplen(CIGAR[i]) + 1;
                end = start-1;
            }
        }
        if(end >= start) {
            cnt += countOverlaps(t,
                hdr->target_name[b->core.tid],
                start,
                end+1,
                (b->core.flag&16)?1:0,
                matchType,
                strandType,
                ffunc);
        }

        printf("%s\t%"PRId32"\n", bam_get_qname(b), cnt);
        cnt = 0;
    }
    bam_destroy1(b);
    free(ks->s);
    free(ks);
}

void usage() {
    fprintf(stderr, "Usage: testCountOverlaps [OPTIONS] <annotation.gtf> <alignments.bam>\n");
    fprintf(stderr, "\n"
"This program demostrates how to use libGTF to count the number of overlaps of\n"
"each alignment in a BAM file. Currently, paired-end reads aren't treated in a\n"
"particularly useful way. It does, however, handle spliced alignments. Note that\n"
"there's no way to only count the number of, for example, exon overlaps. For\n"
"that, one should simply findOverlaps() and then filter the resulting overlapSet.\n"
"\nOPTIONS\n"
"-m STR  Match type. Possible values are 'any', 'exact', 'contain', 'within',\n"
"        'start' and 'end'. These values are equivalent to the 'type' parameter\n"
"        in the findOverlaps() function in GenomicRanges and also Allen's\n"
"        Interval Algebra. 'contain' is simply the opposite of 'within'. The\n"
"        default is 'any'.\n"
"-s STR  Strand type. Possible values are 'ignore', 'same', 'opposite', and\n"
"        'exact'. 'exact' differs from 'same' in how '*' strands are handled.\n"
"        Normally, a subject and query will overlaps if either of them has a '*'\n"
"        strand. The 'exact' option indicates that strands must exactly match.\n"
"-q INT  Minimum MAPQ value. Default is [0].\n"
"-f INT  Denotes a filter function to use. 0: include all entries, 1: exclude\n"
"        entries with gene_biotype set to 'processed_pseudogene', 2: include\n"
"        only entries with gene_biotype set to 'processed_pseudogene'\n"
);
}

int main(int argc, char *argv[]) {
    int matchType = 0;
    int strandType = 0;
    int minMapq = 0;
    FILTER_ENTRY_FUNC ffunc = NULL;
    char c;
    htsFile *fp = NULL;
    bam_hdr_t *hdr = NULL;
    GTFtree *t = NULL;

    opterr = 0; //Disable error messages
    while((c = getopt(argc, argv, "m:s:q:f:h")) >= 0) {
        switch(c) {
        case 'm' :
            if(strcmp(optarg, "any") == 0) matchType = 0;
            else if(strcmp(optarg, "exact") == 0) matchType = 1;
            else if(strcmp(optarg, "contain") == 0) matchType = 2;
            else if(strcmp(optarg, "within") == 0) matchType = 3;
            else if(strcmp(optarg, "start") == 0) matchType = 4;
            else if(strcmp(optarg, "end") == 0) matchType = 5;
            else fprintf(stderr, "Unknown -m option '%s', ignoring\n", optarg);
            break;
        case 's' :
            if(strcmp(optarg, "ignore") == 0) strandType = 0;
            else if(strcmp(optarg, "same") == 0) strandType = 1;
            else if(strcmp(optarg, "opposite") == 0) strandType = 2;
            else if(strcmp(optarg, "exact") == 0) strandType = 3;
            else fprintf(stderr, "Unknown -s option '%s', ignoring\n", optarg);
            break;
        case 'f' :
            if(atoi(optarg) == 1) ffunc = (FILTER_ENTRY_FUNC) noPseudo;
            else if(atoi(optarg) == 2) ffunc = (FILTER_ENTRY_FUNC) onlyPseudo;
            break;
        case 'h' :
            usage();
            return 0;
            break;
        case 'q' :
            minMapq = atoi(optarg);
            if(minMapq < 0) minMapq = 0;
            break;
        case '?' :
        default :
            fprintf(stderr, "Invalid option '%s'\n", argv[optind-1]);
            usage();
            return 1;
            break;
        }
    }

    if(argc == 1) {
        usage();
        return 0;
    }
    if(argc-optind != 2) {
        fprintf(stderr, "Missing either the GTF or BAM file!\n");
        usage();
        return 1;
    }

    //Create the GTFtree
    t = GTF2Tree(argv[optind], NULL);
    if(!t) {
        fprintf(stderr, "Couldn't open %s or there was a problem parsing it.\n", argv[1]);
        return 1;
    }
    sortGTF(t);

    //Open the BAM file
    fp = sam_open(argv[optind+1], "r");
    if(!fp) {
        fprintf(stderr, "Couldn't open %s for reading!\n", argv[optind+1]);
        destroyGTFtree(t);
    }
    hdr = sam_hdr_read(fp);
    if(!hdr) {
        fprintf(stderr, "Couldn't read the header from %s!\n", argv[optind+1]);
        destroyGTFtree(t);
        sam_close(fp);
    }

    countOverlapsBAM2(t, fp, hdr, matchType, strandType, minMapq, ffunc);

    destroyGTFtree(t);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    return 0;
}
