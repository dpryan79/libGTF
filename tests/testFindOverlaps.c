#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include "../gtf.h"
#include "../htslib/htslib/sam.h"

void findOverlapsBAM(GTFtree *t, htsFile *fp, bam_hdr_t *hdr, int matchType, int strandType) {
    bam1_t *b = bam_init1();
    overlapSet *os = os_init();
    kstring_t *ks = calloc(1, sizeof(kstring_t));
    int32_t i;
    assert(b);
    assert(os);
    assert(ks);

    while(sam_read1(fp, hdr, b) > 0) {
        if(b->core.tid < 0) continue;
        if(b->core.flag & BAM_FUNMAP) continue;
        os = findOverlaps(os,
                          t,
                          hdr->target_name[b->core.tid],
                          b->core.pos,
                          bam_endpos(b),
                          (b->core.flag&16)?1:0,
                          matchType,
                          strandType);
        if(os->l) {
            for(i=0; i<os->l; i++) {
                GTFEntry2GTF(ks, t, os->overlaps[i]);
                printf("%s overlaps %s\n", bam_get_qname(b), ks->s);
            }
        }
        os_reset(os);
    }
    bam_destroy1(b);
    os_destroy(os);
    free(ks->s);
    free(ks);
}

void usage() {
    fprintf(stderr, "Usage: testFindOverlaps [OPTIONS] <annotation.gtf> <alignments.bam>\n");
    fprintf(stderr, "\n"
"This program demostrates how to use libGTF to find the overlaps of each\n"
"alignment in a BAM file. Currently, paired-end reads aren't treated in a\n"
"particularly useful way.\n"
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
);
}

int main(int argc, char *argv[]) {
    int matchType = 0;
    int strandType = 0;
    char c;
    htsFile *fp = NULL;
    bam_hdr_t *hdr = NULL;
    GTFtree *t = NULL;

    opterr = 0; //Disable error messages
    while((c = getopt(argc, argv, "m:s:h")) >= 0) {
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
        case 'h' :
            usage();
            return 0;
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
    t = GTF2Tree(argv[optind]);
    if(!t) {
        fprintf(stderr, "Couldn't open %s or there was a problem parsing it.\n", argv[1]);
        return 1;
    }
    sortGTF(t);

    //Open the BAM file
    fp = sam_open(argv[optind+1], "rb");
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

    findOverlapsBAM(t, fp, hdr, matchType, strandType);

    destroyGTFtree(t);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    return 0;
}
