#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "htslib/sam.h"
#include "gtf.h"

/*******************************************************************************
*
* Comparison functions
*
* These are according to Allen's Interval Algebra
*
*******************************************************************************/
static inline int rangeAny(uint32_t start, uint32_t end, GTFentry *e) {
    if(end <= e->start) return -1;
    if(start >= e->end) return 1;
    return 0;
}

static inline int rangeContains(uint32_t start, uint32_t end, GTFentry *e) {
    if(e->start >= start && e->end <= end) return 0;
    if(e->end < end) return -1;
    return 1;
}

static inline int rangeWithin(uint32_t start, uint32_t end, GTFentry *e) {
    if(start >= e->start && end <= e->end) return 0;
    if(start < e->start) return -1;
    return 1;
}

static inline int rangeExact(uint32_t start, uint32_t end, GTFentry *e) {
    if(start == e->start && end == e->end) return 0;
    if(start < e->start) return -1;
    if(end < e->end) return -1;
    return 1;
}

static inline int rangeStart(uint32_t start, uint32_t end, GTFentry *e) {
    if(start == e->start) return 0;
    if(start < e->start) return -1;
    return 1;
}

static inline int rangeEnd(uint32_t start, uint32_t end, GTFentry *e) {
    if(end == e->end) return 0;
    if(end < e->end) return -1;
    return 1;
}

static inline int exactSameStrand(int strand, GTFentry *e) {
    return strand == e->strand;
}

static inline int sameStrand(int strand, GTFentry *e) {
    if(strand == 3 || e->strand == 3) return 1;
    if(strand == e->strand) return 1;
    return 0;
}

static inline int oppositeStrand(int strand, GTFentry *e) {
    if(strand == 3 || e->strand == 3) return 1;
    if(strand != e->strand) return 1;
    return 0;
}

/*******************************************************************************
*
* OverlapSet functions
*
*******************************************************************************/
overlapSet *os_init(GTFtree *t) {
    overlapSet *os = calloc(1, sizeof(overlapSet));
    assert(os);
    os->tree = t;
    return os;
}

void os_reset(overlapSet *os) {
    int i;
    for(i=0; i<os->l; i++) os->overlaps[i] = NULL;
    os->l = 0;
}

void os_destroy(overlapSet *os) {
    if(os->overlaps) free(os->overlaps);
    free(os);
}

overlapSet *os_grow(overlapSet *os) {
    int i;
    os->m++;
    kroundup32(os->m);
    os->overlaps = realloc(os->overlaps, os->m * sizeof(GTFentry*));
    assert(os->overlaps);
    for(i=os->l; i<os->m; i++) os->overlaps[i] = NULL;

    return os;
}

static void os_push(overlapSet *os, GTFentry *e) {
    if(os->l+1 >= os->m) os = os_grow(os);
    os->overlaps[os->l++] = e;
}

void os_exclude(overlapSet *os, int i) {
    int j;
    for(j=i; j<os->l-1; j++) os->overlaps[j] = os->overlaps[j+1];
    os->overlaps[--os->l] = NULL;
}

static int os_sortFunc(const void *a, const void *b) {
    GTFentry *pa = (GTFentry*) a;
    GTFentry *pb = (GTFentry*) b;

    if(pa->start < pb->start) return -1;
    if(pb->start < pa->start) return 1;
    if(pa->end < pb->end) return -1;
    if(pb->end < pa->end) return 1;
    return 0;
}

static void os_sort(overlapSet *os) {
    qsort((void *) os->overlaps, os->l, sizeof(GTFentry*), os_sortFunc);
}

/*******************************************************************************
*
* uniqueSet functions
*
*******************************************************************************/
static uniqueSet *us_init(hashTable *ht) {
    uniqueSet *us = calloc(1, sizeof(uniqueSet));
    assert(us);
    us->ht = ht;
    return us;
}

void us_destroy(uniqueSet *us) {
    if(!us) return;
    if(us->IDs) free(us->IDs);
    free(us);
}

static uniqueSet *us_grow(uniqueSet *us) {
    int i;
    us->m++;
    kroundup32(us->m);
    us->IDs = realloc(us->IDs, us->m * sizeof(int32_t));
    assert(us->IDs);
    for(i=us->l; i<us->m; i++) us->IDs[i] = -1;

    return us;
}

static void us_push(uniqueSet *us, int32_t ID) {
    if(us->l+1 >= us->m) us = us_grow(us);
    us->IDs[us->l++] = ID;
}

char *us_val(uniqueSet *us, int32_t i) {
    if(i>=us->l) return NULL;
    return val2strHT(us->ht, us->IDs[i]);
}

/*******************************************************************************
*
* Overlap set count/unique functions
*
*******************************************************************************/
static int int32_t_cmp(const void *a, const void *b) {
    int32_t ia = *((int32_t*) a);
    int32_t ib = *((int32_t*) b);
    return ia-ib;
}

int32_t cntGeneIDs(overlapSet *os) {
    int32_t i, last, n = 0;
    int32_t IDs[os->l];
    if(os->l == 0) return 0;
    if(os->l == 1) return 1;

    for(i = 0; i<os->l; i++) IDs[i] = os->overlaps[i]->gene_id;
    qsort((void*) IDs, os->l, sizeof(int32_t), int32_t_cmp);

    last = IDs[0];
    n = 1;
    for(i = 1; i<os->l; i++) {
        if(IDs[i] != last) {
            n++;
            last = IDs[i];
        }
    }
    return n;
}

int32_t cntTranscriptIDs(overlapSet *os) {
    int32_t i, last, n = 0;
    int32_t IDs[os->l];
    if(os->l == 0) return 0;
    if(os->l == 1) return 1;

    for(i = 0; i<os->l; i++) IDs[i] = os->overlaps[i]->transcript_id;
    qsort((void*) IDs, os->l, sizeof(int32_t), int32_t_cmp);

    last = IDs[0];
    n = 1;
    for(i = 1; i<os->l; i++) {
        if(IDs[i] != last) {
            n++;
            last = IDs[i];
        }
    }
    return n;
}

int32_t cntAttributes(overlapSet *os, char *attributeName) {
    int32_t IDs[os->l], i, j, key, last, n = 0;
    if(!strExistsHT(os->tree->htAttributes, attributeName)) return n;

    key = str2valHT(os->tree->htAttributes, attributeName);
    for(i=0; i<os->l; i++) {
        IDs[i] = -1;
        for(j=0; j<os->overlaps[i]->nAttributes; j++) {
            if(os->overlaps[i]->attrib[j]->key == key) {
                IDs[i] = os->overlaps[i]->attrib[j]->val;
                break;
            }
        }
    }
    qsort((void*) IDs, os->l, sizeof(int32_t), int32_t_cmp);

    last = IDs[0];
    n = (last >= 0) ? 1 : 0;
    for(i = 1; i<os->l; i++) {
        if(IDs[i] != last) {
            n++;
            last = IDs[i];
        }
    }
    return n;
}

uniqueSet *uniqueGeneIDs(overlapSet *os) {
    if(os->l == 0) return NULL;
    int32_t i, last;
    int32_t IDs[os->l];
    uniqueSet *us = us_init(os->tree->htGenes);

    for(i = 0; i<os->l; i++) IDs[i] = os->overlaps[i]->gene_id;
    qsort((void*) IDs, os->l, sizeof(int32_t), int32_t_cmp);

    last = IDs[0];
    if(IDs[0] >= 0) us_push(us, last);
    for(i=1; i<os->l; i++) {
        if(IDs[i] != last) {
            us_push(us, IDs[i]);
            last = IDs[i];
        }
    }

    if(us->l) return us;
    us_destroy(us);
    return NULL;
}

uniqueSet *uniqueTranscriptIDs(overlapSet *os) {
    if(os->l == 0) return NULL;
    int32_t i, last;
    int32_t IDs[os->l];
    uniqueSet *us = us_init(os->tree->htTranscripts);

    for(i = 0; i<os->l; i++) IDs[i] = os->overlaps[i]->transcript_id;
    qsort((void*) IDs, os->l, sizeof(int32_t), int32_t_cmp);

    last = IDs[0];
    if(IDs[0] >= 0) us_push(us, last);
    for(i=1; i<os->l; i++) {
        if(IDs[i] != last) {
            us_push(us, IDs[i]);
            last = IDs[i];
        }
    }

    if(us->l) return us;
    us_destroy(us);
    return NULL;
}

uniqueSet *uniqueAttributes(overlapSet *os, char *attributeName) {
    if(os->l == 0) return NULL;
    int32_t IDs[os->l], i, j, key, last;
    if(!strExistsHT(os->tree->htAttributes, attributeName)) return NULL;
    uniqueSet *us = us_init(os->tree->htAttributes);

    key = str2valHT(os->tree->htAttributes, attributeName);
    for(i=0; i<os->l; i++) {
        IDs[i] = -1;
        for(j=0; j<os->overlaps[i]->nAttributes; j++) {
            if(os->overlaps[i]->attrib[j]->key == key) {
                IDs[i] = os->overlaps[i]->attrib[j]->val;
                break;
            }
        }
    }
    qsort((void*) IDs, os->l, sizeof(int32_t), int32_t_cmp);

    last = IDs[0];
    if(IDs[0] >= 0) us_push(us, last);
    for(i=1; i<os->l; i++) {
        if(IDs[i] != last) {
            us_push(us, IDs[i]);
            last = IDs[i];
        }
    }

    if(us->l) return us;
    us_destroy(us);
    return NULL;
}

/*******************************************************************************
*
* Node iterator functions
*
*******************************************************************************/
//bit 1: go left, bit 2: go right (a value of 3 is then "do both"
static int centerDirection(uint32_t start, uint32_t end, GTFnode *n) {
    if(n->center >= start && n->center < end) return 3;
    if(n->center < start) return 2;
    return 1;
}

static int matchingStrand(GTFentry *e, int strand, int strandType) {
    if(strandType == GTF_IGNORE_STRAND) return 1;

    if(strandType == GTF_SAME_STRAND) {
        return sameStrand(strand, e);
    } else if(strandType == GTF_OPPOSITE_STRAND) {
        return oppositeStrand(strand, e);
    } else if(strandType == GTF_EXACT_SAME_STRAND) {
        return exactSameStrand(strand, e);
    }

    fprintf(stderr, "[matchingStrand] Unknown strand type %i. Assuming a match.\n", strandType);
    return 1;
}

static void filterStrand(overlapSet *os, int strand, int strandType) {
    int i;

    if(strandType == GTF_IGNORE_STRAND) return;

    for(i=os->l-1; i>=0; i--) {
        if(strandType == GTF_SAME_STRAND) {
            if(!sameStrand(strand, os->overlaps[i])) os_exclude(os, i);
        } else if(strandType == GTF_OPPOSITE_STRAND) {
            if(!oppositeStrand(strand, os->overlaps[i])) os_exclude(os, i);
        } else if(strandType == GTF_EXACT_SAME_STRAND) {
            if(!exactSameStrand(strand, os->overlaps[i])) os_exclude(os, i);
        }
    }
}

static void pushOverlaps(overlapSet *os, GTFentry *e, uint32_t start, uint32_t end, int comparisonType, int direction) {
    int dir;
    if(!e) return;

    switch(comparisonType) {
    case GTF_MATCH_EXACT :
        if((dir = rangeExact(start, end, e)) == 0) {
            os_push(os, e);
        }
        break;
    case GTF_MATCH_WITHIN :
        if((dir = rangeAny(start, end, e)) == 0) {
            if(rangeWithin(start, end ,e) == 0) os_push(os, e);
        }
        break;
    case GTF_MATCH_CONTAIN :
        if((dir = rangeAny(start, end, e)) == 0) {
            if(rangeContains(start, end, e) == 0) os_push(os, e);
        }
        break;
    case GTF_MATCH_START :
        if((dir = rangeStart(start, end, e)) == 0) {
            os_push(os, e);
        }
        break;
    case GTF_MATCH_END :
        if((dir = rangeEnd(start, end, e)) == 0) {
            os_push(os, e);
        }
        break;
    default :
        if((dir = rangeAny(start, end, e)) == 0) {
            os_push(os, e);
        }
        break;
    }

    if(direction) {
        if(dir > 0) return;
        pushOverlaps(os, e->right, start, end, comparisonType, direction);
    } else {
        if(dir < 0) return;
        pushOverlaps(os, e->left, start, end, comparisonType, direction);
    }
}

static int32_t countOverlapsEntry(GTFentry *e, uint32_t start, uint32_t end, int strand, int matchType, int strandType, int direction, int32_t max) {
    int dir;
    int32_t cnt = 0;
    if(!e) return cnt;
    
    switch(matchType) {
    case GTF_MATCH_EXACT :
        if((dir = rangeExact(start, end, e)) == 0) {
            cnt = 1;
        }
        break;
    case GTF_MATCH_WITHIN :
        if((dir = rangeAny(start, end, e)) == 0) {
            if(rangeWithin(start, end, e) == 0) cnt = 1;
        }
        break;
    case GTF_MATCH_CONTAIN :
        if((dir = rangeAny(start, end, e)) == 0) {
            if(rangeContains(start, end, e) == 0) cnt = 1;
        }
        break;
    case GTF_MATCH_START :
        if((dir = rangeStart(start, end, e)) == 0) {
            cnt = 1;
        }
        break;
    case GTF_MATCH_END :
        if((dir = rangeEnd(start, end, e)) == 0) {
            cnt = 1;
        }
        break;
    default :
        if((dir = rangeAny(start, end, e)) == 0) {
            cnt = 1;
        }
        break;
    }

    if(cnt) {
        if(!matchingStrand(e, strand, strandType)) cnt = 0;
    }

    if(cnt >= max) return max;

    if(direction) {
        if(dir > 0) return cnt;
        return cnt + countOverlapsEntry(e->right, start, end, strand, matchType, strandType, direction, max);
    } else {
        if(dir < 0) return cnt;
        return cnt + countOverlapsEntry(e->left, start, end, strand, matchType, strandType, direction, max);
    }
}

static void pushOverlapsNode(overlapSet *os, GTFnode *n, uint32_t start, uint32_t end, int matchType) {
    int dir;
    if(!n) return;
    dir = centerDirection(start, end, n);

    if(dir&1) {
        pushOverlaps(os, n->starts, start, end, matchType, 1);
        pushOverlapsNode(os, n->left, start, end, matchType);
    } 
    if(dir&2) {
        if(dir!=3) pushOverlaps(os, n->ends, start, end, matchType, 0);
        pushOverlapsNode(os, n->right, start, end, matchType);
    }
}

static int32_t countOverlapsNode(GTFnode *n, uint32_t start, uint32_t end, int strand, int matchType, int strandType, int32_t max) {
    int32_t cnt = 0;
    int dir;
    if(!n) return cnt;
    dir = centerDirection(start, end, n);

    if(dir&1) {
        cnt += countOverlapsEntry(n->starts, start, end, strand, matchType, strandType, 1, max);
        if(max && cnt >= max) return max;
        cnt += countOverlapsNode(n->left, start, end, strand, matchType, strandType, max);
        if(max && cnt >= max) return max;
    } 
    if(dir&2) {
        if(dir!=3) cnt += countOverlapsEntry(n->starts, start, end, strand, matchType, strandType, 0, max);
        if(max && cnt >= max) return max;
        cnt += countOverlapsNode(n->right, start, end, strand, matchType, strandType, max);
        if(max && cnt >= max) return max;
    }
    return cnt;
}

/*******************************************************************************
*
* Driver functions for end use.
*
*******************************************************************************/
overlapSet * findOverlaps(overlapSet *os, GTFtree *t, char *chrom, uint32_t start, uint32_t end, int strand, int matchType, int strandType, int keepOS) {
    int32_t tid = str2valHT(t->htChroms, chrom);
    overlapSet *out = os;

    if(out && !keepOS) os_reset(out);
    else out = os_init(t);

    if(tid<0) return out;
    if(!t->balanced) {
        fprintf(stderr, "[findOverlaps] The tree has not been balanced! No overlaps will be returned.\n");
        return out;
    }

    pushOverlapsNode(out, (GTFnode*) t->chroms[tid]->tree, start, end, matchType);
    if(out->l) filterStrand(out, strand, strandType);
    if(out->l) os_sort(out);

    return out;
}

int32_t countOverlaps(GTFtree *t, char *chrom, uint32_t start, uint32_t end, int strand, int matchType, int strandType) {
    int32_t tid = str2valHT(t->htChroms, chrom);
    if(tid<0) return 0;

    if(!t->balanced) {
        fprintf(stderr, "[countOverlaps] The tree has not been balanced! No overlaps will be returned.\n");
        return 0;
    }

    return countOverlapsNode((GTFnode*) t->chroms[tid]->tree, start, end, strand, matchType, strandType, 0);
}

int overlapsAny(GTFtree *t, char *chrom, uint32_t start, uint32_t end, int strand, int matchType, int strandType) {
    int32_t tid = str2valHT(t->htChroms, chrom);
    if(tid<0) return 0;

    if(!t->balanced) {
        fprintf(stderr, "[overlapsAny] The tree has not been balanced! No overlaps will be returned.\n");
        return 0;
    }

    return countOverlapsNode((GTFnode*) t->chroms[tid]->tree, start, end, strand, matchType, strandType, 1);
}

/*******************************************************************************
*
* Convenience functions for alignments
*
*******************************************************************************/
overlapSet *findOverlapsBAM(overlapSet *os, GTFtree *t, bam1_t *b, bam_hdr_t *hdr, int strand, int matchType, int strandType) {
    int32_t i, start, end;
    uint32_t *CIGAR, op;
    char *chrom = NULL;
    overlapSet *out = os;

    if(out) os_reset(out);
    else out = os_init(t);

    if(b->core.tid < 0 || (b->core.flag & BAM_FUNMAP)) return out;
    chrom = hdr->target_name[b->core.tid];

    CIGAR = bam_get_cigar(b);
    start = b->core.pos;
    end = b->core.pos-1;
    for(i=0; i<b->core.n_cigar; i++) {
        op = bam_cigar_op(CIGAR[i]);
        if(bam_cigar_type(op) == 3) { //M, = or X
            end += bam_cigar_oplen(CIGAR[i]);
        } else if(bam_cigar_type(op) == 2) { //D or N
            if(end >= start) {
                os = findOverlaps(os, t, chrom, start, end+1, (b->core.flag&16)?1:0, matchType, strandType, 1);
            }
            start = end + bam_cigar_oplen(CIGAR[i]) + 1;
            end = start-1;
        }
    }
    if(end >= start) {
        os = findOverlaps(os, t, chrom, start, end+1, (b->core.flag&16)?1:0, matchType, strandType, 1);
    }

    return os;
}
