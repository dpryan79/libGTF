#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "gtf.h"

//Nodes for the interval tree
typedef struct {
    int32_t center;
    int32_t l;
    GTFentry *start;
    GTFentry *end;
    struct treeNode *left;
    struct treeNode *right;
} treeNode;

//The sizes shouldn't be preset...
GTFtree * initGTFtree() {
    GTFtree *t = calloc(1, sizeof(GTFtree));
    assert(t);

    //Initialize the hash tables
    t->htChroms = initHT(128);
    t->htSources = initHT(128);
    t->htFeatures = initHT(128);
    t->htGenes = initHT(128);
    t->htTranscripts = initHT(128);
    t->htAttributes = initHT(128);

    return t;
}

void destroyGTFentry(GTFentry *e) {
    int32_t i;
    if(e->right) destroyGTFentry(e->right);
    for(i=0; i<e->nAttributes; i++) free(e->attrib[i]);
    if(e->attrib) free(e->attrib);
    free(e);
}

void destroyGTFnode(GTFnode *n) {
    if(n->left) destroyGTFnode(n->left);
    if(n->starts) destroyGTFentry(n->starts);
    if(n->right) destroyGTFnode(n->right);
    free(n);
}

void destroyGTFchrom(GTFchrom *c, int balanced) {
    if(balanced) destroyGTFnode((GTFnode*) c->tree);
    else destroyGTFentry((GTFentry*) c->tree);
    free(c);
}

//This need to handle htTargets and htGenes still
void destroyGTFtree(GTFtree *t) {
    uint32_t i;
    for(i=0; i<t->n_targets; i++) {
        destroyGTFchrom(t->chroms[i], t->balanced);
    }

    destroyHT(t->htChroms);
    destroyHT(t->htSources);
    destroyHT(t->htGenes);
    destroyHT(t->htTranscripts);
    destroyHT(t->htFeatures);
    destroyHT(t->htAttributes);

    free(t->chroms);
    free(t);
}

void addChrom(GTFtree *t) {
    int i;

    t->n_targets++;
    //Grow if needed
    if(t->n_targets >= t->m) {
        t->m++;
        kroundup32(t->m);
        t->chroms = realloc(t->chroms, t->m * sizeof(GTFchrom *));
        assert(t->chroms);
        for(i=t->n_targets-1; i<t->m; i++) t->chroms[i] = NULL;
    }

    assert(!t->chroms[t->n_targets-1]); //We shouldn't be adding over anything...
    t->chroms[t->n_targets-1] = calloc(1,sizeof(GTFchrom));
    assert(t->chroms[t->n_targets-1]);
}

void addGTFentry(GTFtree *t, GTFline *l) {
    int32_t IDchrom, IDgene, IDtranscript, IDfeature, IDsource;
    assert(t->balanced==0); //Should just switch to insertGTFentry(), which remains to be written

    //Chromosome
    if(!strExistsHT(t->htChroms, l->chrom.s)) {
        addChrom(t);
        IDchrom = addHTelement(t->htChroms, l->chrom.s);
    } else {
        IDchrom = str2valHT(t->htChroms, l->chrom.s);
    }
    //Source
    if(!strExistsHT(t->htSources, l->source.s)) {
        IDsource = addHTelement(t->htSources, l->source.s);
    } else {
        IDsource = str2valHT(t->htSources, l->source.s);
    }
    //gene
    if(!strExistsHT(t->htGenes, l->gene.s)) {
        IDgene = addHTelement(t->htGenes, l->gene.s);
    } else {
        IDgene = str2valHT(t->htGenes, l->gene.s);
    }
    //transcript
    if(!strExistsHT(t->htTranscripts, l->transcript.s)) {
        IDtranscript = addHTelement(t->htTranscripts, l->transcript.s);
    } else {
        IDtranscript = str2valHT(t->htTranscripts, l->transcript.s);
    }
    //feature
    if(!strExistsHT(t->htFeatures, l->feature.s)) {
        IDfeature = addHTelement(t->htFeatures, l->feature.s);
    } else {
        IDfeature = str2valHT(t->htFeatures, l->feature.s);
    }

    //Initialize the entry
    GTFentry *e = malloc(sizeof(GTFentry));
    assert(e);
    e->right = NULL;

    //Copy over the values
    e->chrom = IDchrom;
    e->feature = IDfeature;
    e->source = IDsource;
    e->start = l->start;
    e->end = l->end;
    e->strand = l->strand;
    e->frame = l->frame;
    e->score = l->score;
    e->gene_id = IDgene;
    e->transcript_id = IDtranscript;
    e->nAttributes = l->nAttributes;
    e->attrib = l->attrib;
    assert(l->end > l->start);

    if(t->chroms[IDchrom]->tree) {
        e->left = ((GTFentry*) t->chroms[IDchrom]->tree)->left;
        e->left->right = e;
        ((GTFentry*) t->chroms[IDchrom]->tree)->left = e;
    } else {
        t->chroms[IDchrom]->tree = (void *) e;
        e->left = e;
    }
    t->chroms[IDchrom]->n_entries++;

    GTFline_reset(l);
}

/*******************************************************************************
*
* Sorting functions
*
*******************************************************************************/

//If a list is circular to the right, it won't be afterward
/*
uint32_t getVineLength(GTFentry *e) {
    uint32_t l = 0;
    GTFentry *p = e;
    while(p) {
        l++;
        if(p->right == e) p->right = NULL;
        p = p->right;
    }
    return l;
}
*/

GTFentry *getMiddleR(GTFentry *e, uint32_t pos) {
    uint32_t i;
    GTFentry *tmp, *o = e;

    if(!o->right) return o;
    for(i=1; i<pos; i++) {
        assert(o->right);
        o = o->right;
    }
    tmp = o;
    assert(o->right);
    o = o->right;
    tmp->right = NULL;
    return o;
}
GTFentry *getMiddleL(GTFentry *e, uint32_t pos) {
    uint32_t i;
    GTFentry *tmp, *o = e;

    if(!o->left) {
        return o;
    }
    for(i=1; i<pos; i++) {
        assert(o->left);
        o = o->left;
    }
    tmp = o;
    assert(o->left);
    o = o->left;
    tmp->left = NULL;
    return o;
}

int cmpRangesStart(GTFentry *a, GTFentry *b) {
    if(!b && !a) return 0;
    if(!b) return -1;
    if(!a) return 1;
    if(a->start < b->start) return -1;
    if(b->start < a->start) return 1;
    if(b->end < a->end) return 1;
    return -1;
}
int cmpRangesEnd(GTFentry *a, GTFentry *b) {
    if(!b && !a) return 0;
    if(!a) return 1;
    if(!b) return -1;
    if(a->end > b->end) return -1;
    if(b->end > a->end) return 1;
    if(a->start > b->start) return -1;
    return 1;
}

GTFentry *mergeSortStart(GTFentry *a, GTFentry *b) {
    GTFentry *o = a, *last;
    int i = cmpRangesStart(a,b);

    if(i<0) {
        o = a;
        a = a->right;
    } else if(i>0) {
        o = b;
        b = b->right;
    } else{
         return NULL;
    }
    last = o;
    last->right = NULL;

    while((i=cmpRangesStart(a,b))) {
        if(i>0) {
            last->right= b;
            last = b;
            b = b->right;
        } else {
            last->right= a;
            last = a;
            a = a->right;
        }
    }
    last->right = NULL;
    return o;
}

GTFentry *mergeSortEnd(GTFentry *a, GTFentry *b) {
    GTFentry *o = a, *last;
    int i = cmpRangesEnd(a,b);

    if(i<0) {
        o = a;
        a = a->left;
    } else if(i>0) {
        o = b;
        b = b->left;
    } else {
        return NULL;
    }
    last = o;
    last->left = NULL;

    while((i=cmpRangesEnd(a,b))) {
        if(i<0) {
            assert(a != last);
            last->left = a;
            last = a;
            a = a->left;
        } else {
            assert(b != last);
            last->left = b;
            last = b;
            b = b->left;
        }
    }
    last->left = NULL;
    return o;
}

GTFentry *sortTreeStart(GTFentry *e, uint32_t l) {
    if(l==1) return e;
    uint32_t half = l/2;
    GTFentry *middle = getMiddleR(e, half);

    return mergeSortStart(sortTreeStart(e,half), sortTreeStart(middle,half+(l&1)));
}

GTFentry *sortTreeEnd(GTFentry *e, uint32_t l) {
    if(l==1) {
        e->left = NULL; //The list is circular, so...
        return e;
    }
    uint32_t half = l/2;
    assert(e->left);
    assert(e != e->left);
    GTFentry *middle = getMiddleL(e, half);
    assert(e != middle);
    assert(e != e->left);

    return mergeSortEnd(sortTreeEnd(e,half), sortTreeEnd(middle,half+(l&1)));
}

/*******************************************************************************
*
* Functions for interval tree construction
*
*******************************************************************************/

//Note the returned object is the rightmost interval sorted by end position
GTFentry *sortChrom(GTFchrom *c) {
    GTFentry *e = ((GTFentry *)c->tree)->left;
    ((GTFentry*) c->tree)->left = NULL;
    c->tree = (void *) sortTreeStart((GTFentry *) c->tree, c->n_entries);
    e = sortTreeEnd(e, c->n_entries);
    return e;
}

uint32_t getCenter(GTFentry *ends) {
    GTFentry *slow = ends;
    GTFentry *fast = ends;

    while(fast->left && fast->left->left) {
        slow = slow->left;
        fast = fast->left->left;
    }

    return slow->end-1;
}

/*
GTFentry *getRStarts(GTFentry *starts, uint32_t pos) {
    GTFentry *o, *tmp;
    if(!starts->right) return NULL;
    o = starts;
    while(o->right) {
        if(o->right->start>=pos) {
            tmp = o;
            o = o->right;
            tmp->right = NULL;
            return o;
        }
        o = o->right;
    }
    return NULL;
}
GTFentry *getLEnds(GTFentry *ends, uint32_t pos) {
    GTFentry *o, *tmp;
    if(!ends->left) return NULL;
    o = ends;
    while(o->left) {
        if(o->left->end < pos) {
            tmp = o;
            o = o->left;
            tmp->left = NULL;
            return o;
        }
        o = o->left;
    }
    return NULL;
}
*/

GTFentry *getMembers(GTFentry **members, GTFentry **rStarts, GTFentry *starts, uint32_t pos) {
    GTFentry *tmp, *newStarts = NULL;
    GTFentry *last = NULL, *lastMember = NULL;

    *members = NULL, *rStarts = NULL;

    while(starts && starts->start <= pos) {
        if(starts->end > pos) {
            tmp = starts->right;
            if(!*members) {
                lastMember = starts;
                *members = starts;
            } else {
                lastMember->right = starts;
                lastMember = starts;
            }
            starts->right = NULL;
            starts = tmp;
        } else {
            if(!newStarts) {
                newStarts = starts;
                last = starts;
            } else {
                last->right = starts;
                last = starts;
            }
            starts = starts->right;
        }
    }
    *rStarts = starts;
    if(lastMember) lastMember->right = NULL;
    if(last) last->right = NULL;
    assert(*members);
    return newStarts;
}
GTFentry *getRMembers(GTFentry **members, GTFentry **lEnds, GTFentry *ends, uint32_t pos) {
    GTFentry *tmp, *newEnds = NULL;
    GTFentry *last = NULL, *lastMember = NULL;

    *members = NULL, *lEnds = NULL;

    while(ends && ends->end > pos) {
        tmp = ends->left;
        if(ends->start <= pos) {
            if(!*members) {
                *members = ends;
                lastMember = ends;
            } else {
                lastMember->left = ends;
                lastMember = ends;
            }
        } else {
            if(!newEnds) {
                newEnds = ends;
                last = ends;
            } else {
                last->left = ends;
                last = ends;
            }
        }
        ends->left = NULL;
        ends = tmp;
    }
    *lEnds = ends;
    assert(*members);
    lastMember->left = NULL;
    if(newEnds) last->left = NULL;
    return newEnds;
}

/*
GTFentry *removeMembers(GTFentry *ends, uint32_t pos) {
    GTFentry *newEnds = NULL, *lastEnd = NULL;
    GTFentry *members = NULL, *lastMember = NULL;

    while(ends) {
        if(ends->start > pos) {
            if(!newEnds) {
                newEnds = ends;
                lastEnd = newEnds;
            } else {
                lastEnd->left = ends;
                lastEnd = lastEnd->left;
            }
        } else {
            if(!members) {
                members = ends;
                lastMember = members;
            } else {
                lastMember->left = ends;
            }
            lastMember = ends;
        }
        ends = ends->left;
        if(lastEnd) lastEnd->left = NULL;
        if(lastMember) lastMember->left = NULL;
    }
    if(lastMember) lastMember->left = members; //Make it circular so resorting works
    return newEnds;
}

static GTFentry *fixEnds(GTFentry *starts) {
    GTFentry *curr = starts->right;
    GTFentry *prev = starts;
    while(curr) {
        curr->left = prev;
        prev = curr;
        curr = curr->right;
    }
    starts->left = NULL;
    if(curr) return curr;
    return prev;
}
*/

GTFnode *makeIntervalTree(GTFentry *starts, GTFentry *ends) {
    uint32_t center = getCenter(ends);//, nMembers;
    GTFentry *rStarts = NULL; //getRStarts(starts, center);
    GTFentry *lEnds = NULL; //getLEnds(ends, center);
    GTFentry *memberStarts = NULL, *memberEnds = NULL;
    GTFnode *out = calloc(1, sizeof(GTFnode));
    assert(out);

    starts = getMembers(&memberStarts, &rStarts, starts, center);
    ends = getRMembers(&memberEnds, &lEnds, ends, center);

    out->center = center;
    out->starts = memberStarts;
    out->ends = memberEnds;
    if(lEnds && starts) {
        out->left = makeIntervalTree(starts, lEnds);
    } else {
        out->left = NULL;
    }
    if(rStarts && ends) {
        out->right = makeIntervalTree(rStarts, ends);
    } else {
        out->right = NULL;
    }

    return out;
}

void sortGTF(GTFtree *t) {
    int32_t i;
    GTFentry *ends;

    for(i=0; i<t->n_targets; i++) {
        ends = sortChrom(t->chroms[i]);
        t->chroms[i]->tree = (void*) makeIntervalTree((GTFentry*) t->chroms[i]->tree, ends);
    }
    t->balanced = 1;
}

/*******************************************************************************
*
* GTFline-specific functions
*
*******************************************************************************/
void GTFline_reset(GTFline *l) {
    if(l->chrom.s) {
        l->chrom.l = 0;
        l->chrom.s[0] = '\0';
    }
    if(l->feature.s) {
        l->feature.l = 0;
        l->chrom.s[0] = '\0';
    }
    if(l->source.s) {
        l->source.l = 0;
        l->source.s[0] = '\0';
    }
    if(l->gene.s) {
        l->gene.l = 0;
        l->gene.s[0] = '\0';
    }
    if(l->transcript.s) {
        l->transcript.l = 0;
        l->transcript.s[0] = '\0';
    }
    l->strand = 3;
    l->frame = 3;
    l->nAttributes = 0;
    l->attrib = NULL;
}

void destroyGTFline(GTFline *l) {
    int32_t i;

    if(l->chrom.s) free(l->chrom.s);
    if(l->feature.s) free(l->feature.s);
    if(l->source.s) free(l->source.s);
    if(l->gene.s) free(l->gene.s);
    if(l->transcript.s) free(l->transcript.s);
    if(l->nAttributes) {
        for(i=0; i<l->nAttributes; i++) {
            free(l->attrib[i]);
        }
        free(l->attrib);
    }
    free(l);
}

GTFline *initGTFline() {
    GTFline *o = calloc(1, sizeof(GTFline));
    if(!o) return NULL;

    o->strand = 3;
    o->frame = 3;
    return o;
}

/*******************************************************************************
*
* Accessors
*
*******************************************************************************/
char *GTFgetGeneID(GTFtree *t, GTFentry *e) {
    return val2strHT(t->htGenes, e->gene_id);
}

/*******************************************************************************
*
* Misc. functions
*
*******************************************************************************/
void printBalancedGTF(GTFnode *n, const char *chrom) {
    kstring_t ks, ks2;
    ks.s = NULL; ks.l = ks.m = 0;
    ks2.s = NULL; ks2.l = ks2.m = 0;
    kputs(chrom, &ks);
    kputc(':', &ks);
    kputuw(n->center, &ks);
    if(n->left) {
        kputs(chrom, &ks2);
        kputc(':', &ks2);
        kputuw(n->left->center, &ks2);
        printf("\t\"%s\" -> \"%s\";\n", ks.s, ks2.s);
        printBalancedGTF(n->left, chrom);
    }

    printf("\t\"%s:%"PRIu32"\" [shape=box];\n", chrom, n->center);
    GTFentry *e = n->starts;
    if(e) printGTFvineStart(e, chrom, ks.s);
    if(n->ends) printGTFvineStartR(n->ends, chrom, ks.s);

    if(n->right) {
        ks2.l = 0;
        kputs(chrom, &ks2);
        kputc(':', &ks2);
        kputuw(n->right->center, &ks2);
        printf("\t\"%s\" -> \"%s\";\n", ks.s, ks2.s);
        printBalancedGTF(n->right, chrom);
    }
    free(ks.s);
    if(ks2.s) free(ks2.s);
}

void printGTFvineR(GTFentry *e, const char* chrom) {
    if(e->left == e) return;
    printf("\t\"%s:%"PRIu32"-%"PRIu32"\" -> \"%s:%"PRIu32"-%"PRIu32"\" [color=red];\n", chrom, e->start, e->end, chrom, e->left->start, e->left->end);
    if(!e->left) return;
    printGTFvineR(e->left, chrom);
}
void printGTFvineStartR(GTFentry *e, const char *chrom, const char *str) {
    printf("\t\"%s\" -> \"%s:%"PRIu32"-%"PRIu32"\" [color=red];\n", str, chrom, e->start, e->end);
    if(e->left) printGTFvineR(e, chrom);
}

void printGTFvine(GTFentry *e, const char* chrom) {
    printf("\t\"%s:%"PRIu32"-%"PRIu32"\" -> \"%s:%"PRIu32"-%"PRIu32"\";\n", chrom, e->start, e->end, chrom, e->right->start, e->right->end);
    if(!e->right) return;
    printGTFvine(e->right, chrom);
}
void printGTFvineStart(GTFentry *e, const char *chrom, const char *str) {
    printf("\t\"%s\" -> \"%s:%"PRIu32"-%"PRIu32"\";\n", str, chrom, e->start, e->end);
    if(e->right) printGTFvine(e, chrom);
}

void printGTFtree(GTFtree *t) {
    int32_t i;
    const char *chromName;

    if(t->balanced) printf("digraph balancedTree {\n");
    else printf("digraph unbalancedTree {\n");

    for(i=0; i<t->n_targets; i++) {
        chromName = val2strHT(t->htChroms, i);
        if(t->balanced) {
            printBalancedGTF((GTFnode*) t->chroms[i]->tree, chromName);
        } else {
            printGTFvineStart((GTFentry*) t->chroms[i]->tree, chromName, chromName);
        }
    }
    printf("}\n");
}
