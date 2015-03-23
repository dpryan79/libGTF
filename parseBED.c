#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <zlib.h>
#include <float.h>
#include <math.h>
#include "gtf.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"

KSTREAM_INIT(gzFile, gzread, 16384)

GTFline * parseBEDline(GTFline *o, GTFtree *t, kstring_t ks) {
    char *p, *end;

    GTFline_reset(o);

    //Chromosome
    p = nextField(ks.s);
    assert(kputs(p, &o->chrom));

    //Start
    p = nextField(NULL);
    if(!p) goto err;
    o->start = strtoull(p, &end, 10);
    if(*end) goto err;

    //End
    p = nextField(NULL);
    if(!p) goto err;
    o->end = strtoull(p, &end, 10);
    if(*end) goto err;

    //name, stored as "gene"
    p = nextField(NULL);
    if(!p) return o;
    assert(kputs(p, &o->gene));

    //Score, not handled
    p = nextField(NULL);
    if(!p) return o;
    if(strcmp(p, ".") == 0) o->score = DBL_MAX;
    else {
        o->score = strtod(p, &end);
        if(*end) goto err;
    }

    //Strand
    p = nextField(NULL);
    if(!p) return o;
    if(*p == '+') {
        o->strand = 0;
    } else if(*p == '-') { 
        o->strand = 1;
    }

    //thickStart
    p = nextField(NULL);
    if(!p) return o;
    addAttribute(o, t, "thickStart", p);

    //thickEnd
    p = nextField(NULL);
    if(!p) return o;
    addAttribute(o, t, "thickEnd", p);

    //itemRgb
    p = nextField(NULL);
    if(!p) return o;
    addAttribute(o, t, "itemRgb", p);

    //blockCount
    p = nextField(NULL);
    if(!p) return o;
    addAttribute(o, t, "blockCount", p);

    //blockSizes
    p = nextField(NULL);
    if(!p) return o;
    addAttribute(o, t, "blockSizes", p);

    //blockStarts
    p = nextField(NULL);
    if(!p) return o;
    addAttribute(o, t, "blockStarts", p);

    return o;

err :
    destroyGTFline(o);
    return NULL;
}

GTFtree *BED2Tree(char *fname) {
    gzFile fp = gzopen(fname, "r");
    GTFtree *o = NULL;
    GTFline *line = initGTFline();
    assert(line);
    int dret;
    kstream_t *ks = ks_init(fp);
    kstring_t str;
    str.s = NULL;
    str.l = str.m = 0;
    if(!fp) return NULL;

    o = initGTFtree();
    while(ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
        line = parseBEDline(line, o, str);
        if(!line) break;
        addGTFentry(o, line);
    }

    if(line) destroyGTFline(line);
    ks_destroy(ks);
    gzclose(fp);
    if(str.s) free(str.s);

    return o;
}

void GTFEntry2BED(kstring_t *ks, GTFtree *t, GTFentry *e, int ncols) {
    char *val;
    if(ks->l) ks->l = 0;

    //chrom
    assert(kputs(val2strHT(t->htChroms, e->chrom), ks));
    assert(kputc('\t', ks));

    //Start
    assert(kputuw(e->start, ks));
    assert(kputc('\t', ks));

    //End
    assert(kputuw(e->end, ks));
    if(ncols <= 3) return;

    //Name (in gene_id)
    assert(e->gene_id);
    assert(kputc('\t', ks));
    assert(kputs(val2strHT(t->htGenes, e->gene_id), ks));
    if(ncols <= 4) return;
    
    //Score
    assert(kputc('\t', ks));
    if(e->score == DBL_MAX) assert(kputc('.', ks));
    else assert(kputw(lrint(e->score), ks));
    if(ncols <= 5) return;

    //Strand
    assert(kputc('\t', ks));
    if(e->strand == 0) assert(kputc('+', ks));
    else if(e->strand == 1) assert(kputc('-', ks));
    else assert(kputc('.', ks));
    if(ncols <= 6) return;

    //thickStart, defaults to start
    assert(kputc('\t', ks));
    val = getAttribute(t, e, "thickStart");
    if(val) assert(kputs(val, ks));
    else assert(kputuw(e->start, ks));
    if(ncols <= 7) return;

    //thickEnd
    assert(kputc('\t', ks));
    val = getAttribute(t, e, "thickEnd");
    if(val) assert(kputs(val, ks));
    else assert(kputw(e->end, ks));
    if(ncols <= 8) return;

    //itemRgb
    assert(kputc('\t', ks));
    val = getAttribute(t, e, "itemRgb");
    if(val) assert(kputs(val, ks));
    else assert(kputc('.', ks));
    if(ncols <= 9) return;

    //blockCount
    assert(kputc('\t', ks));
    val = getAttribute(t, e, "blockCount");
    if(val) assert(kputs(val, ks));
    else assert(kputw(1, ks));
    if(ncols <= 10) return;

    //blockSizes
    assert(kputc('\t', ks));
    val = getAttribute(t, e, "blockSizes");
    if(val) assert(kputs(val, ks));
    else assert(kputuw(e->end-e->start, ks));
    if(ncols <= 11) return;

    //blockStarts
    assert(kputc('\t', ks));
    val = getAttribute(t, e, "blockStarts");
    if(val) assert(kputs(val, ks));
    else assert(kputuw(e->start, ks));
}
