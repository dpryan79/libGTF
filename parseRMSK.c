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

GTFline * parseRMSKline(GTFline *o, GTFtree *t, kstring_t ks) {
    char *p, *end;

    if(o) GTFline_reset(o);
    else o = initGTFline();

    //bin
    p = nextField(ks.s);
    if(!p) return o;

    //Score
    p = nextField(NULL);
    if(!p) return o;
    if(strcmp(p, ".") == 0) o->score = DBL_MAX;
    else {
        o->score = strtod(p, &end);
        if(*end) goto err;
    }

    //milliDiv
    p = nextField(NULL);
    if(!p) return o;

    //milliDel
    p = nextField(NULL);
    if(!p) return o;

    //milliIns
    p = nextField(NULL);
    if(!p) return o;

    //Chromosome
    p = nextField(NULL);
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

    //genoLeft
    p = nextField(NULL);
    if(!p) return o;

    //Strand
    p = nextField(NULL);
    if(!p) return o;
    if(*p == '+') {
        o->strand = 0;
    } else if(*p == '-') { 
        o->strand = 1;
    }

    //repName
    p = nextField(NULL);
    if(!p) return o;
    addAttribute(o, t, "repName", p);

    //repClass
    p = nextField(NULL);
    if(!p) return o;
    addAttribute(o, t, "repClass", p);

    //repFamily
    p = nextField(NULL);
    if(!p) return o;
    addAttribute(o, t, "repFamily", p);

    //repClass
    

    return o;

err :
    destroyGTFline(o);
    return NULL;
}

GTFtree *RMSK2Tree(char *fname, FILTER_FUNC ffunc) {
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
        line = parseRMSKline(line, o, str);
        if(*str.s == '#') continue;
        if(!line) break;
        if(ffunc == NULL) {
            addGTFentry(o, line);
        } else if(ffunc((void *) line)) {
            addGTFentry(o, line);
        }
    }

    if(line) destroyGTFline(line);
    ks_destroy(ks);
    gzclose(fp);
    if(str.s) free(str.s);

    return o;
}
