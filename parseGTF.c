#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <zlib.h>
#include <float.h>
#include "gtf.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"


KSTREAM_INIT(gzFile, gzread, 16384)

GTFline * parseGTFline(GTFline *o, GTFtree *t, kstring_t ks) {
    char *p, *end;

    GTFline_reset(o);

    //Chromosome
    p = nextField(ks.s);
    assert(kputs(p, &o->chrom));

    //Source
    p = nextField(NULL);
    assert(kputs(p, &o->source));

    //Feature
    p = nextField(NULL);
    assert(kputs(p, &o->feature));

    //Start
    p = nextField(NULL);
    if(!p) goto err;
    o->start = strtoull(p, &end, 10)-1;
    if(*end) goto err;

    //End
    p = nextField(NULL);
    if(!p) goto err;
    o->end = strtoull(p, &end, 10);
    if(*end) goto err;

    //Score
    p = nextField(NULL);
    if(!p) goto err;
    if(strcmp(p, ".") == 0) o->score = DBL_MAX;
    else {
        o->score = strtod(p, &end);
        if(*end) goto err;
    }

    //Strand (default is 3)
    p = nextField(NULL);
    if(!p) goto err;
    if(*p == '+') {
        o->strand = 0;
    } else if(*p == '-') {
        o->strand = 1;
    }

    //Frame (default is 3)
    p = nextField(NULL);
    if(!p) goto err;
    if(*p == '0') o->frame = 0;
    else if(*p == '1') o->frame = 1;
    else if(*p == '2') o->frame = 2;

    //Attributes
    p = nextField(NULL);
    if(!p) goto err;
    if(!addGTFAttributes(o, t, p)) goto err;

    return o;

err:
    destroyGTFline(o);
    return NULL;
}

GTFtree *GTF2Tree(char *fname, FILTER_FUNC ffunc) {
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
    initGTFre();
    while(ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
        if(strncmp("##FASTA", str.s, 7) == 0) break;
        if(str.s[0] == '#') continue;
        line = parseGTFline(line, o, str);
        if(!line) break;
        if(ffunc != NULL) {
            if(ffunc((void*) line)) addGTFentry(o, line);
        } else {
            addGTFentry(o, line);
        }
    }

    if(line) destroyGTFline(line);
    ks_destroy(ks);
    gzclose(fp);
    if(str.s) free(str.s);
    destroyGTFre();

    return o;
}

void GTFEntry2GTF(kstring_t *ks, GTFtree *t, GTFentry *e) {
    int first = 1;
    char buf[1024];
    int32_t i;
    if(ks->l) ks->l = 0;

    //chrom
    assert(kputs(val2strHT(t->htChroms, e->chrom), ks));
    assert(kputc('\t', ks));

    //Source
    assert(kputs(val2strHT(t->htSources, e->source), ks));
    assert(kputc('\t', ks));

    //Feature
    assert(kputs(val2strHT(t->htFeatures, e->feature), ks));
    assert(kputc('\t', ks));

    //Start
    assert(kputw(e->start+1, ks) != EOF);
    assert(kputc('\t', ks));

    //End
    assert(kputw(e->end, ks) != EOF);
    assert(kputc('\t', ks));

    //Score
    if(e->score == DBL_MAX) assert(kputc('.', ks));
    else {
        snprintf(buf, 1024, "%400.2f", e->score);
        assert(kputs(buf, ks));
    }
    assert(kputc('\t', ks));

    //Strand
    if(e->strand == 0) assert(kputc('+', ks));
    else if(e->strand == 1) assert(kputc('-', ks));
    else assert(kputc('.', ks));
    assert(kputc('\t', ks));

    //Frame
    if(e->frame == 0) assert(kputc('0', ks));
    else if(e->frame == 1) assert(kputc('1', ks));
    else if(e->frame == 2) assert(kputc('2', ks));
    else assert(kputc('.', ks));
    assert(kputc('\t', ks));

    //Attributes
    if(e->gene_id) {
        first = 0;
        assert(kputs("gene_id \"", ks));
        assert(kputs(val2strHT(t->htGenes, e->gene_id), ks));
        assert(kputs("\"", ks));
    }
    if(e->transcript_id) {
        if(!first) kputs("; ", ks);
        first = 0;
        assert(kputs("transcript_id \"", ks));
        assert(kputs(val2strHT(t->htTranscripts, e->transcript_id), ks));
        assert(kputs("\"", ks));
    }
    for(i=0; i<e->nAttributes; i++) {
        if(!first) kputs("; ", ks);
        first = 0;
        assert(kputs(val2strHT(t->htAttributes, e->attrib[i]->key), ks));
        assert(kputs(" \"", ks));
        assert(kputs(val2strHT(t->htAttributes, e->attrib[i]->val), ks));
        assert(kputs("\"", ks));
    }
}    
