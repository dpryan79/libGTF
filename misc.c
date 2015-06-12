#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <pcre.h>
#include "gtf.h"

static pcre *GTFre;

char *nextField(char *str) {
    char *p = strtok(str, "\t");
    if(p) return p;
    p = strtok(str, "\r");
    if(p) return p;
    p = strtok(str, "\n");
    return p;
}

int initGTFre() {
    const char *errPtr = NULL;
    int errorOffset = 0;
    GTFre = pcre_compile("(?:([\\w][\\w_]*) (?|\"([^\"]+)\"|([^\"]+))[;|\r|\n])+", 0, &errPtr, &errorOffset, NULL);
    if(!GTFre) {
        fprintf(stderr, "[initGTFre] Error while compiling regular expression @%d\n", errorOffset);
        fprintf(stderr, "The error was: %s\n", errPtr);
        return 0;
    }
    return 1;
}

void destroyGTFre() {
    pcre_free(GTFre);
}

static Attribute *makeAttribute(GTFline *l, GTFtree *t, char *key, char *value) {
    int32_t idx;
    Attribute *a = malloc(sizeof(Attribute));
    assert(a);

    if(!strExistsHT(t->htAttributes, key)) {
        idx = addHTelement(t->htAttributes, key);
    } else {
        idx = str2valHT(t->htAttributes, key);
    }
    a->key = idx;
    if(!strExistsHT(t->htAttributes, value)) {
        idx = addHTelement(t->htAttributes, value);
    } else {
        idx = str2valHT(t->htAttributes, value);
    }
    a->val = idx;

    return a;
}

void addAttribute(GTFline *l, GTFtree *t, char *key, char *value) {
    Attribute *a = makeAttribute(l, t, key, value);
    l->attrib = realloc(l->attrib, sizeof(Attribute*)*(++l->nAttributes));
    assert(l->attrib);
    l->attrib[l->nAttributes-1] = a;
}

void destroyAttributes(GTFline *l) {
    int i;
    for(i=0; i<l->nAttributes; i++) free(l->attrib[i]);
    l->nAttributes = 0;
}

int addGTFAttributes(GTFline *l, GTFtree *t, char *str) {
    int rc, len = strlen(str);
    int offset = 0, ovector[12];
    char *key, *value;

    while(offset < len && (rc = pcre_exec(GTFre, NULL, str, strlen(str), offset, 0, ovector, 12))) {
        if(rc<0) break;
        assert(rc==3);
        key = strndup(str + ovector[2], ovector[3]-ovector[2]);
        assert(key);
        value = strndup(str + ovector[4], ovector[5]-ovector[4]);
        assert(value);
        addAttribute(l, t, key, value);
        free(key);
        free(value);
        offset = ovector[2*rc-1];
    }
    return 1;
}
