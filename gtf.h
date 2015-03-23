#include <inttypes.h>
#include "htslib/kstring.h"

/*****************
 * Strand macros *
 *****************/
#define GTF_IGNORE_STRAND     0
#define GTF_SAME_STRAND       1
#define GTF_OPPOSITE_STRAND   2
#define GTF_EXACT_SAME_STRAND 3

/***********************
 * Overlap type macros *
 ***********************/
#define GTF_MATCH_ANY     0
#define GTF_MATCH_EXACT   1
#define GTF_MATCH_CONTAIN 2
#define GTF_MATCH_WITHIN  3
#define GTF_MATCH_START   4
#define GTF_MATCH_END     5

typedef struct {
    int32_t key;
    int32_t val;
} Attribute;

/*! @typedef
 @abstract Structure for a single GTF line
 @field	 chrom         Index into the chrom hash table
 @field  source        Index into the source hash table
 @field  feature       Index into the feature hash table
 @field  start         0-based starting position
 @field  end           1-based end position
 @field  score         The score field. A value of DBL_MAX indicates a "."
 @field  strand        0: '+'; 1: '-'; 3: '.'
 @field  frame         0: '0'; 1: '1'; 2: '2'; 3: '.'
 @field  gene_id       Index into the gene_id hash table
 @field  transcript_id Index into the transcript_id hash table
 @discussion The following fields are not yet included:
    * score
    * Attributes other than gene_id and transcript_id
 Positions are 0-based half open ([start, end)), like BED files.
*/

typedef struct GTFentry {
    int32_t chrom;
    int32_t source;
    int32_t feature;
    uint32_t start;
    uint32_t end;
    double score;
    uint8_t strand:4, frame:4;
    int32_t gene_id;
    int32_t transcript_id;
    int nAttributes;
    Attribute **attrib;
    struct GTFentry *left, *right;
} GTFentry;

typedef struct {
    kstring_t chrom;
    kstring_t feature;
    kstring_t source;
    uint32_t start;
    uint32_t end;
    double score;
    uint8_t strand:4, frame: 4;
    kstring_t gene;
    kstring_t transcript;
    int nAttributes;
    Attribute **attrib;
} GTFline;

typedef struct GTFnode {
    uint32_t center;
    GTFentry *starts, *ends;
    struct GTFnode *left, *right;
} GTFnode;

typedef struct {
    int32_t chrom;
    uint32_t n_entries;
    void **tree;
} GTFchrom;

typedef struct hashTableElement {
    int32_t val;
    struct hashTableElement *next;
} hashTableElement;

typedef struct {
    uint64_t l, m;
    hashTableElement **elements;
    char **str;
} hashTable;

typedef struct {
    int32_t n_targets, m;
    int balanced;
    hashTable *htChroms;
    hashTable *htSources;
    hashTable *htFeatures;
    hashTable *htGenes;
    hashTable *htTranscripts;
    hashTable *htAttributes;
    GTFchrom **chroms;
} GTFtree;

typedef struct {
    int32_t l, m;
    GTFentry **overlaps;
} overlapSet;

//gtf.c
GTFtree * initGTFtree();
void destroyGTFtree(GTFtree *t);
void addGTFentry(GTFtree *t, GTFline *l);
void sortGTF(GTFtree *o);
void printGTFtree(GTFtree *t);
void printGTFvineStart(GTFentry *e, const char *chrom, const char *str);
void printGTFvineStartR(GTFentry *e, const char *chrom, const char *str);
void GTFline_reset(GTFline *l);
void destroyGTFline(GTFline *l);
GTFline *initGTFline();

//parseBED
GTFtree *BED2Tree(char *fname);
void GTFEntry2BED(kstring_t *ks, GTFtree *t, GTFentry *e, int ncols);

//parseGTF
GTFtree *GTF2Tree(char *fname);
void GTFEntry2GTF(kstring_t *ks, GTFtree *t, GTFentry *e);

//hashTable.c
hashTable *initHT(uint64_t size);
void destroyHTelement(hashTableElement *e);
void destroyHT(hashTable *ht);
int32_t addHTelement(hashTable *ht, char *s);
uint64_t hashString(char *s);
int strExistsHT(hashTable *ht, char *s);
int32_t str2valHT(hashTable *ht, char *s);
char *val2strHT(hashTable *ht, int32_t val);
int hasAttribute(GTFtree *t, GTFentry *e, char *str);
char *getAttribute(GTFtree *t, GTFentry *e, char *str); //NULL if the attribute isn't there

//findOverlaps.c
overlapSet *os_init();
void os_reset(overlapSet *os);
void os_destroy(overlapSet *os);
overlapSet *os_grow(overlapSet *os);
void os_exclude(overlapSet *os, int i);
overlapSet * findOverlaps(overlapSet *os, GTFtree *t, char *chrom, uint32_t start, uint32_t end, int strand, int matchType, int strandType);
int32_t countOverlaps(GTFtree *t, char *chrom, uint32_t start, uint32_t end, int strand, int matchType, int strandType);
int overlapsAny(GTFtree *t, char *chrom, uint32_t start, uint32_t end, int strand, int matchType, int strandType);

//misc.c
char *nextField(char *str);
//The following MUST be called before addGTFAttributes() or addAttribute()!
int initGTFre();
void addAttribute(GTFline *l, GTFtree *t, char *key, char *value);
int addGTFAttributes(GTFline *l, GTFtree *t, char *str);
//Call this after the above are done
void destroyGTFre();
