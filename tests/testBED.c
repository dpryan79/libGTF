#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include "../gtf.h"

void usage() {
    fprintf(stderr, "Usage: testBED <file.bed>\n");
    fprintf(stderr, "\n"
"This program will write the resulting interval tree representation of a BED\n"
"file in DOT format to stdout\n"
);
}

int main(int argc, char *argv[]) {
    GTFtree *t = NULL;
    char c;

    opterr = 0; //Disable error messages
    while((c = getopt(argc, argv, "h")) >= 0) {
        switch(c) {
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

    if(argc-optind != 1) {
        usage();
        return 0;
    }

    t = BED2Tree(argv[optind], NULL);
    if(!t) {
        fprintf(stderr, "Couldn't open %s or there was a problem parsing it.\n", argv[1]);
        return 1;
    }

    printGTFtree(t);
    sortGTF(t);
    printGTFtree(t);
    destroyGTFtree(t);

    return 0;
}
