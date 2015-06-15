#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "../gtf.h"

void usage() {
    fprintf(stderr, "Usage: testGTF <file.gtf>\n");
    fprintf(stderr, "\n"
"This program will write the resulting interval tree representation of an RMSK\n"
"file in DOT format to stdout\n");
}

int main(int argc, char *argv[]) {
    if(argc != 2) {
        usage();
        return 0;
    }
    if(strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0) {
        usage();
        return 0;
    }

    GTFtree *t = RMSK2Tree(argv[1],NULL);
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
