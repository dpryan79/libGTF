path= /usr/local/bin #This can be changed
CC = gcc
AR = ar
RANLIB = ranlib
OPTS = -Wall -g

.PHONY: all clean htslib install clean-all

.SUFFIXES:.c .o

all: lib testBED testGTF

OBJS = murmur3.o hashTable.o gtf.o findOverlaps.o misc.o parseBED.o parseGTF.o
#VERSION = 0.1.0

##If we're building from a git repo, then append the most recent tag
#ifneq "$(wildcard .git)" ""
#VERSION := $(VERSION)-$(shell git describe --always --dirty)
#endif

#version.h:
#	echo '#define VERSION "$(VERSION)"' > $@

.c.o:
	$(CC) -c $(OPTS) -I../htslib $< -o $@

htslib: 
	$(MAKE) -C ../htslib

libGTF.a: $(OBJS)
	-@rm -f $@
	$(AR) -rcs $@ $(OBJS)

lib: libGTF.a

testBED: lib
	$(CC) $(OPTS) -I../htslib -o tests/testBED tests/testBED.c libGTF.a ../htslib/libhts.a -lz -lpthread -lpcre -lm

testGTF: lib
	$(CC) $(OPTS) -I../htslib -o tests/testGTF tests/testGTF.c libGTF.a ../htslib/libhts.a -lz -lpthread -lpcre -lm

clean:
	rm -f *.o *.a tests/testBED tests/testGTF

clean-all: clean
	make --directory=../htslib clean
