path= /usr/local/bin #This can be changed
CC = gcc
AR = ar
RANLIB = ranlib
OPTS = -Wall -g

.PHONY: all clean htslib install clean-all check

.SUFFIXES:.c .o

all: lib

OBJS = murmur3.o hashTable.o gtf.o findOverlaps.o misc.o parseBED.o parseGTF.o parseRMSK.o cntTable.o
VERSION = 0.1.0

#If we're building from a git repo, then append the most recent tag
ifneq "$(wildcard .git)" ""
VERSION := $(VERSION)-$(shell git describe --always --dirty)
endif

version.h:
	echo '#define libGTF_VERSION "$(VERSION)"' > $@

.c.o:
	$(CC) -c $(OPTS) -Ihtslib $< -o $@

htslib: 
	$(MAKE) -C htslib

libGTF.a: $(OBJS)
	-@rm -f $@
	$(AR) -rcs $@ $(OBJS)

lib: htslib libGTF.a

tests/testBED: lib
	$(CC) $(OPTS) -Ihtslib -o tests/testBED tests/testBED.c libGTF.a htslib/libhts.a -lz -lpthread -lpcre -lm

tests/testGTF: lib
	$(CC) $(OPTS) -Ihtslib -o tests/testGTF tests/testGTF.c libGTF.a htslib/libhts.a -lz -lpthread -lpcre -lm

tests/testFindOverlaps: lib
	$(CC) $(OPTS) -Ihtslib -o tests/testFindOverlaps tests/testFindOverlaps.c libGTF.a htslib/libhts.a -lz -lpthread -lpcre -lm

tests/testCountOverlaps: lib
	$(CC) $(OPTS) -Ihtslib -o tests/testCountOverlaps tests/testCountOverlaps.c libGTF.a htslib/libhts.a -lz -lpthread -lpcre -lm

tests/testOverlapsAny: lib
	$(CC) $(OPTS) -Ihtslib -o tests/testOverlapsAny tests/testOverlapsAny.c libGTF.a htslib/libhts.a -lz -lpthread -lpcre -lm

check: tests/testBED tests/testGTF tests/testFindOverlaps tests/testCountOverlaps tests/testOverlapsAny
	cd tests && ./test.sh && cd ..

clean:
	rm -f *.o *.a tests/testBED tests/testGTF tests/testFindOverlaps tests/testCountOverlaps tests/testOverlapsAny

clean-all: clean
	make --directory=htslib clean
