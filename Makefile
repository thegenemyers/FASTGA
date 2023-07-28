DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

CC = gcc

ALL = Gindex Gshow FastGA Hitter

all: $(ALL)

libfastk.c: gene_core.c gene_core.h
libfastk.h: gene_core.h

DB.c: gene_core.c gene_core.h
DB.h: gene_core.h

Gindex : Gindex.c MSDsort.c libfastk.c libfastk.h DB.c DB.h
	$(CC) $(CFLAGS) -DLCPs -o Gindex Gindex.c MSDsort.c libfastk.c DB.c gene_core.c -lpthread -lm

FastGA: FastGA.c libfastk.c libfastk.h DB.c DB.h RSDsort.c align.c align.h
	$(CC) $(CFLAGS) -o FastGA FastGA.c RSDsort.c libfastk.c align.c DB.c gene_core.c -lpthread -lm

Gshow: Gshow.c libfastk.c libfastk.h gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o Gshow Gshow.c libfastk.c gene_core.c -lpthread -lm

Hitter: Hitter.c align.h align.c DB.c DB.h
	$(CC) $(CFLAGS) -o Hitter Hitter.c align.c DB.c gene_core.c -lpthread -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f FastMGA.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf FastG.tar.gz LICENSE README.md Makefile *.h *.c
