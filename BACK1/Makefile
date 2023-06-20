DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

CC = gcc

ALL = Gindex Amerge Ishow

all: $(ALL)

libfastk.c: gene_core.c gene_core.h
libfastk.h: gene_core.h

DB.c: gene_core.c gene_core.h
DB.h: gene_core.h

Gindex : Gindex.c MSDsort.c libfastk.c libfastk.h DB.c DB.h
	$(CC) $(CFLAGS) -DLCPs -o Gindex Gindex.c MSDsort.c libfastk.c DB.c gene_core.c -lpthread -lm

Amerge: Amerge.c libfastk.c libfastk.h DB.c DB.h RSDsort.c
	$(CC) $(CFLAGS) -o Amerge Amerge.c RSDsort.c libfastk.c DB.c gene_core.c -lpthread -lm

Ishow: Ishow.c libfastk.c libfastk.h gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o Ishow Ishow.c libfastk.c gene_core.c -lpthread -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f FastMGA.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf FastMGA.tar.gz LICENSE README.md Makefile *.h *.c
