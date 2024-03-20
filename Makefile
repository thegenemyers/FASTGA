DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
CFLAGS = -g -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

CC = gcc

ALL = FAtoGDB GDBtoFA GDBstat GDBshow GIXmake GIXshow GIXrm GIXmv GIXcp FastGA ALNshow ALNtoPAF ALNtoPSL ALNreset

all: $(ALL)

libfastk.c: gene_core.c gene_core.h
libfastk.h: gene_core.h

DB.c: gene_core.c gene_core.h
DB.h: gene_core.h

FAtoGDB: FAtoGDB.c DB.c DB.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o FAtoGDB FAtoGDB.c DB.c gene_core.c alncode.c ONElib.c -lm -lz

GDBtoFA: GDBtoFA.c DB.c DB.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o GDBtoFA GDBtoFA.c DB.c gene_core.c alncode.c ONElib.c -lm -lz

GDBstat: GDBstat.c DB.c DB.h
	$(CC) $(CFLAGS) -o GDBstat GDBstat.c DB.c gene_core.c -lpthread -lm

GDBshow: GDBshow.c libfastk.c libfastk.h DB.c DB.h
	$(CC) $(CFLAGS) -o GDBshow GDBshow.c libfastk.c DB.c gene_core.c -lpthread -lm

GIXmake: GIXmake.c MSDsort.c libfastk.c libfastk.h DB.c DB.h
	$(CC) $(CFLAGS) -DLCPs -o GIXmake GIXmake.c MSDsort.c libfastk.c DB.c gene_core.c -lpthread -lm

GIXshow: GIXshow.c libfastk.c libfastk.h gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o GIXshow GIXshow.c libfastk.c gene_core.c -lpthread -lm

GIXrm: GIXrm.c gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o GIXrm GIXrm.c gene_core.c -lm

GIXmv: GIXxfer.c gene_core.c gene_core.h
	$(CC) $(CFLAGS) -DMOVE -o GIXmv GIXxfer.c gene_core.c -lm

GIXcp: GIXxfer.c gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o GIXcp GIXxfer.c gene_core.c -lm

FastGA: FastGA.c libfastk.c libfastk.h DB.c DB.h RSDsort.c align.c align.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o FastGA FastGA.c RSDsort.c libfastk.c align.c DB.c alncode.c gene_core.c ONElib.c -lpthread -lm

ALNshow: ALNshow.c align.h align.c DB.c DB.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ALNshow ALNshow.c align.c DB.c alncode.c gene_core.c ONElib.c -lpthread -lm

ALNtoPAF: ALNtoPAF.c align.h align.c DB.c DB.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ALNtoPAF ALNtoPAF.c align.c DB.c alncode.c gene_core.c ONElib.c -lpthread -lm

ALNtoPSL: ALNtoPSL.c align.h align.c DB.c DB.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ALNtoPSL ALNtoPSL.c align.c DB.c alncode.c gene_core.c ONElib.c -lpthread -lm

ALNreset: ALNreset.c DB.c DB.h ONElib.c ONElib.h alncode.c alncode.h
	$(CC) $(CFLAGS) -o ALNreset ALNreset.c DB.c alncode.c gene_core.c ONElib.c -lpthread -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f FastGA.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf FastGA.tar.gz LICENSE README.md Makefile *.h *.c
