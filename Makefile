DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

CC = gcc

ALL = FAtoGDB GDBtoFA GDBstat GDBshow GIXmake GIXshow GIXrm GIXmv GIXcp FastGA ALNshow ALNtoPAF ALNtoPSL ALNreset ALNplot ONEview

all: $(ALL)

libfastk.c: gene_core.c gene_core.h
libfastk.h: gene_core.h

GDB.c: gene_core.c gene_core.h
GDB.h: gene_core.h

FAtoGDB: FAtoGDB.c GDB.c GDB.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o FAtoGDB FAtoGDB.c GDB.c gene_core.c ONElib.c -lm -lz

GDBtoFA: GDBtoFA.c GDB.c GDB.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o GDBtoFA GDBtoFA.c GDB.c gene_core.c ONElib.c -lm -lz

GDBstat: GDBstat.c GDB.c GDB.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o GDBstat GDBstat.c GDB.c gene_core.c ONElib.c -lpthread -lm -lz

GDBshow: GDBshow.c GDB.h GDB.c select.c select.h hash.c hash.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o GDBshow GDBshow.c ONElib.c GDB.c select.c hash.c gene_core.c -lpthread -lm -lz

GIXmake: GIXmake.c MSDsort.c libfastk.c libfastk.h ONElib.c ONElib.h GDB.c GDB.h
	$(CC) $(CFLAGS) -DLCPs -o GIXmake GIXmake.c MSDsort.c libfastk.c ONElib.c GDB.c gene_core.c -lpthread -lm -lz

GIXshow: GIXshow.c libfastk.c libfastk.h gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o GIXshow GIXshow.c libfastk.c gene_core.c -lpthread -lm

GIXrm: GIXrm.c gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o GIXrm GIXrm.c gene_core.c -lm

GIXmv: GIXxfer.c GDB.c GDB.h gene_core.c ONElib.c ONElib.h gene_core.h
	$(CC) $(CFLAGS) -DMOVE -o GIXmv GIXxfer.c GDB.c ONElib.c gene_core.c -lm -lz

GIXcp: GIXxfer.c GDB.c GDB.h ONElib.c ONElib.h gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o GIXcp GIXxfer.c GDB.c ONElib.c gene_core.c -lm -lz

FastGA: FastGA.c libfastk.c libfastk.h GDB.c GDB.h RSDsort.c align.c align.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o FastGA FastGA.c RSDsort.c libfastk.c align.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

ALNshow: ALNshow.c align.h align.c GDB.c GDB.h select.c select.h hash.c hash.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ALNshow ALNshow.c align.c GDB.c alncode.c select.c hash.c gene_core.c ONElib.c -lpthread -lm -lz

ALNtoPAF: ALNtoPAF.c align.h align.c GDB.c GDB.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ALNtoPAF ALNtoPAF.c align.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

ALNtoPSL: ALNtoPSL.c align.h align.c GDB.c GDB.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ALNtoPSL ALNtoPSL.c align.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

ALNreset: ALNreset.c GDB.c GDB.h ONElib.c ONElib.h alncode.c alncode.h
	$(CC) $(CFLAGS) -o ALNreset ALNreset.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

ALNplot: ALNplot.c hash.c hash.h select.c select.h GDB.c GDB.h ONElib.c ONElib.h alncode.c alncode.h
	$(CC) $(CFLAGS) -o ALNplot ALNplot.c GDB.c alncode.c select.c hash.c gene_core.c ONElib.c -lpthread -lm -lz

ONEview: ONEview.c ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ONEview ONEview.c ONElib.c -lm -lz

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f FastGA.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf FastGA.tar.gz LICENSE README.md Makefile *.h *.c
