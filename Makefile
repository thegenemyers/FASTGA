DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

CC = gcc

ALL = FAtoGDB GDBtoFA GDBstat GDBshow GIXmake GIXshow GIXrm GIXmv GIXcp FastGA LASshow LAStoPSL LAStoPAF LAStoONE LASreset

all: $(ALL)

libfastk.c: gene_core.c gene_core.h
libfastk.h: gene_core.h

DB.c: gene_core.c gene_core.h
DB.h: gene_core.h

FAtoGDB: FAtoGDB.c DB.c DB.h
	$(CC) $(CFLAGS) -o FAtoGDB FAtoGDB.c DB.c gene_core.c -lm -lz

GDBtoFA: GDBtoFA.c DB.c DB.h
	$(CC) $(CFLAGS) -o GDBtoFA GDBtoFA.c DB.c gene_core.c -lm -lz

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

FastGA: FastGA.c libfastk.c libfastk.h DB.c DB.h RSDsort.c align.c align.h
	$(CC) $(CFLAGS) -o FastGA FastGA.c RSDsort.c libfastk.c align.c DB.c gene_core.c -lpthread -lm

LASshow: LASshow.c align.h align.c DB.c DB.h
	$(CC) $(CFLAGS) -o LASshow LASshow.c align.c DB.c gene_core.c -lpthread -lm

LAStoPSL: LAStoPSL.c align.h align.c DB.c DB.h
	$(CC) $(CFLAGS) -o LAStoPSL LAStoPSL.c align.c DB.c gene_core.c -lpthread -lm

LAStoPAF: LAStoPAF.c align.h align.c DB.c DB.h
	$(CC) $(CFLAGS) -o LAStoPAF LAStoPAF.c align.c DB.c gene_core.c -lpthread -lm

LAStoONE: LAStoONE.c align.c align.h DB.c DB.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o LAStoONE LAStoONE.c align.c DB.c gene_core.c ONElib.c -lm

LASreset: LASreset.c DB.c DB.h
	$(CC) $(CFLAGS) -o LASreset LASreset.c DB.c gene_core.c -lpthread -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f FastGA.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf FastGA.tar.gz LICENSE README.md Makefile *.h *.c
