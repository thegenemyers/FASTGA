/*  File: alncode.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: IO for ONEcode .1aln files for Myers FASTGA package
 * Exported functions:
 * HISTORY:
 * Last edited: Feb 25 15:28 2024 (rd109)
 * Created: Sat Feb 24 12:19:16 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include "ONElib.h"
#include "align.h"

OneSchema *make_Aln_Schema();

// open the .1aln file for reading and read the header

OneFile *open_Aln_Read (char *filename, int nThreads,
			int64 *nOverlaps, int *tspace,
			char **db1_name, char **db2_name, char **cpath) ;

// next two routines read the records from the file

void Read_Aln_Overlap(OneFile *of, Overlap *ovl);
int  Read_Aln_Trace  (OneFile *of, uint8 *trace);
void Skip_Aln_Trace  (OneFile *of);

// and equivalents for writing

OneFile *open_Aln_Write (char *filename, int nThreads,
			 char *progname, char *version, char *commandLine,
			 int tspace, char *db1_name, char *db2_name, char *cpath);

void Write_Aln_Overlap(OneFile *of, Overlap *ovl);
void Write_Aln_Trace  (OneFile *of, uint8 *trace, int tlen, int64 *trace64);

// end of file
