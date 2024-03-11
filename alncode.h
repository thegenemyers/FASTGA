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

static char *alnSchemaText =
  "1 3 def 1 0                 schema for aln and FastGA\n"
  ".\n"
  "P 3 seq                     SEQUENCE\n"
  "G s 3 3 INT 3 INT 6 STRING  count, length, and id for group of sequences = a scaffold\n"
  "D n 2 4 CHAR 3 INT          non-acgt chars outside (between) sequences within scaffold\n"
  "O S 1 3 DNA                 sequence\n"
  "D I 1 6 STRING              identifier of sequence\n"
  ".\n"
  "P 3 aln                     ALIGNMENTS\n"
  "D t 1 3 INT                 trace point spacing in a - global\n"
  "G a 0                       a colinear group of alignments (chain)\n"
  "D p 2 3 INT 3 INT           spacing in a,b between end of previous alignment and start of next\n"
  ".                           alignment: a_read[beg..end] b_read[beg..end], 0-indexing\n"
  "O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT\n"
  "D L 2 3 INT 3 INT           lengths of sequences a and b\n" 
  "D R 0                       flag: reverse-complement sequence b\n"
  "D Q 1 3 INT                 quality: alignment confidence in phred units\n"
  "D M 1 3 INT                 match: number of matching bases\n"
  "D D 1 3 INT                 differences: number of diffs = substitions + indels\n"
  "D C 1 6 STRING              cigar string: encodes precise alignment\n"
  "D T 1 8 INT_LIST            trace points in b\n"
  "D X 1 8 INT_LIST            number of differences in alignment per trace interval\n"
;

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
void Write_Aln_Trace  (OneFile *of, uint8 *trace, int tlen);

// end of file
