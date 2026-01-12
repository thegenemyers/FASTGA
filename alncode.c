/*  File: alncode.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2024
 *-------------------------------------------------------------------
 * Description: IO for ONEcode .1aln files for Myers FASTGA package
 * Exported functions:
 * HISTORY:
 * Last edited: Jul 15 13:10 2024 (rd109)
 * Created: Sat Feb 24 12:19:16 2024 (rd109)
 *-------------------------------------------------------------------
 */

#include <string.h>
#include <stdlib.h>

#include "alncode.h"
#include "GDB.h"

static char *alnSchemaText =
  "1 3 def 2 1                 schema for aln and FastGA\n"
  ".\n"
  "P 3 seq                     SEQUENCE\n"
  "O s 2 3 INT 6 STRING        length and id for group of sequences = a scaffold\n"
  "G S                         scaffolds (s) group sequence objects (S)\n"
  "D n 2 4 CHAR 3 INT          non-acgt chars outside (between) sequences within scaffold\n"
  "O S 1 3 DNA                 sequence\n"
  "D I 1 6 STRING              identifier of sequence\n"
  ".\n"
  "P 3 aln                     ALIGNMENTS\n"
  "D t 1 3 INT                 trace point spacing in a - global\n"
  ".                           GDB skeleton (may not be presend)\n"
  "O g 0                       groups scaffolds into a GDB skeleton\n"
  "G S                         collection of scaffolds constituting a GDB\n"
  "O S 1 6 STRING              id for a scaffold\n"
  "D G 1 3 INT                 gap of given length\n"
  "D C 1 3 INT                 contig of given length\n"
  ".\n"
  "O a 0                       groups A's into a colinear chain\n"
  "G A                         chains (a) group alignment objects (A)\n"
  "D p 2 3 INT 3 INT           spacing in a,b between end of previous alignment and start of next\n"
  ".                           alignment: a_read[beg..end] b_read[beg..end], 0-indexing\n"
  "O A 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT\n"
  "D L 2 3 INT 3 INT           lengths of sequences a and b\n"
  "D R 0                       flag: reverse-complement sequence b\n"
  "D D 1 3 INT                 differences: number of diffs = substitions + indels\n"
  "D T 1 8 INT_LIST            trace points in b\n"
  "D X 1 8 INT_LIST            number of differences in alignment per trace interval\n"
  "D Q 1 3 INT                 quality: alignment confidence in phred units (currently unused)\n"
  "D E 1 3 INT                 match: number of equal bases (currently unused)\n"
  "D Z 1 6 STRING              cigar string: encodes precise alignment (currently unused)\n"
  "D U 1 3 INT                 putative unit size of a TR alignment (FASTAN)\n"
;

OneSchema *make_Aln_Schema ()
{ return (oneSchemaCreateFromText(alnSchemaText)); }

  // Open the .1aln file for reading and read the header

OneFile *open_Aln_Read (char *filename, int nThreads,
			int64 *nOverlaps, int *tspace,
			char **db1_name, char **db2_name, char **cpath)
{ OneSchema *schema;
  OneFile   *of;
  OneInfo   *refInfo;
  int        i;

  schema = oneSchemaCreateFromText(alnSchemaText);
  if (schema == NULL) 
    { EPRINTF("Failed to create 1aln schema");
      EXIT(NULL);
    }

  of = oneFileOpenRead(filename,schema,"aln",nThreads);
  if (of == NULL)
    { EPRINTF("Failed to open .1aln file %s",filename);
      oneSchemaDestroy(schema);
      EXIT(NULL);
    }

  if (of->info['A'] == NULL)
    { EPRINTF("No alignments found in aln file");
      goto clean_up;
    }
  *nOverlaps = of->info['A']->given.count;
    
  *db1_name = NULL;
  *db2_name = NULL;
  *cpath    = NULL;
  refInfo = of->info['<'];
  if (refInfo == NULL)
    { EPRINTF("No references in aln file");
      goto clean_up;
    }      
  for (i = 0; i < refInfo->accum.count; ++i)
    if (of->reference[i].count == 1)
      { if (*db1_name != NULL)
          free(*db1_name);
        *db1_name = strdup(of->reference[i].filename);
      }
    else if (of->reference[i].count == 2)
      { if (*db2_name != NULL)
          free(*db2_name);
        *db2_name = strdup(of->reference[i].filename);
      }
    else if (of->reference[i].count == 3)
      { if (*cpath != NULL)
          free(*cpath);
        *cpath = strdup(of->reference[i].filename);
      }
  if (cpath == NULL)
    *cpath = "";

  *tspace = 0;
  while (oneReadLine(of))             // advance to first alignment record
    if (of->lineType == 'A' || of->lineType == 'g')
      break;
    else if (of->lineType == 't')
      *tspace = oneInt(of,0);

  if (*tspace == 0)
    { EPRINTF("Did not find a t-line before first alignment or GDB skeleton");
      goto clean_up;
    }      

  oneSchemaDestroy(schema);
  return (of);

clean_up:
  oneFileClose(of);
  oneSchemaDestroy(schema);
  EXIT(NULL);
}

  // Next two routines read the records from the file

int Read_Aln_Overlap(OneFile *of, Overlap *ovl)
{ if (of->lineType != 'A')
    { EPRINTF("Failed to be at start of alignment in Read_Aln_Overlap()");
      EXIT(1);
    }
    
  ovl->flags = 0;
  ovl->aread = oneInt(of,0);
  ovl->path.abpos = oneInt(of,1);
  ovl->path.aepos = oneInt(of,2);
  ovl->bread = oneInt(of,3);
  ovl->path.bbpos = oneInt(of,4);
  ovl->path.bepos = oneInt(of,5);

  while (oneReadLine(of))
    if (of->lineType == 'T')
       break;
    else if (of->lineType == 'R')
      ovl->flags |= COMP_FLAG;
    else if (of->lineType == 'D')
      ovl->path.diffs = oneInt(of,0);
    else if (of->lineType == 'A')
       break;

  if (of->lineType != 'T')
    { EPRINTF("Failed to find trace record in .1aln object %lld",
                     of->info['A']->accum.count);
      EXIT(1);
    }
  return (0);
}

int Read_Aln_Trace(OneFile *of, uint8 *trace, int *period)
{ int64 *trace64;
  int    tlen;
  int    j, x;
  
  if (of->lineType != 'T')
    { EPRINTF("Failed to be at start of trace in Read_Aln_Trace()");
      EXIT(1);
    }
    
  tlen    = 2*oneLen(of);
  trace64 = oneIntList(of);
  j = 0;
  for (x = 1; x < tlen; x += 2)
    trace[x] = trace64[j++];

  oneReadLine(of);
  if (of->lineType != 'X')
    { EPRINTF("No X-line following a T-line in 1aln file");
      EXIT(1);
    }
  if (2*oneLen(of) != tlen)
    { EPRINTF("X-line and T-lines should have the same length");
      EXIT(1);
    }

  trace64 = oneIntList(of);
  j = 0;
  for (x = 0; x < tlen; x += 2)
    trace[x] = trace64[j++];

  if (period != NULL)
    *period = 0;
  while (oneReadLine(of))       // move to start of next alignment
    if (of->lineType == 'A')
      break;
    else if (of->lineType == 'U' && period != NULL)
      *period = oneInt(of,0);

  return (tlen);
}

int Skip_Aln_Trace(OneFile *of)
{ int    tlen;
  
  if (of->lineType != 'T')
    { EPRINTF("Failed to be at start of trace in Read_Aln_Trace()");
      EXIT(1);
    }
    
  tlen = oneLen(of);

  oneReadLine(of);
  if (of->lineType != 'X')
    { EPRINTF("No X-line following a T-line in 1aln file");
      EXIT(1);
    }
  if (oneLen(of) != tlen)
    { EPRINTF("X-line and T-lines should have the same length");
      EXIT(1);
    }

  while (oneReadLine(of))       // move to start of next alignment
    if (of->lineType == 'A')
      break;

  return (0);
}

  // And these routines write an alignment

OneFile *open_Aln_Write (char *filename, int nThreads,
			 char *progname, char *version, char *commandLine, int tspace,
			 char *db1_name, char *db2_name, char *cpath)
{ OneSchema *schema;
  OneFile   *of;

  schema = oneSchemaCreateFromText(alnSchemaText);
  if (schema == NULL) 
    { EPRINTF("Failed to create 1aln schema");
      EXIT(NULL);
    }

  of = oneFileOpenWriteNew(filename,schema,"aln",true,nThreads);
  if (of == NULL)
    { EPRINTF("Failed to open .1aln file %s",filename);
      EXIT(NULL);
    }

  oneAddProvenance(of,progname,version,commandLine);
  
  oneAddReference(of,db1_name,1);
  if (db2_name != NULL)
    oneAddReference(of,db2_name,2);
  if (cpath)
    oneAddReference(of,cpath,3);

  oneInt(of,0) = tspace;
  oneWriteLine (of,'t',0,0);

  oneSchemaDestroy (schema);
  return of;
}

void Write_Aln_Overlap (OneFile *of, Overlap *ovl)
{ oneInt(of,0) = ovl->aread;
  oneInt(of,1) = ovl->path.abpos;
  oneInt(of,2) = ovl->path.aepos;
  oneInt(of,3) = ovl->bread;
  oneInt(of,4) = ovl->path.bbpos;
  oneInt(of,5) = ovl->path.bepos;
  oneWriteLine (of, 'A', 0, 0);

  if (COMP(ovl->flags))
    oneWriteLine(of,'R',0,0);

  oneInt(of,0) = ovl->path.diffs;
  oneWriteLine (of,'D',0,0);
}

void Write_Aln_Trace (OneFile *of, uint8 *trace, int tlen, int64 *trace64, int period)
{ int j, x;

  j = 0;
  for (x = 1; x < tlen; x += 2)
    trace64[j++] = trace[x];
  oneWriteLine (of,'T',j,trace64);

  j = 0;
  for (x = 0; x < tlen; x += 2)
    trace64[j++] = trace[x];
  oneWriteLine(of,'X',j,trace64);

  if (period != 0)
    { oneInt(of,0) = period;
      oneWriteLine(of,'U',0,NULL);
    }
}

int Copy_Aln_Trace(OneFile *in, OneFile *out)
{ int    tlen;
  
  if (in->lineType != 'T')
    { EPRINTF("Failed to be at start of trace in Read_Aln_Trace()");
      EXIT(1);
    }
    
  tlen = oneLen(in);
  oneWriteLine(out,'T',tlen,oneIntList(in));

  oneReadLine(in);
  if (in->lineType != 'X')
    { EPRINTF("No X-line following a T-line in 1aln file");
      EXIT(1);
    }
  if (oneLen(in) != tlen)
    { EPRINTF("X-line and T-lines should have the same length");
      EXIT(1);
    }
  oneWriteLine(out,'X',tlen,oneIntList(in));

  while (oneReadLine(in))       // move to start of next alignment
    if (in->lineType == 'U')
      { oneInt(out,0) = oneInt(in,0);
        oneWriteLine(out,'U',0,NULL);
      }
    else if (in->lineType == 'A')
      break;

  return (0);
}
