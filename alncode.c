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
  "O a 1 3 INT                 number of alignments for colinear groups = a chain\n"
  "G A                         chains (a) group alignment objects (A)\n"
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
    { fprintf (stderr, "failed to create 1aln schema\n");
      return (NULL);
    }

  of = oneFileOpenRead(filename,schema,"aln",nThreads);
  if (of == NULL)
    { fprintf (stderr,"%s: Failed to open .1aln file %s\n",Prog_Name,filename);
      oneSchemaDestroy(schema);
      return (NULL);
    }

  if (of->info['A'] == NULL)
    { fprintf(stderr,"%s: No alignments found in aln file\n",Prog_Name);
      goto clean_up;
    }
  *nOverlaps = of->info['A']->given.count;
    
  *db1_name = NULL;
  *db2_name = NULL;
  *cpath    = "";
  refInfo = of->info['<'];
  if (refInfo == NULL)
    { fprintf(stderr,"%s: No references in aln file",Prog_Name);
      goto clean_up;
    }      
  for (i = 0; i < refInfo->accum.count; ++i)
    if (of->reference[i].count == 1)
      *db1_name = strdup(of->reference[i].filename);
    else if (of->reference[i].count == 2)
      *db2_name = strdup(of->reference[i].filename);
    else if (of->reference[i].count == 3)
      *cpath = strdup(of->reference[i].filename);

  *tspace = 0;
  while (oneReadLine(of))             // advance to first alignment record
    if (of->lineType == 'A')
      break;
    else if (of->lineType == 't')
      *tspace = oneInt(of,0);

  if (*tspace == 0)
    { fprintf(stderr,"%s: Did not find a t-line before first alignment\n",Prog_Name);
      goto clean_up;
    }      

  oneSchemaDestroy(schema);
  return (of);

clean_up:
  oneFileClose(of);
  oneSchemaDestroy(schema);
  return (NULL);
}

  // Next two routines read the records from the file

void Read_Aln_Overlap(OneFile *of, Overlap *ovl)
{ if (of->lineType != 'A')
    { fprintf(stderr,"%s: Failed to be at start of alignment in Read_Aln_Overlap()\n",Prog_Name);
      exit (1);
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
    { fprintf(stderr,"%s: Failed to find trace record in .1aln object %lld\n",
                     Prog_Name,of->info['A']->accum.count);
      exit (1);
    }
}

int Read_Aln_Trace(OneFile *of, uint8 *trace)
{ int64 *trace64;
  int    tlen;
  int    j, x;
  
  if (of->lineType != 'T')
    { fprintf(stderr,"%s: Failed to be at start of trace in Read_Aln_Trace()\n",Prog_Name);
      exit (1);
    }
    
  tlen    = 2*oneLen(of);
  trace64 = oneIntList(of);
  j = 0;
  for (x = 1; x < tlen; x += 2)
    trace[x] = trace64[j++];

  oneReadLine(of);
  if (of->lineType != 'X')
    { fprintf(stderr,"%s: No X-line following a T-line in 1aln file\n",Prog_Name);
      exit (1);
    }
  if (2*oneLen(of) != tlen)
    { fprintf(stderr,"%s: X-line and T-lines should have the same length\n",Prog_Name);
      exit (1);
    }

  trace64 = oneIntList(of);
  j = 0;
  for (x = 0; x < tlen; x += 2)
    trace[x] = trace64[j++];

  while (oneReadLine(of))       // move to start of next alignment
    if (of->lineType == 'A')
      break;

  return (tlen);
}

void Skip_Aln_Trace(OneFile *of)
{ int    tlen;
  
  if (of->lineType != 'T')
    { fprintf(stderr,"%s: Failed to be at start of trace in Read_Aln_Trace()\n",Prog_Name);
      exit (1);
    }
    
  tlen = oneLen(of);

  oneReadLine(of);
  if (of->lineType != 'X')
    { fprintf(stderr,"%s: No X-line following a T-line in 1aln file\n",Prog_Name);
      exit (1);
    }
  if (oneLen(of) != tlen)
    { fprintf(stderr,"%s: X-line and T-lines should have the same length\n",Prog_Name);
      exit (1);
    }

  while (oneReadLine(of))       // move to start of next alignment
    if (of->lineType == 'A')
      break;
}

  // And these routines write an alignment

OneFile *open_Aln_Write (char *filename, int nThreads,
			 char *progname, char *version, char *commandLine,
			 int tspace, char *db1_name, char *db2_name, char *cpath)
{ OneSchema *schema;
  OneFile   *of;

  schema = oneSchemaCreateFromText(alnSchemaText);
  if (schema == NULL) 
    { fprintf(stderr,"%s: Failed to create 1aln schema\n",Prog_Name);
      return 0;
    }

  of = oneFileOpenWriteNew(filename,schema,"aln",true,nThreads);
  if (of == NULL)
    { fprintf(stderr,"%s: Failed to open .1aln file %s\n",Prog_Name,filename);
      return 0;
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

void Write_Aln_Trace (OneFile *of, uint8 *trace, int tlen)
{ static int    tmax = 0;
  static int64 *trace64 = NULL;
  int    j, x;

  if (tlen > tmax)
    { tmax    = tlen*1.2 + 1024;
      trace64 = (int64 *) Realloc(trace64,(tmax/2)*sizeof(int64),"Reallocating trace vector");
    }

  j = 0;
  for (x = 1; x < tlen; x += 2)
    trace64[j++] = trace[x];
  oneWriteLine (of,'T',j,trace64);

  j = 0;
  for (x = 0; x < tlen; x += 2)
    trace64[j++] = trace[x];
  oneWriteLine(of,'X',j,trace64);
}
