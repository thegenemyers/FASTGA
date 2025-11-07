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
  "D M 1 8 INT_LIST            mask pair list for a contig\n"
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

int Read_Aln_Skeleton(OneFile *of, char *source, GDB *gdb)
{ GDB_SCAFFOLD  *scf;
  GDB_CONTIG    *ctg;
  GDB_MASK      *msk;
  OneProvenance *prov;
  char          *hdr;
  int64          nscaff, ncontig, nprov, nmasks;
  int64          len, seqtot, hdrtot, maxctg, boff, spos, iscaps;

  oneStats(of,'S',&nscaff,NULL,&hdrtot);
  oneStats(of,'C',&ncontig,NULL,NULL);
  oneStats(of,'M',NULL,NULL,&nmasks);
  hdrtot += nscaff;
  nmasks /= 2;

  scf     = Malloc(sizeof(GDB_SCAFFOLD)*nscaff,"Reading GDB Skeleton");
  ctg     = Malloc(sizeof(GDB_CONTIG)*(ncontig+1),"Reading GDB Skeleton");
  if (nmasks > 0)
    msk = Malloc(sizeof(GDB_MASK)*(nmasks+1),"Reading GDB Skeleton");
  else
    msk = NULL;
  hdr   = Malloc(hdrtot,"Reading GDB Skeleton");
  nprov = 0;
  prov  = NULL;
  if (scf == NULL || ctg == NULL || hdr == NULL || (nmasks > 0 && msk == NULL))
    { free(hdr);
      free(ctg);
      free(scf);
      free(msk);
      return (1);
    }

  gdb->freq[0] = gdb->freq[1] = gdb->freq[2] = gdb->freq[3] = .25;

  nscaff  = -1;
  ncontig = 0;
  nmasks  = 0;
  hdrtot  = 0;
  seqtot  = 0;
  maxctg  = 0;
  iscaps  = 0;
  boff = 0;
  spos = 0;
  while (oneReadLine(of))
    switch (of->lineType)
    { case 'g':
      case 'A':
        goto endofsketch;
      case 'S':
        if (nscaff >= 0)
          { scf[nscaff].ectg = ncontig;
            scf[nscaff].slen = spos;
            spos = 0;
          }
        nscaff += 1;
        scf[nscaff].hoff = hdrtot;
        scf[nscaff].fctg = ncontig;
        len = oneLen(of);
        memcpy(hdr+hdrtot,oneString(of),len);
        hdrtot += len;
        hdr[hdrtot++] = '\0';
        break;
      case 'G':
        spos += oneInt(of,0);
        break;
      case 'C':
        len = oneInt(of,0);
        ctg[ncontig].boff = boff;
        ctg[ncontig].moff = nmasks;
        ctg[ncontig].sbeg = spos;
        ctg[ncontig].clen = len;
        ctg[ncontig].scaf = nscaff;
        ncontig += 1;
        if (len > maxctg)
          maxctg = len;
        seqtot  += len;
        boff    += COMPRESSED_LEN(len);
        spos    += len;
        break;
      case 'M':
        { int    i;
          int64 *list;

          len  = oneInt(of,0);
          list = oneIntList(of);
          for (i = 0; i < len; i += 2)
            { msk[nmasks].beg = list[i];
              msk[nmasks].end = list[i+1];
              nmasks += 1;
            }
          iscaps = 1;
        }
        break;
    }
endofsketch:
  scf[nscaff].ectg = ncontig;
  scf[nscaff].slen = spos;
  nscaff += 1;

  ctg[ncontig].moff = nmasks;
  ctg[ncontig].boff = ctg[ncontig-1].boff + COMPRESSED_LEN(ctg[ncontig-1].clen);

  gdb->nprov = nprov;
  gdb->prov  = prov;

  gdb->nscaff    = nscaff;
  gdb->scaffolds = scf;

  gdb->ncontig = ncontig;
  gdb->maxctg  = maxctg;
  gdb->contigs = ctg;

  gdb->iscaps  = iscaps;
  gdb->nmasks  = nmasks;
  gdb->masks   = msk;

  gdb->srcpath  = strdup(source);
  gdb->seqpath  = NULL;

  gdb->hdrtot  = hdrtot;
  gdb->headers = hdr;

  gdb->seqtot   = seqtot;
  gdb->seqstate = EXTERNAL;
  gdb->seqs     = NULL;

  source += strlen(source);
  if (strcmp(source-3,".gz") == 0)
    gdb->seqsrc = IS_FA_GZ;
  else if (strcmp(source-3,".fa") == 0)
    gdb->seqsrc = IS_FA;
  else if (strcmp(source-4,".fna") == 0)
    gdb->seqsrc = IS_FA;
  else if (strcmp(source-6,".fasta") == 0)
    gdb->seqsrc = IS_FA;
  else
    gdb->seqsrc = IS_ONE;

  return (0);
}

void Skip_Aln_Skeletons(OneFile *of)
{ if (of->lineType == 'g')
    { while (oneReadLine(of))
        if (of->lineType == 'A')
          break;
    }
}

void Write_Aln_Skeleton(OneFile *of, GDB *gdb)
{ int64      spos, len;
  char      *head;
  int        s, c;

  oneWriteLine(of,'g',0,0);

  for (s = 0; s < gdb->nscaff; s++)
    { head = gdb->headers + gdb->scaffolds[s].hoff;
      oneWriteLine(of,'S',strlen(head),head);

      spos = 0;
      for (c = gdb->scaffolds[s].fctg; c < gdb->scaffolds[s].ectg; c++)
        { if (gdb->contigs[c].sbeg > spos)
            { oneInt(of,0) = gdb->contigs[c].sbeg - spos;
              oneWriteLine(of,'G',0,0);
            }
          len = gdb->contigs[c].clen;
          oneInt(of,0) = len;
          oneWriteLine(of,'C',0,0);
          spos = gdb->contigs[c].sbeg + len;
          if (c == gdb->ncontig-1)
            len = 2*(gdb->nmasks - gdb->contigs[c].moff);
          else
            len = 2*(gdb->contigs[c+1].moff - gdb->contigs[c].moff);
          if (len > 0)
            oneWriteLine(of,'M',len,(int64 *) (gdb->masks + gdb->contigs[c].moff));
        }
      if (gdb->scaffolds[s].slen > spos)
        { oneInt(of,0) = gdb->scaffolds[s].slen - spos;
          oneWriteLine(of,'G',0,0);
        }
    }
}

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

int Read_Aln_Trace(OneFile *of, uint8 *trace)
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

  while (oneReadLine(of))       // move to start of next alignment
    if (of->lineType == 'A')
      break;

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

void Write_Aln_Trace (OneFile *of, uint8 *trace, int tlen, int64 *trace64)
{ int j, x;

  j = 0;
  for (x = 1; x < tlen; x += 2)
    trace64[j++] = trace[x];
  oneWriteLine (of,'T',j,trace64);

  j = 0;
  for (x = 0; x < tlen; x += 2)
    trace64[j++] = trace[x];
  oneWriteLine(of,'X',j,trace64);
}
