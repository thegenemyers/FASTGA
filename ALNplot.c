/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2024 Chenxi Zhou <chnx.zhou@gmail.com>                          *
 *                                                                               *
 * Permission is hereby granted, free of charge, to any person obtaining a copy  *
 * of this software and associated documentation files (the "Software"), to deal *
 * in the Software without restriction, including without limitation the rights  *
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     *
 * copies of the Software, and to permit persons to whom the Software is         *
 * furnished to do so, subject to the following conditions:                      *
 *                                                                               *
 * The above copyright notice and this permission notice shall be included in    *
 * all copies or substantial portions of the Software.                           *
 *                                                                               *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, *
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE *
 * SOFTWARE.                                                                     *
 *********************************************************************************/

/*****************************************************************************************\
*                                                                                         *
*  Genome Alignment Plotter                                                               *
*                                                                                         *
*  Author:  Chenxi Zhou (with an assist from Gene Myers)                                  *
*  Date  :  June 2024                                                                     *
*                                                                                         *
\*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <zlib.h>

#include "GDB.h"
#include "hash.h"
#include "select.h"
#include "alncode.h"

#undef DEBUG_MAKE_HASH
#undef DEBUG_READ_1ALN
#undef DEBUG_AXIS_CONF

#define MAX_XY_LEN 10000   // Max/Min plotting width/height
#define MIN_XY_LEN 50

#define MAX_LAB_LEN 20     // Max number characters for sequence names
#define MAX_LAB_FRC .2     // Max ratio of label to line plot panel

  //  Command line syntax and global parameter variables

static char *Usage[] =
  { "[-vSL] [-T<int(4)>] [-p[:<output:path>[.pdf]]]",
    "[-l<int(100)>] [-i<float(.7)>] [-n<int(100000)>]",
    "[-H<int(600)>] [-W<int>] [-f<int>] [-t<float>]",
    "<alignment:path>[.1aln|.paf[.gz]]> [<selection>|<FILE> [<selection>|<FILE>]]",
  };

static int    VERBOSE;            // -v
// static int    HIGHLIGHT;          // -h
// static int    TRYADIAG;           // -d
static int    PRINTSID;           // -S
static int    LABELS;             // ! -L

static int    NTHREADS = 4;       // -T
static char  *OUTEPS   = NULL;    // -p

static int    MINALEN  = 100;     // -a
static double MINAIDNT = 0.7;     // -e
static int    MAXALIGN = 100000;  // -n

static int    IMGHEIGH = 0;       // -H
static int    IMGWIDTH = 0;       // -W
static int    FONTSIZE = 0;       // -f
static double LINESIZE = 0;       // -t

  //  Array of segments to plot

typedef struct
  { uint8 flag;
    int   aread, bread;
    int   abpos, bbpos;
    int   aepos, bepos;
    // int   group; // segment group - e.g. chain
  } Segment;

#define DEL_FLAG 0x1
#define COL_GRAY 0x2
#define COL_RED  0x4
#define COL_BLUE 0x8

#define IS_DEL(flag) ((flag)&DEL_FLAG)
#define SET_DEL(flag) ((flag)|=DEL_FLAG)

#define IS_COLOR(flag) ((flag)&0xE)
#define IS_RED(flag)  ((flag)&COL_RED)
#define IS_BLUE(flag) ((flag)&COL_BLUE)
#define IS_GRAY(flag) ((flag)&COL_GRAY)
#define UNSET_COLOR(flag) ((flag)&=0xF1)
#define SET_RED(flag)  (UNSET_COLOR(flag), (flag)|=COL_RED)
#define SET_BLUE(flag) (UNSET_COLOR(flag), (flag)|=COL_BLUE)
#define SET_GRAY(flag) (UNSET_COLOR(flag), (flag)|=COL_GRAY)

static Segment  *segments = NULL;
static int64     nSegment = 0;

  //  Genome skeletons and auxiliary info

static int ISTWO;                           // If two GDBs are different

static GDB           _AGDB, *AGDB = &_AGDB;  //  1st genome skeleton with parts ...
static GDB           _BGDB, *BGDB = &_BGDB;  //  2nd genome skeleton with parts ...

static int           NASCAFF;     // Number scaffolds
static int           NACONTIG;    // Number contigs
static GDB_SCAFFOLD *ASCAFFS;     // Scaffold record array
static GDB_CONTIG   *ACONTIG;     // Contig record array

static int           NBSCAFF;     // Number scaffolds
static int           NBCONTIG;    // Number contigs
static GDB_SCAFFOLD *BSCAFFS;     // Scaffold record array
static GDB_CONTIG   *BCONTIG;     // Contig record array

static Hash_Table   *AHASH;       // Scaffold names hash table
static Hash_Table   *BHASH;       // Scaffold names hash table

static Contig_Range *ACHORD;      // [0..NACONTIG) Portion of each contig to plot (or not)
static Contig_Range *BCHORD;      // [0..NBCONTIG) Portion of each contig to plot (or not)

  //  Threading communication packet (read_1aln_block & aln_clip)

typedef struct
  { int64       beg;
    int64       end;
    OneFile    *in;
    Segment    *segs;
    int64       nseg;
    GDB_CONTIG *bctg;
  } Packet;


  // Helvetica font character width
static double HELVETICA[128] = 
  { 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
    0.278, 0.278, 0.355, 0.556, 0.556, 0.889, 0.667, 0.222, 0.333, 0.333, 0.389, 0.584, 0.278, 0.333, 0.278, 0.278, 
    0.556, 0.556, 0.556, 0.556, 0.556, 0.556, 0.556, 0.556, 0.556, 0.556, 0.278, 0.278, 0.584, 0.584, 0.584, 0.556, 
    1.015, 0.667, 0.667, 0.722, 0.722, 0.667, 0.611, 0.778, 0.722, 0.278, 0.500, 0.667, 0.556, 0.833, 0.722, 0.778, 
    0.667, 0.778, 0.722, 0.667, 0.611, 0.722, 0.667, 0.944, 0.667, 0.667, 0.611, 0.278, 0.278, 0.278, 0.469, 0.556, 
    0.222, 0.556, 0.556, 0.500, 0.556, 0.556, 0.278, 0.556, 0.556, 0.222, 0.222, 0.500, 0.222, 0.833, 0.556, 0.556, 
    0.556, 0.556, 0.333, 0.500, 0.278, 0.556, 0.500, 0.722, 0.500, 0.500, 0.500, 0.334, 0.260, 0.334, 0.584, 0.000
  };

/*******************************************************************************************
 *
 *  Utilities
 *
 *******************************************************************************************/

static char *PDFTOOLS[4] = {"pstopdf", "epstopdf", "ps2pdf", "eps2pdf"};

int run_system_cmd(char *cmd, int retry)
{ int exit_code = system(cmd);
  --retry;
  if ((exit_code != -1 && !WEXITSTATUS(exit_code)) || !retry)
    return (exit_code);
  return (run_system_cmd(cmd, retry));
}

int check_executable(char *exe)
{ char cmd[4096];
  int  exit_code;

  sprintf(cmd, "which %s 1>/dev/null 2>/dev/null", exe);
  exit_code = run_system_cmd(cmd, 1);
  if (exit_code == -1 || WEXITSTATUS(exit_code))
    return (0);
  return (1);
}

char *findPDFtool()
{ uint64 i;
  for (i = 0; i < sizeof(PDFTOOLS) / sizeof(char *); i++)
    if (check_executable(PDFTOOLS[i]))
      return (PDFTOOLS[i]);
  return (NULL);
}

static inline void expandBuffer(char **names, uint32 *n, uint32 *m, uint32 l)
{ if (*n + l > *m)
    { *m <<= 1;
      *names = (char *)Realloc(*names,sizeof(char)*(*m),"Reallocating name array");
      if (*names == NULL)
        exit (1);
    }
}

//  Read next line into a buffer and return a pointer to the buffer
//    the length of the line.  NB: replaces '\n' with '\0'.

static char *read_line(void *input, int gzipd, int nline, char *spath)
{ static char *buffer;
  static int   bmax = 0;
  int len;

  if (bmax == 0)
    { bmax = 500;
      buffer = (char *) malloc(bmax);
      if (buffer == NULL)
        { fprintf(stderr,"%s: Out of memory reading %s\n",Prog_Name,spath);
          exit (1);
        }
    }

  if (gzipd)
    { if (gzgets(input,buffer,bmax) == NULL)
        { if (gzeof(input))
            { free(buffer);
              bmax = 0;
              return (NULL);
            }
          fprintf(stderr,"%s: Could not read next line, %d, of file %s\n",Prog_Name,nline,spath);
          exit (1);
        }
    }
  else
    { if (fgets(buffer,bmax,input) == NULL)
        { if (feof(input))
            { free(buffer);
              bmax = 0;
              return (NULL);
            }
          fprintf(stderr,"%s: Could not read next line, %d, of file %s\n",Prog_Name,nline,spath);
          exit (1);
        }
    }

  len = strlen(buffer);
  while (buffer[len-1] != '\n')
    { bmax = ((int) (1.4*bmax)) + 100;
      buffer = (char *) realloc(buffer,bmax);
      if (buffer == NULL)
        { fprintf(stderr,"%s: Out of memory reading %s\n",Prog_Name,spath);
          exit (1);
        }
      if (gzipd)
        { if (gzgets(input,buffer+len,bmax-len) == NULL)
            { if (gzeof(input))
                fprintf(stderr,"%s: Last line %d of file %s does not end with new-line\n",
                               Prog_Name,nline,spath);
              else
                fprintf(stderr,"%s: Could not read next line, %d, of file %s\n",
                               Prog_Name,nline,spath);
              exit (1);
            }
        }
      else
        { if (fgets(buffer+len,bmax-len,input) == NULL)
            { if (feof(input))
                fprintf(stderr,"%s: Last line %d of file %s does not end with new-line\n",
                               Prog_Name,nline,spath);
              else
                fprintf(stderr,"%s: Could not read next line, %d, of file %s\n",
                               Prog_Name,nline,spath);
              exit (1);
            }
        }
      len += strlen(buffer+len);
    }
  buffer[--len] = '\0';

  return (buffer);
}


/*******************************************************************************************
 *
 *  Read .1aln file creating array of segments and two GDB skeletons, AGDB & BGDB,
 *    and hash tables, AHASH & BHASH, of the scaffold names for each.
 *
 *******************************************************************************************/

void *read_1aln_block(void *args)
{ Packet     *parm  = (Packet *) args;
  int64       beg  = parm->beg;
  int64       end  = parm->end;
  GDB_CONTIG *bctg = parm->bctg;
  OneFile    *in   = parm->in;
  Segment    *segs = parm->segs;
  int64       nseg = 0;
  // int         ngroup;

  //  For each alignment do

  if (!oneGoto (parm->in, 'A', beg+1))
    { fprintf(stderr,"%s: Could not locate to object %lld in 1aln file\n",Prog_Name,beg+1);
      exit (1);
    }
  oneReadLine(parm->in);

  //  Get the group id of the first alignment
  // ngroup = -1; // oneGetGroup (parm->in, beg);
  // if (ngroup < 0) ngroup = 0;

  { int64 i;
    uint32   flags;
    uint8    flag;
    int      aread, bread;
    int      abpos, bbpos;
    int      aepos, bepos;
    // int      group;
    int      diffs;
    int      blocksum, iid;

    //  Set group to be the first
    // group = ngroup;

    for (i = beg; i < end; i++)
      { // read i-th alignment record 

        if (in->lineType != 'A')
          { fprintf(stderr,"%s: Failed to be at start of alignment\n",Prog_Name);
            exit (1);
          }

        flags = 0;
        aread = oneInt(in,0);
        abpos = oneInt(in,1);
        aepos = oneInt(in,2);
        bread = oneInt(in,3);
        bbpos = oneInt(in,4);
        bepos = oneInt(in,5);

        diffs = 0;
        while (oneReadLine(in))
          if (in->lineType == 'R')
            flags |= COMP_FLAG;
          else if (in->lineType == 'D')
            diffs = oneInt(in,0);
          else if (in->lineType == 'A')
            break; // stop at next A-line
          else
            continue; // skip trace lines
      
        if (aepos - abpos < MINALEN || bepos - bbpos < MINALEN)
          continue; // filter by segment size

        blocksum = (aepos-abpos) + (bepos-bbpos);
        iid      = (blocksum - diffs) / 2;
        if (2.*iid / blocksum < MINAIDNT)
          continue; // filter by segment identity

        flag = 0;
        if (COMP(flags))
          { bbpos = bctg[bread].clen - bbpos;
            bepos = bctg[bread].clen - bepos;
          }

        // add to output
        segs->flag  = flag;
        segs->aread = aread;
        segs->abpos = abpos;
        segs->aepos = aepos;
        segs->bread = bread;
        segs->bbpos = bbpos;
        segs->bepos = bepos;
        // segs->group = group;
        segs += 1;
        nseg += 1;

        //  change group after this alignment
        // group = ngroup;
      }
  }

  parm->nseg = nseg;

  return (NULL);
}

void read_1aln(char *oneAlnFile)
{ OneFile   *input;
  int64      novl;

  //  Initiate .1aln file reading and read header information

  { char      *cpath, *tmp;
    FILE      *test;
    char      *src1_name, *src2_name;
    char      *spath, *tpath;
    char      *sptr, *eptr;
    char      *head;
    int        s, type, TSPACE;
    
    input = open_Aln_Read(oneAlnFile,NTHREADS,&novl,&TSPACE,&src1_name,&src2_name,&cpath);
    if (input == NULL)
      { fprintf(stderr,"%s: Could not open .1aln file: %s\n",Prog_Name,oneAlnFile);
        exit (1);
      }

    test = fopen(src1_name,"r");
    if (test == NULL)
      { if (*src1_name != '/')
          test = fopen(Catenate(cpath,"/",src1_name,""),"r");
        if (test == NULL)
          { fprintf(stderr,"%s: Could not find GDB %s\n",Prog_Name,src1_name);
            exit (1);
          }
        tmp = Strdup(Catenate(cpath,"/",src1_name,""),"Allocating expanded name");
        free(src1_name);
        src1_name = tmp;
      }
    fclose(test);

    if (src2_name != NULL)
      { test = fopen(src2_name,"r");
        if (test == NULL)
          { if (*src2_name != '/')
              test = fopen(Catenate(cpath,"/",src2_name,""),"r");
            if (test == NULL)
              { fprintf(stderr,"%s: Could not find GDB %s\n",Prog_Name,src2_name);
                exit (1);
              }
            tmp = Strdup(Catenate(cpath,"/",src2_name,""),"Allocating expanded name");
            free(src2_name);
            src2_name = tmp;
          }
        fclose(test);
      }

    free(cpath);

    ISTWO = 0;
    type  = Get_GDB_Paths(src1_name,NULL,&spath,&tpath,0);
    if (type != IS_GDB)
      Create_GDB(AGDB,spath,type,0,NULL);
    else
      Read_GDB(AGDB,tpath);
    free(spath);
    free(tpath);

    if (src2_name != NULL)
      { type = Get_GDB_Paths(src2_name,NULL,&spath,&tpath,0);
        if (type != IS_GDB)
          Create_GDB(BGDB,spath,type,0,NULL);
        else
          Read_GDB(BGDB,tpath);
        free(spath);
        free(tpath);
        ISTWO = 1;
      }
    else
      BGDB = AGDB;
    free(src1_name);
    free(src2_name);

    NASCAFF  = AGDB->nscaff;
    NACONTIG = AGDB->ncontig;
    ASCAFFS  = AGDB->scaffolds;
    ACONTIG  = AGDB->contigs;

    NBSCAFF  = BGDB->nscaff;
    NBCONTIG = BGDB->ncontig;
    BSCAFFS  = BGDB->scaffolds;
    BCONTIG  = BGDB->contigs;

    AHASH = New_Hash_Table(NASCAFF,0);
    head  = AGDB->headers;
    for (s = 0; s < NASCAFF; s++)
      { sptr = head + ASCAFFS[s].hoff;
        for (eptr = sptr; *eptr != '\0'; eptr++)
          if (isspace(*eptr))
            break;
        *eptr = '\0';
        if (Hash_Lookup(AHASH,sptr) < 0)
          Hash_Add(AHASH,sptr);
        else
          { fprintf(stderr,"%s: Duplicate scaffold name: %s\n",Prog_Name,sptr);
            exit (1);
          }
      }

    if (ISTWO)
      { BHASH = New_Hash_Table(NBSCAFF,0);
        head  = BGDB->headers;
        for (s = 0; s < NBSCAFF; s++)
          { sptr = head + BSCAFFS[s].hoff;
            for (eptr = sptr; *eptr != '\0'; eptr++)
              if (isspace(*eptr))
                break;
            *eptr = '\0';
            if (Hash_Lookup(BHASH,sptr) < 0)
              Hash_Add(BHASH,sptr);
            else
              { fprintf(stderr,"%s: Duplicate scaffold name: %s\n",Prog_Name,sptr);
                exit (1);
              }
          }
      }
    else
      BHASH = AHASH;
  }

#ifdef DEBUG_MAKE_HASH
  { int i;

    fprintf(stderr,"%s: DB_A\n",Prog_Name);
    fprintf(stderr,"%s: List of contigs (NACTG=%d)\n",Prog_Name,NACONTIG);
    fprintf(stderr,"%s:  INDEX   SCAF     OFFSET     LENGTH\n",Prog_Name);
    for (i = 0; i < NACONTIG; i++)
      fprintf(stderr,"%s: %6d %6d %10lld %10lld\n",Prog_Name,i,
                     ACONTIG[i].scaf,ACONTIG[i].sbeg,ACONTIG[i].clen);
    fprintf(stderr,"%s: List of scaffolds (NASCAFF=%d)\n",Prog_Name,NASCAFF);
    for (i = 0; i < NASCAFF; i++)
      fprintf(stderr,"%s: %6d %s\n",Prog_Name,i,Get_Hash_String(AHASH,i));
    if (ISTWO)
      { fprintf(stderr,"%s: DB_B\n",Prog_Name);
        fprintf(stderr,"%s: List of contigs (NBCTG=%d)\n",Prog_Name,NBCONTIG);
        fprintf(stderr,"%s:  INDEX   SCAF     OFFSET     LENGTH\n",Prog_Name);
        for (i = 0; i < NBCONTIG; i++)
          fprintf(stderr,"%s: %6d %6d %10lld %10lld\n",Prog_Name,i,
                         BCONTIG[i].scaf,BCONTIG[i].sbeg,BCONTIG[i].clen);
        fprintf(stderr,"%s: List of scaffolds (NBSCAFF=%d)\n",Prog_Name,NBSCAFF);
        for (i = 0; i < NBSCAFF; i++)
          fprintf(stderr,"%s: %6d %s\n",Prog_Name,i,Get_Hash_String(BHASH,i));
      }
  }
#endif

  // Read alignment segments

  { int p;
    Packet    parm[NTHREADS];
    pthread_t threads[NTHREADS];

    //  Divide .1aln into NTHREADS parts

    segments = (Segment *) Malloc(sizeof(Segment)*novl, "Allocating segment array");
    for (p = 0; p < NTHREADS ; p++)
      { parm[p].beg = (p * novl) / NTHREADS;
        if (p > 0)
          parm[p-1].end = parm[p].beg;
        parm[p].segs = segments + parm[p].beg;
        parm[p].nseg = 0;
        parm[p].in   = input + p;
        parm[p].bctg = BGDB->contigs;
      }
    parm[NTHREADS-1].end = novl;

    // Use NTHREADS to produce alignment segments for each part

    for (p = 1; p < NTHREADS; p++)
      pthread_create(threads+p,NULL,read_1aln_block,parm+p);
    read_1aln_block(parm);
    for (p = 1; p < NTHREADS; p++)
      pthread_join(threads[p],NULL);
  
    // Collect results from different part

    nSegment = parm[0].nseg;
    for (p = 1; p < NTHREADS; p++)
      { if (nSegment < parm[p].beg)
          memmove(segments+nSegment,parm[p].segs,parm[p].nseg);
        nSegment += parm[p].nseg;
      }
    if (nSegment < novl)
      segments = Realloc(segments,sizeof(Segment)*nSegment,"Compacting segment arrary");

#ifdef DEBUG_READ_1ALN
    fprintf(stderr, "%s: %9lld segments loaded\n",Prog_Name,novl);
    fprintf(stderr, "%s: %9lld after filtering\n",Prog_Name,nSegment);
#endif
  }

  oneFileClose(input);
}


/*******************************************************************************************
 *
 *  Read .paf file creating array of segments and two GDB skeletons, AGDB & BGDB
 *    and hash tables, AHASH & BHASH, of the scaffold names for each.
 *
 *******************************************************************************************/

void read_paf(char *pafAlnFile, int gzipd)
{ void *input;
  int   nline;
  int i, naseq, maseq, nbseq, mbseq;
  int     index;
  uint64  nsegs, msegs;
  Segment *segs;
  char *fptrs[11], *fptr, *eptr;
  char *head;
  int alen, blen, *alens, *blens;;
  int aread, bread;
  int abpos, bbpos;
  int aepos, bepos;
  // int group;
  int blocksum, iid;
  uint8 flag;

  if (gzipd)
    input = gzopen(pafAlnFile,"r");
  else
    input = fopen(pafAlnFile,"r");

  AHASH = New_Hash_Table(1024,1);
  BHASH = New_Hash_Table(1024,1);
  
  naseq = 0;
  maseq = 1024;
  nbseq = 0;
  mbseq = 1024;
  nsegs = 0;
  msegs = 4096;

  alens = (int *) Malloc(sizeof(int) * maseq, "Allocating length map");
  blens = (int *) Malloc(sizeof(int) * mbseq, "Allocating length map");
  segments = (Segment *) Malloc(sizeof(Segment) * msegs, "Allocating segment array");
  segs = segments;

  nline = 1;
  while ((eptr = read_line(input,gzipd,nline++,pafAlnFile)) != NULL)
    {
      fptr = eptr;                // parse PAF line
      for (i = 0; i < 11; i++)
        { while (*eptr != '\t' && *eptr != '\n' && *eptr != '\0')
            eptr++;
          if (eptr > fptr)
            fptrs[i] = fptr;

          if (*eptr == '\t')
            { *eptr = '\0';
              fptr = ++eptr;
            }
          else
            { *eptr = '\0';
              *fptr = '\0';
              break; 
            }
        }
      
      if (i != 11)
        continue;

      index = Hash_Lookup(AHASH, fptrs[0]);
      alen  = strtol(fptrs[1], &eptr, 10);
      if (index < 0)
        { alens[naseq++] = alen;
          if (naseq == maseq)
            { maseq <<= 1;
              alens = (int *) Realloc(alens, sizeof(int) * maseq, "Reallocating length map");
            }
          index = Hash_Add(AHASH, fptrs[0]);
        }
      aread  = index;
      abpos  = strtol(fptrs[2], &eptr, 10);
      aepos  = strtol(fptrs[3], &eptr, 10);
      
      index = Hash_Lookup(BHASH, fptrs[5]);
      blen   = strtol(fptrs[6], &eptr, 10);
      if (index < 0)
        { blens[nbseq++] = blen;
          if (nbseq == mbseq)
            { mbseq <<= 1;
              blens = (int *) Realloc(blens, sizeof(int) * mbseq, "Reallocating length map");
            }
          index = Hash_Add(BHASH, fptrs[5]);
        }
      bread  = index;
      bbpos  = strtol(fptrs[7], &eptr, 10);
      bepos  = strtol(fptrs[8], &eptr, 10);
        
      if (aepos - abpos < MINALEN || bepos - bbpos < MINALEN)     // filter by segment size
        continue;

      blocksum = (aepos-abpos) + (bepos-bbpos);
      iid      = strtol(fptrs[9], &eptr, 10);
      if (2.*iid / blocksum < MINAIDNT)                    // filter by segment identity
        continue;

      /***
      group = 0;                     // parse 'cn:i' tag to find segment group if exists
      while (fptr)
        { eptr = fptr;
          while (*eptr != '\t' && *eptr != '\n' && *eptr != '\0')
            eptr++;
          if (eptr - fptr > 5 && !strncmp(fptr, "cn:i:", 5))
            { group = strtol(fptr + 5, 0, 10);
              break;
            }
          if (*eptr != '\t')
            break;
          fptr = ++eptr;
        }
      **/

      flag = 0;                      // Add to segment array
      if (*fptrs[4] == '-')
        { int tmp = bbpos;
          bbpos = bepos;
          bepos = tmp;
        }

      segs->flag  = flag;
      segs->aread = aread;
      segs->bread = bread;
      segs->abpos = abpos;
      segs->aepos = aepos;
      segs->bbpos = bbpos;
      segs->bepos = bepos;
      // segs->group = group;
    
      nsegs += 1;
      if (nsegs == msegs)
        { msegs <<= 1;
          segments = (Segment *) Realloc(segments, sizeof(Segment) * msegs,
                                         "Reallocating segment array");
        }
      segs = segments + nsegs;
    }
 
  if (gzipd)
    gzclose(input);
  else
    fclose(input);

  ISTWO = 1;
  nSegment = nsegs;

  //  Contigs = Scaffolds so the underlying GDB is ...

  AGDB->nscaff = AGDB->ncontig = NASCAFF = NACONTIG = naseq;
  BGDB->nscaff = BGDB->ncontig = NBSCAFF = NBCONTIG = nbseq;

  AGDB->scaffolds = ASCAFFS = Malloc(sizeof(GDB_SCAFFOLD)*(NASCAFF+NBSCAFF),"Allocating scaffolds");
  AGDB->contigs   = ACONTIG = Malloc(sizeof(GDB_CONTIG)*(NACONTIG+NBCONTIG),"Allocating contigs");
  AGDB->headers   = Get_Hash_String(AHASH,0);

  BGDB->scaffolds = BSCAFFS = ASCAFFS + NASCAFF;
  BGDB->contigs   = BCONTIG = ACONTIG + NACONTIG;
  BGDB->headers   = Get_Hash_String(BHASH,0);

  head = AGDB->headers;
  for (i = 0; i < NASCAFF; i++)
    { ASCAFFS[i].slen = alens[i];
      ASCAFFS[i].fctg = i;
      ASCAFFS[i].ectg = i+1;
      ASCAFFS[i].hoff = Get_Hash_String(AHASH,i) - head;
    }
  for (i = 0; i < NACONTIG; i++)
    { ACONTIG[i].clen = alens[i];
      ACONTIG[i].sbeg = 0;
      ACONTIG[i].boff = 0;
      ACONTIG[i].scaf = i;
    }

  head = BGDB->headers;
  for (i = 0; i < NBSCAFF; i++)
    { BSCAFFS[i].slen = blens[i];
      BSCAFFS[i].fctg = i;
      BSCAFFS[i].ectg = i+1;
      BSCAFFS[i].hoff = Get_Hash_String(BHASH,i) - head;
    }
  for (i = 0; i < NBCONTIG; i++)
    { BCONTIG[i].clen = blens[i];
      BCONTIG[i].sbeg = 0;
      BCONTIG[i].boff = 0;
      BCONTIG[i].scaf = i;
    }

#ifdef DEBUG_MAKE_HASH
  { int i;

    fprintf(stderr,"%s: DB_A\n",Prog_Name);
    fprintf(stderr,"%s: List of contigs (NACTG=%d)\n",Prog_Name,NACONTIG);
    fprintf(stderr,"%s:  INDEX   SCAF     OFFSET     LENGTH\n",Prog_Name);
    for (i = 0; i < NACONTIG; i++)
      fprintf(stderr,"%s: %6d %6d %10lld %10lld\n",Prog_Name,i,
                     ACONTIG[i].scaf,ACONTIG[i].sbeg,ACONTIG[i].clen);
    fprintf(stderr,"%s: List of scaffolds (NASCAFF=%d)\n",Prog_Name,NASCAFF);
    for (i = 0; i < NASCAFF; i++)
      fprintf(stderr,"%s: %6d %s\n",Prog_Name,i,Get_Hash_String(AHASH,i));
    if (ISTWO)
      { fprintf(stderr,"%s: DB_B\n",Prog_Name);
        fprintf(stderr,"%s: List of contigs (NBCTG=%d)\n",Prog_Name,NBCONTIG);
        fprintf(stderr,"%s:  INDEX   SCAF     OFFSET     LENGTH\n",Prog_Name);
        for (i = 0; i < NBCONTIG; i++)
          fprintf(stderr,"%s: %6d %6d %10lld %10lld\n",Prog_Name,i,
                         BCONTIG[i].scaf,BCONTIG[i].sbeg,BCONTIG[i].clen);
        fprintf(stderr,"%s: List of scaffolds (NBSCAFF=%d)\n",Prog_Name,NBSCAFF);
        for (i = 0; i < NBSCAFF; i++)
          fprintf(stderr,"%s: %6d %s\n",Prog_Name,i,Get_Hash_String(BHASH,i));
      }
  }
#endif

  free(alens);
  free(blens);
}


/*******************************************************************************************
 *
 *  Chenxi's plotting code
 *
 *******************************************************************************************/

static int USORT(const void *l, const void *r)
{ uint64 x = *((uint64 *) l);
  uint64 y = *((uint64 *) r);

  if (x > y)
    return (1);
  return (-1);
}

/*  Currently not needed

static void axisReverse(uint64 *sarray, int64 *caxis, int64 soff, int beg, int end,
        int *cbeg, int *cend)
{ int i, c, clen;
  int64 coff;
  coff = caxis[(uint32) sarray[beg]];
  for (i = beg; i < end; i++)
    { c = (uint32) sarray[i];
      soff -= caxis[c] - coff; // gap
      clen  = cend[c] - cbeg[c];
      soff -= clen; // ctg
      coff  = caxis[c] + clen;
      caxis[c] = soff;
    }
}

*/

// Now string length of a name is capped by MAX_LAB_LEN
// If the length N exceeds the limit then 
// substring name[MAX_LAB_LEN-2,N-2] is replaced by a "*"
static void addSeqName(char *names, uint32 n, uint32 m, int c0, int c1, 
                       uint32 s, Hash_Table *hash, GDB *gdb, Contig_Range *chord, int rank,
                       char **_names, uint32 *_n, uint32 *_m)
{ uint32 l, n0;
  int64 p;
  char *name;

  n0 = n;

  if (PRINTSID)
    { l = Number_Digits(s+1);
      expandBuffer(&names, &n, &m, l+1);
      sprintf(names+n,"%u",s+1);
      n += l;
    }
  else
    { name = Get_Hash_String(hash,s);
      l = strlen(name);
      expandBuffer(&names, &n, &m, l+1);
      sprintf(names+n,"%s",name);
      n += l;
    }
  if (chord[c0].beg > 0 ||
          gdb->scaffolds[s].fctg != c0 ||      // start from first
          gdb->scaffolds[s].ectg != c1+1 ||
          chord[c1].end != gdb->contigs[c1].clen) // end at last
    { // partial scaffold
      p = gdb->contigs[c0].sbeg + chord[c0].beg;
      p += 1; // make it 1-based
      l = Number_Digits(p);
      l += 1;
      expandBuffer(&names, &n, &m, l+1);
      sprintf(names+n,"_%lld",p);
      n += l;
      p = gdb->contigs[c1].sbeg + chord[c1].end;
      l = Number_Digits(p);
      l += 1;
      expandBuffer(&names, &n, &m, l+1);
      sprintf(names+n,"-%lld",p);
      n += l;
    }
  if (rank < 0)
    { // add a prime
      expandBuffer(&names, &n, &m, 2);
      names[n++] = '\'';
      names[n]   = '\0';
    }
  if (n > n0 + MAX_LAB_LEN) {
    // string length exceed limits
    l = n - n0;
    if (rank < 0) {
      names[n0+MAX_LAB_LEN-3] = '*';
      names[n0+MAX_LAB_LEN-2] = names[n-2]; // the last character
      names[n0+MAX_LAB_LEN-1] = names[n-1]; // the prime
      names[n0+MAX_LAB_LEN]   = '\0';
    } else {
      names[n0+MAX_LAB_LEN-2] = '*';
      names[n0+MAX_LAB_LEN-1] = names[n-1]; // the last character
      names[n0+MAX_LAB_LEN]   = '\0';
    }
    n = n0 + MAX_LAB_LEN;
  }
  n++; // the null terminator

// assign_vars:
  *_names = names;
  *_n = n;
  *_m = m;
}

static double seqNameWidth(char *names, int n)
{ int i, c;
  double l, w;
  w = 0;
  for (i = 0; i < n; i++)
  { l = 0;
    while (*names != '\0')
      { c = (uint8) (*names);
        if (c < 33 || c > 126)
          { fprintf(stderr,"%s: unsupported character '%c'\n",Prog_Name,*names);
            exit (1);
          }
        l += HELVETICA[c];
        names++;
      }
    if (l > w) w = l;
    names++;
  }
  return (w);
}

static double seqNameRenderWidth(char *names, int n, int64 *soff, double unit, double space)
{ int i, c;
  double l, w, s;
  w = 0;
  for (i = 0; i < n; i++)
  { s = (i == 0 ? soff[0] : soff[i] - soff[i-1]);
    if (s * unit < space)
      { names += strlen(names) + 1;
        continue;
      }
    l = 0;
    while (*names != '\0')
      { c = (uint8) (*names);
        if (c < 33 || c > 126)
          { fprintf(stderr,"%s: unsupported character '%c'\n",Prog_Name,*names);
            exit (1);
          }
        l += HELVETICA[c];
        names++;
      }
    if (l > w) w = l;
    names++;
  }
  return (w);
}

static double chooseFontSizeByHeight(int64 *soff, int n, double unit, double minf, double maxf)
{ int i;
  double s, f;
  f = maxf;
  for (i = 0; i < n; i++)
  { s  = (i == 0 ? soff[0] : soff[i] - soff[i-1]);
    s *= unit;
    if (s >= minf && s < f) f = s;
  }
  return (f);
}

int axisConfig(Hash_Table *hash, GDB *gdb, Contig_Range *chord,
               int64 **_caxis, int64 **_saxis, char **_names, int64 *_tseq)
{ int i, j, s, r1, r2, mseq, nseq;
  uint32 c0, c1, c2, nstr, mstr;
  uint64 *sarray;
  int64 *caxis, *saxis, tseq;
  char *names;

  mseq = 0;
  for (i = 0; i < gdb->ncontig; i++)
    if (chord[i].order >= 0)
      mseq++;

  sarray = (uint64 *) Malloc(sizeof(uint64)*mseq,"Allocating seq array");
  saxis  = (int64 *) Malloc(sizeof(int64)*mseq,"Allocating offset array");
  caxis  = (int64 *) Malloc(sizeof(int64)*gdb->ncontig,"Allocating offset array");

  if (LABELS)
    { mstr = 2048;
      names = (char *) Malloc(sizeof(char)*mstr,"Allocating name array");
    }
  else
    names = NULL;

  mseq = 0;
  for (i = 0; i < gdb->ncontig; i++)
    if (chord[i].order)
      sarray[mseq++] = (uint64) abs(chord[i].order) << 32 | i;

  qsort(sarray, mseq, sizeof(uint64), USORT);

  c1 = (uint32) sarray[0];
  r1 = chord[c1].order;
  tseq = 0;
  nseq = 0;
  nstr = 0;
  for (j = 0, i = 1; i < mseq; i++)
    { caxis[c1] = tseq - chord[c1].beg;
      tseq += chord[c1].end - chord[c1].beg;
      c2 = (uint32) sarray[i];
      r2 = chord[c2].order;
      if (chord[c1].end < gdb->contigs[c1].clen ||      // end-clipping
              c1 + 1 < c2 ||          // contigs skipped
              gdb->contigs[c1].scaf != gdb->contigs[c2].scaf || // scaffolds switched
              r1 != r2 ||             // ranks changed
              chord[c2].beg > 0)           // start-clipping
        { // different rank group or gap
          // add new seq [j,i)
          c0 = (uint32) sarray[j];
          s = gdb->contigs[c0].scaf; // get scaffold id
          if (LABELS)
            addSeqName(names,nstr,mstr,c0,c1,s,hash,gdb,chord,r1,&names,&nstr,&mstr);
          saxis[nseq++] = tseq;
          // reverse ctg [j,i)
          // if (r1 < 0) axisReverse(sarray, caxis, tseq, j, i, cbeg, cend);
          j = i;
        }
      else
        tseq += gdb->contigs[c2].sbeg - gdb->contigs[c1].sbeg - gdb->contigs[c1].clen;
      c1 = c2;
      r1 = r2;
    }
  // add new seq [j,i)
  caxis[c1] = tseq - chord[c1].beg;
  tseq += chord[c1].end - chord[c1].beg;
  c0 = (uint32) sarray[j];
  s = gdb->contigs[c0].scaf; // get scaffold id
  if (LABELS)
    addSeqName(names,nstr,mstr,c0,c1,s,hash,gdb,chord,r1,&names,&nstr,&mstr);
  saxis[nseq++] = tseq;
  // if (r1 < 0) axisReverse(sarray, caxis, tseq, j, i, cbeg, cend);

  if (_caxis) *_caxis = caxis;
  if (_saxis) *_saxis = saxis;
  if (_names) *_names = names;
  if (_tseq)  *_tseq = tseq;

#ifdef DEBUG_AXIS_CONF
  fprintf(stderr,"%s: sequence number: %d\n",Prog_Name,nseq);
  fprintf(stderr,"%s: sequence length: %lld\n",Prog_Name,tseq);
  for (i = 0; i < nctg; i++)
    if (cseq[i])
      fprintf(stderr,"%s: sequence [ctg] offset: %-8d %-12lld\n",Prog_Name,i,caxis[i]);
  for (i = 0; i < nseq; i++)
    { fprintf(stderr,"%s: sequence [scf] offset: %-8d %-12lld %s\n",Prog_Name,i,saxis[i],LABELS? names : "");
      if (LABELS)
        names += strlen(names) + 1;
    }
#endif

  free(sarray);

  return (nseq);
}

#define INT_SIGN(x) (((x)>0)-((x)<0))

static void alnConfig()
{ int64 i;
  int abpos, aepos, bbpos, bepos, a, b;
  Segment *segment;

  for (i = 0; i < nSegment; i++)
    { segment = &segments[i];
      if (IS_DEL(segment->flag))
        continue;

      abpos = segment->abpos;
      aepos = segment->aepos;
      bbpos = segment->bbpos;
      bepos = segment->bepos;

      // if (ACHORD[aread].order < 0)         //  orientation not currently available
        // { l = AEND[aread] - ABEG[aread];
          // abpos = l - abpos;
          // aepos = l - aepos;
        // }
      // if (BSEQ[bread].order < 0)
        // { l = BEND[bread] - BBEG[bread];
          // bbpos = l - bbpos;
          // bepos = l - bepos;
        // }

      segment->abpos = abpos;
      segment->aepos = aepos;
      segment->bbpos = bbpos;
      segment->bepos = bepos;

      a = abpos - aepos;
      b = bbpos - bepos;

      // if (HIGHLIGHT && !segment->group)
      //   SET_GRAY(segment->flag);
      // else 
      if (INT_SIGN(a) == INT_SIGN(b))
        SET_RED(segment->flag);
      else
        SET_BLUE(segment->flag);
    }
}

int myers_clip(int *nx1, int *ny1, int *nx2, int *ny2,
               double xmin, double xmax, double ymin, double ymax)
{ double t, x1, x2, y1, y2;
  int    flipx, flipy, inter;

  inter = 0;
  flipx = (*nx1 > *nx2);
  if (flipx)
    { x1 = *nx2; x2 = *nx1;
      y1 = *ny2; y2 = *ny1; 
    }
  else
    { x1 = *nx1; x2 = *nx2;
      y1 = *ny1; y2 = *ny2; 
    }
  if (x2 <= xmin || x1 >= xmax)
    return (-1);
  flipy = (y1 > y2);
  if (flipy)
    { t = x1; x1 = x2; x2 = t;
      t = y1; y1 = y2; y2 = t;
    }
  if (y2 <= ymin || y1 >= ymax)
    return (-1);
  if (y2 > ymax)
    { x2 = x1 + (x2 - x1) * (ymax - y1) / (y2 - y1);
      y2 = ymax;
      inter = 1;
    }
  if (y1 < ymin)
    { x1 = x1 + (x2 - x1) * (ymin - y1) / (y2 - y1);
      y1 = ymin;
      inter = 1;
    }
  if (flipy)
    { t = x1; x1 = x2; x2 = t;
      t = y1; y1 = y2; y2 = t;
    }
  if (x2 > xmax)
    { if (x1 >= xmax)
        return (-1);
      y2 = y1 + (y2 - y1) * (xmax - x1) / (x2 - x1);
      x2 = xmax;
      inter = 1;
    }
  if (x1 < xmin)
    { if (x2 <= xmin)
        return (-1);
      y1 = y1 + (y2 - y1) * (xmin - x1) / (x2 - x1);
      x1 = xmin;
      inter = 1;
    }
  if (inter)
    { if (flipx)
        { *nx1 = (int) (x2 + .499); *nx2 = (int) (x1 + .499);
          *ny1 = (int) (y2 + .499); *ny2 = (int) (y1 + .499);
        }
      else
        { *nx1 = (int) (x1 + .499); *nx2 = (int) (x2 + .499);
          *ny1 = (int) (y1 + .499); *ny2 = (int) (y2 + .499);
        }
    }
  return (0);
}

void *aln_clip(void *args)
{ Packet *parm = (Packet *) args;
  int64 beg    = parm->beg;
  int64 end    = parm->end;

  int64    i;
  Segment *seg;
  int      x, aread, bread, nseg;

  nseg = 0;
  for (i = beg; i < end; i++)
    { seg = segments + i;
      // if (IS_DEL(segment->flag))   //  Not possible
        // continue;
      aread = seg->aread;
      bread = seg->bread;
      if (ACHORD[aread].order == 0 || BCHORD[bread].order == 0)
        { SET_DEL(seg->flag);
          continue;
        }

      x = myers_clip(&(seg->abpos), &(seg->bbpos), &(seg->aepos), &(seg->bepos),
                     ACHORD[aread].beg, ACHORD[aread].end, BCHORD[bread].beg, BCHORD[bread].end);

      if (x < 0)
        SET_DEL(seg->flag);
      else
        nseg += 1;
    }

  parm->nseg = nseg;
  return (NULL);
}

static int DSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);
  
  return (y-x);
}

void aln_filter()
{ int      i, digits;
  int64    nseg;
  int     *sarray, alen;
  Segment *s;

  // recalculate alignment coordinates

  { int       p;
    Packet   *parm;
    pthread_t threads[NTHREADS];

    parm = Malloc(sizeof(Packet)*NTHREADS,"Allocating thread records");
    for (p = 0; p < NTHREADS ; p++)
      { parm[p].beg = (p * nSegment) / NTHREADS;
        if (p > 0)
          parm[p-1].end = parm[p].beg;
      }
    parm[NTHREADS-1].end = nSegment;

    for (p = 1; p < NTHREADS; p++)
      pthread_create(threads+p,NULL,aln_clip,parm+p);
    aln_clip(parm);

    nseg = parm->nseg;
    for (p = 1; p < NTHREADS; p++)
      { pthread_join(threads[p],NULL);
        nseg += parm[p].nseg;
      }

    free(parm);
  }

  if (MAXALIGN == 0 || nseg <= MAXALIGN) return;

  sarray = (int *) Malloc(sizeof(int)*nseg,"Allocating seq array");
  
  nseg = 0;
  for (i = 0, s = segments; i < nSegment; i++, s++)
    if (!IS_DEL(s->flag))
      sarray[nseg++] = s->aepos-s->abpos;

  qsort(sarray, nseg, sizeof(int), DSORT);

  alen   = sarray[MAXALIGN-1];   // replace low digits with 0 up to loosing 10%
  digits = 1;
  while (1)
    { if ((alen/digits)*digits < .9*alen)  
        break;
      digits *= 10;
    }
  digits /= 10;
  alen    = (alen/digits)*digits;

  nseg = 0;
  for (i = 0, s = segments; i < nSegment; i++, s++)
    { if (IS_DEL(s->flag))
        continue;
      if (s->aepos-s->abpos < alen)
        SET_DEL(s->flag);
      else
        nseg++;
    }

  free(sarray);

  if (VERBOSE)
    { fprintf(stderr, "  Using length filter threshold %d\n",alen);
      fprintf(stderr, "  Selected %lld alignments to plot\n",nseg); 
    }
}

#define eps_header(fp,x,y,linewidth) { \
    fprintf(fp,"%%!PS-Adobe-3.0 EPSF-3.0\n"); \
    fprintf(fp,"%%%%BoundingBox:"); \
    fprintf(fp," 1 1 %g %g\n\n",(float)(x),(float)(y)); \
    fprintf(fp,"/C { dup 255 and 255 div exch dup -8 bitshift 255 and 255 div 3"); \
    fprintf(fp," 1 roll -16 bitshift 255 and 255 div 3 1 roll setrgbcolor } bind def\n"); \
    fprintf(fp,"/L { 4 2 roll moveto lineto } bind def\n"); \
    fprintf(fp,"/LX { dup 4 -1 roll exch moveto lineto } bind def\n"); \
    fprintf(fp,"/LY { dup 4 -1 roll moveto exch lineto } bind def\n"); \
    fprintf(fp,"/LS { 3 1 roll moveto show } bind def\n"); \
    fprintf(fp,"/MS { dup stringwidth pop 2 div 4 -1 roll exch sub 3 -1"); \
    fprintf(fp," roll moveto show } bind def\n"); \
    fprintf(fp,"/RS { dup stringwidth pop 4 -1 roll exch sub 3 -1 roll moveto show } bind def\n"); \
    fprintf(fp,"/B { 4 copy 3 1 roll exch 6 2 roll 8 -2 roll moveto lineto"); \
    fprintf(fp," lineto lineto closepath } bind def\n"); \
    fprintf(fp,"%g setlinewidth\n\n",linewidth);\
}

#define eps_font(fp,f,s) do { \
    fprintf(fp,"/FS %d def\n",s); \
    fprintf(fp,"/FS4 FS 4 div def\n"); \
    fprintf(fp,"/%s findfont FS scalefont setfont\n\n",f); \
} while (0)

#define eps_Ralign(fp) do { \
    fprintf(fp,"/RightAlignedText {\n"); \
    fprintf(fp,"  /str exch def\n"); \
    fprintf(fp,"  /y exch def\n"); \
    fprintf(fp,"  /x exch def\n"); \
    fprintf(fp,"  str stringwidth pop\n"); \
    fprintf(fp,"  x exch sub\n"); \
    fprintf(fp,"  y moveto\n"); \
    fprintf(fp,"  str show\n"); \
    fprintf(fp,"} def\n\n"); \
} while (0)

#define eps_RAlignedStr(fp,x,y,s) do { \
    fprintf(fp,"%g %g (%s) RightAlignedText\n",(float)(x),(float)(y),s); \
} while (0)

#define eps_RotatedStr(fp,x,y,r,s) do { \
    fprintf(fp,"/str (%s) def\n",s); \
    fprintf(fp,"gsave\n"); \
    fprintf(fp,"%g %g moveto\n",(float)(x),(float)(y)); \
    fprintf(fp,"%g rotate\n",(float)(r)); \
    fprintf(fp,"str show\n"); \
    fprintf(fp,"grestore\n"); \
} while (0)

#define eps_bottom(fp) fprintf(fp,"stroke showpage\n")
#define eps_color(fp,col) fprintf(fp,"stroke %d C\n",col)
#define eps_gray(fp,gray) fprintf(fp, "%g setgray\n",(float)gray)
#define eps_linewidth(fp, lw) fprintf(fp, "%g setlinewidth\n", (float)(lw))
#define eps_line(fp,x1,y1,x2,y2) \
               fprintf(fp,"%g %g %g %g L\n",(float)(x1),(float)(y1),(float)(x2),(float)(y2))
#define eps_linex(fp,x1,x2,y) fprintf(fp,"%g %g %g LX\n",(float)(x1),(float)(x2),(float)(y))
#define eps_liney(fp,y1,y2,x) fprintf(fp,"%g %g %g LY\n",(float)(y1),(float)(y2),(float)(x))
#define eps_Mstr(fp,x,y,s) fprintf(fp,"%g %g (%s) MS\n",(float)(x),(float)(y),s)
#define eps_Mint(fp,x,y,i) fprintf(fp,"%g %g (%d) MS\n",(float)(x),(float)(y),i)
#define eps_stroke(fp) fprintf(fp,"stroke\n")

#define G_COLOR 0x808080
#define N_COLOR 0xFF0000
#define C_COLOR 0x0080FF

static int SEG_COLOR[3] = { G_COLOR, N_COLOR, C_COLOR };

  // generate eps file

void make_plot(FILE *fo)
{ int    width, height, fsize, maxis, xmargin, ymargin;
  int    nxseq, nyseq;
  char  *xnames, *ynames;
  int64  txseq, tyseq, *cxoff, *sxoff, *cyoff, *syoff;
  double sx, sy;
  double lsize, bsize, gsize; // line size, border size, grid size
  int    i, c;

  // find total length of x- and y-axis and order of plotting

  txseq = tyseq = 0;
  nxseq = axisConfig(BHASH,BGDB,BCHORD,&cxoff,&sxoff,&xnames,&txseq);
  nyseq = axisConfig(AHASH,AGDB,ACHORD,&cyoff,&syoff,&ynames,&tyseq);

  alnConfig();

  width  = IMGWIDTH;
  height = IMGHEIGH;
  if (!height)
    height = (int) ((((double) width) / txseq) * tyseq + .499);
  if (!width)
    width = (int)((((double) height) / tyseq) * txseq + .499);
  
  maxis = (width > height ? width : height);
  if (maxis > MAX_XY_LEN)
    { double scale = (double) MAX_XY_LEN / maxis;
      int    plen  = strlen(Prog_Name);

      width  = (int)(width  * scale + 0.499);
      height = (int)(height * scale + 0.499);
      fprintf(stderr,"%s: Image size too large [%d]x[%d]\n",Prog_Name,width,height);
      fprintf(stderr,"%*s  Shrinking the size to [%d]x[%d]\n",plen,"",width,height);
      
      if (width < MIN_XY_LEN)
        { fprintf(stderr,"%s: Image width too small [%d]\n",Prog_Name,width);
          fprintf(stderr,"%*s  Reseting image width to [%d]\n",plen,"",MAX_XY_LEN);
          fprintf(stderr,"%*s  Axes are no longer at same scale\n",plen,"");
          width = MIN_XY_LEN;  
        }
      if (height < MIN_XY_LEN)
        { fprintf(stderr,"%s: Image height too small [%d]\n",Prog_Name,height);
          fprintf(stderr,"%*s  Reseting image height to [%d]\n",plen,"",MAX_XY_LEN);
          fprintf(stderr,"%*s  Axes are no longer at same scale\n",plen,"");
          height = MIN_XY_LEN;
        }
    }
  
  maxis = (width < height ? width : height);
  if (maxis < MIN_XY_LEN)
    { double scale = (double) MIN_XY_LEN / maxis;
      int    plen  = strlen(Prog_Name);

      width  = (int)(width  * scale + 0.499);
      height = (int)(height * scale + 0.499);
      fprintf(stderr,"%s: Image size too small [%d]x[%d]\n",Prog_Name,width,height);
      fprintf(stderr,"%*s  Expanding the size to [%d]x[%d]\n",plen,"",width,height);

      if (width > MAX_XY_LEN)
        { fprintf(stderr,"%s: Image width too large [%d]\n",Prog_Name,width);
          fprintf(stderr,"%*s  Reseting image width to [%d]\n",plen,"", MAX_XY_LEN);
          fprintf(stderr,"%*s  Axes are no longer at same scale\n",plen,"");
          width = MAX_XY_LEN;
        }
      if (height > MAX_XY_LEN)
        { fprintf(stderr,"%s: Image height too large [%d]\n",Prog_Name,height);
          fprintf(stderr,"%*s  Reseting image height to [%d]\n",plen,"",MAX_XY_LEN);
          fprintf(stderr,"%*s  Axes are no longer at same scale\n",plen,"");
          height = MAX_XY_LEN;
        }
    }

  maxis = (width < height ? width : height);
  
  lsize = LINESIZE;
  if (lsize < 1e-6)
    lsize = (double) maxis / 500;
  bsize = lsize * 2;
  gsize = lsize / 2;

  sx = (double)  width / txseq;
  sy = (double) height / tyseq;

  xmargin = bsize*2;
  ymargin = bsize*2;

  fsize = FONTSIZE;
  if (!fsize)
    { if (LABELS)
        { double xlabw, ylabw; // x/y lab width
          double xfsize, yfsize; // adjusted x/y font size

          // choose a font size by height
          // the smallest that is between 1 and 2 percent of axis size
          xfsize = chooseFontSizeByHeight(sxoff, nxseq, sx, (double) maxis / 100, (double) maxis / 50);
          yfsize = chooseFontSizeByHeight(syoff, nyseq, sy, (double) maxis / 100, (double) maxis / 50);
          fsize  = (xfsize < yfsize? xfsize : yfsize);

          // adjust font size by sequence name width
          xlabw = seqNameWidth(xnames, nxseq);
          ylabw = seqNameWidth(ynames, nyseq);

          // reduce font size if label panel is too big
          if (xlabw * fsize > height * MAX_LAB_FRC)
            fsize = height * MAX_LAB_FRC / xlabw;
          if (ylabw * fsize > width * MAX_LAB_FRC)
            fsize = width * MAX_LAB_FRC / ylabw;

          // adjust x/y lab width
          // some labels will be excluded since no enough space
          xlabw = seqNameRenderWidth(xnames, nxseq, sxoff, sx, fsize);
          ylabw = seqNameRenderWidth(ynames, nyseq, syoff, sy, fsize);

          xmargin += fsize * ylabw;
          ymargin += fsize * xlabw;
        }
      else
        fsize = 10; // this has no effect
    }

  eps_header(fo, width+xmargin+bsize*3, height+ymargin+bsize*3, lsize);
  eps_font(fo, "Helvetica-Narrow", fsize);
  eps_Ralign(fo);

  if (LABELS)
    { char *names;
      double aoff;
      aoff  = xmargin > ymargin ? ymargin * 0.1 : xmargin * 0.1;
      // write x labels
      names = xnames;
      if (sxoff[0] * sx >= fsize)
        eps_RotatedStr(fo, xmargin + bsize + .5 * sxoff[0] * sx - fsize / 2, ymargin - aoff, 270, names);
      names += strlen(names) + 1;
      for (i = 1; i < nxseq; i++)
        { if ((sxoff[i] - sxoff[i-1]) * sx >= fsize)
            eps_RotatedStr(fo, xmargin + bsize + .5 * (sxoff[i-1] + sxoff[i]) * sx - fsize / 2, ymargin - aoff, 270, names);
          names += strlen(names) + 1;
        }
      // write y labels
      names = ynames;
      if (syoff[0] * sy >= fsize)
        eps_RAlignedStr(fo, xmargin - aoff, ymargin + bsize + .5 * syoff[0] * sy - fsize / 2, names);
      names += strlen(names) + 1;
      for (i = 1; i < nyseq; i++)
        { if ((syoff[i] - syoff[i-1]) * sy >= fsize)
            eps_RAlignedStr(fo, xmargin - aoff, ymargin + bsize + .5 * (syoff[i-1] + syoff[i]) * sy - fsize / 2, names);
          names += strlen(names) + 1;
        }
    }

  // write grid lines
  eps_gray(fo, .6);
  eps_linewidth(fo, gsize);
  for (i = 0; i < nyseq-1; i++)
    eps_linex(fo, xmargin, xmargin+bsize*2+width,  ymargin+bsize+syoff[i]*sy-gsize/2);
  for (i = 0; i < nxseq-1; i++)
    eps_liney(fo, ymargin, ymargin+bsize*2+height, xmargin+bsize+sxoff[i]*sx-gsize/2);
  eps_stroke(fo);
  eps_gray(fo, 0);

  // write border lines
  eps_linewidth(fo, bsize);
  eps_linex(fo, xmargin, xmargin+bsize*2+width,  ymargin+bsize/2);
  eps_linex(fo, xmargin, xmargin+bsize*2+width,  ymargin+height+bsize*3/2);
  eps_liney(fo, ymargin, ymargin+bsize*2+height, xmargin+bsize/2);
  eps_liney(fo, ymargin, ymargin+bsize*2+height, xmargin+width+bsize*3/2);
  eps_stroke(fo);

  // write segments
  { int      aread, bread, iflag;
    double   x0, y0, x1, y1, xo, yo, xoff, yoff;
    Segment *seg;

    xoff = xmargin + bsize;
    yoff = ymargin + bsize;

    eps_linewidth(fo, lsize);
    for (c = 0; c < 3; c++)
      { eps_color(fo, SEG_COLOR[c]);
        iflag = (1 << (c+1));
        for (i = 0; i < nSegment; i++)
          { seg = segments + i;
            if (seg->flag != iflag)
              continue;

            aread = seg->aread;
            bread = seg->bread;

            xo = cxoff[bread];
            yo = cyoff[aread];
            x0 = seg->bbpos;
            x1 = seg->bepos;
            y0 = seg->abpos;
            y1 = seg->aepos;

            x0 = xoff + (x0 + xo) * sx;
            x1 = xoff + (x1 + xo) * sx;
            y0 = yoff + (y0 + yo) * sy;
            y1 = yoff + (y1 + yo) * sy;

            eps_line(fo, x0, y0, x1, y1);
          }
        eps_stroke(fo);
      }
    eps_bottom(fo);
  }

  free(cxoff);
  free(cyoff);
  free(sxoff);
  free(syoff);
  free(xnames);
  free(ynames);
}

int main(int argc, char *argv[])
{ char  *xseq, *yseq, *pdf;
  FILE  *foeps;
  char  *pdftool;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("ALNplot")
    
    pdf  = NULL;
    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
          { default:
              ARG_FLAGS("vSL")
              // ARG_FLAGS("vhdSL")
              break;
            case 'f':
              ARG_POSITIVE(FONTSIZE,"Label font size")
              break;
            case 't':
              ARG_REAL(LINESIZE)
              break;
            case 'n':
              ARG_NON_NEGATIVE(MAXALIGN,"Maximium number of lines")
              break;
            case 'i':
              ARG_REAL(MINAIDNT)
              break;
            case 'l':
              ARG_NON_NEGATIVE(MINALEN,"Minimum alignment length")
              break;
            case 'p':
              if (argv[i][2] == ':' && argv[i][3] != '\0')
                pdf = argv[i]+3;
              else
                pdf = "";
              break;
            case 'H':
              ARG_POSITIVE(IMGHEIGH,"Image height")
              break;
            case 'T':
              ARG_POSITIVE(NTHREADS,"Number of threads")
              break;
            case 'W':
              ARG_POSITIVE(IMGWIDTH,"Image width")
              break;
          }
        else
          argv[j++] = argv[i];
    argc = j;
  
    VERBOSE   = flags['v'];
    // HIGHLIGHT = flags['h'];
    // TRYADIAG  = flags['d'];
    PRINTSID  = flags['S'];
    LABELS    = 1-flags['L'];
  
    if (argc < 2 || argc > 4)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[3]);
        fprintf(stderr,"\n");
        fprintf(stderr,"     <selection> = <range> [ , <range> ]*\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"     <range> =     <contig>[-<contig>]    |  <contig>_<int>-(<int>|#)\n");
        fprintf(stderr,"             | @[<scaffold>[-<scaffold>]] | @<scaffold>_<int>-(<int>|#)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"        <contig>   = (<int>|#)[.(<int>|#)]\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"        <scaffold> =  <int>|<string>|#\n");
        fprintf(stderr,"\n");
        // fprintf(stderr,"      -d: try to put alignments along the diagonal line\n");
        fprintf(stderr,"      -S: print sequence IDs as labels instead of names\n");
        fprintf(stderr,"      -L: do not print labels\n");
        fprintf(stderr,"      -T: use -T threads\n");
        fprintf(stderr,"      -p: make PDF output (requires \'[e]ps[to|2]pdf\')\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -l: minimum alignment length\n");
        fprintf(stderr,"      -i: minimum alignment identity\n");
        fprintf(stderr,"      -n: maximum number of lines to display (set '0' to force all)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -H: image height\n");
        fprintf(stderr,"      -W: image width\n");
        fprintf(stderr,"      -f: label font size\n");
        fprintf(stderr,"      -t: line thickness\n");
        // fprintf(stderr,"      -h: highlight sequences in groups\n");
        exit (1);
      }
  
    xseq = NULL;
    yseq = NULL;
    if (argc == 3)
      xseq = argv[2];
    else if (argc == 4)
      { if (strcmp(argv[2],"-") != 0)
          xseq = argv[2];
        yseq = argv[3];
      }
  
    if (pdf != NULL && !(pdftool = findPDFtool()))
      { fprintf(stderr,"%s: Cannot find [e]ps[to|2]pdf needed to produce .pdf output\n",Prog_Name);
        exit (1);
      }
  
    // if (TRYADIAG)
    //   { fprintf(stderr,"%s: diagonalisation (-d) is not supported yet\n",Prog_Name);
    //     exit (1);
    //   }
  
    if (IMGWIDTH && IMGHEIGH)
      fprintf(stderr,"%s: setting both image width and height is not recommended\n",Prog_Name);
  
    if (!IMGWIDTH && !IMGHEIGH)
      IMGHEIGH = 600;
  }
  
  { char *pwd, *root;
    int   ispaf, gzipd;
    FILE *input;
    char *name;

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".1aln");
    name  = Catenate(pwd,"/",root,".1aln");
    input = fopen(name,"r");
    if (input == NULL)
      { free(root);
        root  = Root(argv[1],".paf");
        name  = Catenate(pwd,"/",root,".paf");
        input = fopen(name,"r");
        if (input == NULL)
          { free(root);
            root  = Root(argv[1],".paf.gz");
            name  = Catenate(pwd,"/",root,".paf.gz");
            input = fopen(name,"r");
            if (input == NULL)
              { fprintf(stderr,"%s: Cannot open %s as a .1aln or .paf file\n",Prog_Name,argv[1]);
                exit (1);
              }
            gzipd = 1;
          }
        else
          gzipd = 0;
        ispaf = 1;
      }
    else
      ispaf = 0;
    fclose(input);
    name = strdup(name);

    if (pdf != NULL)
      { if (*pdf == '\0')
          OUTEPS = strdup(Catenate(pwd,"/",root,".eps"));
        else
          { free(pwd);
            free(root);
            pwd    = PathTo(pdf);
            root   = Root(pdf,".pdf");
            OUTEPS = strdup(Catenate(pwd,"/",root,".eps"));
          }
      }
    free(pwd);
    free(root);

    if (ispaf)
      read_paf(name,gzipd);
    else
      read_1aln(name);

    free(name);
  }

  ACHORD = get_selection_contigs(xseq, AGDB, AHASH, 1);
  BCHORD = get_selection_contigs(yseq, BGDB, BHASH, 1);

  aln_filter();

  if (OUTEPS != NULL)
    { foeps = fopen(OUTEPS,"w");
      if (foeps == NULL)
        { fprintf(stderr,"%s: Could not open file %s for writing\n",Prog_Name,OUTEPS);
          exit (1);
        }
    }
  else
    foeps = stdout;

  make_plot(foeps);

  if (foeps != stdout)
    fclose(foeps);
  
  if (OUTEPS != NULL)
    { char cmd[4096];

      sprintf(cmd,"%s %s",pdftool,OUTEPS);
      run_system_cmd(cmd, 1);
      sprintf(cmd,"rm -f %s",OUTEPS);
      run_system_cmd(cmd, 1);
    }
  
  free(BCHORD);
  free(ACHORD);

  Free_Hash_Table(AHASH);
  if (ISTWO)
    Free_Hash_Table(BHASH);

  free(segments);
  free(OUTEPS);
  free(Command_Line);
  free(Prog_Name);

  Catenate(NULL,NULL,NULL,NULL);

  exit (0);
}
