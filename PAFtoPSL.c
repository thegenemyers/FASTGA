/*******************************************************************************************
 *
 *  Utility for converting PAF to PSL format 
 *
 *  Author:    Chenxi Zhou
 *  Creation:  Nov 2024
 *  Last Mod:  Nov 2024
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pthread.h>

#include "gene_core.h"

#undef DEBUG_THREADS

#define VERSION "0.1"

static char *Usage[] = 
            { "[-T<int(8)>] [-C<str(cg:Z:)>] <alignments:path>[.paf]" };

static char *CIGAR_TAG = "cg:Z:";

typedef struct {
  int64  matches;
  int64  misMatches;
  int64  repMatches;
  int64  nCount;
  int64  qNumInsert;
  int64  qBaseInsert;
  int64  tNumInsert;
  int64  tBaseInsert;
  int    strand;
  char  *qname;
  int64  qSize;
  int64  qStart;
  int64  qEnd;
  char  *tname;
  int64  tSize;
  int64  tStart;
  int64  tEnd;
  int64  blockCount;
  int64 *blockSizes;
  int64 *blockStartQ;
  int64 *blockStartT;
  int64  capacity;
} PSL_Bundle;

static int push_align_block(PSL_Bundle *bundle, int64 qpos, int64 tpos, int64 clen)
{ if (bundle->capacity <= bundle->blockCount)
    { bundle->capacity = bundle->capacity * 1.2 + 100;
      bundle->blockSizes  = (int64 *) Realloc(bundle->blockSizes,  sizeof(int64)*bundle->capacity, "Reallocating block size array");
      bundle->blockStartQ = (int64 *) Realloc(bundle->blockStartQ, sizeof(int64)*bundle->capacity, "Reallocating block qPos array");
      bundle->blockStartT = (int64 *) Realloc(bundle->blockStartT, sizeof(int64)*bundle->capacity, "Reallocating block tPos array");
      if (bundle->blockSizes == NULL || bundle->blockStartQ == NULL || bundle->blockStartT == NULL)
        return (1);
    }
  bundle->blockStartQ[bundle->blockCount]  = qpos;
  bundle->blockStartT[bundle->blockCount]  = tpos;
  bundle->blockSizes[bundle->blockCount++] = clen; 
  return (0);
}

static int cigar2psl(char *cigar, PSL_Bundle *bundle)
{ int64 qNumInsert, qBaseInsert, tNumInsert, tBaseInsert;
  int64 qpos, tpos, insl, insr, lens, clen;
  char *c, p;
  int i;

  bundle->blockCount = 0;

  qNumInsert = 0;
  qBaseInsert = 0;
  tNumInsert = 0;
  tBaseInsert = 0;
  qpos = 0;
  tpos = 0;
  insl = 0;
  insr = 0;
  lens = 0;
  c = cigar;
  p = '\0';
  while (*c != '\0')
    { clen = 0;
      while (isdigit(*c))
        clen = 10*clen + (*c++-'0');
      if (clen == 0)
        { fprintf(stderr,"%s: CIGAR operator length is zero '%c'\n",Prog_Name,c[-1]);
          return (1);
        }
      switch ((*c++))
        { case 'M':
          case 'X':
          case '=':
            qpos += clen;
            tpos += clen;
            lens += clen;
            break;
          case 'I':
            /***
            if (p == 'I' || p == 'D')
              { fprintf(stderr,"%s: CIGAR error: operator '%c' followed by '%c'\n",Prog_Name,p,c[-1]);
                return (1);
              }
            **/
            if (p == '\0') // leading insertion
              insl = clen;
            else
              { if (push_align_block(bundle, qpos-lens, tpos-lens, lens))
                  return (1);
                lens = 0;
              }
            qNumInsert += 1;
            qBaseInsert += clen;
            qpos += clen;
            break;
          case 'D':
            /***
            if (p == 'I' || p == 'D')
              { fprintf(stderr,"%s: CIGAR error: operator '%c' followed by '%c'\n",Prog_Name,p,c[-1]);
                return (1);
              }
            **/
            if (p == '\0') // leading deletion
              insl = -clen;
            else
              { if (push_align_block(bundle, qpos-lens, tpos-lens, lens))
                  return (1);
                lens = 0;
              }
            tNumInsert += 1;
            tBaseInsert += clen;
            tpos += clen;
            break;
          default:
            fprintf(stderr,"%s: Invalid CIGAR operator '%c'\n",Prog_Name,c[-1]);
            return (1);
        }
      p = c[-1];
    }
  
  // the last CIGAR block
  if (p == 'I') // trailing insertion
    insr = clen;
  else if (p == 'D') // trailing deletion
    insr = -clen;
  else if (push_align_block(bundle, qpos-lens, tpos-lens, lens))
    return (1);
      
  if (qpos != bundle->qEnd - bundle->qStart)
    { fprintf(stderr,"%s: CIGAR length does not match aignment length (query): %lld!=%lld\n",
              Prog_Name,bundle->qEnd-bundle->qStart,qpos);
      return (1);
    }

  if (tpos != bundle->tEnd - bundle->tStart)
    { fprintf(stderr,"%s: CIGAR length does not match aignment length (target): %lld!=%lld\n",
              Prog_Name,bundle->tEnd-bundle->tStart,tpos);
      return (1);
    }
  
  // handle leading and trailing INDELs
  if (insl > 0) 
    { qNumInsert -= 1;
      qBaseInsert -= insl;
      bundle->qStart += insl;
    }
  else if (insl < 0)
    { tNumInsert -= 1;
      tBaseInsert += insl;
      bundle->tStart -= insl;
    } 

  if (insr > 0) 
    { qNumInsert -= 1;
      qBaseInsert -= insr;
      bundle->qEnd -= insr;
    }
  else if (insr < 0)
    { tNumInsert -= 1;
      tBaseInsert += insr;
      bundle->tEnd += insr;
    }
  
  // update target alignment block positions
  for (i = 0; i < bundle->blockCount; i++)
    bundle->blockStartT[i] += bundle->tStart;

  // update query alignment block positions
  // adjust positions for negative-strand matches
  if (bundle->strand) 
    for (i = 0; i < bundle->blockCount; i++)
      bundle->blockStartQ[i] = bundle->qSize - bundle->qEnd + bundle->blockStartQ[i];
  else
    for (i = 0; i < bundle->blockCount; i++)
      bundle->blockStartQ[i] += bundle->qStart;

  // insertions
  bundle->qNumInsert  = qNumInsert;
  bundle->qBaseInsert = qBaseInsert;
  bundle->tNumInsert  = tNumInsert;
  bundle->tBaseInsert = tBaseInsert;

  // mismatches
  bundle->misMatches = bundle->qEnd - bundle->qStart - qBaseInsert - bundle->matches;
  if (bundle->misMatches < 0)
    { fprintf(stderr,"%s: negative value misMatches: %lld\n",Prog_Name,bundle->misMatches);
      return (1);
    }
  
  // nCount
  lens = 0;
  for (i = 0; i < bundle->blockCount; i++)
    lens += bundle->blockSizes[i];
  bundle->nCount = lens - bundle->matches - bundle->misMatches - bundle->repMatches;
  if (bundle->nCount < 0)
    { fprintf(stderr,"%s: negative value nCount: %lld\n",Prog_Name,bundle->nCount);
      return (1);
    }

  return (0);
}

typedef struct
  { int64    beg;
    int64    end;
    char    *iname;
    FILE    *out;
    int      ret;
  } Packet;

//  Read next line into a buffer and return a pointer to the buffer
//    the length of the line.  NB: replaces '\n' with '\0'.

typedef struct
  { char *buffer;
    int   bmax;
  } Line_Bundle;

static char *read_line(void *input, char *name, Line_Bundle *lb)
{ char *buffer = lb->buffer;
  int   bmax   = lb->bmax;
  int len;

  if (bmax == 0)
    { bmax = 500;
      buffer = (char *) malloc(bmax);
      if (buffer == NULL)
        { fprintf(stderr,"%s: Out of memory reading %s\n",Prog_Name,name);
          return (NULL);
        }
    }

  if (fgets(buffer,bmax,input) == NULL)
    { if (feof(input))
        return (NULL);
      fprintf(stderr,"%s: Could not read next line of file %s (offset %lld)\n",
                     Prog_Name,name,(int64)ftello(input));
      return (NULL);
    }

  len = strlen(buffer);
  while (buffer[len-1] != '\n')
    { bmax = ((int) (1.4*bmax)) + 100;
      buffer = (char *) realloc(buffer,bmax);
      if (buffer == NULL)
        { fprintf(stderr,"%s: Out of memory reading %s\n",Prog_Name,name);
          return (NULL);
        }
      if (fgets(buffer+len,bmax-len,input) == NULL)
        { if (feof(input))
            fprintf(stderr,"%s: Last line of file %s does not end with new-line\n",
                           Prog_Name,name);
          else
            fprintf(stderr,"%s: Could not read next line of file %s (offset %lld)\n",
                           Prog_Name,name,(int64)ftello(input));
          return (NULL);
        }
      len += strlen(buffer+len);
    }
  buffer[--len] = '\0';

  lb->bmax   = bmax;
  lb->buffer = buffer;

  return (buffer);
}

static inline char *parse_cigar(char *s)
{ int l;
	while (*s != '\0') {
    while (isspace(*s)) s++;
    l = 0;
    while (s[l] != '\0' && !isspace(s[l]))
      l++;
    if (l > 5 && !strncmp(CIGAR_TAG,s,5))
      { s[l] = '\0';
        return (s+5);
      }
    s += l;
  }
  return (NULL);
}

void *gen_psl(void *args)
{ Packet  *parm  = (Packet *) args;
  char    *iname = parm->iname;
  FILE    *out = parm->out;
  
  Line_Bundle  _lb,  *lb = &_lb;
  PSL_Bundle   _psl, *psl = &_psl;

  int64  span, amnt;
  FILE  *input;
  int    nfields, linelen;
  char  *fptrs[11], *fptr, *eptr;

  int    i;

  input = fopen(iname,"r");
  if (input == NULL)
    { fprintf(stderr,"%s: Could not open %s for reading\n",Prog_Name,iname);
      return ((void *)1);
    }
  
  fseeko(input,parm->beg,SEEK_SET);
  span = parm->end-parm->beg;


  memset(lb, 0, sizeof(Line_Bundle));
  if (parm->beg > 0)
    amnt = strlen(read_line(input,iname,lb)) + 1;
  else
    amnt = 0;
  
  memset(psl, 0, sizeof(PSL_Bundle));
  psl->capacity    = 1<<16;
  psl->blockSizes  = Malloc(sizeof(int64)*psl->capacity,"Allocating blockSizes array");
  psl->blockStartQ = Malloc(sizeof(int64)*psl->capacity,"Allocating blockStartQ array");
  psl->blockStartT = Malloc(sizeof(int64)*psl->capacity,"Allocating blockStartT array");

  while (amnt <= span)
    { eptr = read_line(input,iname,lb);
      if (eptr == NULL)
        break;
      linelen = strlen(eptr) + 1;
      amnt   += linelen;
      
      fptr = eptr;
      while (*eptr != '\0' && isspace(*eptr))
        eptr += 1;
      nfields = 0;
      while (*eptr != '\0')
        { fptrs[nfields++] = eptr;
          while (*eptr != '\0' && !isspace(*eptr))
            eptr += 1;
          *eptr++ = '\0';
          while (isspace(*eptr))
            eptr += 1;
          if (nfields == 11)
            break;
        }
      if (nfields < 11)
        { fprintf(stderr,"%s: Line of paf has fewer than 11 fields (offset %lld)\n",
                         Prog_Name,(int64)ftello(input)-linelen);
          return ((void *)1);
        }
      
      psl->matches    = strtoll(fptrs[9], NULL, 10);
      //psl->misMatches = strtoll(fptrs[10], NULL, 10); // this is accutually alignment length
      psl->strand     = fptrs[4][0] == '+'? 0 : 1;
      psl->qname      = fptrs[0];
      psl->qSize      = strtoll(fptrs[1], NULL, 10);
      psl->qStart     = strtoll(fptrs[2], NULL, 10);
      psl->qEnd       = strtoll(fptrs[3], NULL, 10);
      psl->tname      = fptrs[5];
      psl->tSize      = strtoll(fptrs[6], NULL, 10);
      psl->tStart     = strtoll(fptrs[7], NULL, 10);
      psl->tEnd       = strtoll(fptrs[8], NULL, 10);

      eptr = parse_cigar(eptr);
      if (eptr == NULL)
        { fprintf(stderr,"%s: PAF line is missing a CIGAR string (offset %lld)\n",
                  Prog_Name,(int64)ftello(input)-linelen);
          return ((void *)1);
        }
      
      if (cigar2psl(eptr, psl))
        { fprintf(stderr,"%s: PAF record parsing error: ", Prog_Name);
          for (i = 0; i < linelen; i++)
            if (fptr[i] == '\0')
              fptr[i] = '\t';
          fprintf(stderr,"%s\n",fptr);
          // return ((void *)1);
          continue;
        }
      
      fprintf(out,"%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%c\t%s\t%lld\t%lld\t%lld\t%s\t%lld\t%lld\t%lld\t%lld\t",
                      psl->matches,
                      psl->misMatches,
                      psl->repMatches,
                      psl->nCount,
                      psl->qNumInsert,
                      psl->qBaseInsert,
                      psl->tNumInsert,
                      psl->tBaseInsert,
                      psl->strand? '-' : '+',
                      psl->qname,
                      psl->qSize,
                      psl->qStart,
                      psl->qEnd,
                      psl->tname,
                      psl->tSize,
                      psl->tStart,
                      psl->tEnd,
                      psl->blockCount);
      
      for (i = 0; i < psl->blockCount; i++)
        fprintf(out,"%lld,",psl->blockSizes[i]);
      fprintf(out,"\t");

      for (i = 0; i < psl->blockCount; i++)
        fprintf(out,"%lld,",psl->blockStartQ[i]);
      fprintf(out,"\t");

      for (i = 0; i < psl->blockCount; i++)
        fprintf(out,"%lld,",psl->blockStartT[i]);
      fprintf(out,"\n");
    }

  free(psl->blockSizes);
  free(psl->blockStartQ);
  free(psl->blockStartT);
  free(lb->buffer);
  fclose(input);

  return (NULL);
}

int main(int argc, char *argv[])
{ Packet    *parm;
  char      *APATH, *AROOT;

  int   NTHREADS;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char*  eptr;

    ARG_INIT("PAFtoPSL")

    (void) flags;

    NTHREADS = 8;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
          case 'C':
            CIGAR_TAG = argv[i]+2;
            if (strlen(CIGAR_TAG) != 5 ||CIGAR_TAG[2] != ':' || CIGAR_TAG[4] != ':')
              { fprintf(stderr,"%s: Not a valid tag name: %s\n",Prog_Name,CIGAR_TAG);
                exit (1);
              }
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 2)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        fprintf(stderr,"      -C: Cigar tag in the PAF file.\n");
        exit (1);
      }

    parm = Malloc(sizeof(Packet)*NTHREADS,"Allocating thread records");
    memset(parm, 0, sizeof(Packet)*NTHREADS);
    if (parm == NULL)
      exit (1);
  }

  { char       *inName, *oprefix, *buffer;
    struct stat state;
#ifndef DEBUG_THREADS
    pthread_t   threads[NTHREADS];
#endif
    int         p;

    APATH  = PathTo(argv[1]);
    AROOT  = Root(argv[1],".paf");
    inName = Strdup(Catenate(APATH,"/",AROOT,".paf"),"Allocating input name");
    if (inName == NULL)
      exit (1);

    if (stat(inName,&state) < 0)
      { fprintf(stderr,"%s: Cannot open %s to get size\n",Prog_Name,inName);
        exit (1);
      }

    oprefix = Strdup(Numbered_Suffix("_psl.",getpid(),"."),"Allocating temp file prefix");

    for (p = 0; p < NTHREADS; p++)
      { parm[p].beg   = (state.st_size * p) / NTHREADS;
        parm[p].end   = (state.st_size * (p+1)) / NTHREADS;
        parm[p].iname = inName;
        parm[p].out = fopen(Numbered_Suffix(oprefix,p,".psl"),"w+");
        if (parm[p].out == NULL)
          { fprintf(stderr,"%s: Cannot open %s.%d.paf for reading & writing\n",Prog_Name,oprefix,p);
            exit (1);
          }
        unlink(Numbered_Suffix(oprefix,p,".psl"));
      }
    
    Numbered_Suffix(NULL,0,NULL);
    free(oprefix);

    //  Launch and then gather threads

#ifdef DEBUG_THREADS
    for (p = 0; p < NTHREADS; p++)
      parm[p].ret = gen_psl(parm+p);
#else
    for (p = 1; p < NTHREADS; p++)
      pthread_create(threads+p,NULL,gen_psl,parm+p);
    parm[0].ret = (int)(intptr_t)gen_psl(parm);
    for (p = 1; p < NTHREADS; p++)
      parm[p].ret = pthread_join(threads[p],NULL);
#endif

    for (p = 0; p < NTHREADS; p++)
      if (parm[p].ret)
        { fprintf(stderr,"%s: conversion unsuccessful\n",Prog_Name);
          exit (1);
        }

    //  Concatenate all result files
#define PUSH_BLOCK 0x100000

    buffer = Malloc(PUSH_BLOCK,"Allocating seek block");

    for (p = 0; p < NTHREADS; p++)
      { FILE *o = parm[p].out;

        rewind(o);
        int x = PUSH_BLOCK;
        while (x >= PUSH_BLOCK)
          { x = fread(buffer,1,PUSH_BLOCK,o);
            if (x > 0)
              fwrite(buffer,x,1,stdout);
          }
        fclose(o);
      }

    free(buffer);
    free(parm);
    free(inName);
    free(AROOT);
    free(APATH);
  }

  free(Command_Line);
  free(Prog_Name);
  Catenate(NULL,NULL,NULL,NULL);
 
  exit (0);
}
