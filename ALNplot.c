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

/********************************** Revision History *****************************
 *                                                                               *
 * 16/04/24 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/

#ifdef __linux__
#define _GNU_SOURCE  // needed for vasprintf() on Linux
#endif

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

#include "DB.h"
#include "align.h"
#include "alncode.h"

#undef DEBUG_MAKE_DICT
#undef DEBUG_READ_1ALN
#undef DEBUG_PARSE_SEQ
#undef DEBUG_AXIS_CONF

static char *Usage[] =
  { "[-dSL] [-T<int(1)>] [-a<int(100)>] [-e<float(0.7)>] [-f<int>] [-t<float>]",
    "[-n<int(100000)>] [-H<int(600)>] [-W<int>] [-p[:<output:path>[.pdf]]]",
    "[-x<range>[,<range>[,...]]|:<FILE>] [-y<range>[,<range>[,...]]|:<FILE>]",
    "<alignment:path>[.1aln|.paf[.gz]]>"
  };

static int    MINALEN  = 100;
static int    IMGWIDTH = 0;
static int    IMGHEIGH = 0;
static int    FONTSIZE = 0;
static double LINESIZE = 0;
static int    MAXALIGN = 100000;
static int    NOLABEL  = 0;
static int    PRINTSID = 0;
static int    TRYADIAG = 0;
static double MINAIDNT = 0.7;
static int    NTHREADS = 1;
static char  *OUTEPS   = NULL;

#define SCAF_SYMBOL '@'
#define LAST_SYMBOL '#'
#define FILE_SYMBOL ':'
#define DLIM_SYMBOL ','

#define MAX_XY_LEN 10000
#define MIN_XY_LEN 50

typedef struct
  { uint8 flag;
    int   aread, bread;
    int   abpos, bbpos;
    int   aepos, bepos;
  } Segment;

#define DEL_FLAG 0x1
#define COL_RED  0x2
#define COL_BLUE 0x4
#define COL_GRAY 0x8

#define IS_DELT(flag) ((flag)&DEL_FLAG)
#define SET_DELT(flag) ((flag)|=DEL_FLAG)

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

static char *PDFTOOLS[4] = {"pstopdf", "epstopdf", "ps2pdf", "eps2pdf"};

typedef struct
  { int64    beg;
    int64    end;
    OneFile *in;
    Segment *segs;
    int64    nseg;
  } Packet;

/***********************************************************************************
 *
 *    DICTIONARY: string hash map
 *    adapted from Richard's implementation
 *
 **********************************************************************************/

typedef struct
  { char **names;
    uint64 *table;
    uint64 max;         /* current number of entries */
    uint64 dim;
    uint64 size;        /* 2^dim = size of tables */
    uint64 new;         /* communication between dictFind() and dictAdd() */
  } DICT;

static int       ISTWO;      // If two DBs are different

//  for 1aln contig to scaffold mapping
static int        NASCF;     // Number scaffolds
static int        NACTG;     // Number contigs
static DICT      *ADIC;      // Scaffold dictionary
static int       *ASCF;      // Scaffold to contig map - first contig
static int       *AMAP;      // Contig to scaffold map
static int       *ALEN;      // Contig length
static int       *AOFF;      // Contig offset map
static int       *ASEQ;      // Contig to plot
static int       *ABEG;      // Contig beg position
static int       *AEND;      // Contig end position
static char      *ANAME;     // User-defined sequence names

static int        NBSCF;     // Number scaffolds
static int        NBCTG;     // Number contigs
static DICT      *BDIC;      // Scaffold dictionary
static int       *BSCF;      // Scaffold to contig map - first contig
static int       *BMAP;      // Contig to scaffold map
static int       *BLEN;      // Contig length
static int       *BOFF;      // Contig offset map
static int       *BSEQ;      // Contig to plot
static int       *BBEG;      // Contig beg position
static int       *BEND;      // Contig end position
static char      *BNAME;     // User-defined sequence names

//  Parse read range grabbing up to 4 values from it as follows:
//    type 1
//          val[0][.val[1]]
//    type 2
//          val[0][.val[1]] - val[2][.val[3]]
//    type 3
//          val[0][.val[1]] _ val[2] - val[3]
//  Return 0 if there is a syntactic error, otherwise return the type.
//  0 values indicate $, -1 indicates not present, -2 (only in val[1] or val[3])
//      indicates val[0] or val[2], respectively, is a scaffold index

int dictFind(DICT *dict, char *s, uint64 *index);

static char *white(char *x)
{ while (isspace(*x))
    x += 1;
  return (x);
}

static char *getint(char *x, int *v)
{ *v = 0;
  while (isdigit(*x))
    *v = 10*(*v) + (*x++ - '0');
  return (x);
}

static char *address(char *x, int *v, int sep)
{ int a;

  x = white(x);
  a = *x++;
  if (a == LAST_SYMBOL)
    v[0] = 0;
  else if (isdigit(a))
    x = getint(x-1,v);
  else
    return (NULL);
  x = white(x);
  if (*x == sep)
    { x = white(x+1);
      a = *x++;
      if (a == LAST_SYMBOL)
        v[1] = 0;
      else if (isdigit(a))
        x = getint(x-1,v+1);
      else
        return (NULL);
      x = white(x);
    }
  else if (sep == ' ')
    v[1] = -2;
  else if (sep == '.')
    v[1] = -1;
  else  // sep == '-'
    return (NULL);
  return (x);
}

static int range(char *src, int *v, DICT *sdict)
{ int   t, sep;
  char *x, *y, *z;
  uint64 s, r;

  v[2] = -1;
  x = white(src);
  if (*x == SCAF_SYMBOL)
    { x += 1;
      sep = ' ';
    }
  else
    sep = '.';
  y = address(x,v,sep);
  if (y == NULL)
    { if (sep == ' ')
        goto try_string;
      else
        return (0);
    }
  t = 1;
  if (*y == '-')
    { y = address(y+1,v+2,sep);
      t = 2;
    }
  else if (*y == '_')
    { y = address(y+1,v+2,'-');
      t = 3;
    }
  if (y == NULL || *y != '\0')
    { if (sep == ' ')
        goto try_string;
      else
        return (0);
    }
  return (t);

try_string:
  if (dictFind(sdict,x,&s))
    { v[1] = -2;
      v[0] = s+1;
      return (1);
    }
  y = rindex(x,'_');
  if (y != NULL)
    { *y = '\0';
      if (dictFind(sdict,x,&s))
        { z = address(y+1,v+2,'-');
          if (z != NULL && *z == '\0')
            { v[1] = -2;
              v[0] = s+1;
              *y = '_';
              return (3);
            }
        }
      *y = '_';
    }
  for (y = index(x,'-'); y != NULL; y = index(y+1,'-'))
    { *y = '\0';
      if (dictFind(sdict,x,&s) && dictFind(sdict,y+1,&r))
        { v[1] = v[3] = -2;
          v[0] = s+1;
          v[2] = r+1;
          *y = '-';
          return (2);
        }
      *y = '-';
    }
  return (0);
}

 //  Interpret value pair s.c into an absolute contig range (1-based)

static int interpret_address(int s, int c, int ncontig, int nscaff, int *map)
{ if (c == -1)
    { if (s == 0)
        return (ncontig);
      if (s > ncontig)
        { fprintf(stderr,"%s: Absolute contig index '%d' is out of bounds\n",Prog_Name,s);
          exit (1);
        }
      return (s);
    }
  if (s == 0)
    s = nscaff;
  else if (s > nscaff)
    { fprintf(stderr,"%s: Scaffold index '%d' is out of bounds\n",Prog_Name,s);
      exit (1);
    }
  if (c == 0)
    return (map[s]);
  if (c == -2)
    return (s);
  if (c > map[s]-map[s-1])
    { fprintf(stderr,"%s: Relative contig index '%d' is out of bounds\n",Prog_Name,c);
      exit (1);
    }
  return (map[s-1] + c);
}

  //  Interpret argument x as a range and map to contig level ranges

static void contig_range(char *x, int ncontig, int nscaff, DICT *sdict, int *smap,
        int *clen, int *coff, int *pbeg, int *pend, int *pfst, int *plst)
{ int t, vals[4], len, sunit;
  int beg, end, fst, lst;

  memset(vals,0,sizeof(int)*4);
  t = range(x,vals,sdict);
  if (t == 0)
    { fprintf(stderr,"%s: Command line argument %s not a valid read range\n",Prog_Name,x);
      exit (1);
    }

  sunit = (vals[1] <= -2);

  beg = end = interpret_address(vals[0],vals[1],ncontig,nscaff,smap) - 1;
  fst = 0;
  if (t == 2)
    { end = interpret_address(vals[2],vals[3],ncontig,nscaff,smap) - 1;
      if (beg > end)
        { fprintf(stderr,"%s: Read range in '%s' is empty\n",Prog_Name,x);
          exit (1);
        }
    }

  if (sunit)
    { beg = smap[beg];
      end = smap[end+1]-1;
    }
  if (t <= 2)
    lst = clen[end];
  else // t == 3
    { fst = vals[2]; //  Type 3: intepret vals[2] and vals[3] as the substring interval
      lst = vals[3];
      if (sunit)
        len = coff[end] + clen[end];
      else
        len = clen[beg];
      if (lst == 0)
        lst = len;
      else if (lst > len)
        { fprintf(stderr,"%s: Substring interval in '%s' is out of bounds\n",Prog_Name,x);
          exit (1);
        }
      if (fst >= lst)
        { fprintf(stderr,"%s: Substring interval in '%s' is empty\n",Prog_Name,x);
          exit (1);
        }
      if (sunit)
        { for ( ; beg <= end; beg++)
            if (fst < coff[beg] + clen[beg])
              break;
          fst -= coff[beg];
          if (fst < 0)
            fst = 0;
          for (; end >= beg; end--)
            if (lst > coff[end])
              break;
          lst -= coff[end];
          if (lst > clen[end])
            lst = clen[end];
        }
      if (fst > 0) fst -= 1; // to support both 0-base and 1-base
    }
  *pbeg = beg;
  *pend = end;
  *pfst = fst;
  *plst = lst;
}

#define dictMax(dict)  ((dict)->max)

static void *remap(void *old, uint64 oldSize, uint64 newSize)
{ void *new = Malloc(newSize, "Allocating name array");
  memset(new, 0, newSize);
  memcpy(new, old, oldSize);
  free(old);
  return (new);
}

static uint64 hashString(char *cp, uint64 n, int isDiff)
{ uint64 i;
  uint64 j, x = 0;
  uint64 rotate = isDiff ? 21 : 13;
  uint64 leftover = 8 * sizeof(uint64) - rotate;

  while (*cp)
    x = (*cp++) ^ ((x >> leftover) | (x << rotate));

  for (j = x, i = n; i < sizeof(uint64); i += n)
    j ^= (x >> i);
  j &= (1 << n) - 1;

  if (isDiff)
    j |= 1;

  return (j);
}

DICT *dictCreate(uint64 size)
{ DICT *dict = (DICT *) Malloc (sizeof(DICT), "Allocating dictionary");
  memset(dict, 0, sizeof(DICT));
  for (dict->dim = 10, dict->size = 1024; dict->size < size; dict->dim += 1, dict->size *= 2);
  dict->table = (uint64 *) Malloc(sizeof(uint64) * dict->size, "Allocating dictionary table");
  memset(dict->table, 0, sizeof(uint64) * dict->size);
  dict->names = (char **) Malloc (sizeof(char *) * dict->size / 2, "Allocating dictionary names");
  memset(dict->names, 0, sizeof(char *) * dict->size / 2);
  return (dict);
}

void dictDestroy(DICT *dict)
{ uint64 i;
  for (i = 1; i <= dict->max; ++i) free(dict->names[i]);
  free(dict->names);
  free(dict->table);
  free(dict);
}

int dictFind(DICT *dict, char *s, uint64 *index)
{ uint64 i, x, d;
  
  if (!dict || !s) return (0);

  x = hashString (s, dict->dim, 0);
  if (!(i = dict->table[x]))
    { dict->new = x;
      return (0);
    }
  else if (!strcmp (s, dict->names[i]))
    { if (index)
        *index = i-1;
      return (1);
    }
  else
    { d = hashString (s, dict->dim, 1);
      while (1)
        { x = (x + d) & ((1 << dict->dim) - 1);
          if (!(i = dict->table[x]))
            { dict->new = x;
              return (0);
            }
          else if (!strcmp (s, dict->names[i]))
          {
            if (index)
              *index = i-1;
            return (1);
          }
        }
    }
  return (0);
}

int dictAdd(DICT *dict, char *s, uint64 *index)
{ uint64 i, x;

  if (dictFind(dict, s, index)) return (0);

  i = ++dict->max;
  dict->table[dict->new] = i;
  dict->names[i] = (char*) Malloc(strlen(s) + 1, "Allocating name space");
  strcpy(dict->names[i], s);
  if (index) *index = i-1;

  if (dict->max > 0.3 * dict->size) /* double table size and remap */
    { uint64 *newTable;
      dict->dim += 1;
      dict->size *= 2;
      dict->names = (char **) remap(dict->names, sizeof(char *) * (dict->max + 1),
                                                 sizeof(char *) * (dict->size / 2));
      newTable = (uint64 *) Malloc(sizeof(uint64) * dict->size, "Allocating new table");
      memset(newTable, 0, sizeof(uint64) * dict->size);
      for (i = 1; i <= dict->max; i++)
        { s = dict->names[i];
          x = hashString (s, dict->dim, 0);
         if (!newTable[x])
           newTable[x] = i;
         else
           { uint64 d = hashString (s, dict->dim, 1);
             while (1)
               { x = (x + d) & ((1 << dict->dim) - 1);
                 if (!newTable[x])
                   { newTable[x] = i;
                     break;
                   }
               }
           }
        }
      free(dict->table);
      dict->table = newTable;
    }
  
  return (1);
}

char *dictName(DICT *dict, uint64 i)
{ return (dict->names[i+1]); }

/**********************************************************************************/

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

// User can  provide a file input for more complicated configuration
// The first three columns of each line will be read
// <range>    +/-    name
// to specify plotting range, orientation and user defined sequence name
// if the sequence name is specified, i. e., three columns in a line
// then the <range> in that line has to be a contiguous scaffold block

static int parseTargetSEQFromFile(char *file, int ncontig, int nscaff, DICT *sdict, int *smap,
        int *cmap, int *clen, int *coff, int **_SEQ, int **_BEG, int **_END, char **_NAME)
{ if (file == NULL || *file == '\0')
    { fprintf(stderr,"%s: Empty range file\n",Prog_Name);
      exit (1);
    }

  FILE *fp;
  char c, *p, *q, *line;
  size_t ln = 0;
  ssize_t read;
  int i, pbeg, pend, pfst, plst, rank;
  int *SEQ, *BEG, *END;
  char *NAME;
  uint32 nstr, mstr, lstr;

  SEQ = *_SEQ;
  BEG = *_BEG;
  END = *_END;

  fp = fopen(file, "r");
  if (fp == NULL)
    { fprintf(stderr,"%s: Cannot open file %s\n",Prog_Name,file);
      exit (1);
    }

  nstr = 0;
  mstr = 2048;
  NAME = (char *)Malloc(sizeof(char)*mstr,"Allocating name array");
  if (NAME == NULL)
    exit (1);
  *NAME = '\0';

  line = NULL;
  rank = 0;
  while ((read = getline(&line, &ln, fp)) != -1)
    { p = q = line;
      while (isspace(*q)) ++q;
      p = q;
      while (*q != '\0' && !isspace(*q)) ++q;
      if (p == q) continue; // empty line
      c = *q;
      *q = '\0';
      contig_range(p, ncontig, nscaff, sdict, smap, clen, coff, &pbeg, &pend, &pfst, &plst);
#ifdef DEBUG_PARSE_SEQ
      fprintf(stderr, "%s: CTG_RANGE pbeg = %-12d pend = %-12d pfst = %-12d plst = %-12d\n",
              Prog_Name,pbeg,pend,pfst,plst);
#endif
      rank++;
      for (i = pbeg; i <= pend; i++)
        { if (SEQ[i])
            { fprintf(stderr,"%s: Overlapped sequence range\n",Prog_Name);
              exit (1);
            }
          SEQ[i] = rank;
          BEG[i] = 0;
          END[i] = clen[i];
        }
      BEG[pbeg] = pfst;
      END[pend] = plst;
      *q = c;

      while (isspace(*q)) ++q;
      p = q;
      while (*q != '\0' && !isspace(*q)) ++q;
      if (p == q) continue; // no second column
      if (*p == '-') // in reverse order
        for (i = pbeg; i <= pend; i++)
          SEQ[i] = -rank;

      while (isspace(*q)) ++q;
      p = q;
      while (*q != '\0' && !isspace(*q)) ++q;
      if (p == q) continue; // no third column
      // check if sequence is contiguous
      if (cmap[pbeg] != cmap[pend])
        { fprintf(stderr,"%s: Cannot set sequence name for noncontiguous range: %s\n",Prog_Name,line);
          exit (1);
        }
      // add rank sequence name
      *q = '\0';
      lstr = Number_Digits(rank) + q - p + 2;
      expandBuffer(&NAME, &nstr, &mstr, lstr+1);
      sprintf(NAME+nstr, "%d %s ", rank, p);
      nstr += lstr;
    }

#ifdef DEBUG_PARSE_SEQ
  if (*NAME == '\0')
    fprintf(stderr,"%s: Customized sequence name buffer is empty\n",Prog_Name);
  else
    fprintf(stderr,"%s: Customized sequence name buffer - %s\n",Prog_Name,NAME);
#endif

  if (*NAME == '\0')
    free(NAME);
  else
    *_NAME = NAME;

  free(line);
  fclose(fp);

  return (0);
}

int parseTargetSEQ(char *x, int ncontig, int nscaff, DICT *sdict, int *smap, int *cmap,
        int *clen, int *coff, int **_SEQ, int **_BEG, int **_END, char **_NAME)
{ char c, *eptr;
  int i, pbeg, pend, pfst, plst, rank;
  int *SEQ, *BEG, *END;

  SEQ = (int *) Malloc(sizeof(int)*3*ncontig,"Allocating sequence array");
  BEG = SEQ + ncontig;
  END = BEG + ncontig;

  *_SEQ  = SEQ;
  *_BEG  = BEG;
  *_END  = END;
  *_NAME = NULL;

  memset(SEQ, 0, sizeof(int)*ncontig);

  if (x == NULL || *x == '-')
    { for (i = 0; i < ncontig; i++)
        { SEQ[i] = 1;
          BEG[i] = 0;
          END[i] = clen[i];
        }
      goto print_range;
    }

  if (*x == '\0')
    { fprintf(stderr, "%s: Empty range\n",Prog_Name);
      exit (1);
    }
  if (*x == FILE_SYMBOL)
    { parseTargetSEQFromFile(x+1,ncontig,nscaff,sdict,smap,cmap,clen,coff,_SEQ,_BEG,_END,_NAME);
      goto print_range;
    }

  rank = 0;
  while (1)
    { eptr = x + 1;
      while (*eptr && *eptr != DLIM_SYMBOL) eptr++;
      c = *eptr;
      *eptr = '\0';
      contig_range(x, ncontig, nscaff, sdict, smap, clen, coff, &pbeg, &pend, &pfst, &plst);
#ifdef DEBUG_PARSE_SEQ
      fprintf(stderr, "%s: CTG_RANGE pbeg = %-12d pend = %-12d pfst = %-12d plst = %-12d\n",
              Prog_Name,pbeg,pend,pfst,plst);
#endif
      rank++;
      for (i = pbeg; i <= pend; i++)
        { if (SEQ[i])
            { fprintf(stderr,"%s: Overlapped sequence range\n",Prog_Name);
              exit (1);
            }
          SEQ[i] = rank;
          BEG[i] = 0;
          END[i] = clen[i];
        }
      BEG[pbeg] = pfst;
      END[pend] = plst;
      *eptr = c;
      if (!c) break;
      x = eptr + 1;
    }

print_range:
#ifdef DEBUG_PARSE_SEQ
  fprintf(stderr,"%s: Sequence range to plot\n",Prog_Name);
  fprintf(stderr,"%s: CONTIG RANK   BEGIN        END\n",Prog_Name);
  for (i = 0; i < ncontig; i++)
    { if (SEQ[i])
        fprintf(stderr,"%s: %-8d %-6d %-12d %-12d\n",Prog_Name,i,SEQ[i],BEG[i],END[i]);
    }
#endif

  return (0);
}

void makeSeqDICTFromDB(char *db1_name, char *db2_name)
{ DAZZ_DB   _db1, *db1 = &_db1;
  DAZZ_DB   _db2, *db2 = &_db2;
  int nascaff, nacontig;
  int nbscaff, nbcontig;

  //  Open DB or DB pair

  { int   s, r, gdb;

    ISTWO  = 0;
    gdb    = Open_DB(db1_name,db1);
    if (gdb < 0)
      { fprintf(stderr,"%s: Could not open DB file: %s\n",Prog_Name,db1_name);
        exit (1);
      }

    if (db1->nreads != db1->treads)
      { fprintf(stderr,"%s: Not the same %d %d\n",Prog_Name,db1->nreads,db1->treads);
        exit (1);
      }

    if (db2_name != NULL)
      { gdb = Open_DB(db2_name,db2);
        if (gdb < 0)
          { fprintf(stderr,"%s: Could not open DB file: %s\n",Prog_Name,db2_name);
            exit (1);
          }

        if (db2->nreads != db2->treads)
          { fprintf(stderr,"%s: Not the same %d %d\n",Prog_Name,db2->nreads,db2->treads);
            exit (1);
          }

        ISTWO = 1;
      }
    else
      db2  = db1;

    //  Build contig to scaffold maps in global vars
    nacontig = db1->nreads;
    nbcontig = db2->nreads;

    if (ISTWO)
      AMAP = (int *) Malloc(sizeof(int)*(4*(nacontig+nbcontig)+2),"Allocating scaffold map");
    else
      AMAP = (int *) Malloc(sizeof(int)*(4*nacontig+1),"Allocating scaffold map");

    if (AMAP == NULL)
      exit (1);

    AOFF   = AMAP + nacontig;
    ALEN   = AOFF + nacontig;
    ASCF   = ALEN + nacontig;

    s = -1;
    for (r = 0; r < nacontig; r++)
      { if (db1->reads[r].origin == 0)
          { s += 1;
            ASCF[s]  = r;
          }
        AMAP[r]  = s;
        AOFF[r]  = db1->reads[r].fpulse;
        ALEN[r]  = db1->reads[r].rlen;
      }
    nascaff = s + 1;
    ASCF[nascaff] = nacontig;

    if (ISTWO)
      { BMAP   = ASCF + nacontig + 1;
        BOFF   = BMAP + nbcontig;
        BLEN   = BOFF + nbcontig;
        BSCF   = BLEN + nbcontig;

        s = -1;
        for (r = 0; r < db2->treads; r++)
          { if (db2->reads[r].origin == 0)
              { s += 1;
                BSCF[s] = r;
              }
            BMAP[r]  = s;
            BOFF[r]  = db2->reads[r].fpulse;
            BLEN[r]  = db2->reads[r].rlen;
          }
        nbscaff = s + 1;
        BSCF[nbscaff] = nbcontig;
      }
    else
      { BMAP   = AMAP;
        BOFF   = AOFF;
        BLEN   = ALEN;
        BSCF   = ASCF;
        nbscaff = nascaff;
      }
  }

  //  Preload all scaffold headers and set header offset

  { int r, hdrs, absent;
    char *HEADER, *eptr;
    struct stat state;

    hdrs = open(Catenate(db1->path,".hdr","",""),O_RDONLY);
    if (hdrs < 0)
      { fprintf(stderr,"%s: Could not open header file of %s\n",Prog_Name,db1->path);
        exit (1);
      }
    if (fstat(hdrs,&state) < 0)
      { fprintf(stderr,"%s: Could not fetch size of %s's header file\n",Prog_Name,db1->path);
        exit (1);
      }

    HEADER = Malloc(state.st_size,"Allocating header table");
    if (HEADER == NULL)
      exit (1);

    if (read(hdrs,HEADER,state.st_size) < 0)
      { fprintf(stderr,"%s: Could not read header file of %s\n",Prog_Name,db1->path);
        exit (1);
      }
    close(hdrs);

    ADIC = dictCreate(nascaff);
    for (r = 0; r < nacontig; r++)
      { if (db1->reads[r].origin == 0)
          { for (eptr = HEADER + db1->reads[r].coff; *eptr != '\n'; eptr++)
              if (isspace(*eptr))
                break;
            *eptr = '\0';

            absent = dictAdd(ADIC, HEADER + db1->reads[r].coff, NULL);
            if (!absent)
              { fprintf(stderr,"%s: Duplicate sequence name: %s\n",Prog_Name,HEADER + db1->reads[r].coff);
                exit (1);
              }
          }
      }
    free(HEADER);

    if (ISTWO)
      { hdrs = open(Catenate(db2->path,".hdr","",""),O_RDONLY);
        if (hdrs < 0)
          { fprintf(stderr,"%s: Could not open header file of %s\n",Prog_Name,db2->path);
            exit (1);
          }
        
        if (fstat(hdrs,&state) < 0)
          { fprintf(stderr,"%s: Could not fetch size of %s's header file\n",Prog_Name,db2->path);
            exit (1);
          }

        HEADER = Malloc(state.st_size,"Allocating header table");
        if (HEADER == NULL)
          exit (1);

        if (read(hdrs,HEADER,state.st_size) < 0)
          { fprintf(stderr,"%s: Could not read header file of %s\n",Prog_Name,db2->path);
            exit (1);
          }
        close(hdrs);

        BDIC = dictCreate(nbscaff);
        for (r = 0; r < nbcontig; r++)
          { if (db2->reads[r].origin == 0)
              { for (eptr = HEADER + db2->reads[r].coff; *eptr != '\n'; eptr++)
                  if (isspace(*eptr))
                    break;
                *eptr = '\0';

                absent = dictAdd(BDIC, HEADER + db2->reads[r].coff, NULL);
                if (!absent)
                  { fprintf(stderr,"%s: Duplicate sequence name: %s\n",Prog_Name,HEADER + db2->reads[r].coff);
                    exit (1);
                  }
              }
          }

        free(HEADER);
      }
    else
      BDIC = ADIC;
  }

  NACTG = nacontig;
  NBCTG = nbcontig;
  NASCF = nascaff;
  NBSCF = nbscaff;

#ifdef DEBUG_MAKE_DICT
  { int i;
    fprintf(stderr,"%s: DB_A\n",Prog_Name);
    fprintf(stderr,"%s: List of contigs (NACTG=%d)\n",Prog_Name,NACTG);
    fprintf(stderr,"%s:  INDEX   SCAF     OFFSET     LENGTH\n",Prog_Name);
    for (i = 0; i < NACTG; i++)
      fprintf(stderr,"%s: %6d %6d %10d %10d\n",Prog_Name,i,AMAP[i],AOFF[i],ALEN[i]);
    fprintf(stderr,"%s: List of scaffolds (NASCF=%d)\n",Prog_Name,NASCF);
    for (i = 0; i < NASCF; i++)
      fprintf(stderr,"%s: %6d %s\n",Prog_Name,i,dictName(ADIC,i));
  }
  if (ISTWO)
    { int i;
      fprintf(stderr,"%s: DB_B\n",Prog_Name);
      fprintf(stderr,"%s: List of contigs (NBCTG=%d)\n",Prog_Name,NBCTG);
      fprintf(stderr,"%s:  INDEX   SCAF     OFFSET     LENGTH\n",Prog_Name);
      for (i = 0; i < NBCTG; i++)
        fprintf(stderr,"%s: %6d %6d %10d %10d\n",Prog_Name,i,BMAP[i],BOFF[i],BLEN[i]);
      fprintf(stderr,"%s: List of scaffolds (NBSCF=%d)\n",Prog_Name,NBSCF);
      for (i = 0; i < NBSCF; i++)
        fprintf(stderr,"%s: %6d %s\n",Prog_Name,i,dictName(BDIC,i));
    }
#endif

  Close_DB(db1);
  if (ISTWO)
    Close_DB(db2);
}

void *read_1aln_block(void *args)
{ Packet *parm  = (Packet *) args;
  int64    beg  = parm->beg;
  int64    end  = parm->end;
  OneFile *in   = parm->in;
  Segment *segs = parm->segs;
  int64    nseg = 0;

  //  For each alignment do

  if (!oneGotoObject (parm->in, beg))
    { fprintf(stderr,"%s: Could not locate to object %lld in 1aln file\n",Prog_Name,beg);
      exit (1);
    }
  
  oneReadLine(parm->in);

  { int64 i;
    uint32   flags;
    uint8    flag;
    int      aread, bread;
    int      abpos, bbpos;
    int      aepos, bepos;
    int      diffs;
    int      blocksum, iid;

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
          { bbpos = BLEN[bread] - bbpos;
            bepos = BLEN[bread] - bepos;
            SET_BLUE(flag);
          }
        else
          SET_RED(flag);

        // add to output
        segs->flag  = flag;
        segs->aread = aread;
        segs->abpos = abpos;
        segs->aepos = aepos;
        segs->bread = bread;
        segs->bbpos = bbpos;
        segs->bepos = bepos;
        segs++;
        nseg++;
      }
  }
  parm->nseg = nseg;

  return (NULL);
}

void read_1aln(char *oneAlnFile)
{ Packet    *parm;
  OneFile   *input;
  int64      novl;

  //  Initiate .1aln file reading and read header information

  { char      *cpath, *tmp;
    FILE      *test;
    char      *db1_name;
    char      *db2_name;
    int        TSPACE;
    
    input = open_Aln_Read(oneAlnFile,NTHREADS,&novl,&TSPACE,&db1_name,&db2_name,&cpath);
    if (input == NULL)
      { fprintf(stderr,"%s: Could not open .1aln file: %s\n",Prog_Name,oneAlnFile);
        exit (1);
      }

    test = fopen(db1_name,"r");
    if (test == NULL)
      { if (*db1_name != '/')
          test = fopen(Catenate(cpath,"/",db1_name,""),"r");

        if (test == NULL)
          { fprintf(stderr,"%s: Could not find .gdb %s\n",Prog_Name,db1_name);
            exit (1);
          }
        
        tmp = Strdup(Catenate(cpath,"/",db1_name,""),"Allocating expanded name");
        free(db1_name);
        db1_name = tmp;
      }
    fclose(test);
  
    if (db2_name != NULL)
      { test = fopen(db2_name,"r");
        if (test == NULL)
          { if (*db2_name != '/')
              test = fopen(Catenate(cpath,"/",db2_name,""),"r");
            
            if (test == NULL)
              { fprintf(stderr,"%s: Could not find .gdb %s\n",Prog_Name,db2_name);
                exit (1);
              }
            
            tmp = Strdup(Catenate(cpath,"/",db2_name,""),"Allocating expanded name");
            free(db2_name);
            db2_name = tmp;
          }
        fclose(test);
      }

    makeSeqDICTFromDB(db1_name, db2_name);

    free(cpath);
    free(db1_name);
    free(db2_name);
  }

  // Read alignment segments

  //  Divide .1aln into NTHREADS parts
  { int p;

    parm = Malloc(sizeof(Packet)*NTHREADS,"Allocating thread records");
    if (parm == NULL)
      exit (1);

    segments = (Segment *) Malloc(sizeof(Segment)*novl, "Allocating segment array");
    if (segments == NULL)
      { fprintf(stderr,"%s: Allocating segment array memory failed\n",Prog_Name);
        exit (1);
      }
    for (p = 0; p < NTHREADS ; p++)
      { parm[p].beg = (p * novl) / NTHREADS;
        if (p > 0)
          parm[p-1].end = parm[p].beg;
        parm[p].segs = segments + parm[p].beg;
        parm[p].nseg = 0;
      }
    parm[NTHREADS-1].end = novl;
  }

  //   Use NTHREADS to produce alignment segments for each part
  {
    int p;
    int64 totSeg;
    pthread_t threads[NTHREADS];
    for (p = 0; p < NTHREADS; p++)
      parm[p].in   = input + p;
    for (p = 1; p < NTHREADS; p++)
      pthread_create(threads+p,NULL,read_1aln_block,parm+p);
    read_1aln_block(parm);
    for (p = 1; p < NTHREADS; p++)
      pthread_join(threads[p],NULL);
  
    // collect results from different part
    nSegment = parm[0].nseg;
    totSeg   = parm[0].end - parm[0].beg;
    for (p = 1; p < NTHREADS; p++)
      { if (nSegment < totSeg)
          memmove(parm[p].segs,segments+nSegment,parm[p].nseg);
        nSegment += parm[p].nseg;
        totSeg   += parm[p].end - parm[p].beg;
      }
#ifdef DEBUG_READ_1ALN
    fprintf(stderr, "%s: %9lld segments loaded\n",Prog_Name,totSeg);
    fprintf(stderr, "%s: %9lld after filtering\n",Prog_Name,nSegment);
#endif
  }

  free(parm);

  oneFileClose(input);
}

#define LINE_SEP '\n'

typedef struct
  { uint64 n, m;
    char *buf;
  } Buffer;

Buffer *newBuffer(uint64 size)
{ Buffer *buffer;
  char *buf;

  if (size < 16) size = 16;

  buffer = (Buffer *) Malloc(sizeof(Buffer), "Allocating buffer memory");
  buf = (char *) Malloc(sizeof(char) * size, "Allocating buffer memory");
  
  if (buffer == NULL || buf == NULL)
    exit (1);

  buffer->n = 0;
  buffer->m = size;
  buffer->buf = buf;
 
  return (buffer);
}

void extendBuffer(Buffer *buffer)
{ buffer->m <<= 1;
  buffer->buf = (char *) Realloc(buffer->buf,sizeof(char)*buffer->m,"Extending buffer memory");
  if (buffer->buf == NULL)
    exit (1);
}

void destroyBuffer(Buffer *buffer)
{ free(buffer->buf);
  free(buffer);
}

int get_until(void *input, int gzipd, Buffer *buffer)
{ uint64 pos, len;
  int eof;
  char *buf;

  eof = 0;
  buffer->n = 0;
  buffer->buf[0] = '\0';
  pos = 0;
  while (1)
    { buf = buffer->buf;
      if (gzipd)
        eof = (gzgets(input,buf+pos,buffer->m-pos) == NULL);
      else
        eof = (fgets(buf+pos,buffer->m-pos,input) == NULL);

      for (len = pos; buf[len] != LINE_SEP && buf[len]; len++) {}

      if (!eof && len == buffer->m-1)
        { extendBuffer(buffer);
          pos = len;
        }
      else
        break;
    }

  return (eof);
}

void read_paf(char *pafAlnFile, int gzipd)
{ void *input;
  Buffer *buffer;
  int i, naseq, maseq, nbseq, mbseq, eof, absent;
  uint64 index, nsegs, msegs;
  Segment *segs;
  char *fptrs[11], *fptr, *eptr;
  int alen, blen, *alens, *blens;;
  int aread, bread;
  int abpos, bbpos;
  int aepos, bepos;
  int blocksum, iid;
  uint8 flag;

  if (gzipd)
    input = gzopen(pafAlnFile,"r");
  else
    input = fopen(pafAlnFile,"r");

  ADIC = dictCreate(1024);
  BDIC = dictCreate(1024);
  
  naseq = 0;
  maseq = 1024;
  nbseq = 0;
  mbseq = 1024;
  nsegs = 0;
  msegs = 4096;

  alens = (int *)Malloc(sizeof(int) * maseq, "Allocating length map");
  blens = (int *)Malloc(sizeof(int) * mbseq, "Allocating length map");
  segments = (Segment *)Malloc(sizeof(Segment) * msegs, "Allocating segment array");
  if (alens == NULL || blens == NULL || segments == NULL)
    exit (1);
  segs = segments;

  buffer = newBuffer(0);
  while (1)
    { eof = get_until(input, gzipd, buffer);
      // parse PAF line
      eptr = buffer->buf;
      fptr = eptr;
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
              break; 
            }
        }
      
      if (i == 11)
      { absent = dictAdd(ADIC, fptrs[0], &index);
        alen   = strtol(fptrs[1], &eptr, 10);
        if (absent)
          { alens[naseq++] = alen;
            if (naseq == maseq)
              { maseq <<= 1;
                alens = (int *) Realloc(alens, sizeof(int) * maseq, "Reallocating length map");
                if (alens == NULL)
                  exit (1);
              }
          }
        aread  = index;
        abpos  = strtol(fptrs[2], &eptr, 10);
        aepos  = strtol(fptrs[3], &eptr, 10);
      
        absent = dictAdd(BDIC, fptrs[5], &index);
        blen   = strtol(fptrs[6], &eptr, 10);
        if (absent)
          { blens[nbseq++] = blen;
            if (nbseq == mbseq)
              { mbseq <<= 1;
                blens = (int *) Realloc(blens, sizeof(int) * mbseq, "Reallocating length map");
                if (blens == NULL)
                  exit (1);
              }
          }
        bread  = index;
        bbpos  = strtol(fptrs[7], &eptr, 10);
        bepos  = strtol(fptrs[8], &eptr, 10);
        
        absent = 0;
        if (aepos - abpos < MINALEN || bepos - bbpos < MINALEN)
          absent = 1; // filter by segment size

        blocksum = (aepos-abpos) + (bepos-bbpos);
        iid      = strtol(fptrs[9], &eptr, 10);
        if (2.*iid / blocksum < MINAIDNT)
          absent = 1; // filter by segment identity

        if (!absent)
          { // Add to segment array
            flag = 0;
            if (*fptrs[4] == '-')
              { int tmp = bbpos;
                bbpos = bepos;
                bepos = tmp;
                SET_BLUE(flag);
              }
            else
              SET_RED(flag);

            segs->flag  = flag;
            segs->aread = aread;
            segs->bread = bread;
            segs->abpos = abpos;
            segs->aepos = aepos;
            segs->bbpos = bbpos;
            segs->bepos = bepos;
        
            nsegs++;
            if (nsegs == msegs)
              { msegs <<= 1;
                segments = (Segment *) Realloc(segments, sizeof(Segment) * msegs,
                                               "Reallocating segment array");
                if (segments == NULL)
                  exit (1);
              }
            segs = segments + nsegs;
          }
      }
      if (eof) break;
    }
 
  if (gzipd)
    gzclose(input);
  else
    fclose(input);

  ISTWO = 1;
  nSegment = nsegs;

  NASCF = NACTG = naseq;
  NBSCF = NBCTG = nbseq;

  AMAP = (int *) Malloc(sizeof(int)*(4*(naseq+nbseq)+2),"Allocating scaffold map");
  AOFF = AMAP + naseq;
  ALEN = AOFF + naseq;
  ASCF = ALEN + naseq;

  memcpy(ALEN, alens, sizeof(int)*naseq);
  for (i = 0; i < naseq; i++)
    { AMAP[i] = i;
      ASCF[i] = i;
    }
  ASCF[naseq] = naseq;
  AOFF[0] = 0;
  for (i = 1; i < naseq; i++)
    AOFF[i] = AOFF[i-1] + ALEN[i-1];

  BMAP = ASCF + naseq + 1;
  BOFF = BMAP + nbseq;
  BLEN = BOFF + nbseq;
  BSCF = BLEN + nbseq;

  memcpy(BLEN, blens, sizeof(int)*nbseq);
  for (i = 0; i < nbseq; i++)
    { BMAP[i] = i;
      BSCF[i] = i;
    }
  BSCF[nbseq] = nbseq;
  BOFF[0] = 0;
  for (i = 1; i < nbseq; i++)
    BOFF[i] = BOFF[i-1] + BLEN[i-1];

  free(alens);
  free(blens);
  destroyBuffer(buffer);
}

static int USORT(const void *l, const void *r)
{ uint64 *x = (uint64 *) l;
  uint64 *y = (uint64 *) r;

  return ((*x > *y) - (*x < *y));
}

static int DSORT(const void *l, const void *r)
{ uint64 *x = (uint64 *) l;
  uint64 *y = (uint64 *) r;
  
  return ((*x < *y) - (*x > *y));
}

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

static void addSeqName(char *names, uint32 n, uint32 m, int c0, int c1, 
        uint32 s, DICT *sdict, int *smap, int *cbeg, int *cend, int *clen,
        int rank, char **_snames, char **_names, uint32 *_n, uint32 *_m)
{ if (NOLABEL)
    return;

  int i;
  uint32 l;
  int64 p;
  char *name;

  if (*_snames != NULL)
    { // check if this rank has user-defined name
      char *eptr, *q;
      int r;
      q = *_snames;
      r = strtol(q, &eptr, 10);
      if (q != eptr && r == abs(rank))
        { // user-defined name found
          q = eptr;
          while (isspace(*q)) q++;
          eptr = q;
          while (*eptr != '\0' && !isspace(*eptr))
            eptr++;
          if (q != eptr)
            { // name is not empty
              l = eptr - q;
              expandBuffer(&names, &n, &m, l+1);
              snprintf(names+n,l+1,"%s",q);
              n += l;
              n++; // the null terminator
            }
          *_snames = eptr;
          goto assign_vars;
        }
    }

  if (PRINTSID)
    { l = Number_Digits(s+1);
      expandBuffer(&names, &n, &m, l+1);
      sprintf(names+n,"%u",s+1);
      n += l;
    }
  else
    { name = dictName(sdict,s);
      l = strlen(name);
      expandBuffer(&names, &n, &m, l+1);
      sprintf(names+n,"%s",name);
      n += l;
    }
  if (cbeg[c0] > 0 ||
          smap[s] != c0 ||      // start from first
          smap[s+1] != c1+1 ||
          cend[c1] != clen[c1]) // end at last
    { // partial scaffold
      p = 0;
      for (i = smap[s]; i < c0; i++)
        p += clen[i];
      p += cbeg[c0];
      p += 1;
      l = Number_Digits(p);
      l += 1;
      expandBuffer(&names, &n, &m, l+1);
      sprintf(names+n,"_%lld",p);
      n += l;
      p -= 1;
      p += cend[c0] - cbeg[c0];
      for (i = c0 + 1; i <= c1; i++)
        p += cend[i] - cbeg[i];
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
  n++; // the null terminator

assign_vars:
  *_names = names;
  *_n = n;
  *_m = m;
}

int axisConfig(DICT *sdict, int *smap, int nctg, int *cseq, int *cbeg, int *cend, int *cmap, int *clen,
        int *coff, char *snames, int64 **_caxis, int64 **_saxis, char **_names, int64 *_tseq)
{ int i, j, s, r1, r2, mseq, nseq;
  uint32 c0, c1, c2, nstr, mstr;
  uint64 *sarray;
  int64 *caxis, *saxis, tseq;
  char *names;

  mseq = 0;
  for (i = 0; i < nctg; i++)
    if (cseq[i])
      mseq++;

  sarray = (uint64 *)Malloc(sizeof(uint64)*mseq,"Allocating seq array");
  caxis = (int64 *)Malloc(sizeof(int64)*nctg,"Allocating offset array");
  saxis = (int64 *)Malloc(sizeof(int64)*mseq,"Allocating offset array");
  if (sarray == NULL || caxis == NULL || saxis == NULL)
    exit (1);

  if (!NOLABEL)
    { mstr = 2048;
      names = (char *) Malloc(sizeof(char)*mstr,"Allocating name array");
      if (names == NULL)
        exit (1);
    }
  else
    names = NULL;

  mseq = 0;
  for (i = 0; i < nctg; i++)
    if (cseq[i])
      sarray[mseq++] = (uint64) abs(cseq[i]) << 32 | i;

  qsort(sarray, mseq, sizeof(uint64), USORT);

  c1 = (uint32) sarray[0];
  r1 = cseq[c1];
  tseq = 0;
  nseq = 0;
  nstr = 0;
  for (j = 0, i = 1; i < mseq; i++)
    { caxis[c1] = tseq - cbeg[c1];
      tseq += cend[c1] - cbeg[c1];
      c2 = (uint32) sarray[i];
      r2 = cseq[c2];
      if (cend[c1] < clen[c1] ||      // end-clipping
              c1 + 1 < c2 ||          // contigs skipped
              cmap[c1] != cmap[c2] || // scaffolds switched
              r1 != r2 ||             // ranks changed
              cbeg[c2] > 0)           // start-clipping
        { // different rank group or gap
          // add new seq [j,i)
          c0 = (uint32) sarray[j];
          s = cmap[c0]; // get scaffold id
          addSeqName(names,nstr,mstr,c0,c1,s,sdict,smap,cbeg,cend,clen,r1,&snames,&names,&nstr,&mstr);
          saxis[nseq++] = tseq;
          // reverse ctg [j,i)
          if (r1 < 0) axisReverse(sarray, caxis, tseq, j, i, cbeg, cend);
          j = i;
        }
      else
        tseq += coff[c2] - coff[c1] - clen[c1];
      c1 = c2;
      r1 = r2;
    }
  // add new seq [j,i)
  caxis[c1] = tseq - cbeg[c1];
  tseq += cend[c1] - cbeg[c1];
  c0 = (uint32) sarray[j];
  s = cmap[c0]; // get scaffold id
  addSeqName(names,nstr,mstr,c0,c1,s,sdict,smap,cbeg,cend,clen,r1,&snames,&names,&nstr,&mstr);
  saxis[nseq++] = tseq;
  if (r1 < 0) axisReverse(sarray, caxis, tseq, j, i, cbeg, cend);

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
    { fprintf(stderr,"%s: sequence [scf] offset: %-8d %-12lld %s\n",Prog_Name,i,saxis[i],NOLABEL? "" : names);
      if (!NOLABEL) names += strlen(names) + 1;
    }
#endif

  free(sarray);

  return (nseq);
}

#define INT_SIGN(x) (((x)>0)-((x)<0))

static void alnConfig()
{ int64 i;
  int aread, bread, abpos, aepos, bbpos, bepos, l, a, b;
  Segment *segment;

  for (i = 0; i < nSegment; i++)
    { segment = &segments[i];
      if (IS_DELT(segment->flag))
        continue;

      aread = segment->aread;
      bread = segment->bread;

      abpos = segment->abpos;
      aepos = segment->aepos;
      bbpos = segment->bbpos;
      bepos = segment->bepos;

      if (ASEQ[aread] < 0)
        { l = AEND[aread] - ABEG[aread];
          abpos = l - abpos;
          aepos = l - aepos;
        }
      if (BSEQ[bread] < 0)
        { l = BEND[bread] - BBEG[bread];
          bbpos = l - bbpos;
          bepos = l - bepos;
        }

      segment->abpos = abpos;
      segment->aepos = aepos;
      segment->bbpos = bbpos;
      segment->bepos = bepos;

      a = abpos - aepos;
      b = bbpos - bepos;

      if (INT_SIGN(a) == INT_SIGN(b))
        SET_RED(segment->flag);
      else
        SET_BLUE(segment->flag);
    }
}

/*************************************************************************
**************** Cohen–Sutherland line clipping algorithm ****************
**************************************************************************
 * given a rectangle clip window bounded by [xmin, xmax, ymin, ymax]
 * and a line segment spanned by [x1, y1] and [x2, y2]
 * return
 *   -1  if the line segment is completely invisible
 *    0  if the line segment is contained, i.e., no clipping
 *   >0  if the line segment is partially visible
 * for non-negative return values, the new coordinates are computed
**************************************************************************/

#define INSIDE 0 // 0000
#define LEFT   1 // 0001
#define RIGHT  2 // 0010
#define BOTTOM 4 // 0100
#define TOP    8 // 1000

static inline int interCode(double x, double y, double xmin, double xmax, double ymin, double ymax)
{ int code = INSIDE;
  
  if (x < xmin)
    code |= LEFT;
  else if (x > xmax)
    code |= RIGHT;
  if (y < ymin)
    code |= BOTTOM;
  else if (y > ymax)
    code |= TOP;

  return (code);
}

int cohen_sutherland_clip(double x1, double y1, double x2, double y2,
        double xmin, double xmax, double ymin, double ymax,
        double *nx1, double *ny1, double *nx2, double *ny2)
{ int code1, code2, code_out;
  double x, y;
  int inter;

  code1 = interCode(x1, y1, xmin, xmax, ymin, ymax);
  code2 = interCode(x2, y2, xmin, xmax, ymin, ymax);

  inter = 0;
  while (1)
    { if (!(code1 | code2))
        break;
      else if (code1 & code2)
        return (-1);
      else
        { code_out = code1? code1 : code2;
          if (code_out & TOP)
            { x = x1 + (x2 - x1) * (ymax - y1) / (y2 - y1);
              y = ymax;
            }
          else if (code_out & BOTTOM)
            { x = x1 + (x2 - x1) * (ymin - y1) / (y2 - y1);
              y = ymin;
            }
          else if (code_out & RIGHT)
            { y = y1 + (y2 - y1) * (xmax - x1) / (x2 - x1);
              x = xmax;
            }
          else
            { y = y1 + (y2 - y1) * (xmin - x1) / (x2 - x1);
              x = xmin;
            }

          if (code_out == code1)
            { x1 = x;
              y1 = y;
              code1 = interCode(x1, y1, xmin, xmax, ymin, ymax);
            }
          else
            { x2 = x;
              y2 = y;
              code2 = interCode(x2, y2, xmin, xmax, ymin, ymax);
            }
          inter++;
        }
    }

  *nx1 = x1;
  *ny1 = y1;
  *nx2 = x2;
  *ny2 = y2;

  return (inter);
}
/**********************************END************************************/

void *aln_clip(void *args)
{ Packet *parm  = (Packet *) args;
  int64 beg  = parm->beg;
  int64 end  = parm->end;
  int64 i;
  Segment *segment;
  int x, aread, bread;
  double abpos, aepos, bbpos, bepos;

  for (i = beg; i < end; i++)
    { segment = segments + i;
      if (IS_DELT(segment->flag))
        continue;
      aread = segment->aread;
      bread = segment->bread;
      if (!ASEQ[aread] || !BSEQ[bread])
        { SET_DELT(segment->flag);
          continue;
        }
      abpos = segment->abpos;
      aepos = segment->aepos;
      bbpos = segment->bbpos;
      bepos = segment->bepos;

      x = cohen_sutherland_clip(abpos, bbpos, aepos, bepos,
              ABEG[aread], AEND[aread], BBEG[bread], BEND[bread],
              &abpos, &bbpos, &aepos, &bepos);

      if (x < 0)
        SET_DELT(segment->flag);
      else if (x > 0)
        { segment->abpos = (int) (abpos + .499);
          segment->aepos = (int) (aepos + .499);
          segment->bbpos = (int) (bbpos + .499);
          segment->bepos = (int) (bepos + .499);
        }
    }

  return (NULL);
}

void aln_filter()
{ int i, digits;
  int64 nseg;
  double alen;
  uint64 *sarray;
  Segment *s;

  // recalculate alignment coordinates
  { int p;
    Packet *parm;
    pthread_t threads[NTHREADS];

    parm = Malloc(sizeof(Packet)*NTHREADS,"Allocating thread records");
    if (parm == NULL)
      exit (1);

    for (p = 0; p < NTHREADS ; p++)
      { parm[p].beg = (p * nSegment) / NTHREADS;
        if (p > 0)
          parm[p-1].end = parm[p].beg;
      }
    parm[NTHREADS-1].end = nSegment;

    for (p = 1; p < NTHREADS; p++)
      pthread_create(threads+p,NULL,aln_clip,parm+p);
    aln_clip(parm);
    for (p = 1; p < NTHREADS; p++)
      pthread_join(threads[p],NULL);

    free(parm);
  }

  nseg = 0;
  for (i = 0, s = segments; i < nSegment; i++, s++)
    if (!IS_DELT(s->flag))
      nseg++;

  if (MAXALIGN == 0 || nseg <= MAXALIGN) return;

  sarray = (uint64 *)Malloc(sizeof(uint64)*nseg,"Allocating seq array");
  if (sarray == NULL)
    exit (1);
  
  nseg = 0;
  for (i = 0, s = segments; i < nSegment; i++, s++)
    if (!IS_DELT(s->flag))
      sarray[nseg++] = (uint64) (s->aepos-s->abpos) << 32 | i;

  qsort(sarray, nseg, sizeof(uint64), DSORT);

  alen = (double) (sarray[MAXALIGN-1] >> 32);
  digits = 1;
  while (alen > digits) {alen /= 10; digits *= 10;};
  alen = ((int) alen + 1) * digits;

  for (i = 0; i < nseg; i++)
    if ((sarray[i] >> 32) < alen)
        SET_DELT(segments[sarray[i]&0xFFFFFFFFU].flag);

  nseg = 0;
  for (i = 0, s = segments; i < nSegment; i++, s++)
    if (!IS_DELT(s->flag))
      nseg++;

  free(sarray);
  
  fprintf(stderr, "%s: using length filter threshold %.0f\n",Prog_Name,alen);
  fprintf(stderr, "%s: selected %lld alignments to plot\n",Prog_Name,nseg); 
}

#define eps_header(fp,x,y,linewidth) { \
    fprintf(fp,"%%!PS-Adobe-3.0 EPSF-3.0\n"); \
    fprintf(fp,"%%%%BoundingBox:"); \
    fprintf(fp," 1 1 %g %g\n\n",(float)(x),(float)(y)); \
    fprintf(fp,"/C { dup 255 and 255 div exch dup -8 bitshift 255 and 255 div 3 1 roll -16 bitshift 255 and 255 div 3 1 roll setrgbcolor } bind def\n"); \
    fprintf(fp,"/L { 4 2 roll moveto lineto } bind def\n"); \
    fprintf(fp,"/LX { dup 4 -1 roll exch moveto lineto } bind def\n"); \
    fprintf(fp,"/LY { dup 4 -1 roll moveto exch lineto } bind def\n"); \
    fprintf(fp,"/LS { 3 1 roll moveto show } bind def\n"); \
    fprintf(fp,"/MS { dup stringwidth pop 2 div 4 -1 roll exch sub 3 -1 roll moveto show } bind def\n"); \
    fprintf(fp,"/RS { dup stringwidth pop 4 -1 roll exch sub 3 -1 roll moveto show } bind def\n"); \
    fprintf(fp,"/B { 4 copy 3 1 roll exch 6 2 roll 8 -2 roll moveto lineto lineto lineto closepath } bind def\n");\
    fprintf(fp,"%g setlinewidth\n\n",linewidth);\
}
#define eps_font(fp,f,s) do { \
    fprintf(fp,"/FS %d def\n",s); \
    fprintf(fp,"/FS4 FS 4 div def\n"); \
    fprintf(fp,"/%s findfont FS scalefont setfont\n\n",f); \
} while (0)

#define eps_bottom(fp) fprintf(fp,"stroke showpage\n")
#define eps_color(fp,col) fprintf(fp,"stroke %d C\n",col)
#define eps_gray(fp,gray) fprintf(fp, "%g setgray\n",(float)gray)
#define eps_linewidth(fp, lw) fprintf(fp, "%g setlinewidth\n", (float)(lw))
#define eps_line(fp,x1,y1,x2,y2) fprintf(fp,"%g %g %g %g L\n",(float)(x1),(float)(y1),(float)(x2),(float)(y2))
#define eps_linex(fp,x1,x2,y) fprintf(fp,"%g %g %g LX\n",(float)(x1),(float)(x2),(float)(y))
#define eps_liney(fp,y1,y2,x) fprintf(fp,"%g %g %g LY\n",(float)(y1),(float)(y2),(float)(x))
#define eps_Mstr(fp,x,y,s) fprintf(fp,"%g %g (%s) MS\n",(float)(x),(float)(y),s)
#define eps_Mint(fp,x,y,i) fprintf(fp,"%g %g (%d) MS\n",(float)(x),(float)(y),i)
#define eps_stroke(fp) fprintf(fp,"stroke\n")

#define N_COLOR 0xFF0000
#define C_COLOR 0x0080FF

static int SEG_COLOR[2] = {N_COLOR, C_COLOR};

void make_plot(FILE *fo)
{ // generate eps file
  int i, c;
  int width, height, fsize, maxis;
  char *xnames, *ynames;
  int nxseq, nyseq;
  int64 txseq, tyseq, *cxoff, *sxoff, *cyoff, *syoff;
  double sx, sy;
  double lsize;

  // find total length of x- and y-axis and order of plotting
  txseq = tyseq = 0;
  nxseq = axisConfig(BDIC,BSCF,NBCTG,BSEQ,BBEG,BEND,BMAP,BLEN,BOFF,BNAME,&cxoff,&sxoff,&xnames,&txseq);
  nyseq = axisConfig(ADIC,ASCF,NACTG,ASEQ,ABEG,AEND,AMAP,ALEN,AOFF,ANAME,&cyoff,&syoff,&ynames,&tyseq);

  alnConfig();

  width  = IMGWIDTH;
  height = IMGHEIGH;
  if (!height)
    height = (int)((double) width / txseq * tyseq + .499);
  if (!width)
    width = (int)((double) height / tyseq * txseq + .499);
  
  maxis = width > height? width : height;
  if (maxis > MAX_XY_LEN)
    { double scale = (double) MAX_XY_LEN / maxis;
      fprintf(stderr,"%s: image size too large [%d]x[%d]\n",Prog_Name,width,height);
      width  = (int)(width  * scale + 0.499);
      height = (int)(height * scale + 0.499);
      fprintf(stderr,"%*s  shrink the size to [%d]x[%d]\n",(int) strlen(Prog_Name),"",width,height);
      
      if (width < MIN_XY_LEN)
        { fprintf(stderr,"%s: image width too small [%d]\n",Prog_Name,width);
          fprintf(stderr,"%*s  reset image width to [%d]\n",(int) strlen(Prog_Name),"",
                                                            MAX_XY_LEN);
          fprintf(stderr,"%*s  image and sequence size are not in proportion\n",
                         (int) strlen(Prog_Name),"");
          width = MIN_XY_LEN;  
        }
      if (height < MIN_XY_LEN)
        { fprintf(stderr,"%s: image height too small [%d]\n",Prog_Name,height);
          fprintf(stderr,"%*s  reset image height to [%d]\n",(int) strlen(Prog_Name),"",
                                                            MAX_XY_LEN);
          fprintf(stderr,"%*s  image and sequence size are not in proportion\n",
                         (int) strlen(Prog_Name),"");
          height = MIN_XY_LEN;
        }
    }
  
  maxis = width < height? width : height;
  if (maxis < MIN_XY_LEN)
    { double scale = (double) MIN_XY_LEN / maxis;
      fprintf(stderr,"%s: image size too small [%d]x[%d]\n",Prog_Name,width,height);
      width  = (int)(width  * scale + 0.499);
      height = (int)(height * scale + 0.499);
      fprintf(stderr,"%*s  rescale the size to [%d]x[%d]\n",(int) strlen(Prog_Name),"",
                     width,height);

      if (width > MAX_XY_LEN)
        { fprintf(stderr,"%s: image width too large [%d]\n",Prog_Name,width);
          fprintf(stderr,"%*s  reset image width to [%d]\n",(int) strlen(Prog_Name),"",
                                                            MAX_XY_LEN);
          fprintf(stderr,"%*s  image and sequence size are not in proportion\n",
                         (int) strlen(Prog_Name),"");
          width = MAX_XY_LEN;
        }
      if (height > MAX_XY_LEN)
        { fprintf(stderr,"%s: image height too large [%d]\n",Prog_Name,height);
          fprintf(stderr,"%*s  reset image height to [%d]\n",(int) strlen(Prog_Name),"",
                                                            MAX_XY_LEN);
          fprintf(stderr,"%*s  image and sequence size are not in proportion\n",
                         (int) strlen(Prog_Name),"");
          height = MAX_XY_LEN;
        }
    }

  maxis = width < height? width : height;
  
  fsize = FONTSIZE;
  if (!fsize) fsize = maxis / 60;
  
  lsize = LINESIZE;
  if (lsize < 1e-6)
      lsize = (double) maxis / 500;

  sx = (double)  width / txseq;
  sy = (double) height / tyseq;

  eps_header(fo, width, height, lsize);
  eps_font(fo, "Helvetica-Narrow", fsize);
  eps_gray(fo, .8);

  if (!NOLABEL)
    { char *names;

      // write x labels
      names = xnames;
      eps_Mstr(fo, .5 * sxoff[0] * sx, fsize * .5, names);
      names += strlen(names) + 1;
      for (i = 1; i < nxseq; i++)
        { eps_Mstr(fo, .5 * (sxoff[i-1] + sxoff[i]) * sx, fsize * .5, names);
          names += strlen(names) + 1;
        }
      eps_stroke(fo);
      fprintf(fo, "gsave %g 0 translate 90 rotate\n", fsize * 1.25);
      
      // write y labels
      names = ynames;
      eps_Mstr(fo, .5 * syoff[0] * sy, 0, names);
      names += strlen(names) + 1;
      for (i = 1; i < nyseq; i++)
        { eps_Mstr(fo, .5 * (syoff[i-1] + syoff[i]) * sy, 0, names);
          names += strlen(names) + 1;
        }
      fprintf(fo, "grestore\n");
      eps_stroke(fo);
    }

  // write grid lines
  eps_linewidth(fo, lsize/2);
  eps_linex(fo, 1, width,  1);
  for (i = 0; i < nyseq; i++)
    eps_linex(fo, 1, width,  syoff[i] * sy);
  eps_liney(fo, 1, height, 1);
  for (i = 0; i < nxseq; i++)
    eps_liney(fo, 1, height, sxoff[i] * sx);
  eps_stroke(fo);

  // write segments
  int aread, bread;
  double x0, y0, x1, y1, xo, yo;
  Segment *segment;
  eps_linewidth(fo, lsize);
  for (c = 0; c < 2; c++)
    { eps_color(fo, SEG_COLOR[c]);
      for (i = 0; i < nSegment; i++)
        { segment = &segments[i];
          if (IS_DELT(segment->flag) ||
                  !(segment->flag&(1<<(c+1))))
            continue;

          aread = segment->aread;
          bread = segment->bread;

          xo = cxoff[bread];
          yo = cyoff[aread];
          x0 = segment->bbpos;
          x1 = segment->bepos;
          y0 = segment->abpos;
          y1 = segment->aepos;

          x0 = (x0 + xo) * sx;
          x1 = (x1 + xo) * sx;
          y0 = (y0 + yo) * sy;
          y1 = (y1 + yo) * sy;

          eps_line(fo, x0, y0, x1, y1);
        }
      eps_stroke(fo);
    }
  eps_bottom(fo);

  free(cxoff);
  free(cyoff);
  free(sxoff);
  free(syoff);
  free(xnames);
  free(ynames);
}

int main(int argc, char *argv[])
{ //  Process options

  int    i, j, k;
  int    flags[128];
  char  *xseq, *yseq, *pdf;
  char  *eptr;
  FILE  *foeps;
  char  *pdftool;

  ARG_INIT("ALNplot")
    
  xseq = NULL;
  yseq = NULL;
  pdf  = NULL;
  j = 1;
  for (i = 1; i < argc; i++)
    if (argv[i][0] == '-')
      switch (argv[i][1])
        { default:
            ARG_FLAGS("dhSL")
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
          case 'e':
            ARG_REAL(MINAIDNT)
            break;
          case 'a':
            ARG_NON_NEGATIVE(MINALEN,"Minimum alignment length")
            break;
          case 'p':
            if (argv[i][2] == ':' && argv[i][3] != '\0')
              pdf = argv[i]+3;
            else
              pdf = "";
            break;
          case 'x':
            xseq = argv[i]+2;
            break;
          case 'y':
            yseq = argv[i]+2;
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

  if (argc < 2 || argc > 4 || flags['h'])
    { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
      fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
      fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
      fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[3]);
      fprintf(stderr,"\n");
      fprintf(stderr,"     <range> =    <contig>[-<contig>]     |  <contig>_<int>-(<int>|#)\n");
      fprintf(stderr,"             | @[<scaffold>[-<scaffold>]] | @<scaffold>_<int>-(<int>|#)\n");
      fprintf(stderr,"\n");
      fprintf(stderr,"        <contig>   = (<int>|#)[.(<int>|#)]\n");
      fprintf(stderr,"\n");
      fprintf(stderr,"        <scaffold> =  <int>|<string>|#\n");
      fprintf(stderr,"\n");
      fprintf(stderr,"      -a: minimum alignment length\n");
      fprintf(stderr,"      -e: minimum alignment similarity\n");
      fprintf(stderr,"      -x: sequences placed on x-axis\n");
      fprintf(stderr,"      -y: sequences placed on y-axis\n");
      fprintf(stderr,"      -d: try to put alignments along the diagonal line\n");
      fprintf(stderr,"      -S: print sequence IDs as labels instead of names\n");
      fprintf(stderr,"      -L: do not print labels\n");
      fprintf(stderr,"      -H: image height\n");
      fprintf(stderr,"      -W: image width\n");
      fprintf(stderr,"      -f: label font size\n");
      fprintf(stderr,"      -t: line thickness\n");
      fprintf(stderr,"      -n: maximum number of lines to display (set '0' to force all)\n");
      fprintf(stderr,"      -T: use -T threads\n");
      fprintf(stderr,"\n");
      fprintf(stderr,"      -p: make PDF output (requires \'[e]ps[to|2]pdf\')\n");
      fprintf(stderr,"\n");
      if (flags['h'])
        exit (0);
      else
        exit (1);
    }

  TRYADIAG = flags['d'];
  PRINTSID = flags['S'];
  NOLABEL  = flags['L'];

  if (pdf != NULL && !(pdftool = findPDFtool()))
    { fprintf(stderr,"%s: Cannot find [e]ps[to|2]pdf needed to produce .pdf output\n",Prog_Name);
      exit (1);
    }

  if (TRYADIAG)
    { fprintf(stderr,"%s: diagonalisation (-d) is not supported yet\n",Prog_Name);
      exit (1);
    }

  if (IMGWIDTH && IMGHEIGH)
    fprintf(stderr,"%s: setting both image width and height is not recommended\n",Prog_Name);

  if (!IMGWIDTH && !IMGHEIGH) IMGHEIGH = 600;

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

  parseTargetSEQ(yseq, NACTG, NASCF, ADIC, ASCF, AMAP, ALEN, AOFF, &ASEQ, &ABEG, &AEND, &ANAME);
  parseTargetSEQ(xseq, NBCTG, NBSCF, BDIC, BSCF, BMAP, BLEN, BOFF, &BSEQ, &BBEG, &BEND, &BNAME);

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
  
  free(AMAP);
  free(ASEQ);
  free(BSEQ);
  free(ANAME);
  free(BNAME);
  dictDestroy(ADIC);
  if (ISTWO)
    dictDestroy(BDIC);

  free(segments);
  free(OUTEPS);
  free(Prog_Name);

  Catenate(NULL,NULL,NULL,NULL);

  exit (0);
}
