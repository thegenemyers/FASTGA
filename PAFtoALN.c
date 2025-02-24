/*******************************************************************************************
 *
 *  Utility for changing the database names in the header of a .1aln file
 *
 *  Author:    Gene Myers in a version for .las files
 *             Modified by Richard Durbin to use .1aln files
 *  Creation:  July 2013
 *  Last Mod:  March 2024
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
#include "GDB.h"
#include "hash.h"
#include "alncode.h"

#undef  DEBUG_THREADS
#undef  DEBUG

#define TSPACE   100
#define VERSION "0.1"

static char *Usage[] = 
            { "[-T<int(8)>] <alignments:path>[.paf]",
              " <source1:path>[.1gdb|<fa_extn>|<1_extn>] [<source2:path>[.1gdb|<fa_extn>|<1_extn>]]"
            };

static GDB_CONTIG *CONTIG1;
static GDB_CONTIG *CONTIG2;

static GDB_SCAFFOLD *SCAFF1;
static GDB_SCAFFOLD *SCAFF2;

static Hash_Table *AHASH;
static Hash_Table *BHASH;

#undef TEST

// HhPp -> 0, Ii -> 1, DdNn -> 2, = -> 3, Xx -> 4, Mm -> 5

static int interp[128] =
  { -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,

    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1,  3, -1, -1,

    -1, -1, -1, -1,  2, -1, -1, -1,
     0,  1, -1, -1, -1,  5,  2, -1,
     0, -1, -1,  1, -1, -1, -1, -1,
     4, -1, -1, -1, -1, -1, -1, -1,

    -1, -1, -1, -1,  2, -1, -1, -1,
     0,  1, -1, -1, -1,  5,  2, -1,
     0, -1, -1,  1, -1, -1, -1, -1,
     4, -1, -1, -1, -1, -1, -1, -1,
  };

  //  Check that cigar string validly covers span of path and return size of trace point array
  //    needed.  If allow_M is 0, then if an 'M' occurs return 0 (error).

static int cigarCheck(Path *path, char *cigar, int allow_M, int comp)
{ int   len;
  int   apos, bpos;
  char *c;

  apos = path->abpos;
  bpos = path->bbpos;
  for (c = cigar; *c != '\0'; c++)
    { len = 0;
      while (isdigit(*c))
        len = 10*len + (*c++-'0');
      if (len == 0)
        len = 1;
      switch (interp[(int) (*c)])
      { case 5:
          if (!allow_M)
            return (0);
        case 4:
        case 3:
          apos += len;
          bpos += len;
          break;
        case 2:
          bpos += len;
          break;
        case 1:
          apos += len;
          break;
        case 0:
          break;
        default:
          fprintf(stderr,"%s: Invalid Cigar symbol %c(%d)\n",Prog_Name,*c,*c);
          exit (1);
      }
    }

  if (path->aepos != apos || path->bepos != bpos)
    { fprintf(stderr,"%s: Cigar span and alignment intervals do not match\n",Prog_Name);
      exit (1);
    }

  if (comp)               //  Reverse the CIGAR string if complemented
    { char *d, *e, x;

      for (d = cigar, c -= 1; d < c; d++, c--) 
        { x = *d; *d = *c; *c = x; }
      for (c = cigar; *c != '\0'; c = e)
        { d = c+1;
          while (isdigit(*d))
            d += 1; 
          e = d;
          for (d--; c < d; c++, d--)
            { x = *d; *d = *c; *c = x; }
        }
    }

  return (((path->aepos-1)/TSPACE - (path->abpos/TSPACE)) + 1);
}

typedef struct
  { int   apos;   //  At point (apos,bpos) in the d.p. matrix (of current contig pair)
    int   bpos;
    char *cptr;   //  if len > 0, then next command is len,*cptr
    int   len;    //      otherwise it is the length,command pair at cptr.
  } Cigar_Position;

  //  if apos or bpos < 0 then off the matrix, if next command is an indel then also
  //     want to skip as these are part of the gap between alignments.
  //  cigarPrefix moves the cursor C forward until both apos and bpos > 0 and the
  //     next command is a diagonal one.

static void cigarPrefix(Cigar_Position *C)
{ int   len;
  int   apos, bpos;
  char *c;

  apos = C->apos;
  bpos = C->bpos;
  len  = C->len;
  for (c = C->cptr; *c != '\0'; c++)
    { if (len <= 0)
        { while (isdigit(*c))
            len = 10*len + (*c++-'0');
          if (len == 0)
            len = 1;
        }
      switch (interp[(int) (*c)])
      { case 5:
        case 4:
        case 3:
          if (apos >= 0 && bpos > 0)
            goto found;
          if (apos < 0 && apos + len >= 0)
            { len  += apos;
              bpos -= apos;
              apos = 0;
              if (bpos >= 0)
                goto found;
            }
          if (bpos < 0 && bpos + len >= 0)
            { len  += bpos;
              apos -= bpos;
              bpos = 0;
              if (apos >= 0)
                goto found;
            }
          apos += len;
        case 2:
          bpos += len;
          break;
        case 1:
          apos += len;
        case 0:
          break;
      }
      len = 0;
    }
found:
  C->cptr = c;
  C->len  = len;
  C->apos = apos;
  C->bpos = bpos;
}

typedef struct
  { int    diff;
    int    tlen;   // length of current result
    uint8 *trace;
    int    mlen;   // current length of trace array (for realloc purposes)
  } TP_Bundle;     //   set this to 0 for the first call, and free trace after the last call

//  Gvien abeg: start map position in target sequence
//        cigar: 0-terminated CIGAR string
//        tspace: trace point spacing
//        bundle: pointer to a TP "bundle"
//  cigar2tp fills in the TP record for the alignment, but stops short if an insert or delete
//     causes a byte overflow.  It returns a pointer into the cigar string of the last command
//     it interpreted, '\0' if it got to the end.

static char *cigar2tp(Cigar_Position *C, int64 aend, int64 bend,
                      int tspace, TP_Bundle *bundle)
{ int64 apos, anext;
  int64 bpos, blast;
  int64 diff, dlast;
  uint8 *trace;
  int   x, len, inc, slen;
  char *c;

  //  In a pass fill in the tps & dfs vectors, and determine
  //    the end point of the alignment in the target and the length of the read.

  diff  = dlast = 0;
  bpos  = blast = C->bpos;
  apos  = C->apos;
  anext = (apos/tspace+1)*tspace;
  trace = bundle->trace;

  slen = 0;
  len  = C->len;
  for (c = C->cptr; *c != '\0'; c++)
    { if (len <= 0)
        { len = 0;
          while (isdigit(*c))
            len = 10*len + (*c++-'0');
          if (len == 0)
            len = 1;
        }
      if (apos >= aend || bpos >= bend)
        { slen = len;
          break;
        }
      x = interp[(int) (*c)];
      if ((x >= 3 || x == 1) && apos+len > aend)
        { slen = (apos+len)-aend;
          len = aend-apos;
        }
      if (x >= 2 && bpos+len > bend)
        { slen = (bpos+len+slen)-bend;
          len = bend-bpos;
        }
      switch (x)
      { case 4:
          while (apos+len > anext)
            { inc = anext-apos; 
              apos += inc;
              bpos += inc;
              diff += inc;
              len  -= inc;
              anext += tspace;
              *trace++ = diff-dlast;
              *trace++ = bpos-blast;
              blast = bpos;
              dlast = diff;
            }
          apos += len;
          bpos += len;
          diff += len;
          break;
        case 3:
          while (apos+len > anext)
            { inc = anext-apos; 
              apos += inc;
              bpos += inc;
              len  -= inc;
              anext += tspace;
              *trace++ = diff-dlast;
              *trace++ = bpos-blast;
              blast = bpos;
              dlast = diff;
            }
          apos += len;
          bpos += len;
          break;
        case 2:
          if ((bpos-blast) + len + (anext-apos) > 200)
            { slen += len;
              break;
            }
          bpos += len;
          diff += len;
          break;
        case 1:
          if (TSPACE + len > 200)
            { slen += len;
              break;
            }
          while (apos+len > anext)
            { inc = anext-apos; 
              apos += inc;
              diff += inc;
              len  -= inc;
              anext += tspace;
              *trace++ = diff-dlast;
              *trace++ = bpos-blast;
              blast = bpos;
              dlast = diff;
            }
          apos += len;
          diff += len;
        case 0:
          break;
      }
      if (slen > 0)
        break;
      len = 0;
    }
  if (apos > anext-tspace)
    { *trace++ = diff-dlast;
      *trace++ = bpos-blast;
    }
  bundle->diff = diff;
  bundle->tlen = trace - bundle->trace;

  C->apos = apos;
  C->bpos = bpos;
  C->cptr = c;
  C->len  = slen;

  return (c);
}

#ifdef TEST

static char *cigars[5] = {
   "23=1X4=1X10=1X8=2I2=1I32=1X3=1X1=1I1X9=1X9=1X12=1X1=2X2=1X5=1I5=1X2=1X2=1I5=1D3=1X6=1I2="
   "1D2=1X4=1X1=1X16=2X1=1X3=1X19=1X10=1X12=1X14=1X16=1I1=1D8=1X11=1X10=2X2=1X3=1I3=2X2=1X2="
   "4D1=1X1=1D4=8D2=1X17=2X5=1I1=1D8=1X9=1X31=1X1=1X7=1X3=1X1=1I1=1D7=1X9=1X4=1X9=1X2=1X2=2X"
   "12=3X13=2X1=2X3=1X5=1X1=2X4=1X25=1X1=1X11=1I1=1D2=1X2=1I1=1X1=2X1=1X2=1X1=2X13=1I1=1D2=1I"
   "1=1D11=1X1=1X4=1X1=1X1=2X7=1X9=1X24I4=1I8=1X1=1X1=2X5=1D3=1X2=1X4=1X3=2X1=1D4=1I15=2X4=2X"
   "16=2X1=1X4=1X9=1I2=1D2=1X1=1X13=1X6=1X2=1X2=1X3=1X3=1X3=1X2=1X7=1X1=1X3=1X2=1I4=2X3=1X1=2X"
   "2=1X5=2X3=1X3=1X3=1X1=1X4=10D5=2D11=1X9=1D13=1X4=1X1=5I4=1X12=2X2=1X2=3X1=1X1=1X16=1X1=1X"
   "3=1X6=2X2=1X2=1X1D1=1X3=1I1=1X2=1I1X2=1D1=2X3=1D2=2I8=1D2=1I3=1D3=1X1D1=1X2=1X1=7I1=1X4=2X"
   "6=1I1=1D3=1X3=2X2=1X2=1X2=1X2=1X3=1X2=1X6=2X1=1I4=1X7=1D1=1X1=1I8=1X2=1X6=1X2=1D3=2D3=1X5="
   "1X1=1X2=1X2=1X2=1X4=1X11=1X28=1X5=1X1=2X1=1X2=1I1X3=2X2=1X1=2D1=1X5=1X1=3D3=1X8=1X6=11D5="
   "2X14=1X9=1X12=1I4=1D3=",

   "3=1I4=1I3=1X19=2X4=1X2=1X3=1X3=1X4=1X2=1X6=1X5=1X2=1X10=1X10=1X12=1X11=1X5=1I1=1X1D"
   "2=1D1=1D2=1D1=1X1D2=1X1=2X1=1I1=1I1=1D3=1I1X3=1X1D2=1X1=1X2=1D19=2X1=1X5=1X2=1X2=1X4="
   "1X14=1X5=8D2=1X3D35=1X15=1X16=1X1=1X1=1X1=1X2=1D12=1X3=1X7D21=1X11=1X4=1X7=1X2=1X2D"
   "1=2X1=1X3=1X4=2X3D1=1X6=3X1=1I3=2I3=5I1=1X1=6I2=2X1=1D3=2I1=1I2=1X1=3I2=3I24=2X10=1X"
   "2=1X2I1=1X1=1D1=1X1=1X2=1D23=1X4=1X10=1X12=1X19=1X4=1X2=1X12=1X12=1X7=1X6=",

   "5=1X18=1X33=1I3=1I6=1I1=1X13=1X27=1X7=1X3=1X2=1X8=1X3=1X7=1X4=1I1=1D2=3X6=1X7=1X18="
   "1X23=1I29=1X2=1X2=1X1=1I1=1D10=1X6=1X2=1X5=1X6=1I2=1D9=1X23=1X4=1X2=1X11=1X23=1X31="
   "1I1=1D8=1X7=1X3=1X7=",

   "5=1X5=1X5=1X3=1X7=1X11=1X32=7D2=7I3=1X13=1X1=1X1=1X8=1X1=1X2=1X5=1I1=1D3=1D2=1I9=1X"
   "8=18D2=15D2=9D11=1X41=1D2=3D3=1X1=1X2=1X1D10=1X4=1X2D17=1X2=1X40=1X1=1X10=1X2=1X8=1D"
   "2=1X4=1X1=1X4=1X4=1X6=1X2=1X9=",

   "3=1I4=1I3=1X19=2X4=1X2=1X3=1X3=1X4=1X2=1X6=1X5=1X2=1X10=1X10=1X12=1X11=1X5=1I1=1X1D"
   "2=1D1=1D2=1D1=1X1D2=1X1=2X1=1I1=1I1=1D3=1I1X3=1X1D2=1X1=1X2=1D19=2X1=1X5=1X2=1X2=1X4="
   "1X14=1X5=8D2=1X3D35=1X15=1X16=1X1=1X1=1X1=1X2=1D12=1X3=1X7D21=1X11=1X4=1X7=1X2=1X2D"
   "1=2X1=1X3=1X4=2X3D1=1X6=3X1=1I3=2I3=5I1=1X1=6I2=2X1=1D3=2I1=1I2=1X1=3I2=3I24=2X10=1X"
   "2=1X2I1=1X1=1D1=1X1=1X2=1D23=1X4=1X10=1X12=1X19=1X4=1X2=1X12=1X12=1X7=1X6=150D"
   "5=1X18=1X33=1I3=1I6=1I1=1X13=1X27=1X7=1X3=1X2=1X8=1X3=1X7=1X4=1I1=1D2=3X6=1X7=1X18="
   "1X23=1I29=1X2=1X2=1X1=1I1=1D10=1X6=1X2=1X5=1X6=1I2=1D9=1X23=1X4=1X2=1X11=1X23=1X31="
   "1I1=1D8=1X7=1X3=1X7=150I"
   "5=1X5=1X5=1X3=1X7=1X11=1X32=7D2=7I3=1X13=1X1=1X1=1X8=1X1=1X2=1X5=1I1=1D3=1D2=1I9=1X"
   "8=18D2=15D2=9D11=1X41=1D2=3D3=1X1=1X2=1X1D10=1X4=1X2D17=1X2=1X40=1X1=1X10=1X2=1X8=1D"
   "2=1X4=1X1=1X4=1X4=1X6=1X2=1X9="
};

static int64 starts[5] = { 126, 31111, 0, 31726, 31111 };

int main(int argc, char *argv[])
{ TP_Bundle answer;
  int       i, j;
  char     *c;

  (void) argc;
  (void) argv;

  answer.mlen = 0;

  for (j = 0; j < 5; j++)
    { c = cigars[j];
      while (1)
        { c = cigar2tp(starts[j],c,100,&answer);

          printf("T ");
          for (i = 1; i < answer.tlen; i += 2)
            printf(" %d",answer.trace[i]);
          printf("\n");

          printf("D ");
          for (i = 0; i < answer.tlen; i += 2)
            printf(" %d",answer.trace[i]);
          printf("\n");

          if (*c == '\0')
            break;

          printf("%c %d\n",*c++,answer.gap);
        }
    }

  exit(0);
}

#endif

/*******************************************************************************************
 *
 *  Read .paf file creating array of segments and two GDB skeletons, AGDB & BGDB
 *    and hash tables, AHASH & BHASH, of the scaffold names for each.
 *
 *******************************************************************************************/

typedef struct
  { int64    beg;
    int64    end;
    OneFile *of;
    char    *iname;
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
          exit (1);
        }
    }

  if (fgets(buffer,bmax,input) == NULL)
    { if (feof(input))
        { free(buffer);
          bmax = 0;
          return (NULL);
        }
      fprintf(stderr,"%s: Could not read next line of file %s (offset %lld)\n",
                     Prog_Name,name,(int64) ftello(input));
      exit (1);
    }

  len = strlen(buffer);
  while (buffer[len-1] != '\n')
    { bmax = ((int) (1.4*bmax)) + 100;
      buffer = (char *) realloc(buffer,bmax);
      if (buffer == NULL)
        { fprintf(stderr,"%s: Out of memory reading %s\n",Prog_Name,name);
          exit (1);
        }
      if (fgets(buffer+len,bmax-len,input) == NULL)
        { if (feof(input))
            fprintf(stderr,"%s: Last line of file %s does not end with new-line\n",
                           Prog_Name,name);
          else
            fprintf(stderr,"%s: Could not read next line of file %s (offset %lld)\n",
                           Prog_Name,name,(int64) ftello(input));
          exit (1);
        }
      len += strlen(buffer+len);
    }
  buffer[--len] = '\0';

  lb->bmax   = bmax;
  lb->buffer = buffer;

  return (buffer);
}

void *gen_1aln(void *args)
{ Packet  *parm  = (Packet *) args;
  char    *iname = parm->iname;
  OneFile *of    = parm->of;   

  Line_Bundle    _lb, *lb = &_lb;
  Overlap        _ovl, *ovl = &_ovl;
  TP_Bundle      _tps, *tps = &_tps;
  Cigar_Position _C, *C = &_C;

  int64  span, amnt;
  FILE  *input;
  int    index, nfields, linelen;
  char  *fptrs[100], *eptr;

#ifdef DEBUG
  char *d;
#endif

  int64  alen, blen, tlen;
  int64 *trace64;
  int    comp;
  int    adel, bdel;
  int    bpos, epos;
  int    aend, bend;
  int    i;

  input = fopen(iname,"r");
  if (input == NULL)
    { fprintf(stderr,"%s: Could not open %s for reading\n",Prog_Name,iname);
      exit (1);
    }

  fseeko(input,parm->beg,SEEK_SET);
  span = parm->end-parm->beg;

  lb->bmax = 0;
  if (parm->beg > 0)
    amnt = strlen(read_line(input,iname,lb)) + 1;
  else
    amnt = 0;

  tps->mlen  = 0;
  tps->trace = NULL;
  trace64    = NULL;
 
  while (amnt <= span)
    { eptr = read_line(input,iname,lb);
      if (eptr == NULL)
        break;
      linelen = strlen(eptr) + 1;
      amnt   += linelen;
    
      while (*eptr != '\0' && isspace(*eptr))    // parse PAF line
        eptr += 1;
      nfields = 0;
      while (*eptr != '\0')
        { fptrs[nfields++] = eptr;
          while (*eptr != '\0' && !isspace(*eptr))
            eptr += 1;
          *eptr++ = '\0';
          while (isspace(*eptr))
            eptr += 1;
        }
 
      if (nfields < 11)
        { fprintf(stderr,"%s: Line of paf has fewer than 11 fields (offset %lld)\n",
                         Prog_Name,(int64) (ftello(input)-linelen));
          exit (1);
        }

      index = Hash_Lookup(AHASH, fptrs[0]);
      alen  = strtol(fptrs[1], &eptr, 10);
      if (index < 0 || alen != SCAFF1[index].slen)
        { fprintf(stderr,"%s: Scaffold name %s not found in first source (offset %lld)\n", 
                         Prog_Name,fptrs[0],(int64) (ftello(input)-linelen));
          exit (1);
        }
      bpos = strtol(fptrs[2], &eptr, 10);
      epos = strtol(fptrs[3], &eptr, 10);
      for (i = SCAFF1[index].fctg; i < SCAFF1[index].ectg; i++)
        if (bpos < CONTIG1[i].sbeg + CONTIG1[i].clen)
          break;
      ovl->aread = i;
      ovl->path.abpos = bpos - CONTIG1[i].sbeg;
      ovl->path.aepos = epos - CONTIG1[i].sbeg;

      index = Hash_Lookup(BHASH, fptrs[5]);
      blen   = strtol(fptrs[6], &eptr, 10);
      if (index < 0 || blen != SCAFF2[index].slen)
        { fprintf(stderr,"%s: Scaffold name %s not found in second source (offset %lld)\n", 
                         Prog_Name,fptrs[5],(int64) (ftello(input)-linelen));
          exit (1);
        }
      bpos = strtol(fptrs[7], &eptr, 10);
      epos = strtol(fptrs[8], &eptr, 10);
      comp = (*fptrs[4] == '-');
      if (comp)
        { for (i = SCAFF2[index].ectg-1; i >= SCAFF2[index].fctg; i--)
            if (epos > CONTIG2[i].sbeg)
              break;
          ovl->bread = i;
          ovl->path.bbpos = (CONTIG2[i].sbeg + CONTIG2[i].clen) - epos;
          ovl->path.bepos = (CONTIG2[i].sbeg + CONTIG2[i].clen) - bpos;
          ovl->flags = COMP_FLAG;
        }
      else
        { for (i = SCAFF2[index].fctg; i < SCAFF2[index].ectg; i++)
            if (bpos < CONTIG2[i].sbeg + CONTIG2[i].clen)
              break;
          ovl->bread = i;
          ovl->path.bbpos = bpos - CONTIG2[i].sbeg;
          ovl->path.bepos = epos - CONTIG2[i].sbeg;
          ovl->flags = 0;
        }

      for (i = 11; i < nfields; i++)
        if (strncmp(fptrs[i],"cg:Z:",5) == 0)
          break;
      if (i >= nfields)
        { fprintf(stderr,"%s: PAF line is missing a CIGAR string (offset %lld)\n",
                         Prog_Name,(int64) (ftello(input)-linelen));
          exit (1);
        }

      tlen = cigarCheck(&(ovl->path),fptrs[i]+5,0,comp);

      if (tlen == 0)
        { fprintf(stderr,"%s: PAF CIGAR string uses M, should be X & = (offset %lld)\n",
                         Prog_Name,(int64) (ftello(input)-linelen));
          exit (1);
        }

      if (tlen >= tps->mlen)
        { tps->mlen  = 1.2*tlen + 250;
          tps->trace = realloc(tps->trace,2*tps->mlen);
          trace64    = realloc(trace64,sizeof(int64)*tps->mlen);
          if (tps->trace == NULL || trace64 == NULL)
            { fprintf(stderr,"Out of memory interpreting CIGAR string");
              exit (1);
            }
        }

      bend = CONTIG2[ovl->bread].clen;
      aend = CONTIG1[ovl->aread].clen;

      oneWriteLine(of,'a',0,0);

      C->len = 0;
      C->apos = ovl->path.abpos;
      C->bpos = ovl->path.bbpos;
      C->cptr = fptrs[i]+5;

      cigarPrefix(C);   //  If starting in a gap skip it.

#ifdef DEBUG
      printf("Initial skip: %d %d",C->apos-ovl->path.abpos,C->bpos-ovl->path.bbpos);
      printf(" %d %.7s\n",C->len,C->cptr);
      fflush(stdout);
#endif

      while (1)
        { ovl->path.abpos = C->apos;
          ovl->path.bbpos = C->bpos;

#ifdef DEBUG
          d = C->cptr;
#endif

          cigar2tp(C,aend,bend,TSPACE,tps);

#ifdef DEBUG
          printf("Piece: %.*s  %d  %.7s\n",(int) (C->cptr-d),d,C->len,C->cptr);
          fflush(stdout);
#endif

          ovl->path.bepos = C->bpos;
          ovl->path.aepos = C->apos;
          ovl->path.diffs = tps->diff;
          Write_Aln_Overlap(of,ovl);
          Write_Aln_Trace(of,tps->trace,tps->tlen,trace64);

          if (*C->cptr == '\0')
            break;

          adel = bdel = 0;
          if (interp[(int) (*C->cptr)] == 1)
            { adel += C->len;
              C->apos += C->len;
              C->cptr += 1;
              C->len   = 0;
            }
          else if (interp[(int) (*C->cptr)] == 2)
            { bdel += C->len;
              C->bpos += C->len;
              C->cptr += 1;
              C->len   = 0;
            }
          while (C->apos >= aend)
            { C->apos += CONTIG1[ovl->aread].sbeg;
              ovl->aread += 1;
              C->apos -= CONTIG1[ovl->aread].sbeg;
              aend = CONTIG1[ovl->aread].clen;
            }
          while (C->bpos >= bend)
            { if (comp)
                { C->bpos -= CONTIG2[ovl->bread].sbeg+CONTIG2[ovl->bread].clen;
                  ovl->bread -= 1;
                  C->bpos += CONTIG2[ovl->bread].sbeg+CONTIG2[ovl->bread].clen;
                  bend = CONTIG1[ovl->bread].clen;
                }
              else
                { C->bpos += CONTIG2[ovl->bread].sbeg;
                  ovl->bread += 1;
                  C->bpos -= CONTIG2[ovl->bread].sbeg;
                  bend = CONTIG1[ovl->bread].clen;
                }
            }

          adel -= C->apos;
          bdel -= C->bpos;
          cigarPrefix(C);
          adel += C->apos;
          bdel += C->bpos;

#ifdef DEBUG
          printf("Internal skip: %d %d",adel,bdel);
          printf(" %d %.7s\n",C->len,C->cptr);
          fflush(stdout);
#endif

          if (adel+bdel > 0)
            { oneInt(of,0) = adel;
              oneInt(of,1) = bdel;
              oneWriteLine(of,'p',0,0);
            }
        }

#ifdef DEBUG
      printf("Done\n");
      fflush(stdout);
#endif
    }

  free(trace64);
  free(tps->trace);

  fclose(input);

  return (NULL);
}


int main(int argc, char *argv[])
{ Packet    *parm;
  GDB       _gdb1, *gdb1 = &_gdb1;
  GDB       _gdb2, *gdb2 = &_gdb2;
  int        ISTWO;
  char      *SPATH1;
  char      *SPATH2;
  char      *APATH, *AROOT;

  int   NTHREADS;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char*  eptr;

    ARG_INIT("PAFtoALN")

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
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc < 3 || argc > 4)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        exit (1);
      }

    parm = Malloc(sizeof(Packet)*NTHREADS,"Allocating thread records");
    if (parm == NULL)
      exit (1);
  }

  //  Find sources

  { char *tpath;
    int   type;

    ISTWO = 0;
    type  = Get_GDB_Paths(argv[2],NULL,&SPATH1,&tpath,0);
    if (type != IS_GDB)
      Create_GDB(gdb1,SPATH1,type,0,NULL,0);
    else
      { Read_GDB(gdb1,tpath);
        if (gdb1->seqs == NULL)
          { fprintf(stderr,"%s: GDB %s must have sequence data\n",Prog_Name,tpath);
            exit (1);
          }
      }
    free(tpath);

    if (argc == 4)
      { type = Get_GDB_Paths(argv[3],NULL,&SPATH2,&tpath,0);
        if (type != IS_GDB)
          Create_GDB(gdb2,SPATH2,type,0,NULL,0);
        else
          { Read_GDB(gdb2,tpath);
            if (gdb2->seqs == NULL)
              { fprintf(stderr,"%s: GDB %s must have sequence data\n",Prog_Name,tpath);
                exit (1);
              }
          }
        free(tpath);
        ISTWO = 1;
      }
    else
      gdb2 = gdb1;
  }

  //  Set up scaffold->contig maps & scaffold name dictionary

  { int   s;
    char *head, *sptr, *eptr;

    CONTIG1 = gdb1->contigs;
    CONTIG2 = gdb2->contigs;
    SCAFF1  = gdb1->scaffolds;
    SCAFF2  = gdb2->scaffolds;

    AHASH = New_Hash_Table(gdb1->nscaff,0);
    head  = gdb1->headers;
    for (s = 0; s < gdb1->nscaff; s++)
      { sptr = head + SCAFF1[s].hoff;
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
      { BHASH = New_Hash_Table(gdb2->nscaff,0);
        head  = gdb2->headers;
        for (s = 0; s < gdb2->nscaff; s++)
          { sptr = head + SCAFF2[s].hoff;
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

  //  Open .1aln file reading and read header information
  
  { char       *inName, *cpath;
    OneFile    *of;
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

    cpath = getcwd(NULL,0);

    of = open_Aln_Write(Catenate(APATH,"/",AROOT,".1aln"),NTHREADS,Prog_Name,
                   VERSION,Command_Line,TSPACE,SPATH1,SPATH2,cpath);

    free(cpath);
    free(SPATH1);
    if (ISTWO)
      free(SPATH2);

    for (p = 0; p < NTHREADS; p++)
      { parm[p].beg   = (state.st_size * p) / NTHREADS;
        parm[p].end   = (state.st_size * (p+1)) / NTHREADS;
        parm[p].iname = inName;
        parm[p].of    = of+p;
      }

    //  Launch and then gather threads

#ifdef DEBUG_THREADS
    for (p = 0; p < NTHREADS; p++)
      gen_1aln(parm+p);
#else
    for (p = 1; p < NTHREADS; p++)
      pthread_create(threads+p,NULL,gen_1aln,parm+p);
    gen_1aln(parm);
    for (p = 1; p < NTHREADS; p++)
      pthread_join(threads[p],NULL);
#endif

    oneFileClose(of);
    free(inName);
    free(AROOT);
    free(APATH);
  }

  if (argc == 4)
    Close_GDB(gdb2);
  Close_GDB(gdb1);
 
  exit (0);
}
