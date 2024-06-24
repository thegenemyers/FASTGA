/*******************************************************************************************
 *
 *  Display a specified set of reads of a database in fasta format.
 *
 *  Author:  Gene Myers
 *  Origin:  Reworking of DAZZ_DB program DBshow.c from the FASTGA module specialized to "gdb" ("dam")
 *  Date  :  Jan 2024
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <fcntl.h>

#include "GDB.h"
#include "hash.h"
#include "select.h"

#define DEBUG_RANGE

static char *Usage =
            "[-hU] [-w<int(80)>] <source:path>[.1gdb] [ <selection>|<FILE> ]";

#define MAX_BUFFER   10001

static int           NCONTIG;  // # of contigs (across all scaffolds)
static int           NSCAFF;   // # of scaffolds
static char         *HEADERS;   // All headers
static GDB_CONTIG   *CONTIGS;   // contig vector of the database
static GDB_SCAFFOLD *SCAFFS;    // scaffold vector of the database
static Hash_Table   *HASH;

int main(int argc, char *argv[])
{ GDB        _gdb, *gdb = &_gdb;
  Selection *select;
  int        llen;

  int         UPPER;
  int         DOSEQ;
  int         WIDTH;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];
    char *eptr;

    ARG_INIT("GDBshow")

    WIDTH = 80;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("hU")
            break;
          case 'w':
            ARG_NON_NEGATIVE(WIDTH,"Line width")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (flags['U'])
      UPPER = UPPER_CASE;
    else
      UPPER = LOWER_CASE;
    DOSEQ = 1-flags['h'];

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"     <selection> = <range> [ , <range> ]*\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"     <range> =    <contig>[-<contig>]     |  <contig>_<int>-(<int>|#)\n");
        fprintf(stderr,"             | @[<scaffold>[-<scaffold>]] | @<scaffold>_<int>-(<int>|#)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"        <contig>   = (<int>|#)[.(<int>|#)]\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"        <scaffold> =  <int>|<string>|#\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -h: Show only the header lines.\n");
        fprintf(stderr,"      -w: Print -w bp per line (default is 80).\n");
        fprintf(stderr,"      -U: Use upper case for DNA (default is lower case).\n");
        exit (1);
      }
  }

  //  Read GDB and establish sorted order of headers

  { int   s, c;
    char *sptr, *eptr;

    Read_GDB(gdb,argv[1]);

    NSCAFF  = gdb->nscaff;
    NCONTIG = gdb->ncontig;
    CONTIGS = gdb->contigs;
    SCAFFS  = gdb->scaffolds;
    HEADERS = gdb->headers;

    HASH  = New_Hash_Table(NSCAFF,1);
    for (s = 0; s < NSCAFF; s++)
      { sptr = HEADERS + SCAFFS[s].hoff;
        for (eptr = sptr; *eptr != '\n'; eptr++)
          if (isspace(*eptr))
            break;
        c = *eptr;
        *eptr = '\0';
        if (Hash_Lookup(HASH,sptr) < 0)
          Hash_Add(HASH,sptr);
        else
          { fprintf(stderr,"%s: Duplicate scaffold name: %s\n",Prog_Name,sptr);
            exit (1);
          }
        *eptr = c;
      }
  }

  //  Determine if extra arg is a range or a file, and if a file then initialize it

  if (argc == 3)
    select = get_selection_list(argv[2],gdb,HASH,&llen);
  else
    select = get_selection_list(NULL,gdb,HASH,&llen);

  //  Display each read (and/or QV streams) in the active DB according to the
  //    range pairs in pts[0..reps) and according to the display options.

  { GDB_CONTIG   *r;
    GDB_SCAFFOLD *s;
    char         *contig;
    int           c, b, e, k;
    int           substr;
    int           len, fst, lst;
    char         *nstring;

    nstring = Malloc(WIDTH+1,"Allocating write buffer\n");
    if (UPPER == 2)
      for (c = 0; c < WIDTH; c++)
        nstring[c] = 'N';
    else
      for (c = 0; c < WIDTH; c++)
        nstring[c] = 'n';
    nstring[WIDTH] = '\0';

    contig = New_Contig_Buffer(gdb);
    substr = 0;

    for (c = 0; c < llen; c++)
      { if (select[c].type % 2)
          { b = e = select[c].src;
            substr = 1;
            fst = select[c].beg;
            lst = select[c].end;
          }
        else
          { b = select[c].beg;
            e = select[c].end;
            substr = 0;
          }
printf("%d %d: %d %d %d %d :: %d %d\n",c,llen,b,e,substr,select[c].type,gdb->nscaff,gdb->ncontig);

        if (select[c].type >= 2)

          { for (k = b; k <= e; k++)
    
              { int cbeg, cend, wpos;
                int j, u, w, f, l;

                r   = CONTIGS + SCAFFS[k].fctg;
                len = SCAFFS[k].slen;
                if (!substr)
                  { fst = 0;
                    lst = len;
                  }

                printf(">%s",HEADERS+SCAFFS[k].hoff);
                if (substr)
                  printf(" :: [%d,%d]",fst,lst);
                printf("\n");

#define ROLLING_PRINT(f,target)			\
  for (j = f; j+(w = WIDTH-wpos) <= l; j += w)	\
    { printf("%.*s\n",w,target);		\
      wpos = 0;					\
    }						\
  if (j < l)					\
    { printf("%.*s",l-j,target);		\
      wpos += l-j;				\
    }

                if (DOSEQ)
                  { cend = 0;
                    wpos = 0;
                    for (u = SCAFFS[k].fctg; u < SCAFFS[k].ectg; u++, r++)
                      { cbeg = r->sbeg;
                        if (cend < fst)
                          cend = fst;
                        if (cend < lst && cbeg > fst)
                          { if (cbeg <= lst)
                              l = cbeg-cend;
                            else
                              l = lst-cend;

                            ROLLING_PRINT(0,nstring)
                          }

                        cend = cbeg + r->clen;
                        if (cbeg < lst && cend > fst)
                          { Get_Contig(gdb,u,UPPER,contig);
                            f = fst-cbeg;
                            l = lst-cbeg;
                            if (f < 0) f = 0; 
                            if (l > r->clen) l = r->clen;

                            ROLLING_PRINT(f,contig+j)
                          }
                      }

                    cbeg = SCAFFS[k].slen;
                    if (cend < fst)
                      cend = fst;
                    if (cend < lst && cbeg > fst)
                      { if (cbeg <= lst)
                          l = cbeg-cend;
                        else
                          l = lst-cend;

                        ROLLING_PRINT(0,nstring)
                      }
                    printf("\n");
                  }
              }
          }

        else

          { for (k = b; k <= e; k++)
              { r   = CONTIGS + k;
                s   = SCAFFS + r->scaf;
                len = r->clen;
                if (!substr)
                  { fst = 0;
                    lst = len;
                  }

                if (len > 0)
                  printf(">%s :: Contig %d[%lld,%lld]\n",
                         HEADERS+s->hoff,k-(s->fctg)+1,r->sbeg+fst,r->sbeg+lst);

                if (DOSEQ)
                  { int j;
    
                    Get_Contig(gdb,k,UPPER,contig);
                    for (j = fst; j+WIDTH < lst; j += WIDTH)
                      printf("%.*s\n",WIDTH,contig+j);
                    if (j < lst)
                      printf("%.*s\n",lst-j,contig+j);
                  }
                fflush(stdout);
              }
         }
      }
  }

  Close_GDB(gdb);

  exit (0);
}
