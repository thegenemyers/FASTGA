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
#include "ANO.h"
#include "hash.h"
#include "select.h"

#define DEBUG_RANGE

static char *Usage = " <source:path>[.1ano]] [ <selection>|<FILE> ]";

#define MAX_BUFFER   10001

static int           NSCAFF;    // # of scaffolds
static char         *HEADERS;   // All headers
static GDB_CONTIG   *CONTIGS;   // contig vector of the database
static GDB_SCAFFOLD *SCAFFS;    // scaffold vector of the database
static Hash_Table   *HASH;      // Hash table of all scaffold names

#define SOEL '<'    //  Start/End of element (scaffold or contig)
#define EOEL '>'
#define SPOS '['    //  Start/ending pposition (but not of element)
#define EPOS ']'

static int64    *MOFF;
static ANO_PAIR *MASK;

static void reverse_print(int n, int fst, int lst, int64 off)
{ int m;

  m = MOFF[n+1]-1;
  while (m >= MOFF[n] && MASK[m].beg >= lst)
    m -= 1;

  for ( ; m >= MOFF[n]; m--)
    { if (MASK[m].end <= fst)
        continue;
      if (MASK[m].end <= lst)
        printf("[%10lld",MASK[m].end + off);
      else
        printf("<%10lld",lst + off);
      if (MASK[m].beg >= fst)
        printf(" - %10lld]",MASK[m].beg + off);
      else
        printf(" - %10lld>",fst + off);
      if (MASK[m].label != NULL)
        printf(" %s",MASK[m].label);
      if (MASK[m].score > 0)
        printf(" score = %d",MASK[m].score);
      printf("\n");
    }
}

static void forward_print(int n, int fst, int lst, int64 off)
{ int m;

  m = MOFF[n];
  while (m < MOFF[n+1] && MASK[m].end <= fst)
    m += 1;

  for ( ; m < MOFF[n+1] && MASK[m].beg < lst; m++)
    { if (MASK[m].end <= fst)
        continue;
      if (MASK[m].beg >= fst)
        printf("[%10lld",MASK[m].beg + off);
      else
        printf("<%10lld",fst + off);
      if (MASK[m].end <= lst)
        printf(" - %10lld]",MASK[m].end + off);
      else
        printf(" - %10lld>",lst + off);
      if (MASK[m].label != NULL)
        printf(" %s",MASK[m].label);
      if (MASK[m].score > 0)
        printf(" score = %d",MASK[m].score);
      printf("\n");
    }
}

int main(int argc, char *argv[])
{ ANO        _ano, *ano = &_ano;
  GDB       *gdb;
  Selection *select;
  int        llen;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];

    ARG_INIT("ANOshow")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"  <selection> = <range>[+-] [ , <range>[+-] ]*\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"     <range> = <object/position> [ - <object/position> ]  | @ | .\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"        <object/position> = @ <scaffold> [ . <contig>] [ : <position> ]\n");
        fprintf(stderr,"                          |                . <contig>  [ : <position> ]\n");
        fprintf(stderr,"                          |                                <position>\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"           <scaffold> = # | <int> | <identifier>\n");
        fprintf(stderr,"           <contig>   = # | <int>\n");
        fprintf(stderr,"           <position> = # | <int> [ . <int> ] [kMG]\n");
        exit (1);
      }
  }

  //  Read ANO and establish sorted order of headers

  { int   s, c;
    char *sptr, *eptr;

    Read_ANO(ano,argv[1],NULL);
    gdb = ano->gdb;

    NSCAFF  = gdb->nscaff;
    CONTIGS = gdb->contigs;
    SCAFFS  = gdb->scaffolds;
    HEADERS = gdb->headers;

    MOFF    = ano->moff;
    MASK    = ano->masks;

    HASH  = New_Hash_Table(NSCAFF,1);
    for (s = 0; s < NSCAFF; s++)
      { sptr = HEADERS + SCAFFS[s].hoff;
        for (eptr = sptr; *eptr != '\0'; eptr++)
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

  //  Get the selection list (every contig by default)

  if (argc == 3)
    select = get_selection_list(argv[2],gdb,HASH,&llen);
  else
    select = get_selection_list(NULL,gdb,HASH,&llen);

  //  Display each interval in the ANO according to the 
  //     range list of the selection, observing the display option

  { int           c, k, n;
    int           b, e;
    int           fst, lst, ori;

    for (c = 0; c < llen; c++)
      { ori = select[c].orient;

        if (select[c].type == SCAFF_SELECTION)

          { for (k = select[c].s1; k <= select[c].s2; k++)

              { b = select[c].c1;
                e = select[c].c2;
                fst = CONTIGS[b].sbeg + select[c].p1;
                lst = CONTIGS[e].sbeg + select[c].p2;
                if (k > select[c].s1)
                  { b = SCAFFS[k].fctg;
                    fst = 0;
                  }
                if (k < select[c].s2)
                  { e = SCAFFS[k].ectg-1;
                    lst = SCAFFS[k].slen;
                  }
    
                if (ori < 0)

                  { printf(">%s ",HEADERS+SCAFFS[k].hoff);
                    printf("%c%lld,%lld%c\n",fst==0?SOEL:SPOS,SCAFFS[k].slen-fst,
                                             SCAFFS[k].slen-lst,lst==SCAFFS[k].slen?EOEL:EPOS);

                    for (n = e; n >= b; n--)
                      { if (n == select[c].c1)
                          fst = select[c].p1;
                        else
                          fst = 0;
			if (n == select[c].c2)
                          lst = select[c].p2;
                        else
                          lst = CONTIGS[n].clen;
                        reverse_print(n,fst,lst,CONTIGS[n].sbeg);
                      }
                  }

                else // ori > 0

                  { printf(">%s ",HEADERS+SCAFFS[k].hoff);
                    printf("%c%d,%d%c\n",fst==0?SOEL:SPOS,fst,lst,lst==SCAFFS[k].slen?EOEL:EPOS);

                    for (n = b; n <= e; n++)
                      { if (n == select[c].c1)
                          fst = select[c].p1;
                        else
                          fst = 0;
			if (n == select[c].c2)
                          lst = select[c].p2;
                        else
                          lst = CONTIGS[n].clen;
                        forward_print(n,fst,lst,CONTIGS[n].sbeg);
                      }
                  }
              }
          }

        else  //  select[c].type == CONTG_SELECTION

          for (k = select[c].c1; k <= select[c].c2; k++)
            { GDB_CONTIG   *r = CONTIGS + k;
              GDB_SCAFFOLD *s = SCAFFS + r->scaf;
              if (k == select[c].c1)
                fst = select[c].p1;
              else
                fst = 0;
              if (k == select[c].c2)
                lst = select[c].p2;
              else
                lst = r->clen;

              if (ori < 0)
                { printf(">%s %c%lld,%lld%c :: Contig %d %c%d,%d%c\n",
                         HEADERS+s->hoff,r->sbeg+lst==s->slen?SOEL:SPOS,r->sbeg+lst,
                                         r->sbeg+fst,r->sbeg+fst==0?EOEL:EPOS,
                         k-(s->fctg)+1,lst==r->clen?SOEL:SPOS,lst,fst,fst==0?EOEL:EPOS);

                  reverse_print(k,fst,lst,0);
                }
              else
                { printf(">%s %c%lld,%lld%c :: Contig %d %c%d,%d%c\n",
                         HEADERS+s->hoff,r->sbeg+fst==0?SOEL:SPOS,r->sbeg+fst,
                                       r->sbeg+lst,r->sbeg+lst==s->slen?EOEL:EPOS,
                         k-(s->fctg)+1,fst==0?SOEL:SPOS,fst,lst,lst==r->clen?EOEL:EPOS);

                  forward_print(k,fst,lst,0);
                }
            }
      }
  }

  Free_ANO(ano);

  exit (0);
}
