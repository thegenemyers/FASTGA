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
            "[-h] [-w<int(80)>] <source:path>[.1gdb] [ <selection>|<FILE> ]";

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

static int Comp_Table[128] = 
  {   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,

      0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,

      0, 'T',   0, 'G',   0,   0,   0, 'C',
      0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0, 'A',   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,

      0, 't',   0, 'g',   0,   0,   0, 'c',
      0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0, 'a',   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,
   };

static void complement(char *seq, int len)
{ int i, j,  x;

  for (i = 0, j = len-1; i < j; i++, j--)
    { x = seq[i];
      seq[i] = Comp_Table[(int) seq[j]];
      seq[j] = Comp_Table[x];
    }
  if (i <= j)
    seq[i] = Comp_Table[(int) seq[i]];
}

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
            ARG_FLAGS("h")
            break;
          case 'w':
            ARG_NON_NEGATIVE(WIDTH,"Line width")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    DOSEQ = 1-flags['h'];

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
        fprintf(stderr,"\n");
        fprintf(stderr,"      -h: Show only the header lines.\n");
        fprintf(stderr,"      -w: Print -w bp per line (default is 80).\n");
        exit (1);
      }
  }

  //  Read GDB and establish sorted order of headers

  { int   s, c;
    char *sptr, *eptr;

    Read_GDB(gdb,argv[1]);

    NSCAFF  = gdb->nscaff;
    CONTIGS = gdb->contigs;
    SCAFFS  = gdb->scaffolds;
    HEADERS = gdb->headers;

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

    if (gdb->iscaps)
      UPPER = UPPER_CASE;
    else
      UPPER = LOWER_CASE;
  }

  //  Get the selection list (every contig by default)

  if (argc == 3)
    select = get_selection_list(argv[2],gdb,HASH,&llen);
  else
    select = get_selection_list(NULL,gdb,HASH,&llen);

  //  Display each contig or scaffold in the GDB according to the 
  //     range list of the selection, observing the display option

  { GDB_CONTIG   *r;
    GDB_SCAFFOLD *s;
    char         *contig;
    char         *nstring;
    int           c, k, x;
    int           fst, lst, ori;

    //  Setup 'n' string for printing scaffold gaps and buffer for contigs

    nstring = Malloc(WIDTH+1,"Allocating write buffer\n");
    if (UPPER == UPPER_CASE)
      for (c = 0; c < WIDTH; c++)
        nstring[c] = 'N';
    else
      for (c = 0; c < WIDTH; c++)
        nstring[c] = 'n';
    nstring[WIDTH] = '\0';

    contig = New_Contig_Buffer(gdb);

    for (c = 0; c < llen; c++)
      { ori = select[c].orient;

        if (select[c].type == SCAFF_SELECTION)

          { for (k = select[c].s1; k <= select[c].s2; k++)
    
              if (ori < 0)

                { int cbeg, cend, wpos;
                  int j, u, w, f, l;

                  if (k == select[c].s1)
                    fst = CONTIGS[select[c].c1].sbeg + select[c].p1;
                  else
                    fst = 0;
                  if (k == select[c].s2)
                    lst = CONTIGS[select[c].c2].sbeg + select[c].p2;
                  else
                    lst = SCAFFS[k].slen;

                  r   = CONTIGS + SCAFFS[k].ectg - 1;
  
                  printf(">%s ",HEADERS+SCAFFS[k].hoff);
                  printf("%c%lld,%lld%c\n",fst==0?SOEL:SPOS,SCAFFS[k].slen-fst,
                                           SCAFFS[k].slen-lst,lst==SCAFFS[k].slen?EOEL:EPOS);
  
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
                    { cbeg = SCAFFS[k].slen;
                      wpos = 0;
                      for (u = SCAFFS[k].ectg-1; u >= SCAFFS[k].fctg; u--, r--)
                        { cend = r->sbeg + r->clen;
                          if (cbeg > lst)
                            cbeg = lst;
                          if (cend < lst && cbeg > fst)
                            { if (cend >= fst)
                                l = cbeg-cend;
                              else
                                l = cbeg-fst;
  
                              ROLLING_PRINT(0,nstring)
                            }
  
                          cbeg = r->sbeg;
                          if (cbeg < lst && cend > fst)
                            { Get_Contig(gdb,u,UPPER,contig);
                              complement(contig,r->clen);
                              f = fst-cbeg;
                              l = lst-cbeg;
                              if (f < 0) f = 0; 
                              if (l > r->clen) l = r->clen;
                              x = r->clen-f;
                              f = r->clen-l;
                              l = x;
  
                              ROLLING_PRINT(f,contig+j)
                            }
                        }
  
                      cend = 0;
                      if (cbeg > lst)
                        cbeg = lst;
                      if (cend < lst && cbeg > fst)
                        { if (cend >= fst)
                            l = cbeg-cend;
                          else
                            l = cbeg-fst;
  
                          ROLLING_PRINT(0,nstring)
                        }
                      printf("\n");
                    }
                }

              else

                { int cbeg, cend, wpos;
                  int j, u, w, f, l;

                  if (k == select[c].s1)
                    fst = CONTIGS[select[c].c1].sbeg + select[c].p1;
                  else
                    fst = 0;
                  if (k == select[c].s2)
                    lst = CONTIGS[select[c].c2].sbeg + select[c].p2;
                  else
                    lst = SCAFFS[k].slen;
  
                  r = CONTIGS + SCAFFS[k].fctg;
  
                  printf(">%s",HEADERS+SCAFFS[k].hoff);
                  printf("%c%d,%d%c\n",fst==0?SOEL:SPOS,fst,lst,lst==SCAFFS[k].slen?EOEL:EPOS);
  
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

        else  //  select[c].type == CONTG_SELECTION

          for (k = select[c].c1; k <= select[c].c2; k++)
            { r = CONTIGS + k;
              s = SCAFFS + r->scaf;
              if (k == select[c].c1)
                fst = select[c].p1;
              else
                fst = 0;
              if (k == select[c].c2)
                lst = select[c].p2;
              else
                lst = r->clen;

              if (ori < 0)
                printf(">%s %c%lld,%lld%c :: Contig %d %c%d,%d%c\n",
                       HEADERS+s->hoff,r->sbeg+lst==s->slen?SOEL:SPOS,r->sbeg+lst,
                                       r->sbeg+fst,r->sbeg+fst==0?EOEL:EPOS,
                       k-(s->fctg)+1,lst==r->clen?SOEL:SPOS,lst,fst,fst==0?EOEL:EPOS);
              else
                printf(">%s %c%lld,%lld%c :: Contig %d %c%d,%d%c\n",
                       HEADERS+s->hoff,r->sbeg+fst==0?SOEL:SPOS,r->sbeg+fst,
                                       r->sbeg+lst,r->sbeg+lst==s->slen?EOEL:EPOS,
                       k-(s->fctg)+1,fst==0?SOEL:SPOS,fst,lst,lst==r->clen?EOEL:EPOS);

              if (DOSEQ)
                { int j;
  
                  Get_Contig(gdb,k,UPPER,contig);

                  if (ori < 0)
                    { x = r->clen - fst;
                      fst = r->clen - lst;
                      lst = x;
                      complement(contig,r->clen);
                    }

                  for (j = fst; j+WIDTH < lst; j += WIDTH)
                    printf("%.*s\n",WIDTH,contig+j);
                  if (j < lst)
                    printf("%.*s\n",lst-j,contig+j);
                }
              fflush(stdout);
            }
      }
  }

  Close_GDB(gdb);

  exit (0);
}
