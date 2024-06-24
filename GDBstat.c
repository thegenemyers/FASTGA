/*******************************************************************************************
 *
 *  Display statistics about the contents of a .gdb and a histogram of its read lengths.
 *
 *  Author:  Gene Myers
 *  Origin:  Reworking of DAZZ_DB program DBstat.c from the FASTGA module specialized to "gdb" ("dam")
 *  Date  :  Jan 2024
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "GDB.h"

static char *Usage = "[-h[<int>,<int>]] [-hlog] <source:path>[.1gdb]";

static GDB_CONTIG    *CONTIGS;
static GDB_SCAFFOLD  *SCAFFS;
static int           *GAPS;

static int CSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  return (CONTIGS[x].clen - CONTIGS[y].clen);
}

static int SSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  return (SCAFFS[x].slen - SCAFFS[y].slen);
}

static int GSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  return (GAPS[x] - GAPS[y]);
}

#define NBINS 20

static int nice_round(int num, int nbins, int *mod)
{ int buck;

  buck = 1;
  while (buck*nbins <= num)
    buck *= 10;
  if (buck >= 10)
    buck /= 10;
  *mod = 0;
  if (buck*nbins*5 <= num)
    { buck *= 5;
      *mod = 1;
    }
  else if (buck*nbins*2 <= num)
    { buck *= 2;
      *mod = 2;
    }
  return (buck);
}

int main(int argc, char *argv[])
{ GDB        _gdb, *gdb = &_gdb;
  int       *ctgsort, *scfsort, *gapsort;
  int        nctg, nscaff, ngap;
  int64      totbps, totspan, totgap;

  int     HIST_LIN;
  int     HIST_LOG;
  int     CBUCK, SBUCK;

  { int   i, j, k, n;
    int   flags[128];

    ARG_INIT("GDBstat")

    (void) flags;

    HIST_LIN = 0;
    HIST_LOG = 0;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'h':
            if (strcmp(argv[i]+2,"log") == 0)
              { HIST_LOG = 1;
                break;
              }
            HIST_LIN = 1;
            if (argv[i][2] == '\0')
              CBUCK = SBUCK = 0;
            else
              { if (sscanf(argv[i]+2,"%d,%d %n",&CBUCK,&SBUCK,&n) != 2 || argv[i][n+2] != '\0')
                  { fprintf(stderr,"%s: Cannot parse option %s as 2 comma separated int's.\n",
                                   Prog_Name,argv[i]);
                    exit (1);
                  }
                if (CBUCK <= 0 || SBUCK <= 0)
                  { fprintf(stderr,"%s: Bucket sizes must be positive int's in -h option.\n",
                                   Prog_Name);
                    exit (1);
                  }
              }
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -h: Display histograms of scaffold & contig lengths.\n");
        fprintf(stderr,"            int's give bucket sizes for respective histograms if given.\n");
        exit (1);
      }
  }

  { int   s, c;
    int64 spos;

    Read_GDB(gdb,argv[1]);

    nctg   = gdb->ncontig;
    nscaff = gdb->nscaff;

    CONTIGS  = gdb->contigs;
    SCAFFS   = gdb->scaffolds;

    GAPS = (int *) Malloc(sizeof(int)*(nctg+nscaff),"Allocating scaffold map");

    ngap = 0;
    for (s = 0; s < nscaff; s++)
      { spos = 0;
        for (c = SCAFFS[s].fctg; c < SCAFFS[s].ectg; c++)
           { if (CONTIGS[c].sbeg > spos)
               GAPS[ngap++] = CONTIGS[c].sbeg - spos;
             spos = CONTIGS[c].sbeg + CONTIGS[c].clen;
           }
        if (spos < SCAFFS[s].slen) 
          GAPS[ngap++] = SCAFFS[s].slen - spos;
      }

    totbps  = gdb->seqtot;
    totspan = 0;
    for (s = 0; s < nscaff; s++)
      totspan += SCAFFS[s].slen;
    totgap = totspan - totbps;

    scfsort = (int *) Malloc(sizeof(int)*nscaff,"Allocating N50 sort array\n");
    ctgsort = (int *) Malloc(sizeof(int)*nctg,"Allocating N50 sort array\n");
    gapsort = (int *) Malloc(sizeof(int)*ngap,"Allocating N50 sort array\n");

    for (c = 0; c < nctg; c++)
      ctgsort[c] = c;
    for (c = 0; c < nscaff; c++)
      scfsort[c] = c;
    for (c = 0; c < ngap; c++)
      gapsort[c] = c;

    qsort(ctgsort,nctg,sizeof(int),CSORT);
    qsort(gapsort,ngap,sizeof(int),GSORT);
    qsort(scfsort,nscaff,sizeof(int),SSORT);
  }

  //  Output overview

  { int cwide, swide, awide;

    cwide = Number_Digits(nctg);
    swide = Number_Digits(totspan);
    awide = Number_Digits(totspan/nscaff);
    cwide += (cwide-1)/3;
    swide += (swide-1)/3;
    awide += (awide-1)/3;

    printf("\nStatistics for assembly %s:\n",Root(argv[1],".gdb"));

    printf("\n  ");
    Print_Number(nscaff,cwide,stdout);
    printf(" scaffolds spanning ");
    Print_Number(totspan,swide,stdout);
    printf("bp, ave. = ");
    Print_Number(totspan/nscaff,awide,stdout);
    printf("bp\n");

    printf("  ");
    Print_Number(nctg,cwide,stdout);
    printf(" contigs containing ");
    Print_Number(totbps,swide,stdout);
    printf("bp, ave. = ");
    Print_Number(totbps/nctg,awide,stdout);
    printf("bp\n");

    if (ngap == 0)
      printf(" No gaps\n");
    else
      { printf("  ");
        Print_Number(ngap,cwide,stdout);
        printf(" gaps    containing ");
        Print_Number(totgap,swide,stdout);
        printf("bp, ave. = ");
        Print_Number(totgap/ngap,awide,stdout);
        printf("bp\n");
      }
  }

  //  Output N<X> table

  { int64 ss, cs, gs;
    int   sf, cf, gf;
    int   n;
    int   cwide, swide, gwide;

    cwide = Number_Digits(CONTIGS[ctgsort[nctg-1]].clen);
    swide = Number_Digits(SCAFFS[scfsort[nscaff-1]].slen);
    gwide = Number_Digits(GAPS[gapsort[ngap-1]]);
    cwide += (cwide-1)/3;
    swide += (swide-1)/3;
    gwide += (gwide-1)/3;
    if (cwide < (int) strlen("Contigs"))
      cwide = strlen("Contigs");
    if (swide < (int) strlen("Scaffolds"))
      swide = strlen("Scaffolds");

    sf = nscaff-1;
    ss = 0;
    cf = nctg-1;
    cs = 0;
    gf = ngap-1;
    gs = 0;
    printf("\n             Contigs%*sScaffolds%*sGaps\n",cwide-4,"",swide-6,"");
    printf("       MAX:  ");
    Print_Number(CONTIGS[ctgsort[cf]].clen,cwide,stdout);
    printf("   ");
    Print_Number(SCAFFS[scfsort[sf]].slen,swide,stdout);
    printf("   ");
    Print_Number(GAPS[gapsort[gf]],gwide,stdout);
    printf("\n");
    for (n = 10; n < 100; n += 10)
      { while (cf >= 0 && cs < totbps*(n/100.))
          { cs += CONTIGS[ctgsort[cf]].clen;
            cf -= 1;
          }
        while (sf >= 0 && ss < totspan*(n/100.))
          { ss += SCAFFS[scfsort[sf]].slen;
            sf -= 1;
          }
        while (gf >= 0 && gs < totgap*(n/100.))
          { gs += GAPS[gapsort[gf]];
            gf -= 1;
          }
        printf("       N%2d:  ",n);
        Print_Number(CONTIGS[ctgsort[cf+1]].clen,cwide,stdout);
        printf("   ");
        Print_Number(SCAFFS[scfsort[sf+1]].slen,swide,stdout);
        printf("   ");
        Print_Number(GAPS[gapsort[gf+1]],gwide,stdout);
        printf("\n");
      }
    printf("       MIN:  ");
    Print_Number(CONTIGS[ctgsort[0]].clen,cwide,stdout);
    printf("   ");
    Print_Number(SCAFFS[scfsort[0]].slen,swide,stdout);
    printf("   ");
    Print_Number(GAPS[gapsort[0]],gwide,stdout);
    printf("\n");
  }

  //  if HIST_LOG then output logarithmic histograms

  if (HIST_LOG)
    { int cbin, cmin, cmod;
      int sbin, smin, smod;
      int cf, ct;
      int sf, st;
      int64 cs, ss;
      int   ccwide, scwide;
      int   cwide, swide;
  
      cwide = Number_Digits(CONTIGS[ctgsort[nctg-1]].clen);
      swide = Number_Digits(SCAFFS[scfsort[nscaff-1]].slen);
      cwide += (cwide-1)/3;
      swide += (swide-1)/3;
      if (cwide < (int) strlen("Contigs"))
        cwide = strlen("Contigs");
      ccwide = Number_Digits(nctg);
      scwide = Number_Digits(nscaff);

      cmin = nice_round(CONTIGS[ctgsort[0]].clen,1,&cmod);
      cbin = nice_round(CONTIGS[ctgsort[nctg-1]].clen,1,&cmod);

      smin = nice_round(SCAFFS[scfsort[0]].slen,1,&smod);
      sbin = nice_round(SCAFFS[scfsort[nscaff-1]].slen,1,&smod);

      cf   = nctg-1;
      cs = 0;
      sf   = nscaff-1;
      ss = 0;
      printf("\n       Contigs%*sScaffolds\n",cwide+ccwide+13,"");
      while (cf >= 0 || sf >= 0)
        { ct = 0;
          while (cf >= 0 && CONTIGS[ctgsort[cf]].clen >= cbin)
            { ct += 1;
              cs += CONTIGS[ctgsort[cf]].clen;
              cf -= 1;
            }
          st = 0;
          while (sf >= 0 && SCAFFS[scfsort[sf]].slen >= sbin)
            { st += 1;
              ss += SCAFFS[scfsort[sf]].slen;
              sf -= 1;
            }
          printf("       ");
          if (cbin >= cmin)
            { Print_Number(cbin,cwide,stdout);
              printf(":  %*d   %5.1f%%",ccwide,ct,(100.*cs)/totbps);
            }
          else
            printf("%*s",cwide+ccwide+12,"");
          if (sbin >= smin)
            { printf("        ");
              Print_Number(sbin,swide,stdout);
              printf(":  %*d   %5.1f%%",scwide,st,(100.*ss)/totspan);
            }
          printf("\n");

          if (cmod == 1)
            cbin = (cbin*2)/5;
          else
            cbin = cbin/2;
          if (smod == 1)
            sbin = (sbin*2)/5;
          else
            sbin = sbin/2;
          cmod = (cmod+1)%3;
          smod = (smod+1)%3;
        }
    }

  //  if HIST_LIN then output classic histograms

  if (HIST_LIN)
    { int cbuck, cbin, cmin;
      int sbuck, sbin, smin;
      int cf, ct;
      int sf, st;
      int64 cs, ss;
      int   ccwide, scwide;
      int   cwide, swide;
  
      cwide = Number_Digits(CONTIGS[ctgsort[nctg-1]].clen);
      swide = Number_Digits(SCAFFS[scfsort[nscaff-1]].slen);
      cwide += (cwide-1)/3;
      swide += (swide-1)/3;
      if (cwide < (int) strlen("Contigs"))
        cwide = strlen("Contigs");
      ccwide = Number_Digits(nctg);
      scwide = Number_Digits(nscaff);
      
      if (CBUCK == 0)
        { cbin = CONTIGS[ctgsort[nctg-1]].clen - CONTIGS[ctgsort[0]].clen;
          sbin = SCAFFS[scfsort[nscaff-1]].slen - SCAFFS[scfsort[0]].slen;
  
          cbuck = nice_round(cbin,NBINS,&cbin);
          sbuck = nice_round(sbin,NBINS,&sbin);
        }
      else
        { cbuck = CBUCK;
          sbuck = SBUCK;
        }
  
      cf   = nctg-1;
      cbin = (CONTIGS[ctgsort[cf]].clen/cbuck)*cbuck;
      cmin = (CONTIGS[ctgsort[0]].clen/cbuck)*cbuck;
      cs = 0;
      sf   = nscaff-1;
      sbin = (SCAFFS[scfsort[sf]].slen/sbuck)*sbuck;
      smin = (SCAFFS[scfsort[0]].slen/sbuck)*sbuck;
      ss = 0;
      printf("\n       Contigs%*sScaffolds\n",cwide+ccwide+13,"");
      while (cf >= 0 || sf >= 0)
        { ct = 0;
          while (cf >= 0 && CONTIGS[ctgsort[cf]].clen >= cbin)
            { ct += 1;
              cs += CONTIGS[ctgsort[cf]].clen;
              cf -= 1;
            }
          st = 0;
          while (sf >= 0 && SCAFFS[scfsort[sf]].slen >= sbin)
            { st += 1;
              ss += SCAFFS[scfsort[sf]].slen;
              sf -= 1;
            }
          printf("       ");
          if (cbin >= cmin)
            { Print_Number(cbin,cwide,stdout);
              printf(":  %*d   %5.1f%%",ccwide,ct,(100.*cs)/totbps);
            }
          else
            printf("%*s",cwide+ccwide+12,"");
          if (sbin >= smin)
            { printf("        ");
              Print_Number(sbin,swide,stdout);
              printf(":  %*d   %5.1f%%",scwide,st,(100.*ss)/totspan);
            }
          printf("\n");
  
          cbin -= cbuck;
          sbin -= sbuck;
        }
    }

  free(gapsort);
  free(ctgsort);
  free(scfsort);
  free(GAPS);
  Close_GDB(gdb);

  exit (0);
}
