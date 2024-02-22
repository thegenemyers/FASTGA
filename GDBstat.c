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

#include "DB.h"

static char *Usage = "[-h[<int>,<int>]] [-hlog] <source:path>[.gdb]";

static DAZZ_READ *READS;

static int CSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  return (READS[x].rlen - READS[y].rlen);
}

static int GSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  return (READS[x].fpulse - READS[y].fpulse);
}

static int64 *SCAFLEN;

static int SSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  return (SCAFLEN[x] - SCAFLEN[y]);
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
{ DAZZ_DB    _db, *db = &_db;
  DAZZ_READ *reads;
  int64     *slen;
  int       *ctgsort, *scfsort, *gapsort;
  int        nctg, nscaff, ngap;
  int        bgap, egap;
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

  { int status;
    int s, r, g;

    status = Open_DB(argv[1],db);
    if (status < 0)
      exit (1);
    if (status == 0)
      { fprintf(stderr,"%s: Cannot open %s as a .gdb\n",Prog_Name,argv[1]);
        exit (1);
      }

    nctg   = db->nreads;
    reads  = db->reads;

    slen = (int64 *) Malloc(sizeof(int64)*2*nctg,"Allocating scaffold map");
    if (slen == NULL)
      exit (1);

    bgap = egap = 0;
    s = -1;
    for (r = 0; r < nctg; r++)
      { if (reads[r].origin == 0)
          { s += 1;
            if (reads[r].fpulse > 0)
              bgap += 1;
            g = 0;
          }
        else
          g = slen[s];
        slen[s] = reads[r].fpulse + reads[r].rlen;
        reads[r].fpulse -= g;
        if (reads[r].rlen == 0)
          egap += 1;
      }
    nscaff = s+1;

    totbps  = db->totlen;
    totspan = 0;
    for (r = 0; r < nscaff; r++)
      totspan += slen[r];
    totgap = totspan - totbps;

    scfsort = (int *) Malloc(sizeof(int)*nscaff,"Allocating N50 sort array\n");
    ctgsort = (int *) Malloc(sizeof(int)*nctg,"Allocating N50 sort array\n");
    gapsort = (int *) Malloc(sizeof(int)*nctg,"Allocating N50 sort array\n");
    if (scfsort == NULL || ctgsort == NULL || gapsort == NULL)
      exit (1);

    for (r = 0; r < nctg; r++)
      ctgsort[r] = gapsort[r] = r;
    for (r = 0; r < nscaff; r++)
      scfsort[r] = r;

    READS = db->reads;
    qsort(ctgsort,nctg,sizeof(int),CSORT);
    qsort(gapsort,nctg,sizeof(int),GSORT);

    SCAFLEN = slen;
    qsort(scfsort,nscaff,sizeof(int),SSORT);
  }

  //  Output overview

  { int cwide, swide, awide;

    ngap = (nctg - nscaff) + bgap;

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
    Print_Number(nctg-egap,cwide,stdout);
    printf(" contigs containing ");
    Print_Number(totbps,swide,stdout);
    printf("bp, ave. = ");
    Print_Number(totbps/(nctg-egap),awide,stdout);
    printf("bp\n");

    printf("  ");
    Print_Number(ngap,cwide,stdout);
    printf(" gaps    containing ");
    Print_Number(totgap,swide,stdout);
    printf("bp, ave. = ");
    Print_Number(totgap/ngap,awide,stdout);
    printf("bp\n");
  }

  //  Output N<X> table

  { int64 ss, cs, gs;
    int   sf, cf, gf;
    int   n;
    int   cwide, swide, gwide;

    cwide = Number_Digits(reads[ctgsort[nctg-1]].rlen);
    swide = Number_Digits(slen[scfsort[nscaff-1]]);
    gwide = Number_Digits(reads[gapsort[nctg-1]].fpulse);
    cwide += (cwide-1)/3;
    swide += (swide-1)/3;
    gwide += (gwide-1)/3;
    if (cwide < (int) strlen("Contigs"))
      cwide = strlen("Contigs");
    if (swide < (int) strlen("Scaffolds"))
      swide = strlen("Scaffolds");

    ngap = nctg - ngap;

    sf = nscaff-1;
    ss = 0;
    cf = nctg-1;
    cs = 0;
    gf = nctg-1;
    gs = 0;
    printf("\n             Contigs%*sScaffolds%*sGaps\n",cwide-4,"",swide-6,"");
    printf("       MAX:  ");
    Print_Number(reads[ctgsort[cf]].rlen,cwide,stdout);
    printf("   ");
    Print_Number(slen[scfsort[sf]],swide,stdout);
    printf("   ");
    Print_Number(reads[gapsort[gf]].fpulse,gwide,stdout);
    printf("\n");
    for (n = 10; n < 100; n += 10)
      { while (cf >= egap && cs < totbps*(n/100.))
          { cs += reads[ctgsort[cf]].rlen;
            cf -= 1;
          }
        while (sf >= 0 && ss < totspan*(n/100.))
          { ss += slen[scfsort[sf]];
            sf -= 1;
          }
        while (gf >= ngap && gs < totgap*(n/100.))
          { gs += reads[gapsort[cf]].fpulse;
            gf -= 1;
          }
        printf("       N%2d:  ",n);
        Print_Number(reads[ctgsort[cf+1]].rlen,cwide,stdout);
        printf("   ");
        Print_Number(slen[scfsort[sf+1]],swide,stdout);
        printf("   ");
        Print_Number(reads[gapsort[gf+1]].fpulse,gwide,stdout);
        printf("\n");
      }
    printf("       MIN:  ");
    Print_Number(reads[ctgsort[egap]].rlen,cwide,stdout);
    printf("   ");
    Print_Number(slen[scfsort[0]],swide,stdout);
    printf("   ");
    Print_Number(reads[gapsort[ngap]].fpulse,gwide,stdout);
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
  
      cwide = Number_Digits(reads[ctgsort[nctg-1]].rlen);
      swide = Number_Digits(slen[scfsort[nscaff-1]]);
      cwide += (cwide-1)/3;
      swide += (swide-1)/3;
      if (cwide < (int) strlen("Contigs"))
        cwide = strlen("Contigs");
      ccwide = Number_Digits(nctg);
      scwide = Number_Digits(nscaff);

      cmin = nice_round(reads[ctgsort[egap]].rlen,1,&cmod);
      cbin = nice_round(reads[ctgsort[nctg-1]].rlen,1,&cmod);

      smin = nice_round(slen[scfsort[0]],1,&smod);
      sbin = nice_round(slen[scfsort[nscaff-1]],1,&smod);

      cf   = nctg-1;
      cs = 0;
      sf   = nscaff-1;
      ss = 0;
      printf("\n       Contigs%*sScaffolds\n",cwide+ccwide+13,"");
      while (cf >= egap || sf >= 0)
        { ct = 0;
          while (cf >= egap && reads[ctgsort[cf]].rlen >= cbin)
            { ct += 1;
              cs += reads[ctgsort[cf]].rlen;
              cf -= 1;
            }
          st = 0;
          while (sf >= 0 && slen[scfsort[sf]] >= sbin)
            { st += 1;
              ss += slen[scfsort[sf]];
              sf -= 1;
            }
          printf("       ");
          if (cbin >= cmin)
            { Print_Number(cbin,cwide,stdout);
              printf(":  %*d   %5.1f%%",ccwide,ct,(100.*cs)/totbps);
            }
          else
            printf("%*s",ccwide+12,"");
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
  
      cwide = Number_Digits(reads[ctgsort[nctg-1]].rlen);
      swide = Number_Digits(slen[scfsort[nscaff-1]]);
      cwide += (cwide-1)/3;
      swide += (swide-1)/3;
      if (cwide < (int) strlen("Contigs"))
        cwide = strlen("Contigs");
      ccwide = Number_Digits(nctg);
      scwide = Number_Digits(nscaff);
      
      if (CBUCK == 0)
        { cbin = reads[ctgsort[nctg-1]].rlen - reads[ctgsort[egap]].rlen;
          sbin = slen[scfsort[nscaff-1]] - slen[scfsort[0]];
  
          cbuck = nice_round(cbin,NBINS,&cbin);
          sbuck = nice_round(sbin,NBINS,&sbin);
        }
      else
        { cbuck = CBUCK;
          sbuck = SBUCK;
        }
  
      cf   = nctg-1;
      cbin = (reads[ctgsort[cf]].rlen/cbuck)*cbuck;
      cmin = (reads[ctgsort[0]].rlen/cbuck)*cbuck;
      cs = 0;
      sf   = nscaff-1;
      sbin = (slen[scfsort[sf]]/sbuck)*sbuck;
      smin = (slen[scfsort[0]]/sbuck)*sbuck;
      ss = 0;
      printf("\n       Contigs%*sScaffolds\n",cwide+ccwide+13,"");
      while (cf >= egap || sf >= 0)
        { ct = 0;
          while (cf >= egap && reads[ctgsort[cf]].rlen >= cbin)
            { ct += 1;
              cs += reads[ctgsort[cf]].rlen;
              cf -= 1;
            }
          st = 0;
          while (sf >= 0 && slen[scfsort[sf]] >= sbin)
            { st += 1;
              ss += slen[scfsort[sf]];
              sf -= 1;
            }
          printf("       ");
          if (cbin >= cmin)
            { Print_Number(cbin,cwide,stdout);
              printf(":  %*d   %5.1f%%",ccwide,ct,(100.*cs)/totbps);
            }
          else
            printf("%*s",ccwide+12,"");
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

  free(ctgsort);
  free(scfsort);
  free(slen);
  Close_DB(db);

  exit (0);
}
