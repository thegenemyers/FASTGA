/*******************************************************************************************
 *
 *  Display statistics about the contents of a .ano and a histogram of its interval lengths.
 *
 *  Author:  Gene Myers
 *  Origin:  Reworking of GDBstat program for new .1ano data object
 *  Date  :  Dec 2025
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ANO.h"

static char *Usage = "[-h[<int>,<int>]] [-hlog] <source:path>[.1ano]";

static int ISORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  return (x - y);
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
{ ANO        _ano, *ano = &_ano;
  GDB        *gdb;
  int        *region, nreg;
  int        *covered, ncov;
  int        *uncovered, nunc;
  int64       totreg, totcov, totunc, totgap;
  int64       numori, numlab, numscr;
  int64       nints;

  int     HIST_LIN;
  int     HIST_LOG;
  int     RBUCK, CBUCK;

  { int   i, j, k, n;
    int   flags[128];

    ARG_INIT("ANOstat")

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
              RBUCK = CBUCK = 0;
            else
              { if (sscanf(argv[i]+2,"%d,%d %n",&RBUCK,&CBUCK,&n) != 2 || argv[i][n+2] != '\0')
                  { fprintf(stderr,"%s: Cannot parse option %s as 2 comma separated int's.\n",
                                   Prog_Name,argv[i]);
                    exit (1);
                  }
                if (RBUCK <= 0 || CBUCK <= 0)
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
        fprintf(stderr,"      -h: Display histogram of interval lengths.\n");
        fprintf(stderr,"            int's give bucket sizes for respective histograms if given.\n");
        exit (1);
      }
  }

  //  Open .ano and compute stats and sorted size arrays

  { int c;

    Read_ANO(ano,argv[1],NULL);

    gdb   = ano->gdb;
    nints = ano->nints;

    region = malloc(sizeof(int)*3*nints + gdb->ncontig);
    if (region == NULL)
      { fprintf(stderr,"%s: Not enough memory\n",Prog_Name);
        exit (1);
      }
    covered   = region + nints;
    uncovered = covered + nints;

    totreg = totcov = totunc = totgap = 0;
    numori = numlab = numscr = 0;
    nreg = ncov = nunc = 0;

    for (c = 0; c < gdb->ncontig; c++)
      { ANO_PAIR *m, *t;
        int b, e;
        int beg, end;

        t = ano->masks + ano->moff[c+1];
        m = ano->masks + ano->moff[c];
        if (m == t)
          continue;
        b = m->beg;
        e = m->end;
        region[nreg++] = e-b;
        totreg += e-b;
        numori += m->orient;
        numlab += (m->label != NULL);
        numscr += (m->score > 0);
        if (b > 0)
          { uncovered[nunc++] = b;
            totunc += b;
          }
        for (m++; m < t; m++)
          { beg = m->beg;
            end = m->end;
            if (e < beg)
              { covered[ncov++] = e-b;
                totcov += e-b;
                b = beg;
                uncovered[nunc++] = b-e;
                totunc += b-e;
                e = end;
              }
            else
              { if (end > e)
                  e = end;
              }
            region[nreg++] = end-beg;
            totreg += end-beg;
            numori += m->orient;
            numlab += (m->label != NULL);
            numscr += (m->score > 0);
          }
        covered[ncov++] = e-b;
        totcov += e-b;
        end = gdb->contigs[c].clen;
        if (e < end)
          { uncovered[nunc++] = end-e;
            totunc += end-e;
          }
        else if (e > end)
          totgap += e-end;
      }

    qsort(region,nreg,sizeof(int),ISORT);
    qsort(covered,ncov,sizeof(int),ISORT);
    qsort(uncovered,nunc,sizeof(int),ISORT);
  }

  //  Output overview

  { printf("\nStatistics for ano file %s:\n",Root(argv[1],".ano"));

    printf("\n  ");
    printf("There are ");
    Print_Number(nints,0,stdout);
    if (numori != 0)
      printf(" oriented");
    else
      printf(" unoriented");
    if (numlab == nints)
      printf(", labelled");
    else if (numlab == 0)
      printf(", unlabelled");
    if (numscr != 0)
      printf(", scored");
    else
      printf(", unscored");
    printf(" intervals");
    if (numlab != nints && numlab != 0)
      { printf(" of which ");
        Print_Number(numlab,0,stdout);
        printf(" are labelled");
      }
    printf("\n");

    printf("\n  ");
    if (totcov == totreg)
      printf("The intervals are all disjoint\n");
    else
      printf("%.1f%% of the interval regions overlap\n",(100.*(totreg-totcov))/totreg);

    printf("\n  The intervals span ");
    if (totreg >= 1000000)
      { Print_Number(totreg/1000000,0,stdout);
        printf(".%lldM",(totreg % 1000000)/100000);
      }
    else if (totreg >= 1000)
      { Print_Number(totreg/1000,0,stdout);
        printf(".%lldK",(totreg % 1000)/100);
      }
    else
      Print_Number(totreg,0,stdout);
    printf("bp and cover ");
    if (totcov >= 1000000)
      { Print_Number(totcov/1000000,0,stdout);
        printf(".%lldM",(totcov % 1000000)/100000);
      }
    else if (totcov >= 1000)
      { Print_Number(totcov/1000,0,stdout);
        printf(".%lldK",(totcov % 1000)/100);
      }
    else
      Print_Number(totcov,0,stdout);
    printf("bp (%.1f%%) of the genome\n",(100.*totcov)/gdb->seqtot);

    if (totgap != 0)
      { printf("\n  The intervals span ");
        Print_Number(totgap,0,stdout);
        printf("bp of the gaps between contigs\n");
      }
    else
      printf("\n  The intervals do not span gaps between contigs\n");
  }

  //  Output N<X> table

  { int64 rs, cs, us;
    int   nr, nc, nu;
    int   n;
    int   rwide, cwide, uwide;

    rwide = Number_Digits(region[nreg-1]);
    rwide += (rwide-1)/3;
    if (rwide < 9)
      rwide = 9;

    cwide = Number_Digits(covered[ncov-1]);
    cwide += (cwide-1)/3;
    if (cwide < 14)
      cwide = 14;

    uwide = Number_Digits(uncovered[nunc-1]);
    uwide += (uwide-1)/3;
    if (uwide < 16)
      uwide = 16;

    nr = nreg-1;
    rs = 0;
    nc = ncov-1;
    cs = 0;
    nu = nunc-1;
    us = 0;
    printf("\n             Intervals%*sCovered Blocks%*sUncovered Blocks\n",rwide-6,"",cwide-11,"");
    printf("       MAX:  ");
    Print_Number(region[nr],rwide,stdout);
    printf("   ");
    Print_Number(covered[nc],cwide,stdout);
    printf("   ");
    Print_Number(uncovered[nu],uwide,stdout);
    printf("\n");
    for (n = 10; n < 100; n += 10)
      { while (nr >= 0 && rs < totreg*(n/100.))
          { rs += region[nr];
            nr -= 1;
          }
        printf("       N%2d:  ",n);
        Print_Number(region[nr+1],rwide,stdout);

        while (nc >= 0 && cs < totcov*(n/100.))
          { cs += covered[nc];
            nc -= 1;
          }
        printf("   ");
        Print_Number(covered[nc+1],cwide,stdout);

        while (nu >= 0 && us < totunc*(n/100.))
          { us += uncovered[nu];
            nu -= 1;
          }
        printf("   ");
        Print_Number(uncovered[nu+1],uwide,stdout);

        printf("\n");
      }
    printf("       MIN:  ");
    Print_Number(region[0],rwide,stdout);
    printf("   ");
    Print_Number(covered[0],cwide,stdout);
    printf("   ");
    Print_Number(uncovered[0],uwide,stdout);
    printf("\n");
  }

  //  if HIST_LOG then output logarithmic histograms

  if (HIST_LOG)
    { int rbin, rmin, rmod;
      int cbin, cmin, cmod;
      int nr, rt;
      int nc, ct;
      int64 rs, cs;
      int   rcwide, ccwide;
      int   rwide, cwide;
  
      rwide = Number_Digits(region[nreg-1]);
      cwide = Number_Digits(covered[ncov-1]);
      rwide += (rwide-1)/3;
      cwide += (cwide-1)/3;
      if (rwide < (int) strlen("Intervals"))
        rwide = strlen("Intervals");
      rcwide = Number_Digits(nreg);
      ccwide = Number_Digits(ncov);

      rmin = nice_round(region[0],1,&rmod);
      rbin = nice_round(region[nreg-1],1,&rmod);

      cmin = nice_round(covered[0],1,&cmod);
      cbin = nice_round(covered[ncov-1],1,&cmod);

      nr   = nreg-1;
      rs = 0;
      nc   = ncov-1;
      cs = 0;
      printf("\n       Intervlas%*sCovered Blocks\n",rwide+rcwide+13,"");
      while (nr >= 0 || nc >= 0)
        { rt = 0;
          while (nr >= 0 && region[nr] >= rbin)
            { rt += 1;
              rs += region[nr];
              nr -= 1;
            }
          ct = 0;
          while (nc >= 0 && covered[nc] >= cbin)
            { ct += 1;
              cs += covered[nc];
              nc -= 1;
            }
          printf("       ");
          if (rbin >= rmin)
            { Print_Number(rbin,rwide,stdout);
              printf(":  %*d   %5.1f%%",rcwide,rt,(100.*rs)/totreg);
            }
          else
            printf("%*s",cwide+rcwide+12,"");
          if (cbin >= cmin)
            { printf("        ");
              Print_Number(cbin,cwide,stdout);
              printf(":  %*d   %5.1f%%",ccwide,ct,(100.*cs)/totcov);
            }
          printf("\n");

          if (rmod == 1)
            rbin = (rbin*2)/5;
          else
            rbin = rbin/2;
          if (cmod == 1)
            cbin = (cbin*2)/5;
          else
            cbin = cbin/2;
          rmod = (rmod+1)%3;
          cmod = (cmod+1)%3;
        }
    }

  //  if HIST_LIN then output classic histograms

  if (HIST_LIN)
    { int rbuck, rbin, rmin;
      int cbuck, cbin, cmin;
      int nr, rt;
      int nc, ct;
      int64 rs, cs;
      int   rcwide, ccwide;
      int   rwide, cwide;
  
      rwide = Number_Digits(region[nreg-1]);
      cwide = Number_Digits(covered[ncov-1]);
      rwide += (rwide-1)/3;
      cwide += (cwide-1)/3;
      if (rwide < (int) strlen("Intervals"))
        rwide = strlen("Intervals");
      rcwide = Number_Digits(nreg);
      ccwide = Number_Digits(ncov);

      if (CBUCK == 0)
        { rbin = region[nreg-1] - region[0];
          cbin = covered[ncov-1] - covered[0];

          rbuck = nice_round(rbin,NBINS,&rbin);
          cbuck = nice_round(cbin,NBINS,&cbin);
        }
      else
        { rbuck = RBUCK;
          cbuck = CBUCK;
        }
      
      nr   = nreg-1;
      rbin = (region[nr]/rbuck)*rbuck;
      rmin = (region[0]/rbuck)*rbuck;
      rs = 0;
      nc   = ncov-1;
      cbin = (covered[nc]/cbuck)*cbuck;
      cmin = (covered[0]/cbuck)*cbuck;
      cs = 0;
      printf("\n       Intervals%*sCovered_Blocks\n",cwide+ccwide+13,"");
      while (nr >= 0 || nc >= 0)
        { rt = 0;
          while (nr >= 0 && region[nr] >= rbin)
            { rt += 1;
              rs += region[nr];
              nr -= 1;
            }
          ct = 0;
          while (nc >= 0 && covered[nc] >= cbin)
            { ct += 1;
              cs += covered[nc];
              nc -= 1;
            }
          printf("       ");
          if (rbin >= rmin)
            { Print_Number(rbin,rwide,stdout);
              printf(":  %*d   %5.1f%%",rcwide,rt,(100.*rs)/totreg);
            }
          else
            printf("%*s",rwide+rcwide+12,"");
          if (cbin >= cmin)
            { printf("        ");
              Print_Number(cbin,cwide,stdout);
              printf(":  %*d   %5.1f%%",ccwide,ct,(100.*cs)/totcov);
            }
          printf("\n");
  
          rbin -= rbuck;
          cbin -= cbuck;
        }
    }

  free(region);
  Free_ANO(ano);

  exit (0);
}
