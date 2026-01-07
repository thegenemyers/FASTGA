/********************************************************************************************
 *
 *  Recreate the .fasta file that gave rise to a given .gdb
 *
 *  Author:  Gene Myers
 *  Origin:  Reworking of DAZZ_DB program DAM2fasta.c for the FASTGA module
 *  Date  :  Jan 2024
 *
 ********************************************************************************************/
/*  Last edited: Jul 25 19:12 2024 (rd109) */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#include "GDB.h"
#include "ANO.h"

static char *Usage = "[-v] <source:path>[.1ano] [ <target:path>[.bed] ]";

int main(int argc, char *argv[])
{ ANO        _ano, *ano = &_ano;
  char       *TPATH, *TROOT;
  char       *SROOT;
  FILE       *output;

  int   VERBOSE; // -v

  //  Process arguments

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("ANOtoBED")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose output\n");
        exit (1);
      }
  }

  //  Open ano & determine source and target parts

  { struct stat tdesc;
    char       *p, *e;

    Read_ANO(ano,argv[1],NULL);

    SROOT = Root(argv[1],".1ano");

    if (argc == 2)
      output = stdout;
    else
      { TPATH = argv[2];
        if (stat(TPATH,&tdesc) >= 0 && (tdesc.st_mode & S_IFMT) == S_IFDIR)
          TROOT = SROOT;
        else
          { p = rindex(TPATH,'/');
            if (p == NULL)
              { TROOT = TPATH;
                TPATH = ".";
              }
            else
              { *p++ = '\0';
                TROOT = p;
              } 
          }

        e = TROOT + strlen(TROOT);
        if (strcmp(e-4,".bed") == 0)
          e[-4] = '\0';

        if (VERBOSE)
          { if (strcmp(TPATH,".") == 0)
              fprintf(stderr,"\n  Creating bed file %s.bed in current directory\n",TROOT);
            else
              fprintf(stderr,"\n  Creating bed file %s.bed at directory %s\n",TROOT,TPATH);
            fflush(stderr);
          }

        output = fopen(Catenate(TPATH,"/",TROOT,".bed"),"w");
        if (output == NULL)
          { fprintf(stderr,"Could not create/open %s/%s.bed for writing\n",TPATH,TROOT);
            exit (1);
          }
      }
  }

  { GDB_CONTIG   *ctg;
    GDB_SCAFFOLD *scf;
    ANO_PAIR     *mask;
    int           nctg;
    int64        *moff;
    char         *hdr, *h;
    char          date[20];
    int           c, m, t;
    time_t        now;

    scf  = ano->gdb->scaffolds;
    ctg  = ano->gdb->contigs;
    hdr  = ano->gdb->headers;
    nctg = ano->gdb->ncontig;

    mask = ano->masks;
    moff = ano->moff;

    printf("# Provenance:\n");
    for (c = 0; c < ano->nprov; c++)
      printf("#  %s  %s\n",ano->prov[c].command,ano->prov[c].date);
    now = time(NULL);
    strftime(date, 20, "%F_%T", localtime(&now));
    printf("#  %s  %s\n",Command_Line,date);

    for (c = 0; c < nctg; c++)
      { t = moff[c+1];
        h = hdr+scf[ctg[c].scaf].hoff;
        for (m = moff[c]; m < t; m++)
          { fprintf(output,"%s\t%lld\t%lld\t",h,mask[m].beg,mask[m].end);
            if (mask[m].label != NULL)
              fprintf(output,"%s",mask[m].label);
            if (mask[m].beg <= mask[m].end)
              fprintf(output,"\t%d\t+\n",mask[m].score);
            else
              fprintf(output,"\t%d\t-\n",mask[m].score);
          }
      }
  }

  free(SROOT);

  Free_ANO(ano);

  exit (0);
}
