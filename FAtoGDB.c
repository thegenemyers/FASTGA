/*******************************************************************************************
 *
 *  Convert a .fasta or .1seq file to a .gdb:
 *
 *  Author:   Gene Myers
 *            Modified by Richard Durbin to also read .1seq files
 *  Origin:   Complete rework/specialization from fasta2DAM of DAZZ_DB for FASTGA module
 *  Creation: Jan 2024
 *  Last Mod: Mar 2024
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <dirent.h>

#include "GDB.h"
#include "ANO.h"

static char *Usage = "[-v] [-L:<log:path>] [-n<int>] <source:path>[<fa_extn>|<1_extn>] [<target:path>[.1gdb]]";

int main(int argc, char *argv[])
{ char *spath, *tpath;
  char *TPATH, *TROOT;
  int   ftype;
  int   NCUT;
  GDB   gdb;
  ANO   ano;

  int   VERBOSE;
  FILE *LOG_FILE;

  //   Process command line

  { int   i, j, k;
    int   flags[128];
    char *eptr;

    ARG_INIT("FAtoGDB")

    NCUT = 0;
    LOG_FILE = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'n':
            ARG_POSITIVE(NCUT,"n-run cutoff");
            break;
          case 'L':
            if (argv[i][2] != ':')
              { fprintf (stderr,"%s: option -L must be followed by :<filename>\n",Prog_Name);
                exit (1);
              }
            LOG_FILE = fopen(argv[i]+3,"a");
            if (LOG_FILE == NULL)
              { fprintf (stderr,"%s: Cannot open logfile %s for output\n",Prog_Name,argv[i]+3);
                exit(1);
              }
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"       -n: Turn runs of n's of length < # into a's.\n");
        fprintf(stderr,"       -L: Output log to specified file.\n");
        exit (1);
      }
  }

  if (LOG_FILE)
    fprintf(LOG_FILE,"\n%s\n",Command_Line);

  if (VERBOSE || LOG_FILE)
    StartTime();

  //  Determine source and target root names, paths, and extensions

  if (argc == 2)
    ftype = Get_GDB_Paths(argv[1],NULL,&spath,&tpath,1);
  else
    ftype = Get_GDB_Paths(argv[1],argv[2],&spath,&tpath,1);

  TPATH = PathTo(tpath);
  TROOT = Root(tpath,NULL);
  if (VERBOSE)
    { if (strcmp(TPATH,".") == 0)
        fprintf(stderr,"\n  Creating genome data base (GDB) %s.1gdb in the current directory\n",
                       TROOT);  
      else
        fprintf(stderr,"\n  Creating genome data base (GDB) %s.1gdb in directory %s\n",
                       TROOT,TPATH);  
      fflush(stderr);
    }

  Create_GDB(&gdb,spath,ftype,1,tpath,NCUT,&ano);

  Write_GDB(&gdb,tpath);
  if (ano.nints > 0)
    { if (VERBOSE)
        { fprintf(stderr,"\n  Masked sequence detected, also creating .ano file %s.1ano\n",TROOT);
          fflush(stderr);
        }
      Write_ANO(&ano,Catenate(TPATH,"/",TROOT,".1ano"),100);
      Free_ANO(&ano);
    }

  Close_GDB(&gdb);

  free(TROOT);
  free(TPATH);

  free(tpath);
  free(spath);

  if (VERBOSE)
    TimeTo(stderr,1,0);
  if (LOG_FILE)
    { TimeTo(LOG_FILE,1,0);
      fclose(LOG_FILE);
    }

  exit (0);
}
