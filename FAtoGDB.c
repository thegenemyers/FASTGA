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

static char *Usage = "[-v] <source:path>[<fa_extn>|<1_extn>] [<target:path>[.1gdb]]";

int main(int argc, char *argv[])
{ char *spath, *tpath;
  char *TPATH, *TROOT;
  int   ftype;
  GDB   gdb;

  int VERBOSE;

  //   Process command line

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("FAtoGDB")

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
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        exit (1);
      }
  }

  //  Determine source and target root names, paths, and extensions

  if (argc == 2)
    ftype = Get_GDB_Paths(argv[1],NULL,&spath,&tpath,1);
  else
    ftype = Get_GDB_Paths(argv[1],argv[2],&spath,&tpath,1);

  if (VERBOSE)
    { TPATH = PathTo(tpath);
      TROOT = Root(tpath,NULL);
      if (strcmp(TPATH,".") == 0)
        fprintf(stderr,"\n  Creating genome data base (GDB) %s.1gdb in the current directory\n",
                       TROOT);  
      else
        fprintf(stderr,"\n  Creating genome data base (GDB) %s.1gdb in directory %s\n",
                       TROOT,TPATH);  
      fflush(stderr);
      free(TROOT);
      free(TPATH);
    }

  Create_GDB(&gdb,spath,ftype,1,tpath);

  Write_GDB(&gdb,tpath);

  free(tpath);
  free(spath);

  exit (0);
}
