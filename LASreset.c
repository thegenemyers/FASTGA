/*******************************************************************************************
 *
 *  Utility for displaying the overlaps in a .las file in a variety of ways including
 *    a minimal listing of intervals, a cartoon, and a full out alignment.
 *
 *  Author:    Gene Myers
 *  Creation:  July 2013
 *  Last Mod:  Jan 2015
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

#include "DB.h"

static char *Usage = "<alignments:path>[.las] <DB1:path>[.gdb] [<DB2:path>[.gdb]]";

int main(int argc, char *argv[])
{ char *SPATH1, *SROOT1;
  char *SPATH2, *SROOT2;
  char *TEMP;
  char *APATH, *AROOT;
  int   tspace;
  int64 novl;
  FILE *input, *output;

  //  Process options

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("LASreset")

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

    if (argc < 3 || argc > 4)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  //  Find .gdb's and parse their path

  { SPATH1 = PathTo(argv[2]);
    SROOT1 = Root(argv[2],".gdb");
    input = Fopen(Catenate(SPATH1,"/",SROOT1,".gdb"),"r");
    if (input == NULL)
      exit (1);
    fclose(input);

    if (argc == 4)
      { SPATH2 = PathTo(argv[3]);
        SROOT2 = Root(argv[3],".gdb");
        input = Fopen(Catenate(SPATH2,"/",SROOT2,".gdb"),"r");
        if (input == NULL)
          exit (1);
        fclose(input);
      }
  }

  //  Open .las file reading and read header information
  
  { int  nlen, one;

    APATH = PathTo(argv[1]);
    AROOT = Root(argv[1],".las");
    input = Fopen(Catenate(APATH,"/",AROOT,".las"),"r");
    if (input == NULL)
      exit (1);

    if (fread(&novl,sizeof(int64),1,input) != 1)
      SYSTEM_READ_ERROR
    if (fread(&tspace,sizeof(int),1,input) != 1)
      SYSTEM_READ_ERROR
    if (tspace < 0)
      { fprintf(stderr,"%s: Garbage .las file, trace spacing < 0 !\n",Prog_Name);
        exit (1);
      }

    fseek(input,sizeof(int64)+sizeof(int),SEEK_SET);

    if (fread(&nlen,sizeof(int),1,input) != 1)
      SYSTEM_READ_ERROR

    fseek(input,nlen,SEEK_CUR);

    if (fread(&nlen,sizeof(int),1,input) != 1)
      SYSTEM_READ_ERROR

    if (nlen == 0)
      one = 1;
    else
      { one = 0;
        fseek(input,nlen,SEEK_CUR);
      }

    if (fread(&nlen,sizeof(int),1,input) != 1)
      SYSTEM_READ_ERROR
    fseek(input,nlen,SEEK_CUR);

    if (one != (argc == 3))
      { if (one)
          fprintf(stderr,"%s: Only 1 DB path is needed\n",Prog_Name);
        else
          fprintf(stderr,"%s: 2 DB paths are required\n",Prog_Name);
        exit (1);
      }
  }

  //  Open temporary output file buffer, write new header, and transfer data 

  { char *name, *buffer;
    int   nlen;

    TEMP   = Strdup(Numbered_Suffix(".lasreset.",getpid(),".las"),"Allocating temp name");
    output = Fopen(TEMP,"w");
    if (output == NULL)
      exit (1);

    if (fwrite(&novl,sizeof(int64),1,output) != 1)
      goto output_error;
    if (fwrite(&tspace,sizeof(int),1,output) != 1)
      goto output_error;

    name = Catenate(SPATH1,"/",SROOT1,".gdb");
    nlen = strlen(name);
    if (fwrite(&nlen,sizeof(int),1,output) != 1)
      goto output_error;
    if (fwrite(name,nlen,1,output) != 1)
      goto output_error;
  
    if (argc == 3)
      { nlen = 0;
        if (fwrite(&nlen,sizeof(int),1,output) != 1)
          goto output_error;
      }
    else
      { name = Catenate(SPATH2,"/",SROOT2,".gdb");
        nlen = strlen(name);
        if (fwrite(&nlen,sizeof(int),1,output) != 1)
          goto output_error;
        if (fwrite(name,nlen,1,output) != 1)
          goto output_error;
      }

    name = getcwd(NULL,0);
    nlen = strlen(name);
    if (fwrite(&nlen,sizeof(int),1,output) != 1)
      goto output_error;
    if (fwrite(name,nlen,1,output) != 1)
      goto output_error;
    free(name);

    buffer = Malloc(1000000,"Transfer buffer");
    if (buffer == NULL)
      exit (1);

    nlen = 1000000;
    while (nlen >= 1000000)
      { nlen = fread(buffer,1,1000000,input);
        fwrite(buffer,1,nlen,output);
      }

    if (rename(TEMP,Catenate(APATH,"/",AROOT,".las")) < 0)
      { fprintf(stderr,"%s: Could not complete the update, failure in rename\n",Prog_Name);
        unlink(TEMP);
        exit (1);
      }

    free(TEMP);
  }

  free(AROOT);
  free(APATH);
  if (argc == 4)
    { free(SROOT2);
      free(SPATH2);
    }
  free(SROOT1);
  free(SPATH1);
 
  exit (0);
    
output_error:
  fprintf(stderr,"%s: Could not write to %s\n",Prog_Name,TEMP);
  unlink(TEMP);
  exit (1);
}
