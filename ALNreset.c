/*******************************************************************************************
 *
 *  Utility for changing the database names in the header of a .1aln file
 *
 *  Author:    Gene Myers in a version for .las files
 *             Modified by Richard Durbin to use .1aln files
 *  Creation:  July 2013
 *  Last Mod:  March 2024
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
#include <pthread.h>

#include "DB.h"
#include "alncode.h"

static char *Usage = "[-T<int(8)>] <alignments:path>[.1aln] <DB1:path>[.gdb] [<DB2:path>[.gdb]]";

typedef struct
  { OneFile *in, *out;
    int64    beg, end;
  } Copy_Args;

static size_t fieldSize[128];

void *threadCopy(void *args)
{ Copy_Args *copy = (Copy_Args *) args;
  OneFile   *in  = copy->in;
  OneFile   *out = copy->out;
  int64      i, end;

  oneGotoObject(in,copy->beg);

  end = copy->end;
  for (i = copy->beg; i < end; i++)
    { oneReadLine(in);
      memcpy(out->field,in->field,fieldSize[(int) in->lineType]);
      oneWriteLine(out,in->lineType,oneLen(in),oneString(in));
    }
  
  return (NULL);
}

int main(int argc, char *argv[])
{ char *SPATH1, *SROOT1;
  char *SPATH2, *SROOT2;
  char *APATH, *AROOT;
  int   NTHREADS;
  char *command;

  { int   n, i;
    char *c;

    n = 0;
    for (i = 1; i < argc; i++)
      n += strlen(argv[i])+1;

    command = Malloc(n+1,"Allocating command string");
    if (command == NULL)
      exit (1);

    c = command;
    if (argc >= 1)
      { c += sprintf(c,"%s",argv[1]);
        for (i = 2; i < argc; i++)
          c += sprintf(c," %s",argv[i]);
      }
    *c = '\0';
  }
  
  //  Process options

  { int    i, j, k;
    int    flags[128];
    char*  eptr;

    ARG_INIT("ALNreset")

    (void) flags;

    NTHREADS = 8;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
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

  { FILE *input;

    SPATH1 = PathTo(argv[2]);
    SROOT1 = Root(argv[2],".gdb");
    input = Fopen(Catenate(SPATH1,"/",SROOT1,".gdb"),"r");
    if (input == NULL)
      exit (1);
    fclose(input);

    SPATH2 = NULL;
    SROOT2 = NULL;
    if (argc == 4)
      { SPATH2 = PathTo(argv[3]);
        SROOT2 = Root(argv[3],".gdb");
        input = Fopen(Catenate(SPATH2,"/",SROOT2,".gdb"),"r");
        if (input == NULL)
          exit (1);
        fclose(input);
      }
  }

  //  Open .1aln file reading and read header information
  
  { OneSchema *schema;
    char      *inFileName, *tmpFileName, *cpath;
    OneFile   *ofIn, *ofOut;
    int        i;

    APATH = PathTo(argv[1]);
    AROOT = Root(argv[1],".1aln");
    inFileName  = Strdup(Catenate(APATH,"/",AROOT,".1aln"),"Allocating input name");
    tmpFileName = Strdup(Numbered_Suffix(".lasreset.",getpid(),".1aln"),"Allocating temp name");
    if (inFileName == NULL || tmpFileName == NULL)
      exit (1);
    
    schema = oneSchemaCreateFromText(alnSchemaText);
    if (schema == NULL) 
      { fprintf(stderr,"%s: Failed to create 1aln schema\n",Prog_Name);
        exit (1);
      }
    ofIn = oneFileOpenRead(inFileName,schema,"aln",1);
    if (ofIn == NULL)
      { fprintf(stderr,"%s: Failed to open .1aln file %s\n",Prog_Name,inFileName);
        exit (1);
      }
    ofOut = oneFileOpenWriteFrom(tmpFileName,ofIn,true,1);
    if (ofOut == NULL)
      { fprintf(stderr,"%s: Failed to open .1aln file %s\n",Prog_Name,tmpFileName);
        exit (1);
      }
    oneAddProvenance(ofOut,Prog_Name,"0.1",command);

    oneAddReference (ofOut,Catenate(SPATH1,"/",SROOT1,".gdb"), 1);
    if (SROOT2 != NULL)
      oneAddReference (ofOut,Catenate(SPATH2,"/",SROOT1,".gdb"), 2);
    cpath = getcwd(NULL,0);
    oneAddReference(ofOut,cpath,3);
    free(cpath);

    if (ofIn->info['A'])  //  file not empty
      { int64      novl = ofIn->info['A']->given.count;
        Copy_Args  args[NTHREADS];
        pthread_t  threads[NTHREADS];
        
        for (i = 0; i < 128;++i)
          if (ofIn->info[i] != NULL)
            fieldSize[i] = ofIn->info[i]->nField*sizeof(OneField);

        for (i = 0; i < NTHREADS; i++)
          { args[i].in  = ofIn + i;
            args[i].out = ofOut + i;
            args[i].beg = (i * novl) / NTHREADS;
            if (i > 0)
              args[i-1].end = args[i].beg;
          }
        args[NTHREADS-1].end = novl;

        for (i = 1; i < NTHREADS; ++i)
          pthread_create(threads+i,NULL,threadCopy,args+i);
        threadCopy(args);
        for (i = 1; i < NTHREADS; i++)
          pthread_join(threads[i], NULL);
      }

    oneFileClose(ofIn);
    oneFileClose(ofOut);
    oneSchemaDestroy(schema);

    if (unlink(inFileName))
      { fprintf(stderr,"%s: Could not replace source\n",Prog_Name);
        unlink(tmpFileName);
        exit (1);
      }
    if (rename(tmpFileName,Catenate(APATH,"/",AROOT,".las")) < 0)
      { fprintf(stderr,"%s: Could mv the temp file %s to the source name\n",Prog_Name,tmpFileName);
        exit (1);
      }

    free(inFileName);
    free(tmpFileName);
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
}
