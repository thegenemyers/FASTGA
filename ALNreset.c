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

#include "GDB.h"
#include "alncode.h"

static char *Usage[] = 
            { "[-T<int(8)>] <alignments:path>[.1aln]",
              " <source1:path>[.1gdb|<fa_extn>|<1_extn>] [<source2:path>[.1gdb|<fa_extn>|<1_extn>]]"
            };

typedef struct
  { OneFile *in, *out;
    int64    beg, end;
  } Copy_Args;

static size_t fieldSize[128];

void *threadCopy(void *args)
{ Copy_Args *copy = (Copy_Args *) args;
  OneFile   *in  = copy->in;
  OneFile   *out = copy->out;
  int64      i, n;

  oneGoto(in,'A',copy->beg+1);

  n = copy->end - copy->beg;
  i = 0;
  for (i = 1; oneReadLine(in); i++)
    { if (in->lineType == 'A' && i > n)
        break;
      memcpy(out->field,in->field,fieldSize[(int) in->lineType]);
      oneWriteLine(out,in->lineType,oneLen(in),oneString(in));
    }
  
  return (NULL);
}

int main(int argc, char *argv[])
{ char *SPATH1;
  char *SPATH2;
  char *APATH, *AROOT;
  int   NTHREADS;

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
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        exit (1);
      }
  }

  //  Find sources

  { char *tpath;

    Get_GDB_Paths(argv[2],NULL,&SPATH1,&tpath,0);
    free(tpath);
    if (argc == 4)
      { Get_GDB_Paths(argv[3],NULL,&SPATH2,&tpath,0);
        free(tpath);
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
    
    schema = make_Aln_Schema();
    if (schema == NULL) 
      { fprintf(stderr,"%s: Failed to create 1aln schema\n",Prog_Name);
        exit (1);
      }
    ofIn = oneFileOpenRead(inFileName,schema,"aln",NTHREADS);
    if (ofIn == NULL)
      { fprintf(stderr,"%s: Failed to open .1aln file %s\n",Prog_Name,inFileName);
        exit (1);
      }
    ofOut = oneFileOpenWriteFrom(tmpFileName,ofIn,true,NTHREADS);
    if (ofOut == NULL)
      { fprintf(stderr,"%s: Failed to open .1aln file %s\n",Prog_Name,tmpFileName);
        exit (1);
      }
    oneAddProvenance(ofOut,Prog_Name,"0.1",Command_Line);

    oneAddReference (ofOut,SPATH1, 1);
    if (argc == 4)
      oneAddReference (ofOut,SPATH2, 2);
    cpath = getcwd(NULL,0);
    oneAddReference(ofOut,cpath,3);
    free(cpath);

    while (oneReadLine(ofIn))         // Transfer any pre-object lines
      { if (ofIn->lineType == 'A')
          break;
        memcpy(ofOut->field,ofIn->field,fieldSize[(int) ofIn->lineType]);
        oneWriteLine(ofOut,ofIn->lineType,oneLen(ofIn),oneString(ofIn));
      }

    if (ofIn->info['A'])  //  file not empty
      { int64      novl = ofIn->info['A']->given.count;
        Copy_Args  args[NTHREADS];
        pthread_t  threads[NTHREADS];

        for (i = 0; i < 128;++i)
          if (ofIn->info[i] != NULL)
            fieldSize[i] = ofIn->info[i]->nField*sizeof(OneField);

        //  Write global lines before 'A'
        while (oneReadLine(ofIn) && ofIn->lineType != 'A')
          { memcpy(ofOut->field,ofIn->field,fieldSize[(int) ofIn->lineType]);
            oneWriteLine(ofOut,ofIn->lineType,oneLen(ofIn),oneString(ofIn));
          }

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
    if (rename(tmpFileName,Catenate(APATH,"/",AROOT,".1aln")) < 0)
      { fprintf(stderr,"%s: Could mv the temp file %s to the source name\n",Prog_Name,tmpFileName);
        exit (1);
      }

    free(inFileName);
    free(tmpFileName);
  }

  free(AROOT);
  free(APATH);
  if (argc == 4)
    free(SPATH2);
  free(SPATH1);
 
  exit (0);
}
