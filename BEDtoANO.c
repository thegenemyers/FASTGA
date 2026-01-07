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

#include "gene_core.h"
#include "GDB.h"
#include "ANO.h"
#include "hash.h"
#include "alncode.h"

#undef  DEBUG_THREADS

#define TSPACE   100
#define VERSION "0.1"

static char *Usage = "[-T<int(8)>] <bed:path>[.bed] [<genome:path>[.1gdb|<fa_extn>|<1_extn>]]";

static GDB_SCAFFOLD *SCAFF;

static Hash_Table *HASH;


/*******************************************************************************************
 *
 *  Read .bed file creating array of mask records and a GDB skeleton, GDB
 *    and hash tables, HASH, of its' scaffold names.
 *
 *******************************************************************************************/

typedef struct
  { int64    beg;
    int64    end;
    OneFile *of;
    char    *iname;
  } Packet;

//  Read next line into a buffer and return a pointer to the buffer
//    the length of the line.  NB: replaces '\n' with '\0'.

typedef struct
  { char *buffer;
    int   bmax;
  } Line_Bundle;

static char *read_line(void *input, char *name, Line_Bundle *lb)
{ char *buffer = lb->buffer;
  int   bmax   = lb->bmax;
  int len;

  if (bmax == 0)
    { bmax = 500;
      buffer = (char *) malloc(bmax);
      if (buffer == NULL)
        { fprintf(stderr,"%s: Out of memory reading %s\n",Prog_Name,name);
          exit (1);
        }
    }

  if (fgets(buffer,bmax,input) == NULL)
    { if (feof(input))
        { free(buffer);
          bmax = 0;
          return (NULL);
        }
      fprintf(stderr,"%s: Could not read next line of file %s (offset %lld)\n",
                     Prog_Name,name,(int64) ftello(input));
      exit (1);
    }

  len = strlen(buffer);
  while (buffer[len-1] != '\n')
    { bmax = ((int) (1.4*bmax)) + 100;
      buffer = (char *) realloc(buffer,bmax);
      if (buffer == NULL)
        { fprintf(stderr,"%s: Out of memory reading %s\n",Prog_Name,name);
          exit (1);
        }
      if (fgets(buffer+len,bmax-len,input) == NULL)
        { if (feof(input))
            fprintf(stderr,"%s: Last line of file %s does not end with new-line\n",
                           Prog_Name,name);
          else
            fprintf(stderr,"%s: Could not read next line of file %s (offset %lld)\n",
                           Prog_Name,name,(int64) ftello(input));
          exit (1);
        }
      len += strlen(buffer+len);
    }
  buffer[--len] = '\0';

  lb->bmax   = bmax;
  lb->buffer = buffer;

  return (buffer);
}

void *gen_1ano(void *args)
{ Packet  *parm  = (Packet *) args;
  char    *iname = parm->iname;
  OneFile *of    = parm->of;   

  Line_Bundle    _lb, *lb = &_lb;

  int64  span, amnt;
  FILE  *input;
  int    nfields, linelen;
  char  *fptrs[100], *eptr;

  int64  beg, end, mbuffer[2], score;
  int    idx;

  input = fopen(iname,"r");
  if (input == NULL)
    { fprintf(stderr,"%s: Could not open %s for reading\n",Prog_Name,iname);
      exit (1);
    }

  fseeko(input,parm->beg,SEEK_SET);
  span = parm->end-parm->beg;

  lb->bmax = 0;
  if (parm->beg > 0)
    amnt = strlen(read_line(input,iname,lb)) + 1;
  else
    amnt = 0;

  while (amnt <= span)
    { eptr = read_line(input,iname,lb);
      if (eptr == NULL)
        break;
      linelen = strlen(eptr) + 1;
      amnt   += linelen;
    
      while (*eptr != '\0' && isspace(*eptr))    // parse BED line
        eptr += 1;

      if (*eptr == '#' || strncmp(eptr,"track:",6) == 0 || strncmp(eptr,"browser:",8) == 0)
        continue;

      nfields = 0;
      while (*eptr != '\0')
        { fptrs[nfields++] = eptr;
          while (*eptr != '\0' && !isspace(*eptr))
            eptr += 1;
          *eptr++ = '\0';
          while (isspace(*eptr))
            eptr += 1;
        }
 
      if (nfields < 3)
        { fprintf(stderr,"%s: Line of bed has fewer than 3 fields (offset %lld)\n",
                         Prog_Name,(int64) (ftello(input)-linelen));
          exit (1);
        }

      idx = Hash_Lookup(HASH, fptrs[0]);
      if (idx < 0)
        { fprintf(stderr,"%s: Scaffold name %s not found in first source (offset %lld)\n", 
                         Prog_Name,fptrs[0],(int64) (ftello(input)-linelen));
          exit (1);
        }
      beg = strtol(fptrs[1], &eptr, 10);
      end = strtol(fptrs[2], &eptr, 10);
      if (beg > end)
        { fprintf(stderr,"%s: Begin coord %lld > End coord %lld (offset %lld)\n",
                         Prog_Name,beg,end,(int64) (ftello(input)-linelen));
          exit (1);
        }
      if (beg < 0)
        { fprintf(stderr,"%s: Begin coord %lld less than 0 (offset %lld)\n",
                         Prog_Name,beg,(int64) (ftello(input)-linelen));
          exit (1);
        }
      if (end > SCAFF[idx].slen)
        { fprintf(stderr,"%s: End coord %lld greater than scaffold length (offset %lld)\n",
                         Prog_Name,end,(int64) (ftello(input)-linelen));
          exit (1);
        }
      if (nfields >= 5)
        score = strtol(fptrs[5], &eptr, 10);
      else
        score = 0;
      if (nfields >= 6)
        { if (strcmp(fptrs[5],"-") == 0)
            { int64 x = beg; beg = end; end = x; } 
        }

      oneInt(of,0) = idx;
      mbuffer[0] = beg;
      mbuffer[1] = end;
      oneWriteLine(of,'M',2,mbuffer);

      if (nfields >= 4)
        oneWriteLine(of,'L',strlen(fptrs[3]),fptrs[3]);

      if (score > 0)
        oneWriteLine(of,'X',1,&score);
    }

  fclose(input);

  return (NULL);
}


int main(int argc, char *argv[])
{ Packet    *parm;
  GDB       _gdb, *gdb = &_gdb;
  char      *SPATH;
  char      *APATH, *AROOT;

  int   NTHREADS;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char*  eptr;

    ARG_INIT("BEDtoANO")

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

    if (argc != 3)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        exit (1);
      }

    parm = Malloc(sizeof(Packet)*NTHREADS,"Allocating thread records");
    if (parm == NULL)
      exit (1);
  }

  //  Find sources

  { char *tpath;
    int   type;

    type = Get_GDB_Paths(argv[2],NULL,&SPATH,&tpath,0);
    if (type != IS_GDB)
      Create_GDB(gdb,SPATH,type,0,NULL,0,NULL);
    else
      Read_GDB(gdb,tpath);
    free(tpath);
  }

  //  Set up scaffold->contig maps & scaffold name dictionary

  { int   s;
    char *head, *sptr, *eptr;

    SCAFF = gdb->scaffolds;
    HASH  = New_Hash_Table(gdb->nscaff,0);
    head  = gdb->headers;
    for (s = 0; s < gdb->nscaff; s++)
      { sptr = head + SCAFF[s].hoff;
        for (eptr = sptr; *eptr != '\0'; eptr++)
          if (isspace(*eptr))
            break;
        *eptr = '\0';
        if (Hash_Lookup(HASH,sptr) < 0)
          Hash_Add(HASH,sptr);
        else
          { fprintf(stderr,"%s: Duplicate scaffold name: %s\n",Prog_Name,sptr);
            exit (1);
          }
      }
  }

  //  Open .1aln file reading and read header information
  
  { char       *inName;
    OneSchema  *schema;
    OneFile    *of;
    struct stat state;
#ifndef DEBUG_THREADS
    pthread_t   threads[NTHREADS];
#endif
    int         p;

    APATH  = PathTo(argv[1]);
    AROOT  = Root(argv[1],".bed");
    inName = Strdup(Catenate(APATH,"/",AROOT,".bed"),"Allocating input name");
    if (inName == NULL)
      exit (1);

    if (stat(inName,&state) < 0)
      { fprintf(stderr,"%s: Cannot open %s to get size\n",Prog_Name,inName);
        exit (1);
      }
      
    schema = make_ANO_Schema();
    if (schema == NULL)
      { fprintf(stderr,"%s: Failed to create ANO schema\n",Prog_Name);
        exit (1);
      }
  
    of = oneFileOpenWriteNew(Catenate(APATH,"/",AROOT,".1ano"),schema,"ano",1,NTHREADS);
    if (of == NULL)
      { fprintf(stderr,"%s: Failed to open file %s/%s.1ano\n",Prog_Name,APATH,AROOT);
        oneSchemaDestroy(schema);
        exit (1);
      } 
    
    oneAddProvenance(of,Prog_Name,"0.1",Command_Line);
      
    oneAddReference(of,gdb->srcpath,1);

    Write_Skeleton(of,gdb);

    for (p = 0; p < NTHREADS; p++)
      { parm[p].beg   = (state.st_size * p) / NTHREADS;
        parm[p].end   = (state.st_size * (p+1)) / NTHREADS;
        parm[p].iname = inName;
        parm[p].of    = of+p;
      }

    //  Launch and then gather threads

#ifdef DEBUG_THREADS
    for (p = 0; p < NTHREADS; p++)
      gen_1ano(parm+p);
#else
    for (p = 1; p < NTHREADS; p++)
      pthread_create(threads+p,NULL,gen_1ano,parm+p);
    gen_1ano(parm);
    for (p = 1; p < NTHREADS; p++)
      pthread_join(threads[p],NULL);
#endif

    oneFileClose(of);
    oneSchemaDestroy(schema);
    free(inName);
    free(AROOT);
    free(APATH);
  }

  Close_GDB(gdb);
 
  exit (0);
}
