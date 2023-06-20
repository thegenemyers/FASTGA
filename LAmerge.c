/*******************************************************************************************
 *
 *  Given a list of sorted .las files, merge them into a single sorted .las file.
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "DB.h"
#include "align.h"

#undef   DEBUG

static char *Usage = "[-va] [-P<dir(/tmp)>] <merge:las> <parts:las> ...";

#define MEMORY 4000   // in Mb

#define MAX_FILES 250

  //  Heap sort of records according to (aread,bread,COMP(flags),abpos) order

#define COMPARE(lp,rp)				\
  if (lp->aread > rp->aread)			\
    bigger = 1;					\
  else if (lp->aread < rp->aread)		\
    bigger = 0;					\
  else if (lp->bread > rp->bread)		\
    bigger = 1;					\
  else if (lp->bread < rp->bread)		\
    bigger = 0;					\
  else if (COMP(lp->flags) > COMP(rp->flags))	\
    bigger = 1;					\
  else if (COMP(lp->flags) < COMP(rp->flags))	\
    bigger = 0;					\
  else if (lp->path.abpos > rp->path.abpos)	\
    bigger = 1;					\
  else if (lp->path.abpos < rp->path.abpos)	\
    bigger = 0;					\
  else if (lp > rp)				\
    bigger = 1;					\
  else						\
    bigger = 0;

static void reheap(int s, Overlap **heap, int hsize)
{ int      c, l, r;
  int      bigger;
  Overlap *hs, *hr, *hl;

  c  = s;
  hs = heap[s];
  while ((l = 2*c) <= hsize)
    { r  = l+1;
      hl = heap[l];
      if (r > hsize)
        bigger = 1;
      else
        { hr = heap[r];
          COMPARE(hr,hl)
        }
      if (bigger)
        { COMPARE(hs,hl)
          if (bigger)
            { heap[c] = hl;
              c = l;
            }
          else
            break;
        }
      else
        { COMPARE(hs,hr)
          if (bigger)
            { heap[c] = hr;
              c = r;
            }
          else
            break;
        }
    }
  if (c != s)
    heap[c] = hs;
}

  //  Heap sort of records according to (aread,abpos) order

#define MAPARE(lp,rp)				\
  if (lp->aread > rp->aread)			\
    bigger = 1;					\
  else if (lp->aread < rp->aread)		\
    bigger = 0;					\
  else if (lp->path.abpos > rp->path.abpos)	\
    bigger = 1;					\
  else if (lp->path.abpos < rp->path.abpos)	\
    bigger = 0;					\
  else if (lp > rp)				\
    bigger = 1;					\
  else						\
    bigger = 0;

static void maheap(int s, Overlap **heap, int hsize)
{ int      c, l, r;
  int      bigger;
  Overlap *hs, *hr, *hl;

  c  = s;
  hs = heap[s];
  while ((l = 2*c) <= hsize)
    { r  = l+1;
      hl = heap[l];
      if (r > hsize)
        bigger = 1;
      else
        { hr = heap[r];
          MAPARE(hr,hl)
        }
      if (bigger)
        { MAPARE(hs,hl)
          if (bigger)
            { heap[c] = hl;
              c = l;
            }
          else
            break;
        }
      else
        { MAPARE(hs,hr)
          if (bigger)
            { heap[c] = hr;
              c = r;
            }
          else
            break;
        }
    }
  if (c != s)
    heap[c] = hs;
}

#ifdef DEBUG

static void showheap(Overlap **heap, int hsize)
{ int i;
  printf("\n");
  for (i = 1; i <= hsize; i++)
    printf(" %3d: %5d, %5d\n",i,heap[i]->aread,heap[i]->bread);
}

#endif

  //  Input block data structure and block fetcher

typedef struct
  { FILE   *stream;
    char   *block;
    char   *ptr;
    char   *top;
    int64   count;
  } IO_block;

static void ovl_reload(IO_block *in, int64 bsize)
{ int64 remains;

  remains = in->top - in->ptr;
  if (remains > 0)
    memmove(in->block, in->ptr, remains);
  in->ptr  = in->block;
  in->top  = in->block + remains;
  in->top += fread(in->top,1,bsize-remains,in->stream);
}

  //  The program

int main(int argc, char *argv[])
{ IO_block *in;
  int64     bsize, osize, psize;
  char     *block, *oblock;
  int       i, c, fway, clen, nfile[argc];
  Overlap **heap;
  int       hsize;
  Overlap  *ovls;
  int64     totl;
  int       tspace, tbytes;
  FILE     *output;
  char     *optr, *otop;

  int       VERBOSE;
  int       MAP_SORT;
  char     *TEMP_PATH;

  //  Process command line

  { int  j, k;
    int  flags[128];
    DIR *dirp;

    ARG_INIT("LAmerge")

    TEMP_PATH = "/tmp";

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("va")
            break;
          case 'P':
            TEMP_PATH = argv[i]+2;
            if ((dirp = opendir(TEMP_PATH)) == NULL)
              { fprintf(stderr,"%s: -P option: cannot open directory %s\n",Prog_Name,TEMP_PATH);
                exit (1);
              }
            closedir(dirp);
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE  = flags['v'];
    MAP_SORT = flags['a'];

    if (argc < 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -a: sort .las by A-read,A-position pairs for map usecase\n");
        fprintf(stderr,"          off => sort .las by A,B-read pairs for overlap piles\n");
        fprintf(stderr,"      -P: Do any intermediate merging in directory -P.\n");
        exit (1);
      }
  }

  //  Determine the number of files and check they are all mergeable

  clen   = 2*strlen(TEMP_PATH) + 50;
  fway   = 0;
  totl   = 0;
  tspace = -1;
  for (c = 2; c < argc; c++)
    { Block_Looper *parse;
      FILE *input;
      char *root, *path;

      parse = Parse_Block_LAS_Arg(argv[c]);

      path = Block_Arg_Path(parse);
      root = Block_Arg_Root(parse);

      clen += strlen(path) + strlen(root) + 30;

      free(root);
      free(path);

      nfile[c] = 0;
      while ((input = Next_Block_Arg(parse)) != NULL)
        { int64 povl;
          int   mspace;

          if (fread(&povl,sizeof(int64),1,input) != 1)
            SYSTEM_READ_ERROR
          totl += povl;
          if (fread(&mspace,sizeof(int),1,input) != 1)
            SYSTEM_READ_ERROR
          if (tspace < 0)
            tspace = mspace;
          else if (tspace != mspace)
            { fprintf(stderr,"%s: trace-point spacing conflict between %s and earlier files",
                             Prog_Name,Block_Arg_Root(parse));
              fprintf(stderr," (%d vs %d)\n",tspace,mspace);
              exit (1);
            }

          fclose(input);
          nfile[c] += 1;
        }

      Free_Block_Arg(parse);
      fway += nfile[c];
    }

  if (VERBOSE)
    { printf("  Merging %d files totaling ",fway);
      Print_Number(totl,0,stdout);
      printf(" records\n");
      fflush(stdout);
    }

  //  Must recursively merge, emit sub-merges, then merge their results

  if (fway > MAX_FILES)
    { Block_Looper *parse;
      int   mul, dim, fsum, cut;
      char  command[clen], *com;
      int   pid;

      mul = 1;
      for (c = 0; mul < fway; c++)
        mul *= MAX_FILES;
      dim = pow(1.*fway,1./c)+1;

      fsum = 0;
      c = 2;

      parse = Parse_Block_LAS_Arg(argv[c]);

      pid = getpid();
      for (i = 1; i <= dim; i++)
        { com = command;
          com += sprintf(com,"LAmerge");
          if (MAP_SORT)
            com += sprintf(com," -a");
          if (mul > 2)
            com += sprintf(com," -P%s",TEMP_PATH);
          com += sprintf(com," %s/LM%d.P%d",TEMP_PATH,pid,i);

          cut = (fway * i) / dim;
          while (fsum + nfile[c] <= cut)
            { com  += sprintf(com," %s",Next_Block_Slice(parse,nfile[c]));
              fsum += nfile[c];

              c += 1;
              if (c >= argc)
                break;

              Free_Block_Arg(parse);

              parse = Parse_Block_LAS_Arg(argv[c]);
            }
          if (c < argc && fsum < cut)
            { int n = cut-fsum;
              com += sprintf(com," %s",Next_Block_Slice(parse,n));
              nfile[c] -= n;
              fsum     += n;
            }
          system(command);
        }

      Free_Block_Arg(parse);

      com = command;
      com += sprintf(com,"LAmerge");
      if (MAP_SORT)
        com += sprintf(com," -a");
      com += sprintf(com," %s %s/LM%d.P%c",argv[1],TEMP_PATH,pid,BLOCK_SYMBOL);
      system(command);

      sprintf(command,"rm %s/LM%d.P*.las",TEMP_PATH,pid);
      system(command);

      exit (0);
    }

  //  Base level merge: Open all the input files and initialize their buffers

  psize  = sizeof(void *);
  osize  = sizeof(Overlap) - psize;
  bsize  = (MEMORY*1000000ll)/(fway + 1);
  block  = (char *) Malloc(bsize*(fway+1)+psize,"Allocating LAmerge blocks");
  in     = (IO_block *) Malloc(sizeof(IO_block)*fway,"Allocating LAmerge IO-reacords");
  if (block == NULL || in == NULL)
    exit (1);
  block += psize;

  fway = 0;
  for (c = 2; c < argc; c++)
    { Block_Looper *parse;
      FILE  *input;

      parse = Parse_Block_LAS_Arg(argv[c]);

      while ((input = Next_Block_Arg(parse)) != NULL)
        { int64  novl;
          int    mspace;
          char  *iblock;

          if (fread(&novl,sizeof(int64),1,input) != 1)
            SYSTEM_READ_ERROR
          if (fread(&mspace,sizeof(int),1,input) != 1)
            SYSTEM_READ_ERROR

          in[fway].stream = input;
          in[fway].block  = iblock = block+fway*bsize;
          in[fway].ptr    = iblock;
          in[fway].top    = iblock + fread(in[fway].block,1,bsize,input);
          in[fway].count  = 0;
          fway += 1;
        }

      Free_Block_Arg(parse);
    }
  if (tspace <= TRACE_XOVR && tspace != 0)
    tbytes = sizeof(uint8);
  else
    tbytes = sizeof(uint16);

  //  Open the output file buffer and write (novl,tspace) header

  { char *pwd, *root;

    pwd    = PathTo(argv[1]);
    root   = Root(argv[1],".las");
    output = Fopen(Catenate(pwd,"/",root,".las"),"w");
    if (output == NULL)
      exit (1);
    free(pwd);
    free(root);

    if (fwrite(&totl,sizeof(int64),1,output) != 1)
      SYSTEM_READ_ERROR
    if (fwrite(&tspace,sizeof(int),1,output) != 1)
      SYSTEM_READ_ERROR

    oblock = block+fway*bsize;
    optr   = oblock;
    otop   = oblock + bsize;
  }

  //  Initialize the heap

  heap = (Overlap **) Malloc(sizeof(Overlap *)*(fway+1),"Allocating heap");
  ovls = (Overlap *) Malloc(sizeof(Overlap)*fway,"Allocating heap");
  if (heap == NULL || ovls == NULL)
    exit (1);

  hsize = 0;
  for (i = 0; i < fway; i++)
    { if (in[i].ptr < in[i].top)
        { ovls[i]     = *((Overlap *) (in[i].ptr - psize));
          in[i].ptr  += osize;
          hsize      += 1;
          heap[hsize] = ovls + i;
        }
    }

  if (hsize > 3)
    { if (MAP_SORT)
        for (i = hsize/2; i > 1; i--)
          maheap(i,heap,hsize);
      else
        for (i = hsize/2; i > 1; i--)
          reheap(i,heap,hsize);
    }

  //  While the heap is not empty do

  while (hsize > 0)
    { Overlap  *ov;
      IO_block *src;
      int64     tsize, span;

      if (MAP_SORT)
        maheap(1,heap,hsize);
      else
        reheap(1,heap,hsize);

      ov  = heap[1];
      src = in + (ov - ovls);

      do
        { src->count += 1;

          tsize = ov->path.tlen*tbytes;
          span  = osize + tsize;
          if (src->ptr + span > src->top)
            ovl_reload(src,bsize);
          if (optr + span > otop)
            { if (fwrite(oblock,1,optr-oblock,output) != (size_t) (optr-oblock))
                SYSTEM_READ_ERROR
              optr = oblock;
            }

          memmove(optr,((char *) ov) + psize,osize);
          optr += osize;
          memmove(optr,src->ptr,tsize);
          optr += tsize;

          src->ptr += tsize;
          if (src->ptr >= src->top)
            { heap[1] = heap[hsize];
              hsize  -= 1;
              break;
            }
          *ov       = *((Overlap *) (src->ptr - psize));
          src->ptr += osize;
        }
      while (CHAIN_NEXT(ov->flags));
    }

  //  Flush output buffer and wind up

  if (optr > oblock)
    { if (fwrite(oblock,1,optr-oblock,output) != (size_t) (optr-oblock))
        SYSTEM_READ_ERROR
    }
  fclose(output);

  for (i = 0; i < fway; i++)
    fclose(in[i].stream);

  for (i = 0; i < fway; i++)
    totl -= in[i].count;
  if (totl != 0)
    { fprintf(stderr,"%s: Did not write all records to %s (%lld)\n",argv[0],argv[1],totl);
      exit (1);
    }

  free(ovls);
  free(heap);
  free(in);
  free(block-psize);

  exit (0);
}
