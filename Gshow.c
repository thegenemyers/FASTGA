#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include "libfastk.h"
#include "DB.h"

#undef START_AT_X

static char *Usage = " <source>[.dam]";

/***********************************************************************************************
 *
 *   POSITION LIST ABSTRACTION:
 *       Routines to manipuate the position or "post" list associated with a k-mer table.
 *
 **********************************************************************************************/

typedef struct
  { int     pbyte;      //  # of bytes for each post (including sign & contig)
    int     cbyte;      //  # of bytes for each sign & contig
    int64   nels;       //  # of posts in the index
    int64   maxp;
    int     freq;
    int     nctg;
    int    *perm;
    int64   cidx;
    uint8  *cache;
    uint8  *cptr;

    int     copn;       //  File currently open
    int     part;       //  Thread # of file currently open
    int     nthr;       //  # of file parts
    int     nsqrt;      //  # of threads/slices (= sqrt(nthr))
    int     nlen;       //  length of path name
    char   *name;       //  Path name for table parts (only # missing)
    uint8  *ctop;       //  Ptr top of current table block in buffer
    int64  *neps;       //  Size of each thread part in elements
  } Post_List;

#define POST_BLOCK 1024

//  Load up the table buffer with the next STREAM_BLOCK suffixes (if possible)

static void More_Post_List(Post_List *P)
{ int    pbyte = P->pbyte;
  uint8 *cache = P->cache;
  int    copn  = P->copn;
  uint8 *ctop;

  if (P->part > P->nthr)
    return;
  while (1)
    { ctop = cache + read(copn,cache,POST_BLOCK*pbyte);
      if (ctop > cache)
        break;
      close(copn);
      P->part += 1;
      if (P->part > P->nthr)
        { P->cptr = NULL;
          return;
        }
      sprintf(P->name+P->nlen,"%d",P->part);
      copn = open(P->name,O_RDONLY);
      lseek(copn,2*sizeof(int)+sizeof(int64),SEEK_SET);
    }
  P->cptr = cache;
  P->ctop = ctop;
  P->copn = copn;
}

static Post_List *Open_Post_List(char *name)
{ Post_List *P;
  int        pbyte, cbyte, nctg;
  int64      nels, maxp, n;
  int        copn;

  int    f, p, flen;
  char  *dir, *root, *full;
  int    pb, cb, nfile, nthreads, freq;

  dir  = PathTo(name);
  root = Root(name,".ktab");
  full = Malloc(strlen(dir)+strlen(root)+20,"Post list name allocation");
  if (full == NULL)
    exit (1);
  sprintf(full,"%s/%s.post",dir,root);
  f = open(full,O_RDONLY);
  sprintf(full,"%s/.%s.post.",dir,root);
  flen = strlen(full);
  free(root);
  free(dir);
  if (f < 0)
    { free(full);
      return (NULL);
    }

  read(f,&pbyte,sizeof(int));
  read(f,&cbyte,sizeof(int));
  pbyte += cbyte;

  read(f,&nfile,sizeof(int));
  read(f,&maxp,sizeof(int64));
  read(f,&freq,sizeof(int));
  nthreads = nfile;
  nfile    = nfile*nfile;
  
  read(f,&nctg,sizeof(int));

  P = Malloc(sizeof(Post_List),"Allocating post record");
  if (P == NULL)
    exit (1);
  P->name   = full;
  P->nlen   = strlen(full);
  P->maxp   = maxp;
  P->cache  = Malloc(POST_BLOCK*pbyte,"Allocating post list buffer\n");
  P->neps   = Malloc(nfile*sizeof(int64),"Allocating parts table of Post_List");
  P->perm   = Malloc(nctg*sizeof(int),"Allocating sort permutation");
  if (P->cache == NULL || P->neps == NULL || P->perm == NULL)
    exit (1);

  read(f,P->perm,sizeof(int)*nctg);

  nels = 0;
  for (p = 1; p <= nfile; p++)
    { sprintf(P->name+P->nlen,"%d",p);
      copn = open(P->name,O_RDONLY);
      if (copn < 0)
        { fprintf(stderr,"%s: Table part %s is missing ?\n",Prog_Name,P->name);
          exit (1);
        }
      read(copn,&pb,sizeof(int));
      read(copn,&cb,sizeof(int));
      pb += cb;
      read(copn,&n,sizeof(int64));
      nels += n;
      P->neps[p-1] = nels;
      if (pbyte != pb)
        { fprintf(stderr,"%s: Post list part %s does not have post size matching stub ?\n",
                         Prog_Name,P->name);
          exit (1);
        }
      close(copn);
    }

  P->pbyte = pbyte;
  P->cbyte = cbyte;
  P->nels  = nels;
  P->nthr  = nfile;
  P->nsqrt = nthreads;
  P->freq  = freq;
  P->nctg  = nctg;

  sprintf(P->name+P->nlen,"%d",1);
  copn = open(P->name,O_RDONLY);
  lseek(copn,2*sizeof(int)+sizeof(int64),SEEK_SET);

  P->copn  = copn;
  P->part  = 1;

  More_Post_List(P);
  P->cidx = 0;

  return (P);
}

static void Free_Post_List(Post_List *P)
{ free(P->neps);
  free(P->cache);
  free(P->perm);
  if (P->copn >= 0)
    close(P->copn);
  free(P->name);
  free(P);
}

static inline void First_Post_Entry(Post_List *P)
{ if (P->cidx != 0)
    { if (P->part != 1)
        { if (P->part <= P->nthr)
            close(P->copn);
          sprintf(P->name+P->nlen,"%d",1);
          P->copn = open(P->name,O_RDONLY);
          P->part = 1;
        }

      lseek(P->copn,sizeof(int)+sizeof(int64),SEEK_SET);

      More_Post_List(P);
      P->cidx = 0;
    }
}

static inline void Next_Post_Entry(Post_List *P)
{ P->cptr += P->pbyte;
  P->cidx += 1;
  if (P->cptr >= P->ctop)
    { if (P->cidx >= P->nels)
        { P->cptr = NULL;
          P->part = P->nthr+1;
          return;
        }
      More_Post_List(P);
    }
}

static inline int64 Current_Post(Post_List *P)
{ int64  post;

  post = 0;
  memcpy((uint8 *) (&post),P->cptr,P->pbyte-P->cbyte);
  return (post);
}

static inline int64 Current_Contig(Post_List *P)
{ int64  cont;

  cont = 0;
  memcpy((uint8 *) (&cont),P->cptr+(P->pbyte-P->cbyte),P->cbyte);
  return (cont);
}

static inline int Current_Sign(Post_List *P)
{ return ((P->cptr[P->pbyte-1] & 0x80) != 0); }


#ifdef START_AT_X

static inline void JumpTo_Post_Index(Post_List *P, int64 del)
{ int   p;
  int64 i;

  P->cptr += del*P->pbyte;
  P->cidx += del;
  if (P->cptr < P->ctop)
    return;

  i = P->cidx;
  p = P->part-1;
  while (i >= P->neps[p])
    p += 1;

  if (p > 0)
    i -= P->neps[p-1];
  p += 1;

  if (P->part != p)
    { if (P->part <= P->nthr)
        close(P->copn);
      sprintf(P->name+P->nlen,"%d",p);
      P->copn = open(P->name,O_RDONLY);
      P->part = p;
    }

  lseek(P->copn,2*sizeof(int) + sizeof(int64) + i*P->pbyte,SEEK_SET);

  More_Post_List(P);
}

#endif


/***********************************************************************************************
 *
 *   LIST INDEX
 *
 **********************************************************************************************/

static void Print_Index(Kmer_Stream *T, Post_List *P)
{ char  *buffer;
  int    invert;
  int64  post, cont, flag;
  int    count, cbyte, pbyte, p;
  uint8 *cptr  = (uint8 *) (&cont);
  uint8 *pptr  = (uint8 *) (&post);
  uint8 *split = (uint8 *) (&count);
  int   *perm  = P->perm;

  (void) Current_Post;   //  directly code memcpy below
  (void) Current_Contig;
  (void) Current_Sign;

  buffer = Current_Kmer(T,NULL);
  cbyte  = P->cbyte;
  pbyte  = P->pbyte - cbyte;
  flag   = (0x1ll << (8*cbyte-1));

  post = cont = 0;

#ifdef START_AT_X
  int64 bidx;

  bidx = 0;
  First_Post_Entry(P);
  for (First_Kmer_Entry(T); T->cidx < X; Next_Kmer_Entry(T))
    { count = Current_Count(T);
      bidx += split[0];
    }
  JumpTo_Post_Index(P,bidx);

  for (; T->cidx < T->nels; Next_Kmer_Entry(T))
#else
  First_Post_Entry(P);
  for (First_Kmer_Entry(T); T->cidx < T->nels; Next_Kmer_Entry(T))
#endif
    { count = Current_Count(T);
      printf(" %6lld: %s %2d\n",T->cidx,Current_Kmer(T,buffer),split[1]);
      fflush(stdout);
      for (p = split[0]; p > 0; p--)
        { memcpy(pptr,P->cptr,pbyte);
          memcpy(cptr,P->cptr+pbyte,cbyte);
          invert = ((cont & flag) != 0);
          if (invert)
            cont -= flag;
          printf("    %11lld:   %c %4d:%10lld\n",P->cidx,invert?'-':'+',perm[cont],post);
          fflush(stdout);
          Next_Post_Entry(P);
        }
    }
}

int main(int argc, char *argv[])
{ Kmer_Stream *T;
  Post_List   *P;

  //  Process options

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("Gshow");

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

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  T = Open_Kmer_Stream(argv[1]);
  P = Open_Post_List(argv[1]);

  Print_Index(T,P); 

  Free_Post_List(P);
  Free_Kmer_Stream(T);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
