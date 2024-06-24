#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include "libfastk.h"
#include "GDB.h"

static char *Usage = "<source>[.gix] [ <address>[-<address>] ] ";

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
    int64  *index;

    int     copn;       //  File currently open
    int     part;       //  Thread # of file currently open
    int     nthr;       //  # of file parts
    int     nlen;       //  length of path name
    char   *name;       //  Path name for table parts (only # missing)
    uint8  *ctop;       //  Ptr top of current table block in buffer
    int64  *neps;       //  Size of each thread part in elements
  } Post_List;

#define POST_BLOCK 0x20000

//  Load up the table buffer with the next STREAM_BLOCK suffixes (if possible)

static void More_Post_List(Post_List *P)
{ int    pbyte = P->pbyte;
  uint8 *cache = P->cache;
  int    copn  = P->copn;
  int    len;
  uint8 *ctop;

  if (P->part > P->nthr)
    return;
  while (1)
    { len  = read(copn,cache,POST_BLOCK*pbyte);
      if (len < 0)
        { fprintf(stderr,"%s: Error reading post file %s\n",Prog_Name,P->name);
          exit (1);
        }
      ctop = cache + len;
      if (len > 0)
        break;
      close(copn);
      P->part += 1;
      if (P->part > P->nthr)
        { P->cptr = NULL;
          return;
        }
      sprintf(P->name+P->nlen,"%d",P->part);
      copn = open(P->name,O_RDONLY);
      if (copn < 0)
        { fprintf(stderr,"%s: Cannot open post file %s for reading\n",Prog_Name,P->name);
          exit (1);
        }
      if (lseek(copn,2*sizeof(int)+sizeof(int64),SEEK_SET) < 0)
        { fprintf(stderr,"%s: Cannot advance post file %s to data part\n",Prog_Name,P->name);
          exit (1);
        }
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

  int    f, p;
  char  *dir, *root, *full;
  int    pb, cb, nfile, freq;

  dir  = PathTo(name);
  root = Root(name,".gix");
  full = Malloc(strlen(dir)+strlen(root)+20,"Post list name allocation");
  if (full == NULL)
    exit (1);
  sprintf(full,"%s/%s.gix",dir,root);
  f = open(full,O_RDONLY);
  if (f < 0)
    { fprintf(stderr,"%s: Cannot open post stub file %s/%s.post\n",Prog_Name,dir,root);
      exit (1);
    }
  sprintf(full,"%s/.%s.post.",dir,root);
  free(root);
  free(dir);

  if (lseek(f,4*sizeof(int)+0x1000000*sizeof(int64),SEEK_SET) < 0) goto open_io_error;

  if (read(f,&pbyte,sizeof(int)) < 0) goto open_io_error;
  if (read(f,&cbyte,sizeof(int)) < 0) goto open_io_error;
  pbyte += cbyte;

  if (read(f,&nfile,sizeof(int)) < 0) goto open_io_error;
  if (read(f,&maxp,sizeof(int64)) < 0) goto open_io_error;
  if (read(f,&freq,sizeof(int)) < 0) goto open_io_error;

  if (read(f,&nctg,sizeof(int)) < 0) goto open_io_error;

  P = Malloc(sizeof(Post_List),"Allocating post record");
  if (P == NULL)
    exit (1);
  P->name   = full;
  P->nlen   = strlen(full);
  P->maxp   = maxp;
  P->cache  = Malloc(POST_BLOCK*pbyte,"Allocating post list buffer\n");
  P->neps   = Malloc(nfile*sizeof(int64),"Allocating parts table of Post_List");
  P->perm   = Malloc(nctg*sizeof(int),"Allocating sort permutation");
  P->index  = Malloc(0x10000*sizeof(int64),"Allocating index array");
  if (P->cache == NULL || P->neps == NULL || P->perm == NULL || P->index == NULL)
    exit (1);

  if (read(f,P->perm,sizeof(int)*nctg) < 0) goto open_io_error;
  if (read(f,P->index,sizeof(int64)*0x10000) < 0) goto open_io_error;
  close(f);

  nels = 0;
  for (p = 1; p <= nfile; p++)
    { sprintf(P->name+P->nlen,"%d",p);
      copn = open(P->name,O_RDONLY);
      if (copn < 0)
        { fprintf(stderr,"%s: Table part %s is missing ?\n",Prog_Name,P->name);
          exit (1);
        }
      if (read(copn,&pb,sizeof(int)) < 0) goto part_io_error;
      if (read(copn,&cb,sizeof(int)) < 0) goto part_io_error;
      pb += cb;
      if (read(copn,&n,sizeof(int64)) < 0) goto part_io_error;
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
  P->freq  = freq;
  P->nctg  = nctg;

  sprintf(P->name+P->nlen,"%d",1);
  copn = open(P->name,O_RDONLY);
  if (copn < 0) goto part_io_error;

  if (lseek(copn,2*sizeof(int)+sizeof(int64),SEEK_SET) < 0)
    { fprintf(stderr,"\n%s: Could not seek file %s\n",Prog_Name,P->name);
      exit (1);
    }

  P->copn  = copn;
  P->part  = 1;

  More_Post_List(P);
  P->cidx = 0;

  return (P);

open_io_error:
  fprintf(stderr,"\n%s: IO error reading file %s/%s.post\n",Prog_Name,dir,root);
  exit (1);

part_io_error:
  fprintf(stderr,"\n%s: IO error reading post part file %s\n",Prog_Name,P->name);
  exit (1);
}

static void Free_Post_List(Post_List *P)
{ free(P->index);
  free(P->perm);
  free(P->neps);
  free(P->cache);
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
          if (P->copn < 0)
            { fprintf(stderr,"\n%s: Could not open post part file %s\n",Prog_Name,P->name);
              exit (1);
            }
          P->part = 1;
        }

      if (lseek(P->copn,sizeof(int)+sizeof(int64),SEEK_SET) < 0)
        { fprintf(stderr,"\n%s: Could not seek file %s\n",Prog_Name,P->name);
          exit (1);
        }

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
      if (P->cidx >= P->nels)
        { P->cptr = NULL;
          P->part = P->nthr+1;
          return;
        }
      sprintf(P->name+P->nlen,"%d",p);
      P->copn = open(P->name,O_RDONLY);
      if (P->copn < 0)
        { fprintf(stderr,"\n%s: Could not open post part file %s\n",Prog_Name,P->name);
          exit (1);
        }
      P->part = p;
    }

  if (lseek(P->copn,2*sizeof(int) + sizeof(int64) + i*P->pbyte,SEEK_SET) < 0)
    { fprintf(stderr,"\n%s: Could not seek file %s\n",Prog_Name,P->name);
      exit (1);
    }

  More_Post_List(P);
}


/***********************************************************************************************
 *
 *   LIST INDEX
 *
 **********************************************************************************************/

void JumpTo_Ktab_Index(Kmer_Stream *T, Post_List *P, int64 idx)
{ int    kpre;
  int64  bidx;
  int    count;
  uint8 *split = (uint8 *) (&count);

  (void) First_Post_Entry;

  GoTo_Kmer_Index(T,idx);
  kpre = (T->cpre >> 8);
  if (kpre > 0)
    bidx = P->index[kpre-1];
  else
    bidx = 0;
  kpre <<= 8;
  if (kpre > 0)
    GoTo_Kmer_Index(T,T->index[kpre-1]);
  else
    First_Kmer_Entry(T);
  while (T->cidx < idx)
    { count = Current_Count(T);
      bidx += split[0];
      Next_Kmer_Entry(T);
    }
  JumpTo_Post_Index(P,bidx);
}

static void Print_Index(Kmer_Stream *T, Post_List *P, int64 bidx, int64 eidx)
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

  JumpTo_Ktab_Index(T,P,bidx);
  for ( ; T->cidx < eidx; Next_Kmer_Entry(T))
    { count = Current_Count(T);
      printf(" %6lld: %s %2d\n",T->cidx,Current_Kmer(T,buffer),split[1]);
      fflush(stdout);
      for (p = split[0]; p > 0; p--)
        { memcpy(pptr,P->cptr,pbyte);
          memcpy(cptr,P->cptr+pbyte,cbyte);
          invert = ((cont & flag) != 0);
          if (invert)
            cont -= flag;
          printf("    %11lld.   %c %4d | %10lld\n",P->cidx,invert?'-':'+',perm[cont],post);
          fflush(stdout);
          Next_Post_Entry(P);
        }
    }
}

static int dna[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,

    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,

    0, 1, 0, 1, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,

    0, 1, 0, 1, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
  };

static int shiftup[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,

    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,

    0, 'C', 0, 'G', 0, 0, 0, 'T',
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,

    0, 'c', 0, 'g', 0, 0, 0, 't',
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
  };

static int64 Interpret(Kmer_Stream *T, char *x, int beg)
{ int   d, n;
  char *u;

  if (sscanf(x,"%d%n",&d,&n) == 1)
    { if (x[n] != 0)
        { fprintf(stderr,"%s: Indx %s is not an integer\n",Prog_Name,x);
          exit (1);
        }
      if (d >= T->nels)
        { fprintf(stderr,"%s: Index %s is out of bounds\n",Prog_Name,x);
          exit (1);
        }
      if (beg)
        return ((int64) d);
      else
        return ((int64) d+1);
    }
  for (n = 0; x[n] != '\0'; n++)
    if (!dna[(int) x[n]])
      { fprintf(stderr,"%s: String %s is not dna (acgt)\n",Prog_Name,x);
        exit (1);
      }
  if (n > T->kmer)
    { fprintf(stderr,"%s: String %s is longer than k-mer size (%d)\n",Prog_Name,x,T->kmer);
      exit (1);
    }
  u = Current_Kmer(T,NULL);
  strcpy(u,x);
  if (!beg)
    { n -= 1;
      while (n >= 0 && (u[n] == 't' || u[n] == 'T'))
        n -= 1;
      if (n < 0)
        return (T->nels);
      else
        u[n] = shiftup[(int) u[n]];
      n += 1;
    }
  while (n < T->kmer)
    u[n++] = 'a';
  GoTo_Kmer_String(T,u);
  free(u);
  return (T->cidx);
}

int main(int argc, char *argv[])
{ Kmer_Stream *T;
  Post_List   *P;
  int64        bidx, eidx;

  //  Process options

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("GIXshow");

    (void) flags;

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

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"          <address> = <int> | <dna:string>\n");
        exit (1);
      }
  }

  T = Open_Kmer_Stream(argv[1]);
  if (T == NULL)
    { fprintf(stderr,"%s: Cannot open k-mer table %s\n",Prog_Name,argv[1]);
      exit (1);
    }
  P = Open_Post_List(argv[1]);
  if (P == NULL)
    exit (1);

  if (argc == 2)
    { bidx = 0;
      eidx = T->nels;
    }
  else
    { char *x = argv[2];
      char *p = index(x,'-');
      if (p != NULL)
        { *p++ = '\0';
          bidx = Interpret(T,x,1);
          eidx = Interpret(T,p,0);
        }
      else
        { bidx = Interpret(T,x,1);
          eidx = Interpret(T,x,0);
        }
    }

  Print_Index(T,P,bidx,eidx); 

  Free_Post_List(P);
  Free_Kmer_Stream(T);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
