/*  Last edited: Feb 21 16:55 2025 (rd109) */
 /*******************************************************************************************
 *
 *  Adaptamer merge phase of a WGA.  Given indices of two genomes (at a specified frequency
 *   cutoff), the adaptemer matches between the k-mers are found in a novel, cache-coherent
 *   merge of the sorted k-mer tables for each genome and seed position pairs are output for
 *   each adaptemer match.
 *
 *  Author  :  Gene Myers
 *  Date    :  February 2023
 *  Last Mod: March 2023 - RD replace .las with .1aln
 *
 *******************************************************************************************/

#define VERSION "0.1"

#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include <sys/resource.h>

#include "libfastk.h"
#include "GDB.h"
#include "align.h"
#include "alncode.h"

#undef    DEBUG_SPLIT
#undef    DEBUG_MERGE
#undef    DEBUG_SORT
#undef    DEBUG_SEARCH
#undef    DEBUG_HIT
#undef    DEBUG_TUBE
#undef    DEBUG_ALIGN
#undef    DEBUG_ENTWINE

#define   CALL_ALIGNER
#define   BOX_ELIM

#undef    LOAD_SEQS

#define   MAX_INT64    0x7fffffffffffffffll

#define   TSPACE       100

static int PTR_SIZE = sizeof(void *);
static int OVL_SIZE = sizeof(Overlap);
static int EXO_SIZE = sizeof(Overlap) - sizeof(void *);

#define    BUCK_SHIFT     6
#define    BUCK_WIDTH    64  //  2^BUCK_SHIFT
#define    BUCK_ANTI    128  //  2*BUCK_WIDTH
#define    BOX_FUZZ      10

static char *Usage[] = { "[-vk] [-T<int(8)>] [-P<dir($TMPDIR)>] [<format(-paf)>]",
                         "[-f<int(10)>] [-c<int(85)> [-s<int(1000)>] [-l<int(100)>] [-i<float(.7)]",
                         "<source1:path>[<precursor>] [<source2:path>[<precursor>]]"
                       };

static int    FREQ;        //  -f: Adaptemer frequence cutoff parameter
static int    VERBOSE;     //  -v: Verbose output
static int    CHAIN_BREAK; //  -s
static int    CHAIN_MIN;   //  -c
static int    ALIGN_MIN;   //  -a
static double ALIGN_RATE;  //  1.-e
static int    NTHREADS;    //  -T
static char  *SORT_PATH;   //  -P
static int    KEEP;        //  -k
static int    SYMMETRIC;   //  -S
static int    SELF;        //  Comparing A to A, or A to B?
static int    OUT_TYPE;    //  -paf = 0; -psl = 1; -one = 2
static int    OUT_OPT;     //  -pafm = 1; -pafx = 2; all others = 0
static char  *ONE_PATH;    //  -one option path
static char  *ONE_ROOT;    //  -one option path

static char *PATH1, *PATH2;   //  GDB & GIX are PATHx/ROOTx[GEXTNx|.gix]
static char *ROOT1, *ROOT2;
static char *GEXTN1, *GEXTN2;

static char *SPATH1, *SPATH2; //  Path name of source if TYPEx <= IS_GDB
static int   TYPE1,  TYPE2;   //  Type of source file (see DNAsource.h)

static int   KMER;         //  K-mer length and # of threads from genome indices

static char *PAIR_NAME;    //  Prefixes for temporary files
static char *ALGN_UNIQ;
static char *ALGN_PAIR;

static int IBYTE;         // # of bytes for an entry in P1
static int IPOST;         // # of bytes in post of a P1 entry
static int ICONT;         // # of bytes in contig of a P1 entry
static int ISIGN;         // byte index in a P1 entry of the sign flag

static int JBYTE;         // # of bytes for an entry in P2
static int JPOST;         // # of bytes in post of a P2 entry
static int JCONT;         // # of bytes in contig of a P2 entry
static int JSIGN;         // byte index in a P2 entry of the sign flag

static int KBYTE;         // # of bytes for a k-mer table entry (both T1 & T2)
static int CBYTE;         // byte of k-mer table entry containing post count
static int LBYTE;         // byte of k-mer table entry containing lcp

static int64 AMXPOS;     //  longest contig in gdb1
static int64 BMXPOS;     //  longest contig in gdb1
static int64 MAXDAG;     //  AMXPOS + BMXPOS
static int   DBYTE;      // # of bytes for a pair diagonal or anti-diagonal

static int    NCONTS;    //  # of A contigs
static int    NPARTS;    //  # of panels A-contigs divided into
static int    ESHIFT;    //  shift to extract P1-contig # from a post

static int   *Select;    //  Select[bucket] = thread file for bucket
static int   *IDBsplit;  //  GDB split: contigs [DBsplit[p],DBsplit[p+1])
static int   *Perm1;     //  Sorted contig permutation of gdb1
static int   *Perm2;     //  Sorted contig permutation of gdb2

typedef struct
  { uint8 *bufr;
    uint8 *btop;
    uint8 *bend;
    int64 *buck;
    int    file;
    int    inum;
  } IOBuffer;

static IOBuffer *N_Units;  //  NTHREADS^2 IO units for + pair temporary files
static IOBuffer *C_Units;  //  NTHREADS^2 IO units for - pair temporary files

typedef struct
  { int   beg;
    int   end;
    int64 off;
  } Range;

extern int rmsd_sort(uint8 *array, int64 nelem, int rsize, int ksize,
                     int nparts, int64 *part, int nthreads, Range *range);

static void Clean_Exit(int status)
{ char *command;
  int   fail;

  if (status == 0 && KEEP)
    exit (0);

  fail = 0;
  if (TYPE2 <= IS_GDB)
    command = Malloc(strlen(PATH1)+strlen(ROOT1)+
                     strlen(PATH2)+strlen(ROOT2)+100,"Allocating command string");
  else
    command = Malloc(strlen(PATH1)+strlen(ROOT1)+100,"Allocating command string");

  if (command == NULL)
    fail = 1;
  else
    { if (TYPE1 <= IS_GDB)
        { if (TYPE1 < IS_GDB)
            sprintf(command,"GIXrm -fg %s/%s%s",PATH1,ROOT1,GEXTN1);
          else
            sprintf(command,"GIXrm -f %s/%s.gix",PATH1,ROOT1);
          if (system(command) != 0)
            fail = 1;
          free(ROOT1);
        }

      if (TYPE2 <= IS_GDB)
        { if (TYPE2 < IS_GDB)
            sprintf(command,"GIXrm -fg %s/%s%s",PATH2,ROOT2,GEXTN2);
          else
            sprintf(command,"GIXrm -f %s/%s.gix",PATH2,ROOT2);
          if (system(command) != 0)
            fail = 1;
          free(ROOT2);
        }

      free(command);
    }

  if (fail)
    fprintf(stderr,"\n%s: Warning: Could not successfully remove .1gdb/.gix\n",Prog_Name);

  exit (status);
}


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
    int     clone;      //  Is this a clone?
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
          Clean_Exit(1);
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
          Clean_Exit(1);
        }
      if (lseek(copn,2*sizeof(int)+sizeof(int64),SEEK_SET) < 0)
        { fprintf(stderr,"%s: Cannot advance post file %s to data part\n",Prog_Name,P->name);
          Clean_Exit(1);
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
    { Clean_Exit(1);
      exit (1);
    }
  sprintf(full,"%s/%s.gix",dir,root);
  f = open(full,O_RDONLY);
  if (f < 0)
    { fprintf(stderr,"%s: Cannot open post stub file %s/%s.post\n",Prog_Name,dir,root);
      Clean_Exit(1);
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
    Clean_Exit(1);
  P->name   = full;
  P->nlen   = strlen(full);
  P->maxp   = maxp;
  P->cache  = Malloc(POST_BLOCK*pbyte,"Allocating post list buffer\n");
  P->neps   = Malloc(nfile*sizeof(int64),"Allocating parts table of Post_List");
  P->perm   = Malloc(nctg*sizeof(int),"Allocating sort permutation");
  P->index  = Malloc(0x10000*sizeof(int64),"Allocating index array");
  if (P->cache == NULL || P->neps == NULL || P->perm == NULL || P->index == NULL)
    { Clean_Exit(1);
      exit (1);
    }

  if (read(f,P->perm,sizeof(int)*nctg) < 0) goto open_io_error;
  if (read(f,P->index,sizeof(int64)*0x10000) < 0) goto open_io_error;
  close(f);

  nels = 0;
  for (p = 1; p <= nfile; p++)
    { sprintf(P->name+P->nlen,"%d",p);
      copn = open(P->name,O_RDONLY);
      if (copn < 0)
        { fprintf(stderr,"%s: Table part %s is missing ?\n",Prog_Name,P->name);
          Clean_Exit(1);
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
          Clean_Exit(1);
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
  P->clone = 0;

  sprintf(P->name+P->nlen,"%d",1);
  copn = open(P->name,O_RDONLY);
  if (copn < 0) goto part_io_error;

  if (lseek(copn,2*sizeof(int)+sizeof(int64),SEEK_SET) < 0)
    { fprintf(stderr,"\n%s: Could not seek file %s\n",Prog_Name,P->name);
      Clean_Exit(1);
      exit (1);
    }

  P->copn  = copn;
  P->part  = 1;

  More_Post_List(P);
  P->cidx = 0;

  return (P);

open_io_error:
  fprintf(stderr,"\n%s: IO error reading file %s/%s.post\n",Prog_Name,dir,root);
  Clean_Exit(1);
  exit (1);

part_io_error:
  fprintf(stderr,"\n%s: IO error reading post part file %s\n",Prog_Name,P->name);
  Clean_Exit(1);
  exit (1);
}

Post_List *Clone_Post_List(Post_List *O)
{ Post_List *P;
  int copn;

  P = Malloc(sizeof(Post_List),"Allocating post record");
  if (P == NULL)
    Clean_Exit(1);

  *P = *O;
  P->clone = 1;

  P->cache = Malloc(POST_BLOCK*O->pbyte,"Allocating post list buffer\n");
  P->name  = Malloc(P->nlen+20,"Allocating post list buffer\n");
  if (P->cache == NULL || P->name == NULL)
    Clean_Exit(1);
  strncpy(P->name,O->name,P->nlen);

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
{ if (!P->clone)
    { free(P->index);
      free(P->perm);
      free(P->neps);
    }
  free(P->name);
  free(P->cache);
  if (P->copn >= 0)
    close(P->copn);
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
              Clean_Exit(1);
            }
          P->part = 1;
        }

      if (lseek(P->copn,2*sizeof(int)+sizeof(int64),SEEK_SET) < 0)
        { fprintf(stderr,"\n%s: Could not seek file %s\n",Prog_Name,P->name);
          Clean_Exit(1);
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

static inline void Current_Post(Post_List *P, uint8 *here)
{ memcpy(here,P->cptr,P->pbyte); }

static inline void GoTo_Post_Index(Post_List *P, int64 i)
{ int    p;

  if (P->cidx == i)
    return;
  P->cidx = i;

  p = 0;
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
          Clean_Exit(1);
        }

      P->part = p;
    }

  if (lseek(P->copn,2*sizeof(int) + sizeof(int64) + i*P->pbyte,SEEK_SET) < 0)
    { fprintf(stderr,"\n%s: Could not seek file %s\n",Prog_Name,P->name);
      Clean_Exit(1);
    }

  More_Post_List(P);
}

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
          Clean_Exit(1);
        }
      P->part = p;
    }

  if (lseek(P->copn,2*sizeof(int) + sizeof(int64) + i*P->pbyte,SEEK_SET) < 0)
    { fprintf(stderr,"\n%s: Could not seek file %s\n",Prog_Name,P->name);
      Clean_Exit(1);
    }

  More_Post_List(P);
}

/***********************************************************************************************
 *
 *   ADAPTAMER MERGE THREAD:  
 *     For each k-mer in T1
 *       Find the longest prefix match to one or more k-mers in T2.
 *       If the total # of positions of the k-mers in T2 <= FREQ, then output the pairs of positions
 *         from P1 & P2 to a file dependent on the slice of the P1 post and the sign of the match.
 *
 **********************************************************************************************/

static int cbyte[41] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
                         3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5,
                         6, 6, 6, 6, 7 };

static int mbyte[41] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0xc0, 0x30, 0x0c, 0x03, 0xc0, 0x30, 0x0c, 0x03, 0xc0, 0x30, 0x0c, 0x03,
                         0xc0, 0x30, 0x0c, 0x03, 0xc0, 0x30, 0x0c, 0x03, 0xc0, 0x30, 0x0c, 0x03,
                         0xc0, 0x30, 0x0c, 0x03, 0xc0 };

#define POST_BUF_LEN  0x1000
#define POST_BUF_MASK 0x0fff

typedef struct
  { Kmer_Stream *T1;
    Kmer_Stream *T2;
    Post_List   *P1;
    Post_List   *P2;
    int          flip;
    int          tid;
    int          pbeg, pend;
    uint8       *cache;
    IOBuffer    *nunit;
    IOBuffer    *cunit;
    int64        nhits;
    int64        g1len;
    int64        tseed;
  } SP;

static void *merge_thread(void *args)
{ SP *parm = (SP *) args;
  int tid          = parm->tid;
  int flip         = parm->flip;
  uint8 *cache     = parm->cache;
  IOBuffer  *nunit = parm->nunit;
  IOBuffer  *cunit = parm->cunit;
  Kmer_Stream *T1  = parm->T1;
  Post_List   *P1  = parm->P1;
  Kmer_Stream *T2  = parm->T2;
  Post_List   *P2  = parm->P2;

  int64   tbeg, tend;

  int     cpre;
  uint8  *ctop, *suf1;

  int     eorun, plen;
  uint8  *rcur, *rend;
  uint8  *vlcp[KMER+1];

  int64   apost;
  uint8  *aptr = (uint8 *) (&apost);

  uint8  *vlow, *vhgh;
  int     pdx, cdx;
  int64   post[POST_BUF_LEN + FREQ];

  int     qcnt, pcnt;
  int64   nhits, g1len, tseed;

#ifdef DEBUG_MERGE
  int64   Tdp;
  char   *tbuffer;
#endif

  { int j;

    for (j = 0; j < NPARTS; j++)
      { nunit[j].bend = nunit[j].bufr + (1000000-(IBYTE+JBYTE+1));
        nunit[j].btop = nunit[j].bufr;
        cunit[j].bend = cunit[j].bufr + (1000000-(IBYTE+JBYTE+1));
        cunit[j].btop = cunit[j].bufr;
      }
  }

  cpre  = -1;
  ctop  = cache;
  vhgh  = cache;
#ifdef DEBUG_MERGE
  tbuffer = Current_Kmer(T1,NULL);
#endif

  nhits = 0;
  g1len = 0;
  tseed = 0;
  apost = 0;
  bzero(post,(POST_BUF_LEN+FREQ)*sizeof(int64));

  First_Post_Entry(P1);
  First_Post_Entry(P2);
  First_Kmer_Entry(T1);
  First_Kmer_Entry(T2);

  if (tid != 0)
    { GoTo_Post_Index(P1,P1->index[parm->pbeg]);
      GoTo_Post_Index(P2,P2->index[parm->pbeg]);
      GoTo_Kmer_Index(T1,T1->index[(parm->pbeg << 8) | 0xff]);
      GoTo_Kmer_Index(T2,T2->index[(parm->pbeg << 8) | 0xff]);
    }
  tend = T1->index[(parm->pend << 8) | 0xff];

  plen = 12;                         //  Keep the dumb compiler happy
  vlcp[plen] = rcur = rend = cache;
  vlow  = cache-KBYTE;
  vhgh  = cache;
  pdx   = POST_BUF_MASK;
  cdx   = 0;
  eorun = 0;

  qcnt = -1;
  tbeg = T1->cidx;
  while (T1->cidx < tend)
    { suf1 = T1->csuf;
#ifdef DEBUG_MERGE
      printf("Doing %s (%lld)\n",Current_Kmer(T1,tbuffer),T1->cidx); fflush(stdout);
#endif

      if (T1->cpre != cpre)  //  New prefix panel
        { int64  bidx;
          uint8 *cp;

          if (VERBOSE && tid == 0)
            { if (tbeg == tend)
                pcnt = 100;
              else
                pcnt = ((T1->cidx - tbeg) * 100) / (tend-tbeg); 
              if (pcnt > qcnt)
                { fprintf(stderr,"\r    Completed %3d%%",pcnt);
                  fflush(stderr);
                }
              qcnt = pcnt;
            }

          //  skip remainder of cache and T2 entries < T1->cpre

          bidx = 0;
          for (cp = vhgh; cp < ctop; cp += KBYTE)
            bidx += cp[CBYTE];
          for (cpre = T1->cpre; T2->cpre < cpre; Next_Kmer_Entry(T2))
            bidx += T2->csuf[CBYTE];
          JumpTo_Post_Index(P2,bidx);

#ifdef DEBUG_MERGE
          Tdp = T2->cidx;
          printf("Loading %lld %06x ...",Tdp,cpre); fflush(stdout);
#endif

          //  load cahce with T2 entries = T1->cpre

          for (cp = cache; T2->cpre == cpre; cp += KBYTE)
            { memcpy(cp,T2->csuf,KBYTE);
              Next_Kmer_Entry(T2);
            }
          ctop = cp;
          ctop[LBYTE] = 11;

          //  if cache is empty then skip to next T1 entry with a greater prefix than cpre

          if (ctop == cache)
            { bidx = 0;
              while (T1->cpre == cpre)
                { bidx += T1->csuf[CBYTE];
                  Next_Kmer_Entry(T1);
                }
              JumpTo_Post_Index(P1,bidx);
#ifdef DEBUG_MERGE
              printf(" ... Empty => to %06x in T1\n",T1->cpre);
#endif
              continue;
            }

          //  start adpatermer merge for prefix cpre

          plen = 12;
          vlcp[plen] = rcur = rend = cache;
          vlow  = cache-KBYTE;
          vhgh  = cache;
          pdx   = POST_BUF_MASK;
          cdx   = 0;
          eorun = 0;

#ifdef DEBUG_MERGE
          printf("... to %lld rcur = rend = %lld, eorun = 0, plen = 12\n",
                 T2->cidx,Tdp+(rcur-cache)/KBYTE);
          fflush(stdout);
#endif
        }

#define ADVANCE(l)					\
{ if (l >= vhgh)					\
    { int n;						\
							\
      for (n = l[CBYTE]; n > 0; n--)			\
        { pdx = (pdx+1) & POST_BUF_MASK;		\
          Current_Post(P2,(uint8 *) (post+pdx));	\
          Next_Post_Entry(P2);				\
        }						\
      vhgh = l+KBYTE;					\
    }							\
  cdx = (cdx + l[CBYTE]) & POST_BUF_MASK;		\
  l += KBYTE;						\
}

      // eorun = 0: rcur <= rend, n[1..plen] = rcur..rend, n[plen+1] < rend[plen+1]
      // eorun = 1: rcur <  rend, n[1..plen] = rcur..rend-1, lcp(rend) < plen 

      else
         { int nlcp;

           nlcp = suf1[LBYTE];
           if (nlcp > plen)
             goto pairs;
           else if (nlcp == plen)
             { if (eorun)
                 goto pairs;
             }
           else
             { if ( ! eorun)
                 ADVANCE(rend)
               while (rend[LBYTE] > nlcp)
                 ADVANCE(rend)
               plen = rend[LBYTE];
               if (plen < nlcp)
                 { eorun = 1;
                   plen  = nlcp;
                   goto pairs;
                 }
               eorun = 0;
               rcur  = rend;
             }
         }

       while (plen < KMER)
         { int h, m, c, d;

           h = cbyte[plen];
           m = mbyte[plen];
           c = suf1[h] & m;
           for (d = rend[h]&m; d < c; d = rend[h]&m)
             { ADVANCE(rend)
               if (rend[LBYTE] < plen)
                 { eorun = 1;
                   goto pairs;
                 }
             }
           if (d > c)
             goto pairs;
           plen += 1;
           vlcp[plen] = rcur = rend;
         }
       ADVANCE(rend)
       eorun = 1;

       //  Get pairs;

    pairs:

#ifdef DEBUG_MERGE
      printf("-> %d[%lld,%lld] %d",plen,Tdp+(vlcp[plen]-cache)/KBYTE,Tdp+(rend-cache)/KBYTE,eorun);
      printf("  [%lld,%lld] %d",Tdp+(vlow-cache)/KBYTE,Tdp+(vhgh-cache)/KBYTE,pdx);
      fflush(stdout);
#endif

      { int       freq, lcs, udx;
        int       asign, acont, adest, bsign;
        IOBuffer *ou;
        uint8    *l, *vcp, *jptr, *btop;
        int       m, n, k, b;
        
        freq = 0;
        vcp = vlcp[plen];
        if (vcp <= vlow)
          {
#ifdef DEBUG_MERGE
            printf("   vlow <= vcp\n");
            fflush(stdout);
#endif
            goto empty;
          }
      
        for (l = rend-KBYTE; l >= vcp; l -= KBYTE)
          { freq += l[CBYTE];
            if (freq >= FREQ)
              { vlow = l;
#ifdef DEBUG_MERGE
                printf("   %d vlow = %lld\n",freq,Tdp+(l-cache)/KBYTE);
                fflush(stdout);
#endif
                goto empty;
              }
          }
        lcs = freq;
        l   = rend;
        if ( ! eorun)
          { udx = cdx;
            l = rend;
            freq += l[CBYTE];
            if (freq >= FREQ)
              {
#ifdef DEBUG_MERGE
                printf("   %d too high at %lld\n",freq,Tdp+(l-cache)/KBYTE);
                fflush(stdout);
#endif
                goto empty;
              }
            ADVANCE(l)
            while (l[LBYTE] >= plen)
              { freq += l[CBYTE];
                if (freq >= FREQ)
                  {
#ifdef DEBUG_MERGE
                    printf("   %d too high at %lld\n",freq,Tdp+(l-cache)/KBYTE);
                    fflush(stdout);
#endif
                    cdx = udx;
                    goto empty;
                  }
                ADVANCE(l)
              }
            cdx = udx;
          }
#ifdef DEBUG_MERGE
        printf("    [%lld,%lld) w %d posts\n",Tdp+(vcp-cache)/KBYTE,Tdp+(l-cache)/KBYTE,freq);
        fflush(stdout);
#endif

        if (cdx >= lcs)
          b = cdx-lcs;
        else
          b = (cdx+POST_BUF_LEN) - lcs;
        if (b + freq > POST_BUF_LEN)
          { m = (b+freq) & POST_BUF_MASK;
            for (m--; m >= 0; m--)
              post[POST_BUF_LEN+m] = post[m];
          }

        if (flip)

          for (n = suf1[CBYTE]; n > 0; n--)
            { Current_Post(P1,aptr);
              bsign = (aptr[JSIGN] & 0x80);
              if (bsign == 0)
                g1len += 1;
              jptr  = (uint8 *) (post+b);
              for (k = 0; k < freq; k++)
                {
#ifdef DEBUG_MERGE
                  if (n == suf1[CBYTE])
                    { int64 ip;
                      int   ss;
  
                      ss = (jptr[ISIGN] & 0x80);
                      jptr[ISIGN] &= 0x7f;
                      printf("      %ld: %c",((int64 *) jptr)-post,ss?'-':'+');
                      ip = 0;
                      memcpy((uint8 *) (&ip),jptr+IPOST,ICONT);
                      // printf(" %4d",Perm1[ip]);
                      printf(" %4lld",ip);
                      ip = 0;
                      memcpy((uint8 *) (&ip),jptr,IPOST);
                      printf(" %9lld\n",ip);
                      fflush(stdout);
                      if (ss)
                        jptr[ISIGN] |= 0x80;
                    }
#endif
                  if (jptr[ISIGN] & 0x80)
                    { jptr += sizeof(int64);
                      continue;
                    }
                  acont = (*((uint64 *)jptr) >> ESHIFT);
                  adest = Select[acont];
                  if (bsign)
                    ou = cunit + adest;
                  else
                    ou = nunit + adest;
                  btop = ou->btop;
                  *btop++ = plen;
                  memcpy(btop,jptr,IBYTE);
                  btop += IBYTE;
                  memcpy(btop,aptr,JBYTE);
                  btop += JBYTE;
  
                  nhits += 1;
                  tseed += plen;
  
                  ou->buck[acont] += 1;
  
                  if (btop >= ou->bend)
                    { if (write(ou->file,ou->bufr,btop-ou->bufr) < 0)
                        { fprintf(stderr,"%s: IO write to file %s/%s.%d.%c failed\n",
                                         Prog_Name,SORT_PATH,PAIR_NAME,ou->inum,
                                         (bsign == (jptr[JSIGN] & 0x80)) ? 'N' : 'C');
                          Clean_Exit(1);
                        }
                      ou->btop = ou->bufr;
                    }
                  else
                    ou->btop = btop;
  
                  jptr += sizeof(int64);
                }
#ifdef DEBUG_MERGE
              { int64 ip;
  
                if (n == suf1[CBYTE])
                  printf("   vs\n");
                aptr[JSIGN] &= 0x7f;
                printf("      %lld: %c",P1->cidx,bsign?'-':'+');
                ip = 0;
                memcpy((uint8 *) (&ip),aptr+JPOST,JCONT);
                // printf(" %4d",Perm2[ip]);
                printf(" %4lld",ip);
                ip = 0;
                memcpy((uint8 *) (&ip),aptr,JPOST);
                printf(" %9lld\n",ip);
                fflush(stdout);
                if (bsign)
                  aptr[JSIGN] |= 0x80;
              }
#endif
              Next_Post_Entry(P1);
            }

        else
  
          for (n = suf1[CBYTE]; n > 0; n--)
            { Current_Post(P1,aptr);
              asign = (aptr[ISIGN] & 0x80);
              if (asign)
                { Next_Post_Entry(P1);
                  continue;
                }
              nhits += freq;
              g1len += 1;
              tseed += freq * plen;
              acont = (apost >> ESHIFT);
              adest = Select[acont];
              jptr  = (uint8 *) (post+b);
              for (k = 0; k < freq; k++)
                { if (jptr[JSIGN] & 0x80)
                    ou = cunit + adest;
                  else
                    ou = nunit + adest;
                  btop = ou->btop;
                  *btop++ = plen;
                  memcpy(btop,aptr,IBYTE);
                  btop += IBYTE;
                  memcpy(btop,jptr,JBYTE);
                  btop += JBYTE;
  
                  ou->buck[acont] += 1;
  
#ifdef DEBUG_MERGE
                  if (n == suf1[CBYTE])
                    { int64 ip;
                      int   ss;
  
                      ss = (jptr[JSIGN] & 0x80);
                      jptr[JSIGN] &= 0x7f;
                      printf("      %ld: %c",((int64 *) jptr)-post,ss?'-':'+');
                      ip = 0;
                      memcpy((uint8 *) (&ip),jptr+JPOST,JCONT);
                      printf(" %4d",Perm2[ip]);
                      ip = 0;
                      memcpy((uint8 *) (&ip),jptr,JPOST);
                      printf(" %9lld\n",ip);
                      fflush(stdout);
                      if (ss)
                        jptr[JSIGN] |= 0x80;
                    }
#endif
  
                  if (btop >= ou->bend)
                    { if (write(ou->file,ou->bufr,btop-ou->bufr) < 0)
                        { fprintf(stderr,"%s: IO write to file %s/%s.%d.%c failed\n",
                                         Prog_Name,SORT_PATH,PAIR_NAME,ou->inum,
                                         (asign == (jptr[JSIGN] & 0x80)) ? 'N' : 'C');
                          Clean_Exit(1);
                        }
                      ou->btop = ou->bufr;
                    }
                  else
                    ou->btop = btop;
  
                  jptr += sizeof(int64);
                }
#ifdef DEBUG_MERGE
              { int64 ip;
  
                if (n == suf1[CBYTE])
                  printf("   vs\n");
                printf("      %lld: +",P1->cidx);
                ip = 0;
                memcpy((uint8 *) (&ip),aptr+IPOST,ICONT);
                printf(" %4d",Perm1[ip]);
                ip = 0;
                memcpy((uint8 *) (&ip),aptr,IPOST);
                printf(" %9lld\n",ip);
                fflush(stdout);
              }
#endif
              Next_Post_Entry(P1);
            }

        Next_Kmer_Entry(T1);
        continue;
      }

    empty:
      JumpTo_Post_Index(P1,(int64) suf1[CBYTE]);
      Next_Kmer_Entry(T1);
    }

  { int j;

    for (j = 0; j < NPARTS; j++)
      { if (nunit[j].btop > nunit[j].bufr)
          if (write(nunit[j].file,nunit[j].bufr,nunit[j].btop-nunit[j].bufr) < 0)
            { fprintf(stderr,"%s: IO write to file %s/%s.%d.N failed\n",
                             Prog_Name,SORT_PATH,PAIR_NAME,nunit[j].inum);
              Clean_Exit(1);
            }
        if (cunit[j].btop > cunit[j].bufr)
          if (write(cunit[j].file,cunit[j].bufr,cunit[j].btop-cunit[j].bufr) < 0)
            { fprintf(stderr,"%s: IO write to file %s/%s.%d.C failed\n",
                             Prog_Name,SORT_PATH,PAIR_NAME,cunit[j].inum);
              Clean_Exit(1);
            }
      }
  }

  parm->nhits = nhits;
  parm->g1len = g1len;
  parm->tseed = tseed;
  return (NULL);
}

#ifdef DEBUG_MERGE

static char dna[4] = { 'a', 'c', 'g', 't' };

static char *fmer[256], _fmer[1280];

static void setup_fmer_table()
{ char *t;
  int   i, l3, l2, l1, l0;
  static int done = 0;

  if (done) return;
  done = 1;

  i = 0;
  t = _fmer;
  for (l3 = 0; l3 < 4; l3++)
   for (l2 = 0; l2 < 4; l2++)
    for (l1 = 0; l1 < 4; l1++)
     for (l0 = 0; l0 < 4; l0++)
       { fmer[i] = t;
         *t++ = dna[l3];
         *t++ = dna[l2];
         *t++ = dna[l1];
         *t++ = dna[l0];
         *t++ = 0;
         i += 1;
       }
}

char *Current_Cachmer(Kmer_Stream *S, uint8 *csuf, char *seq)
{ int    cpre  = S->cpre-1;
  int    hbyte = S->hbyte;

  setup_fmer_table();

  { int    j;
    uint8 *a;
    char  *s;

    s = seq;
    switch (S->ibyte)
    { case 3:
        memcpy(s,fmer[cpre>>16],4);
        s += 4;
        memcpy(s,fmer[cpre>>8 & 0xff],4);
        s += 4;
        memcpy(s,fmer[cpre&0xff],4);
        s += 4;
        break;
      case 2:
        memcpy(s,fmer[cpre>>8],4);
        s += 4;
        memcpy(s,fmer[cpre&0xff],4);
        s += 4;
        break;
      case 1:
        memcpy(s,fmer[cpre],4);
        s += 4;
        break;
    }

    a = csuf;
    for (j = 0; j < hbyte; j++, s += 4)
      memcpy(s,fmer[a[j]],4);
    seq[S->kmer] = '\0';
  }

  return (seq);
}

#endif


static void *self_merge_thread(void *args)
{ SP *parm = (SP *) args;
  int tid          = parm->tid;
  uint8 *cache     = parm->cache;
  IOBuffer  *nunit = parm->nunit;
  IOBuffer  *cunit = parm->cunit;
  Kmer_Stream *T1  = parm->T1;
  Post_List   *P2  = parm->P1;

  int64   tbeg, tend;

  int     cpre;
  uint8  *ctop, *suf1;

  int     eorun, plen;
  uint8  *rcur, *rend;
  uint8  *vlcp[KMER+1];

  uint8  *vlow, *vhgh;
  int     pdx, cdx, odx;
  int64   post[POST_BUF_LEN + FREQ];

  int     qcnt, pcnt;
  int64   nhits, g1len, tseed;

#ifdef DEBUG_MERGE
  int64   Tdp;
  char   *tbuffer;
#endif

  { int j;

    for (j = 0; j < NPARTS; j++)
      { nunit[j].bend = nunit[j].bufr + (1000000-(IBYTE+JBYTE+1));
        nunit[j].btop = nunit[j].bufr;
        cunit[j].bend = cunit[j].bufr + (1000000-(IBYTE+JBYTE+1));
        cunit[j].btop = cunit[j].bufr;
      }
  }

  memset(post,0,sizeof(int64)*(POST_BUF_LEN+FREQ));

  ctop  = cache;
  vhgh  = cache;
#ifdef DEBUG_MERGE
  tbuffer = Current_Kmer(T1,NULL);
#endif

  nhits = 0;
  g1len = 0;
  tseed = 0;

  First_Post_Entry(P2);
  First_Kmer_Entry(T1);

  if (tid != 0)
    { GoTo_Post_Index(P2,P2->index[parm->pbeg]);
      GoTo_Kmer_Index(T1,T1->index[(parm->pbeg<<8) | 0xff]);
    }

  vlcp[11] = rcur = rend = suf1 = cache;   //  Keep the dumb compiler happy
  vlow  = cache-KBYTE;
  vhgh  = cache;
  pdx   = POST_BUF_MASK;
  cdx   = 0;
  odx   = 0;
  plen  = 0;
  eorun = 0;

  qcnt = -1;
  tend = T1->index[(parm->pend<<8) | 0xff];
  tbeg = T1->cidx;
  for (suf1 = ctop; 1; suf1 += KBYTE)
    { if (suf1 >= ctop)
        { uint8 *cp;
          int    i;

          if (VERBOSE && tid == 0)
            { if (tbeg == tend)
                pcnt = 100;
              else
                pcnt = ((T1->cidx - tbeg) * 100) / (tend-tbeg); 
              if (pcnt > qcnt)
                { fprintf(stderr,"\r    Completed %3d%%",pcnt);
                  fflush(stderr);
                }
              qcnt = pcnt;
            }

          if (T1->cidx >= tend)
            break;

#ifdef DEBUG_MERGE
          Tdp = T1->cidx;
          printf("Loading %lld %06x ...",Tdp,cpre); fflush(stdout);
#endif

          //  load cahce with all T1 entries with T1->cpre

          cpre = T1->cpre;
          for (cp = cache; T1->cpre == cpre; cp += KBYTE)
            { memcpy(cp,T1->csuf,KBYTE);
              Next_Kmer_Entry(T1);
            }
          ctop = cp;
          ctop[LBYTE] = 11;

          //  start adpatermer merge for prefix cpre

          vlcp[11] = rcur = rend = suf1 = cache;
          vlow  = cache-KBYTE;
          vhgh  = cache;
          pdx   = POST_BUF_MASK;
          cdx   = 0;
          odx   = 0;

          ADVANCE(rend);
          plen = rend[LBYTE];
          for (i = rcur[LBYTE]; i <= plen; i++)
            vlcp[i] = rcur;
          eorun = (plen <= 11);

#ifdef DEBUG_MERGE
          printf("... to %lld\n",T1->cidx);
          fflush(stdout);

          printf("Doing %s (%lld)\n",Current_Cachmer(T1,suf1,tbuffer),(suf1-cache)/KBYTE+Tdp);
          fflush(stdout);
#endif
        }

      // suf1 = rend-1
      // eorun = 0: rcur <= rend, suf1[1..plen] = rcur..rend, suf1[plen+1] < rend[plen+1]
      // eorun = 1: rcur <  rend, suf1[1..plen] = rcur..rend-1, lcp(rend) < plen 

      else
        { int i;

#ifdef DEBUG_MERGE
          printf("Doing %s (%lld)\n",Current_Cachmer(T1,suf1,tbuffer),(suf1-cache)/KBYTE+Tdp);
          fflush(stdout);
#endif

          odx = cdx;
          if (eorun)
            plen = suf1[LBYTE];
          ADVANCE(rend);
          if (rend[LBYTE] < plen)
            eorun = 1;
          else if (rend[LBYTE] == plen)
            eorun = 0;
          else
            { rcur = rend-KBYTE;
              for (i = plen+1; i <= rend[LBYTE]; i++)
                vlcp[i] = rcur;
              eorun = 0;
              plen = rend[LBYTE];
            }
        }

#ifdef DEBUG_MERGE
      printf("-> %d[%lld,%lld] %d",plen,Tdp+(vlcp[plen]-cache)/KBYTE,Tdp+(rend-cache)/KBYTE,eorun);
      printf("  [%lld,%lld] %d\n",Tdp+(vlow-cache)/KBYTE,Tdp+(vhgh-cache)/KBYTE,pdx);
      fflush(stdout);
#endif

       //  Get pairs;

      { int       freq, lcs, mlen, udx;
        int       isign, jsign;
        int64     ipost;
        int       icont, idest;
        uint8    *iptr, *jptr;
        IOBuffer *ou;
        uint8    *l, *vcp, *btop;
        int       m, n, k, b;

        if (suf1[CBYTE] > 1)
          { if (suf1[CBYTE] >= FREQ)
              goto empty;
            lcs = freq = suf1[CBYTE];
            mlen = KMER;
          }

        else
        
          { freq = 0;
            vcp = vlcp[plen];
            if (vcp <= vlow)
              {
#ifdef DEBUG_MERGE
                printf("   vlow <= vcp\n");
                fflush(stdout);
#endif
                goto empty;
              }
          
            for (l = rend-KBYTE; l >= vcp; l -= KBYTE)
              { freq += l[CBYTE];
                if (freq >= FREQ)
                  { vlow = l;
#ifdef DEBUG_MERGE
                    printf("   %d vlow = %lld\n",freq,Tdp+(l-cache)/KBYTE);
                    fflush(stdout);
#endif
                    goto empty;
                  }
              }
            lcs = freq;
            l   = rend;
            if ( ! eorun)
              { udx = cdx;
                l = rend;
                freq += l[CBYTE];
                if (freq >= FREQ)
                  {
#ifdef DEBUG_MERGE
                    printf("   %d too high at %lld\n",freq,Tdp+(l-cache)/KBYTE);
                    fflush(stdout);
#endif
                    goto empty;
                  }
                ADVANCE(l)
                while (l[LBYTE] >= plen)
                  { freq += l[CBYTE];
                    if (freq >= FREQ)
                      {
#ifdef DEBUG_MERGE
                        printf("   %d too high at %lld\n",freq,Tdp+(l-cache)/KBYTE);
                        fflush(stdout);
#endif
                        cdx = udx;
                        goto empty;
                      }
                    ADVANCE(l)
                  }
                cdx = udx;
              }
#ifdef DEBUG_MERGE
            printf("    [%lld,%lld) w %d posts\n",Tdp+(vcp-cache)/KBYTE,Tdp+(l-cache)/KBYTE,freq);
            fflush(stdout);
#endif
          mlen = plen;
        }

        if (cdx >= lcs)
          b = cdx-lcs;
        else
          b = (cdx+POST_BUF_LEN) - lcs;
        if (b + freq > POST_BUF_LEN)
          { m = (b+freq) & POST_BUF_MASK;
            for (m--; m >= 0; m--)
              post[POST_BUF_LEN+m] = post[m];
          }

        nhits += suf1[CBYTE] * (freq-1);
        g1len += suf1[CBYTE];
        tseed += suf1[CBYTE] * (freq-1) * plen;

        iptr = (uint8 *) (post+odx);
        for (n = suf1[CBYTE]; n > 0; n--, iptr += sizeof(int64))
          { isign = (iptr[ISIGN] & 0x80);
            if (isign)
              iptr[ISIGN] &= 0x7f;
            ipost = *((int64 *) iptr);
            icont = (ipost >> ESHIFT);
            idest = Select[icont];

            jptr  = (uint8 *) (post+b);
            for (k = 0; k < freq; k++, jptr += sizeof(int64))
              {
#ifdef DEBUG_MERGE
                if (n == suf1[CBYTE])
                  { int64 ip;
                    int   ss;

                    if (iptr == jptr)
                      ss = isign;
                    else
                      { ss = (jptr[JSIGN] & 0x80);
                        jptr[JSIGN] &= 0x7f;
                      }
                    printf("      %ld: %c",((int64 *) jptr)-post,ss?'-':'+');
                    ip = 0;
                    memcpy((uint8 *) (&ip),jptr+JPOST,JCONT);
                    printf(" %4d",Perm2[ip]);
                    ip = 0;
                    memcpy((uint8 *) (&ip),jptr,JPOST);
                    printf(" %9lld\n",ip);
                    fflush(stdout);
                    if (ss && iptr != jptr)
                      jptr[JSIGN] |= 0x80;
                  }
#endif
                if (jptr == iptr)
                  continue;

                jsign = (jptr[JSIGN] & 0x80);
                if (isign == jsign)
                  ou = nunit + idest;
                else
                  ou = cunit + idest;
                btop = ou->btop;
                *btop++ = mlen;
                memcpy(btop,iptr,IBYTE);
                btop += IBYTE;
                memcpy(btop,jptr,JBYTE);
                btop += JBYTE;

                ou->buck[icont] += 1;

                if (btop >= ou->bend)
                  { if (write(ou->file,ou->bufr,btop-ou->bufr) < 0)
                      { fprintf(stderr,"%s: IO write to file %s/%s.%d.%c failed\n",
                                       Prog_Name,SORT_PATH,PAIR_NAME,
                                       ou->inum,(isign == jsign) ? 'N' : 'C');
                        Clean_Exit(1);
                      }
                    ou->btop = ou->bufr;
                  }
                else
                  ou->btop = btop;
              }

#ifdef DEBUG_MERGE
            { int64 ip;

              if (n == suf1[CBYTE])
                printf("   vs\n");
              printf("      %ld: %c",((int64 *) iptr)-post,isign?'-':'+');
              ip = 0;
              memcpy((uint8 *) (&ip),iptr+IPOST,ICONT);
              printf(" %4d",Perm1[ip]);
              ip = 0;
              memcpy((uint8 *) (&ip),iptr,IPOST);
              printf(" %9lld\n",ip);
              fflush(stdout);
            }
#endif

            if (isign)
              iptr[ISIGN] |= 0x80;
          }
 empty: continue;
      }
    }

  { int j;

    for (j = 0; j < NPARTS; j++)
      { if (nunit[j].btop > nunit[j].bufr)
          if (write(nunit[j].file,nunit[j].bufr,nunit[j].btop-nunit[j].bufr) < 0)
            { fprintf(stderr,"%s: IO write to file %s/%s.%d.N failed\n",
                             Prog_Name,SORT_PATH,PAIR_NAME,nunit[j].inum);
              Clean_Exit(1);
            }
        if (cunit[j].btop > cunit[j].bufr)
          if (write(cunit[j].file,cunit[j].bufr,cunit[j].btop-cunit[j].bufr) < 0)
            { fprintf(stderr,"%s: IO write to file %s/%s.%d.C failed\n",
                             Prog_Name,SORT_PATH,PAIR_NAME,cunit[j].inum);
              Clean_Exit(1);
            }
      }
  }

  parm->nhits = nhits/2;
  parm->g1len = g1len;
  parm->tseed = tseed/2;
  return (NULL);
}


static void adaptamer_merge(Kmer_Stream *T1, Kmer_Stream *T2,
                            Post_List *P1,   Post_List *P2)
{ SP         parm[NTHREADS];
#ifndef DEBUG_MERGE
  pthread_t  threads[NTHREADS];
#endif
  uint8     *cache;
  int64      nhits, g1len, tseed;
  int        i;

  { Kmer_Stream *tp;
    uint8       *ent;
#ifdef DEBUG_SPLIT
    char        *seq;
#endif
    int   t;
    int64 p;

    if (T1->nels > T2->nels)
      tp = T1;
    else
      tp = T2;
    ent = Current_Entry(tp,NULL);
#ifdef DEBUG_SPLIT
    seq = Current_Kmer(tp,NULL);
#endif
    parm[0].pbeg = 0;
    for (t = 1; t < NTHREADS; t++)
      { p = (tp->nels * t) / NTHREADS;
        GoTo_Kmer_Index(tp,p);
        if (p >= tp->nels)
          parm[t].pbeg = 0xffff;
        else
          { ent = Current_Entry(tp,ent);
            parm[t].pbeg = (tp->cpre >> 8);
          }
      }
    for (t = 0; t < NTHREADS-1; t++)
      parm[t].pend = parm[t+1].pbeg;
    parm[NTHREADS-1].pend = 0xffff;
  }

  parm[0].T1 = T1;
  parm[0].T2 = T2;
  parm[0].P1 = P1;
  parm[0].P2 = P2;
  for (i = 1; i < NTHREADS; i++)
    { parm[i].T1 = Clone_Kmer_Stream(T1);
      parm[i].T2 = Clone_Kmer_Stream(T2);
      parm[i].P1 = Clone_Post_List(P1);
      parm[i].P2 = Clone_Post_List(P2);
    }

  if (SYMMETRIC && P1->maxp > P2->maxp)
    cache = Malloc(NTHREADS*(P1->maxp+1)*KBYTE,"Allocating cache");
  else
    cache = Malloc(NTHREADS*(P2->maxp+1)*KBYTE,"Allocating cache");
  if (cache == NULL)
    Clean_Exit(1);

  for (i = 0; i < NTHREADS; i++)
    { IOBuffer *nu, *cu;

      parm[i].tid   = i;
      if (SYMMETRIC && P1->maxp > P2->maxp)
        parm[i].cache = cache + i * (P1->maxp+1) * KBYTE;
      else
        parm[i].cache = cache + i * (P2->maxp+1) * KBYTE;
      parm[i].nunit = nu = N_Units + i * NPARTS;
      parm[i].cunit = cu = C_Units + i * NPARTS;
      bzero(nu[0].buck,sizeof(int64)*NCONTS);
      bzero(cu[0].buck,sizeof(int64)*NCONTS);
      parm[i].flip = 0;
    }

  nhits = g1len = tseed = 0;

  if (VERBOSE)
    { fprintf(stderr,"\n  Starting adaptive seed merge for G1\n");
      fflush(stderr);
    }

#ifdef DEBUG_MERGE
  for (i = 0; i < NTHREADS; i++)
    merge_thread(parm+i);
#else
  for (i = 1; i < NTHREADS; i++)
    pthread_create(threads+i,NULL,merge_thread,parm+i);
  merge_thread(parm);

  for (i = 1; i < NTHREADS; i++)
    pthread_join(threads[i],NULL);
#endif

  if (VERBOSE)
    { fprintf(stderr,"\r    Completed 100%%\n");
      fflush(stderr);
    }

  for (i = 0; i < NTHREADS; i++)
    { nhits += parm[i].nhits;
      g1len += parm[i].g1len;
      tseed += parm[i].tseed;
    }

  if (SYMMETRIC)
    { if (VERBOSE)
        { fprintf(stderr,"\n  Starting adaptive seed merge for G2\n");
          fflush(stderr);
        }

      for (i = 0; i < NTHREADS; i++)
        { Kmer_Stream *u;
          Post_List   *p;

          u = parm[i].T1;
          parm[i].T1 = parm[i].T2;
          parm[i].T2 = u;
          p = parm[i].P1;
          parm[i].P1 = parm[i].P2;
          parm[i].P2 = p;

          parm[i].flip = 1;
        }

#ifdef DEBUG_MERGE
      for (i = 0; i < NTHREADS; i++)
        merge_thread(parm+i);
#else
      for (i = 1; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,merge_thread,parm+i);
      merge_thread(parm);

      for (i = 1; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);
#endif

      if (VERBOSE)
        { fprintf(stderr,"\r    Completed 100%%\n");
          fflush(stderr);
        }

      for (i = 0; i < NTHREADS; i++)
        { nhits += parm[i].nhits;
          g1len += parm[i].g1len;
          tseed += parm[i].tseed;
        }
    }
   
  free(cache);

  for (i = NTHREADS-1; i >= 1; i--)
    { Free_Kmer_Stream(parm[i].T1);
      Free_Kmer_Stream(parm[i].T2);
      Free_Post_List(parm[i].P1);
      Free_Post_List(parm[i].P2);
    }
  Free_Kmer_Stream(T1);
  Free_Kmer_Stream(T2);

  if (VERBOSE)
    { fprintf(stderr,"\n  Total seeds = %lld, ave. len = %.1f, seeds per genome position = %.1f\n",
                     nhits,(1.*tseed)/nhits,(1.*nhits)/g1len);
      fflush(stderr);
    }
}


static void self_adaptamer_merge(Kmer_Stream *T1, Post_List *P1)
{ SP         parm[NTHREADS];
#ifndef DEBUG_MERGE
  pthread_t  threads[NTHREADS];
#endif
  uint8     *cache;
  int64      nhits, g1len, tseed;
  int        i;

  if (VERBOSE)
    { fprintf(stderr,"  Starting adaptive seed merge\n");
      fflush(stderr);
    }

  { uint8       *ent;
#ifdef DEBUG_SPLIT
    char        *seq;
#endif
    int   t;
    int64 p;

    ent = Current_Entry(T1,NULL);
#ifdef DEBUG_SPLIT
    seq = Current_Kmer(T1,NULL);
#endif
    parm[0].pbeg = 0;
    for (t = 1; t < NTHREADS; t++)
      { p = (T1->nels * t) / NTHREADS;
        GoTo_Kmer_Index(T1,p);
        if (p >= T1->nels)
          parm[t].pbeg = 0xffff;
        else
          { ent = Current_Entry(T1,ent);
            parm[t].pbeg = (T1->cpre >> 8);
          }
      }
    for (t = 0; t < NTHREADS-1; t++)
      parm[t].pend = parm[t+1].pbeg;
    parm[NTHREADS-1].pend = 0xffff;
  }

  parm[0].T1 = T1;
  parm[0].P1 = P1;
  for (i = 1; i < NTHREADS; i++)
    { parm[i].T1 = Clone_Kmer_Stream(T1);
      parm[i].P1 = Clone_Post_List(P1);
    }

  cache = Malloc(NTHREADS*(P1->maxp+1)*KBYTE,"Allocating cache");
  if (cache == NULL)
    Clean_Exit(1);

  for (i = 0; i < NTHREADS; i++)
    { IOBuffer *nu, *cu;

      parm[i].tid   = i;
      parm[i].cache = cache + i * (P1->maxp+1) * KBYTE;
      parm[i].nunit = nu = N_Units + i * NPARTS;
      parm[i].cunit = cu = C_Units + i * NPARTS;
      bzero(nu[0].buck,sizeof(int64)*NCONTS);
      bzero(cu[0].buck,sizeof(int64)*NCONTS);
    }

#ifdef DEBUG_MERGE
  for (i = 0; i < NTHREADS; i++)
    self_merge_thread(parm+i);
#else
  for (i = 1; i < NTHREADS; i++)
    pthread_create(threads+i,NULL,self_merge_thread,parm+i);
  self_merge_thread(parm);

  for (i = 1; i < NTHREADS; i++)
    pthread_join(threads[i],NULL);
#endif

  if (VERBOSE)
    { fprintf(stderr,"\r    Completed 100%%\n");
      fflush(stderr);
    }
   
  free(cache);

  for (i = NTHREADS-1; i >= 1; i--)
    { Free_Kmer_Stream(parm[i].T1);
      Free_Post_List(parm[i].P1);
    }
  Free_Kmer_Stream(T1);

  nhits = g1len = tseed = 0;
  for (i = 0; i < NTHREADS; i++)
    { nhits += parm[i].nhits;
      g1len += parm[i].g1len;
      tseed += parm[i].tseed;
    }

  if (VERBOSE)
    { fprintf(stderr,"\n  Total seeds = %lld, ave. len = %.1f, seeds per G1 position = %.1f\n",
                     nhits,(1.*tseed)/nhits,(1.*nhits)/g1len);
      fflush(stderr);
    }
}


/***********************************************************************************************
 *
 *   PAIR RE-IMPORT AND SORT
 *
 **********************************************************************************************/

typedef struct
  { int       in;
    int       swide;
    int       comp;
    int       inum;
    GDB      *gdb1;
    GDB      *gdb2;
    int64    *buck;
    uint8    *buffer;
    uint8    *sarr;
    Range    *range;
  } RP;

static void *reimport_thread(void *args)
{ RP *parm = (RP *) args;
  int    swide  = parm->swide;
  int    in     = parm->in;
  int    comp   = parm->comp;
  uint8 *sarr   = parm->sarr;
  uint8 *bufr   = parm->buffer;
  int64 *buck   = parm->buck;

  int64  ipost, jpost, icont, jcont, band, anti;
  uint8 *_ipost = (uint8 *) (&ipost);
  uint8 *_jpost = (uint8 *) (&jpost);
  uint8 *_icont = (uint8 *) (&icont);
  uint8 *_jcont = (uint8 *) (&jcont);
  uint8 *_band  = (uint8 *) (&band);
  uint8 *_anti  = (uint8 *) (&anti);

  int    iamt;
  uint8 *x;
  int    iolen, iunit, lcp;
  int64  diag, flag, mask;
  uint8 *bend, *btop, *b;

  iolen = 2*NPARTS*1000000;
  iunit = IBYTE + JBYTE + 1;

  iamt = read(in,bufr,iolen);
  if (iamt < 0)
    { fprintf(stderr,"%s: IO read error for file %s/%s.%d.%c\n",
                     Prog_Name,SORT_PATH,PAIR_NAME,parm->inum,comp?'C':'N');
      Clean_Exit(1);
    }
  bend = bufr + iamt;

  if (bend-bufr < iolen)
    btop = bend;
  else
    btop = bend-iunit;
  b = bufr;

  ipost = 0;
  jpost = 0;
  jcont = 0;
  icont = 0;

  flag = (0x1ll << (8*JCONT-1));
  mask = flag-1;

  if (bend > bufr)

  while (1)
    { lcp = *b++;
      memcpy(_ipost,b,IPOST);
      b += IPOST;
      memcpy(_icont,b,ICONT);
      b += ICONT;
      memcpy(_jpost,b,JPOST);
      b += JPOST;
      memcpy(_jcont,b,JCONT);
      b += JCONT;
      jcont &= mask;

      x = sarr + swide * buck[icont]++;
      *x++ = lcp;
      if (comp)
        { diag = MAXDAG - (ipost + jpost);
          anti = AMXPOS - (ipost - jpost);
        }
      else
        { diag = BMXPOS + (ipost - jpost);
          anti = ipost + jpost;
        }
      band = (diag >> BUCK_SHIFT);
      *x++ = diag-(band<<BUCK_SHIFT);

      memcpy(x,_anti,DBYTE);
      x += DBYTE;
      memcpy(x,_band,DBYTE);
      x += DBYTE;
      memcpy(x,_jcont,JCONT);
      x += JCONT;

      if (b >= btop)
        { int ex = bend-b;
          memcpy(bufr,b,ex);
          bend = bufr+ex;
          iamt = read(in,bend,iolen-ex);
          if (iamt < 0)
            { fprintf(stderr,"%s: IO read error for file %s/%s.%d.%c\n",
                             Prog_Name,SORT_PATH,PAIR_NAME,parm->inum,comp?'C':'N');
              Clean_Exit(1);
            }
          bend += iamt;
          if (bend == bufr)
            break;
          if (bend-bufr < iolen)
            btop = bend;
          else
            btop = bend-iunit;
          b = bufr;
        }
    }

  close(in);

  return (NULL);
}

void print_seeds(uint8 *sarray, int swide, Range *range, int64 *panel,
                 GDB *gdb1, GDB *gdb2, int comp)
{ uint8 *e, *x;
  int    n, p;
  int    lcp, drm;

  int64  anti, dbuck, jcont;
  int64  ipost, jpost, diag;
  uint8 *_anti  = (uint8 *) (&anti);
  uint8 *_dbuck = (uint8 *) (&dbuck);
  uint8 *_jcont = (uint8 *) (&jcont);

  anti  = 0;
  dbuck = 0;
  jcont = 0;

  (void) gdb1;
  (void) gdb2;

  for (n = 0; n < NTHREADS; n++)
    { x = sarray + range[n].off;
      for (p = range[n].beg; p < range[n].end; p++)
        { e = x + panel[p];
          printf("\nContig %d:\n",p);
          while (x < e)
            { lcp = *x++;
              drm = *x++;
              memcpy(_anti,x,DBYTE);
              x += DBYTE;
              memcpy(_dbuck,x,DBYTE);
              x += DBYTE;
              memcpy(_jcont,x,JCONT);
              x += JCONT;

              diag = (dbuck<<BUCK_SHIFT)+drm;
              if (comp)
                { anti += gdb1->contigs[Perm1[p]].clen - AMXPOS;
                  diag += gdb1->contigs[Perm1[p]].clen - MAXDAG;
                }
              else
                diag -= BMXPOS;
             jpost = (anti - diag) >> 1;
             ipost = (anti + diag) >> 1;

             // if (jpost < 0 || jpost > gdb2->contigs[Perm2[jcont]].clen)
               // printf("J index out of bounds\n");
             // if (ipost < 0 || ipost > gdb1->contigs[Perm1[jcont]].clen)
               // printf("I index out of bounds\n");

             printf("  %10ld:  %5d %5lld: %8lld  %10lld x d=%10lld i=%10lld j=%10lld  (%2d)  %2d\n",
                     (x-sarray)/swide,p,jcont,dbuck,anti,diag,ipost,jpost,drm,lcp);
            }
        }
    }
}

#ifdef PRINT_SEQ

static void print_seq(char *seq, int b, int e)
{ static int toA[5] = { 'a', 'c', 'g', 't', '*' };
  int j;

  for (j = b; j < e; j++)
    printf("%c",toA[(int) seq[j]]);
  printf("\n");
}

#endif

static int entwine(Path *jpath, uint8 *jtrace, Path *kpath, uint8 *ktrace, int *where, int show)
{ int   ac, b2, y2, yp, ae;
  int   i, j, k;
  int   num, den, min;

  (void) show;

#ifdef DEBUG_ENTWINE
  if (show)
    printf("\n");
#endif

  *where = -1;

  y2 = jpath->bbpos;
  b2 = kpath->bbpos;
  j  = jpath->abpos/TSPACE;
  k  = kpath->abpos/TSPACE;

  ac = k*TSPACE;

  j = 1 + 2*(k-j);
  k = 1;
  for (i = 1; i < j; i += 2)
    y2 += jtrace[i];

  if (j == 1)
    yp = y2 + (jtrace[j] * (kpath->abpos - jpath->abpos)) / (ac+TSPACE - jpath->abpos);
  else
    yp = y2 + (jtrace[j] * (kpath->abpos - ac)) / TSPACE;

#ifdef DEBUG_ENTWINE
  if (show)
    printf("   @ %5d : %5d %5d = %4d\n",kpath->abpos,yp,b2,b2-yp);
#endif

  num = b2-yp;
  den = 1;
  min = num;

  ae = jpath->aepos;
  if (ae > kpath->aepos)
    ae = kpath->aepos;

  for (ac += TSPACE; ac < ae; ac += TSPACE)
    { y2 += jtrace[j];
      b2 += ktrace[k];
      j += 2;
      k += 2;

#ifdef DEBUG_ENTWINE
      if (show)
        printf("   @ %5d : %5d %5d = %4d\n",ac,y2,b2,b2-y2);
#endif

      i = b2-y2;
      num += i;
      den += 1;
      if (min < 0 && min < i)
        { if (i >= 0)
            min = 0; 
          else
            min = i;
        }
      else if (min > 0 && min > i)
        { if (i <= 0)
            min = 0;
          else
            min = i;
        }
      if (i == 0)
        *where = ac;
    }

  ac -= TSPACE;
  if (ae == jpath->aepos)
    { y2 = jpath->bepos;
      if (kpath->aepos >= ac)
        b2 += (ktrace[k] * (ae - ac)) / TSPACE;
      else
        b2 += (ktrace[k] * (ae - ac)) / (kpath->aepos - ac);
    }
  else
    { b2 = kpath->bepos;
      if (jpath->aepos >= ac)
        y2 += (jtrace[j] * (ae - ac)) / TSPACE;
      else
        y2 += (jtrace[j] * (ae - ac)) / (jpath->aepos - ac);
    }

#ifdef DEBUG_ENTWINE
  if (show)
    printf("   @ %5d : %5d %5d = %4d\n",ae,y2,b2,b2-y2);
#endif

  i = b2-y2;
  num += i;
  den += 1;
  if (min < 0 && min < i)
    { if (i >= 0)
        min = 0; 
      else
        min = i;
    }
  else if (min > 0 && min > i)
    { if (i <= 0)
        min = 0;
      else
        min = i;
    }

#ifdef DEBUG_ENTWINE
  if (show)
    { printf("MINIM = %d AVERAGE = %d",min,num/den);
      if (*where >= 0)
        printf(" WHERE = %d",*where);
      printf("\n");
    }
#endif

  (void) den;

  return (min);
}

typedef struct

  { int         tid;
    GDB        *gdb1, *gdb2;
    FILE       *ofile;
    FILE       *tfile;
    int64       nhits;
    int64       nlass;
    int64       nlive;
    int64       nlcov;
    int64       nmemo;
                            //  See align.h for doc on the following:
    Work_Data  *work;           //  work storage for alignment module
    Align_Spec *spec;           //  alignment spec
    Alignment   align;          //  alignment record
    Overlap     ovl;            //  overlap record

  } Contig_Bundle;

static int ALIGN_SORT(const void *l, const void *r)
{ Overlap *ol = *((Overlap **) l);
  Overlap *or = *((Overlap **) r);

  return (ol->path.abpos - or->path.abpos);
}

//  [beg,end) in the sorted array of width swide elements contain all the adaptive seeds between
//    the contigs in the parameter pair.  Look for seed chains in each pair of diagaonl buckets
//    of sufficient score, and when found search for an alignment, outputing it if found.

static void align_contigs(uint8 *beg, uint8 *end, int swide, int ctg1, int ctg2,
                          Contig_Bundle *pair)
{ Overlap    *ovl   = &(pair->ovl);
  int         comp  = (ovl->flags != 0);
#ifdef CALL_ALIGNER
  Work_Data  *work  = pair->work;
  Align_Spec *spec  = pair->spec;
  Alignment  *align = &(pair->align);
  Path       *path  = align->path;
  FILE       *ofile = pair->ofile;
  FILE       *tfile = pair->tfile;
#endif
#if defined(DEBUG_SEARCH) || defined(DEBUG_HIT)
  int         repgo;
#endif

  uint8 *b, *m, *e;

  int    alnMin;
  double alnRate;
  int64  nhit, nlas, nmem, nliv, ncov;
  int64  alen, blen, mlen;
  int64  aoffset, doffset;

  int    new, aux;
  int64  ndiag, cdiag;
  uint8 *_ndiag = (uint8 *) (&ndiag);

  int    self;

  int64  ipost, apost;
  uint8 *_ipost = (uint8 *) (&ipost);
  uint8 *_apost = (uint8 *) (&apost);

  ctg1 = Perm1[ctg1];
  ctg2 = Perm2[ctg2];

  if (pair->gdb1->contigs[ctg1].boff < 0 || pair->gdb2->contigs[ctg2].boff < 0)
    return;

  alnMin  = ALIGN_MIN - 50;
  alnRate = ALIGN_RATE + .05;

  ndiag = 0;
  ipost = 0;
  apost = 0;

  blen   = pair->gdb2->contigs[ctg2].clen;
  alen   = pair->gdb1->contigs[ctg1].clen;
  mlen   = alen+blen;

  nhit   = 0;
  nlas   = 0;
  nmem   = 0;
  nliv   = 0;
  ncov   = 0;

  if (SELF && ctg1 == ctg2 && !comp)
    self = 1;
  else
    self = 0;

  doffset = alen - MAXDAG;
  aoffset = alen - AMXPOS;

  //  Find segments b,m,e such that b > m and diag is cdiag for elements [b,m) and
  //    cdiag+1 for elements [m,e) (if m < e)
  //  If m == e (i.e. !aux) and b,m = m',e' of previous find (i.e. !new) then don't examine
  //    as the chain for this triple is subset of the chain for the previous triple.

  b = e = beg + (DBYTE+2);
  memcpy(_ndiag,e,DBYTE);
  cdiag = ndiag;
  while (ndiag == cdiag && e < end)
    { e += swide;
      memcpy(_ndiag,e,DBYTE);
    }
  new = 1;

#if defined(DEBUG_SEARCH) || defined(DEBUG_HIT)
  repgo = (ctg1 == 53) && (ctg2 == 11);
  if (repgo)
    printf("\n  Contig %d vs Contig %d\n",ctg1,ctg2);
#endif

  while (1)
    { m = e;
      aux = 0;
      while (ndiag == cdiag+1 && e < end)
        { e += swide;
          memcpy(_ndiag,e,DBYTE);
          aux = 1;
        }

      if (new || aux)
        { int    go, lcp, wch, mix, cov;
          int64  ahgh, alow, amid, alast;
          int64  anti, eant;
          int    dgmin, dgmax, dg;
          uint8 *s, *t;
     
#ifdef DEBUG_SEARCH
          if (repgo)
            { printf("Diag %lld",cdiag);
              if (aux)
                printf("+1");
              printf("\n");
            }
#endif

          //  Have triple b,m,e, b > m, to examine.  Capture the ipost-ordered merge of [b,m) and
          //    [m,e) in list[0..len) and process any above-threshold chains encountered while
          //    doing the merge.

          alast = -1;

          e -= DBYTE;
          m -= DBYTE;

          s = b-DBYTE;
          memcpy(_ipost,s,DBYTE);
          t = m; 
          if (aux)
            memcpy(_apost,t,DBYTE);
          else
            apost = MAX_INT64;

          dgmin = 2*BUCK_WIDTH;
          dgmax = 0;
          ahgh = -CHAIN_BREAK;
          if (apost < ipost)
            alow = apost;
          else
            alow = ipost;
          cov  = 0;
          go   = 1;
          mix  = 0;
          while (go)
            { if (apost < ipost)
                { lcp   = t[-2];
                  dg    = t[-1] + BUCK_WIDTH;
                  anti  = apost;
                  t += swide;
                  if (t >= e)
                    apost = MAX_INT64;
                  else
                    memcpy(_apost,t,DBYTE);
                  wch = 0x2;
                }
              else
                { lcp   = s[-2];
                  dg    = s[-1];
                  anti  = ipost;
                  s += swide;
                  if (s >= m)
                    { if (s > m)
                        go = 0;
                      else
                        ipost = MAX_INT64;
                    }
                  else
                    memcpy(_ipost,s,DBYTE);
                  wch = 0x1;
                }
              lcp <<= 1;

              if (anti < ahgh + CHAIN_BREAK)
                { int64 cps;

                  cps = anti + lcp;
                  if (cps > ahgh)
                    { if (anti >= ahgh)
                        cov += lcp;
                      else
                        cov += cps-ahgh;
                      ahgh = cps;
                    }
                  mix |= wch;
                  if (dg < dgmin)
                    dgmin = dg;
                  else if (dg > dgmax)
                    dgmax = dg;
                }
              else
                { if (cov >= CHAIN_MIN && (mix != 1 || new))

                    //  Have a chain that covers CHAIN_MIN or more anti-diagonals:
                    //    in the "tube" (alow..ahgh,dgmin..dgmax)
                    //    Search for local alignments within it.

                    { nhit += 1;
#ifdef DEBUG_SEARCH
                      if (repgo)
                        printf("                  Process\n");
#endif
#ifdef DEBUG_HIT
                      if (repgo)
                        { printf("Hit on bucket %lld",cdiag);
                          if (aux)
                            printf("+1");
                          printf(" Coverage = %d\n",cov);
                        }
#endif

                      //  Fetch contig sequences if not already loaded
#ifdef CALL_ALIGNER
                      if (ctg1 != ovl->aread)
                        { if (Get_Contig(pair->gdb1,ctg1,NUMERIC,align->aseq) == NULL)
                            Clean_Exit(1);
                          align->alen = alen;
                          ovl->aread  = ctg1;
                          if (comp)
                            Complement_Seq(align->aseq,align->alen);
#ifdef DEBUG_HIT
                          if (repgo)
                            printf("Loading A = %d%c\n",ctg1,comp?'c':'n');
                          fflush(stdout);
#endif
                        }
                      if (ctg2 != ovl->bread)
                        { if (Get_Contig(pair->gdb2,ctg2,NUMERIC,align->bseq) == NULL)
                            Clean_Exit(1);
                          align->blen = blen;
                          ovl->bread = ctg2;
#ifdef DEBUG_HIT
                          if (repgo)
                            printf("Loading B = %dn\n",ctg2);
                          fflush(stdout);
#endif
                        }
#endif

                      dgmin += (cdiag<<BUCK_SHIFT);
                      dgmax += (cdiag<<BUCK_SHIFT);
                      if (comp)
                        { dgmin += doffset;
                          dgmax += doffset;
                          alow  += aoffset;
                          ahgh  += aoffset;
                        }
                      else
                        { dgmin -= BMXPOS;
                          dgmax -= BMXPOS;
                        }

#ifdef DEBUG_HIT
                      if (repgo)
                        { printf("  Box:   Diag = %d:%d  Anti = %lld:%lld:%lld :: %d x %lld\n",
                                 dgmin,dgmax,alow,(alow+ahgh)>>1,ahgh,dgmax-dgmin,ahgh-alow);
                         
                          fflush(stdout);
                        }
#endif

                      if (ahgh > alast)
                        { int rlen;

                          if (alow < alast)
                            alow = alast;
                          ahgh -= BUCK_ANTI;
#ifdef CALL_ALIGNER
                          do {
                            amid = alow + BUCK_ANTI;
                            if (amid > ahgh)
                              { amid = ahgh;
                                if (amid + dgmin < 0)
                                  { dgmin = -amid;
                                    if (dgmin > dgmax)
                                      break;
                                  }
                              }

                            if (self)
                              { if (dgmin > 0)
                                  { if (Local_Alignment(align,work,spec,
                                                        dgmin,dgmax,amid,dgmin-1,-1) == NULL)
                                      Clean_Exit(1);
                                  }
                                else if (dgmax < 0)
                                  { if (Local_Alignment(align,work,spec,
                                                        dgmin,dgmax,amid,-1,-(dgmax+1)) == NULL)
                                      Clean_Exit(1);
                                  }
                                else
                                  path->abpos = path->aepos = 0;
                              }
                            else
                              { if (Local_Alignment(align,work,spec,dgmin,dgmax,amid,-1,-1) == NULL)
                                  Clean_Exit(1);
                              }

                            rlen = path->aepos - path->abpos;
                            if (rlen >= alnMin && alnRate*rlen >= path->diffs)
                              { Compress_TraceTo8(ovl,0);
                                if (fwrite(ovl,OVL_SIZE,1,tfile) != 1)
                                  { fprintf(stderr,
                                           "%s: Cannot write overlap gather file %s/%s.%d.las\n",
                                           Prog_Name,SORT_PATH,ALGN_PAIR,pair->tid);
                                    Clean_Exit(1);
                                  }
                                if (fwrite(ovl->path.trace,ovl->path.tlen,1,tfile) != 1)
                                  { fprintf(stderr,
                                            "%s: Cannot write overlap gather file %s/%s.%d.las\n",
                                            Prog_Name,SORT_PATH,ALGN_PAIR,pair->tid);
                                    Clean_Exit(1);
                                  }
                                nlas += 1;
                                nmem += path->tlen + OVL_SIZE;
                              }

#ifdef DEBUG_ALIGN
                            if (rlen >= ALIGN_MIN && ALIGN_RATE*rlen >= path->diffs)
                              { Decompress_TraceTo16(ovl);
                                printf("\nLocal %lld: %d-%d vs %d %d (%d)\n",nlas+1,
                                       path->abpos,path->aepos,path->bbpos,path->bepos,path->diffs);
                                if (comp)
                                  { Complement_Seq(align->aseq+(alen-path->aepos),
                                                   path->aepos-path->abpos);
                                    Complement_Seq(align->bseq+(blen-path->bepos),
                                                   path->bepos-path->bbpos);
                                  }
                                Compute_Trace_PTS(align,work,TSPACE,GREEDIEST);
                                Print_Alignment(stdout,align,work,4,100,10,0,8,0);
                                fflush(stdout);
                                if (comp)
                                  { Complement_Seq(align->aseq+(alen-path->aepos),
                                                   path->aepos-path->abpos);
                                    Complement_Seq(align->bseq+(blen-path->bepos),
                                                   path->bepos-path->bbpos);
                                  }
                              }
                            else
                              printf("Not found, len = %d\n",
                                     path->aepos-path->abpos);
#endif

                            if (comp)
                              eant = mlen-(path->abpos+path->bbpos);
                            else
                              eant = path->aepos+path->bepos;
    
#ifdef DEBUG_TUBE
                            { int64 bant;

                              if (repgo)
                                { if (comp)
                                    bant = mlen-(path->aepos+path->bepos);
                                  else
                                    bant = path->abpos+path->bbpos;
                                  if (eant < ahgh)
                                    printf("Did not reach top %lld %lld %lld\n",
                                           ahgh-eant,eant-bant,ahgh-alow);
                                  else if (bant > alow)
                                    printf("Did not reach bottom %lld %lld %lld\n",
                                           bant-alow,eant-bant,ahgh-alow);
                                  else
                                    printf("Good\n");
                                }
                            }
#endif

                            if (eant <= alow)
                              alow = amid;
                            else
                              alow = eant;
                          } while (alow < ahgh);
#endif
                          alast = alow;
                        }
#ifdef DEBUG_HIT
                      else if (repgo)
                        printf("BLOCKED %lld\n",alast);
#endif
                    }

#ifdef DEBUG_SEARCH
                  else if (repgo)
                    printf("                  Break\n");
#endif

                  if (go)
                    { cov  = lcp;
                      ahgh = anti + lcp;
                      mix  = wch;
                      alow = anti;
                      dgmin = dgmax = dg;
                    }
                }

#ifdef DEBUG_SEARCH
              if (go && repgo)
                { uint8 *n;

                  dg += (cdiag<<BUCK_SHIFT);
                  if (wch == 0x1)
                    n = s - swide;
                  else
                    n = t - swide;
                  if (comp)
                    printf("   %c %10ld: c %10lld x %10lld %2d %4d (%d)\n",
                           wch==0x1?'.':'+',(n-beg)/swide,anti+aoffset,dg+doffset,n[-2],cov,n[-1]);
                  else
                    printf("   %c %10ld: n %10lld x %10lld %2d %4d (%d)\n",
                           wch==0x1?'.':'+',(n-beg)/swide,anti,dg-BMXPOS,n[-2],cov,n[-1]);
                }
#endif
            }

          e += DBYTE;
          m += DBYTE;

          ipost = apost = 0;
        }

      if (e >= end) break;

      if (aux)
        { b = m;
          cdiag += 1;
          new = 0;
        }
      else
        { b = e;
          cdiag = ndiag;
          while (ndiag == cdiag && e < end)
            { e += swide;
              memcpy(_ndiag,e,DBYTE);
            }
          new = 1;
        }
    }

  //  Detect and remove redundant alignments

  if (nlas > 0)
    { void    *oblock;
      Overlap **perm;
      int      j, k, where, dist;
      Path     tpath;
      void    *tcopy;

      oblock = Malloc(nmem,"Allocating overlap block");
      perm   = Malloc(nlas*sizeof(Overlap *),"Allocating permutation array");
      if (oblock == NULL || perm == NULL)
        Clean_Exit(1);

      rewind(tfile);
      if (fread(oblock,nmem,1,tfile) != 1) 
        { fprintf(stderr,"\n%s: Cannot read overlap gather file %s/%s.%d.las\n",
                         Prog_Name,SORT_PATH,ALGN_PAIR,pair->tid);
          Clean_Exit(1);
        }

      { void *off;

        off = oblock;
        for (j = 0; j < nlas; j++)
          { perm[j] = (Overlap *) off;
            off += OVL_SIZE + ((Overlap *) off)->path.tlen;
          }
      }

      qsort(perm,nlas,sizeof(Overlap *),ALIGN_SORT);

#define ELIMINATED   0x4
#define OWNS_MEMORY  0x8
#define RESET_FLAGS  0x3

      for (j = nlas-1; j >= 0; j--)
        { Overlap *o  = perm[j];
          Path    *op = &(o->path);

          for (k = j+1; k < nlas; k++)
            { Overlap *w  = perm[k];
              Path    *wp = &(w->path);

              if (op->aepos <= wp->abpos)
                break;
              if (w->flags & ELIMINATED)
                continue;

              if (op->abpos == wp->abpos && op->bbpos == wp->bbpos)
                if (op->aepos == wp->aepos && op->bepos == wp->bepos)
                  { if (op->diffs < wp->aepos)    //  = endpoints
                      { w->flags |= ELIMINATED;
                        continue;
                      }
                    else
                      { o->flags |= ELIMINATED;
                        break;
                      }
                  }
                else
                  { if (op->aepos > wp->aepos)    //  a-interval of one contained in the other
                      { w->flags |= ELIMINATED;
                        // printf("  START - %d %d\n",op->aepos-wp->aepos,op->diffs-wp->diffs);
                        continue;
                      }
                    else
                      { o->flags |= ELIMINATED;
                        // printf("  START . %d %d\n",wp->aepos-op->aepos,wp->diffs-op->diffs);
                        break;
                      }
                  }
              else
                if (op->aepos == wp->aepos && op->bepos == wp->bepos) 
                  { if (op->abpos < wp->abpos)   //  a-interval of one contained in the other
                      { w->flags |= ELIMINATED;
                        // printf("  END - %d %d\n",wp->abpos-op->abpos,op->diffs-wp->diffs);
                        continue;
                      }
                    else
                      { o->flags |= ELIMINATED;
                        // printf("  END . %d %d\n",op->abpos-wp->abpos,wp->diffs-op->diffs);
                        break;
                      }
                  }
            }
        }

      for (j = nlas-1; j >= 0; j--)
        { Overlap *o  = perm[j];
          Path    *op = &(o->path);

          if (o->flags & ELIMINATED)
            continue;

          //  Both endpoints of o are distinct

          for (k = j+1; k < nlas; k++)
            { Overlap *w  = perm[k];
              Path    *wp = &(w->path);
              uint8   *otrace, *wtrace;

              if (op->aepos <= wp->abpos)   //  No further a-interval overlap
                break;
              if (w->flags & ELIMINATED)    //  Ignore w as already eliminated
                continue;
              if (op->bepos <= wp->bbpos || op->bbpos >= wp->bepos)  //  If b-intervals disjoint
                continue;

              if (o->flags & OWNS_MEMORY)
                otrace = (uint8 *) op->trace;
              else
                otrace = (uint8 *) (o+1);
              if (w->flags & OWNS_MEMORY)
                wtrace = (uint8 *) wp->trace;
              else
                wtrace = (uint8 *) (w+1);

              dist = entwine(op,otrace,wp,wtrace,&where,0);
              if (where != -1)   // The paths meet at a trace point given by where
                { uint8 *ntrace;
                  int    ocut, wcut;
                  int    d, h, g;

                  //  Fuse here

                  // printf("FUSE %d %d\n",(wp->abpos-op->abpos)+(wp->aepos-op->aepos),
                                       //  op->aepos-wp->abpos);
                  // printf(" %3d: %d-%d vs %d-%d\n            %d-%d vs %d-%d\n",
                         // where,o->path.abpos,o->path.aepos,w->path.abpos,w->path.aepos,
                         // o->path.bbpos,o->path.bepos,w->path.bbpos,w->path.bepos);
                  // dist = entwine(op,otrace,wp,wtrace,&where,1);

                  ocut = 2 * (((where-op->abpos)-1)/TSPACE+1);
                  wcut = 2 * (((where-wp->abpos)-1)/TSPACE+1);
                  op->tlen  = ocut + (wp->tlen-wcut);

                  ntrace = (uint8 *) Malloc(op->tlen,"Allocating new trace");
                  if (ntrace == NULL)
                    Clean_Exit(1);

                  d = 0;
                  h = 0;
                  for (g = 0; g < ocut; g += 2)
                    { d += (ntrace[h] = otrace[g]);
                      ntrace[h+1] = otrace[g+1];
                      h += 2;
                    }
                  for (g = wcut; g < wp->tlen; g += 2)
                    { d += (ntrace[h] = wtrace[g]);
                      ntrace[h+1] = wtrace[g+1];
                      h += 2;
                    }
                  
                  if (o->flags & OWNS_MEMORY)
                    free(otrace);
                  if (w->flags & OWNS_MEMORY)
                    free(wtrace);
                  op->diffs = d;
                  op->aepos = wp->aepos;
                  op->bepos = wp->bepos;
                  w->flags |= ELIMINATED;
                  o->flags |= OWNS_MEMORY;
                  op->trace = ntrace;
                  continue;
                }
#ifdef BOX_ELIM
              if (dist != 0)   // The paths do not cross
                { if ((op->aepos - op->abpos) + BOX_FUZZ >= wp->aepos - wp->abpos)
                    { if (wp->aepos <= op->aepos+BOX_FUZZ && wp->bbpos >= op->bbpos-BOX_FUZZ &&
                          wp->bepos <= op->bepos+BOX_FUZZ)
                        { w->flags |= ELIMINATED;
                          continue;
                        }
                    }
                  else
                    { if (op->aepos <= wp->aepos+BOX_FUZZ && op->bbpos >= wp->bbpos-BOX_FUZZ &&
                          op->bepos <= wp->bepos+BOX_FUZZ && op->abpos >= wp->abpos-BOX_FUZZ)
                        { o->flags |= ELIMINATED;
                          continue;
                        }
                    }
                }
#endif
              continue;

              // printf("OTHER\n");
              // printf(" %3d x %3d: %d-%d vs %d-%d\n            %d-%d vs %d-%d\n",
                     // j,k,o->path.abpos,o->path.aepos,w->path.abpos,w->path.aepos,
                     // o->path.bbpos,o->path.bepos,w->path.bbpos,w->path.bepos);
              dist = entwine(&(o->path),(uint8 *) (o+1),&(w->path),(uint8 *) (w+1),&where,1);

              if (op->abpos <= wp->abpos && op->aepos >= wp->aepos)
                { w->flags |= ELIMINATED;
                  printf("  CONTAIN - %d %d\n",
                         (wp->abpos-op->abpos)+(op->aepos-wp->aepos),op->diffs-wp->diffs);
                  // continue;
                }
              if (wp->abpos <= op->abpos && wp->aepos >= op->aepos)
                { o->flags |= ELIMINATED;
                  printf("  CONTAIN . %d %d\n",
                         (op->abpos-wp->abpos)+(wp->aepos-op->aepos),wp->diffs-op->diffs);
                  // break;
                }

              printf("\nAlign 1\n");
              tpath = o->path;
              align->path = &tpath;
              tpath.trace = tcopy = Malloc(sizeof(uint16)*tpath.tlen,"Trace");
              memcpy(tcopy,o+1,tpath.tlen);
              { uint16 *t16 = (uint16 *) tcopy;
                uint8  *t8  = (uint8  *) tcopy;
                int     nn;

                for (nn = tpath.tlen-1; nn >= 0; nn--)
                  t16[nn] = t8[nn];
              }
              Compute_Trace_PTS(align,work,TSPACE,GREEDIEST);
              Print_Reference(stdout,align,work,4,100,10,0,8,0);
              fflush(stdout);
              free(tcopy);

              printf("\nAlign 2\n");
              tpath = w->path;
              align->path = &tpath;
              tpath.trace = tcopy = Malloc(sizeof(uint16)*tpath.tlen,"Trace");
              memcpy(tcopy,w+1,tpath.tlen);
              { uint16 *t16 = (uint16 *) tcopy;
                uint8  *t8  = (uint8  *) tcopy;
                int     nn;

                for (nn = tpath.tlen-1; nn >= 0; nn--)
                  t16[nn] = t8[nn];
              }
              Compute_Trace_PTS(align,work,TSPACE,GREEDIEST);
              Print_Reference(stdout,align,work,4,100,10,0,8,0);
              fflush(stdout);
              free(tcopy);

              align->path = &(ovl->path);
            }
        }

      nmem = 0;
      for (j = 0; j < nlas; j++)
        { Overlap *o = perm[j];
          int      hasmem;

          if (o->flags & ELIMINATED)
            continue;
          hasmem = (o->flags & OWNS_MEMORY);
          o->flags &= RESET_FLAGS;
          if (fwrite( ((char *) o)+PTR_SIZE, EXO_SIZE, 1, ofile) != 1)
            { fprintf(stderr,"%s: Could not write to overlap block file %s/%s.%d.las\n",
                             Prog_Name,SORT_PATH,ALGN_UNIQ,pair->tid);
              Clean_Exit(1);
            }
          if (hasmem)
            { if (fwrite(o->path.trace, o->path.tlen, 1, ofile) != 1)
                { fprintf(stderr,"%s: Could not write to overlap block file %s/%s.%d.las\n",
                                 Prog_Name,SORT_PATH,ALGN_UNIQ,pair->tid);
                  Clean_Exit(1);
                }
              free(o->path.trace);
            }
          else
            { if (fwrite( (char *) (o+1), o->path.tlen, 1, ofile) != 1)
                { fprintf(stderr,"%s: Could not write to overlap block file %s/%s.%d.las\n",
                                 Prog_Name,SORT_PATH,ALGN_UNIQ,pair->tid);
                  Clean_Exit(1);
                }
            }
          nliv += 1;
          ncov += o->path.aepos - o->path.abpos;
          nmem += EXO_SIZE + o->path.tlen;
        }

      rewind (tfile);
      free(perm);
      free(oblock);
    }
  else
    nmem = 0;

  pair->nhits += nhit;
  pair->nlass += nlas;
  pair->nlive += nliv;
  pair->nlcov += ncov;
  pair->nmemo += nmem;
}


typedef struct
  { int       tid;
    int       swide;
    int       comp;
    int64    *panel;
    uint8    *sarr;
    Range    *range;
    GDB       gdb1;
    GDB       gdb2;
    FILE     *ofile;
    FILE     *tfile;
    int64     nhits;
    int64     nlass;
    int64     nlive;
    int64     nlcov;
    int64     nmemo;
    int       tmaxl;
  } TP;

static void *search_seeds(void *args)
{ TP *parm = (TP *) args;
  int      swide  = parm->swide;
  int      comp   = parm->comp;
  int64   *panel  = parm->panel;
  uint8   *sarray = parm->sarr;
  Range   *range  = parm->range;
  int      beg    = range->beg;
  int      end    = range->end;
  GDB     *gdb1   = &(parm->gdb1);
  GDB     *gdb2   = &(parm->gdb2);
  int      foffs  = swide-JCONT;
  FILE    *ofile  = parm->ofile;
  FILE    *tfile  = parm->tfile;

  int    icrnt;
  int64  jcrnt;
  uint8 *_jcrnt = (uint8 *) (&jcrnt);

  Contig_Bundle _pair, *pair = &_pair;

  uint8 *x, *e, *b;

  jcrnt = 0;

  pair->tid  = parm->tid;
  pair->gdb1 = gdb1;
  pair->gdb2 = gdb2;
  pair->align.aseq = New_Contig_Buffer(gdb1);
  pair->align.bseq = New_Contig_Buffer(gdb2);
  if (pair->align.bseq == NULL || pair->align.bseq == NULL)
    Clean_Exit(1);
  pair->align.path = &(pair->ovl.path);
  if (comp)
    { pair->ovl.flags   = COMP_FLAG;
      pair->align.flags = ACOMP_FLAG;
    }
  else
    { pair->ovl.flags   = 0;
      pair->align.flags = 0;
    }
  pair->ovl.aread = -1;
  pair->ovl.bread = -1;
  pair->work = New_Work_Data();
  pair->spec = New_Align_Spec(1.-ALIGN_RATE,100,gdb1->freq,0);
  if (pair->work == NULL || pair->spec == NULL)
    Clean_Exit(1);
  pair->ofile = ofile;
  pair->tfile = tfile;
  pair->nhits = 0;
  pair->nlass = 0;
  pair->nlive = 0;
  pair->nlcov = 0;
  pair->nmemo = 0;

  x = sarray + range->off;
  for (icrnt = beg; icrnt < end; icrnt++)
    { e = x + panel[icrnt];
      if (e > x)
        { memcpy(_jcrnt,x+foffs,JCONT);
          b = x;
          for (x += swide; x < e; x += swide)
            if (memcmp(_jcrnt,x+foffs,JCONT))
              { align_contigs(b,x,swide,icrnt,(int) jcrnt,pair);
                memcpy(_jcrnt,x+foffs,JCONT);
                b = x;
              }
          align_contigs(b,x,swide,icrnt,jcrnt,pair);
        }
    }

  Free_Align_Spec(pair->spec);
  Free_Work_Data(pair->work);
  free(pair->align.aseq-1);
  free(pair->align.bseq-1);

  parm->nhits += pair->nhits;
  parm->nlass += pair->nlass;
  parm->nlive += pair->nlive;
  parm->nlcov += pair->nlcov;
  parm->nmemo += pair->nmemo;
  return (NULL);
}

static int SORT_MAP(const void *x, const void *y)
{ Overlap *ol = *((Overlap **) x);
  Overlap *or = *((Overlap **) y);

  int      al, ar;
  int      bl, br;
  int      cl, cr;
  int      pl, pr;

  al = ol->aread;
  ar = or->aread;
  if (al != ar)
    return (al-ar);

  pl = ol->path.abpos;
  pr = or->path.abpos;
  if (pl != pr)
    return (pl-pr);

  bl = ol->bread;
  br = or->bread;
  if (bl != br)
    return (bl-br);

  cl = COMP(ol->flags);
  cr = COMP(or->flags);
  if (cl != cr)
    return (cl-cr);

  if (ol < or)
    return (-1);
  else if (ol > or)
    return (1);
  else
    return (0);
}

static void *la_sort(void *args)
{ TP *parm = (TP *) args;

  FILE *fid  = parm->ofile;
  int64 novl = parm->nlive;
  int64 size = parm->nmemo;
  
  void     *iblock, *off;
  Overlap **perm;
  int       j, tmaxl;
  
  if (novl == 0)
    return (NULL);

  iblock = Malloc(size+PTR_SIZE,"Allocating overlap block");
  perm   = Malloc(sizeof(Overlap *)*novl,"Allocating permutation array");
  if (iblock == NULL || perm == NULL)
    Clean_Exit(1);
  iblock += PTR_SIZE;

  rewind(fid);

  if (fread(iblock,size,1,fid) != 1)
    { fprintf(stderr,"\n%s: Cannot not read overlap block file %s/%s.%d.las\n",
                     Prog_Name,SORT_PATH,ALGN_UNIQ,parm->tid);
      Clean_Exit(1);
    }
  
  rewind(fid);
    
  off = iblock-PTR_SIZE;
  for (j = 0; j < novl; j++)
    { perm[j] = (Overlap *) off;
      off += EXO_SIZE + ((Overlap *) off)->path.tlen;
    }

  qsort(perm,novl,sizeof(Overlap *),SORT_MAP);

  tmaxl = 0;
  for (j = 0; j < novl; j++)
    { Overlap *o = perm[j];

      if (fwrite( ((void *) o)+PTR_SIZE, EXO_SIZE, 1, fid) != 1)
        { fprintf(stderr,"\n%s: Cannot not write sorted overlap block file %s/%s.%d.las\n",
                         Prog_Name,SORT_PATH,ALGN_UNIQ,parm->tid);
          Clean_Exit(1);
        }
      if (fwrite( (void *) (o+1), o->path.tlen, 1, fid) != 1)
        { fprintf(stderr,"\n%s: Cannot not write sorted overlap block file %s/%s.%d.las\n",
                         Prog_Name,SORT_PATH,ALGN_UNIQ,parm->tid);
          Clean_Exit(1);
        }
      if (o->path.tlen > tmaxl)
        tmaxl = o->path.tlen;
    }

  rewind(fid);

  free(perm);
  free(iblock-PTR_SIZE);

  parm->tmaxl = tmaxl;
  return (NULL);
}

#define MEMORY 4000   // in Mb

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
    void   *block;
    void   *ptr;
    void   *top;
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

static int la_merge(TP *parm)
{ IO_block *in;
  int64     bsize;
  char     *block;
  int       i, c;
  Overlap **heap;
  int       hsize, tmaxl;
  Overlap  *ovls;
  int64     totl;
  OneFile  *of;
  int64    *trace64;

  //  Base level merge: Open all the input files and initialize their buffers

  bsize  = (MEMORY*1000000ll)/NTHREADS;
  block  = (char *) Malloc(bsize*NTHREADS+PTR_SIZE,"Allocating LAmerge blocks");
  in     = (IO_block *) Malloc(sizeof(IO_block)*NTHREADS,"Allocating LAmerge IO-reacords");
  if (block == NULL || in == NULL)
    return (1);
  block += PTR_SIZE;

  tmaxl = 0;
  totl = 0;
  for (c = 0; c < NTHREADS; c++)
    { void *iblock;

      in[c].stream = parm[c].ofile;
      in[c].block  = iblock = block+c*bsize;
      in[c].ptr    = iblock;
      in[c].top    = iblock + fread(iblock,1,bsize,parm[c].ofile);
      in[c].count  = 0;
      totl += parm[c].nlive;
      if (parm[c].tmaxl > tmaxl)
        tmaxl = parm[c].tmaxl;
    }

  //  Initialize the heap

  heap = (Overlap **) Malloc(sizeof(Overlap *)*(NTHREADS+1),"Allocating heap");
  ovls = (Overlap *) Malloc(sizeof(Overlap)*NTHREADS,"Allocating heap");
  trace64 = (int64 *) Malloc(sizeof(int64)*(tmaxl/2),"Allocating int64 trace vector");
  if (heap == NULL || ovls == NULL)
    return (1);

  hsize = 0;
  for (i = 0; i < NTHREADS; i++)
    { if (in[i].ptr < in[i].top)
        { ovls[i]     = *((Overlap *) (in[i].ptr - PTR_SIZE));
          in[i].ptr  += EXO_SIZE;
          hsize      += 1;
          heap[hsize] = ovls + i;
        }
    }

  if (hsize > 3)
    for (i = hsize/2; i > 1; i--)
      maheap(i,heap,hsize);

  //  Open the output file buffer and write (novl,tspace) header

  { char *db1_name;
    char *db2_name;
    char *cpath;

    if (TYPE1 < IS_GDB && !KEEP)
      db1_name = Strdup(SPATH1,"db1_name");
    else
      db1_name = Strdup(Catenate(PATH1,"/",ROOT1,GEXTN1),"db1_name");
    if (SELF)
      db2_name = NULL;
    else
      { if (TYPE2 < IS_GDB && !KEEP)
          db2_name = Strdup(SPATH2, "db2_name");
        else
          db2_name = Strdup(Catenate(PATH2,"/",ROOT2,GEXTN2), "db2_name");
      }
    cpath = getcwd(NULL,0);

    of = open_Aln_Write(Catenate(ONE_PATH,"/",ONE_ROOT,".1aln"), 1,
                        Prog_Name, VERSION, Command_Line,
                        TSPACE, db1_name, db2_name, cpath);

    free(cpath);
    if (db2_name != NULL)
      free(db2_name);
    free(db1_name);
  }

  //  While the heap is not empty do

  while (hsize > 0)
    { Overlap  *ov;
      IO_block *src;
      int64     tsize, span;

      maheap(1,heap,hsize);

      ov  = heap[1];
      src = in + (ov - ovls);

      src->count += 1;

      tsize = ov->path.tlen;
      span  = EXO_SIZE + tsize;
      if (src->ptr + span > src->top)
        ovl_reload(src,bsize);

      Write_Aln_Overlap (of, ov);
      Write_Aln_Trace (of, src->ptr, tsize, trace64);

      src->ptr += tsize;
      if (src->ptr >= src->top)
        { heap[1] = heap[hsize];
          hsize  -= 1;
          continue;
        }
      *ov       = *((Overlap *) (src->ptr - PTR_SIZE));
      src->ptr += EXO_SIZE;
    }

  oneFileClose(of);

  for (i = 0; i < NTHREADS; i++)
    fclose(parm[i].ofile);

  for (i = 0; i < NTHREADS; i++)
    totl -= in[i].count;
  if (totl != 0)
    { fprintf(stderr,"%s: Did not write all records to %s/%s.1aln (%lld)\n",
                     Prog_Name,ONE_PATH,ONE_ROOT,totl);
      return (1);
    }

  free(trace64);
  free(ovls);
  free(heap);
  free(in);
  free(block-PTR_SIZE);

  return (0);
}

static void pair_sort_search(GDB *gdb1, GDB *gdb2)
{ uint8 *sarray;
  int    swide;
  int64  nels;

  RP     rarm[NTHREADS];
  TP     tarm[NTHREADS];
  pthread_t threads[NTHREADS];
  int64    *panel;
  Range     range[NTHREADS];

  IOBuffer *unit[2], *nu;
  int       nused;
  int       i, p, j, u;

  if (VERBOSE)
    { fprintf(stderr,"\n  Starting seed sort and alignment search, %d parts\n",2*NPARTS);
      fflush(stderr);
    }

  unit[0] = N_Units;
  unit[1] = C_Units;

  { int64 cum, nelmax;

    nelmax = 0;
    for (u = 0; u < 2; u++)
      { cum = 0;
        nu = unit[u];

        for (j = 0; j < NCONTS; j++)
          { for (i = 0; i < NTHREADS; i++)
              { cum += nu[i].buck[j];
                nu[i].buck[j] = cum;
              }
            if (j+1 == NCONTS || Select[j] != Select[j+1])
              { if (cum > nelmax)
                  nelmax = cum;
                cum = 0;
              }
          }

        for (j = NCONTS-1; j >= 0; j--)
          { for (i = NTHREADS-1; i >= 1; i--)
              nu[i].buck[j] = nu[i-1].buck[j];
            if (j == 0 || Select[j] != Select[j-1])
              nu[0].buck[j] = 0;
            else
              nu[0].buck[j] = nu[NTHREADS-1].buck[j-1];
          }
      }

    swide  = 2*DBYTE + JCONT + 2;
    sarray = Malloc((nelmax+1)*swide,"Sort Array");
    panel  = Malloc(NCONTS*sizeof(int64),"Bucket Array");
    if (sarray == NULL || panel == NULL)
      Clean_Exit(1);
  }

#ifdef LOAD_SEQS
  Load_Sequences(gdb2,COMPRESSED);
#endif

  for (p = 0; p < NTHREADS; p++)
    { rarm[p].swide  = swide;
      rarm[p].sarr   = sarray;
      rarm[p].buffer = N_Units[p].bufr;   //  NB: Units have been transposed
      rarm[p].range  = range+p;
      rarm[p].gdb1   = gdb1;
      rarm[p].gdb2   = gdb2;

      tarm[p].tid    = p;
      tarm[p].swide  = swide;
      tarm[p].sarr   = sarray;
      tarm[p].panel  = panel;
      tarm[p].range  = range+p;

      tarm[p].gdb1   = *gdb1;
      tarm[p].gdb2   = *gdb2;
      if (p > 0)
        { if (gdb1->seqstate == EXTERNAL)
            { tarm[p].gdb1.seqs = fopen(gdb1->seqpath,"r");
              if (tarm[p].gdb1.seqs == NULL)
                { fprintf(stderr,"%s: Cannot open another copy of GDB\n",Prog_Name);
                  Clean_Exit(1);
                }
            }
          if (gdb2->seqstate == EXTERNAL)
            { tarm[p].gdb2.seqs = fopen(gdb2->seqpath,"r");
              if (tarm[p].gdb2.seqs == NULL)
                { fprintf(stderr,"%s: Cannot open another copy of GDB\n",Prog_Name);
                  Clean_Exit(1);
                }
            }
        }

      tarm[p].nhits = 0;
      tarm[p].nlass = 0;
      tarm[p].nlive = 0;
      tarm[p].nlcov = 0;
      tarm[p].nmemo = 0;

      tarm[p].ofile = fopen(Catenate(SORT_PATH,"/",ALGN_UNIQ,Numbered_Suffix(".",p,".las")),"w+");
      if (tarm[p].ofile == NULL)
        { fprintf(stderr,"%s: Cannot open %s/%s.%d.las for writing\n",
                         Prog_Name,SORT_PATH,ALGN_UNIQ,p);
          Clean_Exit(1);
        }
      unlink(Catenate(SORT_PATH,"/",ALGN_UNIQ,Numbered_Suffix(".",p,".las")));

      tarm[p].tfile = fopen(Catenate(SORT_PATH,"/",ALGN_PAIR,Numbered_Suffix(".",p,".las")),"w+");
      if (tarm[p].tfile == NULL)
        { fprintf(stderr,"%s: Cannot open %s/%s.%d.las for reading & writing\n",
                         Prog_Name,SORT_PATH,ALGN_PAIR,p);
          Clean_Exit(1);
        }
      unlink(Catenate(SORT_PATH,"/",ALGN_PAIR,Numbered_Suffix(".",p,".las")));
    }

  for (u = 0; u < 2; u++)
   for (i = 0; i < NPARTS; i++)
    { nu = unit[u] + i*NTHREADS;

      if (VERBOSE)
        { fprintf(stderr,"\r    Loading seeds for part %d  ",u*NPARTS+i+1);
          fflush(stderr);
        }

      for (p = 0; p < NTHREADS; p++)
        { rarm[p].in = nu[p].file;
          lseek(nu[p].file,0,SEEK_SET);
          rarm[p].buck = nu[p].buck;
          rarm[p].comp = u;
          rarm[p].inum = nu[p].inum;
        }

#ifdef DEBUG_SORT
      for (p = 0; p < NTHREADS; p++)
        reimport_thread(rarm+p);
#else
      for (p = 1; p < NTHREADS; p++)
        pthread_create(threads+p,NULL,reimport_thread,rarm+p);
      reimport_thread(rarm);
      for (p = 1; p < NTHREADS; p++)
        pthread_join(threads[p],NULL);
#endif

#ifdef DEBUG_SORT
      for (p = 0; p < NTHREADS; p++)
        printf("  %d",nu[p].file);
      printf("\n");
      for (j = 0; j < NCONTS; j++)
        { printf(" %4d:",j);
          for (p = 0; p < NTHREADS; p++)
            printf(" %10lld",unit[u][p].buck[j]);
          printf("\n");
        }
     fflush(stdout);
#endif

      { int64 prev, next;

        bzero(panel,sizeof(int64)*NCONTS);
        prev = 0;
        next = 0;
        for (j = IDBsplit[i]; j < IDBsplit[i+1]; j++)
          { next = nu[NTHREADS-1].buck[j];
            panel[j] = (next - prev)*swide;
            prev = next;
          }
        nels = next;

#ifdef DEBUG_SORT
        for (p = 0; p < NCONTS; p++)
          if (panel[p] > 0)
            printf(" %2d(%2d): %10lld %10lld\n",p,Perm1[p],panel[p],panel[p]/swide);
#endif

        if (VERBOSE)
          { fprintf(stderr,"\r    Sorting seeds for part %d  ",u*NPARTS+i+1);
            fflush(stderr);
          }

        nused = rmsd_sort(sarray,nels,swide,swide-2,NCONTS,panel,NTHREADS,range);

#ifdef DEBUG_SORT
        print_seeds(sarray,swide,range,panel,gdb1,gdb2,u);
#endif
      }

      if (VERBOSE)
        { fprintf(stderr,"\r    Searching seeds for part %d",u*NPARTS+i+1);
          fflush(stderr);
        }

      for (p = 0; p < nused; p++)
        tarm[p].comp = u;

#if defined(DEBUG_SORT) || defined(DEBUG_SEARCH) || defined(DEBUG_HIT) || defined(DEBUG_ALIGN)
      for (p = 0; p < nused; p++)
        search_seeds(tarm+p);
#else
      for (p = 1; p < nused; p++)
        pthread_create(threads+p,NULL,search_seeds,tarm+p);
      search_seeds(tarm);
      for (p = 1; p < nused; p++)
        pthread_join(threads[p],NULL);
#endif
    }

  free(panel);
  free(sarray);
  for (p = 0; p < NTHREADS; p++)
    fclose(tarm[p].tfile);
  for (p = 1; p < NTHREADS; p++)
    { if (gdb2->seqstate == EXTERNAL)
        fclose(tarm[p].gdb2.seqs);
      if (gdb1->seqstate == EXTERNAL)
        fclose(tarm[p].gdb1.seqs);
    }

  if (VERBOSE)
    { int64 nhit, nlas, nliv, ncov;

      fprintf(stderr,"\r    Done                        \n");

      nhit = nlas = nliv = ncov = 0;
      for (p = 0; p < NTHREADS; p++)
        { nhit += tarm[p].nhits;
          nlas += tarm[p].nlass;
          nliv += tarm[p].nlive;
          ncov += tarm[p].nlcov;
        }
      if (nliv == 0)
        fprintf(stderr,
           "\n  Total hits over %dbp = %lld, %lld aln's, 0 non-redundant aln's of ave len 0\n",
                       CHAIN_MIN/2,nhit,nlas);
      else
        fprintf(stderr,
          "\n  Total hits over %dbp = %lld, %lld aln's, %lld non-redundant aln's of ave len %lld\n",
                       CHAIN_MIN/2,nhit,nlas,nliv,ncov/nliv);
      fflush(stderr);
    }

  if (VERBOSE)
    { fprintf(stderr,"\n  Sorting and merging alignments\n");
      fflush(stderr);
    }

#ifdef DEBUG_LASORT
  for (p = 0; p < NTHREADS; p++)
    la_sort(tarm+p);
#else
  for (p = 1; p < NTHREADS; p++)
    pthread_create(threads+p,NULL,la_sort,tarm+p);
  la_sort(tarm);
  for (p = 1; p < NTHREADS; p++)
    pthread_join(threads[p],NULL);
#endif

  if (la_merge(tarm))
    Clean_Exit(1);
}

static void short_GDB_fix(GDB *gdb)
{ int i;

  if (gdb->ncontig >= NTHREADS)
    return;

  //  Add additional conitgs of length KMER that are the first bit of the 0th read/contig
  //    Mark as "fake" with -1 in origin field.

  gdb->contigs = Realloc(gdb->contigs,sizeof(GDB_CONTIG)*NTHREADS,"Reallocating GDB contig vector");
  for (i = gdb->ncontig; i < NTHREADS; i++)
    { gdb->contigs[i] = gdb->contigs[0];
      gdb->contigs[i].clen = KMER;
      gdb->contigs[i].boff = -1;
    }
  gdb->seqtot += (NTHREADS-gdb->ncontig)*KMER;
  if (gdb->maxctg < KMER)
    gdb->maxctg = KMER;
  gdb->ncontig = NTHREADS;
}

int main(int argc, char *argv[])
{ Kmer_Stream *T1, *T2;
  Post_List   *P1, *P2;
  GDB _gdb1, *gdb1 = &_gdb1;
  GDB _gdb2, *gdb2 = &_gdb2;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;
    FILE  *test;

    ARG_INIT("FastGA");

    FREQ = 10;
    CHAIN_BREAK = 2000;   //  2x in anti-diagonal space
    CHAIN_MIN   =  170;
    ALIGN_MIN   =  100;
    ALIGN_RATE  = .3;
    SORT_PATH = getenv("TMPDIR");
    if (SORT_PATH == NULL)
      SORT_PATH = ".";
    NTHREADS    = 8;

    OUT_TYPE    = 0;
    OUT_OPT     = 0;
    ONE_PATH    = NULL;
    ONE_ROOT    = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vkS")
            break;
          case '1':
            if (strncmp(argv[i]+1,"1:",2) == 0)
              { OUT_TYPE = 2;
                ONE_PATH = PathTo(argv[i]+3);
                ONE_ROOT = Root(argv[i]+3,".1aln");
                test = fopen(Catenate(ONE_PATH,"/",ONE_ROOT,".1aln"),"w");
                if (test == NULL)
                  { fprintf(stderr,"%s: Cannot open %s/%s.1aln for output\n",
                                   Prog_Name,ONE_PATH,ONE_ROOT);
                    exit (1);
                  }
                fclose(test);
                break;
              }
            fprintf(stderr,"%s: Do not recognize option %s\n",Prog_Name,argv[i]);
            exit (1);
          case 'c':
            ARG_NON_NEGATIVE(CHAIN_MIN,"minimum seed cover");
            CHAIN_MIN <<= 1;
            break;
          case 'f':
            ARG_NON_NEGATIVE(FREQ,"maximum seed frequency");
            break;
          case 'i':
            ARG_REAL(ALIGN_RATE);
            if (ALIGN_RATE < .55 || ALIGN_RATE >= 1.)
              { fprintf(stderr,"%s: '-i' minimum alignment similarity must be in [0.55,1.0)\n",
                               Prog_Name);
                exit (1);
              }
            ALIGN_RATE = 1.-ALIGN_RATE;
            break;
          case 'l':
            ARG_NON_NEGATIVE(ALIGN_MIN,"minimum alignment length");
            break;
          case 'p':
            if (strncmp(argv[i]+1,"paf",3) == 0)
              { OUT_TYPE = 0;
                if (argv[1][4] == '\0')
                  { OUT_OPT = 0;
                    break;
                  }
                else if (strcmp(argv[i]+4,"m") == 0)
                  { OUT_OPT = 1;
                    break;
                  }
                else if (strcmp(argv[i]+4,"x") == 0)
                  { OUT_OPT = 2;
                    break;
                  }
              }
            else if (strcmp(argv[i]+1,"psl") == 0)
              { OUT_TYPE = 1;
                break;
              }
            fprintf(stderr,"%s: Do not recognize option %s\n",Prog_Name,argv[i]);
            exit (1);
          case 's':
            ARG_NON_NEGATIVE(CHAIN_BREAK,"seed chain break threshold");
            CHAIN_BREAK <<= 1;
            break;
          case 'P':
            SORT_PATH = argv[i]+2;
            break;
          case 'T':
            ARG_NON_NEGATIVE(NTHREADS,"number of threads to use");
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE   = flags['v'];
    KEEP      = flags['k'];
    SYMMETRIC = flags['S'];

    if (argc != 3 && argc != 2)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"\n");
        fprintf(stderr,"         <format> = -paf[mx] | -psl | -1:<align:path>[.1aln]\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"         <precursor> = .gix | .1gdb | <fa_extn> | <1_extn>\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"             <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"             <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -k: Keep any generated .1gdb's and .gix's.\n");
        // fprintf(stderr,"      -S: Use both A- and B adaptamer seeds\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        fprintf(stderr,"      -P: Directory to use for temporary files.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -paf: Stream PAF output\n");
        fprintf(stderr,"        -pafx: Stream PAF output with CIGAR sring with X's\n");
        fprintf(stderr,"        -pafm: Stream PAF output with CIGAR sring with ='s\n");
        fprintf(stderr,"      -psl: Stream PSL output\n");
        fprintf(stderr,"      -1: Generate 1-code output to specified file\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -f: adaptive seed count cutoff\n");
        fprintf(stderr,"      -c: minimum seed chain coverage in both genomes\n");
        fprintf(stderr,"      -s: threshold for starting a new seed chain\n");
        fprintf(stderr,"      -l: minimum alignment length\n");
        fprintf(stderr,"      -i: minimum alignment identity\n");
        fprintf(stderr,"\n");
        exit (1);
      }

    if (FREQ > 255)
      { fprintf(stderr,"%s: The maximum allowable frequency cutoff is 255\n",Prog_Name);
        exit (1);
      }
  }

  //  Parse source names and make precursors if necessary

  { char *p;
    FILE *input;
    char *tpath1;
    char *tpath2;

    PATH1 = argv[1];
    p = rindex(PATH1,'/');
    if (p == NULL)
      { ROOT1 = PATH1;
        PATH1 = ".";
      }
    else
      { *p++ = '\0';
        ROOT1 = p;
      }

    input = fopen(Catenate(PATH1,"/",ROOT1,".gix"),"r");
    if (input != NULL)
      { fclose(input);
        TYPE1 = IS_GDB+1;
      }
    else if (strcmp(ROOT1+(strlen(ROOT1)-4),".gix") == 0)
      { ROOT1[strlen(ROOT1)-4] = '\0';
        TYPE1 = IS_GDB+1;
      }
    else
      { TYPE1 = Get_GDB_Paths(Catenate(PATH1,"/",ROOT1,""),NULL,&SPATH1,&tpath1,0);
        ROOT1 = Root(tpath1,NULL);

        input = fopen(Catenate(PATH1,"/",ROOT1,".gix"),"r");
        if (input != NULL)
          { TYPE1 = IS_GDB+1;
            fclose(input);
          }
        else
          { input = fopen(Catenate(PATH1,"/",ROOT1,".1gdb"),"r");
            if (input != NULL)
              { fclose(input);
                TYPE1 = IS_GDB;
              }
            else
              { input = fopen(Catenate(PATH1,"/",ROOT1,".gdb"),"r");
                if (input != NULL)
                  { fclose(input);
                    TYPE1 = IS_GDB;
                  }
              }
          }
      }

    TYPE2 = IS_GDB+1;
    if (argc == 3)
      { PATH2 = argv[2];
        p = rindex(PATH2,'/');
        if (p == NULL)
          { ROOT2 = PATH2;
            PATH2 = ".";
          }
        else
          { *p++ = '\0';
            ROOT2 = p;
          }

        input = fopen(Catenate(PATH2,"/",ROOT2,".gix"),"r");
        if (input != NULL)
          { fclose(input);
            TYPE2 = IS_GDB+1;
          }
        else if (strcmp(ROOT2+(strlen(ROOT2)-4),".gix") == 0)
          { ROOT2[strlen(ROOT2)-4] = '\0';
            TYPE2 = IS_GDB+1;
          }
        else
          { TYPE2 = Get_GDB_Paths(Catenate(PATH2,"/",ROOT2,""),NULL,&SPATH2,&tpath2,0);
            ROOT2 = Root(tpath2,NULL);

            input = fopen(Catenate(PATH2,"/",ROOT2,".gix"),"r");
            if (input != NULL)
              { TYPE2 = IS_GDB+1;
                fclose(input);
              }
            else
              { input = fopen(Catenate(PATH2,"/",ROOT2,".1gdb"),"r");
                if (input != NULL)
                  { fclose(input);
                    TYPE2 = IS_GDB;
                  }
                else
                  { input = fopen(Catenate(PATH2,"/",ROOT2,".gdb"),"r");
                    if (input != NULL)
                      { fclose(input);
                        TYPE2 = IS_GDB;
                      }
                  }
              }
          }
      }

    if (TYPE1 <= IS_GDB)
      { char *command;

        command = Malloc(strlen(SPATH1)+100,"Allocating command string");
        if (command == NULL)
          exit (1);

        if (TYPE1 < IS_GDB)
          { sprintf(command,"FAtoGDB%s %s",VERBOSE?" -v":"",SPATH1);
            if (system(command) != 0)
              { fprintf(stderr,"\n%s: Call to FAtoGDB failed\n",Prog_Name);
                exit (1);
              }
          }

        sprintf(command,"GIXmake%s -T%d -P%s -f%d %s",
                        VERBOSE?" -v":"",NTHREADS,SORT_PATH,FREQ,tpath1);
        if (system(command) != 0)
          { fprintf(stderr,"\n%s: Call to GIXmake failed\n",Prog_Name);
            Clean_Exit(1);
          }

        free(command);
        free(tpath1);
      }

    if (TYPE2 <= IS_GDB)
      { char *command;

        command = Malloc(strlen(SPATH2)+100,"Allocating command string");
        if (command == NULL)
          exit (1);

        if (TYPE2 < IS_GDB)
          { sprintf(command,"FAtoGDB%s %s",VERBOSE?" -v":"",SPATH2);
            if (system(command) != 0)
              { fprintf(stderr,"\n%s: Call to FAtoGDB failed\n",Prog_Name);
                Clean_Exit(1);
              }
          }
        sprintf(command,"GIXmake%s -T%d -P%s -f%d %s",
                        VERBOSE?" -v":"",NTHREADS,SORT_PATH,FREQ,tpath2);
        if (system(command) != 0)
          { fprintf(stderr,"\n%s: Call to GIXmake failed\n",Prog_Name);
            Clean_Exit(1);
          }
        free(command);
        free(tpath2);
      }
  }

  if (VERBOSE)
    StartTime();

  //  Get full path string for sorting subdirectory (in variable SORT_PATH)

  { char  *cpath, *spath;
    DIR   *dirp;

    if (SORT_PATH[0] != '/')
      { cpath = getcwd(NULL,0);
        if (SORT_PATH[0] == '.')
          { if (SORT_PATH[1] == '/')
              spath = Catenate(cpath,SORT_PATH+1,"","");
            else if (SORT_PATH[1] == '\0')
              spath = cpath;
            else
              { fprintf(stderr,"\n%s: -P option: . not followed by /\n",Prog_Name);
                Clean_Exit(1);
                exit (1);
              }
          }
        else
          spath = Catenate(cpath,"/",SORT_PATH,"");
        SORT_PATH = Strdup(spath,"Allocating path");
        free(cpath);
      }
    else
      SORT_PATH = Strdup(SORT_PATH,"Allocating path");

    if ((dirp = opendir(SORT_PATH)) == NULL)
      { fprintf(stderr,"\n%s: -P option: cannot open directory %s\n",Prog_Name,SORT_PATH);
        Clean_Exit(1);
      }
    closedir(dirp);
  }

  SELF = (argc == 2);

  T1 = Open_Kmer_Stream(Catenate(PATH1,"/",ROOT1,".gix"));
  if (SELF)
    T2 = T1;
  else
    T2 = Open_Kmer_Stream(Catenate(PATH2,"/",ROOT2,".gix"));
  if (T1 == NULL)
    { fprintf(stderr,"%s: Cannot find genome index for %s/%s.gix\n",Prog_Name,PATH1,ROOT1);
      Clean_Exit(1);
    }
  if (T2 == NULL)
    { fprintf(stderr,"%s: Cannot find genome index for %s/%s.gix\n",Prog_Name,PATH2,ROOT2);
      Clean_Exit(1);
    }
  
  P1 = Open_Post_List(Catenate(PATH1,"/",ROOT1,".gix"));
  if (SELF)
    P2 = P1;
  else
    P2 = Open_Post_List(Catenate(PATH2,"/",ROOT2,".gix"));
  if (P1 == NULL)
    { fprintf(stderr,"%s: Cannot find genome index for %s/%s.gix\n",Prog_Name,PATH1,ROOT1);
      Clean_Exit(1);
    }
  if (P2 == NULL)
    { fprintf(stderr,"%s: Cannot find genome index for %s/%s.gix\n",Prog_Name,PATH2,ROOT2);
      Clean_Exit(1);
    }

  Perm1  = P1->perm;
  Perm2  = P2->perm;
  KMER   = T1->kmer;

  { FILE *file;
    char *fname;

    GEXTN1 = ".gdb";
    fname = Catenate(PATH1,"/",ROOT1,GEXTN1);
    if ((file = fopen(fname,"r")) == NULL)
      GEXTN1 = ".1gdb";
    else
      fclose(file);

    if (Read_GDB(gdb1,Catenate(PATH1,"/",ROOT1,GEXTN1)) < 0)
      Clean_Exit(1);
    short_GDB_fix(gdb1);

    if (SELF)
      gdb2 = gdb1;
    else
      { GEXTN2 = ".gdb";
        fname = Catenate(PATH1,"/",ROOT1,GEXTN2);
        if ((file = fopen(fname,"r")) == NULL)
          GEXTN2 = ".1gdb";
        else
          fclose(file);

        if (Read_GDB(gdb2,Catenate(PATH2,"/",ROOT2,GEXTN2)) < 0)
          Clean_Exit(1);
        short_GDB_fix(gdb2);
      }
  }

  if (OUT_TYPE != 2)
    { ONE_ROOT = Strdup(Numbered_Suffix("_oaln.",getpid(),""),"Allocating temp name");
      ONE_PATH = SORT_PATH;
      if (ONE_ROOT == NULL)
        Clean_Exit(1);
    }
  ALGN_UNIQ = Strdup(Numbered_Suffix("_uniq.",getpid(),""),"Allocating temp name");
  PAIR_NAME = Strdup(Numbered_Suffix("_pair.",getpid(),""),"Allocating temp name");
  ALGN_PAIR = Strdup(Numbered_Suffix("_algn.",getpid(),""),"Allocating temp name");
  if (ALGN_UNIQ == NULL || PAIR_NAME == NULL || ALGN_PAIR == NULL)
    Clean_Exit(1);

  if (P1->freq < FREQ)
    { fprintf(stderr,"%s: Genome index %s/%s.gix cutoff of %d < requested cutoff\n",
                     Prog_Name,PATH1,ROOT1,P1->freq);
      Clean_Exit(1);
    }
  if (P2->freq < FREQ)
    { fprintf(stderr,"%s: Genome index %s/%s.gix cutoff %d < requested cutoff\n",
                     Prog_Name,PATH2,ROOT2,P2->freq);
      Clean_Exit(1);
    }
  if (T1->kmer != T2->kmer)
    { fprintf(stderr,"%s: Indices not made with the same k-mer size (%d vs %d)\n",
                     Prog_Name,T1->kmer,T2->kmer);
      Clean_Exit(1);
    }

    //  Make sure you can open (NTHREADS + 3) * 2 * NTHREADS + tid files at one time.
    //    tid is typically 3 unless using valgrind or other instrumentation.

    { struct rlimit rlp;
      int           tid;
      uint64        nfiles;

      tid = open(".xxx",O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
      close(tid);
      unlink(".xxx");

      nfiles = (NTHREADS+3)*2*NTHREADS + tid + 100 ; // RD 250221 add 100 to allow for a few wrapping processes
      getrlimit(RLIMIT_NOFILE,&rlp);
      if (nfiles > rlp.rlim_max)
        { fprintf(stderr,"\n%s: Cannot open %lld files simultaneously\n",Prog_Name,nfiles);
          Clean_Exit(1);
        }
      rlp.rlim_cur = nfiles;
      setrlimit(RLIMIT_NOFILE,&rlp);
    }

  IBYTE = P1->pbyte;
  ICONT = P1->cbyte;
  IPOST = IBYTE-ICONT;
  ISIGN = IBYTE-1;

  JBYTE = P2->pbyte;
  JCONT = P2->cbyte;
  JPOST = JBYTE-JCONT;
  JSIGN = JBYTE-1;

  KBYTE = T2->pbyte;
  CBYTE = T2->hbyte;
  LBYTE = CBYTE+1;

  ESHIFT = 8*IPOST;

  { int64 cum;      // DBYTE accommodates the sum of the largest contig positions in each GDB !
    int   r, len;

    AMXPOS = 0;
    for (r = 0; r < gdb1->ncontig; r++)
      { len = gdb1->contigs[r].clen;
        if (len > AMXPOS)
          AMXPOS = len;
      }

    if (SELF)
      BMXPOS = AMXPOS;
    else
      { BMXPOS = 0;
        for (r = 0; r < gdb2->ncontig; r++)
          { len = gdb2->contigs[r].clen;
            if (len > BMXPOS)
              BMXPOS = len;
          }
      }

    MAXDAG = AMXPOS+BMXPOS;
    DBYTE = 0;
    cum   = 1;
    while (cum < MAXDAG)
      { cum   *= 256;
        DBYTE += 1;
      }
  }

  if (VERBOSE)
    { fprintf(stderr,"\n  Using %d threads\n",NTHREADS);
      fflush(stderr);
    }

  { int64 npost, cum, t;   //  Compute GDB split into NTHREADS parts
    int   p, r, x;

    NCONTS = gdb1->ncontig;

    IDBsplit = Malloc((NTHREADS+1)*sizeof(int),"Allocating GDB1 partitions");
    Select   = Malloc(NCONTS*sizeof(int),"Allocating GDB1 partition");
    if (IDBsplit == NULL || Select == NULL)
      Clean_Exit(1);

    npost = gdb1->seqtot;
    IDBsplit[0] = 0;
    Select[0] = 0;
    p = 0;
    r = NTHREADS;
    t = npost/NTHREADS;
    cum = gdb1->contigs[Perm1[0]].clen;
    for (x = 1; x < NCONTS; x++)
      { if (cum >= t && x >= r)
          { p += 1;
            IDBsplit[p] = x;
            t = (npost*(p+1))/NTHREADS;
            r += NTHREADS;
          }
        Select[x] = p;
        cum += gdb1->contigs[Perm1[x]].clen;
      }
    NPARTS = p+1;
    IDBsplit[NPARTS] = NCONTS;

#ifdef DEBUG_SPLIT
    for (x = 0; x < NPARTS; x++)
      { printf(" %2d: %4d - %4d\n",x,IDBsplit[x],IDBsplit[x+1]);
        for (r = IDBsplit[x]; r < IDBsplit[x+1]; r++)
          if (Select[r] != x)
            printf("  Not OK: %d->%d\n",r,Select[x]);
      }
#endif
  }

  { int    i, j, k, x;   // Setup temporary pair file IO buffers
    uint8 *buffer;
    int64 *bucks;
    char  *name;
    int   *nfile, *cfile;

    N_Units = Malloc(NPARTS*NTHREADS*sizeof(IOBuffer),"IO buffers");
    C_Units = Malloc(NPARTS*NTHREADS*sizeof(IOBuffer),"IO buffers");
    buffer  = Malloc(2*NPARTS*NTHREADS*1000000,"IO buffers");
    bucks   = Malloc(2*NTHREADS*NCONTS*sizeof(int64),"IO buffers");
    if (N_Units == NULL || C_Units == NULL || buffer == NULL || bucks == NULL)
      Clean_Exit(1);

    k = 0;
    for (i = 0; i < NTHREADS; i++)
      for (j = 0; j < NPARTS; j++)
        { N_Units[k].bufr = buffer + (2*k) * 1000000; 
          C_Units[k].bufr = buffer + (2*k+1) * 1000000; 
          N_Units[k].buck = bucks + (2*i) * NCONTS; 
          C_Units[k].buck = bucks + (2*i+1) * NCONTS; 
          N_Units[k].inum = k;
          C_Units[k].inum = k;
          name = Catenate(SORT_PATH,"/",PAIR_NAME,Numbered_Suffix(".",k,".N"));
          N_Units[k].file = open(name,O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
          if (N_Units[k].file < 0)
            { fprintf(stderr,"%s: Cannot open %s for reading & writing\n",Prog_Name,name);
              Clean_Exit(1);
            }
          unlink(name);
          name = Catenate(SORT_PATH,"/",PAIR_NAME,Numbered_Suffix(".",k,".C"));
          C_Units[k].file = open(name,O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
          if (C_Units[j].file < 0)
            { fprintf(stderr,"%s: Cannot open %s for reading & writing\n",Prog_Name,name);
              Clean_Exit(1);
            }
          unlink(name);
          k += 1;
        }

#ifdef DEBUG_SPLIT
    for (i = 0; i < NTHREADS; i++)
      { for (j = 0; j < NPARTS; j++)
          printf("%d ",N_Units[i*NPARTS+j].file);
        printf("\n");
        for (j = 0; j < NPARTS; j++)
          printf("%ld ",N_Units[i*NPARTS+j].buck-bucks);
        printf("\n");
        for (j = 0; j < NPARTS; j++)
          printf("%ld ",N_Units[i*NPARTS+j].bufr-buffer);
        printf("\n");
      }
#endif

    if (SELF)
      self_adaptamer_merge(T1,P1);
    else
      adaptamer_merge(T1,T2,P1,P2);

    if (VERBOSE)
      TimeTo(stderr,0);

    //  Transpose N_unit & C_unit matrices

    nfile = (int *) buffer;
    cfile = nfile + NPARTS*NTHREADS;
    k = 0;
    for (j = 0; j < NPARTS; j++)
      for (i = 0; i < NTHREADS; i++)
        { N_Units[k].bufr = 
          C_Units[k].bufr = buffer + i * (2*NPARTS*1000000); 
          x = i*NPARTS+j;
          nfile[k] = N_Units[x].file;
          cfile[k] = C_Units[x].file;
          N_Units[k].buck = bucks + (2*i) * NCONTS; 
          C_Units[k].buck = bucks + (2*i+1) * NCONTS; 
          N_Units[k].inum = x;
          C_Units[k].inum = x;
          k += 1;
        }
    k = 0;
    for (j = 0; j < NPARTS; j++)
      for (i = 0; i < NTHREADS; i++)
        { N_Units[k].file = nfile[k];
          C_Units[k].file = cfile[k];
          k += 1;
        }

#ifdef DEBUG_SPLIT
    for (j = 0; j < NPARTS; j++)
      { for (i = 0; i < NTHREADS; i++)
          printf("%d ",N_Units[j*NTHREADS+i].file);
        printf("\n");
        for (i = 0; i < NTHREADS; i++)
          printf("%ld ",N_Units[j*NTHREADS+i].buck-bucks);
        printf("\n");
        for (i = 0; i < NTHREADS; i++)
          printf("%ld ",N_Units[j*NTHREADS+i].bufr-buffer);
        printf("\n");
      }
#endif
  
    pair_sort_search(gdb1,gdb2);

    if (VERBOSE)
      TimeTo(stderr,0);

    free(N_Units->buck);
    free(N_Units->bufr);
    free(C_Units);
    free(N_Units);

    if (OUT_TYPE != 2)
      { char *command;

        if (VERBOSE)
          { fprintf(stderr,"\n  Converting aln's to %s-format\n",OUT_TYPE==0?"PAF":"PSL");
            fflush(stderr);
          }

        command = Malloc(strlen(ONE_ROOT)+strlen(ONE_PATH)+100,"Allocating command buffer");
        if (command == NULL)
          { unlink(Catenate(ONE_PATH,"/",ONE_ROOT,".1aln"));
            Clean_Exit(1);
          }

        switch (OUT_TYPE)
        { case 0: // PAF
            sprintf(command,"ALNtoPAF%s -T%d %s/%s",
                            OUT_OPT==2?" -x":(OUT_OPT==1?" -m":""),NTHREADS,ONE_PATH,ONE_ROOT);
            break;
          case 1: // PSL
            sprintf(command,"ALNtoPSL -T%d %s/%s",NTHREADS,ONE_PATH,ONE_ROOT);
            break;
        }

        if (system(command) != 0)
          { switch (OUT_TYPE)
            { case 0:
                fprintf(stderr,"\n%s: Call to ALNtoPAF failed\n",Prog_Name);
                break;
              case 1:
                fprintf(stderr,"\n%s: Call to ALNtoPSL failed\n",Prog_Name);
                break;
            }
            unlink(Catenate(ONE_PATH,"/",ONE_ROOT,".1aln"));
            Clean_Exit(1);
          }

        free(command);
        unlink(Catenate(ONE_PATH,"/",ONE_ROOT,".1aln"));

        if (VERBOSE)
          TimeTo(stderr,0);
      }
  }

  if (VERBOSE)
    TimeTo(stderr,1);

  free(Select);
  free(IDBsplit);

  free(ALGN_PAIR);
  free(PAIR_NAME);
  free(ALGN_UNIQ);
  free(ONE_ROOT);

  if ( ! SELF)
    Close_GDB(gdb2);
  Close_GDB(gdb1);

  if ( ! SELF)
    Free_Post_List(P2);
  Free_Post_List(P1);

  free(SORT_PATH);

  if (TYPE2 <= IS_GDB)
    free(SPATH2);
  if (TYPE1 <= IS_GDB)
    free(SPATH1);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  Clean_Exit(0);
}
