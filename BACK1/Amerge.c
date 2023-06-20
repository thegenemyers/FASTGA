/*******************************************************************************************
 *
 *  Adaptamer merge phase of a WGA.  Given indices of two genomes (at a specified frequency
 *   cutoff), the adaptemer matches between the k-mers are found in a novel, cache-coherent
 *   merge of the sorted k-mer tables for each genome and seed position pairs are output for
 *   each adaptemer match.
 *
 *  Author:  Gene Myers
 *  Date  :  February 2023
 *
 *******************************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include "libfastk.h"
#include "DB.h"

#undef    DEBUG_MERGE
#undef    DEBUG_SORT
#undef    DEBUG_SEARCH
#undef    DEBUG_HIT

#define   CHAIN_BREAK  500
#define   CHAIN_MIN    100

static char *Usage = "[-v] -f<int> <source1>[.dam] <source2>[.dam]";

static int   FREQ;     //  Adaptemer frequence cutoff parameter
static int   VERBOSE;  //  Verbose output
static int   KMER;
static int   NTHREADS;
static char *PAIR_NAME;

static int IBYTE;   // # of bytes for a post in G1/P1 (excluding sign/slice)
static int JBYTE;   // # of bytes for a post in G2/P2 (excluding sign/slice)
static int KBYTE;   // # of bytes for a k-mer table entry (both T1 & T2)
static int CBYTE;   // byte of k-mer table entry containing post count
static int LBYTE;   // byte of k-mer table entry containing lcp
static int DBYTE;   // # of bytes for a pair diagonal

static int   *IDBsplit;     //  DB split: contig [DBsplit[p],DBsplit[p+1])
static int64 *IDBpost;      //  DB post: post of start of contig DBsplit[p] is DBpost[p]

static int   *JDBsplit;     //  DB split: contig [DBsplit[p],DBsplit[p+1])
static int64 *JDBpost;      //  DB post: post of start of contig DBsplit[p] is DBpost[p]

typedef struct
  { char  *name;
    uint8 *bufr;
    uint8 *btop;
    uint8 *bend;
    int64 *buck;
    int    file;
  } IOBuffer;

static IOBuffer *N_Units;  //  NTHREADS^2 IO units for + pair temporary files
static IOBuffer *C_Units;  //  NTHREADS^2 IO units for - pair temporary files

typedef struct
  { int   beg;
    int   end;
    int64 off;
  } Range;

extern void rmsd_sort(uint8 *array, int64 nelem, int rsize, int ksize,
                      int64 *part, int nthreads, Range *range);


/***********************************************************************************************
 *
 *   POSITION LIST ABSTRACTION:
 *       Routines to manipuate the position or "post" list associated with a k-mer table.
 *
 **********************************************************************************************/

typedef struct
  { int     pbyte;      //  # of bytes for each post (including sign & slice)
    int64   nels;       //  # of posts in the index
    int64   maxp;
    int     freq;
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
      lseek(copn,sizeof(int)+sizeof(int64),SEEK_SET);
    }
  P->cptr = cache;
  P->ctop = ctop;
  P->copn = copn;
}

static Post_List *Open_Post_List(char *name)
{ Post_List *P;
  int        pbyte;
  int64      nels, maxp, n;
  int        copn;

  int    f, p, flen;
  char  *dir, *root, *full;
  int    pb, nfile, nthreads, freq;

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
  read(f,&nfile,sizeof(int));
  read(f,&maxp,sizeof(int64));
  read(f,&freq,sizeof(int));
  nthreads = nfile;
  nfile    = nfile*nfile;

  P = Malloc(sizeof(Post_List),"Allocating post record");
  if (P == NULL)
    exit (1);
  P->name   = full;
  P->nlen   = strlen(full);
  P->maxp   = maxp;
  P->cache  = Malloc(POST_BLOCK*pbyte,"Allocating post list buffer\n");
  P->neps   = Malloc(nfile*sizeof(int64),"Allocating parts table of Post_List");
  if (P->cache == NULL || P->neps == NULL)
    exit (1);

  nels = 0;
  for (p = 1; p <= nfile; p++)
    { sprintf(P->name+P->nlen,"%d",p);
      copn = open(P->name,O_RDONLY);
      if (copn < 0)
        { fprintf(stderr,"%s: Table part %s is missing ?\n",Prog_Name,P->name);
          exit (1);
        }
      read(copn,&pb,sizeof(int));
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
  P->nels  = nels;
  P->nthr  = nfile;
  P->nsqrt = nthreads;
  P->freq  = freq;

  sprintf(P->name+P->nlen,"%d",1);
  copn = open(P->name,O_RDONLY);
  lseek(copn,sizeof(int)+sizeof(int64),SEEK_SET);

  P->copn  = copn;
  P->part  = 1;

  More_Post_List(P);
  P->cidx = 0;

  return (P);
}

static void Free_Post_List(Post_List *P)
{ free(P->cache);
  free(P->neps);
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
  memcpy((uint8 *) (&post),P->cptr,P->pbyte);
  return (post);
}

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
      sprintf(P->name+P->nlen,"%d",p);
      P->copn = open(P->name,O_RDONLY);
      P->part = p;
    }

  lseek(P->copn,sizeof(int) + sizeof(int64) + i*P->pbyte,SEEK_SET);

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
      sprintf(P->name+P->nlen,"%d",p);
      P->copn = open(P->name,O_RDONLY);
      P->part = p;
    }

  lseek(P->copn,sizeof(int) + sizeof(int64) + i*P->pbyte,SEEK_SET);

  More_Post_List(P);
}


/***********************************************************************************************
 *
 *   The internal data structure for a table (taken from libfastk.c) needs to be visible so
 *     that the "neps" array can be accessed in order to synchronize thread starts.
 *
 **********************************************************************************************/

typedef struct
  { int    kmer;       //  Kmer length
    int    minval;     //  The minimum count of a k-mer in the stream
    int64  nels;       //  # of elements in entire table
                   //  Current position (visible part)
    int64  cidx;       //  current element index
    uint8 *csuf;       //  current element suffix
    int    cpre;       //  current element prefix
                   //  Other useful parameters
    int    ibyte;      //  # of bytes in prefix
    int    kbyte;      //  Kmer encoding in bytes
    int    tbyte;      //  Kmer+count entry in bytes
    int    hbyte;      //  Kmer suffix in bytes (= kbyte - ibyte)
    int    pbyte;      //  Kmer,count suffix in bytes (= tbyte - ibyte)
                   //  Hidden parts
    int    ixlen;      //  length of prefix index (= 4^(4*ibyte))
    int    shift;      //  shift for inverse mapping
    uint8 *table;      //  The (huge) table in memory
    int64 *index;      //  Prefix compression index
    int   *inver;      //  inverse prefix index
    int    copn;       //  File currently open
    int    part;       //  Thread # of file currently open
    int    nthr;       //  # of thread parts
    int    nlen;       //  length of path name
    char  *name;       //  Path name for table parts (only # missing)
    uint8 *ctop;       //  Ptr top of current table block in buffer
    int64 *neps;       //  Size of each thread part in elements
    int    clone;      //  Is this a clone?
  } _Kmer_Stream;


/***********************************************************************************************
 *
 *   ADAPTAMER MERGE THREAD:  
 *     For each k-mer in T1
 *       Find the longest prefix match to one or more k-mers in T2.
 *       If the totwl # of positions of the k-mers in T2, then output the pairs of positions
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
    int          tid;
    uint8       *cache;
    IOBuffer    *nunit;
    IOBuffer    *cunit;
    int64        nhits;
    int64        g1len;
  } SP;

static void *merge_thread(void *args)
{ SP *parm = (SP *) args;
  int tid          = parm->tid;
  uint8 *cache     = parm->cache;
  IOBuffer  *nunit = parm->nunit;
  IOBuffer  *cunit = parm->cunit;
  Kmer_Stream *T1  = parm->T1;
  Post_List   *P1  = parm->P1;
  Kmer_Stream *T2  = parm->T2;
  Post_List   *P2  = parm->P2;

  int     spart;
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
  int64   nhits, g1len;

#ifdef DEBUG_MERGE
  int64   Tdp;
  char   *tbuffer;
#endif

  { int j;

    for (j = 0; j < NTHREADS; j++)
      { nunit[j].bend = nunit[j].bufr + (1000000-(IBYTE+JBYTE+2));
        nunit[j].btop = nunit[j].bufr;
        cunit[j].bend = cunit[j].bufr + (1000000-(IBYTE+JBYTE+2));
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

  spart = P1->nsqrt * tid - 1;
  First_Post_Entry(P1);
  First_Post_Entry(P2);
  First_Kmer_Entry(T1);
  First_Kmer_Entry(T2);

  if (tid != 0)
    { GoTo_Post_Index(P1,P1->neps[spart]);
      GoTo_Post_Index(P2,P2->neps[spart]);
      GoTo_Kmer_Index(T1,((_Kmer_Stream *) T1)->neps[spart]);
      GoTo_Kmer_Index(T2,((_Kmer_Stream *) T2)->neps[spart]);
    }

  qcnt = -1;
  tend = ((_Kmer_Stream *) T1)->neps[spart+P1->nsqrt];
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
            { pcnt = ((T1->cidx - tbeg) * 100) / (tend-tbeg); 
              if (pcnt > qcnt)
                { printf("\r    Completed %3d%%",pcnt);
                  fflush(stdout);
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
          printf("... to %lld rcur = rend = %lld, eorunn = 0, plen = 12\n",
                 T2->cidx,Tdp+(rcur-cache)/KBYTE);
          fflush(stdout);
#endif
        }

#define ADVANCE(l) 	 					\
{ if (l >= vhgh)						\
    { int n;							\
								\
      for (n = l[CBYTE]; n > 0; n--)				\
        { pdx = (pdx+1) & POST_BUF_MASK;			\
          post[pdx] = Current_Post(P2);				\
          Next_Post_Entry(P2);					\
        }							\
      vhgh = l+KBYTE;						\
    }								\
  cdx = (cdx + l[CBYTE]) & POST_BUF_MASK;			\
  l += KBYTE;							\
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
printf(" P2 = %10lld\n",P2->cidx);
      fflush(stdout);
#endif

      { int       freq, lcs, udx;
        int       apane, asign;
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

        nhits += suf1[CBYTE] * freq;
        g1len += suf1[CBYTE];

        for (n = suf1[CBYTE]; n > 0; n--)
          { apost = Current_Post(P1);
            apane = aptr[IBYTE];
            asign = apane & 0x80;
            apane = apane & 0x7f;
            jptr  = (uint8 *) (post+b);
            for (k = 0; k < freq; k++)
              { if (asign == (jptr[JBYTE] & 0x80))
                  ou = nunit + apane;
                else
                  ou = cunit + apane;
                btop  = ou->btop;
                *btop++ = plen;
                memcpy(btop,aptr,IBYTE);
                btop += IBYTE;
                memcpy(btop,jptr,JBYTE);
                btop += JBYTE;
                ou->buck[*btop++ = (jptr[JBYTE] & 0x7f)] += 1;

#ifdef DEBUG_MERGE
                if (n == suf1[CBYTE])
                  { int64 ip;

                    printf("      %ld: %c (%d)",((int64 *) jptr)-post,(jptr[JBYTE] & 0x80)?'-':'+',
                                                jptr[JBYTE]&0x7f);
                    ip = 0;
                    memcpy((uint8 *) (&ip),jptr,JBYTE);
                    printf(" %lld\n",ip);
                    fflush(stdout);
                  }
#endif

                if (btop >= ou->bend)
                  { write(ou->file,ou->bufr,btop-ou->bufr);
                    ou->btop = ou->bufr;
                  }
                else
                  ou->btop = btop;

                jptr += sizeof(int64);
              }
#ifdef DEBUG_MERGE
            if (n == suf1[CBYTE])
              printf("   vs\n");
            printf("      %lld: %c (%d)",P1->cidx,asign?'-':'+',apane);
            aptr[IBYTE] = 0;
            printf(" %lld\n",apost);
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

    for (j = 0; j < NTHREADS; j++)
      { if (nunit[j].btop > nunit[j].bufr)
          write(nunit[j].file,nunit[j].bufr,nunit[j].btop-nunit[j].bufr);
        close(nunit[j].file);
        if (cunit[j].btop > cunit[j].bufr)
          write(cunit[j].file,cunit[j].bufr,cunit[j].btop-cunit[j].bufr);
        close(cunit[j].file);
      }
  }

  parm->nhits = nhits;
  parm->g1len = g1len;
  return (NULL);
}
  

static void adaptamer_merge(char *g1,        char *g2,
                            Kmer_Stream *T1, Kmer_Stream *T2,
                            Post_List *P1,   Post_List *P2)
{ SP         parm[NTHREADS];
#ifndef DEBUG_MERGE
  pthread_t  threads[NTHREADS];
#endif
  uint8     *cache;
  int64      nhits, g1len;
  int        i, k;

  if (VERBOSE)
    { fprintf(stdout,"  Starting adaptamer merge\n");
      fflush(stdout);
    }

  parm[0].T1 = T1;
  parm[0].T2 = T2;
  parm[0].P1 = P1;
  parm[0].P2 = P2;
  for (i = 1; i < NTHREADS; i++)
    { parm[i].T1 = Clone_Kmer_Stream(T1);
      parm[i].T2 = Clone_Kmer_Stream(T2);
      parm[i].P1 = Open_Post_List(g1);
      parm[i].P2 = Open_Post_List(g2);
    }

  cache = Malloc(NTHREADS*(P2->maxp+1)*KBYTE,"Allocating cache");
  if (cache == NULL)
    exit (1);

  for (i = 0; i < NTHREADS; i++)
    { IOBuffer *nu, *cu;

      parm[i].tid   = i;
      parm[i].cache = cache + i * (P2->maxp+1) * KBYTE;
      parm[i].nunit = nu = N_Units + i * NTHREADS;
      parm[i].cunit = cu = C_Units + i * NTHREADS;
      for (k = 0; k < NTHREADS; k++)
        { nu[k].file = open(nu[k].name,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU);
          if (nu[k].file < 0)
            { fprintf(stderr,"%s: Cannot open %s for reading\n",Prog_Name,nu[k].name);
              exit (1);
            }
          cu[k].file = open(cu[k].name,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU);
          if (cu[k].file < 0)
            { fprintf(stderr,"%s: Cannot open %s for reading\n",Prog_Name,cu[k].name);
              exit (1);
            }
          bzero(nu[k].buck,sizeof(int64)*NTHREADS);
          bzero(cu[k].buck,sizeof(int64)*NTHREADS);
        }
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
    { printf("\r    Completed 100%%\n");
      fflush(stdout);
    }
   
  free(cache);

  for (i = NTHREADS-1; i >= 0; i--)
    { Free_Kmer_Stream(parm[i].T1);
      Free_Kmer_Stream(parm[i].T2);
      Free_Post_List(parm[i].P1);
      Free_Post_List(parm[i].P2);
    }

  nhits = g1len = 0;
  for (i = 0; i < NTHREADS; i++)
    { nhits += parm[i].nhits;
      g1len += parm[i].g1len;
    }

  if (VERBOSE)
    printf("\n  Total seeds = %lld,  per G1 position = %.1f\n\n",nhits,(1.*nhits)/g1len);
}



/***********************************************************************************************
 *
 *   PAIR RE-IMPORT AND SORT
 *
 **********************************************************************************************/

typedef struct
  { int       in;
    int       swide;
    int64     ioff;
    int64    *buck;
    uint8    *buffer;
    uint8    *sarr;
    Range    *range;
  } RP;

static void *reimport_N_thread(void *args)
{ RP *parm = (RP *) args;
  int    swide  = parm->swide;
  int    in     = parm->in;
  int64  ioff   = parm->ioff;
  uint8 *sarr   = parm->sarr;
  uint8 *bufr   = parm->buffer;
  int64 *buck   = parm->buck;

  int64  ipost, jpost, pdiag;
  uint8 *_ipost = (uint8 *) (&ipost);
  uint8 *_jpost = (uint8 *) (&jpost);
  uint8 *_pdiag = (uint8 *) (&pdiag);

  uint8 *x;
  int    iolen, iunit, pane, lcp;
  int64  drem;
  uint8 *bend, *btop, *b;

  iolen = 2*NTHREADS*1000000;
  iunit = IBYTE + JBYTE + 2;
  bend = bufr + read(in,bufr,iolen);
  if (bend-bufr < iolen)
    btop = bend;
  else
    btop = bend-iunit;
  b = bufr;

  ipost = jpost = 0;

  while (1)
    { lcp = *b++;
      memcpy(_ipost,b,IBYTE);
      b += IBYTE;
      memcpy(_jpost,b,JBYTE);
      b += JBYTE;
      pane = *b++;

      x = sarr + swide * buck[pane]++;
      *x++ = lcp;
      drem  = (jpost-ipost)+ioff;
      pdiag = (drem >> 6);
      *x++ = drem-(pdiag<<6);
      memcpy(x,_ipost,IBYTE);
      x += IBYTE;
      memcpy(x,_pdiag,DBYTE);
      x += DBYTE;
      *x++ = 0;

      if (b >= btop)
        { int ex = bend-b;
          memcpy(bufr,b,ex);
          bend = bufr+ex;
          bend += read(in,bend,iolen-ex);
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

static void *reimport_C_thread(void *args)
{ RP *parm = (RP *) args;
  int    swide  = parm->swide;
  int    in     = parm->in;
  uint8 *sarr   = parm->sarr;
  uint8 *bufr   = parm->buffer;
  int64 *buck   = parm->buck;

  int64  ipost, jpost, panti;
  uint8 *_ipost = (uint8 *) (&ipost);
  uint8 *_jpost = (uint8 *) (&jpost);
  uint8 *_panti = (uint8 *) (&panti);

  uint8 *x;
  int    iolen, iunit, pane, lcp;
  int64  arem;
  uint8 *bend, *btop, *b;

  iolen = 2*NTHREADS*1000000;
  iunit = IBYTE + JBYTE + 2;
  bend = bufr + read(in,bufr,iolen);
  if (bend-bufr < iolen)
    btop = bend;
  else
    btop = bend-iunit;
  b = bufr;

  ipost = jpost = 0;

  while (1)
    { lcp = *b++;
      memcpy(_ipost,b,IBYTE);
      b += IBYTE;
      memcpy(_jpost,b,JBYTE);
      b += JBYTE;
      pane = *b++;

      x = sarr + swide * buck[pane]++;
      *x++ = lcp;
      arem  = jpost+ipost;
      panti = (arem >> 6);
      *x++ = arem-(panti<<6);
      memcpy(x,_ipost,IBYTE);
      x += IBYTE;
      memcpy(x,_panti,DBYTE);
      x += DBYTE;
      *x++ = 0;

      if (b >= btop)
        { int ex = bend-b;
          memcpy(bufr,b,ex);
          bend = bufr+ex;
          bend += read(in,bend,iolen-ex);
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

void print_seeds(uint8 *sarray, int64 nels, int swide)
{ uint8 *e;
  int64  i;

  int64  ipost, diag;
  uint8 *_ipost = (uint8 *) (&ipost);
  uint8 *_diag = (uint8 *) (&diag);

  ipost = 0;
  diag  = 0;

  e = sarray;
  for (i = 0; i < nels; i++)
    { memcpy(_diag,e+(IBYTE+2),DBYTE);
      memcpy(_ipost,e+2,IBYTE);
      printf("  %10lld:  %10lld  %10lld  (%2d)  %2d\n",i,diag,ipost,e[1],e[0]);
      e += swide;
    }
}

typedef struct
  { int       pane;
    int       swide;
    int64     ioff;
    uint8    *sbeg;
    uint8    *send;
  } TP;

static void *search_N_seeds(void *args)
{ TP *parm = (TP *) args;
  int    pane   = parm->pane;
  int    swide  = parm->swide;
  uint8 *sbeg   = parm->sbeg;
  uint8 *send   = parm->send;

  uint8 **list;
  int     lmax;

  uint8 *b, *m, *e;

int64 nhit;

  int    new, aux;
  int64  ndiag, cdiag;
  uint8 *_ndiag = (uint8 *) (&ndiag);

  int64  ipost, apost, mpost;
  uint8 *_ipost = (uint8 *) (&ipost);
  uint8 *_apost = (uint8 *) (&apost);
  uint8 *_mpost = (uint8 *) (&mpost);

  ndiag = 0;
  ipost = 0;
  apost = 0;
  mpost = 0;

  lmax = 1000;
  list = Malloc(lmax*sizeof(uint8 *),"Allocating merge index");
  if (list == NULL)
    exit (1);

nhit = 0;

  (void) pane;

  b = e = sbeg;
  memcpy(_ndiag,e+(IBYTE+2),DBYTE);
  cdiag = ndiag;
  while (ndiag == cdiag && e < send)
    { e += swide;
      memcpy(_ndiag,e+(IBYTE+2),DBYTE);
    }

  while (1)
    { m = e;
      aux = 0;
      while (ndiag == cdiag+1 && e < send)
        { e += swide;
          memcpy(_ndiag,e+(IBYTE+2),DBYTE);
          aux = 1;
        }

      if (new || aux)
        { int    go, lcp, mix;
          int64  lps, cov, npost;
          uint8 *s, *t;
          int    len, lst;
       
#ifdef DEBUG_SEARCH
          printf("Diag %lld",cdiag);
          if (aux)
            printf("+1");
          printf("\n");
#endif

          if (e-b > lmax)
            { lmax = (e-b)*1.2+100;
              list = Realloc(list,lmax*sizeof(uint8 *),"Expanding merge index");
              if (list == NULL)
                exit (1);
            }

          s = b;
          memcpy(_ipost,s+2,IBYTE);
          t = m; 
          if (aux)
            memcpy(_apost,t+2,IBYTE);
          else
            apost = 0x7fffffffffffffffll;
          lps = 0x8000000000000000ll;
 
          mix = 0;
          go  = 1;
          for (len = 0; go; len++)
            { if (ipost < apost)
                { lcp = s[0];
                  npost = ipost;
                  list[len] = s;
                  s += swide;
                  if (s >= m)
                    ipost = 0x7fffffffffffffffll;
                  else
                    memcpy(_ipost,s+2,IBYTE);
                  mix |= 0x1;
                }
              else
                { lcp = t[0];
                  npost = apost;
                  list[len] = t;
                  t += swide;
                  if (t >= e)
                    { if (t > e)
                        go = 0;
                      else
                        apost = 0x7fffffffffffffffll;
                    }
                  else
                    memcpy(_apost,t+2,IBYTE);
                  mix |= 0x2;
                }

              if (npost >= lps + CHAIN_BREAK)
                { if (lps >= 0)
                    { if (cov > CHAIN_MIN && (mix != 1 || new))
                        { int64  dg, lg;
                          int64  fp, lp;
                          uint8 *n;
                          int    k;

nhit += 1;

#ifdef DEBUG_SEARCH
                          printf("                  Process\n");
#endif
#ifdef DEBUG_HIT
                          printf("Hit %lld",cdiag);
                          if (aux)
                            printf("+1");
                          printf(" %lld %d\n",cov,lst);
#endif

                          n = list[lst];
                          memcpy(_mpost,n+2,IBYTE);
                          fp = mpost;
                          lp = mpost + n[0];
                          if (n < m)
                            lg = (cdiag << 6) + n[1];
                          else
                            lg = (cdiag << 6) + n[1] + 64;
                          for (k = lst+1; k <= len; k++)
                            { if (k == len)
                                dg = -1;
                              else
                                { n = list[k];
                                  memcpy(_mpost,n+2,IBYTE);
                                  if (n < m)
                                    dg = (cdiag << 6) + n[1];
                                  else
                                    dg = (cdiag << 6) + n[1] + 64;
                                }
                              if (dg == lg && mpost <= lp)
                                { if (mpost + n[0] > lp)
                                    lp = mpost + n[0];
                                }
                              else
                                {
#ifdef DEBUG_HIT
                                  if (dg < 0)
                                    printf("        TERM:  %10lld %2lld [%lld]\n",fp,lp-fp,lg);
                                  else
                                    printf("  %10ld:  %10lld %2lld [%lld]\n",
                                           (n-sbeg)/swide,fp,lp-fp,lg);
#endif
                                  fp = mpost;
                                  lg = dg;
                                  lp = mpost + n[0];
                                }
                            }
                          fflush(stdout);
                        }
#ifdef DEBUG_SEARCH
                      else
                        printf("                  Break\n");
#endif
                    }
                  cov = lcp;
                  lps = npost + lcp;
                  lst = len;
                  mix = 0;
                }
              else
                { int64 cps;

                  cps = npost + lcp;
                  if (cps > lps)
                    { if (npost >= lps)
                        cov += lcp;
                      else
                        cov += cps-lps;
                      lps = cps;
                    }
                }
#ifdef DEBUG_SEARCH
              if (go)
                { uint8 *n = list[len];
                  printf("   %c %10ld: %10lld %2d %4lld (%d)\n",
                         n<m?'.':'+',(n-sbeg)/swide,npost,n[0],cov,n[1]);
                }
#endif
            }
          ipost = apost = 0;
        }

      if (e >= send) break;

      if (aux)
        { b = m;
          cdiag += 1;
          new = 0;
        }
      else
        { b = e;
          cdiag = ndiag;
          while (ndiag == cdiag && e < send)
            { e += swide;
              memcpy(_ndiag,e+(IBYTE+2),DBYTE);
            }
          new = 1;
        }
    }

  free(list);

printf("Nhits = %lld\n",nhit);

  return (NULL);
}

static void *search_C_seeds(void *args)
{ TP *parm = (TP *) args;
  int    pane   = parm->pane;
  int    swide  = parm->swide;
  uint8 *sbeg   = parm->sbeg;
  uint8 *send   = parm->send;

  (void) pane;
  (void) swide;
  (void) sbeg;
  (void) send;

  return (NULL);
}


static void pair_sort_search()
{ uint8 *sarray;
  int    swide;
  RP     rarm[NTHREADS];
  TP     tarm[NTHREADS];
#ifndef DEBUG_SORT
  pthread_t threads[NTHREADS];
#endif
  int64     panel[NTHREADS];
  Range     range[NTHREADS];
  int64     ioff;
  IOBuffer *unit[2], *nu;
  int       i, p, j, u;

  (void) reimport_C_thread;

  if (VERBOSE)
    { fprintf(stdout,"  Starting seed sort and alignment search, %d parts\n",2*NTHREADS);
      fflush(stdout);
    }

  unit[0] = N_Units;
  unit[1] = C_Units;

  { int64 x, z, nelmax;

    nelmax = 0;
    for (u = 0; u < 2; u++)
      for (i = 0; i < NTHREADS; i++)
        { nu = unit[u] + i*NTHREADS;
          x = 0;
          for (j = 0; j < NTHREADS; j++) 
            for (p = 0; p < NTHREADS; p++)
              { z = nu[p].buck[j];
                nu[p].buck[j] = x;
                x += z;
              }
          if (x > nelmax)
            nelmax = x;
        }

    swide  = IBYTE + DBYTE + 3;
    sarray = Malloc((nelmax+1)*swide,"Sort Array");
    if (sarray == NULL)
      exit (1);
  }

  for (p = 0; p < NTHREADS; p++)
    { rarm[p].swide  = swide;
      rarm[p].sarr   = sarray;
      rarm[p].buffer = N_Units[p].bufr;   //  NB: Units have been transposed
      rarm[p].range  = range+p;

      tarm[p].pane   = p;
      tarm[p].swide  = swide;
    }

  for (u = 0; u < 2; u++)
   for (i = 0; i < NTHREADS; i++)
    { nu = unit[u] + i*NTHREADS;

      if (VERBOSE)
        { fprintf(stdout,"\r    Loading seeds for part %d  ",u*NTHREADS+i+1);
          fflush(stdout);
        }

      ioff = IDBpost[i+1] - IDBpost[i];
      for (p = 0; p < NTHREADS; p++)
        { rarm[p].in = open(nu[p].name,O_RDONLY);
          if (rarm[p].in < 0)
            { fprintf(stderr,"%s: Cannot open %s for reading\n",Prog_Name,nu[p].name);
              exit (1);
            }
          rarm[p].buck = nu[p].buck;
          rarm[p].ioff = ioff;
        }

#ifdef DEBUG_SORT
      if (u == 0)
        for (p = 0; p < NTHREADS; p++)
          reimport_N_thread(rarm+p);
      else
        for (p = 0; p < NTHREADS; p++)
          reimport_C_thread(rarm+p);
#else
      if (u == 0)
        { for (p = 1; p < NTHREADS; p++)
            pthread_create(threads+p,NULL,reimport_N_thread,rarm+p);
          reimport_N_thread(rarm);
        }
      else
        { for (p = 1; p < NTHREADS; p++)
            pthread_create(threads+p,NULL,reimport_C_thread,rarm+p);
          reimport_C_thread(rarm);
        }
      for (p = 1; p < NTHREADS; p++)
        pthread_join(threads[p],NULL);
#endif

      for (p = 0; p < NTHREADS; p++)
        unlink(nu[p].name);

#ifdef DEBUG_SORT
      for (p = 0; p < NTHREADS; p++)
        printf("  %s",nu[p].name);
      printf("\n");
      for (j = 0; j < NTHREADS; j++)
        { for (p = 0; p < NTHREADS; p++)
            printf(" %10lld",nu[p].buck[j]);
          printf("\n");
        }
#endif

      { int64 prev, next, nels;

        bzero(panel,sizeof(int64)*NTHREADS);
        prev = 0;
        for (j = 0; j < NTHREADS; j++)
          { next = nu[NTHREADS-1].buck[j];
            panel[j] = (next - prev)*swide;
            prev = next;
          }
        nels = next;

#ifdef DEBUG_SORT
        for (p = 0; p < NTHREADS; p++)
          printf(" %2d: %10lld %10lld\n",p,panel[p],panel[p]/swide);
#endif

        if (VERBOSE)
          { fprintf(stdout,"\r    Sorting seeds for part %d  ",u*NTHREADS+i+1);
            fflush(stdout);
          }

        rmsd_sort(sarray,nels,swide,swide-2,panel,NTHREADS,range);

        for (p = 0; p < NTHREADS; p++)
          { tarm[p].ioff = ioff;
            tarm[p].sbeg = sarray + range[p].off;
          }
        for (p = 1; p < NTHREADS; p++)
          tarm[p-1].send = tarm[p].sbeg;
        tarm[NTHREADS-1].send = sarray + swide *nels;
      }

    if (VERBOSE)
      { fprintf(stdout,"\r    Searching seeds for part %d",u*NTHREADS+i+1);
        fflush(stdout);
      }

#if defined(DEBUG_SEARCH) || defined(DEBUG_HIT)
      if (u == 0)
        for (p = 0; p < NTHREADS; p++)
          search_N_seeds(tarm+p);
      else
        for (p = 0; p < NTHREADS; p++)
          search_C_seeds(tarm+p);
#else
      if (u == 0)
        { for (p = 1; p < NTHREADS; p++)
            pthread_create(threads+p,NULL,search_N_seeds,tarm+p);
          search_N_seeds(tarm);
        }
      else
        { for (p = 1; p < NTHREADS; p++)
            pthread_create(threads+p,NULL,search_C_seeds,tarm+p);
          search_C_seeds(tarm);
        }
      for (p = 1; p < NTHREADS; p++)
        pthread_join(threads[p],NULL);
#endif

#ifdef DEBUG_SORT
        print_seeds(sarray,nels,swide);
#endif
    }

  if (VERBOSE)
    fprintf(stdout,"\r    Done                        \n");
}


int main(int argc, char *argv[])
{ Kmer_Stream *T1, *T2;
  Post_List   *P1, *P2;
  DAZZ_DB _DB1, *DB1 = &_DB1;
  DAZZ_DB _DB2, *DB2 = &_DB2;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("Amerge");

    FREQ = -1;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'f':
            ARG_NON_NEGATIVE(FREQ,"maximum seed frequency");
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc != 3 || FREQ < 0)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  T1 = Open_Kmer_Stream(argv[1]);
  T2 = Open_Kmer_Stream(argv[2]);
  
  P1 = Open_Post_List(argv[1]);
  P2 = Open_Post_List(argv[2]);

  if (Open_DB(argv[1],DB1) < 0)
    { fprintf(stderr,"%s: Cannot open Dazzler DB %s\n",Prog_Name,argv[1]);
      exit (1);
    }
  Trim_DB(DB1);

  if (Open_DB(argv[1],DB2) < 0)
    { fprintf(stderr,"%s: Cannot open Dazzler DB %s\n",Prog_Name,argv[1]);
      exit (1);
    }
  Trim_DB(DB2);

  KMER      = T1->kmer;
  NTHREADS  = P1->nsqrt;
  PAIR_NAME = Strdup(Numbered_Suffix("._pair.",getpid(),"."),"Allocating temp name");

  if (P1->freq < FREQ)
    { fprintf(stderr,"%s: Genome index for %s cutoff %d < requested cutoff\n",
                     Prog_Name,argv[1],P1->freq);
      exit (1);
    }
  if (P2->freq < FREQ)
    { fprintf(stderr,"%s: Genome index for %s cutoff %d < requested cutoff\n",
                     Prog_Name,argv[2],P2->freq);
      exit (1);
    }
  if (T1->kmer != T2->kmer)
    { fprintf(stderr,"%s: Indices not made with the same k-mer size (%d vs %d)\n",
                     Prog_Name,T1->kmer,T2->kmer);
      exit (1);
    }

  IBYTE = P1->pbyte-1;
  JBYTE = P2->pbyte-1;
  KBYTE = T2->pbyte;
  CBYTE = T2->hbyte;
  LBYTE = CBYTE+1;
  if (IBYTE > JBYTE)
    DBYTE = IBYTE;
  else
    DBYTE = JBYTE;

  { int    k;   // Setup temporary pair file IO buffers
    uint8 *buffer;
    int64 *bucks;

    N_Units = Malloc(NTHREADS*NTHREADS*sizeof(IOBuffer),"IO buffers");
    C_Units = Malloc(NTHREADS*NTHREADS*sizeof(IOBuffer),"IO buffers");
    buffer  = Malloc(2*NTHREADS*NTHREADS*1000000,"IO buffers");
    bucks   = Malloc(2*NTHREADS*NTHREADS*NTHREADS*sizeof(int64),"IO buffers");
    if (N_Units == NULL || C_Units == NULL || buffer == NULL || bucks == NULL)
      exit (1);

    for (k = 0; k < NTHREADS*NTHREADS; k++)
      { N_Units[k].bufr = buffer + (2 * k) * 1000000; 
        C_Units[k].bufr = buffer + (2 * k + 1) * 1000000; 
        N_Units[k].name = Strdup(Numbered_Suffix(PAIR_NAME,k,".N"),"Temp file name");
        C_Units[k].name = Strdup(Numbered_Suffix(PAIR_NAME,k,".C"),"Temp file name");
        N_Units[k].buck = bucks + (2 * k) * NTHREADS; 
        C_Units[k].buck = bucks + (2 * k + 1) * NTHREADS; 
      }
  }

  { int64 npost, cum, t;   //  Compute DB split into NTHREADS parts
    int   p, r;

    IDBsplit = Malloc(2*(NTHREADS+1)*sizeof(int),"Allocating DB partitions");
    IDBpost  = Malloc(2*(NTHREADS+1)*sizeof(int64),"Allocating DB partitions");
    if (IDBsplit == NULL || IDBpost == NULL)
      exit (1);
    JDBsplit = IDBsplit + (NTHREADS+1);
    JDBpost  = IDBpost + (NTHREADS+1);

    npost = DB1->totlen;
    cum   = 0;

    IDBsplit[0] = 0;
    IDBpost [0] = 0;
    p = 1;
    t = (npost*p)/NTHREADS;
    for (r = 0; r < DB1->treads; r++)
      { cum += DB1->reads[r].rlen;
        if (cum >= t)
          { IDBsplit[p] = r+1;
            IDBpost [p] = cum;
            p += 1;
            t = (npost*p)/NTHREADS;
          }
      }
    IDBsplit[NTHREADS] = DB1->treads;
    IDBpost [NTHREADS] = npost;

    npost = DB2->totlen;
    cum   = 0;

    JDBsplit[0] = 0;
    JDBpost [0] = 0;
    p = 1;
    t = (npost*p)/NTHREADS;
    for (r = 0; r < DB2->treads; r++)
      { cum += DB2->reads[r].rlen;
        if (cum >= t)
          { JDBsplit[p] = r+1;
            JDBpost [p] = cum;
            p += 1;
            t = (npost*p)/NTHREADS;
          }
      }
    JDBsplit[NTHREADS] = DB2->treads;
    JDBpost [NTHREADS] = npost;
  }

  adaptamer_merge(argv[1],argv[2],T1,T2,P1,P2);

  { IOBuffer u;      //  Transpose unit matrix
    int      k, j;

    for (k = 0; k < NTHREADS; k++)
      for (j = k+1; j < NTHREADS; j++)
        { u = N_Units[j*NTHREADS+k];
          N_Units[j*NTHREADS+k] = N_Units[k*NTHREADS+j];
          N_Units[k*NTHREADS+j] = u;
          u = C_Units[j*NTHREADS+k];
          C_Units[j*NTHREADS+k] = C_Units[k*NTHREADS+j];
          C_Units[k*NTHREADS+j] = u;
        }
  }

  pair_sort_search();

  free(IDBpost);
  free(IDBsplit);

  { int k;

    for (k = NTHREADS*NTHREADS-1; k >= 0; k--)
      { free(C_Units[k].name);
        free(N_Units[k].name);
      }
    free(N_Units->buck);
    free(N_Units->bufr);
    free(C_Units);
    free(N_Units);
  }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
