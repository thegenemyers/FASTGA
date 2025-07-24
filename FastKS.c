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

#undef DEBUG_MERGE

static char *Usage[] = { "[-vk] [-b:<name>] [-T<int(8)>] [-P<dir($TMPDIR)>]",
                         "<source1:path>[<precursor>] <source2:path>[<precursor>]"
                       };

static int    VERBOSE;     //  -v: Verbose output
static int    NTHREADS;    //  -T
static char  *SORT_PATH;   //  -P
static int    KEEP;        //  -k
static FILE  *HFILE;       //  -b

static char *PATH1, *PATH2;   //  GDB & GIX are PATHx/ROOTx[GEXTNx|.gix]
static char *ROOT1, *ROOT2;
static char *GEXTN1, *GEXTN2;

static char *SPATH1, *SPATH2; //  Path name of source if TYPEx <= IS_GDB
static int   TYPE1,  TYPE2;   //  Type of source file (see DNAsource.h)

static int   KMER;         //  K-mer length and # of threads from genome indices

static int KBYTE;         // # of bytes for a k-mer table entry (both T1 & T2)
static int CBYTE;         // byte of k-mer table entry containing post count
static int LBYTE;         // byte of k-mer table entry containing lcp

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

typedef struct
  { Kmer_Stream *T1;
    Kmer_Stream *T2;
    int          tid;
    int          pbeg, pend;
    int64       *histu;
    int64       *histl;
  } SP;

static void *merge_thread(void *args)
{ SP *parm = (SP *) args;
  int tid          = parm->tid;
  int64 *histu     = parm->histu;
  int64 *histl     = parm->histl;
  Kmer_Stream *T1  = parm->T1;
  Kmer_Stream *T2  = parm->T2;

  int64   tbeg, tend, maxp;

  int     cpre;
  uint8  *ctop, *suf1;
  uint8  *cache;

  uint8   puint8;
  int     eorun, plen;
  uint8  *rcur, *rend;
  uint8  *vlcp[KMER+1];

  int     qcnt, pcnt;

#ifdef DEBUG_MERGE
  int64   Tdp;
  char   *tbuffer;
#endif

#ifdef DEBUG_MERGE
  tbuffer = Current_Kmer(T1,NULL);
#endif

  First_Kmer_Entry(T1);
  First_Kmer_Entry(T2);

  if (tid != 0)
    { GoTo_Kmer_Index(T1,T1->index[(parm->pbeg << 8) | 0xff]);
      GoTo_Kmer_Index(T2,T2->index[(parm->pbeg << 8) | 0xff]);
    }
  tend = T1->index[(parm->pend << 8) | 0xff];

  { int t, u;

    t = T2->cpre;
    if (t == 0)
      maxp = T2->index[t];
    else
      maxp = T2->index[t] - T2->index[t-1];
    u = (parm->pend << 8) | 0xff;
    for (t++; t <= u; t++)
      if (T2->index[t]-T2->index[t-1] > maxp)
        maxp = T2->index[t]-T2->index[t-1];
  }

  cache = Malloc((maxp+1)*KBYTE,"Allocating cache");
  if (cache == NULL)
    Clean_Exit(1);

  cpre  = -1;
  ctop  = cache;

  plen = 12;                         //  Keep the dumb compiler happy
  vlcp[plen] = rcur = rend = cache;
  eorun = 0;

  qcnt = -1;
  for (tbeg = T1->cidx; T1->cidx < tend; Next_Kmer_Entry(T1))
    { suf1 = T1->csuf;
#ifdef DEBUG_MERGE
      printf("Doing %s (%lld)\n",Current_Kmer(T1,tbuffer),T1->cidx); fflush(stdout);
#endif

      if (T1->cpre != cpre)  //  New prefix panel
        { uint8 *cp;

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

          for (cpre = T1->cpre; T2->cpre < cpre; Next_Kmer_Entry(T2))
            ;

#ifdef DEBUG_MERGE
          Tdp = T2->cidx;
          printf("Loading %lld %06x ...",Tdp,cpre); fflush(stdout);
#endif

          //  load cache with T2 entries = T1->cpre

          for (cp = cache; T2->cpre == cpre; cp += KBYTE)
            { memcpy(cp,T2->csuf,KBYTE);
              Next_Kmer_Entry(T2);
            }
          ctop = cp;
          ctop[LBYTE] = 11;

          //  if cache is empty then skip to next T1 entry with a greater prefix than cpre

          if (ctop == cache)
            { while (T1->cpre == cpre)
                Next_Kmer_Entry(T1);
#ifdef DEBUG_MERGE
              printf(" ... Empty => to %06x in T1\n",T1->cpre);
#endif
              continue;
            }

          //  start adpatermer merge for prefix cpre

          plen = 12;
          vlcp[plen] = rcur = rend = cache;
          eorun = 0;

#ifdef DEBUG_MERGE
          printf("... to %lld rcur = rend = %lld, eorun = 0, plen = 12\n",
                 T2->cidx,Tdp+(rcur-cache)/KBYTE);
          fflush(stdout);
#endif
        }

      // eorun = 0: rcur <= rend, n[1..plen] = rcur..rend, n[plen+1] < rend[plen+1]
      // eorun = 1: rcur <  rend, n[1..plen] = rcur..rend-1, lcp(rend) < plen 

      else
         { int nlcp;

           nlcp = suf1[LBYTE];
           if (nlcp > plen)
             goto pairs;                   // voids last
           else if (nlcp == plen)
             { if (eorun)                  // voids last
                 goto pairs;
             }
           else
             { if ( ! eorun)               // barrier
                 rend += KBYTE;
               while (rend[LBYTE] > nlcp)
                 rend += KBYTE;
               plen = rend[LBYTE];
               if (plen < nlcp)
                 { eorun = 1;
                   plen  = nlcp;
                   goto pairs;               //  not any good, shorter than last
                 }
               eorun = 0;                  //  good only if extended
               rcur  = rend;
             }
         }

       while (plen < KMER)
         { int h, m, c, d;

           h = cbyte[plen];
           m = mbyte[plen];
           c = suf1[h] & m;
           for (d = rend[h]&m; d < c; d = rend[h]&m)
             { rend += KBYTE;
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
       rend += KBYTE;
       eorun = 1;

       //  Get pairs;

    pairs:

#ifdef DEBUG_MERGE
      printf("-> %d[%lld,%lld,%lld] %d\n",
             plen,Tdp+(vlcp[plen]-cache)/KBYTE,Tdp+(rcur-cache)/KBYTE,Tdp+(rend-cache)/KBYTE,eorun);
      fflush(stdout);
#endif

      { int    freq;
        uint8 *vcp;

        histl[plen] += 1;
        if (HFILE != NULL)
          { puint8 = (uint8) plen;
            fwrite(&puint8,sizeof(uint8),1,HFILE);
          }
        
        freq = 0;
        vcp = vlcp[plen];
        if (vcp < rend)
          { if (rend-vcp > KBYTE)
              continue;
            freq = vcp[CBYTE];
          }
        if ( ! eorun)
          { if (rend[KBYTE+LBYTE] >= plen)
              continue;
            freq += rend[CBYTE];
          }
#ifdef DEBUG_MERGE
        printf("    w %d hits\n",freq);
        fflush(stdout);
#endif

        if (freq == 1 && suf1[CBYTE] == 1)
          { if (suf1[LBYTE] < plen && suf1[KBYTE+LBYTE] < plen)
              histu[plen] += 1;
          }
      }
    }

  free(cache);
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

static void adaptamer_merge(Kmer_Stream *T1, Kmer_Stream *T2, int64 g1len)
{ SP         parm[NTHREADS];
#ifndef DEBUG_MERGE
  pthread_t  threads[NTHREADS];
#endif
  int64     *histu, *histl;
  int        i, t;

  (void) g1len;

  { Kmer_Stream *tp;
    uint8       *ent;
    int   t;
    int64 p;

    if (T1->nels > T2->nels)
      tp = T1;
    else
      tp = T2;
    ent = Current_Entry(tp,NULL);
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

  histu = Malloc(sizeof(int64)*(KMER+1)*NTHREADS,"Unique k-mer histogram");
  histl = Malloc(sizeof(int64)*(KMER+1)*NTHREADS,"Adapatmer length histogram");
  bzero(histu,sizeof(int64)*(KMER+1)*NTHREADS);
  bzero(histl,sizeof(int64)*(KMER+1)*NTHREADS);
  for (i = 0; i < NTHREADS; i++)
    { parm[i].histu = histu + i*(KMER+1);
      parm[i].histl = histl + i*(KMER+1);
    }

  parm[0].T1 = T1;
  parm[0].T2 = T2;
  for (i = 1; i < NTHREADS; i++)
    { parm[i].T1 = Clone_Kmer_Stream(T1);
      parm[i].T2 = Clone_Kmer_Stream(T2);
    }

  for (i = 0; i < NTHREADS; i++)
    parm[i].tid   = i;

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

  for (i = 1; i < NTHREADS; i++)
    { for (t = 1; t <= KMER; t++)
        { histu[t] += parm[i].histu[t];
          histl[t] += parm[i].histl[t];
        }
    }

  if (VERBOSE)
    { fprintf(stderr,"\r    Completed 100%%\n");
      fflush(stderr);
    }

  printf("   K:  unique-mers   adapt-mers\n");
  for (t = 1; t <= KMER; t++)
    printf(" %2d: %10lld %10lld\n",t,histu[t],histl[t]);

  for (i = NTHREADS-1; i >= 1; i--)
    { Free_Kmer_Stream(parm[i].T1);
      Free_Kmer_Stream(parm[i].T2);
    }
  Free_Kmer_Stream(T1);
  Free_Kmer_Stream(T2);
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
  GDB _gdb1, *gdb1 = &_gdb1;
  GDB _gdb2, *gdb2 = &_gdb2;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("FastKS");

    SORT_PATH = getenv("TMPDIR");
    if (SORT_PATH == NULL)
      SORT_PATH = ".";
    NTHREADS = 8;
    HFILE    = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vk")
            break;
          case 'b':
            if (strncmp(argv[i]+1,"b:",2) == 0)
              { HFILE = fopen(argv[i]+3,"w");
                if (HFILE == NULL)
                  { fprintf(stderr,"%s: Cannot open %s for output\n",
                                   Prog_Name,argv[i]+3);
                    exit (1);
                  }
                break;
              }
            fprintf(stderr,"%s: Do not recognize option %s\n",Prog_Name,argv[i]);
            exit (1);
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

    if (argc != 3)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"         <precursor> = .gix | .1gdb | <fa_extn> | <1_extn>\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"             <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"             <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -b: output adaptamer lengths for each A entry.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -k: Keep any generated .1gdb's and .gix's.\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        fprintf(stderr,"      -P: Directory to use for temporary files.\n");
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
    PATH2 = argv[2];
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
                        VERBOSE?" -v":"",NTHREADS,SORT_PATH,10,tpath1);
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
                        VERBOSE?" -v":"",NTHREADS,SORT_PATH,10,tpath2);
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

  T1 = Open_Kmer_Stream(Catenate(PATH1,"/",ROOT1,".gix"),2);
  T2 = Open_Kmer_Stream(Catenate(PATH2,"/",ROOT2,".gix"),2);
  if (T1 == NULL)
    { fprintf(stderr,"%s: Cannot find genome index for %s/%s.gix\n",Prog_Name,PATH1,ROOT1);
      Clean_Exit(1);
    }
  if (T2 == NULL)
    { fprintf(stderr,"%s: Cannot find genome index for %s/%s.gix\n",Prog_Name,PATH2,ROOT2);
      Clean_Exit(1);
    }
  
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

    GEXTN2 = ".gdb";
    fname = Catenate(PATH1,"/",ROOT1,GEXTN2);
    if ((file = fopen(fname,"r")) == NULL)
      GEXTN2 = ".1gdb";
    else
      fclose(file);

    if (Read_GDB(gdb2,Catenate(PATH2,"/",ROOT2,GEXTN2)) < 0)
      Clean_Exit(1);
    short_GDB_fix(gdb2);
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

      nfiles = (NTHREADS+3)*2*NTHREADS + tid + 100 ; // 100 allows for a few wrapping processes
      getrlimit(RLIMIT_NOFILE,&rlp);
      if (nfiles > rlp.rlim_max)
        { fprintf(stderr,"\n%s: Cannot open %lld files simultaneously\n",Prog_Name,nfiles);
          Clean_Exit(1);
        }
      rlp.rlim_cur = nfiles;
      setrlimit(RLIMIT_NOFILE,&rlp);
    }

  KBYTE = T2->pbyte;
  CBYTE = T2->hbyte;
  LBYTE = CBYTE+1;

  if (VERBOSE)
    { fprintf(stderr,"\n  Using %d threads\n",NTHREADS);
      fflush(stderr);
    }

  adaptamer_merge(T1,T2,gdb1->seqtot);

  if (VERBOSE)
    TimeTo(stderr,1,0);

  Close_GDB(gdb2);
  Close_GDB(gdb1);

  free(SORT_PATH);

  if (TYPE2 <= IS_GDB)
    free(SPATH2);
  if (TYPE1 <= IS_GDB)
    free(SPATH1);

  fclose(HFILE);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  Clean_Exit(0);
}
