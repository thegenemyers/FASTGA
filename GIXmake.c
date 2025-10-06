/*******************************************************************************************
 *
 *  Produce a k-mer index of genome's contigs suitable for finding adaptamer seed matches.
 *    As such only k-mers whose count is less than the adaptamer frequency cutoff are in
 *    the index which consists of a FASTK k-mer table and an associated position list that
 *    contains the positions at which each k-mer occurs in order of the k-mers of the table.
 *    What would normally be the count field for a k-mer contains the lcp of the k-mer with
 *    its predecessor, and the # of positions at which the k-mer occurs, each in a byte.
 *
 *  Author:  Gene Myers
 *  Date  :  February 2023
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include <sys/resource.h>

#include "libfastk.h"
#include "GDB.h"

#undef  DEBUG_THREADS
#undef  DEBUG_MAP
#undef  DEBUG_SETUP
#undef  DEBUG_SORT

  //  These constants determine the syncmer (TMER,SMER) and whether to add a spaced seed (defunct)

#define TMER 12        //  Syncmer len
#define SMER  8        //  S-mer len
#define SOFF  4        //  TMER-SMER

static char *Usage[] =
    { "[-v] [-L:<log:path>] [-T<int(8)>] [-P<dir($TMPDIR)>] [-k<int(40)]",
      "( <source:path>[.1gdb]  |  <source:path>[<fa_extn>|<1_extn>] [<target:path>[.gix]] )"
    };

static int   VERBOSE;    //  -v
static char *LOG_PATH;   //  -L path name for log file
static FILE *LOG_FILE;   //     file handle for log
static char *SORT_PATH;  //  -P
static char *TPATH;
static char *TROOT;
static char *POST_NAME;

static int NTHREADS;   //  by default 8
static int NPARTS;     //  # of parts for sorts (in [8,64], power of 2)
static int KMER;       //  by default 40, must be >= 12 and divisible by 4
static int KBYTES;     //  Bytes for 2-bit compress k-mer (KMER/4)
static int MASK;       //  Set if gdb has masks
static int MBYTES;     //  Bytes for masked k-mer (KBYTES+1)

#define BUFF_MAX   1000000
#define SCAN_MAX  10000000    //  Must be divisible by 4
#define NUM_BUCK      1024

static int   *Perm;        //  Size sorted permutation of contigs
static int   *InvP;        //  Inverse of Perm

static int    PostBytes;   //  # of bytes needed for a position
static int    ContBytes;   //  # of bytes for a contig + sign bit

                           //  For p in [0,NTHREADS):
static int   *DBsplit;     //    DB split: contig [DBsplit[p],DBsplit[p+1])
static int64 *DBpost;      //    DB post: post of start of contig DBsplit[p] is DBpost[p]
                           //  For p in [0,NPARTS):
static int   *Ksplit;      //    Kmer split:  1st bytes [Ksplit[p],Ksplit[p+1])

static int64 **Buckets;    //  NTHREADS x NUM_BUCK array:
                           //  B[i][j] = # of posts whose canonical k-mer's 1st byte
                           //            is j from block DBsplit[i],DBsplit[i+1]
                           //  Also used as "fingers" for loading initial sort array.

static int *Units;    //  NTHREADS*NPARTS IO units for distribution and import & k-mer table parts

static void *MyBlock; //  1 large block for all other work structures save sorting array

static int    Comp[256];   //  DNA complement of packed byte
static int    Select[NUM_BUCK]; //  1st k-mer byte bucket -> block (of NTHREAD)

static int    TMap[256]   //  4-bp code map
 = { 0xff, 0xd4, 0xf5, 0xfd, 0xe4, 0xad, 0x21, 0xa5, 0xed, 0x64, 0xbf, 0xa9, 0xf3, 0x70, 0xd6, 0xf0,
     0xca, 0x89, 0xcb, 0xc9, 0x82, 0x9d, 0x13, 0x79, 0x0a, 0x0f, 0x25, 0x19, 0x3e, 0x47, 0xa3, 0xa8,
     0xf9, 0x5e, 0xe8, 0xa1, 0xb0, 0x71, 0x1d, 0x8c, 0xde, 0x69, 0xe7, 0x7c, 0x56, 0x3f, 0x90, 0xa4,
     0xeb, 0x45, 0x59, 0xf1, 0x97, 0x4c, 0x08, 0xa0, 0xb8, 0x4a, 0x86, 0xc8, 0xcd, 0x98, 0x7d, 0xfc,
     0xef, 0x4d, 0x83, 0x7e, 0xdc, 0x66, 0x2b, 0x8e, 0xe0, 0xa7, 0xd0, 0xa2, 0x88, 0x5f, 0x7f, 0xd9,
     0x9b, 0x78, 0xd1, 0x8b, 0xc3, 0x8f, 0x2d, 0xe6, 0x18, 0x27, 0x2c, 0x24, 0x94, 0xb7, 0xce, 0xbd,
     0x0d, 0x04, 0x1c, 0x09, 0x16, 0x23, 0x00, 0x1e, 0x1a, 0x29, 0x2e, 0x15, 0x01, 0x10, 0x2a, 0x20,
     0xbe, 0x31, 0x43, 0x58, 0xc2, 0xaa, 0x1f, 0xe5, 0xc5, 0x9e, 0xcf, 0xc6, 0x68, 0xb2, 0x80, 0xf4,
     0xf8, 0x53, 0xb6, 0x93, 0x76, 0x37, 0x11, 0x40, 0xda, 0x51, 0xba, 0x46, 0x42, 0x30, 0x60, 0x6d,
     0x5c, 0x39, 0x9f, 0x48, 0x6c, 0x62, 0x28, 0x67, 0x06, 0x12, 0x26, 0x0e, 0x33, 0x50, 0xa6, 0x63,
     0xdd, 0x3b, 0xab, 0x4b, 0x72, 0x5b, 0x22, 0x6f, 0xb4, 0x61, 0x92, 0x99, 0x36, 0x38, 0x65, 0xac,
     0x4f, 0x2f, 0x32, 0x44, 0x54, 0x3c, 0x03, 0x5d, 0x73, 0x3a, 0x77, 0x84, 0x8d, 0x4e, 0x49, 0xd2,
     0xfb, 0x91, 0x6a, 0xcc, 0x8a, 0x35, 0x02, 0x55, 0x7a, 0x34, 0x96, 0x3d, 0xd3, 0x41, 0x85, 0xf2,
     0xb1, 0x75, 0xc4, 0xb5, 0xbb, 0xb3, 0x1b, 0xd5, 0x07, 0x05, 0x17, 0x0b, 0x7b, 0xd7, 0xdf, 0xea,
     0xe3, 0x57, 0xc0, 0x95, 0x9c, 0x6e, 0x14, 0xae, 0xb9, 0x6b, 0xc1, 0x81, 0x87, 0x74, 0xd8, 0xe2,
     0xec, 0x52, 0xbc, 0xe9, 0xe1, 0xdb, 0x0c, 0xf7, 0xaf, 0x5a, 0x9a, 0xc7, 0xfa, 0xf6, 0xee, 0xfe
   };

typedef struct
  { int   beg;
    int   end;
    int64 off;
  } Range;

extern void msd_sort(uint8 *array, int64 nelem, int rsize, int ksize,
                     int64 *part, int beg, int end, int nthreads);


/***********************************************************************************************
 *
 *   DISTRIBUTION PHASE:  Output signed posts for each 1/NTHREAD of the data and the k-mers
 *        At completion file ._post.<jobid>.<i*NTHREADS+j>.idx contains the signed posts from
 *             contigs (DBsplit[j],DBplist[j+1]) whose k-mers 1st byte is in Ksplit[i],Ksplit[i+1].
 *
 **********************************************************************************************/

#ifdef DEBUG_SYNCMERS

static int is_syncmer(int x)
{ int y, z, m, b;
  int q, s;
 
  b = 256;
  q = 16;
  for (s = 16; s >= 0; s -= 2)
    { y = (x>>s) & 0xff;
      z = Comp[y];
      if (TMap[y] < TMap[z])
        m = TMap[y];
      else
        m = TMap[z];
      if (m < b)
        { q = s;
          b = m;
        }
    }
  return (q == 16 || m <= b);
}

#endif

typedef struct
  { int     tid;
    GDB     gdb;
    int64  *buck;
    uint8  *seq;
    uint8  *buff;
  } DP;

  //  Determine distribution of 1st 4-mers of the 40-mers going into index

static void *sample_thread(void *args)
{ DP    *parm  = (DP *) args;
  int    tid  = parm->tid;
  int64 *buck = parm->buck;
  GDB   *gdb  = &parm->gdb;
  char  *SEQ  = (char *) (parm->seq + 1);

  int     r, ren;

  bzero(buck,sizeof(int64)*NUM_BUCK);

  ren  = DBsplit[tid+1];
  for (r = DBsplit[tid]; r < ren; r++)
    { int   len;
      int    beg, end;

      int    i, n;
      uint8 *seq, *sn;
      int    ne[8];
      int    nq[8], cq[8];
      int    mzr[8];
      int    min4, pos4;

      if (gdb->contigs[r].boff < 0)
        continue;

      len = gdb->contigs[r].clen;
      if (len > SCAN_MAX)
        end = SCAN_MAX;
      else
        end = len;
   
      seq = (uint8 *) Get_Contig_Piece(gdb,r,0,end,NUMERIC,SEQ);   //  Load the contig part

      //  Set up for scan

#ifdef DEBUG_MAP
      for (i = 0; i < 3; i++)
        printf(" %4d: %02x\n",i,SEQ[i]);
#endif

      sn = seq+3;
      n = (seq[0] << 4) | (seq[1] << 2) | seq[2];
      for (i = 0; i < 4; i++)
        { int c;

          ne[i+4] = n  = ((n<<2) | sn[i]) & 0xff;
          c = Comp[n];
          nq[i] = TMap[n];
          cq[i] = TMap[c];
#ifdef DEBUG_MAP
          printf(" %4d: %02x | %02x %02x | %02x %02x\n",
                 i,sn[i],n,Comp[n],nq[i],cq[i]);
#endif
        }
      sn += 4;
      min4 = 0x10000;
      pos4 = 0;
      for (i = 0; i < SOFF; i++)
        { int nh, ch, c;
          int mn, mc, mz;

          ne[i] = n  = ((n<<2) | sn[i]) & 0xff;
          c  = Comp[n];
          nh = TMap[n];
          ch = TMap[c];
          mn = (nq[i] << 8) | nh;
          mc = cq[i] | (ch << 8);
          nq[i] = nh;
          cq[i] = ch;
          if (mn < mc)
            mzr[i] = mz = mn;
          else
            mzr[i] = mz = mc;
#ifdef DEBUG_MAP
          printf(" %4d: %02x | %02x %02x | %02x %02x | %04x %04x %04x\n",
                 i,sn[i],n,c,nh,ch,mn,mc,mz);
#endif

          if (mz < min4)
            { min4 = mz;
              pos4 = i;
            }
        }

      beg = SOFF;
      while (1)
        { end -= SMER;

#ifdef DEBUG_MAP
          printf("------- %d-%d\n",beg,end);
#endif
          for (i = beg; i <= end; i++)
            { int iq, jq;
              int nh, ch;
              int mn, mc, mz;
              int j, w, c, p;

              iq = (i&0x7);
              w = ne[iq];
              ne[iq] = n = ((n<<2) | sn[i]) & 0xff;
              c  = Comp[n];
              nh = TMap[n];
              ch = TMap[c];

              jq = i&0x3;
              mn = (nq[jq] << 8) | nh;
              mc = cq[jq] | (ch << 8);
              nq[jq] = nh;
              cq[jq] = ch;
              if (mn < mc)
                mz = mzr[jq] = mn;
              else
                mz = mzr[jq] = mc;

#ifdef DEBUG_MAP
              printf(" %4d: %02x | %02x %02x | %02x %02x | %04x %04x %04x :: %04x(%d)\n",
                      i,sn[i],n,c,nh,ch,mn,mc,mz,min4,pos4);
#endif

              if (mz < min4)              //  right-end of 12-syncmer
                { min4 = mz;
                  pos4 = i;
#ifdef DEBUG_MAP
                  printf("      Hit R");
#endif
                }
              else if (pos4 == i-SOFF)          //   left-end of 12-syncmer
                { min4 = mzr[(++pos4)&0x3];
                  for (j = pos4+1; j <= i; j++)
                    if (mzr[j&0x3] < min4)
                      { min4 = mzr[j&0x3];
                        pos4 = j;
                      }
#ifdef DEBUG_MAP
                  printf("      Hit L");
#endif
                }
              else if (mz > min4)
                continue;
#ifdef DEBUG_MAP
              else
                printf("      Hit RE");
#endif

              p = ne[(i+4)&0x7];
              buck[w<<2|p>>6] += 1;
              buck[c<<2|(3-(p&0x3))] += 1;
            }
          end += SMER;

          if (end == len)
            break;
          beg = end;
          end += SCAN_MAX;
          if (end > len)
            end = len;

          beg -= (SMER-1);
          sn = (uint8 *) Get_Contig_Piece(gdb,r,beg+(SMER-1),end,NUMERIC,(char *) SEQ) - beg;
        }
    }

  return (NULL);
}

typedef struct
  { uint8 *ptr;
    uint8 *end;
    uint8 *buffer;
    uint64 last;
    int    out;
    int    inum;
  } Packet;

//  The thread scans seq and sends those posts assigned to file tid*Nthreads + ? to their
//    designated file relative to the last emission in compressed form:
//       x0   -> byte,   6-bit post delta with sign x
//       x10  -> short, 13-bit ...
//       x110 -> short, 28-bit ...
//       x111 -> 0x10000000 spacer

void push(Packet *pack, int64 post, int comp)
{ uint8 *b = pack->ptr;
  uint8 *p;
  uint32 x;
  uint8 *xbyte = (uint8 *) &x;

  x = post - pack->last;
  if (x < 0x3f)
    { if (comp)
        *b++ = 0x80 | x;
      else
        *b++ = x;
    }
  else if (x < 0x1fff)
    { if (comp)
        x |= 0xc000;
      else
        x |= 0x4000;
      *b++ = xbyte[1];
      *b++ = xbyte[0];
    }
  else
    { while (x >= 0x10000000)
        { if (comp)
            *b++ = 0xf0;
          else
            *b++ = 0x70;
          if (b >= pack->end)
            { p = pack->buffer;
              if (write(pack->out,p,b-p) < 0)
                { fprintf(stderr,"%s: IO write to file %s%d.idx failed\n",
                                  Prog_Name,POST_NAME,pack->inum);
                   exit (1);
                 }
              b = p;
            }
          x -= 0x10000000;
        }
      if (comp)
        x |= 0xe0000000;
      else
        x |= 0x60000000;
      *b++ = xbyte[3];
      *b++ = xbyte[2];
      *b++ = xbyte[1];
      *b++ = xbyte[0];
    }
  if (b >= pack->end)
    { p = pack->buffer;
      if (write(pack->out,p,b-p) < 0)
        { fprintf(stderr,"%s: IO write to file %s%d.idx failed\n",
                         Prog_Name,POST_NAME,pack->inum);
          exit (1);
        }
      b = p;
    }
  pack->ptr  = b;
  pack->last = post;
}

static void *scan_thread(void *args)
{ DP    *parm  = (DP *) args;
  int    tid  = parm->tid;
  int64 *buck = parm->buck;
  GDB   *gdb  = &parm->gdb;
  char  *SEQ  = (char *) (parm->seq + 1);
  uint8 *buff = parm->buff;

  Packet packs[NPARTS];

  int     lnK, KMT;
  int     r, ren;
  int64   post;

  bzero(buck,sizeof(int64)*NUM_BUCK);

  for (r = 0; r < NPARTS; r++)
    { packs[r].buffer = buff+r*BUFF_MAX;
      packs[r].ptr    = packs[r].buffer;
      packs[r].end    = packs[r].buffer + (BUFF_MAX-4);
      packs[r].out    = Units[r*NTHREADS+tid];
      packs[r].inum   = NTHREADS*tid+r;
      packs[r].last   = DBpost[tid];
    }

  KMT  = KMER-TMER;
  post = DBpost[tid];
  ren  = DBsplit[tid+1];
  for (r = DBsplit[tid]; r < ren; r++)
    { int   len;
      int    beg, end;

      int    i, n;
      uint8 *seq, *sn;
      int    ne[8];
      int    nq[8], cq[8];
      int    mzr[8];
      int    min4, pos4;

      if (gdb->contigs[r].boff < 0)
        continue;

      len = gdb->contigs[r].clen;
      if (len > SCAN_MAX)
        end = SCAN_MAX;
      else
        end = len;
   
      seq = (uint8 *) Get_Contig_Piece(gdb,r,0,end,NUMERIC,SEQ);   //  Load the contig part

      //  Set up for scan

      lnK = len-KMER;

#ifdef DEBUG_MAP
      for (i = 0; i < 3; i++)
        printf(" %4d: %02x\n",i,SEQ[i]);
#endif

      sn = seq+3;
      n = (seq[0] << 4) | (seq[1] << 2) | seq[2];
      for (i = 0; i < 4; i++)
        { int c;

          ne[i+4] = n  = ((n<<2) | sn[i]) & 0xff;
          c = Comp[n];
          nq[i] = TMap[n];
          cq[i] = TMap[c];
#ifdef DEBUG_MAP
          printf(" %4d: %02x | %02x %02x | %02x %02x\n",
                 i,sn[i],n,Comp[n],nq[i],cq[i]);
#endif
        }
      sn += 4;
      min4 = 0x10000;
      pos4 = 0;
      for (i = 0; i < SOFF; i++)
        { int nh, ch, c;
          int mn, mc, mz;

          ne[i] = n  = ((n<<2) | sn[i]) & 0xff;
          c  = Comp[n];
          nh = TMap[n];
          ch = TMap[c];
          mn = (nq[i] << 8) | nh;
          mc = cq[i] | (ch << 8);
          nq[i] = nh;
          cq[i] = ch;
          if (mn < mc)
            mzr[i] = mz = mn;
          else
            mzr[i] = mz = mc;
#ifdef DEBUG_MAP
          printf(" %4d: %02x | %02x %02x | %02x %02x | %04x %04x %04x\n",
                 i,sn[i],n,c,nh,ch,mn,mc,mz);
#endif

          if (mz < min4)
            { min4 = mz;
              pos4 = i;
            }
        }

      beg = SOFF;
      while (1)
        { end -= SMER;

#ifdef DEBUG_MAP
          printf("------- %d-%d\n",beg,end);
#endif
          for (i = beg; i <= end; i++)
            { int iq, jq;
              int nh, ch;
              int mn, mc, mz;
              int j, w, c, p;

              iq = (i&0x7);
              w = ne[iq];
              ne[iq] = n = ((n<<2) | sn[i]) & 0xff;
              c  = Comp[n];
              nh = TMap[n];
              ch = TMap[c];

              jq = i&0x3;
              mn = (nq[jq] << 8) | nh;
              mc = cq[jq] | (ch << 8);
              nq[jq] = nh;
              cq[jq] = ch;
              if (mn < mc)
                mz = mzr[jq] = mn;
              else
                mz = mzr[jq] = mc;

#ifdef DEBUG_MAP
              printf(" %4d: %02x | %02x %02x | %02x %02x | %04x %04x %04x :: %04x(%d)\n",
                      i,sn[i],n,c,nh,ch,mn,mc,mz,min4,pos4);
#endif

              if (mz < min4)              //  right-end of 12-syncmer
                { min4 = mz;
                  pos4 = i;
#ifdef DEBUG_MAP
                  printf("      Hit R");
#endif
                }
              else if (pos4 == i-SOFF)          //   left-end of 12-syncmer
                { min4 = mzr[(++pos4)&0x3];
                  for (j = pos4+1; j <= i; j++)
                    if (mzr[j&0x3] < min4)
                      { min4 = mzr[j&0x3];
                        pos4 = j;
                      }
#ifdef DEBUG_MAP
                  printf("      Hit L");
#endif
                }
              else if (mz > min4)
                continue;
#ifdef DEBUG_MAP
              else
                printf("      Hit RE");
#endif

	      j = i-SOFF;
              p = ne[(i+4)&0x7];
              if (j <= lnK)
                { w = (w<<2|p>>6);
                  buck[w] += 1;
#ifdef DEBUG_MAP
                  printf(" @ %d >> %02x %d ",j,w,Select[w]);
#endif
                  push(packs+Select[w],post+j,0);
                }
              if (j >= KMT)
                { c = (c<<2|(3-(p&0x3)));
                  buck[c] += 1;
#ifdef DEBUG_MAP
                  printf(" >> %02x %d",c,Select[c]);
#endif
                  push(packs+Select[c],post+j,1);
                }
#ifdef DEBUG_MAP
              printf("\n");
#endif
            }
          end += SMER;

          if (end == len)
            break;
          beg = end;
          end += SCAN_MAX;
          if (end > len)
            end = len;

          beg -= (SMER-1);
          sn = (uint8 *) Get_Contig_Piece(gdb,r,beg+(SMER-1),end,NUMERIC,(char *) SEQ) - beg;
        }
      post += len;
    }

  for (r = 0; r < NPARTS; r++)
    if (packs[r].ptr > packs[r].buffer)
      write(packs[r].out,packs[r].buffer,packs[r].ptr-packs[r].buffer);

  return (NULL);
}

//  Create a thread for each DB section, and have it distribute syncmer posts to NTHREAD files
//     according to the 1st byte of its k-mer.
  
void distribute(GDB *gdb)
{ DP      parm[NTHREADS];
  uint8  *buff, *seq;
#ifndef DEBUG_THREADS
  pthread_t threads[NTHREADS];
#endif
  int i;

  //  60MB per thread for work

  seq  = (uint8 *) MyBlock;;
  buff = seq + (SCAN_MAX+8)*NTHREADS;

  for (i = 0; i < NTHREADS; i++)
    { parm[i].tid  = i;
      parm[i].seq  = seq + i*(SCAN_MAX+8);
      parm[i].buff = buff + i*NPARTS*BUFF_MAX;
      parm[i].buck = Buckets[i];
      parm[i].gdb  = *gdb;
      if (i > 0)
        { parm[i].gdb.seqs = fopen(gdb->seqpath,"r");
          if (parm[i].gdb.seqs == NULL)
            { fprintf(stderr,"%s: Cannot open another copy of GDB\n",Prog_Name);
              exit (1);
            }
        }
    }

#ifdef DEBUG_THREADS
  for (i = 0; i < NTHREADS; i++)
    sample_thread(parm+i);
#else
  for (i = 1; i < NTHREADS; i++)
    pthread_create(threads+i,NULL,sample_thread,parm+i);
  sample_thread(parm);
  for (i = 1; i < NTHREADS; i++)
    pthread_join(threads[i],NULL);
#endif

  { int64 t, *buck;
    int   n;

    buck = Buckets[0];
    for (i = 1; i < NTHREADS; i++)
      for (n = 0; n < NUM_BUCK; n++)
        buck[n] += Buckets[i][n];

#ifdef DEBUG_SETUP
    printf("\nBuckets\n");
    for (i = 0; i < NUM_BUCK; i++)
      printf(" %3d: %10lld\n",i,buck[i]);
#endif

    for (i = 1; i < NUM_BUCK; i++)
      buck[i] = buck[i-1] + buck[i];

    Ksplit[0] = 0;
    n = 1;
    t = buck[NUM_BUCK-1]/NPARTS;
    for (i = 0; i < NUM_BUCK; i++)
      { if (buck[i] >= t)
          { if (buck[i]-t > t-buck[i-1])
              { Select[i] = n;
                Ksplit[n] = i;
              }
            else
              { Select[i] = n-1;
                Ksplit[n] = i+1;
              }
            n += 1;
            t = (n*buck[NUM_BUCK-1])/NPARTS;
          }
        else
          Select[i] = n-1;
      }
    Ksplit[NPARTS] = NUM_BUCK;

#ifdef DEBUG_SETUP
    printf("Prepatory %lld\n",buck[NUM_BUCK-1]);
    for (i = 0; i <= NPARTS; i++)
      printf(" %3d: %3d %lld\n",i,Ksplit[i],buck[Ksplit[i+1]]-buck[Ksplit[i]]);
    for (i = 0; i < NUM_BUCK; i++)
      printf(" %3d: %8lld  %d\n",i,buck[i],Select[i]);
    fflush(stdout);
#endif
  }

#ifdef DEBUG_THREADS
  for (i = 0; i < NTHREADS; i++)
    scan_thread(parm+i);
#else
  for (i = 1; i < NTHREADS; i++)
    pthread_create(threads+i,NULL,scan_thread,parm+i);
  scan_thread(parm);
  for (i = 1; i < NTHREADS; i++)
    pthread_join(threads[i],NULL);
#endif

  for (i = 1; i < NTHREADS; i++)
    fclose(parm[i].gdb.seqs);
}


/***********************************************************************************************
 *
 *   SORT PHASE:  Read in posts for a given k-mer 1st byte range, recompute their 2-bit compressed
 *        canonical k-mer, and place them and their k-mer in a sort array in order of 1st byte
 *        using the Bucket array as a set of "fingers".  Then complete the sort with a min-first
 *        radix sort.  Finally output the k-mer table and index for each k-mer range.
 *
 **********************************************************************************************/


#ifdef DEBUG_SORT

static void print_table(uint8 *array, int swide, int64 nelem)
{ static char  DNA[4] = { 'a', 'c', 'g', 't' };
  static char *fmer[256], _fmer[1280];
  static int   first = 1;

  int64 x, post, cont, flag, end;
  int   k;
  uint8 *pust = (uint8 *) &post;
  uint8 *cust = (uint8 *) &cont;

  if (first)
    { char *t;
      int   i, l3, l2, l1, l0;

      first = 0;

      i = 0;
      t = _fmer;
      for (l3 = 0; l3 < 4; l3++)
       for (l2 = 0; l2 < 4; l2++)
        for (l1 = 0; l1 < 4; l1++)
         for (l0 = 0; l0 < 4; l0++)
           { fmer[i] = t;
             *t++ = DNA[l3];
             *t++ = DNA[l2];
             *t++ = DNA[l1];
             *t++ = DNA[l0];
             *t++ = 0;
             i += 1;
           }
    }

  flag = (0x1ll << (8*ContBytes-1));
  post = 0;
  cont = 0;
  end  = nelem * swide;
  for (x = 0; x < end; x += swide)
    { printf(" %10lld %2d: ",x,array[x]);
      for (k = 1; k < KBYTES; k++)
        printf("%s",fmer[array[x+k]]);
      if (MASK)
        printf("(%d)",array[x+KBYTES});
      for (k = 0; k < PostBytes; k++)
        pust[k] = array[x+MBYTES+k];
      for (k = 0; k < ContBytes; k++)
        cust[k] = array[x+MBYTES+PostBytes+k];
      if (cont & flag)
        printf(" - %11lld (%d)\n",post,cont-flag);
      else
        printf(" + %11lld (%d)\n",post,cont); 
      fflush(stdout);
    }
}

#endif

typedef struct
  { GDB      gdb;
    int      tid;
    int      inum;
    int      in;
    int      swide;
    char    *seq;
    uint8   *sarr;
    uint8   *buff;
  } SP;

//  Read post file, uncompressing it, recomputing the canonical k-mer at its absolute
//   position and loading the k-mer and the signed post into the initial soring array
//   for the current 1st byte panel.

static void *setup_thread_plain(void *args)
{ SP *parm = (SP *) args;
  int in        = parm->in;
  int tid       = parm->tid;
  int swide     = parm->swide;
  char  *SEQ    = parm->seq+1;
  uint8 *sarr   = parm->sarr;
  uint8 *buffer = parm->buff;
  GDB   *gdb    = &(parm->gdb);
  int64 *buck   = Buckets[tid];

  int    iamt;
  int    len;
  uint8 *neq, *keq;
  int    ncntg, inv;
  int64  post, bost, last, lost;
  int64  cont, nont, flag;
  uint8 *lust = (uint8 *) (&last);
  uint8 *bust = (uint8 *) (&bost);
  uint8 *cust = (uint8 *) (&cont);
  uint8 *nust = (uint8 *) (&nont);
  int64  nextpost, basepost;
  uint8 *bend, *btop, *b;

  flag = (0x1ll << (8*ContBytes-1));

  if (lseek(in,0,SEEK_SET) < 0)
    { fprintf(stderr,"%s: Rewind of file %s%d.idx failed\n",
                     Prog_Name,POST_NAME,parm->inum);
      exit (1);
    }

  iamt = read(in,buffer,BUFF_MAX);
  if (iamt < 0)
    { fprintf(stderr,"%s: IO read to file %s%d.idx failed\n",
                     Prog_Name,POST_NAME,parm->inum);
      exit (1);
    }
  bend = buffer + iamt;
  if (bend-buffer < BUFF_MAX)
    btop = bend;
  else
    btop = bend-4;
  b = buffer;
  if (btop <= buffer)   //  Nothing to do
    { close(in);
      return (NULL);
    }

  ncntg = DBsplit[tid];
  post  = DBpost[tid];
  nextpost = post;
  basepost = post;
  while (1)
    { { int   v;
        int64 extra;

        v = *b++;                    //  Get next encoded post, uncompress it and "un-delta"
        inv = ((v & 0x80) != 0);
        if (v & 0x40)
          if (v & 0x20)
            { extra = 0;
              while (v & 0x1)
                { extra -= 0x10000000;
                  v = *b++;
                }
              last = 0;
              lust[3] = (v & 0x0f); 
              lust[2] = *b++;
              lust[1] = *b++;
              lust[0] = *b++;
              last += extra;
            }
          else
            { last = 0;
              lust[1] = (v & 0x1f); 
              lust[0] = *b++;
            }
        else
          last = (v & 0x3f);
        post += last;
      }

      if (post >= nextpost)    //  if position is not in current contig, then fetch the next one
        {
          do
            { len = gdb->contigs[ncntg].clen;
              basepost = nextpost;
              if (gdb->contigs[ncntg].boff >= 0)
                nextpost += len;
              ncntg += 1;
            }
          while (post >= nextpost);

          cont = InvP[ncntg-1];
          nont = cont | flag;
          lost = 0;
        }

      { int    i;        //  load the k-mer / post pair into the next available slot in the
        uint8 *n, *x;    //    initial sort array using the bucket index

        bost = post-basepost;

        if (bost >= lost)   //  Get another buffer full of the current contig
          { int64  top, bot;
            int    n;

            lost = bost + (SCAN_MAX - (2*KMER-8));
            if (lost > len)
              lost = len;
            bot = bost-(KMER-8);
            top = lost+KMER;
            if (top > len)
              top = len;
            if (bot < 0)
              bot = 0;
            neq = (uint8 *) Get_Contig_Piece(gdb,ncntg-1,bot,top,NUMERIC,SEQ) - bot;

            keq = neq+3;
            top -= 3;
            n = (neq[bot] << 4) | (neq[bot+1] << 2) | neq[bot+2];
            for (i = bot; i < top; i++)
              neq[i] = n = (((n << 2) | keq[i]) & 0xff);
            keq = neq-4;
          }

        if (inv)
          { bost += TMER;
            n = keq+bost;
            x = sarr + swide * buck[(Comp[*n]<<2)|(3-(n[-4]&0x3))]++;
            *x++ = 0;
            for (i = 4; i < KMER; i += 4)
              *x++ = Comp[n[-i]];
            *x++ = 0;
            for (i = 0; i < PostBytes; i++)
              *x++ = bust[i];
            for (i = 0; i < ContBytes; i++)
              *x++ = nust[i];
          }
        else
          { n = neq+bost;
            x = sarr + swide * buck[(*n<<2)|(n[4]>>6)]++;
            *x++ = 0;
            for (i = 4; i < KMER; i += 4)
              *x++ = n[i];
            *x++ = 0;
            for (i = 0; i < PostBytes; i++)
              *x++ = bust[i];
            for (i = 0; i < ContBytes; i++)
              *x++ = cust[i];
          }
      }

      if (b >= btop)           //  refill input buffer if needed
        { int ex = bend-b;
          memcpy(buffer,b,ex); 
          bend = buffer+ex;
          iamt = read(in,bend,BUFF_MAX-ex);
          if (iamt < 0)
            { fprintf(stderr,"%s: IO read to file %s%d.idx failed\n",
                             Prog_Name,POST_NAME,parm->inum);
              exit (1);
            }
          bend += iamt;
          if (bend == buffer)
            break;
          if (bend-buffer < BUFF_MAX)
            btop = bend;
          else
            btop = bend-4;
          b = buffer;
        }
    }

  close(in);

  return (NULL);
}

static void *setup_thread_with_masks(void *args)
{ SP *parm = (SP *) args;
  int in        = parm->in;
  int tid       = parm->tid;
  int swide     = parm->swide;
  char  *SEQ    = parm->seq+1;
  uint8 *sarr   = parm->sarr;
  uint8 *buffer = parm->buff;
  GDB   *gdb    = &(parm->gdb);
  int64 *buck   = Buckets[tid];

  int    iamt;
  int    len;
  uint8 *neq, *keq;
  int    ncntg, inv;
  int64  post, bost, last, lost;
  int64  cont, nont, flag;
  uint8 *lust = (uint8 *) (&last);
  uint8 *bust = (uint8 *) (&bost);
  uint8 *cust = (uint8 *) (&cont);
  uint8 *nust = (uint8 *) (&nont);
  int64  nextpost, basepost;
  uint8 *bend, *btop, *b;

  GDB_CONTIG *contigs = gdb->contigs;
  GDB_MASK   *masks   = gdb->masks, tempm;
  int         nextm, lastm, pbg;

  flag = (0x1ll << (8*ContBytes-1));

  if (lseek(in,0,SEEK_SET) < 0)
    { fprintf(stderr,"%s: Rewind of file %s%d.idx failed\n",
                     Prog_Name,POST_NAME,parm->inum);
      exit (1);
    }

  iamt = read(in,buffer,BUFF_MAX);
  if (iamt < 0)
    { fprintf(stderr,"%s: IO read to file %s%d.idx failed\n",
                     Prog_Name,POST_NAME,parm->inum);
      exit (1);
    }
  bend = buffer + iamt;
  if (bend-buffer < BUFF_MAX)
    btop = bend;
  else
    btop = bend-4;
  b = buffer;
  if (btop <= buffer)   //  Nothing to do
    { close(in);
      return (NULL);
    }

  lastm = -1;
  ncntg = DBsplit[tid];
  post  = DBpost[tid];
  nextpost = post;
  basepost = post;
  while (1)
    { { int   v;
        int64 extra;

        v = *b++;                    //  Get next encoded post, uncompress it and "un-delta"
        inv = ((v & 0x80) != 0);
        if (v & 0x40)
          if (v & 0x20)
            { extra = 0;
              while (v & 0x1)
                { extra -= 0x10000000;
                  v = *b++;
                }
              last = 0;
              lust[3] = (v & 0x0f); 
              lust[2] = *b++;
              lust[1] = *b++;
              lust[0] = *b++;
              last += extra;
            }
          else
            { last = 0;
              lust[1] = (v & 0x1f); 
              lust[0] = *b++;
            }
        else
          last = (v & 0x3f);
        post += last;
      }

      if (post >= nextpost)    //  if position is not in current contig, then fetch the next one
        {
          if (lastm >= 0)
            masks[lastm] = tempm;

          do
            { len = contigs[ncntg].clen;
              basepost = nextpost;
              if (contigs[ncntg].boff >= 0)
                nextpost += len;
              ncntg += 1;
            }
          while (post >= nextpost);

	  nextm = contigs[ncntg-1].moff;
          lastm = contigs[ncntg].moff;
          tempm = masks[lastm];
          masks[lastm].beg = masks[lastm].end = contigs[ncntg-1].clen+1;

          cont = InvP[ncntg-1];
          nont = cont | flag;
          lost = 0;
        }

      { int    i;        //  load the k-mer / post pair into the next available slot in the
        uint8 *n, *x;    //    initial sort array using the bucket index

        bost = post-basepost;

        while (bost >= masks[nextm].end)
          nextm += 1;
        if (bost < masks[nextm].beg)
          pbg = 0;
        else 
          { pbg = masks[nextm].end - bost;
            if (pbg > KMER)
              pbg = KMER;
          }

        if (bost >= lost)   //  Get another buffer full of the current contig
          { int64  top, bot;
            int    n;

            lost = bost + (SCAN_MAX - (2*KMER-8));
            if (lost > len)
              lost = len;
            bot = bost-(KMER-8);
            top = lost+KMER;
            if (top > len)
              top = len;
            if (bot < 0)
              bot = 0;
            neq = (uint8 *) Get_Contig_Piece(gdb,ncntg-1,bot,top,NUMERIC,SEQ) - bot;

            keq = neq+3;
            top -= 3;
            n = (neq[bot] << 4) | (neq[bot+1] << 2) | neq[bot+2];
            for (i = bot; i < top; i++)
              neq[i] = n = (((n << 2) | keq[i]) & 0xff);
            keq = neq-4;
          }

        if (inv)
          { bost += TMER;
            n = keq+bost;
            x = sarr + swide * buck[(Comp[*n]<<2)|(3-(n[-4]&0x3))]++;
            *x++ = 0;
            for (i = 4; i < KMER; i += 4)
              *x++ = Comp[n[-i]];
            *x++ = pbg;
            for (i = 0; i < PostBytes; i++)
              *x++ = bust[i];
            for (i = 0; i < ContBytes; i++)
              *x++ = nust[i];
          }
        else
          { n = neq+bost;
            x = sarr + swide * buck[(*n<<2)|(n[4]>>6)]++;
            *x++ = 0;
            for (i = 4; i < KMER; i += 4)
              *x++ = n[i];
            *x++ = pbg;
            for (i = 0; i < PostBytes; i++)
              *x++ = bust[i];
            for (i = 0; i < ContBytes; i++)
              *x++ = cust[i];
          }
      }

      if (b >= btop)           //  refill input buffer if needed
        { int ex = bend-b;
          memcpy(buffer,b,ex); 
          bend = buffer+ex;
          iamt = read(in,bend,BUFF_MAX-ex);
          if (iamt < 0)
            { fprintf(stderr,"%s: IO read to file %s%d.idx failed\n",
                             Prog_Name,POST_NAME,parm->inum);
              exit (1);
            }
          bend += iamt;
          if (bend == buffer)
            break;
          if (bend-buffer < BUFF_MAX)
            btop = bend;
          else
            btop = bend-4;
          b = buffer;
        }
    }

  close(in);

  return (NULL);
}

typedef struct
  { int      tid;
    int      swide;
    uint8   *sarr;
    int64   *prefix;
    int      part;
    int64    off;
    int64    span;
  } RP;

static int64 NKMER[NUM_BUCK];

static pthread_mutex_t TMUTEX;
static pthread_cond_t  TCOND;

//  Tstack[0..Tavail-1] is a stack of available threads at any moment.
//  It is always manipulated inside the mutex TMUTEX

static int *Tstack;
static int  Tavail;

//  for 1st k-mer byte range [beg,end), find each group of equal k-mers, then
//    overwrite to the bottom of the range.  K-mer payloads contain the prefix mask, count,
//    and then the contig/position location.

static void *compress_thread(void *args)
{ RP     *parm   = (RP *) args;
  int     swide  = parm->swide;
  uint8  *sarray = parm->sarr;
  int64  *prefix = parm->prefix;
  int     beg    = parm->part;
  int64   off    = parm->off;
  int64   span   = parm->span;

  int64   nkmer;
  int64   x, e, y;
  int     w, k, z;
  int     lcp, idx, flc;
  uint8  *b;

  nkmer = 0;
  x = off;
  b = (sarray + x) + 1;
  lcp = sarray[x];
  prefix += ((beg>>2) << 16);
  e = x + span;
  z = (sarray[e] == 0);   //  Caution: end of panel can have 0 lcp
  if (z)
    sarray[e] = 12;
  while (x < e)
    { y = x+swide;
      w = 1;
      while (sarray[y] == 0)
        { y += swide;
          w += 1;
        }

      //  sorted w entries in [x,y) are all equal
      //  output to k-mer table

      idx = ((sarray[x+1] << 8) | sarray[x+2]);
      prefix[idx] += w;
      nkmer += w;
      flc = lcp;
      while (x < y)
        { for (k = 3; k < MBYTES; k++)
            *b++ = sarray[x+k];
          *b++ = flc;
          flc  = 40;
          for (k = MBYTES; k < swide; k++)
            *b++ = sarray[x+k];
          x += swide;
        }

      lcp = sarray[x];
    }
  if (z)
    lcp = sarray[e] = 0;

  NKMER[beg] = nkmer;

#ifndef DEBUG_THREADS

  pthread_mutex_lock(&TMUTEX);   //  Put this thread back on the avail stack
    Tstack[Tavail++] = parm->tid;
  pthread_mutex_unlock(&TMUTEX);

  pthread_cond_signal(&TCOND);   //  Signal a thread is available

#endif

  return (NULL);
}

  //  Writes over 2GB don't work on some systems, patch to overcome said

static inline int64 big_write(int f, uint8 *buffer, int64 bytes)
{ int64 v, x;

  v = 0;
  while (bytes > 0x70000000)
    { x = write(f,buffer,0x70000000);
      if (x < 0)
        return (-1);
      v += x;
      bytes  -= x;
      buffer += x;
    }
  x = write(f,buffer,bytes);
  if (x < 0)
    return (-1);
  return (v+x);
}

void k_sort(GDB *gdb)
{ uint8 *sarray;
  uint8 *buffer;
  char  *seqbuf;
  int    p, part, swide;
  int64  panel[NUM_BUCK];
  int    tstack[NTHREADS];
  SP     sarm[NTHREADS];
  RP     rarm[NTHREADS];
  int    kbeg, kend;
  int64  nkmer;
  int64 *prefix;
#ifndef DEBUG_THREADS
  pthread_t threads[NTHREADS];
#endif

  //  Complete the bucket array by accumulating counts across k-mer 1st bytes

  { int   i, j; 
    int64 x, cum, nelmax;

    nelmax = 0;
    cum = 0;
    for (i = 0; i < NUM_BUCK; i++)
      { for (j = 0; j < NTHREADS; j++)
          { x = Buckets[j][i];
            Buckets[j][i] = cum;
            cum += x;
          }
        if (i >= NUM_BUCK-1 || Select[i] != Select[i+1])
          { if (cum > nelmax)
              nelmax = cum;
            cum = 0;
          }
      }

#ifdef SHOW_FINGERS
    for (i = 0; i < NUM_BUCK; i++)
      { printf(" %3d:",i);
        for (j = 0; j < NTHREADS; j++)
          printf(" %9lld",Buckets[j][i]);
        printf("  :: %d\n",Select[i]);
      }
#endif

    prefix = (int64 *) MyBlock;
    seqbuf = (char *) (prefix + 0x1000000);
    buffer = (uint8 *) (seqbuf + (SCAN_MAX+8)*NTHREADS);
    bzero(prefix,sizeof(int64)*0x1000000);

    swide  = MBYTES + PostBytes + ContBytes;
    sarray = Malloc(nelmax*swide+1,"Allocating sort array");
  }

  //  Must open a separate file descriptor to DB bases for each setup thread !!!

  for (p = 0; p < NTHREADS; p++)
    { sarm[p].tid   = p;
      sarm[p].swide = swide;
      sarm[p].sarr  = sarray;
      sarm[p].seq   = seqbuf + p*(SCAN_MAX+8);
      sarm[p].buff  = buffer + p*BUFF_MAX;
      sarm[p].gdb   = *gdb;
      if (p > 0)
        { sarm[p].gdb.seqs = fopen(gdb->seqpath,"r");
          if (sarm[p].gdb.seqs == NULL)
            { fprintf(stderr,"%s: Cannot open another copy of GDB\n",Prog_Name);
              exit (1);
            }
        }
    }

  for (p = 0; p < NTHREADS; p++)
    { rarm[p].tid    = p;
      rarm[p].swide  = swide;
      rarm[p].sarr   = sarray;
      rarm[p].prefix = prefix;
    }

  nkmer = 0;

  Tstack = tstack;
  for (p = 0; p < NTHREADS; p++)
    Tstack[p] = p;
  Tavail = NTHREADS;

  pthread_mutex_init(&TMUTEX,NULL);
  pthread_cond_init(&TCOND,NULL);

  for (part = 0; part < NPARTS; part++)

    { kbeg = Ksplit[part];
      kend = Ksplit[part+1];

      for (p = 0; p < NTHREADS; p++)
        { sarm[p].in = Units[part*NTHREADS+p];
          sarm[p].inum = part*NTHREADS+p;
        }

#ifdef DEBUG_THREADS
      for (p = 0; p < NTHREADS; p++)
        { if (MASK)
            setup_thread_with_masks(sarm+p);
          else
            setup_thread_plain(sarm+p);
        }
#else
      for (p = 1; p < NTHREADS; p++)
        { if (MASK)
            pthread_create(threads+p,NULL,setup_thread_with_masks,sarm+p);
          else
            pthread_create(threads+p,NULL,setup_thread_plain,sarm+p);
        }
      if (MASK)
        setup_thread_with_masks(sarm);
      else
        setup_thread_plain(sarm);
      for (p = 1; p < NTHREADS; p++)
        pthread_join(threads[p],NULL);
#endif

      { int64 prev, next;

        prev = 0;
        for (p = kbeg; p < kend; p++)
          { next = Buckets[NTHREADS-1][p];
            panel[p] = (next - prev)*swide;
            prev = next;
          }

        if (VERBOSE)
          { fprintf(stderr,"\r    Sorting part %d  ",part+1);
            fflush(stderr);
          }

        if (prev > 0)
          { msd_sort(sarray,next,swide,KBYTES,panel,kbeg,kend,NTHREADS);

#ifdef DEBUG_SORT
            print_table(sarray,swide,next);
#endif
          }
      }

      if (VERBOSE)
        { fprintf(stderr,"\r    Compressing part %d",part+1);
          fflush(stderr);
        }

      { int64 off;
        int   tid;

        off = 0;
        for (p = kbeg; p < kend; p++)
          {
#ifdef DEBUG_THREADS
            tid = 0;
#else
            pthread_mutex_lock(&TMUTEX);

            if (Tavail <= 0)                       //  all threads are busy, wait
              pthread_cond_wait(&TCOND,&TMUTEX);

            tid = Tstack[--Tavail];                //  thread tid is available

            pthread_mutex_unlock(&TMUTEX);
#endif

            // Launching job on thread tid
        
            rarm[tid].part = p;
            rarm[tid].off  = off;
            rarm[tid].span = panel[p];
            off += panel[p];

#ifdef DEBUG_THREADS
            compress_thread(rarm);
#else
            pthread_create(threads+tid,NULL,compress_thread,rarm+tid);
#endif
          }

#ifndef DEBUG_THREADS
        pthread_mutex_lock(&TMUTEX);   //  Wait for all the jobs to complete
        while (Tavail < NTHREADS)
          pthread_cond_wait(&TCOND,&TMUTEX);
        pthread_mutex_unlock(&TMUTEX);
#endif
      }

      { int   tout;
        int64 nents, xp, off;

        nents = 0;
        for (p = kbeg; p < kend; p++)
          nents += NKMER[p];
        nkmer += nents;

        if (VERBOSE)
          { fprintf(stderr,"\r    Outputing part %d  ",part+1);
            fflush(stderr);
          }

        tout = open(Catenate(TPATH,"/.",TROOT,Numbered_Suffix(".ktab.",part+1,"")),
                            O_WRONLY|O_CREAT|O_TRUNC,0666);
        if (tout < 0)
          { fprintf(stderr,"%s: Cannot open part file %s/.%s.ktab.%d for writing\n",
                           Prog_Name,TPATH,TROOT,p+1);
            goto remove_parts;
          }
        if (write(tout,&KMER,sizeof(int)) < 0) goto part_error;
        if (write(tout,&nents,sizeof(int64)) < 0) goto part_error;

        off = 0;
        for (p = kbeg; p < kend; p++)
          { xp = big_write(tout,sarray+off+1,NKMER[p]*(swide-2));
            if (xp != NKMER[p]*(swide-2))
              goto part_error;
            off += panel[p];
          }
        close(tout);
      }
    }

  for (p = 1; p < NTHREADS; p++)
    fclose(sarm[p].gdb.seqs);

  if (VERBOSE)
    { int64 ktot = gdb->seqtot - (KMER-1)*gdb->ncontig;

      fprintf(stderr,"\r    Done                                           \n");
      fprintf(stderr,"\n  Sampled:   %11lld (%5.1f%%) kmers/positions\n",
                     nkmer,(100.*nkmer)/ktot);
      fflush(stderr);
    }
  if (LOG_FILE)
    { int64 ktot = gdb->seqtot - (KMER-1)*gdb->ncontig;

      fprintf(LOG_FILE,"\n  Sampled:   %11lld (%5.1f%%) kmers/positions\n",
                       nkmer,(100.*nkmer)/ktot);
    }

  { int   tab;
    int   x, freq;
    int64 y;
    int64 maxpre;

    tab = open(Catenate(TPATH,"/",TROOT,".gix"),O_WRONLY|O_CREAT|O_TRUNC,0666);
    if (tab < 0)
      { fprintf(stderr,"%s: Cannot open %s/%s.gix for writing\n",Prog_Name,TPATH,TROOT);
        goto remove_parts;
      }
    if (write(tab,&KMER,sizeof(int)) < 0) goto gix_error;
    if (write(tab,&NPARTS,sizeof(int)) < 0) goto gix_error;
    x = 1;
    if (write(tab,&x,sizeof(int)) < 0) goto gix_error;
    x = 3;
    if (write(tab,&x,sizeof(int)) < 0) goto gix_error;

    maxpre = prefix[0];
    for (x = 1; x < 0x1000000; x++)
      { if (prefix[x] > maxpre)
          maxpre = prefix[x];
        prefix[x] += prefix[x-1];
      }
    if (write(tab,prefix,sizeof(int64)*0x1000000) < 0) goto gix_error;

    freq = 0;
    if (write(tab,&PostBytes,sizeof(int)) < 0) goto gix_error;
    if (write(tab,&ContBytes,sizeof(int)) < 0) goto gix_error;
    if (write(tab,&NPARTS,sizeof(int)) < 0) goto gix_error;
    if (write(tab,&maxpre,sizeof(int64)) < 0) goto gix_error;
    if (write(tab,&freq,sizeof(int)) < 0) goto gix_error;
    if (write(tab,&(gdb->ncontig),sizeof(int)) < 0) goto gix_error;
    if (write(tab,Perm,sizeof(int)*gdb->ncontig) < 0) goto gix_error;

    y = -1;
    if (write(tab,&y,sizeof(int64)) < 0) goto gix_error;

    close(tab);
  }
 
  free(sarray);
  return;

gix_error:
  unlink(Catenate(TPATH,"/",TROOT,".gix"));
  fprintf(stderr,"%s: IO error while writing %s/%s.gix\n",Prog_Name,TPATH,TROOT);
  goto remove_parts;

part_error:
  fprintf(stderr,"%s: IO error writing to part file %s/.%s.ktab.%d\n",Prog_Name,TPATH,TROOT,part+1);
remove_parts:
  for (p = 1; p <= NTHREADS; p++)
    unlink(Catenate(TPATH,"/.",TROOT,Numbered_Suffix(".ktab.",p,"")));
  exit (1);
}


/***********************************************************************************************
 *
 *   MAIN
 *
 **********************************************************************************************/

static void short_GDB_fix(GDB *gdb)
{ int i;

  if (gdb->ncontig >= NTHREADS)
    return;

  //  Add additional conitgs of length KMER that are the first bit of the 0th read/contig
  //    Mark as "fake" with -1 in origin field.

  gdb->contigs = Realloc(gdb->contigs,sizeof(GDB_CONTIG)*NTHREADS,"Reallocating DB read vector");
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

static GDB_CONTIG *CONTIGS;

static int LSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  return (CONTIGS[y].clen - CONTIGS[x].clen);
}

int main(int argc, char *argv[])
{ GDB    _gdb, *gdb = &_gdb;
  int     ftype;
  char   *spath, *tpath;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("GIXmake");

    KMER = 40;
    LOG_PATH = NULL;
    LOG_FILE = NULL;
    NTHREADS = 8;
    SORT_PATH = getenv("TMPDIR");
    if (SORT_PATH == NULL)
      SORT_PATH = ".";

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'k':
            ARG_NON_NEGATIVE(KMER,"index k-mer size");
            break;
          case 'L':
            if (argv[i][2] != ':')
              { fprintf (stderr,"%s: option -L must be followed by :<filename>\n",Prog_Name);
                exit (1);
              }
            LOG_PATH = argv[i]+3;
            LOG_FILE = fopen (LOG_PATH,"a");
            if (LOG_FILE == NULL)
              { fprintf (stderr,"%s: Cannot open logfile %s for output\n",Prog_Name,LOG_PATH);
                exit(1);
              }
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

    VERBOSE = flags['v'];

    KBYTES  = (KMER>>2);
    if (argc < 2 || argc > 3)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -L: Output log to specified file.\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        fprintf(stderr,"      -P: Directory to use for temporary files.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -k: index k-mer size\n");
        exit (1);
      }

    if ((KMER & 0x3) != 0)
      { fprintf(stderr,"%s: K-mer size must be a multiple of 4\n",Prog_Name);
        exit (1);
      }
    if (NTHREADS > 32)
      { fprintf(stderr,"%s: # of threads can be at most 32, more doesn't help.\n",Prog_Name);
        exit (1);
      }
  }

  //  Determine source and target root names, paths, and extensions

  if (argc == 3)
    ftype = Get_GDB_Paths(argv[1],argv[2],&spath,&tpath,0);
  else
    ftype = Get_GDB_Paths(argv[1],NULL,&spath,&tpath,0);
  TPATH = PathTo(tpath);
  TROOT = Root(tpath,NULL);

  if (ftype == IS_GDB && argc == 3)
    { fprintf(stderr,
              "%s: Cannot create a .gix with a different location and root name than its .gdb\n",
              Prog_Name);
      exit (1);
    }

  if (ftype != IS_GDB)
    { char *command;

      command = Malloc(strlen(spath)+strlen(tpath)+100,"Allocating command string");
      if (LOG_FILE)
        { fclose(LOG_FILE);
          sprintf(command,"FAtoGDB%s -L:%s %s %s",VERBOSE?" -v":"",LOG_PATH,spath,tpath);
        }
      else
        sprintf(command,"FAtoGDB%s %s %s",VERBOSE?" -v":"",spath,tpath);
      if (system(command) != 0)
        { fprintf(stderr,"\n%s: Call to FAtoGDB failed\n",Prog_Name);
          exit (1);
        }
      if (LOG_PATH)
        LOG_FILE = fopen(LOG_PATH,"a");
      free(command);
    }

  if (VERBOSE || LOG_FILE)
    StartTime();

  if (LOG_FILE)
    fprintf(LOG_FILE,"\n%s\n", Command_Line);

  if (VERBOSE)
    { if (strcmp(TPATH,".") == 0)
        fprintf(stderr,"\n  Creating genome index (GIX) %s.gix in the current directory\n",TROOT);
      else
        fprintf(stderr,"\n  Creating genome genome index (GIX) %s.gix in directory %s\n",
                       TROOT,TPATH);
      fflush(stderr);
    }

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
                exit (1);
              }
          }
        else
          spath = Catenate(cpath,"/",SORT_PATH,"");
        SORT_PATH = Strdup(spath,"Allocating temporary sort path");
        free(cpath);
      }
    else
      SORT_PATH = Strdup(SORT_PATH,"Allocating temporary sort path");

    if ((dirp = opendir(SORT_PATH)) == NULL)
      { fprintf(stderr,"\n%s: -P option: cannot open directory %s\n",Prog_Name,SORT_PATH);
        exit (1);
      }
    closedir(dirp);

    if (SORT_PATH[strlen(SORT_PATH)-1] == '/')
      SORT_PATH[strlen(SORT_PATH)-1] = '\0';
  }

  //  Open GDB

  POST_NAME = Strdup(Catenate(SORT_PATH,"/.",Numbered_Suffix("post.",getpid(),"."),""),
                     "Allocating post index name");

  Read_GDB(gdb,tpath);
  short_GDB_fix(gdb);

  MASK   = (gdb->nmasks > 0);
  MBYTES = KBYTES+1;

  { int i, l0, l1, l2, l3;   //  Compute byte complement table

    i = 0;
    for (l0 = 3; l0 >= 0; l0 -= 1)
     for (l1 = 12; l1 >= 0; l1 -= 4)
      for (l2 = 48; l2 >= 0; l2 -= 16)
       for (l3 = 192; l3 >= 0; l3 -= 64)
         Comp[i++] = (l3 | l2 | l1 | l0);
  }

  //  Compute NTHREADS GDB partition

  { int64 npost, range, cum, t;   //  Compute DB split into NTHREADS parts
    int   p, r, len;

    DBsplit = Malloc((NTHREADS+1)*sizeof(int),"Allocating DB split arrays");
    DBpost  = Malloc((NTHREADS+1)*sizeof(int64),"Allocating DB split arrays");

    npost = gdb->seqtot;
    cum   = 0;
    range = 0;

    DBsplit[0] = 0;
    DBpost [0] = 0;
    p = 1;
    t = (npost*p)/NTHREADS;
    for (r = 0; r < gdb->ncontig; r++)
      { len = gdb->contigs[r].clen;
        cum += len;
        while (cum >= t)
          { DBsplit[p] = r+1;
            DBpost [p] = cum;
            p += 1;
            t = (npost*p)/NTHREADS;
          }
        if (range < len)
          range = len;
      }
    DBsplit[NTHREADS] = gdb->ncontig;
    DBpost [NTHREADS] = npost;

    PostBytes = 0;                 //  # of bytes for encoding a post
    cum = 1;
    while (cum < range)
      { cum *= 256;
        PostBytes += 1;
      }

    range = 2*gdb->ncontig;
    ContBytes = 0;                 //  # of bytes for encoding a contig + sign bit
    cum = 1;
    while (cum < range)
      { cum *= 256;
        ContBytes += 1;
      }
  }

  //  Determine how many parts for 4GB sorts and allocate 1st byte split array
  //    (to be determined during training in the first phase

  { int64 nels;
    int   nbit;

    nels = 0x100000000ll / (ContBytes + PostBytes + KBYTES + 2);
    nbit = (.81 * (gdb->seqtot - (KMER-1)*gdb->ncontig)) / nels;

    NPARTS = ((nbit-1)/NTHREADS+1)*NTHREADS;
    if (NPARTS < 8)
      NPARTS = 8;
    else if (NPARTS > 64)
      NPARTS = 64;

    Ksplit = Malloc((NPARTS+1)*sizeof(int),"Allocating split vector");
  }

  //  Make sure you can open (NPARTS + 2) * NTHREADS + 1 + tid files at one time.
  //    tid is typically 3 unless using valgrind or other instrumentation.

  { struct rlimit rlp;
    int           tid;
    uint64        nfiles;

    tid = open(".xxx",O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
    close(tid);
    unlink(".xxx");

    nfiles = (NPARTS+2)*NTHREADS + 3 + tid + 100;    // RD 250221 add 100 to allow for a
                                                         //   few wrapping processes
    getrlimit(RLIMIT_NOFILE,&rlp);
    if (nfiles > rlp.rlim_max)
      { fprintf(stderr,"\n%s: Cannot open %lld files simultaneously\n",Prog_Name,nfiles);
        exit (1);
      } 
    rlp.rlim_cur = nfiles;
    if (setrlimit(RLIMIT_NOFILE,&rlp) < 0)
      { fprintf(stderr,"%s: Could not increase IO unit resourc to %d\n",
                       Prog_Name,(2*NTHREADS+4)*NTHREADS);
        exit (1);
      }
  } 

  //  Produce perms for length sorted order of contigs

  { int i;
 
    Perm = Malloc(2*gdb->ncontig*sizeof(int),"Allocating sort permutation arrays");
    InvP = Perm + gdb->ncontig;

    for (i = 0; i < gdb->ncontig; i++)
      Perm[i] = i;
  
    CONTIGS = gdb->contigs;
    qsort(Perm,gdb->ncontig,sizeof(int),LSORT);

    for (i = 0; i < gdb->ncontig; i++)
      InvP[Perm[i]] = i;
  }

  if (VERBOSE)
    { fprintf(stderr,"\n  Partitioning K-mers via pos-lists into %d parts\n",NPARTS);
      fflush(stderr);
    }
  if (LOG_FILE)
    fprintf(LOG_FILE,"\n  Partitioning K-mers via pos-lists into %d parts\n",NPARTS);

  { int p;   //  Setup distribution bucket array

    Buckets = Malloc(NTHREADS*sizeof(int64 *),"Allocating distribution buckets");
    Buckets[0] = Malloc(NTHREADS*NUM_BUCK*sizeof(int64),"Allocating distribution buckets");
    for (p = 1; p < NTHREADS; p++)
      Buckets[p] = Buckets[p-1] + NUM_BUCK;
  }

  { int   p, i, k;         //  Open IO units for distribution and reimport
    char *name;

    Units = Malloc(NTHREADS*NPARTS*sizeof(int),"Allocating IO Units");

    k = 0;
    for (p = 0; p < NTHREADS; p++)
      for (i = 0; i < NPARTS; i++)
        { name = Numbered_Suffix(POST_NAME,k,".idx");
          Units[k] = open(name,O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
          if (Units[k] < 0)
            { fprintf(stderr,"%s: Cannot open %s for reading & writing\n",Prog_Name,name);
              exit (1);
            }
          unlink(name);
          k += 1;
        }
  }

  { int64 p1, p2;

    p1 = ((SCAN_MAX+8) + NPARTS*BUFF_MAX)*NTHREADS;
    p2 = ((SCAN_MAX+8) + BUFF_MAX)*NTHREADS + sizeof(int64)*0x1000000;
    if (p1 > p2)
      MyBlock = Malloc(p1,"Allocating memory block");
    else
      MyBlock = Malloc(p2,"Allocating memory block");
  }
    
  distribute(gdb);   //  Distribute k-mers to 1st byte partitions, encoded as compressed
                     //    relative positions of the given k-mers
  if (VERBOSE)
    { TimeTo(stderr,0,LOG_FILE==NULL);
      fprintf(stderr,"\n  Starting sort & index output of each part\n");
      fflush(stderr);
    }
  if (LOG_FILE)
    { TimeTo(LOG_FILE,0,1);
      fprintf(LOG_FILE,"\n  Sorting & output of each part\n");
    }

  k_sort(gdb);  //  Reimport the post listings, recreating the k-mers and sorting
                //    them with their posts to produce the final genome index.

  free(MyBlock);
  free(Units);
  free(Buckets[0]);
  free(Buckets);
  free(Perm);
  free(Ksplit);
  free(DBpost);
  free(DBsplit);

  Close_GDB(gdb);

  free(POST_NAME);

  free(TROOT);
  free(TPATH);

  free(tpath);
  free(spath);

  if (VERBOSE)
    { TimeTo(stderr,0,LOG_FILE==NULL);
      TimeTo(stderr,1,0);
    }
  if (LOG_FILE)
    { TimeTo(LOG_FILE,0,1);
      TimeTo(LOG_FILE,1,0);
      fclose(LOG_FILE);
    }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
