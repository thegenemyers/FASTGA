/*  Last edited: Feb 21 16:55 2025 (rd109) */
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

#undef  DEBUG_MAP
#undef  DEBUG_THREADS
#undef  DEBUG_SORT

  //  These constants determine the syncmer (TMER,SMER) and whether to add a spaced seed

#define TMER 12        //  Syncmer len
#define SMER  8        //  S-mer len (4-8 only)
#define SRES  4        //  SMER-4
#define SOFF  4        //  TMER-SMER

#undef  SPACER         //  Single spaced-seed?

static char *Usage[] =
    { "[-v] [-T<int(8)>] [-P<dir($TMPDIR)>] [-k<int(40)] [-f<int(10)>]",
      "( <source:path>[.1gdb]  |  <source:path>[<fa_extn>|<1_extn>] [<target:path>[.gix]] )"
    };

static int   FREQ;       //  -f
static int   VERBOSE;    //  -v
static char *SORT_PATH;  //  -P
static char *TPATH;
static char *TROOT;
static char *POST_NAME;

static int NTHREADS;   //  by default 8
static int KMER;       //  by default 40, must be >= 12 and divisible by 4
static int KBYTES;     //  Bytes for 2-bit compress k-mer (KMER/4)

#define BUFFER_LEN  1000000

static int   *Perm;        //  Size sorted permutation of contigs
static int   *InvP;        //  Inverse of Perm

static int    PostBytes;   //  # of bytes needed for a position
static int    ContBytes;   //  # of bytes for a contig + sign bit

static int    Comp[256];   //  DNA complement of packed byte
static int    Select[256]; //  1st k-mer byte -> block (of NTHREAD)

                           //  For p in [0,NTHREADS):
static int   *DBsplit;     //    DB split: contig [DBsplit[p],DBsplit[p+1])
static int64 *DBpost;      //    DB post: post of start of contig DBsplit[p] is DBpost[p]
static int   *Ksplit;      //    Kmer split:  1st bytes [Ksplit[p],Ksplit[p+1])

static int64 **Buckets;    //  NTHREADS x 256 array:
                           //  B[i][j] = # of posts whose canonical k-mer's 1st byte
                           //            is j from block DBsplit[i],DBsplit[i+1]
                           //  Also used as "fingers" for loading initial sort array.

static int *Units;  //  NTHREADS^2 IO units for distribution and import & k-mer table parts
static int *Pnits;  //  NTHREADS^2 extra IO units for post parts

typedef struct
  { int   beg;
    int   end;
    int64 off;
  } Range;

extern void msd_sort(uint8 *array, int64 nelem, int rsize, int ksize,
                     int64 *part, int nthreads, Range *range);


/***********************************************************************************************
 *
 *   DISTRIBUTION PHASE:  Output signed posts for each 1/NTHREAD of the data and the k-mers
 *        At completion file ._post.<jobid>.<i*NTHREADS+j>.idx contains the signed posts from
 *             contigs (DBsplit[j],DBplist[j+1]) whose k-mers 1st byte is in Ksplit[i],Ksplit[i+1].
 *
 **********************************************************************************************/

#ifdef GEN_HASH     //  Code to generate hash map

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{ int i1, i2, i3, i4;
  int i;
  int hash[256];

  for (i1 = 0; i1 < 4; i1++)
  for (i2 = 0; i2 < 4; i2++)
  for (i3 = 0; i3 < 4; i3++)
  for (i4 = 0; i4 < 4; i4++)
    hash[(i1<<6)|(i2<<4)|(i3<<2)|i4] = (i1<<6)|(((i2+1)&0x3)<<4)|(((i3+2)&0x3)<<2)|((i4+3)&0x3);

  printf("static int hash[256] =\n  { %02x",hash[0]);
  for (i = 1; i < 256; i++)
    if (i%8 == 0)
      printf(",\n    %02x",hash[i]);
    else
      printf(", %02x",hash[i]);
  printf("\n  };\n");

  exit (0);
}

#endif

static int hash[256] =
  { 0x1b, 0x18, 0x19, 0x1a, 0x1f, 0x1c, 0x1d, 0x1e, 0x13, 0x10, 0x11, 0x12, 0x17, 0x14, 0x15, 0x16,
    0x2b, 0x28, 0x29, 0x2a, 0x2f, 0x2c, 0x2d, 0x2e, 0x23, 0x20, 0x21, 0x22, 0x27, 0x24, 0x25, 0x26,
    0x3b, 0x38, 0x39, 0x3a, 0x3f, 0x3c, 0x3d, 0x3e, 0x33, 0x30, 0x31, 0x32, 0x37, 0x34, 0x35, 0x36,
    0x0b, 0x08, 0x09, 0x0a, 0x0f, 0x0c, 0x0d, 0x0e, 0x03, 0x00, 0x01, 0x02, 0x07, 0x04, 0x05, 0x06,
    0x5b, 0x58, 0x59, 0x5a, 0x5f, 0x5c, 0x5d, 0x5e, 0x53, 0x50, 0x51, 0x52, 0x57, 0x54, 0x55, 0x56,
    0x6b, 0x68, 0x69, 0x6a, 0x6f, 0x6c, 0x6d, 0x6e, 0x63, 0x60, 0x61, 0x62, 0x67, 0x64, 0x65, 0x66,
    0x7b, 0x78, 0x79, 0x7a, 0x7f, 0x7c, 0x7d, 0x7e, 0x73, 0x70, 0x71, 0x72, 0x77, 0x74, 0x75, 0x76,
    0x4b, 0x48, 0x49, 0x4a, 0x4f, 0x4c, 0x4d, 0x4e, 0x43, 0x40, 0x41, 0x42, 0x47, 0x44, 0x45, 0x46,
    0x9b, 0x98, 0x99, 0x9a, 0x9f, 0x9c, 0x9d, 0x9e, 0x93, 0x90, 0x91, 0x92, 0x97, 0x94, 0x95, 0x96,
    0xab, 0xa8, 0xa9, 0xaa, 0xaf, 0xac, 0xad, 0xae, 0xa3, 0xa0, 0xa1, 0xa2, 0xa7, 0xa4, 0xa5, 0xa6,
    0xbb, 0xb8, 0xb9, 0xba, 0xbf, 0xbc, 0xbd, 0xbe, 0xb3, 0xb0, 0xb1, 0xb2, 0xb7, 0xb4, 0xb5, 0xb6,
    0x8b, 0x88, 0x89, 0x8a, 0x8f, 0x8c, 0x8d, 0x8e, 0x83, 0x80, 0x81, 0x82, 0x87, 0x84, 0x85, 0x86,
    0xdb, 0xd8, 0xd9, 0xda, 0xdf, 0xdc, 0xdd, 0xde, 0xd3, 0xd0, 0xd1, 0xd2, 0xd7, 0xd4, 0xd5, 0xd6,
    0xeb, 0xe8, 0xe9, 0xea, 0xef, 0xec, 0xed, 0xee, 0xe3, 0xe0, 0xe1, 0xe2, 0xe7, 0xe4, 0xe5, 0xe6,
    0xfb, 0xf8, 0xf9, 0xfa, 0xff, 0xfc, 0xfd, 0xfe, 0xf3, 0xf0, 0xf1, 0xf2, 0xf7, 0xf4, 0xf5, 0xf6,
    0xcb, 0xc8, 0xc9, 0xca, 0xcf, 0xcc, 0xcd, 0xce, 0xc3, 0xc0, 0xc1, 0xc2, 0xc7, 0xc4, 0xc5, 0xc6
  };

#ifdef DEBUG_SYNCMERS

static int is_syncmer(int x)
{ int y, z, m, b;
  int q, s;
 
  b = 256;
  q = 16;
  for (s = 16; s >= 0; s -= 2)
    { y = (x>>s) & 0xff;
      z = Comp[y];
      if (hash[y] < hash[z])
        m = hash[y];
      else
        m = hash[z];
      if (m < b)
        { q = s;
          b = m;
        }
    }
  return (q == 16 || m <= b);
}

#endif

static void split_4bits(GDB *gdb, int64 *buck)
{ int KM4 = KMER-(4+SRES);

  int     counter;
  uint8  *seq;
  int     len;
  uint8  *neq, *ceq;
  uint8  *npl, *cpl;
  uint16 *mzr;
  int     r, i;

  seq = (uint8 *) New_Contig_Buffer(gdb);
  neq = (uint8 *) New_Contig_Buffer(gdb);
  ceq = (uint8 *) New_Contig_Buffer(gdb);
  mzr = (uint16 *) Malloc(sizeof(uint16)*(gdb->maxctg+4),"Allocating minimizer array") + 1;
  if (seq == NULL || neq == NULL || ceq == NULL || mzr == NULL)
    exit (1);
  npl = neq + SRES;
  cpl = ceq + SRES;

  counter = 0;
  for (r = 0; r < gdb->ncontig; r++)
    { if (gdb->contigs[r].boff < 0)
        continue;

      len = gdb->contigs[r].clen;
      if (len < TMER)
        continue;

      Get_Contig(gdb,r,NUMERIC,(char *) seq);   //  Load the contig

      //  compute neq, ceq, and mzr from seq for contig

      len -= 3;

      { uint8 *s1, *s2, *s3, *n1;
        int    n, c, k;

        s1 = seq+1;
        s2 = seq+2;
        s3 = seq+3;
        for (i = 0; i < len; i += 4)
          { neq[i] = (seq[i] << 6) | (s1[i] << 4) | (s2[i] << 2) | s3[i];
            ceq[i] = Comp[neq[i]];
          }
        n1 = neq-1;
        for (k = 1; k < 4; k++)
          for (i = k; i < len; i += 4)
            { neq[i] = (n1[i] << 2) | s3[i];
              ceq[i] = Comp[neq[i]];
            }

        len -= SRES;
  
        for (i = 0; i < len; i++)
#if SRES == 0
          { n = hash[neq[i]];
            c = hash[ceq[i]];
#else
          { n = hash[neq[i]] * 251 + hash[npl[i]];
            c = hash[ceq[i]] * 251 + hash[cpl[i]];
#endif
            if (n < c)
              mzr[i] = n;
            else
              mzr[i] = c;
          }
      }

      //  Compute processed seq array from neq,ceq
    
      { int    min4, pos4;
        int    w, j;

        min4 = mzr[0];
        pos4 = 0;
        for (i = 1; i < SOFF; i++)
          if (mzr[i] < min4)
            { min4 = mzr[i];
              pos4 = i;
            }
        for (i = SOFF; i < len; i++)
          { if (mzr[i] < min4)              //  right-end of 12-syncmer
              { min4 = mzr[i];
                pos4 = i;
              }
            else if (pos4 == i-SOFF)          //   left-end of 12-syncmer
              { min4 = mzr[++pos4];
                for (j = pos4+1; j <= i; j++)
                  if (mzr[j] < min4)
                    { min4 = mzr[j];
                      pos4 = j;
                    }
              }
            else if (mzr[i] > min4)
              continue;

            j = i-SOFF;
            if (j+KM4 < len)
              { w = neq[j];
                buck[w] += 1;
                counter += 1;
              }
            if (i >= KM4)
              { w = cpl[i];
                buck[w] += 1;
                counter += 1;
              }
            if (counter >= 1000000)
              goto enough;
          }
      }
    }
enough:

  { int64 t;
    int   r;

    for (i = 1; i < 256; i++)
      buck[i] = buck[i-1] + buck[i];

    Ksplit[0] = 0;
    r = 1;
    t = buck[255]/NTHREADS;
    for (i = 0; i < 256; i++)
      { if (buck[i] >= t)
          { if (buck[i]-t > t-buck[i-1])
              { Select[i] = r;
                Ksplit[r] = i;
              }
            else
              { Select[i] = r-1;
                Ksplit[r] = i+1;
              }
            r += 1;
            t = (r*buck[255])/NTHREADS;
          }
        else
          Select[i] = r-1;
      }
    Ksplit[NTHREADS] = 256;
  }
    
#ifdef EXAMINE_TRIAL
  printf("Prepatory %lld %d\n",buck[255],counter);

  for (i = 0; i <= NTHREADS; i++)
    printf(" %3d: %3d\n",i,Ksplit[i]);
  for (i = 0; i < 256; i++)
    printf(" %3d: %8lld  %d\n",i,buck[i],Select[i]);
  fflush(stdout);
#endif

  free(seq-1);
  free(neq-1);
  free(ceq-1);
  free(mzr-1);
}

typedef struct
  { int    tid;
    int64 *buck;
    GDB    gdb;
  } DP;

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
{ DP   *parm  = (DP *) args;
  int    tid  = parm->tid;
  int64 *buck = parm->buck;
  GDB   *gdb  = &parm->gdb;
  int    KM4  = KMER-(4+SRES);

  Packet packs[NTHREADS];

  uint8  *seq, *buffer;
  int     len, ren;
  uint8  *neq, *ceq;
  uint8  *npl, *cpl;
  uint16 *mzr;
  int     r, i;
  int64   post;

  buffer = malloc(NTHREADS*BUFFER_LEN);
  for (i = 0; i < NTHREADS; i++)
    { packs[i].buffer = buffer+i*BUFFER_LEN;
      packs[i].ptr    = packs[i].buffer;
      packs[i].end    = packs[i].buffer + (BUFFER_LEN-4);
      packs[i].out    = Units[i*NTHREADS+tid];
      packs[i].inum   = NTHREADS*tid+i;
      packs[i].last   = DBpost[tid];
    }

  seq = (uint8 *) New_Contig_Buffer(gdb);
  neq = (uint8 *) New_Contig_Buffer(gdb);
  ceq = (uint8 *) New_Contig_Buffer(gdb);
  mzr = (uint16 *) Malloc(sizeof(uint16)*(gdb->maxctg+4),"Allocating minimizer array") + 1;
  if (seq == NULL || neq == NULL || ceq == NULL || mzr == NULL)
    exit (1);
  npl = neq + SRES;
  cpl = ceq + SRES;

  post = DBpost[tid];
  ren  = DBsplit[tid+1];
  for (r = DBsplit[tid]; r < ren; r++)
    { if (gdb->contigs[r].boff < 0)
        continue;

      len = gdb->contigs[r].clen;
      if (len < TMER)
        { post += len;
          continue;
        }

      Get_Contig(gdb,r,NUMERIC,(char *) seq);   //  Load the contig
    
#ifdef DEBUG_MAP
      printf("Src:");
      for (i = 0; i < 50; i++)
        printf(" %2d",seq[i]);
      printf("\n");
      fflush(stdout);
#endif

      //  compute neq, ceq, and mzr from seq for contig

      len -= 3;

      { uint8 *s1, *s2, *s3, *n1;
        int    n, c, k;

        s1 = seq+1;
        s2 = seq+2;
        s3 = seq+3;
        for (i = 0; i < len; i += 4)
          { neq[i] = (seq[i] << 6) | (s1[i] << 4) | (s2[i] << 2) | s3[i];
            ceq[i] = Comp[neq[i]];
          }
        n1 = neq-1;
        for (k = 1; k < 4; k++)
          for (i = k; i < len; i += 4)
            { neq[i] = (n1[i] << 2) | s3[i];
              ceq[i] = Comp[neq[i]];
            }

        len -= SRES;
  
        for (i = 0; i < len; i++)
#if SRES == 0
          { n = hash[neq[i]];
            c = hash[ceq[i]];
#else
          { n = hash[neq[i]] * 251 + hash[npl[i]];
            c = hash[ceq[i]] * 251 + hash[cpl[i]];
#endif
            if (n < c)
              mzr[i] = n;
            else
              mzr[i] = c;
          }
      }
    
#ifdef DEBUG_MAP
      printf("Neq:");
      for (i = 0; i < 50; i++)
        printf(" %02x",neq[i]);
      printf("\n");

      printf("Ceq:");
      for (i = 0; i < 50; i++)
        printf(" %02x",ceq[i]);
      printf("\n");

      printf("Mzr:");
      for (i = 0; i < 50; i++)
        printf(" %02x",mzr[i]);
      printf("\n");

      printf("    ");
      for (i = 0; i < 50; i++)
        printf(" %2d",i);
      printf("\n");
      fflush(stdout);
#endif

      //  Compute processed seq array from neq,ceq

#ifdef LINEAR_TIME_CODE
        { int min4[16];
          int pos4[16];

          bot = 0;
          min4[bot] = mzr[0];
          pos4[bot] = 0;
          for (i = 1; i < len; i++)
            { while (mzr[i] <= min4[top])
                { if (top == bot)
                    { min4[top] = mzr[i];
                      pos4[top] = i;
                      goto hit;
                    }
                  top = ((top-1) & 0xf);
                }
              top = ((top+1) & 0xff);
              min4[top] = mzr[i];
              pos4[top] = i;
              if (i-8 == pos4[bot])
                bot = ((bot+1) & 0xf);
              else
                continue;
  
            hit:
              Hit at i-8
            }
        }
#endif

      { int    min4, pos4;
        int    w, j;

        min4 = mzr[0];
        pos4 = 0;
        for (i = 1; i < SOFF; i++)
          if (mzr[i] < min4)
            { min4 = mzr[i];
              pos4 = i;
            }
        for (i = SOFF; i < len; i++)
          {
#ifdef DEBUG_MAP
            printf("%3d: %02x %02x %d\n",i,mzr[i],min4,pos4);
#endif
            if (mzr[i] < min4)              //  right-end of 12-syncmer
              { min4 = mzr[i];
                pos4 = i;
#ifdef DEBUG_MAP
                printf("      Hit R");
#endif
              }
            else if (pos4 == i-SOFF)          //   left-end of 12-syncmer
              { min4 = mzr[++pos4];
                for (j = pos4+1; j <= i; j++)
                  if (mzr[j] < min4)
                    { min4 = mzr[j];
                      pos4 = j;
                    }
#ifdef DEBUG_MAP
                printf("      Hit L");
#endif
              }
            else if (mzr[i] > min4)
              continue;
#ifdef DEBUG_MAP
            else
              printf("      Hit RE");
#endif

	    j = i-SOFF;
            if (j+KM4 < len)
              { w = neq[j];
                buck[w] += 1;
#ifdef DEBUG_MAP
                printf(" %02x %d ::",w,Select[w]);
#endif
                push(packs+Select[w],post+j,0);
              }
            if (i >= KM4)
              { w = cpl[i];
                buck[w] += 1;
#ifdef DEBUG_MAP
                printf(" %02x %d",w,Select[w]);
#endif
                push(packs+Select[w],post+j,1);
              }
#ifdef DEBUG_MAP
             printf("\n");
#endif
          }
      }

      post += len + SRES + 3;
    }

  for (i = 0; i < NTHREADS; i++)
    if (packs[i].ptr > packs[i].buffer)
      write(packs[i].out,packs[i].buffer,packs[i].ptr-packs[i].buffer);

  free(buffer);
  free(seq-1);
  free(neq-1);
  free(ceq-1);
  free(mzr-1);

  return (NULL);
}

//  Create a thread for each DB section, and have it distribute syncmer posts to NTHREAD files
//     according to the 1st byte of its k-mer.

void distribute(GDB *gdb)
{ int i;
  DP  parm[NTHREADS];
#ifndef DEBUG_THREADS
  pthread_t threads[NTHREADS];
#endif

  Buckets = malloc(NTHREADS*sizeof(int64 *));
  Buckets[0] = malloc(sizeof(int64)*256*NTHREADS);

  for (i = 0; i < NTHREADS; i++)
    { parm[i].tid  = i;
      parm[i].buck = Buckets[i] = Buckets[0] + i*256;
      parm[i].gdb  = *gdb;
      if (i > 0)
        { parm[i].gdb.seqs = fopen(gdb->seqpath,"r");
          if (parm[i].gdb.seqs == NULL)
            { fprintf(stderr,"%s: Cannot open another copy of GDB\n",Prog_Name);
              exit (1);
            }
        }
    }

  bzero(Buckets[0],sizeof(int64)*256*NTHREADS);

  split_4bits(gdb,Buckets[0]);   //  Empirically sets up Ksplit,Select based on 1st 1Mbp

  bzero(Buckets[0],sizeof(int64)*256);

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
      for (k = 0; k < PostBytes; k++)
        pust[k] = array[x+KBYTES+k];
      for (k = 0; k < ContBytes; k++)
        cust[k] = array[x+KBYTES+PostBytes+k];
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
    uint8   *sarr;
    uint8   *buff;
  } SP;

//  Read post file, uncompressing it, recomputing the canonical k-mer at its absolute
//   position and loading the k-mer and the signed post into the initial soring array
//   for the current 1st byte panel.

static void *setup_thread(void *args)
{ SP *parm = (SP *) args;
  int in        = parm->in;
  int tid       = parm->tid;
  int swide     = parm->swide;
  uint8 *sarr   = parm->sarr;
  uint8 *buffer = parm->buff;
  GDB   *gdb    = &(parm->gdb);
  int64 *buck   = Buckets[tid];

  int    iamt;
  uint8 *seq;
  int    len;
  uint8 *neq, *ceq, *keq;
  int    ncntg, inv;
  int64  post, bost, last;
  int64  cont, nont, flag;
  uint8 *lust = (uint8 *) (&last);
  uint8 *bust = (uint8 *) (&bost);
  uint8 *cust = (uint8 *) (&cont);
  uint8 *nust = (uint8 *) (&nont);
  int64  nextpost, basepost;
  uint8 *bend, *btop, *b;

  seq = (uint8 *) New_Contig_Buffer(gdb);
  neq = (uint8 *) New_Contig_Buffer(gdb);
  ceq = (uint8 *) New_Contig_Buffer(gdb);
  if (seq == NULL || neq == NULL || ceq == NULL)
    exit (1);
  keq = ceq+8;

  flag = (0x1ll << (8*ContBytes-1));

  if (lseek(in,0,SEEK_SET) < 0)
    { fprintf(stderr,"%s: Rewind of file %s%d.idx failed\n",
                     Prog_Name,POST_NAME,parm->inum);
      exit (1);
    }

  iamt = read(in,buffer,BUFFER_LEN);
  if (iamt < 0)
    { fprintf(stderr,"%s: IO read to file %s%d.idx failed\n",
                     Prog_Name,POST_NAME,parm->inum);
      exit (1);
    }
  bend = buffer + iamt;
  if (bend-buffer < BUFFER_LEN)
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
        { int    i, k, l;
          uint8 *s1, *s2, *s3, *n1;

          do
            { len = gdb->contigs[ncntg].clen;
              basepost = nextpost;
              if (gdb->contigs[ncntg].boff >= 0)
                nextpost += len;
              ncntg += 1;
            }
          while (post >= nextpost);

          Get_Contig(gdb,ncntg-1,NUMERIC,(char *) seq);

          cont = InvP[ncntg-1];
          nont = cont | flag;

          s1 = seq+1;
          s2 = seq+2;
          s3 = seq+3;
          n1 = neq-1;
          l  = len-3;
          for (i = 0; i < l; i += 4)
            { neq[i] = (seq[i] << 6) | (s1[i] << 4) | (s2[i] << 2) | s3[i];
              ceq[i] = Comp[neq[i]];
            }
          for (k = 1; k < 4; k++)
            for (i = k; i < l; i += 4)
              { neq[i] = (n1[i] << 2) | s3[i];
                ceq[i] = Comp[neq[i]];
              }
        }

      { int    i;        //  load the k-mer / post pair into the next available slot in the
        uint8 *n, *x;    //    initial sort array using the bucket index

        bost = post-basepost;
        if (inv)
          { n = keq+bost;
            bost += 12;
            x = sarr + swide * buck[*n]++;
            *x++ = 0;
            for (i = 4; i < KMER; i += 4)
              *x++ = n[-i];
            for (i = 0; i < PostBytes; i++)
              *x++ = bust[i];
            for (i = 0; i < ContBytes; i++)
              *x++ = nust[i];

#ifdef SPACER
            buck[*n] += 1;
            *x++ = 0;
            *x++ = n[-4];
            for (i = 9; i < KMER; i += 4)
              *x++ = n[-i];
            for (i = 0; i < PostBytes; i++)
              *x++ = bust[i];
            for (i = 0; i < ContBytes; i++)
              *x++ = nust[i];
#endif
          }
        else
          { n = neq+bost;
            x = sarr + swide * buck[*n]++;
            *x++ = 0;
            for (i = 4; i < KMER; i += 4)
              *x++ = n[i];
            for (i = 0; i < PostBytes; i++)
              *x++ = bust[i];
            for (i = 0; i < ContBytes; i++)
              *x++ = cust[i];

#ifdef SPACER
            buck[*n] += 1;
            *x++ = 0;
            *x++ = n[4];
            for (i = 9; i < KMER; i += 4)
              *x++ = n[i];
            for (i = 0; i < PostBytes; i++)
              *x++ = bust[i];
            for (i = 0; i < ContBytes; i++)
              *x++ = cust[i];
#endif
          }
      }

      if (b >= btop)           //  refill input buffer if needed
        { int ex = bend-b;
          memcpy(buffer,b,ex); 
          bend = buffer+ex;
          iamt = read(in,bend,BUFFER_LEN-ex);
          if (iamt < 0)
            { fprintf(stderr,"%s: IO read to file %s%d.idx failed\n",
                             Prog_Name,POST_NAME,parm->inum);
              exit (1);
            }
          bend += iamt;
          if (bend == buffer)
            break;
          if (bend-buffer < BUFFER_LEN)
            btop = bend;
          else
            btop = bend-4;
          b = buffer;
        }
    }

  free(seq-1);
  free(neq-1);
  free(ceq-1);
  close(in);

  return (NULL);
}

typedef struct
  { int      inum;
    int      swide;
    int64   *panel;
    uint8   *sarr;
    uint8   *buff;
    Range   *range;
    int      tout;
    int      pout;
    int64   *prefix;
    int64   *posfix;
    int64    nelim;
    int64    nbase;
    int64    nkmer;
    int64    npost;
  } RP;

//  for 1st k-mer byte range [beg,end), find each group of equal k-mers and if < FREQ, then
//    write to FastK table under construtions and output (signed) posts of each equal
//    k-mer to the associated post list.  K-mer 2byte payloads contain the count in the
//    first byte and the lcp with its predecessor in the 2nd byte.

static void *output_thread(void *args)
{ RP     *parm  = (RP *) args;
  int     swide  = parm->swide;
  int64  *panel  = parm->panel;
  uint8  *sarray = parm->sarr;
  Range  *range  = parm->range;
  int     beg    = range->beg;
  int     end    = range->end;
  int     tout   = parm->tout;
  int     pout   = parm->pout;
  int64  *prefix = parm->prefix;
  int64  *posfix = parm->posfix;

  uint8 *buf1 = parm->buff;
  uint8 *buf2 = buf1 + 500000;
  uint8 *bed1 = (buf1 + 500000) - (KBYTES-1);
  uint8 *bed2 = (buf2 + 500000) - (PostBytes+ContBytes);

  int64   nelim, nkmer, nbase;
  int64   x, e, y;
  int     o, w, k, z, lcp, idx;
  uint8 *_w = (uint8 *) &w;
  uint8  *b, *c;

  nelim = nkmer = nbase = 0;
  o = KMER;
  b = buf1;
  c = buf2;
  x = range->off;
  lcp = sarray[x];
  prefix += (beg << 16);
  posfix += (beg << 8);
  for (o = beg; o < end; o++, prefix += 0x10000, posfix += 0x100)
   { e = x + panel[o];
     z = (sarray[e] == 0);   //  Caution: end of panel can have 0 lcp
     if (z)
       sarray[e] = 1;
     for (e = x + panel[o]; x < e; )
      { y = x+swide;
        w = 1;
        while (sarray[y] == 0)
          { y += swide;
            w += 1;
          }

        //  sorted w entries in [x,y) are all equal

        if (w >= FREQ)
          { x = y;
            if (sarray[x] < lcp)
              lcp = sarray[x];
            nelim += 1;
            nbase += w;
            continue;
          }

        //  if < FREQ then output to k-mer table and post list

        idx = sarray[x+1];
        posfix[idx] += (y-x)/swide;
        idx = ((idx << 8) | sarray[x+2]);
        prefix[idx] += 1;

        for (k = 3; k < KBYTES; k++)
          *b++ = sarray[x+k];
        *b++ = _w[0];
        *b++ = lcp;
        if (b >= bed1)
          { if (write(tout,buf1,b-buf1) < 0)
              { fprintf(stderr,"%s: IO write to file %s%d.ktb failed\n",
                               Prog_Name,POST_NAME,parm->inum);
                exit (1);
              }
            b = buf1;
          }
        nkmer += 1;

        while (x < y)
          { for (k = KBYTES; k < swide; k++)
              *c++ = sarray[x+k];
            if (c >= bed2)
              { if (write(pout,buf2,c-buf2) < 0)
                  { fprintf(stderr,"%s: IO write to file %s%d.pst failed\n",
                                   Prog_Name,POST_NAME,parm->inum);
                    exit (1);
                  }
                c = buf2;
              }
            x += swide;
          }

        lcp = sarray[x];
      }
    if (z)
      lcp = sarray[e] = 0;
   }
  if (b > buf1)
    if (write(tout,buf1,b-buf1) < 0)
      { fprintf(stderr,"%s: IO write to file %s%d.ktb failed\n",
                       Prog_Name,POST_NAME,parm->inum);
        exit (1);
      }
  if (c > buf2)
    if (write(pout,buf2,c-buf2) < 0)
      { fprintf(stderr,"%s: IO write to file %s%d.pst failed\n",
                       Prog_Name,POST_NAME,parm->inum);
        exit (1);
      }

  parm->npost = (x-range->off)/swide - nbase;
  parm->nelim = nelim;
  parm->nbase = nbase;
  parm->nkmer = nkmer;
  return (NULL);
}

typedef struct
  { int      tid;
    uint8   *bufr;
    int64    nkmer;
    int64    npost;
    int      fail;
    int      tout, pout;
  } CP;

static void *catenate_thread(void *args)
{ CP     *parm = (CP *) args;
  uint8  *bufr = parm->bufr;
  int     tid  = parm->tid;
  int     tout = parm->tout;
  int     pout = parm->pout;

  int tin, pin;
  int p, x;

  parm->fail = 1;

  if (write(tout,&KMER,sizeof(int)) < 0) goto tout_write_fail;
  if (write(tout,&parm->nkmer,sizeof(int64)) < 0) goto tout_write_fail;

  for (p = tid*NTHREADS; p < (tid+1)*NTHREADS; p++)
    { tin = Units[p];
      if (lseek(tin,0,SEEK_SET) < 0)
        { fprintf(stderr,"%s: Rewind of file %s%d.ktb failed\n",Prog_Name,POST_NAME,p);
          return (NULL);
        }
      while (1)
        { x = read(tin,bufr,BUFFER_LEN);
          if (x < 0)
            { fprintf(stderr,"%s: IO read from file %s%d.ktb failed\n",Prog_Name,POST_NAME,p);
              return (NULL);
            }
          if (x == 0)
            break;
          if (write(tout,bufr,x) < 0) goto tout_write_fail;
        }
      close(tin);
    }
  close(tout);

  if (write(pout,&PostBytes,sizeof(int)) < 0) goto pout_write_fail;
  if (write(pout,&ContBytes,sizeof(int)) < 0) goto pout_write_fail;
  if (write(pout,&parm->npost,sizeof(int64)) < 0) goto pout_write_fail;

  for (p = tid*NTHREADS; p < (tid+1)*NTHREADS; p++)
    { pin = Pnits[p];
      if (lseek(pin,0,SEEK_SET) < 0)
        { fprintf(stderr,"%s: Rewind of file %s%d.pst failed\n",Prog_Name,POST_NAME,p);
          return (NULL);
        }
      while (1)
        { x = read(pin,bufr,BUFFER_LEN);
          if (x < 0)
            { fprintf(stderr,"%s: IO read from file %s%d.pst failed\n",Prog_Name,POST_NAME,p);
              return (NULL);
            }
          if (x == 0)
            break;
          if (write(pout,bufr,x) < 0) goto pout_write_fail;
        }
      close(pin);
    }
  close(pout);

  parm->fail = 0;
  return (NULL);

tout_write_fail:
  fprintf(stderr,"%s: IO error writing to part file %s/.%s.ktab.%d\n",Prog_Name,TPATH,TROOT,tid+1);
  return (NULL);

pout_write_fail:
  fprintf(stderr,"%s: IO error writing to part file %s/.%s.post.%d\n",Prog_Name,TPATH,TROOT,tid+1);
  return (NULL);
}

void k_sort(GDB *gdb)
{ uint8 *sarray;
  uint8 *buffer;
  int    p, part, swide;
  int64  panel[256];
  Range  range[NTHREADS];
  SP     sarm[NTHREADS];
  RP     rarm[NTHREADS];
  CP     carm[NTHREADS];
  int64  nelim, nbase, nkmer, npost;
  int64 *prefix;
  int64 *posfix;
#ifndef DEBUG_THREADS
  pthread_t threads[NTHREADS];
#endif

  //  Complete the bucket array by accumulating counts across k-mer 1st bytes

  { int   i, j; 
    int64 x, cum, nelmax;

#ifdef SPACER
    for (j = 0; j < NTHREADS; j++)
      for (i = 0; i < 256; i++)
        Buckets[j][i] *= 2;
#endif

    nelmax = 0;
    cum = 0;
    for (i = 0; i < 256; i++)
      { for (j = 0; j < NTHREADS; j++)
          { x = Buckets[j][i];
            Buckets[j][i] = cum;
            cum += x;
          }
        if (i >= 255 || Select[i] != Select[i+1])
          { if (cum > nelmax)
              nelmax = cum;
            cum = 0;
          }
      }

#ifdef SHOW_FINGERS
    for (i = 0; i < 256; i++)
      { printf(" %3d:",i);
        for (j = 0; j < NTHREADS; j++)
          printf(" %9lld",Buckets[j][i]);
        printf("  :: %d\n",Select[i]);
      }
#endif

    swide  = KBYTES + PostBytes + ContBytes;
    prefix = Malloc(sizeof(int64)*0x1000000,"Allocating prefix array");
    posfix = Malloc(sizeof(int64)*0x10000,"Allocating postfix array");
    sarray = Malloc(nelmax*swide+1,"Allocating sort array");
    buffer = Malloc(NTHREADS*BUFFER_LEN,"Allocating input buffers");
    bzero(prefix,sizeof(int64)*0x1000000);
    bzero(posfix,sizeof(int64)*0x10000);
  }

  //  Must open a separate file descriptor to DB bases for each setup thread !!!

  for (p = 0; p < NTHREADS; p++)
    { sarm[p].tid   = p;
      sarm[p].swide = swide;
      sarm[p].sarr  = sarray;
      sarm[p].buff  = buffer + p*BUFFER_LEN;
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
    { rarm[p].swide  = swide;
      rarm[p].sarr   = sarray;
      rarm[p].buff   = buffer + p*BUFFER_LEN;
      rarm[p].range  = range + p;
      rarm[p].panel  = panel;
      rarm[p].prefix = prefix;
      rarm[p].posfix = posfix;
    }

  for (p = 0; p < NTHREADS; p++)
    { carm[p].tid   = p;
      carm[p].bufr  = buffer + p*BUFFER_LEN;
      carm[p].nkmer = 0;
      carm[p].npost = 0;
    }

  nelim = nbase = nkmer = npost = 0;

  for (part = 0; part < NTHREADS; part++)

    { for (p = 0; p < NTHREADS; p++)
        { sarm[p].in = Units[part*NTHREADS+p];
          sarm[p].inum = part*NTHREADS+p;
        }

#ifdef DEBUG_THREADS
      for (p = 0; p < NTHREADS; p++)
        setup_thread(sarm+p);
#else
      for (p = 1; p < NTHREADS; p++)
        pthread_create(threads+p,NULL,setup_thread,sarm+p);
      setup_thread(sarm);
      for (p = 1; p < NTHREADS; p++)
        pthread_join(threads[p],NULL);
#endif

      { int64 prev, next;

        bzero(panel,sizeof(int64)*256);
        prev = 0;
        for (p = Ksplit[part]; p < Ksplit[part+1]; p++)
          { next = Buckets[NTHREADS-1][p];
            panel[p] = (next - prev)*swide;
            prev = next;
          }

        if (prev == 0)
          { for (p = 0; p < NTHREADS; p++)
              rarm[p].range->off = rarm[p].range->beg = rarm[p].range->end = 0;
          }
        else
          { if (VERBOSE)
              { fprintf(stderr,"\r    Sorting part %d  ",part+1);
                fflush(stderr);
              }

            msd_sort(sarray,next,swide,KBYTES,panel,NTHREADS,range);

#ifdef DEBUG_SORT
            print_table(sarray,swide,next);
#endif
          }
      }

      for (p = 0; p < NTHREADS; p++)   //  open units for T^2 index parts
        { char *name;
          int   inum;

          inum = part*NTHREADS+p;
          name = Numbered_Suffix(POST_NAME,inum,".ktb");
          rarm[p].tout = Units[inum] = open(name,O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
          if (rarm[p].tout < 0)
            { fprintf(stderr,"%s: Cannot open %s for reading & writing\n",Prog_Name,name);
              exit (1);
            }
          unlink(name);

          name = Numbered_Suffix(POST_NAME,inum,".pst");
          rarm[p].pout = Pnits[inum] = open(name,O_RDWR|O_CREAT|O_TRUNC,S_IRWXU);
          if (rarm[p].pout < 0)
            { fprintf(stderr,"%s: Cannot open %s for reading & writing\n",Prog_Name,name);
              exit (1);
            }
          unlink(name);

          rarm[p].inum = inum;
        }

      if (VERBOSE)
        { fprintf(stderr,"\r    Outputing part %d",part+1);
          fflush(stderr);
        }

#ifdef DEBUG_THREADS
      for (p = 0; p < NTHREADS; p++)
        output_thread(rarm+p);
#else
      for (p = 1; p < NTHREADS; p++)
        pthread_create(threads+p,NULL,output_thread,rarm+p);
      output_thread(rarm);
      for (p = 1; p < NTHREADS; p++)
        pthread_join(threads[p],NULL);
#endif

      for (p = 0; p < NTHREADS; p++)
        { nelim += rarm[p].nelim;
          nbase += rarm[p].nbase;
          nkmer += rarm[p].nkmer;
          npost += rarm[p].npost;
          carm[part].nkmer += rarm[p].nkmer;
          carm[part].npost += rarm[p].npost;
        }
    }

  if (VERBOSE)
    { fprintf(stderr,"\n    Concat'ing parts\n");
      fflush(stderr);
    }

  for (p = 1; p < NTHREADS; p++)
    fclose(sarm[p].gdb.seqs);

  for (p = 0; p < NTHREADS; p++)
    { carm[p].tout = open(Catenate(TPATH,"/.",TROOT,
                          Numbered_Suffix(".ktab.",p+1,"")),
                          O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU);
      if (carm[p].tout < 0)
        { fprintf(stderr,"%s: Cannot open part file %s/.%s.ktab.%d for writing\n",
                         Prog_Name,TPATH,TROOT,p+1);
          goto remove_parts;
        }
      carm[p].pout = open(Catenate(TPATH,"/.",TROOT,
                          Numbered_Suffix(".post.",p+1,"")),
                          O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU);
      if (carm[p].pout < 0)
        { fprintf(stderr,"%s: Cannot open part file %s/.%s.post.%d for writing\n",
                         Prog_Name,TPATH,TROOT,p+1);
          goto remove_parts;
        }
    }

#ifdef DEBUG_THREADS
  for (p = 0; p < NTHREADS; p++)
    catenate_thread(carm+p);
#else
  for (p = 1; p < NTHREADS; p++)
    pthread_create(threads+p,NULL,catenate_thread,carm+p);
  catenate_thread(carm);
  for (p = 1; p < NTHREADS; p++)
    pthread_join(threads[p],NULL);
#endif

  for (p = 0; p < NTHREADS; p++)
    if (carm[p].fail)
      goto remove_parts;

  if (VERBOSE)
    { fprintf(stderr,"\r    Done                                           \n");
      fprintf(stderr,"\n  Kept:    %11lld kmers, %11lld(%5.1f%%) positions\n",
                     nkmer,npost-nbase,((npost-nbase)*100.)/npost);
      fprintf(stderr,"  Dropped: %11lld kmers, %11lld(%5.1f%%) positions\n",
                     nelim,nbase,(nbase*100.)/npost);
      fflush(stderr);
    }

  { int   tab;
    int   x;
    int64 maxpre;

    tab = open(Catenate(TPATH,"/",TROOT,".gix"),O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU);
    if (tab < 0)
      { fprintf(stderr,"%s: Cannot open %s/%s.gix for writing\n",Prog_Name,TPATH,TROOT);
        goto remove_parts;
      }
    if (write(tab,&KMER,sizeof(int)) < 0) goto gix_error;
    if (write(tab,&NTHREADS,sizeof(int)) < 0) goto gix_error;
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

    if (write(tab,&PostBytes,sizeof(int)) < 0) goto gix_error;
    if (write(tab,&ContBytes,sizeof(int)) < 0) goto gix_error;
    if (write(tab,&NTHREADS,sizeof(int)) < 0) goto gix_error;
    if (write(tab,&maxpre,sizeof(int64)) < 0) goto gix_error;
    if (write(tab,&FREQ,sizeof(int)) < 0) goto gix_error;
    if (write(tab,&(gdb->ncontig),sizeof(int)) < 0) goto gix_error;
    if (write(tab,Perm,sizeof(int)*gdb->ncontig) < 0) goto gix_error;

    for (x = 1; x < 0x10000; x++)
      posfix[x] += posfix[x-1];
    if (write(tab,posfix,sizeof(int64)*0x10000) < 0) goto gix_error;

    close(tab);
  }
 
  free(buffer);
  free(sarray);
  free(posfix);
  free(prefix);
  return;

gix_error:
  unlink(Catenate(TPATH,"/",TROOT,".gix"));
  fprintf(stderr,"%s: IO error while writing %s/%s.gix\n",Prog_Name,TPATH,TROOT);
  goto remove_parts;

remove_parts:
  for (p = 1; p <= NTHREADS; p++)
    { unlink(Catenate(TPATH,"/.",TROOT,Numbered_Suffix(".ktab.",p,"")));
      unlink(Catenate(TPATH,"/.",TROOT,Numbered_Suffix(".post.",p,"")));
    }
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

    FREQ = 10;
    KMER = 40;
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
          case 'f':
            ARG_NON_NEGATIVE(FREQ,"maximum seed frequency");
            break;
          case 'k':
            ARG_NON_NEGATIVE(KMER,"index k-mer size");
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
        fprintf(stderr,"      -T: Number of threads to use.\n");
        fprintf(stderr,"      -P: Directory to use for temporary files.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -k: index k-mer size\n");
        fprintf(stderr,"      -f: adaptive seed count cutoff\n");
        exit (1);
      }

    if (FREQ > 255)
      { fprintf(stderr,"%s: The maximum allowable frequency cutoff is 255\n",Prog_Name);
        exit (1);
      }
    if ((KMER & 0x3) != 0)
      { fprintf(stderr,"%s: K-mer size must be a multiple of 4\n",Prog_Name);
        exit (1);
      }
  }

  if (VERBOSE)
    StartTime();

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

  if (VERBOSE)
    { if (ftype != IS_GDB)
        if (strcmp(TPATH,".") == 0)
          { fprintf(stderr,"\n  Creating genome data base and index (GDB/GIX) %s.1gdb/gix",TROOT);
            fprintf(stderr," in the current directory\n");
          }
        else
          { fprintf(stderr,"\n  Creating genome data base and index (GDB/GIX)");
            fprintf(stderr," %s.1gdb/gix in directory %s\n",TROOT,TPATH);
          }
      else
        if (strcmp(TPATH,".") == 0)
          fprintf(stderr,"\n  Creating genome index (GIX) %s.gix in the current directory\n",
                         TROOT);
        else
          fprintf(stderr,"\n  Creating genome genome index (GIX) %s.gix in directory %s\n",
                         TROOT,TPATH);
      fflush(stdout);
    }

  if (ftype != IS_GDB)
    { char *command;

      command = Malloc(strlen(spath)+strlen(tpath)+100,"Allocating command string");
      sprintf(command,"FAtoGDB%s %s %s",VERBOSE?" -v":"",spath,tpath);
      if (system(command) != 0)
        { fprintf(stderr,"\n%s: Call to FAtoGDB failed\n",Prog_Name);
          exit (1);
        }
      free(command);
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
  }

  //  Make sure you can open (2 * NTHREADS + 2) * NTHREADS + 1 + tid files at one time.
  //    tid is typically 3 unless using valgrind or other instrumentation.

  { struct rlimit rlp;
    int           tid;
    uint64        nfiles;

    tid = open(".xxx",O_CREAT|O_TRUNC|O_WRONLY,S_IRWXU);
    close(tid);
    unlink(".xxx");

    nfiles = (2*NTHREADS+2)*NTHREADS + 1 + tid + 100 ; // RD 250221 add 100 to allow for a few wrapping processes
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

  //  Open GDB

  POST_NAME = Strdup(Catenate(SORT_PATH,"/.",Numbered_Suffix("post.",getpid(),"."),""),
                     "Allocating post index name");

  Read_GDB(gdb,tpath);
  short_GDB_fix(gdb);

  { int i, l0, l1, l2, l3;   //  Compute byte complement table

    i = 0;
    for (l0 = 3; l0 >= 0; l0 -= 1)
     for (l1 = 12; l1 >= 0; l1 -= 4)
      for (l2 = 48; l2 >= 0; l2 -= 16)
       for (l3 = 192; l3 >= 0; l3 -= 64)
         Comp[i++] = (l3 | l2 | l1 | l0);
  }

  //  Compute NTHREADS 1st byte partitions based on bp frequency

  //    Now handled empirically by split_4bits called at the start of distribute

  Ksplit = Malloc((NTHREADS+1)*sizeof(int),"Allocating split vector");

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

  { int i;   //  Produce perms for length sorted order of contigs
 
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
    { fprintf(stderr,"\n  Partitioning K-mers via pos-lists into %d parts\n",NTHREADS);
      fflush(stderr);
    }

  { int p;   //  Setup distribution bucket array

    Buckets = Malloc(NTHREADS*sizeof(int64 *),"Allocating distribution buckets");
    Buckets[0] = Malloc(NTHREADS*256*sizeof(int64),"Allocating distribution buckets");
    bzero(Buckets[0],NTHREADS*256*sizeof(int64));
    for (p = 1; p < NTHREADS; p++)
      Buckets[p] = Buckets[p-1] + 256;
  }

  { int   p, i, k;         //  Open IO units for distribution and reimport
    char *name;

    Units = Malloc(2*NTHREADS*NTHREADS*sizeof(int),"Allocating IO Units");
    Pnits = Units + NTHREADS*NTHREADS;

    k = 0;
    for (p = 0; p < NTHREADS; p++)
      for (i = 0; i < NTHREADS; i++)
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

  distribute(gdb);   //  Distribute k-mers to 1st byte partitions, encoded as compressed
                     //    relative positions of the given k-mers

  if (VERBOSE)
    { fprintf(stderr,"  Starting sort & index output of each part\n");
      fflush(stderr);
    }

  k_sort(gdb);  //  Reimport the post listings, recreating the k-mers and sorting
                //    them with their posts to produce the final genome index.

  free(Units);

  free(Buckets[0]);
  free(Buckets);
  free(Perm);
  free(DBpost);
  free(DBsplit);
  free(Ksplit);

  Close_GDB(gdb);

  free(POST_NAME);

  free(TROOT);
  free(TPATH);

  free(tpath);
  free(spath);

  if (VERBOSE)
    TimeTo(stderr,0);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
