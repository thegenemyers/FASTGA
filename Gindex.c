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
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include "libfastk.h"
#include "DB.h"

static char *Usage = "[-v] [-T<int(8)] [-k<int(40)] -f<int> <source>[.dam]";

static int   FREQ;     //  Adaptemer frequence cutoff parameter
static int   VERBOSE;  //  Verbose output
static char *PATH;
static char *ROOT;
static char *POST_NAME;

static int NTHREADS;   //  by default 8
static int KMER;       //  by default 40, must be >= 12 and divisible by 4
static int KBYTES;     //  Bytes for 2-bit compress k-mer (KMER/4)

#undef  DEBUG_MAP
#undef  DEBUG_THREADS
#undef  DEBUG_SORT

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

typedef struct
  { int     beg;
    int     end;
    uint8  *seq;
    uint8  *neq;
    uint8  *ceq;
    int64   buck[256];
  } BP;

// Thread computes neq[beg..end-1] and ceq[beg..end-1] given seq

static void *pack_thread(void *args)
{ BP *parm = (BP *) args;
  int     beg    = parm->beg;
  int     end    = parm->end;
  uint8  *seq    = parm->seq;
  uint8  *neq    = parm->neq;
  uint8  *ceq    = parm->ceq;

  int    i, k;
  uint8 *s1, *s2, *s3, *n1;
  
  s1 = seq+1;
  s2 = seq+2;
  s3 = seq+3;
  for (i = beg; i < end; i += 4)
    { neq[i] = (seq[i] << 6) | (s1[i] << 4) | (s2[i] << 2) | s3[i];
      ceq[i] = Comp[neq[i]];
    }
  n1 = neq-1;
  for (k = 1; k < 4; k++)
    for (i = beg+k; i < end; i += 4)
      { neq[i] = (n1[i] << 2) | s3[i];
        ceq[i] = Comp[neq[i]];
      }

  return (NULL);
}

// Given neq & ceq, upon completion, for i in [beg,end), seq[i] = file k-mer starting at i should
//   go to and whether it should be complemented or not (0x80 flag) in order to be canonical.
//   Also accumulates # of these k-mers with a given 1st byte in buck.

static void *map_thread(void *args)
{ BP *parm = (BP *) args;
  int     beg    = parm->beg;
  int     end    = parm->end;
  uint8  *seq    = parm->seq;
  uint8  *neq    = parm->neq;
  uint8  *ceq    = parm->ceq;
  int64  *buck   = parm->buck;

  int    kspn = KMER-4;
  int    i, u, v;

  for (i = beg; i < end; i++)
    { for (u = i, v = i+kspn; neq[u] == ceq[v]; u += 4, v -= 4)
        if (u >= v)
          break;
      if (ceq[v] < neq[u])
        { u = ceq[i+kspn];
          seq[i] = Select[u] | 0x80;
        }
      else
        { u = neq[i];
          seq[i] = Select[u];
        }
      buck[u] += 1;
    }

  return (NULL);
}

typedef struct
  { int    tid;
    uint8 *seq;
    int    len;
    int    out;
    uint8 *buffer;
    uint8 *bend;
    int64  post;
    int64  last;
  } DP;

//  The thread scans seq and sends those posts assigned to file tid*Nthreads + ? to their
//    designated file relative to the last emission in compressed form:
//       x0   -> byte,   6-bit post delta with sign x
//       x10  -> short, 13-bit ...
//       x110 -> short, 28-bit ...
//       x111 -> 0x10000000 spacer

static void *distribute_thread(void *args)
{ DP *parm = (DP *) args;
  int    tid    = parm->tid;
  uint8 *seq    = parm->seq;
  int    len    = parm->len;
  int    out    = parm->out;
  uint8 *buffer = parm->buffer;
  uint8 *bend   = parm->bend;

  int64  post, last;
  uint8 *lust = (uint8 *) (&last);
  uint8 *b;
  int    i, u;

  b = buffer;
  last = parm->last;
  post = parm->post;
  for (i = 0; i < len; i++, post++)
    { u = seq[i];
      if ((u & 0x7f) == tid)
        { last = post-last;
          if (last < 0x3f)
            { if (u & 0x80)
                *b++ = 0x80 | last;
              else
                *b++ = last;
            }
          else if (last < 0x1fff)
            { if (u & 0x80)
                last |= 0xc000;
              else
                last |= 0x4000;
              *b++ = lust[1];
              *b++ = lust[0];
            }
          else
            { while (last >= 0x10000000)
                { if (u & 0x80)
                    *b++ = 0xf0;
                  else
                    *b++ = 0x70;
                  if (b >= bend)
                    { write(out,buffer,b-buffer);
                      b = buffer;
                    }
                  last -= 0x10000000;
                }
              if (u & 0x80)
                last |= 0xe0000000;
              else
                last |= 0x60000000;
              *b++ = lust[3];
              *b++ = lust[2];
              *b++ = lust[1];
              *b++ = lust[0];
            }
          if (b >= bend)
            { write(out,buffer,b-buffer);
              b = buffer;
            }
          last = post;
        }
    }
  if (b > buffer)
    write(out,buffer,b-buffer);

  parm->last = last;
  parm->post = post;
  return (NULL);
}

// Distribute posts to the appropriate NTHREADS^2 files based on section of the DB and 1st byte
//   of its canonical k-mer

void distribute(DAZZ_DB *DB)
{ uint8 *seq;
  int    len, ren;
  uint8 *neq, *ceq;
  int    p, r, i, j;

  DP        parm[NTHREADS];
  BP        barm[NTHREADS];
#ifndef DEBUG_THREADS
  pthread_t threads[NTHREADS];
#endif

  seq = (uint8 *) New_Read_Buffer(DB);   //  Allocate work vectorss and set up fixed parts of
  neq = (uint8 *) New_Read_Buffer(DB);   //    of the thread records
  ceq = (uint8 *) New_Read_Buffer(DB);

  parm[0].buffer = Malloc(1000000*NTHREADS,"IO Buffer");
  for (i = 0; i < NTHREADS; i++)
    { parm[i].tid    = i;
      parm[i].seq    = seq;
      parm[i].buffer = parm[0].buffer + 1000000*i;
      parm[i].bend   = parm[i].buffer +  999996;
      parm[i].post   = 0;
      parm[i].last   = 0;
    }

  for (i = 0; i < NTHREADS; i++)
    { barm[i].seq  = seq;
      barm[i].neq  = neq;
      barm[i].ceq  = ceq;
      bzero(barm[i].buck,256*sizeof(int64));
    }

  //  For each segment of the input ...

  for (p = 0; p < NTHREADS; p++)
                                     //  Open the NTHREAD files to recieve posts from this segment
    { for (i = 0; i < NTHREADS; i++)
        { parm[i].out = open(Numbered_Suffix(POST_NAME,i*NTHREADS+p,".idx"),
                                O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU);
          if (parm[p].out < 0)
            { fprintf(stderr,"%s: Cannot open %s%d.idx for writing\n",
                             Prog_Name,POST_NAME,i*NTHREADS+p);
              exit (1);
            }
        }

      //  For each contig in the segment ...

      ren = DBsplit[p+1];
      for (r = DBsplit[p]; r < ren; r++)
        { Load_Read(DB,r,(char *) seq,0);   //  Load the contig
    
          len = DB->reads[r].rlen;
    
#ifdef DEBUG_MAP
          printf("Src:");
          for (i = 0; i < len; i++)
            printf(" %d",seq[i]);
          printf("\n");
#endif
          //  In segments each thread computes neq and ceq from seq
    
          len -= 3;
          barm[0].beg = 0;
          for (i = 1; i < NTHREADS; i++)
            barm[i-1].end = barm[i].beg = (((int64) i)*len)/NTHREADS;
          barm[NTHREADS-1].end = len;
    
#ifdef DEBUG_THREADS
          for (i = 0; i < NTHREADS; i++)
            pack_thread(barm+i);
#else
          for (i = 1; i < NTHREADS; i++)
            pthread_create(threads+i,NULL,pack_thread,barm+i);
          pack_thread(barm);
          for (i = 1; i < NTHREADS; i++)
            pthread_join(threads[i],NULL);
#endif
    
#ifdef DEBUG_MAP
          printf("Neq:");
          for (i = 0; i < len; i++)
            printf(" %02x",neq[i]);
          printf("\n");

          printf("Ceq:");
          for (i = 0; i < len; i++)
            printf(" %02x",ceq[i]);
          printf("\n");
#endif
          //  In segments each thread computes processed seq array from neq,ceq
    
          len -= (KMER-4);
          barm[0].beg = 0;
          for (i = 1; i < NTHREADS; i++)
            barm[i-1].end = barm[i].beg = (((int64) i)*len)/NTHREADS;
          barm[NTHREADS-1].end = len;
    
#ifdef DEBUG_THREADS
          for (i = 0; i < NTHREADS; i++)
            map_thread(barm+i);
#else
          for (i = 1; i < NTHREADS; i++)
            pthread_create(threads+i,NULL,map_thread,barm+i);
          map_thread(barm);
          for (i = 1; i < NTHREADS; i++)
            pthread_join(threads[i],NULL);
#endif
    
#ifdef DEBUG_MAP
          printf("Out:");
          for (i = 0; i < len; i++)
            printf(" (%d)%d %02x %02x\n",(seq[i]&0x80)!=0,seq[i]&0x7f,neq[i],ceq[i+KMER-4]);
          printf("\n");
#endif

          //  Each thread scans *all* of seq, writing the k-mers assigned to "its" file

          for (i = 0; i < NTHREADS; i++)
            parm[i].len = len;
    
#ifdef DEBUG_THREADS
          for (i = 0; i < NTHREADS; i++)
            distribute_thread(parm+i);
#else
          for (i = 1; i < NTHREADS; i++)
            pthread_create(threads+i,NULL,distribute_thread,parm+i);
          distribute_thread(parm);
          for (i = 1; i < NTHREADS; i++)
            pthread_join(threads[i],NULL);
#endif
    
          for (i = 0; i < NTHREADS; i++)
            parm[i].post += (KMER-1);
        }

      //  Accumulate bucket counts into Bucket vector for data panel

      { int64 *buck;

        buck = Buckets[p];
        for (i = 0; i < NTHREADS; i++)
          { close(parm[i].out);
            for (j = 0; j < 256; j++)
              buck[j] += barm[i].buck[j];
            parm[i].last = parm[i].post;
          }
      }
    }
 
  free(parm[0].buffer);
  free(seq-1);
  free(neq-1);
  free(ceq-1);

  //  Complete the bucket array by accumulating counts across k-mer 1st bytes

  { int   i, j; 
    int64 cum;

    cum = Buckets[NTHREADS-1][0];
    for (i = 1; i < 256; i++)
      { for (j = 0; j < NTHREADS; j++)
          Buckets[j][i] += cum;
        if (i < 255 && Select[i] == Select[i+1])
          cum = Buckets[NTHREADS-1][i];
        else
          cum = 0;
      }
  }
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
  { DAZZ_DB  DB;
    int      tid;
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
  DAZZ_DB *DB   = &(parm->DB);
  int64   *buck = Buckets[tid];

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

  seq = (uint8 *) New_Read_Buffer(DB);
  neq = (uint8 *) New_Read_Buffer(DB);
  ceq = (uint8 *) New_Read_Buffer(DB);
  keq = ceq+(KMER-4);

  flag = (0x1ll << (8*ContBytes-1));

  bend = buffer + read(in,buffer,1000000);
  if (bend-buffer < 1000000)
    btop = bend;
  else
    btop = bend-4;
  b = buffer;

  ncntg = DBsplit[tid];
  post  = DBpost[tid];
  nextpost = post;
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
            { len = DB->reads[ncntg].rlen;
              basepost = nextpost;
              ncntg    += 1;
              nextpost += len;
            }
          while (post >= nextpost);

          Load_Read(DB,ncntg-1,(char *) seq,0);

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
            x = sarr + swide * buck[*n]++;
            *x++ = 0;
            for (i = 4; i < KMER; i += 4)
              *x++ = n[-i];
            for (i = 0; i < PostBytes; i++)
              *x++ = bust[i];
            for (i = 0; i < ContBytes; i++)
              *x++ = nust[i];
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
          }
      }

      if (b >= btop)           //  refill input buffer if needed
        { int ex = bend-b;
          memcpy(buffer,b,ex); 
          bend = buffer+ex;
          bend += read(in,bend,1000000-ex);
          if (bend == buffer)
            break;
          if (bend-buffer < 1000000)
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
  { int      swide;
    int64   *panel;
    uint8   *sarr;
    uint8   *buff;
    Range   *range;
    int      tout;
    int      pout;
    int64   *prefix;
    int64    nelim;
    int64    nbase;
    int64    nkmer;
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

  uint8 *buf1 = parm->buff;
  uint8 *buf2 = buf1 + 500000;
  uint8 *bed1 = (buf1 + 500000) - (KBYTES-1);
  uint8 *bed2 = (buf2 + 500000) - (PostBytes+ContBytes);

  int64   nelim, nkmer, nbase;
  int64   x, e, y;
  int     o, w, k, z, lcp;
  uint8 *_w = (uint8 *) &w;
  uint8  *b, *c;

  nelim = nkmer = nbase = 0;
  o = KMER;
  write(tout,&o,sizeof(int));       // Skip prolog of each part file
  write(tout,&nkmer,sizeof(int64));
  write(pout,&o,sizeof(int));
  write(pout,&o,sizeof(int));
  write(pout,&nkmer,sizeof(int64));

  b = buf1;
  c = buf2;
  x = range->off;
  lcp = sarray[x];
  prefix += (beg << 16);
  for (o = beg; o < end; o++, prefix += 0x10000)
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

        prefix[(sarray[x+1] << 8) | sarray[x+2]] += 1;

        for (k = 3; k < KBYTES; k++)
          *b++ = sarray[x+k];
        *b++ = _w[0];
        *b++ = lcp;
        if (b >= bed1)
          { write(tout,buf1,b-buf1);
            b = buf1;
          }
        nkmer += 1;

        while (x < y)
          { for (k = KBYTES; k < swide; k++)
              *c++ = sarray[x+k];
            if (c >= bed2)
              { write(pout,buf2,c-buf2);
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
    write(tout,buf1,b-buf1);
  if (c > buf2)
    write(pout,buf2,c-buf2);

  lseek(tout,0,SEEK_SET);        //  Write prologs of each part file
  o = KMER;
  write(tout,&o,sizeof(int));
  write(tout,&nkmer,sizeof(int64));

  lseek(pout,0,SEEK_SET);
  o = PostBytes;
  write(pout,&o,sizeof(int));
  o = ContBytes;
  write(pout,&o,sizeof(int));
  e = (x-range->off)/swide - nbase;
  write(pout,&e,sizeof(int64));

  close(tout);
  close(pout);

  parm->nelim = nelim;
  parm->nbase = nbase;
  parm->nkmer = nkmer;
  return (NULL);
}

//  PLEASE DOCUMENT ME

void k_sort(DAZZ_DB *DB)
{ uint8 *sarray;
  uint8 *buffer;
  int    p, part, swide;
  int64  panel[256];
  Range  range[NTHREADS];
  SP     sarm[NTHREADS];
  RP     rarm[NTHREADS];
  int64  nelim, nbase, nkmer;
  int64 *prefix;
#ifndef DEBUG_THREADS
  pthread_t threads[NTHREADS];
#endif

  { int64 x, nelmax;

    nelmax = 0;
    for (p = 0; p < NTHREADS; p++)
      { x = Buckets[NTHREADS-1][Ksplit[p+1]-1];
        if (nelmax < x)
          nelmax = x;
      }


    swide  = KBYTES + PostBytes + ContBytes;
    prefix = Malloc(sizeof(int64)*0x1000000,"Prefix array");
    sarray = Malloc(nelmax*swide+1,"Sort Array");
    buffer = Malloc(NTHREADS*1000000,"Input Buffers");
    if (prefix == NULL || sarray == NULL || buffer == NULL)
      exit (1);
    bzero(prefix,sizeof(int64)*0x1000000);
  }

  { int i, j;

    for (i = 255; i >= 0; i--)
      { for (j = NTHREADS-1; j >= 1; j--)
          Buckets[j][i] = Buckets[j-1][i];
        if (i == 0 || Select[i] != Select[i-1])
          Buckets[0][i] = 0;
        else
          Buckets[0][i] = Buckets[NTHREADS-1][i-1];
      }
  }

  //  Must open a separate file descriptor to DB bases for each setup thread !!!

  for (p = 0; p < NTHREADS; p++)
    { sarm[p].tid   = p;
      sarm[p].swide = swide;
      sarm[p].sarr  = sarray;
      sarm[p].buff  = buffer + p*1000000;
      sarm[p].DB    = *DB;
      if (p > 0)
        { sarm[p].DB.bases = fopen(Catenate(DB->path,"","",".bps"),"r");
          if (sarm[p].DB.bases == NULL)
            { fprintf(stderr,"%s: Cannot open another copy of DB\n",Prog_Name);
              exit (1);
            }
        }
    }

  for (p = 0; p < NTHREADS; p++)
    { rarm[p].swide  = swide;
      rarm[p].sarr   = sarray;
      rarm[p].buff   = buffer + p*1000000;
      rarm[p].range  = range + p;
      rarm[p].panel  = panel;
      rarm[p].prefix = prefix;
    }

  nelim = nbase = nkmer = 0;

  for (part = 0; part < NTHREADS; part++)

    { for (p = 0; p < NTHREADS; p++)
        { sarm[p].in = open(Numbered_Suffix(POST_NAME,part*NTHREADS+p,".idx"),O_RDONLY);
          if (sarm[p].in < 0)
            { fprintf(stderr,"%s: Cannot open %s%d.idx for reading\n",
                             Prog_Name,POST_NAME,part*NTHREADS+p);
              exit (1);
            }
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

      for (p = 0; p < NTHREADS; p++)
        unlink(Numbered_Suffix(POST_NAME,part*NTHREADS+p,".idx"));

      { int64 prev, next;

        bzero(panel,sizeof(int64)*256);
        prev = 0;
        for (p = Ksplit[part]; p < Ksplit[part+1]; p++)
          { next = Buckets[NTHREADS-1][p];
            panel[p] = (next - prev)*swide;
            prev = next;
          }

        if (VERBOSE)
          { fprintf(stdout,"\r    Sorting part %d  ",part+1);
            fflush(stdout);
          }

        msd_sort(sarray,next,swide,KBYTES,panel,NTHREADS,range);

#ifdef DEBUG_SORT
        print_table(sarray,swide,next);
#endif
      }

      for (p = 0; p < NTHREADS; p++)
        { rarm[p].tout = open(Catenate(PATH,"/.",ROOT,
                                    Numbered_Suffix(".ktab.",part*NTHREADS+p+1,"")),
                                O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU);
          if (rarm[p].tout < 0)
            { fprintf(stderr,"%s: Cannot open %s/.%s.ktab.%d for writing\n",
                             Prog_Name,PATH,ROOT,part*NTHREADS+p+1);
              exit (1);
            }
          rarm[p].pout = open(Catenate(PATH,"/.",ROOT,
                                    Numbered_Suffix(".post.",part*NTHREADS+p+1,"")),
                                O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU);
          if (rarm[p].pout < 0)
            { fprintf(stderr,"%s: Cannot open %s/.%s.post.%d for writing\n",
                             Prog_Name,PATH,ROOT,part*NTHREADS+p+1);
              exit (1);
            }
        }

      if (VERBOSE)
        { fprintf(stdout,"\r    Outputing part %d",part+1);
          fflush(stdout);
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
        }
    }

  for (p = 1; p < NTHREADS; p++)
    fclose(sarm[p].DB.bases);

  if (VERBOSE)
    { int64 npost = DB->totlen - DB->treads*(KMER-1);

      fprintf(stdout,"\r    Done                                           \n");
      fprintf(stdout,"\n  Kept:    %11lld kmers, %11lld(%5.1f%%) positions\n",
                     nkmer,npost-nbase,((npost-nbase)*100.)/npost);
      fprintf(stdout,"  Dropped: %11lld kmers, %11lld(%5.1f%%) positions\n\n",
                     nelim,nbase,(nbase*100.)/npost);
      fflush(stdout);
    }

  { FILE *tab, *idx;
    int   x;
    int64 maxpre;

    tab = fopen(Catenate(PATH,"/",ROOT,".ktab"),"w");
    if (tab == NULL)
      { fprintf(stderr,"%s: Cannot open %s/%s.ktab for writing\n",Prog_Name,PATH,ROOT);
        exit (1);
      }
    x = KMER;
    fwrite(&x,sizeof(int),1,tab);
    x = NTHREADS*NTHREADS;
    fwrite(&x,sizeof(int),1,tab);
    x = 1;
    fwrite(&x,sizeof(int),1,tab);
    x = 3;
    fwrite(&x,sizeof(int),1,tab);

    maxpre = prefix[0];
    for (x = 1; x < 0x1000000; x++)
      { if (prefix[x] > maxpre)
          maxpre = prefix[x];
        prefix[x] += prefix[x-1];
      }
    fwrite(prefix,sizeof(int64),0x1000000,tab);
    fclose(tab);

    idx = fopen(Catenate(PATH,"/",ROOT,".post"),"w");
    if (idx == NULL)
      { fprintf(stderr,"%s: Cannot open %s/%s.post for writing\n",Prog_Name,PATH,ROOT);
        exit (1);
      }
    x = PostBytes;
    fwrite(&x,sizeof(int),1,idx);
    x = ContBytes;
    fwrite(&x,sizeof(int),1,idx);
    x = NTHREADS;                       //  # of files = square, but would like # of threads too
    fwrite(&x,sizeof(int),1,idx);
    fwrite(&maxpre,sizeof(int64),1,idx);
    fwrite(&FREQ,sizeof(int),1,idx);
    fwrite(&(DB->treads),sizeof(int),1,idx);
    fwrite(Perm,sizeof(int),DB->treads,idx);
    fclose(idx);
  }
 
  free(buffer);
  free(sarray);
  free(prefix);
}


/***********************************************************************************************
 *
 *   MAIN
 *
 **********************************************************************************************/

static void short_DB_fix(DAZZ_DB *DB)
{ int i;

  if (DB->treads >= NTHREADS)
    return;

  //  Add additional reads of length KMER that are the first bit of the 0th read/contig
  //    Mark as "fake" with -1 in origin field.

  DB->reads = Realloc(DB->reads,sizeof(DAZZ_READ)*NTHREADS,"Reallocating DB read vector");
  for (i = DB->treads; i < NTHREADS; i++)
    { DB->reads[i] = DB->reads[0];
      DB->reads[i].origin = -1;
      DB->reads[i].rlen   = KMER;
    }
  DB->treads = NTHREADS;
}

static DAZZ_READ *READS;

static int LSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  return (READS[y].rlen - READS[x].rlen);
}

int main(int argc, char *argv[])
{ DAZZ_DB _DB, *DB = &_DB;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;


    ARG_INIT("Gindex");

    FREQ = -1;
    KMER = 40;
    NTHREADS = 8;

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
            ARG_NON_NEGATIVE(KMER,"maximum seed frequency");
            break;
          case 'T':
            ARG_NON_NEGATIVE(NTHREADS,"maximum seed frequency");
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    KBYTES  = (KMER>>2);
    if (argc != 2 || FREQ < 0)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
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

  //  Open and trim DB

  PATH  = PathTo(argv[1]);
  ROOT  = Root(argv[1],".dam");
  POST_NAME = Strdup(Numbered_Suffix("._post.",getpid(),"."),"Allocating temp name");

  if (Open_DB(argv[1],DB) < 0)
    { fprintf(stderr,"%s: Cannot open Dazzler DB %s\n",Prog_Name,argv[1]);
      exit (1);
    }
  Trim_DB(DB);
  short_DB_fix(DB);

  { int i, l0, l1, l2, l3;   //  Compute byte complement table

    i = 0;
    for (l0 = 3; l0 >= 0; l0 -= 1)
     for (l1 = 12; l1 >= 0; l1 -= 4)
      for (l2 = 48; l2 >= 0; l2 -= 16)
       for (l3 = 192; l3 >= 0; l3 -= 64)
         Comp[i++] = (l3 | l2 | l1 | l0);
  }

  { int    i, n;     //  Compute NTHREADS 1st byte partitions based on bp frequency
    double p, t;

    Ksplit = Malloc((NTHREADS+1)*sizeof(int),"Allocating Kmer Split");
    if (Ksplit == NULL)
      exit (1);

    p = 0.;
    n = 0;
    t = 1./NTHREADS;
    Ksplit[0] = 0;
    for (i = 0; i < 256; i++)
      { p += DB->freq[i >> 6] * DB->freq[(i >> 4) & 0x3] 
           * DB->freq[(i >> 2) & 0x3] * DB->freq[i&0x3];
        while (p*(2.-p) > t)
          { n += 1;
            Ksplit[n] = i;
            t = (n+1.)/NTHREADS;
          }
        Select[i] = n;
      }
    Ksplit[NTHREADS] = 256;
  }

  { int64 npost, range, cum, t;   //  Compute DB split into NTHREADS parts
    int   p, r, len;

    DBsplit = Malloc((NTHREADS+1)*sizeof(int),"Allocating DB Split");
    DBpost  = Malloc((NTHREADS+1)*sizeof(int64),"Allocating DB Split");
    if (DBsplit == NULL || DBpost == NULL)
      exit (1);

    npost = DB->totlen;
    cum   = 0;
    range = 0;

    DBsplit[0] = 0;
    DBpost [0] = 0;
    p = 1;
    t = (npost*p)/NTHREADS;
    for (r = 0; r < DB->treads; r++)
      { len = DB->reads[r].rlen;
        cum += len;
        if (cum >= t)
          { DBsplit[p] = r+1;
            DBpost [p] = cum;
            p += 1;
            t = (npost*p)/NTHREADS;
          }
        if (range < len)
          range = len;
      }
    DBsplit[NTHREADS] = DB->treads;
    DBpost [NTHREADS] = npost;

    PostBytes = 0;                 //  # of bytes for encoding a post
    cum = 1;
    while (cum < range)
      { cum *= 256;
        PostBytes += 1;
      }

    range = 2*DB->treads;
    ContBytes = 0;                 //  # of bytes for encoding a contig + sign bit
    cum = 1;
    while (cum < range)
      { cum *= 256;
        ContBytes += 1;
      }
  }

  { int i;   //  Produce perms for length sorted order of contigs
 
    Perm = Malloc(2*DB->treads*sizeof(int),"Allocating sort permutation\n");
    InvP = Perm + DB->treads;

    for (i = 0; i < DB->treads; i++)
      Perm[i] = i;
  
    READS = DB->reads;
    qsort(Perm,DB->treads,sizeof(int),LSORT);

    for (i = 0; i < DB->treads; i++)
      InvP[Perm[i]] = i;
  }

  if (VERBOSE)
    { fprintf(stdout,"  Partitioning K-mers via pos-lists into %d parts\n",NTHREADS);
      fflush(stdout);
    }

  { int p;   //  Setup distribution bucket array

    Buckets = Malloc(NTHREADS*sizeof(int64 *),"Allocating Distribution Buckets");
    if (Buckets == NULL)
      exit (1);
    Buckets[0] = Malloc(NTHREADS*256*sizeof(int64),"Allocating Distribution Buckets");
    if (Buckets[0] == NULL)
      exit (1);
    bzero(Buckets[0],NTHREADS*256*sizeof(int64));
    for (p = 1; p < NTHREADS; p++)
      Buckets[p] = Buckets[p-1] + 256;
  }

  distribute(DB);   //  Distribute k-mers to 1st byte partitions, encoded as compressed
                    //    relative positions of the given k-mers

  if (VERBOSE)
    { fprintf(stdout,"  Starting sort & index output of each part\n");
      fflush(stdout);
    }

  k_sort(DB);       //  Reimport the post listings, recreating the k-mers and sorting
                    //    them with their posts to produce final genome index.

  free(Buckets[0]);
  free(Buckets);
  free(Perm);
  free(DBpost);
  free(DBsplit);
  free(Ksplit);

  Close_DB(DB);

  free(POST_NAME);
  free(ROOT);
  free(PATH);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
