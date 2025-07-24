#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>

#include "gene_core.h"

#undef  IS_SORTED
#undef  SHOW_STUFF
#undef  DEBUG_SORT

#define SMAX  20
#define NMAX   3

#define THR0 15
#define THR1 15
#define THR2  8
#define GAP1  9
#define GAP2  4

static int S_thr0, S_thr1, S_thr2;
static int S_gap1, S_gap2;

static int    RSIZE;
static int    KSIZE;
static uint8 *ARRAY;

#ifdef LCPs

static int LCP_Table[256] =
  { 4, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, };

#endif

#ifdef DEBUG_SORT

static void print_table(uint8 *array, int64 asize)
{ static char  DNA[4] = { 'a', 'c', 'g', 't' };
  static char *fmer[256], _fmer[1280];
  static int64 flag = 0;

  int64 x, post;
  int   k;
  uint8 *pust = (uint8 *) &post;

  if (flag == 0)
    { char *t;
      int   i, l3, l2, l1, l0;

      flag = (0x80llu << (8*(RSIZE-KSIZE)-8));
  
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

  post = 0;
  for (x = 0; x < asize; x += RSIZE)
    { printf(" %10lld %2d: ",x,array[x]);
      for (k = 1; k < KSIZE; k++)
        printf("%s",fmer[array[x+k]]);
      for (k = KSIZE; k < RSIZE; k++)
        pust[k-KSIZE] = array[x+k];
      if (post & flag)
        printf(" %11lld\n",-(post-flag));
      else
        printf(" %11lld\n",post);
      fflush(stdout);
    }
}

#endif

static inline void mycpy(uint8 *a, uint8 *b, int n)
{ while (n--)
    *a++ = *b++;
}

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n--)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}

#ifdef LCPs

static inline int mylcp(uint8 *a, uint8 *b, int n)
{ int i;
  for (i = n; i < KSIZE; i++)
    if (a[i] != b[i])
      return ((i<<2) + LCP_Table[a[i]^b[i]]);
  return (0);
}

#endif

#ifdef IS_SORTED

static inline void sorted(uint8 *array, int64 asize, int digit)
{ int64 p, i;
  int   first = 1;
  int64 beg;

  for (p = RSIZE; p < asize; p += RSIZE)
    if (mycmp(array + (p-RSIZE) + 1,array + p + 1,KSIZE-1) > 0)
      { if (first)
          { first = 0;
            beg = (array-ARRAY)/RSIZE;
            printf("A[%lld-%lld]: %d\n",beg,beg+asize/RSIZE,digit);
          }
        printf("  Not sorted %12lld: ",p/RSIZE);
        for (i = 0; i < KSIZE; i++)
          printf(" %02x",array[(p-RSIZE)+i]);
        printf(" vs ");
        for (i = 0; i < KSIZE; i++)
          printf(" %02x",array[p+i]);
        printf("\n");
      }
}

#endif

static inline void gap_sort(uint8 *array, int asize, int gap, int cmp, int rem)
{ int    i, j;
  uint8  temp[RSIZE];
  uint8 *garray;

  garray = array + gap;
  for (i = gap; i < asize; i += RSIZE)
    { j = i-gap;
      if (mycmp(array+j,array+i,cmp) <= 0)
        continue;
      mycpy(temp,array+i,rem);
      mycpy(array+i,array+j,rem);
      for(j -= gap; j >= 0; j -= gap)
        { if (mycmp(array+j,temp,cmp) <= 0)
            break;
          mycpy(garray+j,array+j,rem);
        }
      mycpy(garray+j,temp,rem);
    }
}

static inline void shell_sort(uint8 *array, int asize, int digit)
{ int    cmp, rem;
  int    i, j;
  uint8 *garray;
#ifdef LCPs
  int    lcp;
#endif
  
  cmp    = KSIZE-digit;
  rem    = RSIZE-digit;
  garray = array+digit;

  // if (asize > S_thr1)
    // gap_sort(garray,asize,S_gap1,cmp,rem);
  if (asize > S_thr2)
    gap_sort(garray,asize,S_gap2,cmp,rem);
  gap_sort(garray,asize,RSIZE,cmp,rem);

  j = 0; 
  for (i = RSIZE; i < asize; i += RSIZE)
#ifdef LCPs
    if ((lcp = mylcp(array+j,array+i,digit)) != 0)
      { array[i] = lcp;
        j = i;
      }
#else
    if (mycmp(garray+j,garray+i,cmp) != 0)
      { array[i] = 1;
        j = i;
      }
#endif
}

static void radix_sort(uint8 *array, int64 asize, int digit, int64 *alive)
{ int64  n, len[256];
  int    y, q, ntop;
  int    nzero[256];
#ifdef LCPs
  int    lcp, p;
#endif

  { uint8 *end[256];
    uint8 *u, *arrow = array + digit;
    int64  o;
    int    rems;
    int    e, x, z;

    uint8 *off[256];
    uint8  temp[RSIZE];
    uint8 *stack[SMAX];

    while (1)
      { e = arrow[0];
        for (o = RSIZE; o < asize; o += RSIZE)
          { if (arrow[o] != e)
              break;
          }
        if (o < asize)
          break;
        digit += 1;
        if (digit >= KSIZE)
          return;
        arrow += 1;
      }

    ntop = 1;
    nzero[0] = e;
    alive[e] = o;
    for (; o < asize; o += RSIZE)
      { x = arrow[o];
        if (alive[x] == 0)
          nzero[ntop++] = x;
        alive[x] += RSIZE;
      }

    u = arrow;
    if (ntop <= NMAX)
      { for (y = 1; y < ntop; y++)
          { x = nzero[y];
            for (z = y-1; z >= 0; z--)
              if (nzero[z] < x)
                break;
              else
                nzero[z+1] = nzero[z];
            nzero[z+1] = x;
          }
        for (y = 0; y < ntop; y++)
          { x = nzero[y];
            len[x] = alive[x];
            alive[x] = 0;
            off[x] = u;
            end[x] = u += len[x];
          }
      }
    else
      { ntop = 0;
        for (x = 0; x < 256; x++)
          if (alive[x] > 0)
            { len[x] = alive[x];
              alive[x] = 0;
              off[x] = u;
              end[x] = u += len[x];
              nzero[ntop++] = x;
            }
      }

    rems = RSIZE-digit;
    for (y = 0; y < ntop; y++)
      { uint8   *p;
        int      t, s;
        int      z;

        x = nzero[y];
        while (off[x] < end[x])
          { t = *off[x];

            if (t == x)
              off[x] += RSIZE;
	    else
              { s = 0;
                stack[s++] = off[x];
                while (s < SMAX)
                  if (t == x)
                    { off[x] += RSIZE;
                      break;
                    }
                  else
                    { u = off[t];
                      while ((z = *u) == t)
                        u += RSIZE;
                      off[t] = u+RSIZE;
                      stack[s++] = u;
                      t = z;
                    }

                u = stack[--s];
                mycpy(temp,u,rems);
	        while (s > 0)
                  { p = stack[--s];
                    mycpy(u,p,rems);
                    u = p;
                  }
                mycpy(u,temp,rems);
              }
          }
      }
  }
  
#ifdef LCPs
  lcp = (digit++ << 2);
  p   = 0;
#else
  digit += 1;
#endif
  if (digit < KSIZE)
    for (y = 0; y < ntop; y++)
      { q = nzero[y];
        n = len[q];
        if (n > S_thr0)
          radix_sort(array, n, digit, alive);
        else if (n > RSIZE)
          shell_sort(array, n, digit);
#ifdef LCPs
        *array = lcp + LCP_Table[p^q];
        p = q;
#else
        *array = 1;
#endif
        array += n;
      }
  else
    for (y = 0; y < ntop; y++)
      { q = nzero[y];
        n = len[q];
#ifdef LCPs
        *array = lcp + LCP_Table[p^q];
        p = q;
#else
        *array = 1;
#endif
        array += n;
      }
}

typedef struct
  { int64 off;
    int64 span;
    int   tid;
  } Packet;

static pthread_mutex_t TMUTEX;
static pthread_cond_t  TCOND;

//  Tstack[0..Tavail-1] is a stack of available threads at any moment.
//  It is always manipulated inside the mutex TMUTEX

static int *Tstack;
static int  Tavail;

static void *sort_thread(void *arg) 
{ Packet *param = (Packet *) arg;
  int64   off   = param->off;
  int64   span  = param->span;
  int     tid   = param->tid;

  int      x;
  int64    alive[256];

  for (x = 0; x < 256; x++)
    alive[x] = 0;

#ifdef SHOW_STUFF
  printf("Bucket %3d: %12lld - %12lld\n",beg,off,off+span);
#endif

  radix_sort(ARRAY + off, span, 1, alive);

  pthread_mutex_lock(&TMUTEX);   //  Put this thread back on the avail stack
    Tstack[Tavail++] = tid;
  pthread_mutex_unlock(&TMUTEX);

  pthread_cond_signal(&TCOND);   //  Signal a thread is available

  return (NULL);
}

void msd_sort(uint8 *array, int64 nelem, int rsize, int ksize,
              int64 *part, int beg, int end, int nthreads) 
{ Packet    packet[nthreads];
  int       tstack[nthreads];
#ifndef SHOW_STUFF
  pthread_t threads[nthreads];
#endif

  int   x, n;
  int64 off;
  int64 asize;

  asize = nelem*rsize;

  ARRAY = array;
  RSIZE = rsize;
  KSIZE = ksize;

  pthread_mutex_init(&TMUTEX,NULL);
  pthread_cond_init(&TCOND,NULL);

  if (KSIZE > 1)

    { S_thr0 = THR0*RSIZE;
      S_thr1 = THR1*RSIZE;
      S_thr2 = THR2*RSIZE;
      S_gap1 = GAP1*RSIZE;
      S_gap2 = GAP2*RSIZE;

      Tstack = tstack;
      for (x = 0; x < nthreads; x++)
        { Tstack[x] = x;
          packet[x].tid = x;
        }
      Tavail = nthreads;

      off = 0;
      for (x = beg; x < end; x++)
        if (part[x] > 0)
          { int tid;

            pthread_mutex_lock(&TMUTEX);
 
            if (Tavail <= 0)                       //  all threads are busy, wait
              pthread_cond_wait(&TCOND,&TMUTEX);

            tid = Tstack[--Tavail];                //  thread tid is available
 
            pthread_mutex_unlock(&TMUTEX);

            // Launching job on thread tid

            packet[tid].off  = off;
            packet[tid].span = part[x];
            off += part[x];

            pthread_create(threads+tid,NULL,sort_thread,packet+tid);
          }

      pthread_mutex_lock(&TMUTEX);   //  Wait for all the jobs to complete
      while (Tavail < nthreads)
        pthread_cond_wait(&TCOND,&TMUTEX);
      pthread_mutex_unlock(&TMUTEX);
    }

#ifdef DEBUG_SORT
  sum = 0;
  for (x = beg; x < end; x++)
    { print_table(array+sum,part[x]);
      sum += part[x];
    }
#endif

#ifdef IS_SORTED
  sum = 0;
  for (x = beg; x < end; x++)
    { sorted(array+sum,part[x],0);
      sum += part[x];
    }
#endif

#ifdef LCPs
  array[0] = 0;
  n = beg;
#else
  array[0] = 1;
#endif
  off = part[beg]; 
  for (x = beg+1; x < end; x++)
    { if (part[x] == 0)
        continue;
#ifdef LCPs
      if ((x&0x300) == (n&0x300))
        array[off] = 1 + LCP_Table[x^n];
      else
        array[off] = 0;
      n = x;
#else
      array[off] = 1;
#endif
      off += part[x];
    }
  array[asize] = 1;
}
