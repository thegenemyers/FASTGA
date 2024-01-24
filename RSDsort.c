#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
static int    PSIZE;   //  = RSIZE+1
static int    KSIZE;
static int64 *PARTS;
static uint8 *ARRAY;

#ifdef DEBUG_SORT

static void print_table(uint8 *array, int64 asize)
{ int64 x;
  int   k;

  for (x = 0; x < asize; x += RSIZE)
    { printf(" %10lld (%10lld) %2d: ",x/RSIZE,(x+(array-ARRAY))/RSIZE,array[x]);
      for (k = 1; k < KSIZE; k++)
        printf(" %02x",array[x-k]);
      printf(" ::");
      for (k = KSIZE; k < RSIZE; k++)
        printf(" %02x",array[x-k]);
      printf("\n");
    }
}

#endif

static inline void mycpy(uint8 *a, uint8 *b, int n)
{ while (n--)
    *a-- = *b--;
}

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n--)
    { if (*a-- != *b--)
        return (a[1] < b[1] ? -1 : 1);
    }
  return (0);
}

#ifdef IS_SORTED

static inline void sorted(uint8 *array, int64 asize)
{ int64 p, i;
  int   first = 1;
  int64 beg;

  for (p = RSIZE; p < asize; p += RSIZE)
    if (mycmp(array + (p-RSIZE),array + p,KSIZE) > 0)
      { if (first)
          { first = 0;
            beg = (array-ARRAY)/RSIZE;
            printf("A[%lld-%lld]: \n",beg,beg+asize/RSIZE);
          }
        printf("  Not sorted %12lld: ",(p-RSIZE)/RSIZE);
        for (i = 0; i < KSIZE; i++)
          printf(" %02x",array[(p-RSIZE)-i]);
        printf(" vs %12lld: ",p/RSIZE);
        for (i = 0; i < KSIZE; i++)
          printf(" %02x",array[p-i]);
        printf("\n");
      }
}

#endif

static inline void gap_sort(uint8 *array, int asize, int gap, int cmp, int rem)
{ int    i, j;
  uint8  _temp[PSIZE], *temp = _temp+RSIZE;
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
  uint8 *garray;
  
  cmp    = KSIZE-digit;
  rem    = RSIZE-digit;
  garray = array-digit;

  // if (asize > S_thr1)
    // gap_sort(garray,asize,S_gap1,cmp,rem);
  if (asize > S_thr2)
    gap_sort(garray,asize,S_gap2,cmp,rem);
  gap_sort(garray,asize,RSIZE,cmp,rem);
}

static void radix_sort(uint8 *array, int64 asize, int digit, int64 *alive)
{ int64  n, len[256];
  int    y, q, ntop;
  int    nzero[256];

  { uint8 *end[256];
    uint8 *u, *arrow = array - digit;
    int64  o;
    int    rems;
    int    e, x, z;

    uint8 *off[256];
    uint8 _temp[PSIZE], *temp = _temp+RSIZE;
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
        arrow -= 1;
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

  digit += 1;
  if (digit <= KSIZE)
    for (y = 0; y < ntop; y++)
      { q = nzero[y];
        n = len[q];
        if (n > S_thr0)
          radix_sort(array, n, digit, alive);
        else if (n > RSIZE)
          shell_sort(array, n, digit);
        array += n;
      }
}

typedef struct
  { int   beg;
    int   end;
    int64 off;
  } Range;

static void *sort_thread(void *arg) 
{ Range *param = (Range *) arg;

  int      beg   = param->beg;
  int      end   = param->end;
  int64    off   = param->off;

  int      x;
  int64    alive[256];

  if (KSIZE <= 1)
    return (NULL);

  for (x = 0; x < 256; x++)
    alive[x] = 0;

  for (x = beg; x < end; x++)
    { if (PARTS[x] == 0)
        continue;

#ifdef SHOW_STUFF
      printf("Bucket %3d: %12lld - %12lld\n",x,off,off+PARTS[x]);
#endif

      radix_sort(ARRAY + off, PARTS[x], 0, alive);

      off += PARTS[x];
    }

  return (NULL);
}

int rmsd_sort(uint8 *array, int64 nelem, int rsize, int ksize,
              int nparts, int64 *part, int nthreads, Range *parms)
{
#ifndef SHOW_STUFF
  pthread_t     threads[nthreads];
#endif

  int   x, n;
  int64 sum, off;
  int64 asize;
  int64 thr;
  int   beg;

  asize = nelem*rsize;
  array += (rsize-1);

  ARRAY = array;
  PARTS = part;
  RSIZE = rsize;
  PSIZE = rsize+1;
  KSIZE = ksize;

  S_thr0 = THR0*RSIZE;
  S_thr1 = THR1*RSIZE;
  S_thr2 = THR2*RSIZE;
  S_gap1 = GAP1*RSIZE;
  S_gap2 = GAP2*RSIZE;

  n   = 0;
  thr = asize / nthreads;
  off = 0;
  sum = 0;
  for (x = 0; x < nparts; x++)
    if (part[x] > 0)
      break;
  beg = x;
  for (; x < nparts; x++)
    if (part[x] > 0)
      { sum += part[x];
        if (sum >= thr)
          { parms[n].end = x+1;
            parms[n].beg = beg;
            parms[n].off = off;
            n  += 1;
            thr = (asize * (n+1))/nthreads;
            beg = x+1;
            off = sum;
          }
      }
  for (x = n; x < nthreads; x++)
    { parms[n].beg = parms[n].end = parms[x-1].end;
      parms[n].off = asize;
    }
  nthreads = n;

#ifdef SHOW_STUFF
  for (x = 0; x < nthreads; x++)
    sort_thread(parms+x);
#else
  for (x = 1; x < nthreads; x++)
    pthread_create(threads+x,NULL,sort_thread,parms+x);

  sort_thread(parms);

  for (x = 1; x < nthreads; x++)
    pthread_join(threads[x],NULL);
#endif

#ifdef DEBUG_SORT
  sum = 0;
  for (x = 0; x < nparts; x++)
    { print_table(array+sum,part[x]);
      sum += part[x];
    }
#endif

#ifdef IS_SORTED
  sum = 0;
  for (x = 0; x < nparts; x++)
    { sorted(array+sum,part[x]);
      sum += part[x];
    }
#endif

  return (nthreads);
}
