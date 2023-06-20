/*******************************************************************************************
 *
 *  Approximate match filter of Rasmussen, Stoye, Myers (RECOMB 2005)
 *
 *  Author:  Gene Myers
 *  Date  :  August 2004
 *
 ********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "filter.h"
#include "filter_params.h"

#undef  TEST_TABLE
#undef  TEST_DIAG
#undef  TEST_CUTOFF

#define OUTER_LOOP  /* When alone, turn on REPORT_SIZES */
#define INNER_LOOP

#undef  DEBUG_FILTER
#undef  DEBUG_HITS

#undef  REPORT_SIZES
#undef  REPORT_MEMORY

/* Shared index and filter arrays */

typedef struct {
  int minim;
  int maxim;
  int count;
} DiagRecord;

static int   Kmask = -1;
static int  *Table;          /* [0..Kmask+1] */
static int  *Tuples = NULL;  /* [0..<Seqlen>-KMERLEN] */
static int   Map[128];

static int   LenMax = -1;
static char *Alast = NULL;

static int  TooFrequent;

static DiagRecord *DiagVec; /* [0..(Alen + BINSIZE - 1) >> BINSHIFT] */


/*** UTILITY ROUTINES ***/

static void OutOfMemory(char *where)
{ fprintf(stderr,"Filter: Out of memory (%s)\n",where);
  exit (1);
}

static void *ckalloc(int size, char *where)
{ void *p;
  p = malloc(size);
  if (p == NULL) OutOfMemory(where);
#ifdef REPORT_MEMORY
  fprintf(stderr,"Alloc:   %10d (%s)\n",size,where);
#endif
  return (p);
}

static void *ckrealloc(void *p, int size, char *where)
{ p = realloc(p,size);
  if (p == NULL) OutOfMemory(where);
#ifdef REPORT_MEMORY
  fprintf(stderr,"Realloc: %10d (%s)\n",size,where);
#endif
  return (p);
}


/*** INDEX MATCH PARAMETERS ***/

static int KMERLEN  = 13;
static int HITMIN   =  9;
static int MAXDIFFS =  4;
static int TRAPMIN  = 60;
static int BINSHIFT =  3;
static int BINMASK  =  7;
static int BINSIZE  = 12;
static int SUPPRESS =  0;

void Reset_Filter()
{ free(Tuples);
  free(Table);
  Kmask  = -1;
  LenMax = -1;
  Tuples = NULL;
  Alast  = NULL;
}

int Set_Filter_Params(double Erate, int MinLen, int KmerLen, int BinShift, int Suppress)
{ int HitThresh;

  if (KmerLen <= 1)
    return (1);
  if (KmerLen > MinLen)
    return (1);
  if (Erate < 0. || Erate >= 1.)
    return (1);
  if (KmerLen > max_gram(Erate))
    { printf("Max_Gram = %d\n",max_gram(Erate));
      return (1);
    }

  HitThresh = max_hits(Erate,KmerLen,MinLen);
  if (HitThresh < 1)
    { printf("Hitthresh = %d\n",HitThresh);
      return (1);
    }
  if (max_diffs(Erate,KmerLen,HitThresh) >= (1 << BinShift))
    { printf("Max_diffs = %d\n",max_diffs(Erate,KmerLen,HitThresh));
      return (1);
    }

  if (KmerLen != KMERLEN)
    Reset_Filter();

  KMERLEN  = KmerLen;
  HITMIN   = HitThresh;
  MAXDIFFS = max_diffs(Erate,KMERLEN,HITMIN);
  TRAPMIN  = max_space(KMERLEN,HITMIN,MAXDIFFS);
  BINSHIFT = BinShift;
  BINMASK  = (1 << BinShift) - 1;
  BINSIZE  = (1 << BinShift) + MAXDIFFS;
  SUPPRESS = Suppress;

  return (0);
}


/*** INDEX CONSTRUCTION */


/* Build index table for sequence S of length Slen. */

static void TableBuild(char *S, int Slen)
{ int   i, c;
  int   x, h;
  char *s;

  s = S+(KMERLEN-1);

  for (c = 0; c <= Kmask; c++)
    Table[c] = 0;

  h = 0;
  c = 0;
  for (i = 0; i < KMERLEN-1; i++)
    { x = Map[(int) (S[i])];
      if (x >= 0)
        c = (c << 2) | x;
      else
        { c = 0; h = i+1; }
    }
  for (i = 0; i <= Slen-KMERLEN; i++)
    { x = Map[(int) (s[i])];
      if (x >= 0)
        c = ((c << 2) | x) & Kmask;
      else
        { c = 0; h = i + KMERLEN; }
      if (i >= h)
        Table[c+1] += 1;
    }

  for (c = 2; c <= Kmask; c++)
    Table[c] += Table[c-1];

  h = 0;
  c = 0;
  for (i = 0; i < KMERLEN-1; i++)
    { x = Map[(int) (S[i])];
      if (x >= 0)
        c = (c << 2) | x;
      else
        { c = 0; h = i+1; }
    }
  for (i = 0; i <= Slen-KMERLEN; i++)
    { x = Map[(int) (s[i])];
      if (x >= 0)
        c = ((c << 2) | x) & Kmask;
      else
        { c = 0; h = i + KMERLEN; }
      if (i >= h)
        Tuples[Table[c]++] = i;
    }

  for (c = Kmask; c >= 0; c--)
    Table[c+1] = Table[c];
  Table[0] = 0;

#ifdef TEST_CUTOFF
  x = 0;
  for (c = 0; c <= Kmask; c++)
    if (Table[c+1] - Table[c] >= TooFrequent)
      x +=1;
  fprintf(stderr,"Table built: %d %d-mers are too frequent (>= %d)\n",x,KMERLEN,TooFrequent);
#endif
  
#ifdef TEST_TABLE
  { int i, c, j;

    static char Convert[] = { 'a', 'c', 'g', 't' };

    for (c = 0; c <= Kmask; c++)
      { printf("Table[%d = ",c);
        for (i = KMERLEN-1; i >= 0; i--)
          printf("%c",Convert[c>>(i*2) & 0x3]);
        printf("]\n");
        for (i = Table[c]; i < Table[c+1]; i++)
          printf("  %d (%.*s)\n",Tuples[i],KMERLEN,S+Tuples[i]);
      }
  }
#endif
}


/*** FILTER ***/

/* Apply index to find filtered hits between sequences, returning pointer to
   array of HitTube of length in the integer pointed at by Hitlen             */

HitTube *Match_Filter(char *A, int Alen, char *B, int Blen, int asym, int *Hitlen)
{ int        HitMax, BinMax;
  HitTube   *HitList;
  int        Divpt;
  int        bintop, hits;
#ifdef REPORT_SIZES
  long long  sum;
#endif

  if (Kmask < 0)
    { int i;
      for (i = 0; i < 128; i++)
        Map[i] = -1;
      Map['a'] = Map['A'] = 0;
      Map['c'] = Map['C'] = 1;
      Map['g'] = Map['G'] = 2;
      Map['t'] = Map['T'] = 3;
      Kmask = (1 << (2*KMERLEN)) - 1;
      Table = ckalloc(sizeof(int)*(Kmask+2),"K-mer index");
    }

  if (Alen > LenMax)
    { LenMax  = 1.2*Alen + 5000;
      Tuples  = ckrealloc(Tuples,sizeof(int)*LenMax,"Tuple table");
    }

  if (SUPPRESS == 0)
    TooFrequent = Alen+1;
  else
    TooFrequent = SUPPRESS;

  if (A != Alast)
    TableBuild(A,Alen);
  Alast = A;

  if (B == NULL)
    { *Hitlen = 0;
      return NULL;
    }

  if (asym)
    Divpt = Alen;
  else
    Divpt = -1;

  BinMax  = ((LenMax + BINSIZE - 1) >> BINSHIFT) + 1;
  DiagVec = (DiagRecord *) ckalloc(sizeof(DiagRecord)*BinMax,"Bin vector");

  HitMax = 100000;
  if (HitMax < Alen/48) HitMax = Alen/48;
  HitList = ckalloc(sizeof(HitTube)*HitMax,"Hit list");

  bintop = ((Alen + BINSIZE - 1) >> BINSHIFT) + 1;

#ifdef OUTER_LOOP

  { int i, c;
    int x, h;
    char *b;
    int id, ia, ib;

    for (i = 0; i < bintop; i++)
      { DiagRecord *dp;
        dp = DiagVec + i;
        dp->count = 0;
        dp->maxim = 0;
      }

    hits = 0;
#ifdef REPORT_SIZES
    sum = 0;
#endif
    h = 0;
    c = 0;
    for (i = 0; i < KMERLEN-1; i++)
      { x = Map[(int) (B[i])];
        if (x >= 0)
          c = (c << 2) | x;
        else
          { c = 0; h = i+1; }
      }

    b  = B + (KMERLEN-1);
    id = -TRAPMIN;
    ia =  Alen;
    ib = -MAXDIFFS;
    for (i = 0; i <= Blen-KMERLEN; i++, id++, ia++)
      { x = Map[(int) (b[i])];
        if (x >= 0)
          c = ((c << 2) | x) & Kmask;
        else
          { c = 0; h = i + KMERLEN; }

        if (i >= h)
          { int j, m;

            j = Table[c];
            m = Table[c+1];
            if (m > j + TooFrequent)
              m = j + TooFrequent;
#ifdef REPORT_SIZES
            sum += m-j;
#endif
#ifdef DEBUG_FILTER
            printf("%11d: %11d: %11d, %11d\n",i,c,j,m);
#endif

#ifdef INNER_LOOP

            for (; j < m; j++)
              { DiagRecord *dp;
                int         k, e;

                e  = ia - Tuples[j];
                if (e <= Divpt) continue;

                k  = e >> BINSHIFT;
                dp = DiagVec + (k % bintop);
#ifdef DEBUG_FILTER
                printf("   %10d -> %10d -> %9d [%3d,%3d]\n",
                       e-Alen,k,k%bintop,dp->maxim,dp->count);
#endif
                if (dp->maxim < id)
                  { if (dp->count >= HITMIN)
                      { HitTube *hp;
                        if (hits >= HitMax)
                          { HitMax = (1.1*Blen*hits)/i + 5000;
                            HitList = ckrealloc(HitList,
                                                sizeof(HitTube)*HitMax,"Hit list");
                          }
                        hp = HitList + hits;
                        hp->diag     = Alen - (k << BINSHIFT);
                        hp->bstart   = dp->minim;
                        hp->bfinish  = dp->maxim + KMERLEN;
                        hits += 1;
#ifdef DEBUG_HITS
                        printf("   %10d -> %10d [%3d,%3d]\n",
                               Tuples[j]-i,hp->diag,hp->bstart,hp->bfinish);
#endif
                      }
                    dp->count = 0;
                  }
                if (dp->count == 0)
                  dp->minim = i;
                if (dp->maxim < i)
                  { dp->count += 1;
                    dp->maxim = i;
                  }

                if ((e & BINMASK) < MAXDIFFS)
                  { if (dp == DiagVec)
                      dp  = DiagVec + (bintop-1);
                    else
                      dp -= 1;
#ifdef DEBUG_FILTER
                    printf(" * %10d -> %10d -> %9d [%3d,%3d]\n",
                           e-Alen,k-1,(k-1)%bintop,dp->maxim,dp->count);
#endif
                    if (dp->maxim < id)
                      { if (dp->count >= HITMIN)
                          { HitTube *hp;
                            if (hits >= HitMax)
                              { HitMax = (1.1*Blen*hits)/i + 5000;
                                HitList = ckrealloc(HitList,
                                                  sizeof(HitTube)*HitMax,"Hit List");
                              }
                            hp = HitList + hits;
                            hp->diag     = Alen - ((k - 1) << BINSHIFT);
                            hp->bstart   = dp->minim;
                            hp->bfinish  = dp->maxim + KMERLEN;
                            hits += 1;
#ifdef DEBUG_HITS
                            printf(" * %10d -> %10d [%3d,%3d]\n",
                                   Tuples[j]-i,hp->diag,hp->bstart,hp->bfinish);
#endif
                          }
                        dp->count = 0;
                      }
                    if (dp->count == 0)
                      dp->minim = i;
                    if (dp->maxim < i)
                      { dp->count += 1;
                        dp->maxim = i;
                      }
                  }
              }

#endif /* INNER_LOOP */

          }

        if (ib++ == BINMASK)
          { DiagRecord *dp;
            int         k;

            ib = 0;
            k  = ((i-MAXDIFFS) >> BINSHIFT);
            dp = DiagVec + (k % bintop);
#ifdef DEBUG_FILTER
            printf(" + %10d -> %10d -> %9d [%3d,%3d]\n",
                   (i-MAXDIFFS)-Alen,k,k%bintop,dp->maxim,dp->count);
#endif

#ifdef INNER_LOOP

            if (dp->count >= HITMIN)
              { HitTube *hp;
                if (hits >= HitMax)
                  { HitMax = (1.1*Blen*hits)/i + 5000;
                    HitList = ckrealloc(HitList,
                                      sizeof(HitTube)*HitMax,"Hit List");
                  }
                hp = HitList + hits;
                hp->diag     = Alen - (k << BINSHIFT);
                hp->bstart   = dp->minim;
                hp->bfinish  = dp->maxim + KMERLEN;
                hits += 1;
#ifdef DEBUG_HITS
                printf(" + %10d -> %10d [%3d,%3d]\n",
                       Alen-(i-MAXDIFFS),hp->diag,hp->bstart,hp->bfinish);
#endif
              }

#endif  /*  INNER_LOOP  */
            dp->count = 0;
            dp->maxim = 0;
          }
      }

    x   = BINMASK + 1;
    ia -= 1;
    for (h = Alen + ((Blen-KMERLEN) & BINMASK); h >= 0; h -= x)
      { DiagRecord *dp;
        int         k;

        k  = (ia - h) >> BINSHIFT;
        dp = DiagVec + (k % bintop);
        if (dp->count >= HITMIN)
          { HitTube *hp;
            if (hits >= HitMax)
              { HitMax = (1.1*Blen*hits)/i + 5000;
                HitList = ckrealloc(HitList,sizeof(HitTube)*HitMax,"Hit list");
              }
            hp = HitList + hits;
            hp->diag     = Alen - (k << BINSHIFT);
            hp->bstart   = dp->minim;
            hp->bfinish  = dp->maxim + KMERLEN;
            hits += 1;
#ifdef DEBUG_HITS
            printf(" > %10d -> %10d [%3d,%3d]\n",h-(Blen-KMERLEN),hp->diag,hp->bstart,hp->bfinish);
#endif
          }
      }
  }

#ifdef REPORT_SIZES
  fprintf(stderr,"\n  %9lld %d-mers (%f%% of matrix)\n",sum,KMERLEN,(100.*sum/Alen)/Blen);
  fprintf(stderr,"  %9lld seed hits (%f%% of matrix)\n",hits,(100.*hits/Alen)/Blen);
  fflush(stderr);
#endif

#ifdef TEST_DIAG
  { int i;

    for (i = 0; i < hits; i++)
      printf("Diag %d [%d,%d]\n",
             HitList[i].diag,HitList[i].bstart,HitList[i].bfinish);
  }
#endif

  free(DiagVec);

#else
  hits = 0;
#endif

  *Hitlen = hits;
  return (HitList);
}
