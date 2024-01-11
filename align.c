/*******************************************************************************************
 *
 *  Fast alignment discovery and trace generation along with utilites for displaying alignments
 *     Based on previously unpublished ideas from 2005, subsequently refined in 2013-14.  Basic
 *     idea is to keep a dynamically selected interval of the f.r. waves from my 1986 O(nd) paper.
 *     A recent cool idea is to not record all the details of an alignment while discovering it
 *     but simply record trace points through which the optimal alignment passes every 100bp,
 *     allowing rapid recomputation of the alignment details between trace points.
 *
 *  Author :  Gene Myers
 *  First  :  June 2013
 *  Current:  June 1, 2014
 *
 ********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>

#include "DB.h"
#include "align.h"

#undef    DEBUG_PASSES     //  Show forward / backward extension termini for Local_Alignment
#undef    DEBUG_POINTS     //  Show trace points
#undef    DEBUG_WAVE       //  Show waves of Local_Alignment
#undef     SHOW_MATCH_WAVE //  For waves of Local_Alignment also show # of matches
#undef    SHOW_TRAIL       //  Show trace at the end of forward and reverse passes
#undef    SHOW_TPS         //  Show trace points as they are encountered in a wave

#undef  DEBUG_EXTEND       //  Show waves of Extend_Until_Overlap

#undef  DEBUG_ALIGN        //  Show division points of Compute_Trace
#undef  DEBUG_TRACE        //  Show trace additions for Compute_Trace
#undef  DEBUG_SCRIPT       //  Show script additions for Compute_Trace
#undef  DEBUG_AWAVE        //  Show F/R waves of Compute_Trace

#undef  SHOW_TRACE         //  Show full trace for Print_Alignment

#undef  WAVE_STATS


/****************************************************************************************\
*                                                                                        *
*  Working Storage Abstraction                                                           *
*                                                                                        *
\****************************************************************************************/

typedef struct            //  Hidden from the user, working space for each thread
  { int     vecmax;
    void   *vector;
    int     celmax;
    void   *cells;
    int     pntmax;
    void   *points;
    int     tramax;
    void   *trace;
    int     alnmax;
    void   *alnpts;
  } _Work_Data;

Work_Data *New_Work_Data()
{ _Work_Data *work;
  
  work = (_Work_Data *) Malloc(sizeof(_Work_Data),"Allocating work data block");
  if (work == NULL)
    EXIT(NULL);
  work->vecmax = 0;
  work->vector = NULL;
  work->pntmax = 0;
  work->points = NULL;
  work->tramax = 0;
  work->trace  = NULL;
  work->alnmax = 0;
  work->alnpts = NULL;
  work->celmax = 0;
  work->cells  = NULL;
  return ((Work_Data *) work);
}

static int enlarge_vector(_Work_Data *work, int newmax)
{ void *vec;
  int   max;

  max = ((int) (newmax*1.2)) + 10000;
  vec = Realloc(work->vector,max,"Enlarging DP vector");
  if (vec == NULL)
    EXIT(1);
  work->vecmax = max;
  work->vector = vec;
  return (0);
}

static int enlarge_points(_Work_Data *work, int newmax)
{ void *vec;
  int   max;

  max = ((int) (newmax*1.2)) + 10000;
  vec = Realloc(work->points,max,"Enlarging point vector");
  if (vec == NULL)
    EXIT(1);
  work->pntmax = max;
  work->points = vec;
  return (0);
}

static int enlarge_alnpts(_Work_Data *work, int newmax)
{ void *vec;
  int   max;

  max = ((int) (newmax*1.2)) + 10000;
  vec = Realloc(work->alnpts,max,"Enlarging point vector");
  if (vec == NULL)
    EXIT(1);
  work->alnmax = max;
  work->alnpts = vec;
  return (0);
}

static int enlarge_trace(_Work_Data *work, int newmax)
{ void *vec;
  int   max;

  max = ((int) (newmax*1.2)) + 10000;
  vec = Realloc(work->trace,max,"Enlarging trace vector");
  if (vec == NULL)
    EXIT(1);
  work->tramax = max;
  work->trace  = vec;
  return (0);
}

void Free_Work_Data(Work_Data *ework)
{ _Work_Data *work = (_Work_Data *) ework;
  if (work->vector != NULL)
    free(work->vector);
  if (work->cells != NULL)
    free(work->cells);
  if (work->trace != NULL)
    free(work->trace);
  if (work->points != NULL)
    free(work->points);
  if (work->alnpts != NULL)
    free(work->alnpts);
  free(work);
}


/****************************************************************************************\
*                                                                                        *
*  ADAPTIVE PATH FINDING                                                                 *
*                                                                                        *
\****************************************************************************************/

  //  Absolute/Fixed Parameters

#define BVEC  uint64     //  Can be uint32 if PATH_LEN <= 32

#define TRIM_LEN    15   //  Report as the tip, the last wave maximum for which the last
                         //     2*TRIM_LEN edits are prefix-positive at rate ave_corr*f(bias)
                         //     (max value is 20)

#define DUB_TRIM    45   //  = 3*TRIM_LEN

#define PATH_LEN    60   //  Follow the last PATH_LEN columns/edges (max value is 63)

  //  Derivative fixed parameters

#define PATH_TOP  0x1000000000000000ll   //  Must be 1 << PATH_LEN
#define PATH_INT  0x0fffffffffffffffll   //  Must be PATH_TOP-1
#define TRIM_MASK 0x7fff                 //  Must be (1 << TRIM_LEN) - 1
#define TRIM_MLAG 250                    //  How far can last trim point be behind best point
#define WAVE_LAG   60                    //  How far can worst point be behind the best point

static double Bias_Factor[10] = { .690, .690, .690, .690, .780,
                                  .850, .900, .933, .966, 1.000 };

  //  Adjustable paramters

typedef struct
  { double ave_corr;
    int    trace_space;
    int    reach;
    float  freq[4];
    int    ave_path;
    int16 *score;
    int16 *table;
  } _Align_Spec;
 
/* Fill in bit table: TABLE[x] = 1 iff the alignment modeled by x (1 = match, 0 = mismatch)
     has a non-negative score for every suffix of the alignment under the scoring scheme
     where match = MATCH and mismatch = -1.  MATCH is set so that an alignment with TRIM_PCT
     matches has zero score ( (1-TRIM_PCT) / TRIM_PCT ).                                     */

#define FRACTION 1000  //  Implicit fractional part of scores, i.e. score = x/FRACTION

typedef struct
  { int    mscore;
    int    dscore;
    int16 *table;
    int16 *score;
  } Table_Bits;

static void set_table(int bit, int prefix, int score, int max, Table_Bits *parms)
{ if (bit >= TRIM_LEN)
    { parms->table[prefix] = (int16) (score-max);
      parms->score[prefix] = (int16) score;
    }
  else
    { if (score > max)
        max = score;
      set_table(bit+1,(prefix<<1),score - parms->dscore,max,parms);
      set_table(bit+1,(prefix<<1) | 1,score + parms->mscore,max,parms);
    }
}

/* Create an alignment specification record including path tip tables & values */

Align_Spec *New_Align_Spec(double ave_corr, int trace_space, float *freq, int reach)
{ _Align_Spec *spec;
  Table_Bits   parms;
  double       match;
  int          bias;

  spec = (_Align_Spec *) Malloc(sizeof(_Align_Spec),"Allocating alignment specification");
  if (spec == NULL)
    EXIT(NULL);

  spec->ave_corr    = ave_corr;
  spec->trace_space = trace_space;
  spec->reach       = reach;
  spec->freq[0]     = freq[0];
  spec->freq[1]     = freq[1];
  spec->freq[2]     = freq[2];
  spec->freq[3]     = freq[3];

  match = freq[0] + freq[3];
  if ((match <= 0.) == (match > 0.))   //  frequencies accidentally undefined?
    match = .5;
  if (match > .5)
    match = 1.-match;
  bias = (int) ((match+.025)*20.-1.);
  if (match < .2)
    { fprintf(stderr,"Warning: Base bias worse than 80/20%% ! (New_Align_Spec)\n");
      fprintf(stderr,"         Capping bias at this ratio.\n");
      bias = 3; 
    }

  spec->ave_path = (int) (PATH_LEN * (1. - Bias_Factor[bias] * (1. - ave_corr)));
  parms.mscore   = (int) (FRACTION * Bias_Factor[bias] * (1. - ave_corr));
  parms.dscore   = FRACTION - parms.mscore;

  parms.score = (int16 *) Malloc(sizeof(int16)*(TRIM_MASK+1)*2,"Allocating trim table");
  if (parms.score == NULL)
    { free(spec);
      EXIT(NULL);
    }
  parms.table = parms.score + (TRIM_MASK+1);

  set_table(0,0,0,0,&parms);

  spec->table = parms.table;
  spec->score = parms.score;

  return ((Align_Spec *) spec);
}

void Free_Align_Spec(Align_Spec *espec)
{ _Align_Spec *spec = (_Align_Spec *) espec;
  free(spec->score);
  free(spec);
}

double Average_Correlation(Align_Spec *espec)
{ return (((_Align_Spec *) espec)->ave_corr); }

int Trace_Spacing(Align_Spec *espec)
{ return (((_Align_Spec *) espec)->trace_space); }

float *Base_Frequencies(Align_Spec *espec)
{ return (((_Align_Spec *) espec)->freq); }

int Overlap_If_Possible(Align_Spec *espec)
{ return (((_Align_Spec *) espec)->reach); }


/****************************************************************************************\
*                                                                                        *
*  LOCAL ALIGNMENT FINDER: forward_/reverse_wave and Local_Alignment                     *
*                                                                                        *
\****************************************************************************************/


#ifdef WAVE_STATS

static int64 MAX, TOT, NWV;
static int64 RESTARTS;

void Init_Stats()
{ MAX = TOT = NWV = 0;
  RESTARTS = 0;
}

void Print_Stats()
{ printf("\nMax = %lld  Ave = %.1f  # = %lld\n",MAX,(1.*TOT)/NWV,NWV);
  printf("\nRestarts = %lld\n",RESTARTS);
}

#endif


#ifdef DEBUG_WAVE

static void print_wave(int *V, int *M, int low, int hgh, int besta)
{ int k, bestk;

  (void) M;
  printf("  [%6d,%6d]: ",low,hgh);
  for (k = low; k <= hgh; k++)
    { if (besta == V[k])
        bestk = k;
      // printf(" %3d",(V[k]+k)/2);
      printf(" %3d",besta-V[k]);
    }
  printf(" : %d (%d,%d)\n",besta,(besta+bestk)/2,(besta-bestk)/2);
#ifdef SHOW_MATCH_WAVE
  printf("                   ");
  for (k = low; k <= hgh; k++)
    printf(" %3d",M[k]);
  printf("\n");
#endif
  fflush(stdout);
}

#endif

/* At each furthest reaching point, keep a-coordinate of point (V), bitvector
     recording the last TRIM_LEN columns of the implied alignment (T), and the
     # of matches (1-bits) in the bitvector (M).                               */

typedef struct
  { int ptr;
    int diag;
    int diff;
    int mark;
  } Pebble;

static int VectorEl = 6*sizeof(int) + sizeof(BVEC);

static int forward_wave(_Work_Data *work, _Align_Spec *spec, Alignment *align, Path *bpath,
                        int *mind, int maxd, int mida, int minp, int maxp, int aoff, int boff)
{ char *aseq  = align->aseq;
  char *bseq  = align->bseq;
  Path *apath = align->path;

  int     hgh, low, dif;
  int     vlen, vmin, vmax;
  int    *V, *M;
  int    *_V, *_M;
  BVEC   *T;
  BVEC   *_T;

  int    *HA, *HB;
  int    *_HA, *_HB;
  int    *NA, *NB;
  int    *_NA, *_NB;
  Pebble *cells;
  int     avail, cmax;

  int     TRACE_SPACE = spec->trace_space;
  int     PATH_AVE    = spec->ave_path;
  int     REACH       = spec->reach;
  int16  *SCORE       = spec->score;
  int16  *TABLE       = spec->table;

  int     besta, besty;
  int     trima, trimy, trimd;
  int     trimha, trimhb;
  int     morea, morey, mored;
  int     moreha, morehb;
  int     more, morem, lasta;
  int     aclip, bclip;

  hgh = maxd;
  low = *mind;
  dif = 0;

  { int span, wing;

    span = (hgh-low)+1;
    vlen = work->vecmax/VectorEl;
    wing = (vlen - span)/2;
    vmin = low - wing;
    vmax = hgh + wing;

    _V  = ((int *) work->vector);
    _M  = _V + vlen;
    _HA = _M + vlen;
    _HB = _HA + vlen;
    _NA = _HB + vlen;
    _NB = _NA + vlen;
    _T  = ((BVEC *) (_NB + vlen));

    V  = _V-vmin;
    M  = _M-vmin;
    HA = _HA-vmin;
    HB = _HB-vmin;
    NA = _NA-vmin;
    NB = _NB-vmin;
    T  = _T-vmin;

    cells = (Pebble *) (work->cells);
    cmax  = work->celmax;
    avail = 0;
  }

  /* Compute 0-wave starting from mid-line */

  more  = 1;
  aclip =  INT32_MAX;
  bclip = -INT32_MAX;

  besta  = trima  = morea = lasta = mida;
  besty  = trimy  = morey = (mida-hgh) >> 1;
  trimd  = mored  = 0;
  trimha = moreha = 0;
  trimhb = morehb = 1;
  morem  = -1;

  { int   k;
    char *a;

    a  = aseq + hgh;
    for (k = hgh; k >= low; k--)
      { int     y, c, d;
        int     ha, hb;
        int     na, nb;
        Pebble *pb;

        y = (mida-k) >> 1;

        if (avail >= cmax-1)
          { cmax  = ((int) (avail*1.2)) + 10000;
            cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
            if (cells == NULL)
              EXIT(1);
            work->celmax = cmax;
            work->cells  = (void *) cells;
          }

        na = (((y+k)+(TRACE_SPACE-aoff))/TRACE_SPACE-1)*TRACE_SPACE+aoff;
#ifdef SHOW_TPS
        printf(" A %d: %d,%d,0,%d\n",avail,-1,k,na); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
        pb->mark = na;
        ha  = avail++;
        na += TRACE_SPACE;

        nb = ((y+(TRACE_SPACE-boff))/TRACE_SPACE-1)*TRACE_SPACE+boff;
#ifdef SHOW_TPS
        printf(" B %d: %d,%d,0,%d\n",avail,-1,k,nb); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
        pb->mark = nb;
        hb  = avail++;
        nb += TRACE_SPACE;

        while (1)
          { c = bseq[y];
            if (c == 4)
              { more = 0;
                if (bclip < k)
                  bclip = k;
                break;
              }
            d = a[y];
            if (c != d)
              { if (d == 4) 
                  { more  = 0;
                    aclip = k;
                  }
                break;
              }
            y += 1;
          }
        c = (y << 1) + k;

        while (y+k >= na)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" A %d: %d,%d,0,%d\n",avail,ha,k,na); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = ha;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = na;
            ha  = avail++;
            na += TRACE_SPACE;
          }
        while (y >= nb)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" B %d: %d,%d,0,%d\n",avail,hb,k,nb); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = hb;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = nb;
            hb  = avail++;
            nb += TRACE_SPACE;
          }

        if (c > besta)
          { besta  = trima = lasta = c;
            besty  = trimy = y;
            trimha = ha;
            trimhb = hb;
          }

        V[k]  = c;
        T[k]  = PATH_INT;
        M[k]  = PATH_LEN;
        HA[k] = ha;
        HB[k] = hb;
        NA[k] = na;
        NB[k] = nb;

        a -= 1;
      }
  }

  if (more == 0)
    { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
        more = 1;
      if (hgh >= aclip)
        { hgh = aclip-1;
          if (morem <= M[aclip])
            { morem  = M[aclip];
              morea  = V[aclip];
              morey  = (morea - aclip)/2;
              moreha = HA[aclip];
              morehb = HB[aclip];
            }
        }
      if (low <= bclip)
        { low = bclip+1;
          if (morem <= M[bclip])
            { morem  = M[bclip];
              morea  = V[bclip];
              morey  = (morea - bclip)/2;
              moreha = HA[bclip];
              morehb = HB[bclip];
            }
        }
      aclip =  INT32_MAX;
      bclip = -INT32_MAX;
    }

#ifdef DEBUG_WAVE
  printf("\nFORWARD WAVE:\n");
  print_wave(V,M,low,hgh,besta);
#endif

  /* Compute successive waves until no furthest reaching points remain */

  while (more && lasta >= besta - TRIM_MLAG)
    { int     k, n;
      int     ua, ub;
      BVEC    t;
      int     am, ac, ap;
      char   *a;

      low -= 1;
      hgh += 1;

      if (low <= vmin || hgh >= vmax)
        { int   span, wing;
          int64 move;
          int64 vd, md, had, hbd, nad, nbd, td;

          span = (hgh-low)+1;
          if (.8*vlen < span)
            { if (enlarge_vector(work,vlen*VectorEl))
                EXIT(1);

              move = ((void *) _V) - work->vector;
              vlen = work->vecmax/VectorEl;

              _V  = (int *) work->vector;
              _M  = _V + vlen;
              _HA = _M + vlen;
              _HB = _HA + vlen;
              _NA = _HB + vlen;
              _NB = _NA + vlen;
              _T  = ((BVEC *) (_NB + vlen));
            }
          else
            move = 0;

          wing = (vlen - span)/2;

          vd  = ((void *) ( _V+wing)) - (((void *) ( V+low)) - move);
          md  = ((void *) ( _M+wing)) - (((void *) ( M+low)) - move);
          had = ((void *) (_HA+wing)) - (((void *) (HA+low)) - move);
          hbd = ((void *) (_HB+wing)) - (((void *) (HB+low)) - move);
          nad = ((void *) (_NA+wing)) - (((void *) (NA+low)) - move);
          nbd = ((void *) (_NB+wing)) - (((void *) (NB+low)) - move);
          td  = ((void *) ( _T+wing)) - (((void *) ( T+low)) - move);

          if (vd < 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));
          if (md < 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (had < 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (hbd < 0)
            memmove(_HB+wing,  ((void *) (HB+low)) - move, span*sizeof(int));
          if (nad < 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (nbd < 0)
            memmove(_NB+wing,  ((void *) (NB+low)) - move, span*sizeof(int));
          if (td < 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));

          if (td > 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));
          if (nbd > 0)
            memmove(_NB+wing,  ((void *) (NB+low)) - move, span*sizeof(int));
          if (nad > 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (hbd > 0)
            memmove(_HB+wing,  ((void *) (HB+low)) - move, span*sizeof(int));
          if (had > 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (md > 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (vd > 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));

          vmin = low-wing;
          vmax = hgh+wing;

          V  =  _V-vmin;
          M  =  _M-vmin;
          HA = _HA-vmin;
          HB = _HB-vmin;
          NA = _NA-vmin;
          NB = _NB-vmin;
          T  =  _T-vmin;
        }

      if (low >= minp)
        { NA[low] = NA[low+1];
          NB[low] = NB[low+1];
          V[low]  = -1;
        }
      else
        low += 1;

      if (hgh <= maxp)
        { NA[hgh] = NA[hgh-1];
          NB[hgh] = NB[hgh-1];
          V[hgh]  = am = -1;
        }
      else
        am = V[--hgh];

      dif += 1;

      ac = V[hgh+1] = V[low-1] = -1;
      a  = aseq + hgh;
      t  = PATH_INT;
      n  = PATH_LEN;
      ua = ub = -1;
      for (k = hgh; k >= low; k--)
        { int     y, m;
          int     ha, hb;
          int     c, d;
          BVEC    b;
          Pebble *pb;

          ap = ac;
          ac = am;
          am = V[d = k-1];

          if (ac < am)
            if (am < ap)
              { c  = ap+1;
                m  = n;
                b  = t;
                ha = ua;
                hb = ub;
              }
            else
              { c  = am+1;
                m  = M[d];
                b  = T[d]; 
                ha = HA[d];
                hb = HB[d];
              }
          else
            if (ac < ap)
              { c  = ap+1;
                m  = n;
                b  = t;
                ha = ua;
                hb = ub;
              }
            else
              { c  = ac+2;
                m  = M[k];
                b  = T[k];
                ha = HA[k];
                hb = HB[k];
              }

          if ((b & PATH_TOP) != 0)
            m -= 1;
          b <<= 1;

          y = (c-k) >> 1;
          while (1)
            { c = bseq[y];
              if (c == 4)
                { more = 0;
                  if (bclip < k)
                    bclip = k;
                  break;
                }
              d = a[y];
              if (c != d)
                { if (d == 4) 
                    { more  = 0;
                      aclip = k;
                    }
                  break;
                }
              y += 1;
              if ((b & PATH_TOP) == 0)
                m += 1;
              b = (b << 1) | 1;
            }
          c = (y << 1) + k;

          while (y+k >= NA[k])
            { if (cells[ha].mark < NA[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" A %d: %d,%d,%d,%d\n",avail,ha,k,dif,NA[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = ha;
                  pb->diag = k;
                  pb->diff = dif;
                  pb->mark = NA[k];
                  ha = avail++;
                }
              NA[k] += TRACE_SPACE;
            }

          while (y >= NB[k])
            { if (cells[hb].mark < NB[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" B %d: %d,%d,%d,%d\n",avail,hb,k,dif,NB[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = hb;
                  pb->diag = k;
                  pb->diff = dif;
                  pb->mark = NB[k];
                  hb = avail++;
                }
              NB[k] += TRACE_SPACE;
            }

          if (c > besta)
            { besta = c;
              besty = y;
              if (m >= PATH_AVE)
                { lasta = c;
                  if (TABLE[b & TRIM_MASK] >= 0)
                    if (TABLE[(b >> TRIM_LEN) & TRIM_MASK] + SCORE[b & TRIM_MASK] >= 0)
                      { trima  = c;
                        trimy  = y;
                        trimd  = dif;
                        trimha = ha;
                        trimhb = hb;
                      }
                }
            }

          t  = T[k];
          n  = M[k];
          ua = HA[k];
          ub = HB[k];
          V[k]  = c;
          T[k]  = b;
          M[k]  = m;
          HA[k] = ha;
          HB[k] = hb;

          a -= 1;
        }

      if (more == 0)
        { if (bseq[besty] != 4 && aseq[besta-besty] != 4)
            more = 1;
          if (hgh >= aclip)
            { hgh = aclip-1;
              if (morem <= M[aclip])
                { morem  = M[aclip];
                  morea  = V[aclip];
                  morey  = (morea - aclip)/2;
                  mored  = dif;
                  moreha = HA[aclip];
                  morehb = HB[aclip];
                }
            }
          if (low <= bclip)
            { low = bclip+1;
              if (morem <= M[bclip])
                { morem  = M[bclip];
                  morea  = V[bclip];
                  morey  = (morea - bclip)/2;
                  mored  = dif;
                  moreha = HA[bclip];
                  morehb = HB[bclip];
                }
            }
          aclip =  INT32_MAX;
          bclip = -INT32_MAX;
        }

      n = besta - WAVE_LAG;
      while (hgh >= low)
        if (V[hgh] < n)
          hgh -= 1;                               
        else
          { while (V[low] < n)
              low += 1;
            break;
          }

#ifdef WAVE_STATS
      k = (hgh-low)+1;
      if (k > MAX)
        MAX = k;
      TOT += k;
      NWV += 1;
#endif

#ifdef DEBUG_WAVE
      print_wave(V,M,low,hgh,besta);
#endif
    }

  { uint16 *atrace = (uint16 *) apath->trace;
    uint16 *btrace = (uint16 *) bpath->trace;
    int     atlen, btlen;
    int     trimx;
    int     a, b, k, h;
    int     d, e;

    if (morem >= 0 && REACH)
      { trimx  = morea-morey;
        trimy  = morey;
        trimd  = mored;
        trimha = moreha;
        trimhb = morehb;
      }
    else
      trimx = trima-trimy;

    atlen = btlen = 0;

    a = -1;
    for (h = trimha; h >= 0; h = b)
      { b = cells[h].ptr; 
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = (mida-k)/2;
    e = 0;
#ifdef SHOW_TRAIL
    printf("  A path = (%5d,%5d)\n",(mida+k)/2,b); fflush(stdout);
#endif
    for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
      { k = cells[h].diag;
        a = cells[h].mark - k;
        d = cells[h].diff;
        atrace[atlen++] = (uint16) (d-e);
        atrace[atlen++] = (uint16) (a-b);
#ifdef SHOW_TRAIL
        printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,a-b); fflush(stdout);
#endif
        b = a;
        e = d;
      }
    if (b+k != trimx)
      { atrace[atlen++] = (uint16) (trimd-e);
        atrace[atlen++] = (uint16) (trimy-b);
#ifdef SHOW_TRAIL
        printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimy-b); fflush(stdout);
#endif
      }
    else if (b != trimy)
      { atrace[atlen-1] = (uint16) (atrace[atlen-1] + (trimy-b));
        atrace[atlen-2] = (uint16) (atrace[atlen-2] + (trimd-e));
#ifdef SHOW_TRAIL
        printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimy-b); fflush(stdout);
#endif
      }

    a = -1;
    for (h = trimhb; h >= 0; h = b)
      { b = cells[h].ptr; 
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = (mida+k)/2;
    e = 0;
    low = k;
#ifdef SHOW_TRAIL
    printf("  B path = (%5d,%5d)\n",b,(mida-k)/2); fflush(stdout);
#endif
    for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
      { k = cells[h].diag;
        a = cells[h].mark + k;
        d = cells[h].diff;
        btrace[btlen++] = (uint16) (d-e);
        btrace[btlen++] = (uint16) (a-b);  
#ifdef SHOW_TRAIL
        printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a,a-k,d-e,a-b); fflush(stdout);
#endif
        b = a;
        e = d;
      }
    if (b-k != trimy)
      { btrace[btlen++] = (uint16) (trimd-e);
        btrace[btlen++] = (uint16) (trimx-b);  
#ifdef SHOW_TRAIL
        printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimx-b); fflush(stdout);
#endif
      }
    else if (b != trimx)
      { btrace[btlen-1] = (uint16) (btrace[btlen-1] + (trimx-b));
        btrace[btlen-2] = (uint16) (btrace[btlen-2] + (trimd-e));
#ifdef SHOW_TRAIL
        printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimx-b); fflush(stdout);
#endif
      }

    apath->aepos = trimx;
    apath->bepos = trimy;
    apath->diffs = trimd;
    apath->tlen  = atlen;
    bpath->tlen  = btlen;
  }

  *mind = low;
  return (0);
}

/*** Reverse Wave ***/

static int reverse_wave(_Work_Data *work, _Align_Spec *spec, Alignment *align, Path *bpath,
                        int mind, int maxd, int mida, int minp, int maxp, int aoff, int boff)
{ char *aseq  = align->aseq - 1;
  char *bseq  = align->bseq - 1;
  Path *apath = align->path;

  int     hgh, low, dif;
  int     vlen, vmin, vmax;
  int    *V, *M;
  int    *_V, *_M;
  BVEC   *T;
  BVEC   *_T;

  int    *HA, *HB;
  int    *_HA, *_HB;
  int    *NA, *NB;
  int    *_NA, *_NB;
  Pebble *cells;
  int     avail, cmax;

  int     TRACE_SPACE = spec->trace_space;
  int     PATH_AVE    = spec->ave_path;
  int     REACH       = spec->reach;
  int16  *SCORE       = spec->score;
  int16  *TABLE       = spec->table;

  int     besta, besty;
  int     trima, trimy, trimd;
  int     trimha, trimhb;
  int     morea, morey, mored;
  int     moreha, morehb;
  int     more, morem, lasta;
  int     aclip, bclip;

  hgh = maxd;
  low = mind;
  dif = 0;

  { int span, wing;

    span = (hgh-low)+1;
    vlen = work->vecmax/VectorEl;
    wing = (vlen - span)/2;
    vmin = low - wing;
    vmax = hgh + wing;

    _V  = ((int *) work->vector);
    _M  = _V + vlen;
    _HA = _M + vlen;
    _HB = _HA + vlen;
    _NA = _HB + vlen;
    _NB = _NA + vlen;
    _T  = ((BVEC *) (_NB + vlen));

    V  = _V-vmin;
    M  = _M-vmin;
    HA = _HA-vmin;
    HB = _HB-vmin;
    NA = _NA-vmin;
    NB = _NB-vmin;
    T  = _T-vmin;

    cells = (Pebble *) (work->cells);
    cmax  = work->celmax;
    avail = 0;
  }

  more  = 1;
  aclip = -INT32_MAX;
  bclip =  INT32_MAX;

  besta  = trima  = morea = lasta = mida;
  besty  = trimy  = morey = (mida-hgh) >> 1;
  trimd  = mored  = 0;
  trimha = moreha = 0;
  trimhb = morehb = 1;
  morem  = -1;

  { int   k;
    char *a;

    a = aseq + low;
    for (k = low; k <= hgh; k++)
      { int     y, c, d;
        int     ha, hb;
        int     na, nb;
        Pebble *pb;

        y = (mida-k) >> 1;

        if (avail >= cmax-1)
          { cmax  = ((int) (avail*1.2)) + 10000;
            cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
            if (cells == NULL)
              EXIT(1);
            work->celmax = cmax;
            work->cells  = (void *) cells;
          }

        na = (((y+k)+(TRACE_SPACE-aoff)-1)/TRACE_SPACE-1)*TRACE_SPACE+aoff;
#ifdef SHOW_TPS
        printf(" A %d: -1,%d,0,%d\n",avail,k,na+TRACE_SPACE); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
        pb->mark = y+k;
        ha  = avail++;

        nb = ((y+(TRACE_SPACE-boff)-1)/TRACE_SPACE-1)*TRACE_SPACE+boff;
#ifdef SHOW_TPS
        printf(" B %d: -1,%d,0,%d\n",avail,k,nb+TRACE_SPACE); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
        pb->mark = y;
        hb  = avail++;

        while (1)
          { c = bseq[y];
            if (c == 4)
              { more = 0;
                if (bclip > k)
                  bclip = k;
                break;
              }
            d = a[y];
            if (c != d)
              { if (d == 4) 
                  { more  = 0;
                    aclip = k;
                  }
                break;
              }
            y -= 1;
          }
        c = (y << 1) + k;

        while (y+k <= na)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" A %d: %d,%d,0,%d\n",avail,ha,k,na); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = ha;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = na;
            ha  = avail++;
            na -= TRACE_SPACE;
          }
        while (y <= nb)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" B %d: %d,%d,0,%d\n",avail,hb,k,nb); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = hb;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = nb;
            hb  = avail++;
            nb -= TRACE_SPACE;
          }

        if (c < besta)
          { besta  = trima = lasta = c;
            besty  = trimy = y;
            trimha = ha;
            trimhb = hb;
          }

        V[k]  = c;
        T[k]  = PATH_INT;
        M[k]  = PATH_LEN;
        HA[k] = ha;
        HB[k] = hb;
        NA[k] = na;
        NB[k] = nb;

        a += 1;
      }
  }

  if (more == 0)
    { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
        more = 1;
      if (low <= aclip)
        { low = aclip+1;
          if (morem <= M[aclip])
            { morem  = M[aclip];
              morea  = V[aclip];
              morey  = (morea - aclip)/2;
              moreha = HA[aclip];
              morehb = HB[aclip];
            }
        }
      if (hgh >= bclip)
        { hgh = bclip-1;
          if (morem <= M[bclip])
            { morem  = M[bclip];
              morea  = V[bclip];
              morey  = (morea - bclip)/2;
              moreha = HA[bclip];
              morehb = HB[bclip];
            }
        }
      aclip = -INT32_MAX;
      bclip =  INT32_MAX;
    }

#ifdef DEBUG_WAVE
  printf("\nREVERSE WAVE:\n");
  print_wave(V,M,low,hgh,besta);
#endif

  while (more && lasta <= besta + TRIM_MLAG)
    { int    k, n;
      int    ua, ub;
      BVEC   t;
      int    am, ac, ap;
      char  *a;

      low -= 1;
      hgh += 1;

      if (low <= vmin || hgh >= vmax)
        { int   span, wing;
          int64 move, vd, md, had, hbd, nad, nbd, td;

          span = (hgh-low)+1;
          if (.8*vlen < span)
            { if (enlarge_vector(work,vlen*VectorEl))
                EXIT(1);

              move = ((void *) _V) - work->vector;
              vlen = work->vecmax/VectorEl;

              _V  = (int *) work->vector;
              _M  = _V + vlen;
              _HA = _M + vlen;
              _HB = _HA + vlen;
              _NA = _HB + vlen;
              _NB = _NA + vlen;
              _T  = ((BVEC *) (_NB + vlen));
            }
          else
            move = 0;

          wing = (vlen - span)/2;

          vd  = ((void *) ( _V+wing)) - (((void *) ( V+low)) - move);
          md  = ((void *) ( _M+wing)) - (((void *) ( M+low)) - move);
          had = ((void *) (_HA+wing)) - (((void *) (HA+low)) - move);
          hbd = ((void *) (_HB+wing)) - (((void *) (HB+low)) - move);
          nad = ((void *) (_NA+wing)) - (((void *) (NA+low)) - move);
          nbd = ((void *) (_NB+wing)) - (((void *) (NB+low)) - move);
          td  = ((void *) ( _T+wing)) - (((void *) ( T+low)) - move);

          if (vd < 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));
          if (md < 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (had < 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (hbd < 0)
            memmove(_HB+wing,  ((void *) (HB+low)) - move, span*sizeof(int));
          if (nad < 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (nbd < 0)
            memmove(_NB+wing,  ((void *) (NB+low)) - move, span*sizeof(int));
          if (td < 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));

          if (td > 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));
          if (nbd > 0)
            memmove(_NB+wing,  ((void *) (NB+low)) - move, span*sizeof(int));
          if (nad > 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (hbd > 0)
            memmove(_HB+wing,  ((void *) (HB+low)) - move, span*sizeof(int));
          if (had > 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (md > 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (vd > 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));

          vmin = low-wing;
          vmax = hgh+wing;

          V  =  _V-vmin;
          M  =  _M-vmin;
          HA = _HA-vmin;
          HB = _HB-vmin;
          NA = _NA-vmin;
          NB = _NB-vmin;
          T  =  _T-vmin;
        }

      if (low >= minp)
        { NA[low] = NA[low+1];
          NB[low] = NB[low+1];
          V[low]  = ap = INT32_MAX;
        }
      else
        ap = V[++low]; 

      if (hgh <= maxp)
        { NA[hgh] = NA[hgh-1];
          NB[hgh] = NB[hgh-1];
          V[hgh] = INT32_MAX;
        }
      else
        hgh -= 1;

      dif += 1;

      ac = V[hgh+1] = V[low-1] = INT32_MAX;
      a  = aseq + low;
      t  = PATH_INT;
      n  = PATH_LEN;
      ua = ub = -1;
      for (k = low; k <= hgh; k++)
        { int     y, m;
          int     ha, hb;
          int     c, d;
          BVEC    b;
          Pebble *pb;

          am = ac;
          ac = ap;
          ap = V[d = k+1];

          if (ac > ap)
            if (ap > am)
              { c = am-1;
                m  = n;
                b  = t;
                ha = ua;
                hb = ub;
              }
            else
              { c  = ap-1;
                m  = M[d];
                b  = T[d];
                ha = HA[d];
                hb = HB[d];
              }
          else
            if (ac > am)
              { c  = am-1;
                m  = n;
                b  = t;
                ha = ua;
                hb = ub;
              }
            else
              { c  = ac-2;
                m  = M[k];
                b  = T[k];
                ha = HA[k];
                hb = HB[k];
              }

          if ((b & PATH_TOP) != 0)
            m -= 1;
          b <<= 1;

          y = (c-k) >> 1;
          while (1)
            { c = bseq[y];
              if (c == 4)
                { more = 0;
                  if (bclip > k)
                    bclip = k;
                  break;
                }
              d = a[y];
              if (c != d)
                { if (d == 4) 
                    { more  = 0;
                      aclip = k;
                    }
                  break;
                }
              y -= 1;
              if ((b & PATH_TOP) == 0)
                m += 1;
              b = (b << 1) | 1;
            }
          c = (y << 1) + k;

          while (y+k <= NA[k])
            { if (cells[ha].mark > NA[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" A %d: %d,%d,%d,%d\n",avail,ha,k,dif,NA[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = ha;
                  pb->diag = k;
                  pb->diff = dif;
                  pb->mark = NA[k];
                  ha = avail++;
                }
              NA[k] -= TRACE_SPACE;
            }
          while (y <= NB[k])
            { if (cells[hb].mark > NB[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" B %d: %d,%d,%d,%d\n",avail,hb,k,dif,NB[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = hb;
                  pb->diag = k;
                  pb->diff = dif;
                  pb->mark = NB[k];
                  hb = avail++;
                }
              NB[k] -= TRACE_SPACE;
            }

          if (c < besta)
            { besta = c;
              besty = y;
              if (m >= PATH_AVE)
                { lasta = c;
                  if (TABLE[b & TRIM_MASK] >= 0)
                    if (TABLE[(b >> TRIM_LEN) & TRIM_MASK] + SCORE[b & TRIM_MASK] >= 0)
                      { trima  = c;
                        trimy  = y;
                        trimd  = dif;
                        trimha = ha;
                        trimhb = hb;
                      }
                }
            }

          t  = T[k];
          n  = M[k];
          ua = HA[k];
          ub = HB[k];
          V[k]  = c;
          T[k]  = b;
          M[k]  = m;
          HA[k] = ha;
          HB[k] = hb;

          a += 1;
        }

      if (more == 0)
        { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
            more = 1;
          if (low <= aclip)
            { low = aclip+1;
              if (morem <= M[aclip])
                { morem  = M[aclip];
                  morea  = V[aclip];
                  morey  = (morea - aclip)/2;
                  mored  = dif;
                  moreha = HA[aclip];
                  morehb = HB[aclip];
                }
            }
          if (hgh >= bclip)
            { hgh = bclip-1;
              if (morem <= M[bclip])
                { morem  = M[bclip];
                  morea  = V[bclip];
                  morey  = (morea - bclip)/2;
                  mored  = dif;
                  moreha = HA[bclip];
                  morehb = HB[bclip];
                }
            }
          aclip = -INT32_MAX;
          bclip =  INT32_MAX;
        }

      n = besta + WAVE_LAG;
      while (hgh >= low)
        if (V[hgh] > n)
          hgh -= 1;                               
        else
          { while (V[low] > n)
              low += 1;
            break;
          }

#ifdef WAVE_STATS
      k = (hgh-low)+1;
      if (k > MAX)
        MAX = k;
      TOT += k;
      NWV += 1;
#endif

#ifdef DEBUG_WAVE
      print_wave(V,M,low,hgh,besta);
#endif
    }

  { uint16 *atrace = (uint16 *) apath->trace;
    uint16 *btrace = (uint16 *) bpath->trace;
    int     atlen, btlen;
    int     trimx;
    int     a, b, k, h;
    int     d, e;

    if (morem >= 0 && REACH)
      { trimx  = morea-morey;
        trimy  = morey;
        trimd  = mored;
        trimha = moreha;
        trimhb = morehb;
      }
    else
      trimx = trima-trimy;

    atlen = btlen = 0;

    a = -1;
    for (h = trimha; h >= 0; h = b)
      { b = cells[h].ptr; 
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = cells[h].mark - k;
    e = 0;
#ifdef SHOW_TRAIL
    printf("  A path = (%5d,%5d)\n",b+k,b); fflush(stdout);
#endif
    if ((b+k)%TRACE_SPACE != aoff)
      { h = cells[h].ptr;
        if (h < 0)
          { a = trimy;
            d = trimd;
          }
        else
          { k = cells[h].diag;
            a = cells[h].mark - k;
            d = cells[h].diff;
          }
#ifdef SHOW_TRAIL
        printf("    +%4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,b-a); fflush(stdout);
#endif
        if (apath->tlen == 0)
          { atrace[--atlen] = (uint16) (b-a);
            atrace[--atlen] = (uint16) (d-e);
          }
        else
          { atrace[1] = (uint16) (atrace[1] + (b-a));
            atrace[0] = (uint16) (atrace[0] + (d-e));
          }
        b = a;
        e = d;
      }
    if (h >= 0)
      { for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
          { k = cells[h].diag;
            a = cells[h].mark - k;
            atrace[--atlen] = (uint16) (b-a);
            d = cells[h].diff;
            atrace[--atlen] = (uint16) (d-e);
#ifdef SHOW_TRAIL
            printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,b-a); fflush(stdout);
#endif
            b = a;
            e = d;
          }
        if (b+k != trimx)
          { atrace[--atlen] = (uint16) (b-trimy);
            atrace[--atlen] = (uint16) (trimd-e);
#ifdef SHOW_TRAIL
            printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimy); fflush(stdout);
#endif
          }
        else if (b != trimy)
          { atrace[atlen+1] = (uint16) (atrace[atlen+1] + (b-trimy));
            atrace[atlen]   = (uint16) (atrace[atlen]   + (trimd-e));
#ifdef SHOW_TRAIL
            printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimy); fflush(stdout);
#endif
          }
      }

    a = -1;
    for (h = trimhb; h >= 0; h = b)
      { b = cells[h].ptr; 
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = cells[h].mark + k;
    e = 0;
#ifdef SHOW_TRAIL
    printf("  B path = (%5d,%5d)\n",b,b-k); fflush(stdout);
#endif
    if ((b-k)%TRACE_SPACE != boff)
      { h = cells[h].ptr;
        if (h < 0)
          { a = trimx;
            d = trimd;
          } 
        else
          { k = cells[h].diag;
            a = cells[h].mark + k;
            d = cells[h].diff;
          }
#ifdef SHOW_TRAIL
        printf("    +%4d: (%5d,%5d): %3d / %3d\n",h,a,a-k,d-e,b-a); fflush(stdout);
#endif
        if (bpath->tlen == 0)
          { btrace[--btlen] = (uint16) (b-a);
            btrace[--btlen] = (uint16) (b-a);
          }
        else
          { btrace[1] = (uint16) (btrace[1] + (b-a));
            btrace[0] = (uint16) (btrace[0] + (d-e));
          }
        b = a;
        e = d;
      }

    if (h >= 0)
      { for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
          { k = cells[h].diag;
            a = cells[h].mark + k;
            btrace[--btlen] = (uint16) (b-a);
            d = cells[h].diff;
            btrace[--btlen] = (uint16) (d-e);
#ifdef SHOW_TRAIL
            printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a,a-k,d-e,b-a); fflush(stdout);
#endif
            b = a;
            e = d;
          }
        if (b-k != trimy)
          { btrace[--btlen] = (uint16) (b-trimx);
            btrace[--btlen] = (uint16) (trimd-e);
#ifdef SHOW_TRAIL
            printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimx); fflush(stdout);
#endif
          }
        else if (b != trimx)
          { btrace[btlen+1] = (uint16) (btrace[btlen+1] + (b-trimx));
            btrace[btlen]   = (uint16) (btrace[btlen]   + (trimd-e));
#ifdef SHOW_TRAIL
            printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimx); fflush(stdout);
#endif
          }
      }

    apath->abpos = trimx;
    apath->bbpos = trimy;
    apath->diffs = apath->diffs + trimd;
    apath->tlen  = apath->tlen  - atlen;
    apath->trace = atrace + atlen;
    bpath->tlen  = bpath->tlen  - btlen;
    bpath->trace = btrace + btlen;
  }

  return (0);
}


/* Find the longest local alignment between aseq and bseq through (xcnt,ycnt)
   See associated .h file for the precise definition of the interface.
*/

Path *Local_Alignment(Alignment *align, Work_Data *ework, Align_Spec *espec,
                      int low, int hgh, int anti, int lbord, int hbord)
{ _Work_Data  *work = ( _Work_Data *) ework;
  _Align_Spec *spec = (_Align_Spec *) espec;

  Path *apath, *bpath;
  int   aoff, boff;
  int   minp, maxp;
  int   selfie;
  int   fshort, rshort;

  { int alen, blen;
    int maxtp, wsize;

    alen = align->alen;
    blen = align->blen;

    if (hgh-low >= 7500)
      wsize = VectorEl*(hgh-low+1);
    else
      wsize = VectorEl*10000;
    if (wsize >= work->vecmax)
      if (enlarge_vector(work,wsize))
        EXIT(NULL);

    if (alen < blen)
      maxtp = 2*(blen/spec->trace_space+2);
    else
      maxtp = 2*(alen/spec->trace_space+2);
    wsize = 4*maxtp*sizeof(uint16) + sizeof(Path);
    if (wsize > work->pntmax)
      if (enlarge_points(work,wsize))
        EXIT(NULL);

    apath = align->path;
    bpath = (Path *) work->points;

    apath->trace = ((uint16 *) (bpath+1)) + maxtp;
    bpath->trace = ((uint16 *) apath->trace) +  2*maxtp;
  }

#ifdef DEBUG_PASSES
  printf("\n");
#endif

  selfie = (align->aseq == align->bseq);

  while (((anti-hgh) >> 1) < 0)
    hgh -= 1;
   
  if (lbord < 0)
    { if (selfie && low >= 0)
        minp = 1;
      else
        minp = -INT32_MAX;
    }
  else
    minp = low-lbord;
  if (hbord < 0)
    { if (selfie && hgh <= 0)
        maxp = -1;
      else
        maxp = INT32_MAX;
    }
  else
    maxp = hgh+hbord;

  if (ACOMP(align->flags))
    { aoff = align->alen % spec->trace_space;
      boff = 0;
    }
  else if (COMP(align->flags))
    { aoff = 0;
      boff = align->blen % spec->trace_space;
    }
  else
    { aoff = 0;
      boff = 0;
    }

  if (forward_wave(work,spec,align,bpath,&low,hgh,anti,minp,maxp,aoff,boff))
    EXIT(NULL);

  fshort = ((apath->aepos + apath->bepos) - anti < DUB_TRIM);

#ifdef DEBUG_PASSES
  printf("F1 (%d,%d) ~ %d => (%d,%d) %d\n",
         (2*anti+(low+hgh))/4,(2*anti-(low+hgh))/4,hgh-low,
         apath->aepos,apath->bepos,apath->diffs);
#endif

  if (reverse_wave(work,spec,align,bpath,low,low,anti,minp,maxp,aoff,boff))
    EXIT(NULL);

#ifdef DEBUG_PASSES
  printf("R1 (%d,%d) => (%d,%d) %d\n",
         (anti+low)/2,(anti-low)/2,apath->abpos,apath->bbpos,apath->diffs);
#endif

  rshort = (anti - (apath->abpos + apath->bbpos) < DUB_TRIM);

  if (fshort)
    { if (rshort)
        { apath->aepos = apath->abpos = (apath->abpos+apath->aepos)/2;
          apath->bepos = apath->bbpos = (apath->bbpos+apath->bepos)/2;
          bpath->aepos = bpath->abpos = (bpath->abpos+bpath->aepos)/2;
          bpath->bepos = bpath->bbpos = (bpath->bbpos+bpath->bepos)/2;
          apath->tlen  = 0;
          bpath->tlen  = 0;
        }
      else
        { low  = apath->abpos - apath->bbpos;
          anti = apath->abpos + apath->bbpos;
          apath->tlen = bpath->tlen = 0;
          if (forward_wave(work,spec,align,bpath,&low,low,anti,minp,maxp,aoff,boff))
            EXIT(NULL);
        }
    }
  else
    { if (rshort)
        { low  = apath->aepos - apath->bepos;
          anti = apath->aepos + apath->bepos;
          apath->tlen = bpath->tlen = 0;
          apath->diffs = 0;
          if (reverse_wave(work,spec,align,bpath,low,low,anti,minp,maxp,aoff,boff))
            EXIT(NULL);
        }
    }
        

  bpath->diffs = apath->diffs;
  if (ACOMP(align->flags))
    { uint16 *trace = (uint16 *) apath->trace;
      uint16  p;
      int     i, j;

      bpath->aepos = apath->bepos;
      bpath->bepos = apath->aepos;
      bpath->abpos = apath->bbpos;
      bpath->bbpos = apath->abpos;

      apath->abpos = align->alen - bpath->bepos;
      apath->bbpos = align->blen - bpath->aepos;
      apath->aepos = align->alen - bpath->bbpos;
      apath->bepos = align->blen - bpath->abpos;
      i = apath->tlen-2;
      j = 0;
      while (j < i)
        { p = trace[i];
          trace[i] = trace[j];
          trace[j] = p;
          p = trace[i+1];
          trace[i+1] = trace[j+1];
          trace[j+1] = p;
          i -= 2;
          j += 2;
        }
    }
  else if (COMP(align->flags))
    { uint16 *trace = (uint16 *) bpath->trace;
      uint16  p;
      int     i, j;

      bpath->abpos = align->blen - apath->bepos;
      bpath->bbpos = align->alen - apath->aepos;
      bpath->aepos = align->blen - apath->bbpos;
      bpath->bepos = align->alen - apath->abpos;
      i = bpath->tlen-2;
      j = 0;
      while (j < i)
        { p = trace[i];
          trace[i] = trace[j];
          trace[j] = p;
          p = trace[i+1];
          trace[i+1] = trace[j+1];
          trace[j+1] = p;
          i -= 2;
          j += 2;
        }
    }
  else
    { bpath->aepos = apath->bepos;
      bpath->bepos = apath->aepos;
      bpath->abpos = apath->bbpos;
      bpath->bbpos = apath->abpos;
    }

#ifdef DEBUG_POINTS
  { uint16 *trace = (uint16 *) apath->trace;
    int     a, h;

    printf("\nA-path (%d,%d)->(%d,%d)",apath->abpos,apath->bbpos,apath->aepos,apath->bepos);
    printf(" %c\n",((COMP(align->flags) || ACOMP(align->flags)) ? 'c' : 'n'));
    a = apath->bbpos;
    for (h = 1; h < apath->tlen; h += 2)
      { int dif = trace[h-1];
        int del = trace[h];
        a += del;
        printf("      %d / %d (%d)\n",dif,del,a);
      }
  }

  { uint16 *trace = (uint16 *) bpath->trace;
    int     a, h;

    printf("\nB-path (%d,%d)->(%d,%d)",bpath->abpos,bpath->bbpos,bpath->aepos,bpath->bepos);
    printf(" %c [%d,%d]\n",((COMP(align->flags) || ACOMP(align->flags)) ? 'c' : 'n'),
                           align->blen,align->alen);
    a = bpath->bbpos;
    for (h = 1; h < bpath->tlen; h += 2)
      { int dif = trace[h-1];
        int del = trace[h];
        a += del;
        printf("      %d / %d (%d)\n",dif,del,a);
      }
  }
#endif

  return (bpath);
}


/****************************************************************************************\
*                                                                                        *
*  EXTENSION VERSION OF LOCAL ALIGNMENT                                                  *
*                                                                                        *
\****************************************************************************************/

static int VectorEn = 4*sizeof(int) + sizeof(BVEC);

static int forward_extend(_Work_Data *work, _Align_Spec *spec, Alignment *align,
                          int midd, int mida, int minp, int maxp)
{ char *aseq  = align->aseq;
  char *bseq  = align->bseq;
  Path *apath = align->path;

  int     hgh, low, dif;
  int     vlen, vmin, vmax;
  int    *V, *M;
  int    *_V, *_M;
  BVEC   *T;
  BVEC   *_T;

  int    *HA, *NA;
  int    *_HA, *_NA;
  Pebble *cells;
  int     avail, cmax;

  int     TRACE_SPACE = spec->trace_space;
  int     PATH_AVE    = spec->ave_path;
  int16  *SCORE       = spec->score;
  int16  *TABLE       = spec->table;

  int     besta, besty;
  int     trima, trimy, trimd;
  int     trimha;
  int     morea, morey, mored;
  int     moreha;
  int     more, morem, lasta;
  int     aclip, bclip;

  hgh = midd;
  low = midd;
  dif = 0;

  { int span, wing;

    span = (hgh-low)+1;
    vlen = work->vecmax/VectorEn;
    wing = (vlen - span)/2;
    vmin = low - wing;
    vmax = hgh + wing;

    _V  = ((int *) work->vector);
    _M  = _V + vlen;
    _HA = _M + vlen;
    _NA = _HA + vlen;
    _T  = ((BVEC *) (_NA + vlen));

    V  = _V-vmin;
    M  = _M-vmin;
    HA = _HA-vmin;
    NA = _NA-vmin;
    T  = _T-vmin;

    cells = (Pebble *) (work->cells);
    cmax  = work->celmax;
    avail = 0;
  }

  /* Compute 0-wave starting from mid-line */

  more  = 1;
  aclip =  INT32_MAX;
  bclip = -INT32_MAX;

  besta  = trima  = morea = lasta = mida;
  besty  = trimy  = morey = (mida-hgh) >> 1;
  trimd  = mored  = 0;
  trimha = moreha = 0;
  morem  = -1;

  { int   k;
    char *a;

    a  = aseq + hgh;
    for (k = hgh; k >= low; k--)
      { int     y, c, d;
        int     ha, na;
        Pebble *pb;

        y = (mida-k) >> 1;

        if (avail >= cmax-1)
          { cmax  = ((int) (avail*1.2)) + 10000;
            cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
            if (cells == NULL)
              EXIT(1);
            work->celmax = cmax;
            work->cells  = (void *) cells;
          }

        na = ((y+k)/TRACE_SPACE)*TRACE_SPACE;
#ifdef SHOW_TPS
        printf(" A %d: %d,%d,0,%d\n",avail,-1,k,na); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
        pb->mark = na;
        ha  = avail++;
        na += TRACE_SPACE;

        while (1)
          { c = bseq[y];
            if (c == 4)
              { more = 0;
                if (bclip < k)
                  bclip = k;
                break;
              }
            d = a[y];
            if (c != d)
              { if (d == 4) 
                  { more  = 0;
                    aclip = k;
                  }
                break;
              }
            y += 1;
          }
        c = (y << 1) + k;

        while (y+k >= na)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" A %d: %d,%d,0,%d\n",avail,ha,k,na); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = ha;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = na;
            ha  = avail++;
            na += TRACE_SPACE;
          }

        if (c > besta)
          { besta  = trima = lasta = c;
            besty  = trimy = y;
            trimha = ha;
          }

        V[k]  = c;
        T[k]  = PATH_INT;
        M[k]  = PATH_LEN;
        HA[k] = ha;
        NA[k] = na;

        a -= 1;
      }
  }

  if (more == 0)
    { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
        more = 1;
      if (hgh >= aclip)
        { hgh = aclip-1;
          if (morem <= M[aclip])
            { morem  = M[aclip];
              morea  = V[aclip];
              morey  = (morea - aclip)/2;
              moreha = HA[aclip];
            }
        }
      if (low <= bclip)
        { low = bclip+1;
          if (morem <= M[bclip])
            { morem  = M[bclip];
              morea  = V[bclip];
              morey  = (morea - bclip)/2;
              moreha = HA[bclip];
            }
        }
      aclip =  INT32_MAX;
      bclip = -INT32_MAX;
    }

#ifdef DEBUG_WAVE
  printf("\nFORWARD WAVE:\n");
  print_wave(V,M,low,hgh,besta);
#endif

  /* Compute successive waves until no furthest reaching points remain */

  while (more && lasta >= besta - TRIM_MLAG)
    { int     k, n;
      int     ua;
      BVEC    t;
      int     am, ac, ap;
      char   *a;

      if (low <= vmin || hgh >= vmax)
        { int   span, wing;
          int64 move;
          int64 vd, md, had, nad, td;

          span = (hgh-low)+1;
          if (.8*vlen < span)
            { if (enlarge_vector(work,vlen*VectorEn))
                EXIT(1);

              move = ((void *) _V) - work->vector;
              vlen = work->vecmax/VectorEn;

              _V  = (int *) work->vector;
              _M  = _V + vlen;
              _HA = _M + vlen;
              _NA = _HA + vlen;
              _T  = ((BVEC *) (_NA + vlen));
            }
          else
            move = 0;

          wing = (vlen - span)/2;

          vd  = ((void *) ( _V+wing)) - (((void *) ( V+low)) - move);
          md  = ((void *) ( _M+wing)) - (((void *) ( M+low)) - move);
          had = ((void *) (_HA+wing)) - (((void *) (HA+low)) - move);
          nad = ((void *) (_NA+wing)) - (((void *) (NA+low)) - move);
          td  = ((void *) ( _T+wing)) - (((void *) ( T+low)) - move);

          if (vd < 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));
          if (md < 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (had < 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (nad < 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (td < 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));

          if (td > 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));
          if (nad > 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (had > 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (md > 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (vd > 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));

          vmin = low-wing;
          vmax = hgh+wing;

          V  =  _V-vmin;
          M  =  _M-vmin;
          HA = _HA-vmin;
          NA = _NA-vmin;
          T  =  _T-vmin;
        }

      if (low > minp)
        { low -= 1;
          NA[low] = NA[low+1];
          V[low]  = -1;
        }
      if (hgh < maxp)
        { hgh += 1;
          NA[hgh] = NA[hgh-1];
          V[hgh]  = am = -1;
        }
      else
        am = V[hgh];
      dif += 1;

      ac = V[hgh+1] = V[low-1] = -1;
      a  = aseq + hgh;
      t  = PATH_INT;
      n  = PATH_LEN;
      ua = -1;
      for (k = hgh; k >= low; k--)
        { int     y, m;
          int     ha;
          int     c, d;
          BVEC    b;
          Pebble *pb;

          ap = ac;
          ac = am;
          am = V[d = k-1];

          if (ac < am)
            if (am < ap)
              { c  = ap+1;
                m  = n;
                b  = t;
                ha = ua;
              }
            else
              { c  = am+1;
                m  = M[d];
                b  = T[d]; 
                ha = HA[d];
              }
          else
            if (ac < ap)
              { c  = ap+1;
                m  = n;
                b  = t;
                ha = ua;
              }
            else
              { c  = ac+2;
                m  = M[k];
                b  = T[k];
                ha = HA[k];
              }

          if ((b & PATH_TOP) != 0)
            m -= 1;
          b <<= 1;

          y = (c-k) >> 1;
          while (1)
            { c = bseq[y];
              if (c == 4)
                { more = 0;
                  if (bclip < k)
                    bclip = k;
                  break;
                }
              d = a[y];
              if (c != d)
                { if (d == 4) 
                    { more  = 0;
                      aclip = k;
                    }
                  break;
                }
              y += 1;
              if ((b & PATH_TOP) == 0)
                m += 1;
              b = (b << 1) | 1;
            }
          c = (y << 1) + k;

          while (y+k >= NA[k])
            { if (cells[ha].mark < NA[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" A %d: %d,%d,%d,%d\n",avail,ha,k,dif,NA[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = ha;
                  pb->diag = k;
                  pb->diff = dif;
                  pb->mark = NA[k];
                  ha = avail++;
                }
              NA[k] += TRACE_SPACE;
            }

          if (c > besta)
            { besta = c;
              besty = y;
              if (m >= PATH_AVE)
                { lasta = c;
                  if (TABLE[b & TRIM_MASK] >= 0)
                    if (TABLE[(b >> TRIM_LEN) & TRIM_MASK] + SCORE[b & TRIM_MASK] >= 0)
                      { trima  = c;
                        trimy  = y;
                        trimd  = dif;
                        trimha = ha;
                      }
                }
            }

          t  = T[k];
          n  = M[k];
          ua = HA[k];
          V[k]  = c;
          T[k]  = b;
          M[k]  = m;
          HA[k] = ha;

          a -= 1;
        }

      if (more == 0)
        { if (bseq[besty] != 4 && aseq[besta-besty] != 4)
            more = 1;
          if (hgh >= aclip)
            { hgh = aclip-1;
              if (morem <= M[aclip])
                { morem  = M[aclip];
                  morea  = V[aclip];
                  morey  = (morea - aclip)/2;
                  mored  = dif;
                  moreha = HA[aclip];
                }
            }
          if (low <= bclip)
            { low = bclip+1;
              if (morem <= M[bclip])
                { morem  = M[bclip];
                  morea  = V[bclip];
                  morey  = (morea - bclip)/2;
                  mored  = dif;
                  moreha = HA[bclip];
                }
            }
          aclip =  INT32_MAX;
          bclip = -INT32_MAX;
        }

      n = besta - WAVE_LAG;
      while (hgh >= low)
        if (V[hgh] < n)
          hgh -= 1;                               
        else
          { while (V[low] < n)
              low += 1;
            break;
          }

#ifdef WAVE_STATS
      k = (hgh-low)+1;
      if (k > MAX)
        MAX = k;
      TOT += k;
      NWV += 1;
#endif

#ifdef DEBUG_WAVE
      print_wave(V,M,low,hgh,besta);
#endif
    }

  { uint16 *atrace = (uint16 *) apath->trace;
    int     atlen;
    int     trimx;
    int     a, b, k, h;
    int     d, e;

    if (morem >= 0)
      { trimx  = morea-morey;
        trimy  = morey;
        trimd  = mored;
        trimha = moreha;
      }
    else
      trimx = trima-trimy;

    atlen = 0;

    a = -1;
    for (h = trimha; h >= 0; h = b)
      { b = cells[h].ptr; 
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = (mida-k)/2;
    e = 0;
#ifdef SHOW_TRAIL
    printf("  A path = (%5d,%5d)\n",(mida+k)/2,b); fflush(stdout);
#endif
    for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
      { k = cells[h].diag;
        a = cells[h].mark - k;
        d = cells[h].diff;
        atrace[atlen++] = (uint16) (d-e);
        atrace[atlen++] = (uint16) (a-b);
#ifdef SHOW_TRAIL
        printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,a-b); fflush(stdout);
#endif
        b = a;
        e = d;
      }
    if (b+k != trimx)
      { atrace[atlen++] = (uint16) (trimd-e);
        atrace[atlen++] = (uint16) (trimy-b);
#ifdef SHOW_TRAIL
        printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimy-b); fflush(stdout);
#endif
      }
    else if (b != trimy)
      { atrace[atlen-1] = (uint16) (atrace[atlen-1] + (trimy-b));
        atrace[atlen-2] = (uint16) (atrace[atlen-2] + (trimd-e));
#ifdef SHOW_TRAIL
        printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,trimy-b); fflush(stdout);
#endif
      }

    apath->aepos = trimx;
    apath->bepos = trimy;
    apath->diffs = trimd;
    apath->tlen  = atlen;
  }

  return (0);
}

static int reverse_extend(_Work_Data *work, _Align_Spec *spec, Alignment *align,
                          int midd, int mida, int minp, int maxp)
{ char *aseq  = align->aseq - 1;
  char *bseq  = align->bseq - 1;
  Path *apath = align->path;

  int     hgh, low, dif;
  int     vlen, vmin, vmax;
  int    *V, *M;
  int    *_V, *_M;
  BVEC   *T;
  BVEC   *_T;

  int    *HA, *NA;
  int    *_HA, *_NA;
  Pebble *cells;
  int     avail, cmax;

  int     TRACE_SPACE = spec->trace_space;
  int     PATH_AVE    = spec->ave_path;
  int16  *SCORE       = spec->score;
  int16  *TABLE       = spec->table;

  int     besta, besty;
  int     trima, trimy, trimd;
  int     trimha;
  int     morea, morey, mored;
  int     moreha;
  int     more, morem, lasta;
  int     aclip, bclip;

  hgh = midd;
  low = midd;
  dif = 0;

  { int span, wing;

    span = (hgh-low)+1;
    vlen = work->vecmax/VectorEn;
    wing = (vlen - span)/2;
    vmin = low - wing;
    vmax = hgh + wing;

    _V  = ((int *) work->vector);
    _M  = _V + vlen;
    _HA = _M + vlen;
    _NA = _HA + vlen;
    _T  = ((BVEC *) (_NA + vlen));

    V  = _V-vmin;
    M  = _M-vmin;
    HA = _HA-vmin;
    NA = _NA-vmin;
    T  = _T-vmin;

    cells = (Pebble *) (work->cells);
    cmax  = work->celmax;
    avail = 0;
  }

  more  = 1;
  aclip = -INT32_MAX;
  bclip =  INT32_MAX;

  besta  = trima  = morea = lasta = mida;
  besty  = trimy  = morey = (mida-hgh) >> 1;
  trimd  = mored  = 0;
  trimha = moreha = 0;
  morem  = -1;

  { int   k;
    char *a;

    a = aseq + low;
    for (k = low; k <= hgh; k++)
      { int     y, c, d;
        int     ha, na;
        Pebble *pb;

        y = (mida-k) >> 1;

        if (avail >= cmax-1)
          { cmax  = ((int) (avail*1.2)) + 10000;
            cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
            if (cells == NULL)
              EXIT(1);
            work->celmax = cmax;
            work->cells  = (void *) cells;
          }

        na = ((y+k+TRACE_SPACE-1)/TRACE_SPACE-1)*TRACE_SPACE;
#ifdef SHOW_TPS
        printf(" A %d: -1,%d,0,%d\n",avail,k,na+TRACE_SPACE); fflush(stdout);
#endif
        pb = cells+avail;
        pb->ptr  = -1;
        pb->diag = k;
        pb->diff = 0;
        pb->mark = y+k;
        ha  = avail++;

        while (1)
          { c = bseq[y];
            if (c == 4)
              { more = 0;
                if (bclip > k)
                  bclip = k;
                break;
              }
            d = a[y];
            if (c != d)
              { if (d == 4) 
                  { more  = 0;
                    aclip = k;
                  }
                break;
              }
            y -= 1;
          }
        c = (y << 1) + k;

        while (y+k <= na)
          { if (avail >= cmax)
              { cmax  = ((int) (avail*1.2)) + 10000;
                cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),"Reallocating trace cells");
                if (cells == NULL)
                  EXIT(1);
                work->celmax = cmax;
                work->cells  = (void *) cells;
              }
#ifdef SHOW_TPS
            printf(" A %d: %d,%d,0,%d\n",avail,ha,k,na); fflush(stdout);
#endif
            pb = cells+avail;
            pb->ptr  = ha;
            pb->diag = k;
            pb->diff = 0;
            pb->mark = na;
            ha  = avail++;
            na -= TRACE_SPACE;
          }

        if (c < besta)
          { besta  = trima = lasta = c;
            besty  = trimy = y;
            trimha = ha;
          }

        V[k]  = c;
        T[k]  = PATH_INT;
        M[k]  = PATH_LEN;
        HA[k] = ha;
        NA[k] = na;

        a += 1;
      }
  }

  if (more == 0)
    { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
        more = 1;
      if (low <= aclip)
        { low = aclip+1;
          if (morem <= M[aclip])
            { morem  = M[aclip];
              morea  = V[aclip];
              morey  = (morea - aclip)/2;
              moreha = HA[aclip];
            }
        }
      if (hgh >= bclip)
        { hgh = bclip-1;
          if (morem <= M[bclip])
            { morem  = M[bclip];
              morea  = V[bclip];
              morey  = (morea - bclip)/2;
              moreha = HA[bclip];
            }
        }
      aclip = -INT32_MAX;
      bclip =  INT32_MAX;
    }

#ifdef DEBUG_WAVE
  printf("\nREVERSE WAVE:\n");
  print_wave(V,M,low,hgh,besta);
#endif

  while (more && lasta <= besta + TRIM_MLAG)
    { int    k, n;
      int    ua;
      BVEC   t;
      int    am, ac, ap;
      char  *a;

      if (low <= vmin || hgh >= vmax)
        { int   span, wing;
          int64 move, vd, md, had, nad, td;

          span = (hgh-low)+1;
          if (.8*vlen < span)
            { if (enlarge_vector(work,vlen*VectorEn))
                EXIT(1);

              move = ((void *) _V) - work->vector;
              vlen = work->vecmax/VectorEn;

              _V  = (int *) work->vector;
              _M  = _V + vlen;
              _HA = _M + vlen;
              _NA = _HA + vlen;
              _T  = ((BVEC *) (_NA + vlen));
            }
          else
            move = 0;

          wing = (vlen - span)/2;

          vd  = ((void *) ( _V+wing)) - (((void *) ( V+low)) - move);
          md  = ((void *) ( _M+wing)) - (((void *) ( M+low)) - move);
          had = ((void *) (_HA+wing)) - (((void *) (HA+low)) - move);
          nad = ((void *) (_NA+wing)) - (((void *) (NA+low)) - move);
          td  = ((void *) ( _T+wing)) - (((void *) ( T+low)) - move);

          if (vd < 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));
          if (md < 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (had < 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (nad < 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (td < 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));

          if (td > 0)
            memmove( _T+wing,  ((void *) ( T+low)) - move, span*sizeof(BVEC));
          if (nad > 0)
            memmove(_NA+wing,  ((void *) (NA+low)) - move, span*sizeof(int));
          if (had > 0)
            memmove(_HA+wing,  ((void *) (HA+low)) - move, span*sizeof(int));
          if (md > 0)
            memmove( _M+wing,  ((void *) ( M+low)) - move, span*sizeof(int));
          if (vd > 0)
            memmove( _V+wing,  ((void *) ( V+low)) - move, span*sizeof(int));

          vmin = low-wing;
          vmax = hgh+wing;

          V  =  _V-vmin;
          M  =  _M-vmin;
          HA = _HA-vmin;
          NA = _NA-vmin;
          T  =  _T-vmin;
        }

      if (low > minp)
        { low -= 1;
          NA[low] = NA[low+1];
          V[low]  = ap = INT32_MAX;
        }
      else
        ap = V[low]; 
      if (hgh < maxp)
        { hgh += 1;
          NA[hgh] = NA[hgh-1];
          V[hgh] = INT32_MAX;
        }
      dif += 1;

      ac = V[hgh+1] = V[low-1] = INT32_MAX;
      a  = aseq + low;
      t  = PATH_INT;
      n  = PATH_LEN;
      ua = -1;
      for (k = low; k <= hgh; k++)
        { int     y, m;
          int     ha;
          int     c, d;
          BVEC    b;
          Pebble *pb;

          am = ac;
          ac = ap;
          ap = V[d = k+1];

          if (ac > ap)
            if (ap > am)
              { c = am-1;
                m  = n;
                b  = t;
                ha = ua;
              }
            else
              { c  = ap-1;
                m  = M[d];
                b  = T[d];
                ha = HA[d];
              }
          else
            if (ac > am)
              { c  = am-1;
                m  = n;
                b  = t;
                ha = ua;
              }
            else
              { c  = ac-2;
                m  = M[k];
                b  = T[k];
                ha = HA[k];
              }

          if ((b & PATH_TOP) != 0)
            m -= 1;
          b <<= 1;

          y = (c-k) >> 1;
          while (1)
            { c = bseq[y];
              if (c == 4)
                { more = 0;
                  if (bclip > k)
                    bclip = k;
                  break;
                }
              d = a[y];
              if (c != d)
                { if (d == 4) 
                    { more  = 0;
                      aclip = k;
                    }
                  break;
                }
              y -= 1;
              if ((b & PATH_TOP) == 0)
                m += 1;
              b = (b << 1) | 1;
            }
          c = (y << 1) + k;

          while (y+k <= NA[k])
            { if (cells[ha].mark > NA[k])
                { if (avail >= cmax)
                    { cmax  = ((int) (avail*1.2)) + 10000;
                      cells = (Pebble *) Realloc(cells,cmax*sizeof(Pebble),
                                                       "Reallocating trace cells");
                      if (cells == NULL)
                        EXIT(1);
                      work->celmax = cmax;
                      work->cells  = (void *) cells;
                    }
#ifdef SHOW_TPS
                  printf(" A %d: %d,%d,%d,%d\n",avail,ha,k,dif,NA[k]); fflush(stdout);
#endif
                  pb = cells+avail;
                  pb->ptr  = ha;
                  pb->diag = k;
                  pb->diff = dif;
                  pb->mark = NA[k];
                  ha = avail++;
                }
              NA[k] -= TRACE_SPACE;
            }

          if (c < besta)
            { besta = c;
              besty = y;
              if (m >= PATH_AVE)
                { lasta = c;
                  if (TABLE[b & TRIM_MASK] >= 0)
                    if (TABLE[(b >> TRIM_LEN) & TRIM_MASK] + SCORE[b & TRIM_MASK] >= 0)
                      { trima  = c;
                        trimy  = y;
                        trimd  = dif;
                        trimha = ha;
                      }
                }
            }

          t  = T[k];
          n  = M[k];
          ua = HA[k];
          V[k]  = c;
          T[k]  = b;
          M[k]  = m;
          HA[k] = ha;

          a += 1;
        }

      if (more == 0)
        { if (bseq[besty] != 4 && aseq[besta - besty] != 4)
            more = 1;
          if (low <= aclip)
            { low = aclip+1;
              if (morem <= M[aclip])
                { morem  = M[aclip];
                  morea  = V[aclip];
                  morey  = (morea - aclip)/2;
                  mored  = dif;
                  moreha = HA[aclip];
                }
            }
          if (hgh >= bclip)
            { hgh = bclip-1;
              if (morem <= M[bclip])
                { morem  = M[bclip];
                  morea  = V[bclip];
                  morey  = (morea - bclip)/2;
                  mored  = dif;
                  moreha = HA[bclip];
                }
            }
          aclip = -INT32_MAX;
          bclip =  INT32_MAX;
        }

      n = besta + WAVE_LAG;
      while (hgh >= low)
        if (V[hgh] > n)
          hgh -= 1;                               
        else
          { while (V[low] > n)
              low += 1;
            break;
          }

#ifdef WAVE_STATS
      k = (hgh-low)+1;
      if (k > MAX)
        MAX = k;
      TOT += k;
      NWV += 1;
#endif

#ifdef DEBUG_WAVE
      print_wave(V,M,low,hgh,besta);
#endif
    }

  { uint16 *atrace = (uint16 *) apath->trace;
    int     atlen;
    int     trimx;
    int     a, b, k, h;
    int     d, e;

    if (morem >= 0)
      { trimx  = morea-morey;
        trimy  = morey;
        trimd  = mored;
        trimha = moreha;
      }
    else
      trimx = trima-trimy;

    atlen = 0;

    a = -1;
    for (h = trimha; h >= 0; h = b)
      { b = cells[h].ptr; 
        cells[h].ptr = a;
        a = h;
      }
    h = a;

    k = cells[h].diag;
    b = cells[h].mark - k;
    e = 0;
#ifdef SHOW_TRAIL
    printf("  A path = (%5d,%5d)\n",b+k,b); fflush(stdout);
#endif
    if ((b+k)%TRACE_SPACE != 0)
      { h = cells[h].ptr;
        if (h < 0)
          { a = trimy;
            d = trimd;
          }
        else
          { k = cells[h].diag;
            a = cells[h].mark - k;
            d = cells[h].diff;
          }
#ifdef SHOW_TRAIL
        printf("    +%4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,b-a); fflush(stdout);
#endif
        atrace[--atlen] = (uint16) (b-a);
        atrace[--atlen] = (uint16) (d-e);
        b = a;
        e = d;
      }
    if (h >= 0)
      { for (h = cells[h].ptr; h >= 0; h = cells[h].ptr)
          { k = cells[h].diag;
            a = cells[h].mark - k;
            atrace[--atlen] = (uint16) (b-a);
            d = cells[h].diff;
            atrace[--atlen] = (uint16) (d-e);
#ifdef SHOW_TRAIL
            printf("     %4d: (%5d,%5d): %3d / %3d\n",h,a+k,a,d-e,b-a); fflush(stdout);
#endif
            b = a;
            e = d;
          }
        if (b+k != trimx)
          { atrace[--atlen] = (uint16) (b-trimy);
            atrace[--atlen] = (uint16) (trimd-e);
#ifdef SHOW_TRAIL
            printf("           (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimy); fflush(stdout);
#endif
          }
        else if (b != trimy)
          { atrace[atlen+1] = (uint16) (atrace[atlen+1] + (b-trimy));
            atrace[atlen]   = (uint16) (atrace[atlen]   + (trimd-e));
#ifdef SHOW_TRAIL
            printf("         @ (%5d,%5d): %3d / %3d\n",trimx,trimy,trimd-e,b-trimy); fflush(stdout);
#endif
          }
      }

    apath->abpos = trimx;
    apath->bbpos = trimy;
    apath->diffs = trimd;
    apath->tlen  = - atlen;
    apath->trace = atrace + atlen;
  }

  return (0);
}


/* Find the longest local alignment between aseq and bseq through (xcnt,ycnt)
   See associated .h file for the precise definition of the interface.
*/

int Find_Extension(Alignment *align, Work_Data *ework, Align_Spec *espec,
                   int diag, int anti, int lbord, int hbord, int prefix)
{ _Work_Data  *work = ( _Work_Data *) ework;
  _Align_Spec *spec = (_Align_Spec *) espec;

  Path *apath;
  int   minp, maxp;

  { int alen, blen;
    int maxtp, wsize;

    alen = align->alen;
    blen = align->blen;

    wsize = VectorEn*10000;
    if (wsize >= work->vecmax)
      if (enlarge_vector(work,wsize))
        EXIT(1);

    if (alen < blen)
      maxtp = 2*(blen/spec->trace_space+2);
    else
      maxtp = 2*(alen/spec->trace_space+2);
    wsize = 2*maxtp*sizeof(uint16);
    if (wsize > work->pntmax)
      if (enlarge_points(work,wsize))
        EXIT(1);

    apath = align->path;
    apath->trace = ((uint16 *) work->points) + maxtp;
  }

#ifdef DEBUG_PASSES
  printf("\n");
#endif

  if (lbord < 0)
    minp = -INT32_MAX;
  else
    minp = diag-lbord;
  if (hbord < 0)
    maxp = INT32_MAX;
  else
    maxp = diag+hbord;

  if (prefix)
    { if (reverse_extend(work,spec,align,diag,anti,minp,maxp))
        EXIT(1);
      apath->aepos = (anti+diag)/2;
      apath->bepos = (anti-diag)/2;
#ifdef DEBUG_PASSES
      printf("E1 (%d,%d) => (%d,%d) %d\n",
             (anti+diag)/2,(anti-diag)/2,apath->abpos,apath->bbpos,apath->diffs);
#endif
    }
  else
    { if (forward_extend(work,spec,align,diag,anti,minp,maxp))
        EXIT(1);
      apath->abpos = (anti+diag)/2;
      apath->bbpos = (anti-diag)/2;
#ifdef DEBUG_PASSES
      printf("F1 (%d,%d) => (%d,%d) %d\n",
             (anti+diag)/2,(anti-diag)/2,apath->aepos,apath->bepos,apath->diffs);
#endif
     }

#ifdef DEBUG_POINTS
  { uint16 *trace = (uint16 *) apath->trace;
    int     a, h;

    printf("\nA-path (%d,%d)->(%d,%d)",apath->abpos,apath->bbpos,apath->aepos,apath->bepos);
    printf(" %c\n",(COMP(align->flags) ? 'c' : 'n'));
    a = apath->bbpos;
    for (h = 1; h < apath->tlen; h += 2)
      { int dif = trace[h-1];
        int del = trace[h];
        a += del;
        printf("      %d / %d (%d)\n",dif,del,a);
      }
  }
#endif

  return (0);
}


/****************************************************************************************\
*                                                                                        *
*  OVERLAP MANIPULATION                                                                  *
*                                                                                        *
\****************************************************************************************/

static int64 PtrSize   = sizeof(void *);
static int64 OvlIOSize = sizeof(Overlap) - sizeof(void *);

int Read_Overlap(FILE *input, Overlap *ovl)
{ if (fread( ((char *) ovl) + PtrSize, OvlIOSize, 1, input) != 1)
    return (1);
  return (0);
}

int Read_Trace(FILE *input, Overlap *ovl, int tbytes)
{ if (tbytes > 0 && ovl->path.tlen > 0)
    { if (fread(ovl->path.trace, tbytes*ovl->path.tlen, 1, input) != 1)
        return (1);
    }
  return (0);
}

int Write_Overlap(FILE *output, Overlap *ovl, int tbytes)
{ if (fwrite( ((char *) ovl) + PtrSize, OvlIOSize, 1, output) != 1)
    return (1);
  if (ovl->path.trace != NULL)
    if (fwrite(ovl->path.trace,tbytes,ovl->path.tlen,output) != (size_t) ovl->path.tlen)
      return (1);
  return (0);
}

int Compress_TraceTo8(Overlap *ovl, int check)
{ uint16 *t16 = (uint16 *) ovl->path.trace;
  uint8  *t8  = (uint8  *) ovl->path.trace;
  int     j, x;

  if (check)
    for (j = 0; j < ovl->path.tlen; j++)
      { x = t16[j];
        if (x > 255)
          { fprintf(stderr,"%s: Compression of trace to bytes fails, value too big\n",Prog_Name);
            EXIT(1);
          }
        t8[j] = (uint8) x;
      }
  else
    for (j = 0; j < ovl->path.tlen; j++)
      t8[j] = (uint8) (t16[j]);
  return (0);
}

void Decompress_TraceTo16(Overlap *ovl)
{ uint16 *t16 = (uint16 *) ovl->path.trace;
  uint8  *t8  = (uint8  *) ovl->path.trace;
  int     j;

  for (j = ovl->path.tlen-1; j >= 0; j--)
    t16[j] = t8[j];
}

void Print_Overlap(FILE *output, Overlap *ovl, int tbytes, int indent)
{ int     i;

  fprintf(output,"%*s%d vs. ",indent,"",ovl->aread);
  if (COMP(ovl->flags))
    fprintf(output,"c(%d)\n",ovl->bread);
  else
    fprintf(output,"%d\n",ovl->bread);
  fprintf(output,"%*s  [%d,%d] vs [%d,%d] w. %d diffs\n",indent,"",
                 ovl->path.abpos,ovl->path.aepos,ovl->path.bbpos,ovl->path.bepos,ovl->path.diffs);

  if (tbytes == 1)
    { uint8 *trace = (uint8 *) (ovl->path.trace);
      if (trace != NULL)
        { int p = ovl->path.bbpos + trace[1];
          fprintf(output,"%*sTrace: %3d/%5d",indent,"",trace[0],p);
          for (i = 3; i < ovl->path.tlen; i += 2)
            { if (i%10 == 0)
                fprintf(output,"\n%*s",indent+6,"");
              p += trace[i];
              fprintf(output," %3d/%5d",trace[i-1],p);
            }
          fprintf(output,"\n");
        }
    }
  else
    { uint16 *trace = (uint16 *) (ovl->path.trace);
      if (trace != NULL)
        { int p = ovl->path.bbpos + trace[1];
          fprintf(output,"%*sTrace: %3d/%5d",indent,"",trace[0],p);
          for (i = 3; i < ovl->path.tlen; i += 2)
            { if (i%10 == 0)
                fprintf(output,"\n%*s",indent+6,"");
              p += trace[i];
              fprintf(output," %3d/%5d",trace[i-1],p);
            }
          fprintf(output,"\n");
        }
    }
}

int Check_Trace_Points(Overlap *ovl, int tspace, int verbose, char *fname)
{ int i, p, q; 

  if (tspace != 0)
    { if (((ovl->path.aepos-1)/tspace - ovl->path.abpos/tspace)*2 != ovl->path.tlen-2)
        { if (verbose) 
            EPRINTF(EPLACE,"  %s: Wrong number of trace points\n",fname);
          return (1);
        }         
      p = ovl->path.bbpos;
      if (tspace <= TRACE_XOVR)
        { uint8 *trace8 = (uint8 *) ovl->path.trace;
          for (i = 1; i < ovl->path.tlen; i += 2)
            p += trace8[i];
        }
      else      
        { uint16 *trace16 = (uint16 *) ovl->path.trace;
          for (i = 1; i < ovl->path.tlen; i += 2)
            p += trace16[i];
        }
      if (p != ovl->path.bepos)
        { if (verbose)
            EPRINTF(EPLACE,"  %s: Trace point sum != aligned interval\n",fname);
          return (1); 
        }         
    }
  else
    { uint16 *trace16 = (uint16 *) ovl->path.trace;

      p = ovl->path.bbpos;
      q = ovl->path.abpos;
      for (i = 1; i < ovl->path.tlen; i += 2)
        { p += trace16[i];
          q += trace16[i-1];
        }
      if (p != ovl->path.bepos || q != ovl->path.aepos)
        { if (verbose)
            EPRINTF(EPLACE,"  %s: Trace point sum != aligned interval\n",fname);
          return (1); 
        }         
    }
  return (0);
}


void Flip_Alignment(Alignment *align, int full)
{ char *aseq  = align->aseq;
  char *bseq  = align->bseq;
  int   alen  = align->alen;
  int   blen  = align->blen;
  Path *path  = align->path;
  int   comp  = COMP(align->flags);

  int  *trace = (int *) path->trace;
  int   tlen  = path->tlen;

  int   i, j, p;

  if (comp)
    { p = path->abpos;
      path->abpos = blen - path->bepos;
      path->bepos = alen - p;
      p = path->aepos;
      path->aepos = blen - path->bbpos;
      path->bbpos = alen - p;

      if (full)
        { alen += 2;
          blen += 2;

          for (i = 0; i < tlen; i++)
            if ((p = trace[i]) < 0)
              trace[i] = alen + p;
            else
              trace[i] = p - blen;

          i = tlen-1;
          j = 0;
          while (j < i)
            { p = trace[i];
              trace[i] = trace[j];
              trace[j] = p;
              i -= 1;
              j += 1;
            }

          alen -= 2;
          blen -= 2;
        }
    }
  else
    { p = path->abpos;
      path->abpos = path->bbpos;
      path->bbpos = p;
      p = path->aepos;
      path->aepos = path->bepos;
      path->bepos = p;

      if (full)
        for (i = 0; i < tlen; i++)
          trace[i] = - (trace[i]);
    }

  align->aseq  = bseq;
  align->bseq  = aseq;
  align->alen  = blen;
  align->blen  = alen;
}


/****************************************************************************************\
*                                                                                        *
*  ALIGNMENT PRINTING                                                                    *
*                                                                                        *
\****************************************************************************************/

/* Complement the sequence in fragment aseq.  The operation does the
   complementation/reversal in place.  Calling it a second time on a
   given fragment restores it to its original state.                */

void Complement_Seq(char *aseq, int len)
{ char *s, *t;
  int   c;

  s = aseq;
  t = aseq + (len-1);
  while (s < t)
    { c    = 3 - *s;
      *s++ = (char) (3 - *t);
      *t-- = (char) c;
    }
  if (s == t)
    *s = (char) (3 - *s);
}


/* Print an alignment to file between a and b given in trace (unpacked).
   Prefix gives the length of the initial prefix of a that is unaligned.  */

static char ToL[8] = { 'a', 'c', 'g', 't', '.', '[', ']', '-' };
static char ToU[8] = { 'A', 'C', 'G', 'T', '.', '[', ']', '-' };

int Print_Alignment(FILE *file, Alignment *align, Work_Data *ework,
                    int indent, int width, int border, int upper, int coord)
{ _Work_Data *work  = (_Work_Data *) ework;
  int        *trace = align->path->trace;
  int         tlen  = align->path->tlen;

  char *Abuf, *Bbuf, *Dbuf;
  int   i, j, o;
  char *a, *b;
  char  mtag, dtag;
  int   prefa, prefb;
  int   aend, bend;
  int   comp, blen;
  int   sa, sb;
  int   match, diff;
  char *N2A;

  if (trace == NULL) return (0);

#ifdef SHOW_TRACE
  fprintf(file,"\nTrace:\n");
  for (i = 0; i < tlen; i++)
    fprintf(file,"  %3d\n",trace[i]);
#endif

  o = sizeof(char)*3*(width+1);
  if (o > work->vecmax)
    if (enlarge_vector(work,o))
      EXIT(1);

  if (upper)
    N2A = ToU;
  else
    N2A = ToL;

  Abuf = (char *) work->vector;
  Bbuf = Abuf + (width+1);
  Dbuf = Bbuf + (width+1);

  aend = align->path->aepos;
  bend = align->path->bepos;

  comp = COMP(align->flags);
  blen = align->blen;

  Abuf[width] = Bbuf[width] = Dbuf[width] = '\0';

  // buffer/output next column

#define COLUMN(x,y)							\
{ int u, v;								\
  if (o >= width)							\
    { fprintf(file,"\n");						\
      fprintf(file,"%*s",indent,"");					\
      if (coord > 0)							\
        { if (sa < aend)						\
            fprintf(file," %*d",coord,sa);				\
          else								\
            fprintf(file," %*s",coord,"");				\
          fprintf(file," %s\n",Abuf);					\
          fprintf(file,"%*s %*s %s\n",indent,"",coord,"",Dbuf);		\
          fprintf(file,"%*s",indent,"");				\
          if (sb < bend)						\
            if (comp)							\
              fprintf(file," %*d",coord,blen-sb);			\
            else							\
              fprintf(file," %*d",coord,sb);				\
          else								\
            fprintf(file," %*s",coord,"");				\
          fprintf(file," %s",Bbuf);					\
        }								\
      else								\
        { fprintf(file," %s\n",Abuf);					\
          fprintf(file,"%*s %s\n",indent,"",Dbuf);			\
          fprintf(file,"%*s %s",indent,"",Bbuf);			\
        }								\
      fprintf(file," %5.1f%%\n",(100.*diff)/(diff+match));		\
      o  = 0;								\
      sa = i-1;								\
      sb = j-1;								\
      match = diff = 0;							\
    }									\
  u = (x);								\
  v = (y);								\
  if (u == 4 || v == 4)							\
    Dbuf[o] = ' ';							\
  else if (u == v)							\
    Dbuf[o] = mtag;							\
  else									\
    Dbuf[o] = dtag;							\
  Abuf[o] = N2A[u];							\
  Bbuf[o] = N2A[v];							\
  o += 1;								\
}

  a = align->aseq - 1;
  b = align->bseq - 1;

  o  = 0;
  i = align->path->abpos;
  j = align->path->bbpos;

  prefa = 0;
  for (prefa = 0; prefa < border && a[i] != 4; prefa++)
    i -= 1;
  i += 1;

  prefb = 0;
  for (prefb = 0; prefb < border && b[j] != 4; prefb++)
    j -= 1;
  j += 1;

  sa   = i-1;
  sb   = j-1;
  mtag = ':';
  dtag = ':';

  while (prefa > prefb)
    { COLUMN(a[i],4)
      i += 1;
      prefa -= 1;
    }
  while (prefb > prefa)
    { COLUMN(4,b[j])
      j += 1;
      prefb -= 1;
    }
  while (prefa > 0)
    { COLUMN(a[i],b[j])
      i += 1;
      j += 1;
      prefa -= 1;
    }

  mtag = '[';
  if (prefb > 0)
    COLUMN(5,5)

  mtag  = '|';
  dtag  = '*';

  match = diff = 0;

  { int p, c;      /* Output columns of alignment til reach trace end */

    for (c = 0; c < tlen; c++)
      if ((p = trace[c]) < 0)
        { p = -p;
          while (i != p)
            { COLUMN(a[i],b[j])
              if (a[i] == b[j])
                match += 1;
              else
                diff += 1;
              i += 1;
              j += 1;
            }
          COLUMN(7,b[j])
          j += 1;
          diff += 1;
        }
      else
        { while (j != p)
            { COLUMN(a[i],b[j])
              if (a[i] == b[j])
                match += 1;
              else
                diff += 1;
              i += 1;
              j += 1;
            }
          COLUMN(a[i],7)
          i += 1;
          diff += 1;
        }
    p = align->path->aepos;
    while (i <= p)
      { COLUMN(a[i],b[j])
        if (a[i] == b[j])
          match += 1;
        else
          diff += 1;
        i += 1;
        j += 1;
      }
  }

  { int c;     /* Output remaining column including unaligned suffix */

    mtag = ']';
    if (a[i] != 4 && b[j] != 4 && border > 0)
      COLUMN(6,6)

    mtag = ':';
    dtag = ':';

    c = 0;
    while (c < border && (a[i] != 4 || b[j] != 4))
      { if (a[i] != 4)
          if (b[j] != 4)
            { COLUMN(a[i],b[j])
              i += 1;
              j += 1;
            }
          else
            { COLUMN(a[i],4)
              i += 1;
            }
        else
          { COLUMN(4,b[j])
            j += 1;
          }
        c += 1;
      }
  }

  /* Print remainder of buffered col.s */

  fprintf(file,"\n");
  fprintf(file,"%*s",indent,"");
  if (coord > 0)
    { if (sa < aend)
        fprintf(file," %*d",coord,sa);
      else
        fprintf(file," %*s",coord,"");
      fprintf(file," %.*s\n",o,Abuf);
      fprintf(file,"%*s %*s %.*s\n",indent,"",coord,"",o,Dbuf);
      fprintf(file,"%*s",indent,"");
      if (sb < bend)
        if (comp)
          fprintf(file," %*d",coord,blen-sb);
        else
          fprintf(file," %*d",coord,sb);
      else
        fprintf(file," %*s",coord,"");
      fprintf(file," %.*s",o,Bbuf);
    }
  else
    { fprintf(file," %.*s\n",o,Abuf);
      fprintf(file,"%*s %.*s\n",indent,"",o,Dbuf);
      fprintf(file,"%*s %.*s",indent,"",o,Bbuf);
    }
  if (diff+match > 0)
    fprintf(file," %5.1f%%\n",(100.*diff)/(diff+match));
  else
    fprintf(file,"\n");

  fflush(file);
  return (0);
}

int Print_Reference(FILE *file, Alignment *align, Work_Data *ework,
                    int indent, int block, int border, int upper, int coord)
{ _Work_Data *work  = (_Work_Data *) ework;
  int        *trace = align->path->trace;
  int         tlen  = align->path->tlen;

  char *Abuf, *Bbuf, *Dbuf;
  int   i, j, o;
  char *a, *b;
  char  mtag, dtag;
  int   prefa, prefb;
  int   aend, bend;
  int   comp, blen;
  int   sa, sb, s0;
  int   match, diff;
  char *N2A;
  int   vmax;

  if (trace == NULL) return (0);

#ifdef SHOW_TRACE
  fprintf(file,"\nTrace:\n");
  for (i = 0; i < tlen; i++)
    fprintf(file,"  %3d\n",trace[i]);
#endif

  vmax = work->vecmax/3;
  o = sizeof(char)*6*(block+1);
  if (o > vmax)
    { if (enlarge_vector(work,3*o))
        EXIT(1);
      vmax = work->vecmax/3;
    }

  Abuf = (char *) work->vector;
  Bbuf = Abuf + vmax;
  Dbuf = Bbuf + vmax;

  if (upper)
    N2A = ToU;
  else
    N2A = ToL;

  aend = align->path->aepos;
  bend = align->path->bepos;

  comp = COMP(align->flags);
  blen = align->blen;

#define BLOCK(x,y)							\
{ int u, v;								\
  if (i%block == 1 && i != s0 && x < 4 && o > 0)			\
    { fprintf(file,"\n");						\
      fprintf(file,"%*s",indent,"");					\
      if (coord > 0)							\
        { if (sa < aend)						\
            fprintf(file," %*d",coord,sa);				\
          else								\
            fprintf(file," %*s",coord,"");				\
          fprintf(file," %.*s\n",o,Abuf);				\
          fprintf(file,"%*s %*s %.*s\n",indent,"",coord,"",o,Dbuf);	\
          fprintf(file,"%*s",indent,"");				\
          if (sb < bend)						\
            if (comp)							\
              fprintf(file," %*d",coord,blen-sb);			\
            else							\
              fprintf(file," %*d",coord,sb);				\
          else								\
            fprintf(file," %*s",coord,"");				\
          fprintf(file," %.*s",o,Bbuf);					\
        }								\
      else								\
        { fprintf(file," %.*s\n",o,Abuf);				\
          fprintf(file,"%*s %.*s\n",indent,"",o,Dbuf);			\
          fprintf(file,"%*s %.*s",indent,"",o,Bbuf);			\
        }								\
      fprintf(file," %5.1f%%\n",(100.*diff)/(diff+match));		\
      o  = 0;								\
      sa = i-1;								\
      sb = j-1;								\
      match = diff = 0;							\
    }									\
  u = (x);								\
  v = (y);								\
  if (u == 4 || v == 4)							\
    Dbuf[o] = ' ';							\
  else if (u == v)							\
    Dbuf[o] = mtag;							\
  else									\
    Dbuf[o] = dtag;							\
  Abuf[o] = N2A[u];							\
  Bbuf[o] = N2A[v];							\
  o += 1;								\
  if (o >= vmax)							\
    { if (enlarge_vector(work,3*o))					\
        EXIT(1);							\
      vmax = work->vecmax/3;						\
      memmove(work->vector+2*vmax,Dbuf,o);				\
      memmove(work->vector+vmax,Bbuf,o);				\
      memmove(work->vector,Abuf,o);					\
      Abuf = (char *) work->vector;					\
      Bbuf = Abuf + vmax;						\
      Dbuf = Bbuf + vmax;						\
    }									\
}

  a = align->aseq - 1;
  b = align->bseq - 1;

  o  = 0;
  i = align->path->abpos;
  j = align->path->bbpos;

  prefa = 0;
  for (prefa = 0; prefa < border && a[i] != 4; prefa++)
    i -= 1;
  i += 1;

  prefb = 0;
  for (prefb = 0; prefb < border && b[j] != 4; prefb++)
    j -= 1;
  j += 1;

  s0   = i;
  sa   = i-1;
  sb   = j-1;
  mtag = ':';
  dtag = ':';

  while (prefa > prefb)
    { BLOCK(a[i],4)
      i += 1;
      prefa -= 1;
    }
  while (prefb > prefa)
    { BLOCK(4,b[j])
      j += 1;
      prefb -= 1;
    }
  while (prefa > 0)
    { BLOCK(a[i],b[j])
      i += 1;
      j += 1;
      prefa -= 1;
    }

  mtag = '[';
  if (prefb > 0)
    BLOCK(5,5)

  mtag  = '|';
  dtag  = '*';

  match = diff = 0;

  { int p, c;      /* Output columns of alignment til reach trace end */

    for (c = 0; c < tlen; c++)
      if ((p = trace[c]) < 0)
        { p = -p;
          while (i != p)
            { BLOCK(a[i],b[j])
              if (a[i] == b[j])
                match += 1;
              else
                diff += 1;
              i += 1;
              j += 1;
            }
          BLOCK(7,b[j])
          j += 1;
          diff += 1;
        }
      else
        { while (j != p)
            { BLOCK(a[i],b[j])
              if (a[i] == b[j])
                match += 1;
              else
                diff += 1;
              i += 1;
              j += 1;
            }
          BLOCK(a[i],7)
          i += 1;
          diff += 1;
        }
    p = align->path->aepos;
    while (i <= p)
      { BLOCK(a[i],b[j])
        if (a[i] == b[j])
		match += 1;
	else
          diff += 1;
        i += 1;
        j += 1;
      }
  }

  { int c;     /* Output remaining column including unaligned suffix */

    mtag = ']';
    if (a[i] != 4 && b[j] != 4 && border > 0)
      BLOCK(6,6)

    mtag = ':';
    dtag = ':';

    c = 0;
    while (c < border && (a[i] != 4 || b[j] != 4))
      { if (a[i] != 4)
          if (b[j] != 4)
            { BLOCK(a[i],b[j])
              i += 1;
              j += 1;
            }
          else
            { BLOCK(a[i],4)
              i += 1;
            }
        else
          { BLOCK(4,b[j])
            j += 1;
          }
        c += 1;
      }
  }

  /* Print remainder of buffered col.s */

  fprintf(file,"\n");
  fprintf(file,"%*s",indent,"");
  if (coord > 0)
    { if (sa < aend)
        fprintf(file," %*d",coord,sa);
      else
        fprintf(file," %*s",coord,"");
      fprintf(file," %.*s\n",o,Abuf);
      fprintf(file,"%*s %*s %.*s\n",indent,"",coord,"",o,Dbuf);
      fprintf(file,"%*s",indent,"");
      if (sb < bend)
        if (comp)
          fprintf(file," %*d",coord,blen-sb);
        else
          fprintf(file," %*d",coord,sb);
      else
        fprintf(file," %*s",coord,"");
      fprintf(file," %.*s",o,Bbuf);
    }
  else
    { fprintf(file," %.*s\n",o,Abuf);
      fprintf(file,"%*s %.*s\n",indent,"",o,Dbuf);
      fprintf(file,"%*s %.*s",indent,"",o,Bbuf);
    }
  if (diff+match > 0)
    fprintf(file," %5.1f%%\n",(100.*diff)/(diff+match));
  else
    fprintf(file,"\n");

  fflush(file);
  return (0);
}

/* Print an ASCII representation of the overlap in align between fragments
   a and b to given file.                                                  */

static inline void repchar(FILE *file, int symbol, int rep)
{ while (rep-- > 0)
    fputc(symbol,file);
}

void Alignment_Cartoon(FILE *file, Alignment *align, int indent, int coord)
{ int   alen = align->alen;
  int   blen = align->blen;
  Path *path = align->path;
  int   comp = COMP(align->flags);
  int   w;

  fprintf(file,"%*s",indent,"");
  if (path->abpos > 0)
    fprintf(file,"    %*d ",coord,path->abpos);
  else
    fprintf(file,"%*s",coord+5,"");
  if (path->aepos < alen)
    fprintf(file,"%*s%d",coord+8,"",alen-path->aepos);
  fprintf(file,"\n");

  fprintf(file,"%*s",indent,"");
  if (path->abpos > 0)
    { fprintf(file,"A ");
      w = Number_Digits((int64) path->abpos);
      repchar(file,' ',coord-w);
      repchar(file,'=',w+3);
      fputc('+',file);
      repchar(file,'-',coord+5);
    }
  else
    { fprintf(file,"A %*s",coord+4,"");
      repchar(file,'-',coord+5);
    }

  if (path->aepos < alen)
    { fputc('+',file);
      w = Number_Digits((int64) (alen-path->aepos));
      repchar(file,'=',w+2);
      fputc('>',file);
      repchar(file,' ',w);
    }
  else
    { fputc('>',file);
      repchar(file,' ',coord+3);
    }

  { int asub, bsub;

    asub = path->aepos - path->abpos;
    bsub = path->bepos - path->bbpos;
    fprintf(file,"   dif/(len1+len2) = %d/(%d+%d) = %5.2f%%\n",
                 path->diffs,asub,bsub,(200.*path->diffs)/(asub+bsub));
  }

  { int   sym1e, sym2e;
    int   sym1p, sym2p;

    if (comp > 0)
      { sym1p = '<'; sym2p = '-'; sym1e = '<'; sym2e = '='; }
    else
      { sym1p = '-'; sym2p = '>'; sym1e = '='; sym2e = '>'; }

    fprintf(file,"%*s",indent,"");
    if (path->bbpos > 0)
      { fprintf(file,"B ");
        w = Number_Digits((int64) path->bbpos);
        repchar(file,' ',coord-w);
        fputc(sym1e,file);
        repchar(file,'=',w+2);
        fputc('+',file);
        repchar(file,'-',coord+5);
      }
    else
      { fprintf(file,"B ");
        repchar(file,' ',coord+3);
        fputc(sym1p,file);
        repchar(file,'-',coord+5);
      }
    if (path->bepos < blen)
      { fprintf(file,"+");
        w = Number_Digits((int64) (blen-path->bepos));
        repchar(file,'=',w+2);
        fprintf(file,"%c\n",sym2e);
      }
    else
      fprintf(file,"%c\n",sym2p);
  }

  fprintf(file,"%*s",indent,"");
  if (path->bbpos > 0)
    fprintf(file,"    %*d ",coord,path->bbpos);
  else
    fprintf(file,"%*s",coord+5,"");
  if (path->bepos < blen)
    fprintf(file,"%*s%d",coord+8,"",blen-path->bepos);
  fprintf(file,"\n");

  fflush(file);
}


/****************************************************************************************\
*                                                                                        *
*  O(ND) trace algorithm                                                                 *
*                                                                                        *
\****************************************************************************************/


#ifdef DEBUG_AWAVE

static void print_awave(int *V, int low, int hgh)
{ int k;

  printf("  [%6d,%6d]: ",low,hgh);
  for (k = low; k <= hgh; k++)
    printf(" %3d",V[k]);
  printf("\n");
  fflush(stdout);
}

#endif

#ifdef DEBUG_ALIGN

static int depth = 0;

#endif

typedef struct
  { int    *Stop;        //  Ongoing stack of alignment indels
    uint16 *Trace;       //  Base of Trace Vector
    char *Aabs, *Babs;   //  Absolute base of A and B sequences

    int  **PVF, **PHF;   //  List of waves for iterative np algorithms
    int   mida,  midb;   //  mid point division for mid-point algorithms

    int   *VF,   *VB;    //  Forward/Reverse waves for nd algorithms
  } Trace_Waves;

static int split_nd(char *A, int M, char *B, int N, Trace_Waves *wave, int *px, int *py)
{ int x, y;
  int D;

  int  *VF = wave->VF;
  int  *VB = wave->VB;
  int   flow;  //  fhgh == D !
  int   blow, bhgh;
  char *a;

  y = 0;
  if (N < M)
    while (y < N && B[y] == A[y])
      y += 1;
  else
    { while (y < M && B[y] == A[y])
        y += 1;
      if (y >= M && N == M)
        { *px = *py = M;
          return (0);
        }
    }

  flow   = 0;
  VF[0]  = y;
  VF[-1] = -2;

  x = N-M;
  a = A-x;
  y = N-1;
  if (N > M)
    while (y >= x && B[y] == a[y])
      y -= 1;
  else
    while (y >= 0 && B[y] == a[y])
      y -= 1;

  blow = bhgh = -x;
  VB += x;
  VB[blow]   = y;
  VB[blow-1] = N+1;

  for (D = 1; 1; D += 1)
    { int   k, r;
      int   am, ac, ap;

      //  Forward wave

      flow -= 1;
      am = ac = VF[flow-1] = -2;

      a = A + D;
      x = M - D;
      for (k = D; k >= flow; k--)
        { ap = ac;
          ac = am+1;
          am = VF[k-1];

          if (ac < am)
            if (ap < am)
              y  = am;
            else
              y = ap;
          else
            if (ap < ac)
              y = ac;
            else
              y = ap;

          if (blow <= k && k <= bhgh)
            { r = VB[k];
              if (y > r)
                { D = (D<<1)-1;
                  if (ap > r)
                    y = ap;
                  else if (ac > r)
                    y = ac;
                  else
                    y = r+1;
                  x = k+y;
                  *px = x;
                  *py = y;
                  return (D);
                }
            }

          if (N < x)
            while (y < N && B[y] == a[y])
              y += 1;
          else
            while (y < x && B[y] == a[y])
              y += 1;
          
          VF[k] = y;
          a -= 1;
          x += 1;
        }

#ifdef DEBUG_AWAVE
      print_awave(VF,flow,D);
#endif

      //  Reverse Wave

      bhgh += 1;
      blow -= 1;
	am = ac = VB[blow-1] = N+1;

      a = A + bhgh;
      x = -bhgh;
      for (k = bhgh; k >= blow; k--)
        { ap = ac+1;
          ac = am;
          am = VB[k-1];

          if (ac > am)
            if (ap > am)
              y  = am;
            else
              y = ap;
          else
            if (ap > ac)
              y = ac;
            else
              y = ap;

          if (flow <= k && k <= D)
            { r = VF[k];
	        if (y <= r)
                { D = (D << 1);
                  if (ap <= r)
                    y = ap;
                  else if (ac <= r)
                    y = ac;
                  else
                    y = r;
                  x = k+y;
                  *px = x;
                  *py = y;
                  return (D);
                }
            }

          y -= 1;
          if (x > 0)
            while (y >= x && B[y] == a[y])
              y -= 1;
          else
            while (y >= 0 && B[y] == a[y])
              y -= 1;

          VB[k] = y;
          a -= 1;
          x += 1;
        }

#ifdef DEBUG_AWAVE
      print_awave(VB,blow,bhgh);
#endif
    }
}

static int trace_nd(char *A, int M, char *B, int N, Trace_Waves *wave, int tspace)
{ int x, y;
  int D, s;

#ifdef DEBUG_ALIGN
  printf("%*s %ld,%ld: %d vs %d\n",depth,"",A-wave->Aabs,B-wave->Babs,M,N);
  fflush(stdout);
#endif

  if (M <= 0)
    { y = (((A-wave->Aabs)/tspace) << 1);
      wave->Trace[y] += N;
      wave->Trace[y+1] += N;
#ifdef DEBUG_TRACE
      printf("%*s Adding1 (%d,%d) to tp %d(%d,%d)\n",depth,"",N,N,y>>1,
                                                     wave->Trace[y+1],wave->Trace[y]);
      fflush(stdout);
#endif
      return (N);
    }

  if (N <= 0)
    { x = A - wave->Aabs;
      y = x/tspace;
      x = (y+1)*tspace - x;
      y <<= 1;
      for (s = M; s > 0; s -= x, x = tspace)
        { if (x > s)
            x = s;
          wave->Trace[y] += x;
#ifdef DEBUG_TRACE
          printf("%*s Adding2 (0,%d) to tp %d(%d,%d)\n",depth,"",x,y>>1,
                                                        wave->Trace[y+1],wave->Trace[y]);
          fflush(stdout);
#endif
          y += 2;
        }
      return (M);
    }

  D = split_nd(A,M,B,N,wave,&x,&y);

  if (D > 1)
    {
#ifdef DEBUG_ALIGN
      printf("%*s (%d,%d) @ %d\n",depth,"",x,y,D);
      fflush(stdout);
      depth += 2;
#endif

      s = A-wave->Aabs;
      if ((s/tspace+1)*tspace - s >= x)
        { s = ((s/tspace)<<1);
          wave->Trace[s] += (D+1)/2;
          wave->Trace[s+1] += y;
#ifdef DEBUG_TRACE
          printf("%*s Adding3 (%d,%d) to tp %d(%d,%d)\n",depth,"",y,(D+1)/2,s>>1,
                                                         wave->Trace[s+1],wave->Trace[s]);
          fflush(stdout);
#endif
        }
      else
        trace_nd(A,x,B,y,wave,tspace);

      s = (A+x)-wave->Aabs;
      if ((s/tspace+1)*tspace - s >= M-x)
        { s = ((s/tspace)<<1);
          wave->Trace[s] += D/2;
          wave->Trace[s+1] += N-y;
#ifdef DEBUG_TRACE
          printf("%*s Adding4 (%d,%d)) to tp %d(%d,%d)\n",depth,"",N-y,D/2,s>>1,
                                                          wave->Trace[s+1],wave->Trace[s]);
          fflush(stdout);
#endif
        }
      else
        trace_nd(A+x,M-x,B+y,N-y,wave,tspace);

#ifdef DEBUG_ALIGN
      depth -= 2;
#endif
    }

  else
    { int u, v;

      if (D == 0 || M < N)
        s = x;
      else
        s = x-1;
      if (s > 0)
        { u = A - wave->Aabs;
          v = u/tspace;
          u = (v+1)*tspace - u;
          for (v <<= 1; s > 0; s -= u, u = tspace)
            { if (u > s)
                u = s;
              wave->Trace[v+1] += u;
#ifdef DEBUG_TRACE
              printf("%*s Adding5 (%d,0)) to tp %d(%d,%d)\n",depth,"",u,v>>1,
                                                             wave->Trace[v+1],wave->Trace[v]);
              fflush(stdout);
#endif
              v += 2;
            }
        }

      if (D == 0)
        return (D);

      if (M < N)
        y = ((((A+x)-wave->Aabs)/tspace)<<1);
      else
        y = ((((A+(x-1))-wave->Aabs)/tspace)<<1);
      wave->Trace[y] += 1;
      if (M <= N)
        wave->Trace[y+1] += 1;
#ifdef DEBUG_TRACE
      printf("%*s Adding5 (%d,1)) to tp %d(%d,%d)\n",depth,"",N>=M,y>>1,
                                                     wave->Trace[y+1],wave->Trace[y]);
      fflush(stdout);
#endif

      s = M-x;
      if (s > 0)
        { u = (A+x) - wave->Aabs;
          v = u/tspace;
          u = (v+1)*tspace - u;
          for (v <<= 1; s > 0; s -= u, u = tspace)
            { if (u > s)
                u = s;
              wave->Trace[v+1] += u;
#ifdef DEBUG_TRACE
              printf("%*s Adding5 (%d,0)) to tp %d(%d,%d)\n",depth,"",u,v>>1,
                                                             wave->Trace[v+1],wave->Trace[v]);
              fflush(stdout);
#endif
              v += 2;
            }
        }
    }

  return (D);
}

static int dandc_nd(char *A, int M, char *B, int N, Trace_Waves *wave)
{ int x, y;
  int D;

#ifdef DEBUG_ALIGN
  printf("%*s %ld,%ld: %d vs %d\n",depth,"",A-wave->Aabs,B-wave->Babs,M,N);
#endif

  if (M <= 0)
    { x = (wave->Aabs-A)-1;
      for (y = 1; y <= N; y++)
        { *wave->Stop++ = x;
#ifdef DEBUG_SCRIPT
          printf("%*s *I %ld(%ld)\n",depth,"",y+(B-wave->Babs),(A-wave->Aabs)+1);
#endif
        }
      return (N);
    }

  if (N <= 0)
    { y = (B-wave->Babs)+1;
      for (x = 1; x <= M; x++)
        { *wave->Stop++ = y;
#ifdef DEBUG_SCRIPT
          printf("%*s *D %ld(%ld)\n",depth,"",x+(A-wave->Aabs),(B-wave->Babs)+1);
#endif
        }
      return (M);
    }

  D = split_nd(A,M,B,N,wave,&x,&y);

  if (D > 1)
    { 
#ifdef DEBUG_ALIGN
      printf("%*s (%d,%d) @ %d\n",depth,"",x,y,D);
      fflush(stdout);
      depth += 2;
#endif

      dandc_nd(A,x,B,y,wave);
      dandc_nd(A+x,M-x,B+y,N-y,wave);

#ifdef DEBUG_ALIGN
      depth -= 2;
#endif
    }

  else if (D == 1)

    { if (M > N)
        { *wave->Stop++ = (B-wave->Babs)+y+1;
#ifdef DEBUG_SCRIPT
          printf("%*s  D %ld(%ld)\n",depth,"",(A-wave->Aabs)+x,(B-wave->Babs)+y+1);
#endif
        }

      else if (M < N)
        { *wave->Stop++ = (wave->Aabs-A)-x-1;
#ifdef DEBUG_SCRIPT
          printf("%*s  I %ld(%ld)\n",depth,"",(B-wave->Babs)+y,(A-wave->Aabs)+x+1);
#endif
        }

#ifdef DEBUG_SCRIPT
      else
        printf("%*s  %ld S %ld\n",depth,"",(wave->Aabs-A)+x,(B-wave->Babs)+y);
#endif
    }

  return (D);
}

int Compute_Alignment(Alignment *align, Work_Data *ework, int task, int tspace)
{ _Work_Data *work = (_Work_Data *) ework;
  Trace_Waves wave;

  int     L, D;
  int     asub, bsub;
  char   *aseq, *bseq;
  Path   *path;
  int    *trace;
  uint16 *strace;

  path = align->path;
  asub = path->aepos-path->abpos;
  bsub = path->bepos-path->bbpos;
  aseq = align->aseq+path->abpos;
  bseq = align->bseq+path->bbpos;

  L = 0;
  if (task != DIFF_ONLY)
    { if (task == DIFF_TRACE || task == PLUS_TRACE)
        L = 2*(((path->aepos + (tspace-1))/tspace - path->abpos/tspace) + 1)*sizeof(uint16);
      else if (asub < bsub)
        L = bsub*sizeof(int);
      else
        L = asub*sizeof(int);
      if (L > work->alnmax)
        if (enlarge_alnpts(work,L))
          EXIT(1);
    }

  trace  = ((int *) work->alnpts);
  strace = ((uint16 *) work->alnpts);

  if (asub > bsub)
    D = (4*asub+6)*sizeof(int);
  else
    D = (4*bsub+6)*sizeof(int);
  if (D > work->vecmax)
    if (enlarge_vector(work,D))
      EXIT(1);
  
  if (asub > bsub)
    { wave.VF = ((int *) work->vector) + (asub+1);
      wave.VB = wave.VF + (2*asub+3);
    }
  else
    { wave.VF = ((int *) work->vector) + (bsub+1);
      wave.VB = wave.VF + (2*bsub+3);
    }

  wave.Aabs = align->aseq;
  wave.Babs = align->bseq;

  if (task == DIFF_ONLY)
    { wave.mida = -1;
      if (asub <= 0)
        path->diffs = bsub;
      else if (bsub <= 0)
        path->diffs = asub;
      else
        path->diffs = split_nd(aseq,asub,bseq,bsub,&wave,&wave.mida,&wave.midb);
      path->trace = NULL;
      path->tlen  = -1;
      return (0);
    }

  else if (task < DIFF_ONLY && wave.mida >= 0)
    { int x = wave.mida;
      int y = wave.midb;

      if (task == PLUS_ALIGN)
        { wave.Stop = trace;
          dandc_nd(aseq,x,bseq,y,&wave);
          dandc_nd(aseq+x,asub-x,bseq+y,bsub-y,&wave);
          path->tlen = wave.Stop - trace;
        }
      else
        { int i, n;

          wave.Trace = strace - 2*(path->abpos/tspace);
          n  = L/sizeof(uint16);
          for (i = 0; i < n; i++)
            strace[i] = 0;

          trace_nd(aseq,x,bseq,y,&wave,tspace);
          trace_nd(aseq+x,asub-x,bseq+y,bsub-y,&wave,tspace);

          if (strace[n-1] != 0)   //  Last element is to capture all inserts on TP boundary
            { strace[n-3] += strace[n-1];
              strace[n-4] += strace[n-2];
            }
          path->tlen = n-2;

#ifdef DEBUG_SCRIPT
          printf("  Trace:\n");
          for (i = 0; i < path->tlen; i += 2)
            printf("    %3d  %3d\n",strace[i],strace[i+1]);
          fflush(stdout);
#endif
        }
    }

  else 
    { if (task == DIFF_ALIGN)
        { wave.Stop   = trace;
          path->diffs = dandc_nd(aseq,asub,bseq,bsub,&wave);
          path->tlen  = wave.Stop - trace;
        }
      else
        { int i, n;

          wave.Trace  = strace - 2*(path->abpos/tspace);
          n  = L/sizeof(uint16);
          for (i = 0; i < n; i++)
            strace[i] = 0;
          path->diffs = trace_nd(aseq,asub,bseq,bsub,&wave,tspace);

          if (strace[n-1] != 0)   //  Last element is to capture all inserts on TP boundary
            { strace[n-3] += strace[n-1];
              strace[n-4] += strace[n-2];
            }
          path->tlen = n-2;

#ifdef DEBUG_SCRIPT
          printf("  Trace:\n");
          for (i = 0; i < path->tlen; i += 2)
            printf("    %3d  %3d\n",strace[i],strace[i+1]);
          fflush(stdout);
#endif
        }
    }

  path->trace = trace;
  return (0);
}


/****************************************************************************************\
*                                                                                        *
*  O(NP) tracing algorithms                                                              *
*                                                                                        *
\****************************************************************************************/

/* Iterative O(np) algorithm for finding the alignment between two substrings (specified
     by a Path record).  The variation includes handling substitutions and guarantees
     to find left-most alignments so that low complexity runs are always aligned in
     the same way.
*/

#ifdef DEBUG_ALIGN

static int ToA[4] = { 'a', 'c', 'g', 't' };

#endif

static char *TP_Align =
         "Bad alignment between trace points (Compute_Trace), source DB likely incorrect";

static int iter_np(char *A, int M, char *B, int N, Trace_Waves *wave, int mode, int dmax)
{ int  **PVF = wave->PVF; 
  int  **PHF = wave->PHF;
  int    D;
  int    del = M-N;

  { int  *F0, *F1, *F2;
    int  *HF;
    int   low, hgh;
    int   posl, posh;

#ifdef DEBUG_ALIGN
    printf("\n    BASE %ld,%ld: %d vs %d\n",A-wave->Aabs,B-wave->Babs,M,N);
    printf("    A = ");
    for (D = 0; D < M; D++)
      printf("%c",ToA[(int) A[D]]);
    printf("\n");
    printf("    B = ");
    for (D = 0; D < N; D++)
      printf("%c",ToA[(int) B[D]]);
    printf("\n");
#endif

    if (del >= 0)
      { low = 0;
        hgh = del;
      }
    else
      { low = del;
        hgh = 0;
      }

    posl = -dmax;
    posh =  dmax;
    if (wave->Aabs == wave->Babs)
      { if (B == A)
          { EPRINTF(EPLACE,"%s: self comparison starts on diagonal 0 (Compute_Trace)\n",Prog_Name);
            EXIT(-1);
          }
        else if (B < A)
          { if ((B-A)+1 > posl)
              posl = (B-A)+1;
          }
        else
          { if ((B-A)-1 < posh)
              posh = (B-A)-1;
          }
      }

    F1 = PVF[-2];
    F0 = PVF[-1];

    for (D = low-1; D <= hgh+1; D++)
      F1[D] = F0[D] = -2;
    F0[0] = -1;

    low += 1;
    hgh -= 1;

    for (D = 0; 1; D += 1)
      { int   k, i, j;
        int   am, ac, ap;
        char *a;

        if (D > dmax)
          { EPRINTF(EPLACE,"%s: %s\n",Prog_Name,TP_Align);
            EXIT(-1);
          }

        F2 = F1;
        F1 = F0;
        F0 = PVF[D];
        HF = PHF[D];

        if ((D & 0x1) == 0)
          { if (low > posl)
              low -= 1;
            if (hgh < posh)
              hgh += 1;
          }
        F0[hgh+1] = F0[low-1] = -2;

#define FS_MOVE(mdir,pdir)			\
  ac = F1[k]+1;					\
  if (ac < am)					\
    if (ap < am)				\
      { HF[k] = mdir;				\
        j = am;					\
      }						\
    else					\
      { HF[k] = pdir;				\
        j = ap;					\
      }						\
  else						\
    if (ap < ac)				\
      { HF[k] = 0;				\
        j = ac;					\
      }						\
    else					\
      { HF[k] = pdir;				\
        j = ap;					\
      }						\
						\
  if (N < i)					\
    while (j < N && B[j] == a[j])		\
      j += 1;					\
  else						\
    while (j < i && B[j] == a[j])		\
      j += 1;					\
  F0[k] = j;

        j = -2;
        a = A + hgh;
        i = M - hgh;
        for (k = hgh; k > del; k--)
          { ap = j+1;
            am = F2[k-1];
            FS_MOVE(-1,4)
            a -= 1;
            i += 1;
          }

        j = -2;
        a = A + low;
        i = M - low;
        for (k = low; k < del; k++)
          { ap = F2[k+1]+1;
            am = j;
            FS_MOVE(2,1)
            a += 1;
            i -= 1;
          }

        ap = F0[del+1]+1;
        am = j;
        FS_MOVE(2,4)

#ifdef DEBUG_AWAVE
        print_awave(F0,low,hgh);
        print_awave(HF,low,hgh);
#endif

        if (F0[del] >= N)
          break;
      }
  }

  { int   k, h, m, e, c;
    int   ap = (wave->Aabs-A)-1;
    int   bp = (B-wave->Babs)+1;

    PHF[0][0] = 3;

    c = N;
    k = del;
    e = PHF[D][k];
    PHF[D][k] = 3;

    if (mode == UPPERMOST)

      while (e != 3)
        { h = k+e;
          if (e > 1)
            h -= 3;
          else if (e == 0)
            D -= 1;
          else
            D -= 2;

          if (h < k)       // => e = -1 or 2,  UPPERMOST
            { char *a;
  
              a = A + k;
              if (k < 0)
                m = -k;
              else
                m = 0;
              if (PVF[D][h] <= c)
                c = PVF[D][h]-1;
              while (c >= m && a[c] == B[c])
                c -= 1;
              if (e == -1)  //  => edge is 2, others are 1, and 0
                { if (c <= PVF[D+2][k+1])
                    { e = 4;
                      h = k+1;
                      D = D+2;
                    }
                  else if (c == PVF[D+1][k])
                    { e = 0;
                      h = k;
                      D = D+1;
                    }
                  else
                    PVF[D][h] = c+1;
                }
              else      //   => edge is 0, others are 1, and 2 (if k != del), 0 (otherwise)
                { if (k == del)
                    m = D;
                  else
                    m = D-2;
                  if (c <= PVF[m][k+1])
                    { if (k == del)
                        e = 4;
                      else
                        e = 1;
                      h = k+1;
                      D = m;
                    }
                  else if (c == PVF[D-1][k])
                    { e = 0;
                      h = k;
                      D = D-1;
                    }
                  else
                    PVF[D][h] = c+1;
                }
            }

          m = PHF[D][h];
          PHF[D][h] = e;
          e = m;
          k = h;
        }

    else if (mode == LOWERMOST)

      while (e != 3)
        { h = k+e;
          if (e > 1)
            h -= 3;
          else if (e == 0)
            D -= 1;
          else
            D -= 2;

          if (h > k)       // => e = 1 or 4,   LOWERMOST
            { char *a;
  
              a = A + k;
              if (k < 0)
                m = -k;
              else
                m = 0;
              if (PVF[D][h] < c)
                c = PVF[D][h];
              while (c >= m && a[c] == B[c])
                c -= 1;
              if (e == 1)  //  => edge is 2, others are 1, and 0
                { if (c < PVF[D+2][k-1])
                    { e = 2;
                      h = k-1;
                      D = D+2;
                    }
                  else if (c == PVF[D+1][k])
                    { e = 0;
                      h = k;
                      D = D+1;
                    }
                  else
                    PVF[D][h] = c--;
                }
              else      //   => edge is 0, others are 1, and 2 (if k != del), 0 (otherwise)
                { if (k == del)
                    m = D;
                  else
                    m = D-2;
                  if (c < PVF[m][k-1])
                    { if (k == del)
                        e = 2;
                      else
                        e = -1;
                      h = k-1;
                      D = m;
                    }
                  else if (c == PVF[D-1][k])
                    { e = 0;
                      h = k;
                      D = D-1;
                    }
                  else
                    PVF[D][h] = c--;
                }
            }

          m = PHF[D][h];
          PHF[D][h] = e;
          e = m;
          k = h;
        }

    else //  mode == GREEDIEST

      while (e != 3)
        { h = k+e;
          if (e > 1)
            h -= 3;
          else if (e == 0)
            D -= 1;
          else
            D -= 2;

          m = PHF[D][h];
          PHF[D][h] = e;
          e = m;
          k = h;
        }

    k = D = 0;
    e = PHF[D][k];
    while (e != 3)
      { h = k-e;
        c = PVF[D][k];
        if (e > 1)
          h += 3;
        else if (e == 0)
          D += 1;
        else
          D += 2;
#ifdef DEBUG_SCRIPT
        if (h > k)
          printf("     D %d(%d)\n",(c-k)-(ap-1),c+bp);
        else if (h < k)
          printf("     I %d(%d)\n",c+(bp-1),(c+k)-ap);
        else
          printf("     %d S %d\n",(c+k)-(ap+1),c+(bp-1));
#endif
        if (h > k)
          *wave->Stop++ = bp+c;
        else if (h < k)
          *wave->Stop++ = ap-(c+k);
        k = h;
        e = PHF[D][h];
      }
  }

  return (D + abs(del));
}

static int middle_np(char *A, int M, char *B, int N, Trace_Waves *wave, int mode, int dmax)
{ int  **PVF = wave->PVF; 
  int  **PHF = wave->PHF;
  int    D;
  int    del = M-N;

  { int  *F0, *F1, *F2;
    int  *HF;
    int   low, hgh;
    int   posl, posh;

#ifdef DEBUG_ALIGN
    printf("\n%*s BASE %ld,%ld: %d vs %d\n",depth,"",A-wave->Aabs,B-wave->Babs,M,N);
    printf("%*s A = ",depth,"");
    for (D = 0; D < M; D++)
      printf("%c",ToA[(int) A[D]]);
    printf("\n");
    printf("%*s B = ",depth,"");
    for (D = 0; D < N; D++)
      printf("%c",ToA[(int) B[D]]);
    printf("\n");
#endif

    if (del >= 0)
      { low = 0;
        hgh = del;
      }
    else
      { low = del;
        hgh = 0;
      }

    posl = -dmax;
    posh =  dmax;
    if (wave->Aabs == wave->Babs)
      { if (B == A)
          { EPRINTF(EPLACE,"%s: self comparison starts on diagonal 0 (Compute_Trace)\n",Prog_Name);
            EXIT(1);
          }
        else if (B < A)
          { if ((B-A)+1 > posl)
              posl = (B-A)+1;
          }
        else
          { if ((B-A)-1 < posh)
              posh = (B-A)-1;
          }
      }

    F1 = PVF[-2];
    F0 = PVF[-1];

    for (D = low-1; D <= hgh+1; D++)
      F1[D] = F0[D] = -2;
    F0[0] = -1;

    low += 1;
    hgh -= 1;

    for (D = 0; 1; D += 1)
      { int   k, i, j;
        int   am, ac, ap;
        char *a;

        if (D > dmax)
          { EPRINTF(EPLACE,"%s: %s\n",Prog_Name,TP_Align);
            EXIT(-1);
          }

        F2 = F1;
        F1 = F0;
        F0 = PVF[D];
        HF = PHF[D];

        if ((D & 0x1) == 0)
          { if (low > posl)
              low -= 1;
            if (hgh < posh)
              hgh += 1;
          }
        F0[hgh+1] = F0[low-1] = -2;

        j = -2;
        a = A + hgh;
        i = M - hgh;
        for (k = hgh; k > del; k--)
          { ap = j+1;
            am = F2[k-1];
            FS_MOVE(-1,4)
            a -= 1;
            i += 1;
          }

        j = -2;
        a = A + low;
        i = M - low;
        for (k = low; k < del; k++)
          { ap = F2[k+1]+1;
            am = j;
            FS_MOVE(2,1)
            a += 1;
            i -= 1;
          }

        ap = F0[del+1]+1;
        am = j;
        FS_MOVE(2,4)

#ifdef DEBUG_AWAVE
        print_awave(F0,low,hgh);
        print_awave(HF,low,hgh);
#endif

        if (F0[del] >= N)
          break;
      }
  }

  { int   k, h, m, e, c;
    int   d, f;

    d = D + abs(del);
    c = N;
    k = del;

    if (mode == UPPERMOST)

      for (f = d/2; d > f; d--)
        { e = PHF[D][k];
          h = k+e;
          if (e > 1)
            h -= 3;
          else if (e == 0)
            D -= 1;
          else
            D -= 2;

          if (h < k)       // => e = -1 or 2,  UPPERMOST
            { char *a;
  
              a = A + k;
              if (k < 0)
                m = -k;
              else
                m = 0;
              if (PVF[D][h] <= c)
                c = PVF[D][h]-1;
              while (c >= m && a[c] == B[c])
                c -= 1;
              if (e == -1)  //  => edge is 2, others are 1, and 0
                { if (c <= PVF[D+2][k+1])
                    { e = 4;
                      h = k+1;
                      D = D+2;
                    }
                  else if (c == PVF[D+1][k])
                    { e = 0;
                      h = k;
                      D = D+1;
                    }
                  else
                    PVF[D][h] = c+1;
                }
              else      //   => edge is 0, others are 1, and 2 (if k != del), 0 (otherwise)
                { if (k == del)
                    m = D;
                  else
                    m = D-2;
                  if (c <= PVF[m][k+1])
                    { if (k == del)
                        e = 4;
                      else
                        e = 1;
                      h = k+1;
                      D = m;
                    }
                  else if (c == PVF[D-1][k])
                    { e = 0;
                      h = k;
                      D = D-1;
                    }
                  else
                    PVF[D][h] = c+1;
                }
            }

          k = h;
        }

    else if (mode == LOWERMOST)

      for (f = d/2; d > f; d--)
        { e = PHF[D][k];
          h = k+e;
          if (e > 1)
            h -= 3;
          else if (e == 0)
            D -= 1;
          else
            D -= 2;

          if (h > k)       // => e = 1 or 4,   LOWERMOST
            { char *a;
  
              a = A + k;
              if (k < 0)
                m = -k;
              else
                m = 0;
              if (PVF[D][h] < c)
                c = PVF[D][h];
              while (c >= m && a[c] == B[c])
                c -= 1;
              if (e == 1)  //  => edge is 2, others are 1, and 0
                { if (c < PVF[D+2][k-1])
                    { e = 2;
                      h = k-1;
                      D = D+2;
                    }
                  else if (c == PVF[D+1][k])
                    { e = 0;
                      h = k;
                      D = D+1;
                    }
                  else
                    PVF[D][h] = c--;
                }
              else      //   => edge is 0, others are 1, and 2 (if k != del), 0 (otherwise)
                { if (k == del)
                    m = D;
                  else
                    m = D-2;
                  if (c < PVF[m][k-1])
                    { if (k == del)
                        e = 2;
                      else
                        e = -1;
                      h = k-1;
                      D = m;
                    }
                  else if (c == PVF[D-1][k])
                    { e = 0;
                      h = k;
                      D = D-1;
                    }
                  else
                    PVF[D][h] = c--;
                }
            }

          k = h;
        }

    else //  mode == GREEDIEST

      for (f = d/2; d > f; d--)
        { e = PHF[D][k];
          h = k+e;
          if (e > 1)
            h -= 3;
          else if (e == 0)
            D -= 1;
          else
            D -= 2;
          k = h;
        }

    wave->midb = (B-wave->Babs) + PVF[D][k];
    wave->mida = (A-wave->Aabs) + k + PVF[D][k];
  }

  return (0);
}


/****************************************************************************************\
*                                                                                        *
*  COMPUTE_TRACE FLAVORS                                                                 *
*                                                                                        *
\****************************************************************************************/

static char *TP_Error = "Trace point out of bounds (Compute_Trace), source DB likely incorrect";

int Compute_Trace_PTS(Alignment *align, Work_Data *ework, int trace_spacing, int mode)
{ _Work_Data *work = (_Work_Data *) ework;
  Trace_Waves wave;

  Path   *path;
  char   *aseq, *bseq;
  int     alen, blen;
  uint16 *points;
  int     tlen;
  int     ab, bb;
  int     ae, be;
  int     diffs, dmax;

  alen   = align->alen;
  blen   = align->blen;
  path   = align->path;
  aseq   = align->aseq;
  bseq   = align->bseq;
  tlen   = path->tlen;
  points = (uint16 *) path->trace;

  { int64 s;
    int   d;
    int   M, N;
    int   nmax;
    int   **PVF, **PHF;

    M = path->aepos-path->abpos;
    N = path->bepos-path->bbpos;
    if (M < N)
      s = N*sizeof(int);
    else
      s = M*sizeof(int);
    if (s > work->tramax)
      if (enlarge_trace(work,s))
        EXIT(1);

    nmax = 0;
    dmax = 0;
    for (d = 1; d < tlen; d += 2)
      { if (points[d-1] > dmax)
          dmax = points[d-1];
        if (points[d] > nmax)
          nmax = points[d];
      }
    if (tlen <= 1)
      nmax = N;
    if (dmax & 0x1)
      dmax += 1;

    s = (dmax+3)*2*((trace_spacing+nmax+3)*sizeof(int) + sizeof(int *));

    if (s > work->vecmax)
      if (enlarge_vector(work,s))
        EXIT(1);

    wave.PVF = PVF = ((int **) (work->vector)) + 2;
    wave.PHF = PHF = PVF + (dmax+3);

    s = trace_spacing+nmax+3;
    PVF[-2] = ((int *) (PHF + (dmax+1))) + (nmax+1);
    for (d = -1; d <= dmax; d++)
      PVF[d] = PVF[d-1] + s;
    PHF[-2] = PVF[dmax] + s;
    for (d = -1; d <= dmax; d++)
      PHF[d] = PHF[d-1] + s;
  }

  wave.Stop = (int *) (work->trace);
  wave.Aabs = aseq;
  wave.Babs = bseq;

  { int i, d;

    diffs = 0;
    ab = path->abpos;
    ae = (ab/trace_spacing)*trace_spacing;
    bb = path->bbpos;
    tlen -= 2;
    for (i = 1; i < tlen; i += 2)
      { ae = ae + trace_spacing;
        be = bb + points[i];
        if (ae > alen || be > blen)
          { EPRINTF(EPLACE,"%s: %s\n",Prog_Name,TP_Error);
            EXIT(1);
          }
        d  = iter_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave,mode,dmax);
        if (d < 0)
          EXIT(1);
        diffs += d;
        ab = ae;
        bb = be;
      }
    ae = path->aepos;
    be = path->bepos;
    if (ae > alen || be > blen)
      { EPRINTF(EPLACE,"%s: %s\n",Prog_Name,TP_Error);
        EXIT(1);
      }
    d  = iter_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave,mode,dmax);
    if (d < 0)
      EXIT(1);
    diffs += d;
  }

  path->trace = work->trace;
  path->tlen  = wave.Stop - ((int *) path->trace);
  path->diffs = diffs;

  return (0);
}

int Compute_Trace_MID(Alignment *align, Work_Data *ework, int trace_spacing, int mode)
{ _Work_Data *work = (_Work_Data *) ework;
  Trace_Waves wave;

  Path   *path;
  char   *aseq, *bseq;
  int     alen, blen;
  uint16 *points;
  int     tlen;
  int     ab, bb;
  int     ae, be;
  int     diffs, dmax;

  alen   = align->alen;
  blen   = align->blen;
  path   = align->path;
  aseq   = align->aseq;
  bseq   = align->bseq;
  tlen   = path->tlen;
  points = (uint16 *) path->trace;

  { int64 s;
    int   d;
    int   M, N;
    int   nmax;
    int   **PVF, **PHF;

    M = path->aepos-path->abpos;
    N = path->bepos-path->bbpos;
    if (M < N)
      s = N*sizeof(int);
    else
      s = M*sizeof(int);
    if (s > work->tramax)
      if (enlarge_trace(work,s))
        EXIT(1);

    nmax = 0;
    dmax = 0;
    for (d = 1; d < tlen; d += 2)
      { if (points[d-1] > dmax)
          dmax = points[d-1];
        if (points[d] > nmax)
          nmax = points[d];
      }
    if (tlen <= 1)
      nmax = N;
    if (dmax & 0x1)
      dmax += 1;

    s = (dmax+3)*4*((trace_spacing+nmax+3)*sizeof(int) + sizeof(int *));

    if (s > work->vecmax)
      if (enlarge_vector(work,s))
        EXIT(1);

    wave.PVF = PVF = ((int **) (work->vector)) + 2;
    wave.PHF = PHF = PVF + (dmax+3);

    s = trace_spacing+nmax+3;
    PVF[-2] = ((int *) (PHF + (dmax+1))) + (nmax+1);
    for (d = -1; d <= dmax; d++)
      PVF[d] = PVF[d-1] + s;
    PHF[-2] = PVF[dmax] + s;
    for (d = -1; d <= dmax; d++)
      PHF[d] = PHF[d-1] + s;
  }

  wave.Stop = ((int *) work->trace);
  wave.Aabs = aseq;
  wave.Babs = bseq;

  { int i, d;
    int as, bs;
    int af, bf;

    diffs = 0;
    ab = as = af = path->abpos;
    ae = (ab/trace_spacing)*trace_spacing;
    bb = bs = bf = path->bbpos;
    tlen -= 2;
    for (i = 1; i < tlen; i += 2) 
      { ae = ae + trace_spacing;
        be = bb + points[i];
        if (ae > alen || be > blen)
          { EPRINTF(EPLACE,"%s: %s\n",Prog_Name,TP_Error);
            EXIT(1);
          }
        if (middle_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave,mode,dmax))
          EXIT(1);
        af = wave.mida;
        bf = wave.midb;
        d  = iter_np(aseq+as,af-as,bseq+bs,bf-bs,&wave,mode,dmax);
        if (d < 0)
          EXIT(1);
        diffs += d;
        ab = ae;
        bb = be;
        as = af;
        bs = bf;
      }

    ae = path->aepos;
    be = path->bepos;

    if (ae > alen || be > blen)
      { EPRINTF(EPLACE,"%s: %s\n",Prog_Name,TP_Error);
        EXIT(1);
      }
    if (middle_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave,mode,dmax))
      EXIT(1);
    af = wave.mida;
    bf = wave.midb;
    d  = iter_np(aseq+as,af-as,bseq+bs,bf-bs,&wave,mode,dmax);
    if (d < 0)
      EXIT(1);
    diffs += d;
    as = af;
    bs = bf;
    
    d += iter_np(aseq+af,ae-as,bseq+bf,be-bs,&wave,mode,dmax);
    if (d < 0)
      EXIT(1);
    diffs += d;
  }

  path->trace = work->trace;
  path->tlen  = wave.Stop - ((int *) path->trace);
  path->diffs = diffs;

  return (0);
}

int Compute_Trace_IRR(Alignment *align, Work_Data *ework, int mode)
{ _Work_Data *work = (_Work_Data *) ework;
  Trace_Waves wave;

  Path   *path;
  char   *aseq, *bseq;
  int     alen, blen;
  uint16 *points;
  int     tlen;
  int     ab, bb;
  int     ae, be;
  int     diffs, dmax;

  alen   = align->alen;
  blen   = align->blen;
  path   = align->path;
  aseq   = align->aseq;
  bseq   = align->bseq;
  tlen   = path->tlen;
  points = (uint16 *) path->trace;

  { int64 s;
    int   d;
    int   M, N;
    int   mmax, nmax;
    int   **PVF, **PHF;

    M = path->aepos-path->abpos;
    N = path->bepos-path->bbpos;
    if (M < N)
      s = N*sizeof(int);
    else
      s = M*sizeof(int);
    if (s > work->tramax)
      if (enlarge_trace(work,s))
        EXIT(1);

    nmax = mmax = 0;
    for (d = 0; d < tlen; d += 2)
      { if (points[d] > mmax)
          mmax = points[d];
        if (points[d+1] > nmax)
          nmax = points[d+1];
      }
    if (tlen <= 1)
      { mmax = M;
        nmax = N;
      }
    if (mmax > nmax)
      dmax = nmax;
    else
      dmax = mmax;

    s = (dmax+3)*2*((mmax+nmax+3)*sizeof(int) + sizeof(int *));

    if (s > work->vecmax)
      if (enlarge_vector(work,s))
        EXIT(1);

    wave.PVF = PVF = ((int **) (work->vector)) + 2;
    wave.PHF = PHF = PVF + (dmax+3);

    s = mmax+nmax+3;
    PVF[-2] = ((int *) (PHF + (dmax+1))) + (nmax+1);
    for (d = -1; d <= dmax; d++)
      PVF[d] = PVF[d-1] + s;
    PHF[-2] = PVF[dmax] + s;
    for (d = -1; d <= dmax; d++)
      PHF[d] = PHF[d-1] + s;
  }

  wave.Stop = (int *) (work->trace);
  wave.Aabs = aseq;
  wave.Babs = bseq;

  { int i, d;

    diffs = 0;
    ab = path->abpos;
    bb = path->bbpos;
    for (i = 0; i < tlen; i += 2)
      { ae = ab + points[i];
        be = bb + points[i+1];
        if (ae > alen || be > blen)
          { EPRINTF(EPLACE,"%s: %s\n",Prog_Name,TP_Error);
            EXIT(1);
          }
        d = iter_np(aseq+ab,ae-ab,bseq+bb,be-bb,&wave,mode,dmax);
        if (d < 0)
          EXIT(1);
        diffs += d;
        ab = ae;
        bb = be;
      }
  }

  path->trace = work->trace;
  path->tlen  = wave.Stop - ((int *) path->trace);
  path->diffs = diffs;

  return (0);
}

#undef DEBUG_BOX
#undef DEBUG_DP
#undef DEBUG_BACK
#undef BOX_STATS

#define LONG_SNAKE 50

#ifdef DEBUG

static int ASCII[5] = { 'a', 'c', 'g', 't', '.' };

static inline void print_string(char *a, int l)
{ int i;

  for (i = 0; i < l; i++)
    printf("%c",ASCII[(int) a[i]]);
}

#endif

static inline int hamming(char *a, char *b, int n)
{ int h, i, x, y;

  h = 0;
  for (i = 0; i < n; i++)
    { x = *a++;
      if (x == 4)
        break;
      y = *b++;
      if (x != y)
        { if (y == 4)
            break;
          else
            h += 1;
        }
    }
  return (h);
}

static inline int snake(char *a, char *b)
{ int i, x;

  for (i = 0; 1; i++)
    { x = *a++;
      if (x == 4)
        break;
      if (x != *b++)
        break;
    }
  return (i);
}

static inline int rsnake(char *a, char *b)
{ int i, x;

  for (i = 0; 1; i++)
    { x = *--a;
      if (x == 4)
        break;
      if (x != *--b)
        break;
    }
  return (i);
}

#ifdef BOX_STATS

static int   MaxBxArea;
static int   MaxBxWidth;
static int   MaxBxHeight;
static int64 SumBx;
static int   NumBx;
static int   BxHist[101];
static int   BxExtend;
static int   BxGaps;

void BeginBoxStats()
{ int i;

  MaxBxArea   = 0;
  MaxBxWidth  = 0;
  MaxBxHeight = 0;
  SumBx       = 0;
  NumBx       = 0;
  for (i = 0; i <= 100; i++)
    BxHist[i] = 0;
  BxExtend = 0;
  BxGaps   = 0;
}

void EndBoxStats()
{ int i;

  printf("\n# of Boxes = %d with average work %lld\n",NumBx,SumBx/NumBx);
  printf("\nMax Work  = %d\n",MaxBxArea);
  printf("Max Diags = %d\n",MaxBxWidth);
  printf("Max Waves = %d\n",MaxBxHeight);
  printf("\nBox extended = %d\n",BxExtend);
  printf("Gaps removed = %d\n",BxGaps);
  printf("\nHistogram of box work:\n");
  for (i = 0; i <= 100; i++)
    if (BxHist[i] > 0)
      printf(" %3d00: %10d\n",i,BxHist[i]);
}

#endif

void Gap_Improver(Alignment *aln, Work_Data *ework)
{ _Work_Data *work = (_Work_Data *) ework;
  int        *F, *H;
  int        *f, *h;

  char  *A, *B;
  int    x;
  int    p, q;
  int    d, m;
  int   *t, T;
  int    Fpos, Lpos, Fdag, Hamm, Gaps, Diag;
  int    passes;
 
  A = aln->aseq-1;
  B = aln->bseq-1;
  t = (int *) aln->path->trace;
  T = aln->path->tlen;
  F = (int *) work->vector;

  d = aln->path->abpos - aln->path->bbpos;
  q = t[0];
  x = 0;
  while (x < T)
    { p = q;
      m = x;
      Fdag = d;
      Fpos = p;
      Hamm = 0;
      Gaps = 1;
      while (1)
        { x += 1;
          q = 0;
          if (x >= T || (q = t[x]) != p)
            { m = x-m;
              if (p < 0)
                { d -= m;
                  if (q >= 0)
                    break;
                  if (p-q >= LONG_SNAKE)
                    break;
                  Hamm += hamming(A-p,B-(d+p),p-q);
                }
              else
                { d += m;
                  if (q <= 0)
                    break;
                  if (q-p >= LONG_SNAKE)
                    break;
                  Hamm += hamming(A+(p+d),B+p,q-p);
                }
              Gaps += 1;
              p = q;
              m = x;
            }
        }
      if (Gaps == 1)
        continue;
      Lpos = p;
      Diag = abs(Fdag-d)+1;

      // Process box

      p = Diag*(Gaps+Hamm+1)*sizeof(int);
      if (p > work->vecmax)
        { if (enlarge_vector(work,p))
            EXIT (1);
          F = (int *) work->vector;
        }
      H = F + Diag;

#ifdef BOX_STATS
      { int hgt  = Gaps+Hamm+1;
        int area = Diag*hgt;
        if (area > MaxBxArea)
          MaxBxArea = area;
        if (Diag > MaxBxWidth)
          MaxBxWidth = Diag;
        if (hgt > MaxBxHeight)
          MaxBxHeight = hgt;
        NumBx += 1;
        SumBx += area;
        if (area >= 10000)
          BxHist[100] += 1;
        else
          BxHist[area/100] += 1;
      }
#endif
#ifdef DEBUG_BOX
      printf("Box:  %5d :: %4d x %3d (%2d+%2d)  :: %5d .. %5d  %6d .. %6d\n",
             (Gaps+Hamm+1)*Diag,abs(Fpos-Lpos),Diag,Hamm,Gaps,Fpos,Lpos,Fdag,d);
      fflush(stdout);
#endif

      if (Fpos < 0)
        { Fpos = -Fpos;
          Lpos = -Lpos;

          while (A[Fpos-1] != B[(Fpos-Fdag)-1] && A[Fpos-1] != 4 && B[(Fpos-Fdag)-1] != 4)
            { Fpos -= 1;
#ifdef BOX_STATS
              BxExtend += 1;
#endif
            }
          while (A[Lpos] != B[Lpos-d] && A[Lpos] != 4 && B[Lpos-d] != 4)
            { Lpos += 1;
#ifdef BOX_STATS
              BxExtend += 1;
#endif
            }

          f = F;
          *f++ = p = Fpos + snake(A+Fpos,B+(Fpos-Fdag));
          for (m = Fdag-1; m >= d; m--)
            *f++ = Fpos-1;
          passes = 0;

#ifdef DEBUG_DP
          printf(" %2d:",passes);
          for (m = Fdag; m >= d; m--)
            printf(" %d",F[Fdag-m]);
          printf("\n");
          fflush(stdout);
#endif

          h = H;
          p = Fpos;
          while (p < Lpos)
            { int b, c;

              b = Fpos;
              c = 0;
              f = F;
              for (m = Fdag; m >= d; m--)
                { p = b;
                  if (*f >= b)
                    { b = *f;
                      c = 0;
                      p = b+1;
                    }
                  else
                    c += 1;
                  *h++ = c;
                  *f++ = p += snake(A+p,B+(p-m));
                }
              passes += 1;

#ifdef DEBUG_DP
              printf(" %2d:",passes);
              for (m = Fdag; m >= d; m--)
                printf(" %d(%2d)",F[Fdag-m],h[(d-m)-1]);
              printf("\n");
              fflush(stdout);
#endif
            }

          if (passes < Gaps+Hamm)
            { int y, k;

              p = Lpos;
              m = d;
              y = x;
#ifdef DEBUG_BACK
              printf("Short cut %d\n",(Gaps+Hamm)-passes);
              printf("Path (%d,%d)",p,m);
#endif
#ifdef BOX_STATS
              BxGaps += (Gaps+Hamm)-passes;
#endif
              while (h > H)
                { p -= rsnake(A+p,B+(p-m));
                  if (p < Fpos)
                    p = Fpos;
                  h -= Diag;
                  k = h[Fdag-m];
                  if (k == 0)
                    p -= 1;
                  else
                    { m += k;
                      for (; k > 0; k--)
                        t[--y] = -p;
                    }
#ifdef DEBUG_BACK
                  printf(" (%d,%d)",p,m);
#endif
                }
#ifdef DEBUG_BACK
              printf("\n");
#endif
            }
        }
      else
        { while (B[Fpos-1] != A[(Fpos+Fdag)-1] && B[Fpos-1] != 4 && A[(Fpos+Fdag)-1] != 4)
            { Fpos -= 1;
#ifdef BOX_STATS
              BxExtend += 1;
#endif
            }
          while (B[Lpos] != A[Lpos+d] && B[Lpos] != 4 && A[Lpos+d] != 4)
            { Lpos += 1;
#ifdef BOX_STATS
              BxExtend += 1;
#endif
            }

          f = F;
          *f++ = p = Fpos + snake(A+(Fpos+Fdag),B+Fpos);
          for (m = Fdag+1; m <= d; m++)
            *f++ = Fpos-1;
          passes = 0;

#ifdef DEBUG_DP
          printf(" %2d:",passes);
          for (m = Fdag; m <= d; m++) 
            printf(" %d",F[m-Fdag]);
          printf("\n");
          fflush(stdout);
#endif

          h = H;
          p = Fpos;
          while (p < Lpos)
            { int b, c;

              b = Fpos;
              c = 0;
              f = F;
              for (m = Fdag; m <= d; m++) 
                { p = b;
                  if (*f >= b)
                    { b = *f;
                      c = 0;
                      p = b+1;
                    }
                  else
                    c += 1;
                  *h++ = c; 
                  *f++ = p += snake(A+(m+p),B+p);
                }
              passes += 1;

#ifdef DEBUG_DP
              printf(" %2d:",passes);
              for (m = Fdag; m <= d; m++) 
                printf(" %d(%2d)",F[m-Fdag],h[(m-d)-1]);
              printf("\n");
              fflush(stdout);
#endif
            }

          if (passes < Gaps+Hamm)
            { int y, k;

              p = Lpos;
              m = d;
              y = x;
#ifdef DEBUG_BACK
              printf("Short cut %d\n",(Gaps+Hamm)-passes);
              printf("Path (%d,%d)",p,m);
#endif
#ifdef BOX_STATS
              BxGaps += (Gaps+Hamm)-passes;
#endif
              while (h > H)
                { p -= rsnake(A+(p+m),B+p);
                  if (p < Fpos)
                    p = Fpos;
                  h -= Diag;
                  k = h[m-Fdag];
                  if (k == 0)
                    p -= 1;
                  else
                    { m -= k;
                      for (; k > 0; k--)
                        t[--y] = p;
                    }
#ifdef DEBUG_BACK
                  printf(" (%d,%d)",p,m);
#endif
                }
#ifdef DEBUG_BACK
              printf("\n");
#endif
            }
        }
    }
}
