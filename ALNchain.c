/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2024 Chenxi Zhou <chnx.zhou@gmail.com>                          *
 *                                                                               *
 * Permission is hereby granted, free of charge, to any person obtaining a copy  *
 * of this software and associated documentation files (the "Software"), to deal *
 * in the Software without restriction, including without limitation the rights  *
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     *
 * copies of the Software, and to permit persons to whom the Software is         *
 * furnished to do so, subject to the following conditions:                      *
 *                                                                               *
 * The above copyright notice and this permission notice shall be included in    *
 * all copies or substantial portions of the Software.                           *
 *                                                                               *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, *
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE *
 * SOFTWARE.                                                                     *
 *********************************************************************************/

/********************************** Revision History *****************************
 *                                                                               *
 * 22/05/24 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/

#define VERSION "0.1"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "GDB.h"
#include "align.h"
#include "alncode.h"

static char *Usage[] = 
  { "[-v] [-g<int(10000)>] [-l<int(10000)>] [-p<float(0.1)>] [-q<float(0.1)>]",
    "[-z<int(1000)>] [-s<int(10000)>] [-n<int(1)>] [-c<float(0.5)>] [-e<0.0>]",
    "[-f<int(1000)>] [-o<output:path>[.1aln]] <alignments:path>[.1aln]"
  };

#undef DEBUG_CHAIN

static int VERBOSE; // -v

  //  Genome skeletons and auxiliary info

static int ISTWO;                           // If two GDBs are different

static GDB           _AGDB, *AGDB = &_AGDB;  //  1st genome skeleton with parts ...
static GDB           _BGDB, *BGDB = &_BGDB;  //  2nd genome skeleton with parts ...

static int           NASCAFF;     // Number scaffolds
static int           NACONTIG;    // Number contigs
static GDB_SCAFFOLD *ASCAFFS;     // Scaffold record array
static GDB_CONTIG   *ACONTIG;     // Contig record array

static int           NBSCAFF;     // Number scaffolds
static int           NBCONTIG;    // Number contigs
static GDB_SCAFFOLD *BSCAFFS;     // Scaffold record array
static GDB_CONTIG   *BCONTIG;     // Contig record array

static size_t fieldSize[128];

// tree node
typedef struct TNODE
  { struct TNODE *L, *R;
    struct TNODE *next;
    int    bread;
    int    abpos, aepos;
    int    bbpos, bepos;
    int    clen;
    int    active;
    double score;
    int    which;
  } TNODE;

#define INTERNAL 1
#define HEAD     2

typedef struct
 { int beg;
   int end;
 } RANGE;

typedef struct
 { double event;
   int which;
 } ORDER;

static inline int BPOSX(TNODE *node)
{ return (node->abpos); }

static inline int BPOSY(TNODE *node)
{ return (node->bbpos); }

static inline int EPOSX(TNODE *node)
{ return (node->aepos); }

static inline int EPOSY(TNODE *node)
{ return (node->bepos); }

static int (*BPFUNCS[2])(TNODE *) = { BPOSX, BPOSY };

static int (*EPFUNCS[2])(TNODE *) = { EPOSX, EPOSY };

static inline void swap(int *a, int *b)
{ int t = *a; *a = *b; *b = t; }

// insertion sort for smaller arrays
static int partition5(TNODE *nodes, int *order, int low, int high, int (*position)(TNODE *))
{ int i, j, k, p;
  for (i = low + 1; i <= high; i++)
    { k = order[i];
      p = position(nodes + k);
      j = i - 1;
      while (j >= low && position(nodes + order[j]) > p)
        { order[j + 1] = order[j];
          j--;
        }
      order[j + 1] = k;
    }
  return ((low + high) / 2);
}

// place k at the right position
// positions lower than k has smaller values
// positions higher than k has larger values
static int partition(TNODE *nodes, int *order, int low, int high, int k, int (*position)(TNODE *))
{ int i, j, p;
  p = position(nodes + order[k]);
  swap(&order[k], &order[high]);
  i = low;
  for (j = low; j < high; j++)
    { if (position(nodes + order[j]) <= p)
        { swap(&order[i], &order[j]);
          i++;
        }
    }
  swap(&order[i], &order[high]);
  return (i);
}

static int quickSelect(TNODE *nodes, int *order, int low, int high, int k, int (*position)(TNODE *));

#define USE_MEDIAN_MEDIAN

// median of medians
// to guarantee O(n) worst-case
static int selectPivot(TNODE *nodes, int *order, int low, int high, int (*position)(TNODE *))
{
#ifndef USE_MEDIAN_MEDIAN
  // it seems even better than using median of medians?
  // may be a good empirical guess since X-Y correlate?
  return ((low + high) / 2);
#endif
  
  if (high - low < 5)
    return (partition5(nodes, order, low, high, position));
  
  int i, n, m, l, h;
  n = (high - low + 5) / 5;
  for (i = 0, l = low; l <= high; i++, l += 5)
    { h = l + 4;
      if (h > high)
        h = high;
      m = partition5(nodes, order, l, h, position);
      swap(&order[low + i], &order[m]);
    }

  return (quickSelect(nodes, order, low, low + n - 1, low + n / 2, position));
}

// partition by median
static int quickSelect(TNODE *nodes, int *order, int low, int high, int k, int (*position)(TNODE *))
{ int i;
  i = selectPivot(nodes, order, low, high, position);
  i = partition(nodes, order, low, high, i, position);
  if (i == k)
    return (k);
  else if (i > k)
    return (quickSelect(nodes, order, low, i - 1, k, position));
  else
    return (quickSelect(nodes, order, i + 1, high, k, position));
}

TNODE *buildKDTree(TNODE *nodes, int *order, int low, int high, int depth)
{ if (low > high) return (NULL);
  
  int i;
  TNODE *root;

  i = (low + high) >> 1;
  quickSelect(nodes, order, low, high, i, EPFUNCS[depth & 1]);

  root = nodes + order[i];
  root->L = buildKDTree(nodes, order, low, i - 1, depth + 1);
  root->R = buildKDTree(nodes, order, i + 1, high, depth + 1);
 
  return (root);
}

int RORDER(const void *l, const void *r)
{ RANGE *x = (RANGE *) l;
  RANGE *y = (RANGE *) r;
  int xm, ym;

  xm = x->beg;
  ym = y->beg;
  if (xm == ym)
    { xm = x->end;
      ym = y->end;
    }
  return ((xm > ym) - (xm < ym));
}

int rangeSetSize(RANGE *ranges, int n)
{ int i, c;
  c = 0;
  for (i = 0; i < n; i++)
    c += ranges[i].end - ranges[i].beg;
  return (c);
}

#define UNSORTED 0
#define SORTED   1

int mergeRangeFuzzy(RANGE *ranges, int n, int fz, int sorted)
{ if (n <= 0)
    return 0;

  if (!sorted)
    qsort(ranges, n, sizeof(RANGE), RORDER);

  int i, r;
  r = 0;
  for (i = 1; i < n; i++)
    { if (ranges[i].beg <= ranges[r].end + fz)
        { if (ranges[i].end > ranges[r].end)
            ranges[r].end = ranges[i].end;
        }
      else
        ranges[++r] = ranges[i];
    }
  return (r + 1);
}

void mergeSorted(RANGE *aranges, int na, RANGE *branges, int nb)
{ if (na <= 0 || nb <= 0)
    return;

  int i, j, k;

  i = na - 1;
  j = nb - 1;
  k = na + nb - 1;
  while (i >= 0 && j >= 0)
    { if (RORDER(aranges+i, branges+j) > 0)
        aranges[k--] = aranges[i--];
      else
        aranges[k--] = branges[j--];
    }

  while (j >= 0)
    aranges[k--] = branges[j--];
}

int sortedRangeOverlap(RANGE *ranges, int n)
{ int i, ovl, end;
  ovl = 0;
  end = ranges[0].end;
  for (i = 1; i < n; i++)
    { if (ranges[i].beg <= end)
        { if (ranges[i].end > end)
            { ovl += end - ranges[i].beg;
              end  = ranges[i].end;
            }
          else
            ovl += ranges[i].end - ranges[i].beg;
        }
      else
        end = ranges[i].end;
    }
  
  return (ovl); 
}

int sortedRange(RANGE *ranges, int n)
{ if (n <= 0)
    return (1);
  int i;
  for (i = 1; i < n; i++)
    if (RORDER(ranges+i-1, ranges+i) > 0)
      return 0;
  return 1;
}

int mergeSortedRangeFuzzy(RANGE *aranges, int na, RANGE *branges, int nb, int fz, int *len)
{ mergeSorted(aranges, na, branges, nb);
  
  if (len)
    *len = sortedRangeOverlap(aranges, na + nb);
  
  return (mergeRangeFuzzy(aranges, na + nb, fz, SORTED));
}

int TORDER(const void *l, const void *r)
{ TNODE *x = (TNODE *) l;
  TNODE *y = (TNODE *) r;
  int xm, ym;

  xm = x->bread;
  ym = y->bread;
  if (xm == ym)
    { xm = x->abpos;
      ym = y->abpos;
    }
  return ((xm > ym) - (xm < ym));
}

static inline double alnSize(TNODE *node)
{ return ((node->aepos - node->abpos) + (node->bepos - node->bbpos)); }

void KDRangeChain(TNODE *root, TNODE *query, int maxGap, int maxOvl, double penGap, double penOvl, int depth)
{ if (root == NULL || query == NULL) return;

  int axis, rpos, qpos, gap[2], ovl[2], ext[2];
  double score;

  axis = depth & 1;
  rpos = EPFUNCS[axis](root);
  qpos = BPFUNCS[axis](query);
  gap[0] = BPOSX(query) - EPOSX(root);
  gap[1] = BPOSY(query) - EPOSY(root);
  ovl[0] = ovl[1] = 0;
  if (gap[0] < 0)
    { ovl[0] = -gap[0];
      gap[0] = 0;
    }
  if (gap[1] < 0)
    { ovl[1] = -gap[1];
      gap[1] = 0;
    }
  ext[0] = EPOSX(query) - (gap[0]>0?BPOSX(query):EPOSX(root));
  ext[1] = EPOSY(query) - (gap[1]>0?BPOSY(query):EPOSY(root));

  if (root->active && root != query &&
          ext[0] > 0 && ext[1] > 0 &&
          gap[0] <= maxGap && gap[1] <= maxGap && 
          ovl[0] <= maxOvl && ovl[1] <= maxOvl &&
          ovl[0] < EPOSX(query) - BPOSX(query) &&
          ovl[1] < EPOSY(query) - BPOSY(query))
    { score = ext[0] + ext[1] 
        - (double) gap[0] * penGap - (double) gap[1] * penGap 
        - (double) ovl[0] * penOvl - (double) ovl[1] * penOvl;
      //if (score > 0 && root->score + score > query->score)
      if (root->score + score > query->score)
        { query->next  = root;
          query->clen  = root->clen + 1;
          query->score = root->score + score;
        }
    }

  if (maxOvl == INT32_MAX || qpos - maxOvl <= rpos)
    KDRangeChain(root->L, query, maxGap, maxOvl, penGap, penOvl, depth + 1);
  if (maxOvl == INT32_MAX || qpos + maxGap >= rpos)
    KDRangeChain(root->R, query, maxGap, maxOvl, penGap, penOvl, depth + 1);
}

static int SORDER(const void *l, const void *r)
{ ORDER *x = (ORDER *) l;
  ORDER *y = (ORDER *) r;
  return ((x->event < y->event) - (x->event > y->event));
}

void backtrackLocal(TNODE *node, int maxDrop, double penGap, double penOvl)
{ if (node->active)
    return;
  
  double minScore, score;
  int clen, gap[2], ovl[2], ext[2];
  TNODE *head, *next;

  head = node;
  minScore = node->score;
  head->active = HEAD;
  next = node->next;
  while (next)
    { if (next->active ||
            next->score > maxDrop + minScore)
        { // break chain
          node->next = NULL;
          break;
        }
      if (next->score < minScore)
        minScore = next->score;
      next->active = INTERNAL;
      node = next;
      next = node->next;
    }

  // recalculate chain score
  node = head;
  score = alnSize(node);
  next = node->next;
  clen = 1;
  while (next)
    { gap[0] = BPOSX(node) - EPOSX(next);
      gap[1] = BPOSY(node) - EPOSY(next);
      ovl[0] = ovl[1] = 0;
      if (gap[0] < 0)
        { ovl[0] = -gap[0];
          gap[0] = 0;
        }
      if (gap[1] < 0)
        { ovl[1] = -gap[1];
          gap[1] = 0;
        }
      ext[0] = (gap[0]>0?EPOSX(next):BPOSX(node)) - BPOSX(next);
      ext[1] = (gap[1]>0?EPOSY(next):BPOSY(node)) - BPOSY(next);

      score += ext[0] + ext[1]
          - (double) gap[0] * penGap - (double) gap[1] * penGap
          - (double) ovl[0] * penOvl - (double) ovl[1] * penOvl;

      node = next;
      next = node->next;
      clen++;
    }

  head->score = score;
  head->clen  = clen;
}

static int popLocalChain(TNODE *nodes, int acnt, double minScore, int maxDrop, int minFrag, 
        double penGap, double penOvl)
{ int i, nchain;
  ORDER *order;
  
  order = (ORDER *)Malloc(sizeof(ORDER)*acnt,"Allocating order array");
  if (order == NULL)
    exit (1);

  for (i = 0; i < acnt; i++)
    { order[i].event = nodes[i].score;
      order[i].which = i;
      nodes[i].active = 0;
    }

  qsort(order, acnt, sizeof(ORDER), SORDER);

  for (i = 0; i < acnt; i++)
    backtrackLocal(nodes + order[i].which, maxDrop, penGap, penOvl);

  minScore *= 2; // chain score is calculated in both X and Y
  nchain = 0;
  for (i = 0; i < acnt; i++)
    { if (nodes[i].active != HEAD)
        continue;
      if (nodes[i].score < minScore || nodes[i].clen < minFrag)
        { nodes[i].active = 1;
          continue;
        }
      nchain++;
    }

#ifdef DEBUG_CHAIN
  fprintf(stderr, "%s: %d local chains were identified\n",Prog_Name,nchain);
#endif

  free(order);

  return (nchain);
}

int localChain(TNODE *nodes, int acnt, int *order, int maxGap, int maxOvl, 
        double penGap, double penOvl, int maxDrop, int minFrag, double minScore)
{ int i, nchain;
  TNODE *root, *node;

  for (i = 0; i < acnt; i++)
    order[i] = i;

  root = buildKDTree(nodes, order, 0, acnt - 1, 0);

  for (i = 0; i < acnt; i++)
    { node = nodes + i;
      KDRangeChain(root, node, maxGap, maxOvl, penGap, penOvl, 0);
      node->active = INTERNAL;
    }

  nchain = popLocalChain(nodes, acnt, minScore, maxDrop, minFrag, penGap, penOvl);

  return (nchain);
}

void reverseRangeStrand(RANGE *ranges, int n, int len)
{ int i, t;
  for (i = 0; i < n; i ++)
    { t = ranges[i].beg;
      ranges[i].beg = len - ranges[i].end;
      ranges[i].end = len - t;
    }
}

int filterChain(TNODE *nodes, int acnt, int alen, int blen, double maxCov, double minExt, int fzMerge)
{ int i, rev, nchain, nfilter, nseg, xmseg, ymseg, xnseg, ynseg, xlen, xcov, ylen, ycov, xext, yext;
  RANGE *xranges, *yranges, *xrange1, *yrange1, *range;
  ORDER *order;
  TNODE *node;

  nchain = 0;
  for (i = 0; i < acnt; i++)
    if (nodes[i].active == HEAD)
      nchain++;
  
  if (nchain == 0)
    return (0);

  order = (ORDER *)Malloc(sizeof(ORDER)*nchain,"Allocating order array");
  xranges = (RANGE *)Malloc(sizeof(RANGE)*acnt*4,"Allocating range space");
  if (order == NULL || xranges == NULL)
    exit (1);
  xrange1 = xranges + acnt;
  yranges = xrange1 + acnt;
  yrange1 = yranges + acnt;

  xext = (double) alen * minExt;
  yext = (double) blen * minExt;
  
  nchain = 0;
  for (i = 0; i < acnt; i++)
    { if (nodes[i].active == HEAD)
        { order[nchain].event = nodes[i].score;
          order[nchain].which = i;
          nchain++;
        }
    }

  qsort(order, nchain, sizeof(ORDER), SORDER);
  
  // primary chain
  node = nodes + order[0].which;
  rev  = node->bread & 1;
  nseg = 0;
  while (node)
    { range = xranges + nseg;
      range->beg = node->abpos;
      range->end = node->aepos;
      range = yranges + nseg;
      range->beg = node->bbpos;
      range->end = node->bepos;
      nseg++;
      node = node->next;
    }
  if (rev)
    reverseRangeStrand(yranges, nseg, blen);
  xmseg = mergeRangeFuzzy(xranges, nseg, fzMerge, UNSORTED);
  ymseg = mergeRangeFuzzy(yranges, nseg, fzMerge, UNSORTED);
  
  nfilter = 0;
  for (i = 1; i < nchain; i++)
    { xcov = ycov = 0;
      
      node = nodes + order[i].which;
      nseg = 0;
      while (node)
        { range = xrange1 + nseg;
          range->beg = node->abpos;
          range->end = node->aepos;
          nseg++;
          node = node->next;
        }
      nseg  = mergeRangeFuzzy(xrange1, nseg, 0, UNSORTED);
      xlen  = rangeSetSize(xrange1, nseg);
      xnseg = mergeSortedRangeFuzzy(xrange1, nseg, xranges, xmseg, fzMerge, &xcov);

      node = nodes + order[i].which;
      rev  = node->bread & 1;
      nseg = 0;
      while (node)
        { range = yrange1 + nseg;
          range->beg = node->bbpos;
          range->end = node->bepos;
          nseg++;
          node = node->next;
        }
      if (rev)
        reverseRangeStrand(yrange1, nseg, blen);
      nseg  = mergeRangeFuzzy(yrange1, nseg, 0, UNSORTED);
      ylen  = rangeSetSize(yrange1, nseg);
      ynseg = mergeSortedRangeFuzzy(yrange1, nseg, yranges, ymseg, fzMerge, &ycov);

      if ((xcov > xlen * maxCov && ycov > ylen * maxCov) || 
          (xlen - xcov < xext && ylen - ycov < yext))
        { nodes[order[i].which].active = INTERNAL;
          nfilter++;

#ifdef DEBUG_CHAIN
          fprintf(stderr, "%s: chain %d filtered out: X=%d/%d=%.3f Y=%d/%d=%.3f\n",
                  Prog_Name,order[i].which,xcov,xlen,(double)xcov/xlen,xcov,ylen,(double)ycov/ylen);
#endif
        }
      else
        { xmseg = xnseg;
          memcpy(xranges, xrange1, sizeof(RANGE) * xmseg);
          ymseg = ynseg;
          memcpy(yranges, yrange1, sizeof(RANGE) * ymseg);
#ifdef DEBUG_CHAIN
          fprintf(stderr, "%s: chain %d kept: X=%d/%d=%.3f Y=%d/%d=%.3f\n",
                  Prog_Name,order[i].which,xcov,xlen,(double)xcov/xlen,ycov,ylen,(double)ycov/ylen);
#endif
        }
    }

#ifdef DEBUG_CHAIN
  fprintf(stderr, "%s: %d non-primary chains were filtered out\n",Prog_Name,nfilter);
#endif

  free(order);
  free(xranges);
  
  return (nfilter);
}

int main(int argc, char *argv[])
{ Overlap   _ovl, *ovl = &_ovl;
  int       amax;
  char     *CommandLine;

  OneFile  *input, *output;
  int64      novl;
  
  int       maxGap   = 10000;
  int       maxOvl   = 10000;
  double    penGap   = .10;
  double    penOvl   = .10;
  double    maxCov   = .50;
  double    minExt   = .0;
  int       minScore = 10000;
  int       maxDrop  = 1000;
  int       fzMerge  = 1000;
  int       minFrag  = 1;
  char     *outfile  = NULL;

  { int   n, i;
    char *c;

    n = 1; // need to be at least one to accommodate final '\0'
    for (i = 1; i < argc; i++)
      n += strlen(argv[i])+1;

    CommandLine = Malloc(n+1,"Allocating command string");
    if (CommandLine == NULL)
      exit (1);

    c = CommandLine;
    if (argc > 1)
      { c += sprintf(c,"%s",argv[1]);
        for (i = 2; i < argc; i++)
          c += sprintf(c," %s",argv[i]);
      }
    *c = '\0';
  }

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("ALNchain")

    (void) flags;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'g':
            ARG_NON_NEGATIVE(maxGap,"Max gap size");
            break;
          case 'l':
            ARG_NON_NEGATIVE(maxOvl,"Max overlap size");
            break;
          case 'p':
            ARG_REAL(penGap);
            break;
          case 'q':
            ARG_REAL(penOvl);
            break;
          case 'c':
            ARG_REAL(maxCov);
            break;
          case 'e':
            ARG_REAL(minExt);
            break;
          case 's':
            ARG_NON_NEGATIVE(minScore,"Min chain score")
            break;
          case 'n':
            ARG_NON_NEGATIVE(minFrag,"Min alignment fragment")
            break;
          case 'z':
            ARG_NON_NEGATIVE(maxDrop,"Max chain score drop")
            break;
          case 'f':
            ARG_NON_NEGATIVE(fzMerge,"Max gap for fuzzy merge")
            break;
          case 'o':
            outfile = argv[i] + 2;
            if (*outfile == '\0')
              { fprintf(stderr,"%s: Empty output file name\n",Prog_Name);
                 exit (1);
              }
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE  = flags['v'];

    if (argc < 2)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -g: maximum gap size\n");
        fprintf(stderr,"      -l: maximum overlap size\n");
        fprintf(stderr,"      -p: a gap of size G cost (-p)*G\n");
        fprintf(stderr,"      -q: an overlap of size O cost (-q)*O\n");
        fprintf(stderr,"      -z: score drop threshold for breaking a chain\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -s: minimum chain score\n");
        fprintf(stderr,"      -n: minimum number of alignment fragments in a chain\n");
        fprintf(stderr,"      -c: maximum coverage as a fraction of chain size\n");
        fprintf(stderr,"      -e: minimum extension as a fraction of sequence size\n");
        fprintf(stderr,"      -f: maximum gap for fuzzy merge\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -o: 1-code output file name\n");
        fprintf(stderr,"      -v: verbose mode\n");
        fprintf(stderr,"\n");

        exit (1);
      }
  }

  //  Initiate .1aln file reading and read header information

  { char      *path, *root, *pwd, *cpath;
    FILE      *test;
    char      *src1_name, *src2_name;
    int        TSPACE;
    
    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".1aln");
    path  = Catenate(pwd,"/",root,".1aln");
    input = open_Aln_Read(path,1,&novl,&TSPACE,&src1_name,&src2_name,&cpath);
    if (input == NULL)
      { fprintf(stderr,"%s: Could not open .1aln file: %s\n",Prog_Name,path);
        exit (1);
      }
    if (novl == 0)
      { fprintf(stderr, "%s: Empty .1aln file: %s\n",Prog_Name,path);
        exit (1);
      }
    if (outfile == NULL)
      path = Catenate(pwd,"/",root,".chain.1aln");
    else
      { free(pwd);
        free(root);
        pwd  = PathTo(outfile);
        root = Root(outfile,".1aln");
        path = Catenate(pwd,"/",root,".1aln");
      }
    test = fopen(path,"w");
    if (test == NULL)
      { fprintf(stderr,"%s: Cannot open %s for output\n",Prog_Name,path);
        exit (1);
      }
    fclose(test);
    output = oneFileOpenWriteFrom(path,input,true,1);
    if (output == NULL)
      { fprintf(stderr,"%s: Could not open .1aln file to write: %s\n",Prog_Name,path);
        exit (1);
      }
    oneAddProvenance(output,Prog_Name,VERSION,CommandLine);
    free(pwd);
    free(root);

    ISTWO = (src2_name != NULL);

    if (input->lineType == 'g')
      Read_Aln_Skeleton(input,src1_name,AGDB);
    else
      Get_GDB(AGDB,src1_name,cpath,0);
    if (ISTWO)
      { if (input->lineType == 'g')
          Read_Aln_Skeleton(input,src2_name,BGDB);
        else
          Get_GDB(BGDB,src2_name,cpath,0);
      }
    else
      BGDB = AGDB;

    free(cpath);
    free(src1_name);
    free(src2_name);

    NASCAFF  = AGDB->nscaff;
    NACONTIG = AGDB->ncontig;
    ASCAFFS  = AGDB->scaffolds;
    ACONTIG  = AGDB->contigs;

    NBSCAFF  = BGDB->nscaff;
    NBCONTIG = BGDB->ncontig;
    BSCAFFS  = BGDB->scaffolds;
    BCONTIG  = BGDB->contigs;
  }

  //  Get field sizes

  { int i;
    for (i = 0; i < 128;++i)
      if (input->info[i] != NULL)
        fieldSize[i] = input->info[i]->nField*sizeof(OneField);
  }

  //  Write global lines before 'A'
  //  Input will be placed at the first 'A' line

  { oneGoto(input,'A',0); // need this as open_Aln_Read advances to the first 'A' line
    while (oneReadLine(input) && input->lineType != 'A')
      { memcpy(output->field,input->field,fieldSize[(int) input->lineType]);
        oneWriteLine(output,input->lineType,oneLen(input),oneString(input));
      }
  }

  //  Read the file to get maximums
  
  { int i;
    int aread, ascaf, acnt;

    amax  = 0;
    acnt  = 0;
    aread = -1;
    ascaf = -1;
    for (i = 0; i < novl; i++)
      { Read_Aln_Overlap(input,ovl);
        Skip_Aln_Trace(input);
        if (aread != ovl->aread && ascaf != ACONTIG[ovl->aread].scaf)
          { if (acnt > amax)
              amax = acnt;
            acnt  = 0;
            ascaf = ACONTIG[ovl->aread].scaf;
          }
        aread = ovl->aread;
        acnt += 1;
      }
    if (acnt > amax)
      amax = acnt;
    if (VERBOSE)
      fprintf(stderr,"%s: There are %d alignments for largest scaffold\n", Prog_Name, amax);
  }

  //  Read the file and chain
  
  { int      i, j, k, b;
    int      aread, acnt, alen, ascaf;
    int      bread, boff, blen;
    int      apulse, bpulse;
    int      fcnt;
    int     *order, nchain, nalign;
    TNODE   *nodes, *node;

    oneGoto(input,'A',1);
    oneReadLine(input);

    nodes = (TNODE *) Malloc(sizeof(TNODE)*amax,"Allocating alignment space");
    order = (int *) Malloc(sizeof(int)*amax*2,"Allocating order arrays");
    if (nodes == NULL || order == NULL)
      exit (1);

    nchain =  0;
    nalign =  0;
    acnt   =  0;
    fcnt   =  0;
    aread  = -1;
    ascaf  = -1;
    for (i = 0; i <= novl; i++)
      { if (i < novl)
          { Read_Aln_Overlap(input,ovl);
            Skip_Aln_Trace(input);
          }

        if (acnt > 0 && (ACONTIG[ovl->aread].scaf != ascaf || i == novl))
          { alen  = ASCAFFS[ascaf].slen;

            qsort(nodes, acnt, sizeof(TNODE), TORDER);
            
            { for (j = 0, k = 0, b = nodes->bread; j <= acnt; j++)
                { if (j == acnt || nodes[j].bread != b)
                    { localChain(nodes + k, j - k, order, maxGap, maxOvl, penGap, penOvl,
                                 maxDrop, minFrag, minScore);
                      if (j == acnt)
                        break;
                      k = j;
                      b = nodes[j].bread;
                    }
                }

              for (j = 0, k = 0, b = nodes->bread >> 1; j <= acnt; j++)
                { if (j == acnt || (nodes[j].bread >> 1) != b)
                    { filterChain(nodes + k, j - k, alen, BSCAFFS[b].slen, maxCov, minExt,
                                  fzMerge);
                      if (j == acnt)
                        break;
                      k = j;
                      b = nodes[j].bread >> 1;
                    }
                }
            }

#ifdef DEBUG_CHAIN
            { // some stats for the chaining results
              int acov = 0;
              for (j = 0; j < acnt; j++)
                { if (nodes[j].active != HEAD)
                    continue;
                  node  = nodes + j;
                  acov += node->aepos - node->abpos;
                  while (node != NULL)
                    { acov += node->aepos - node->abpos;
                      node  = node->next;
                    }
                }
              fprintf(stderr,"%s: scaffold %6d - %d/%d %.3f\n",Prog_Name,ascaf,acov,alen,
                             (double)acov/alen);
            }
#endif

            { // mark aln records
              for (j = 0; j < acnt; j++)
                if (nodes[j].active != HEAD)
                  nodes[j].active = 0;
              for (j = 0; j < acnt; j++)
                { if (nodes[j].active != HEAD)
                    continue;
                  nchain++;
                  nalign++;
                  node = nodes[j].next;
                  while (node != NULL)
                    { node->active = INTERNAL;
                      node = node->next;
                      nalign++;
                    }
                }
              
              // write 1aln file for chained records
              for (j = 0; j < acnt; j++)
                { if (!nodes[j].active)
                    continue;
                  if (!oneGoto(input, 'A', nodes[j].which+1))
                    { fprintf(stderr,"%s: Can't locate to object %d in aln file\n",
                                     Prog_Name,nodes[j].which+1);
                      exit (1);
                    }
                  //  Copy 'A' line
                  oneReadLine(input);
                  memcpy(output->field,input->field,fieldSize[(int) input->lineType]);
                  oneWriteLine(output,input->lineType,oneLen(input),oneString(input));
                  //  Copy all lines before reaching next 'A' line or EOF
                  while (oneReadLine(input) && input->lineType != 'A')
                    { memcpy(output->field,input->field,fieldSize[(int) input->lineType]);
                      oneWriteLine(output,input->lineType,oneLen(input),oneString(input));
                    }
                }
            }

            if (i >= novl)
              break;
            fcnt += acnt;
            acnt  = 0;

            if (!oneGoto(input, 'A', fcnt+1))
              { fprintf(stderr,"%s: Can't locate to object %d in aln file\n",Prog_Name,fcnt);
                exit (1);
              }
            oneReadLine(input);
            Read_Aln_Overlap(input,ovl);
            Skip_Aln_Trace(input);
          }

        aread  = ovl->aread;
        bread  = ovl->bread;
        ascaf  = ACONTIG[aread].scaf;
        apulse = ACONTIG[aread].sbeg;
        bpulse = BCONTIG[bread].sbeg;
        b = BCONTIG[bread].scaf << 1;
        if (COMP(ovl->flags))
          b |= 1;
        node = nodes + acnt++;
        node->bread = b;
        node->abpos = ovl->path.abpos + apulse;
        node->aepos = ovl->path.aepos + apulse;
        if (COMP(ovl->flags))
          { boff = bpulse + BCONTIG[bread].clen;
            blen = BSCAFFS[BCONTIG[bread].scaf].slen;
            node->bbpos = blen - (boff - ovl->path.bbpos);
            node->bepos = blen - (boff - ovl->path.bepos);
          }
        else
          { node->bbpos = ovl->path.bbpos + bpulse;
            node->bepos = ovl->path.bepos + bpulse;
          }
        node->L = NULL;
        node->R = NULL;
        node->next = NULL;
        node->clen = 1;
        node->active = 0;
        node->score = alnSize(node);
        node->which = i;
      }
    
    fprintf(stderr,"%s: retained %d alignments in %d chains\n",Prog_Name,nalign,nchain);

    free(nodes);
    free(order);
  
    oneFileClose(input);
    oneFileClose(output);
  }

  Close_GDB(AGDB);
  if (ISTWO)
    Close_GDB(BGDB);

  free(CommandLine);
  free(Command_Line);
  free(Prog_Name);

  Catenate(NULL,NULL,NULL,NULL);

  exit (0);
}
