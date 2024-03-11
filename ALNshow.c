/*******************************************************************************************
 *
 *  Utility for displaying the overlaps in a .1aln file in a variety of ways including
 *    a minimal listing of intervals, and a full out alignment.
 *
 *  Author:    Gene Myers
 *             Modified by Richard Durbin to work on .1aln files
 *  Creation:  July 2013
 *  Last Mod:  March 2023
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/uio.h>
#include <fcntl.h>

#include "DB.h"
#include "align.h"
#include "alncode.h"

static char *Usage[] =
    { "[-arU] [-i<int(4)>] [-w<int(100)>] [-b<int(10)>] ",
      "    <alignments:path>[.1aln] [ <range> [<range>] ]"
    };

#define LAST_SYMBOL  '#'
#define SCAF_SYMBOL  '@'

static char      *AHEADERS;
static DAZZ_READ *AREADS;
static int       *AMAP;
static int       *HPERMA;

static char      *BHEADERS;
static DAZZ_READ *BREADS;
static int       *BMAP;
static int       *HPERMB;

//  Sorted table of Headers

static int HSORTA(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  return (strcmp(AHEADERS+AREADS[AMAP[x]].coff,AHEADERS+AREADS[AMAP[y]].coff));
}

static int HSORTB(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  return (strcmp(BHEADERS+BREADS[BMAP[x]].coff,BHEADERS+BREADS[BMAP[y]].coff));
}

static char *AHEAD(int x)
{ return (AHEADERS+AREADS[AMAP[HPERMA[x]]].coff); }

static char *BHEAD(int x)
{ return (BHEADERS+BREADS[BMAP[HPERMB[x]]].coff); }

static int lookup(char *x, char *(*head)(int), int nscaff)
{ int l, r, m;

  // smallest l s.t. HPERM(l) >= x (or NSCAFF if does not exist)

  l = 0;
  r = nscaff;
  while (l < r)
    { m = ((l+r) >> 1);
      if (strcmp(head(m),x) < 0)
        l = m+1;
      else
        r = m;
    }

  if (l >= nscaff || strncmp(head(l),x,strlen(x)) != 0)
    return (-1);

  return (l);
}

static int matches(char *x, int *pe, char *(*head)(int), int nscaff)
{ int b, e;

  b = lookup(x,head,nscaff);
  for (e = b+1; b < nscaff; e++)
    if (strncmp(head(e),x,strlen(x)) != 0)
      break;
  *pe = e-1;
  return (b);
}


//  Parse read range grabbing up to 4 values from it as follows:
//    type 1
//          val[0][.val[1]]
//    type 2
//          val[0][.val[1]] - val[2][.val[3]]
//    type 3
//          val[0][.val[1]] _ val[2] - val[3]
//  Return 0 if there is a syntactic error, otherwise return the type.
//  0 values indicate $, -1 indicates not present, -2 (only in val[1] or val[3])
//      indicates val[0] or val[2], respectively, is a scaffold index, -3 (like -2)
//      implies val is index into arg of string value.

static char *white(char *x)
{ while (isspace(*x))
    x += 1;
  return (x);
}

static char *getint(char *x, int *v)
{ *v = 0;
  while (isdigit(*x))
    *v = 10*(*v) + (*x++ - '0');
  return (x);
}

static char *address(char *x, int *v, int sep)
{ int a;

  x = white(x);
  a = *x++;
  if (a == LAST_SYMBOL)
    v[0] = 0;
  else if (isdigit(a))
    x = getint(x-1,v);
  else
    return (NULL);
  x = white(x);
  if (*x == sep)
    { x = white(x+1);
      a = *x++;
      if (a == LAST_SYMBOL)
        v[1] = 0;
      else if (isdigit(a))
        x = getint(x-1,v+1);
      else
        return (NULL);
      x = white(x);
    }
  else if (sep == ' ')
    v[1] = -2;
  else if (sep == '.')
    v[1] = -1;
  else  // sep == '-'
    return (NULL);
  return (x);
}

static int range(char *src, int *v, char *(*head)(int), int nscaff)
{ int   t, sep;
  char *x, *y, *z;

  v[2] = -1;
  x = white(src);
  if (*x == SCAF_SYMBOL)
    { x += 1;
      sep = ' ';
    }
  else
    sep = '.';
  y = address(x,v,sep);
  if (y == NULL)
    { if (sep == ' ')
        goto try_string;
      else
        return (0);
    }
  t = 1;
  if (*y == '-')
    { y = address(y+1,v+2,sep);
      t = 2;
    }
  else if (*y == '_')
    { y = address(y+1,v+2,'-');
      t = 3;
    }
  if (y == NULL || *y != '\0')
    { if (sep == ' ')
        goto try_string;
      else
        return (0);
    }
  return (t);

try_string:
  if (lookup(x,head,nscaff) >= 0)
    { v[1] = -3;
      v[0] = x-src;
      return (1);
    }
  y = rindex(x,'_');
  if (y != NULL)
    { *y = '\0';
      if (lookup(x,head,nscaff) >= 0)
        { z = address(y+1,v+2,'-');
          if (z != NULL && *z == '\0')
            { v[1] = -3;
              v[0] = x-src;  
              return (3);
            }
        }
      *y = '_';
    }
  for (y = index(x,'-'); y != NULL; y = index(y+1,'-'))
    { *y = '\0'; 
      if (lookup(x,head,nscaff) >= 0 && lookup(y+1,head,nscaff) >= 0)
        { v[1] = v[3] = -3;
          v[0] = x-src;
          v[2] = (y+1)-src;
          return (2);
        }
      *y = '-';
    }
  return (0);
}

 //  Interpret value pair s.c into an absolute contig range (1-based)

static int interpret_address(int s, int c, int ncontig, int nscaff, int *map)
{ if (c == -1)
    { if (s == 0)
        return (ncontig);
      if (s > ncontig)
        { fprintf(stderr,"%s: Absolute contig index '%d' is out of bounds\n",Prog_Name,s);
          exit (1);
        }
      return (s);
    }
  if (s == 0)
    s = nscaff;
  else if (s > nscaff)
    { fprintf(stderr,"%s: Scaffold index '%d' is out of bounds\n",Prog_Name,s);
      exit (1);
    }
  if (c == 0)
    return (map[s]);
  if (c == -2)
    return (s);
  if (c > map[s]-map[s-1])
    { fprintf(stderr,"%s: Relative contig index '%d' is out of bounds\n",Prog_Name,c);
      exit (1);
    }
  return (map[s-1] + c);
}

  //  Interpret argument x as a range and map to contig level ranges

static void contig_range(char *x, DAZZ_DB *db, int ncontig, int nscaff, int *map, char *(*head)(int), int *hperm,
                         int *pbeg, int *pend, int *pfst, int *plst)
{ int t, vals[4], len, sunit;
  int beg, end, fst, lst;

  t = range(x,vals,head,nscaff);
  if (t == 0)
    { fprintf(stderr,"%s: Command line argument %s not a valid read range\n",Prog_Name,x);
      exit (1);
    }

  sunit = (vals[1] <= -2);

  if (vals[1] == -3)
    { beg = matches(x+vals[0],&end,head,nscaff);
      if (t == 2)
        { matches(x+vals[2],&end,head,nscaff);
          if (beg > end)
            { fprintf(stderr,"%s: Scaffold range in '%s' is empty\n",Prog_Name,x);
              exit (1);
            }
        }
      if (t == 3)
        { if (end > beg)
            { fprintf(stderr,"%s: String in %s matches multiple scaffolds, so substring ambiguous.\n",Prog_Name,x);
              exit (1);
            }
          beg = end = hperm[beg];
        }
      else
        { *pbeg = beg;
          *pend = end;
          *pfst = -1;
          return;
        }
    }
  else
    { beg = end = interpret_address(vals[0],vals[1],ncontig,nscaff,map) - 1;
      fst = 0;
      if (t == 2)
        { end = interpret_address(vals[2],vals[3],ncontig,nscaff,map) - 1;
          if (beg > end)
            { fprintf(stderr,"%s: Read range in '%s' is empty\n",Prog_Name,x);
              exit (1);
            }
        }
    }

  if (sunit)
    { beg = map[beg];
      end = map[end+1]-1; 
    }
  if (t <= 2)
    lst = db->reads[end].rlen;
  else // t == 3
    { fst = vals[2]; //  Type 3: intepret vals[2] and vals[3] as the substring interval
      lst = vals[3];
      if (sunit)
        len = db->reads[end].fpulse + db->reads[end].rlen;
      else
        len = db->reads[beg].rlen;
      if (lst == 0)
        lst = len;
      else if (lst > len)
        { fprintf(stderr,"%s: Substring interval in '%s' is out of bounds\n",Prog_Name,x);
          exit (1);
        }
      if (fst >= lst)
        { fprintf(stderr,"%s: Substring interval in '%s' is empty\n",Prog_Name,x);
          exit (1);
        }
      if (sunit)
        { for ( ; beg <= end; beg++)
            if (fst < db->reads[beg].fpulse + db->reads[beg].rlen)
              break;
          fst -= db->reads[beg].fpulse; 
          if (fst < 0)
            fst = 0;
          for (; end >= beg; end--)
            if (lst > db->reads[end].fpulse)
              break;
          lst -= db->reads[end].fpulse;
          if (lst > db->reads[end].rlen)
            lst = db->reads[end].rlen;
        }
    }
  *pbeg = beg;
  *pend = end;
  *pfst = fst;
  *plst = lst;
}


int main(int argc, char *argv[])
{ DAZZ_DB   _db1, *db1 = &_db1; 
  DAZZ_DB   _db2, *db2 = &_db2; 
  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;

  char   *db1_name, *db2_name;
  OneFile *input;
  int64   novl;
  int     tspace;

  int     ALIGN, REFERENCE;
  int     INDENT, WIDTH, BORDER, UPPERCASE;
  int     ISTWO;

  int    *amap, *alen, *actg, *ascf, nascaff, nacontig, amaxlen, actgmax;
  int    *bmap, *blen, *bctg, *bscf, nbscaff, nbcontig, bmaxlen, bctgmax;
  int     abeg, aend, afst, alst;
  int     bbeg, bend, bfst, blst;
  char   *abst, *aest;
  char   *bbst, *best;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("ALNshow")

    INDENT    = 4;
    WIDTH     = 100;
    BORDER    = 10;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("arU")
            break;
          case 'i':
            ARG_NON_NEGATIVE(INDENT,"Indent")
            break;
          case 'w':
            ARG_POSITIVE(WIDTH,"Alignment width")
            break;
          case 'b':
            ARG_NON_NEGATIVE(BORDER,"Alignment border")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    ALIGN     = flags['a'];
    REFERENCE = flags['r'];
    UPPERCASE = flags['U'];

    if (argc < 2 || argc > 4)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"     <range> =    <contig>[-<contig>]     |  <contig>_<int>-(<int>|#)\n");
        fprintf(stderr,"             | @[<scaffold>[-<scaffold>]] | @<scaffold>_<int>-(<int>|#)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"        <contig>   = (<int>|#)[.(<int>|#)]\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"        <scaffold> =  <int>|<string>|#\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -a: Show the alignment of each LA.\n");
        fprintf(stderr,"      -r: Show the alignment of each LA with -w bp's of A in each row.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -U: Show alignments in upper case.\n");
        fprintf(stderr,"      -i: Indent alignments and cartoons by -i.\n");
        fprintf(stderr,"      -w: Width of each row of alignment in symbols (-a) or bps (-r).\n");
        fprintf(stderr,"      -b: # of border bp.s to show on each side of LA.\n");
        exit (1);
      }
  }

  //  Initiate .las file reading and read header information
  
  { char  *pwd, *root, *cpath;
    FILE  *test;

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".1aln");
    input = open_Aln_Read(Catenate(pwd,"/",root,".1aln"),1,
			   &novl,&tspace,&db1_name,&db2_name,&cpath) ;
    if (input == NULL)
      exit (1);
    free(root);
    free(pwd);

    test = fopen(db1_name,"r");
    if (test == NULL)
      { test = fopen(Catenate(cpath,db1_name,"",""),"r");
        if (test == NULL)
          { fprintf(stderr,"%s: Could not find GDB %s\n",Prog_Name,db1_name);
            exit (1);
          }
        pwd = Strdup(Catenate(cpath,db1_name,"",""),"Allocating expanded name");
	free(db1_name);
        db1_name = pwd;
      }
    fclose(test);

    if (db2_name != NULL)
      { test = fopen(db2_name,"r");
        if (test == NULL)
          { test = fopen(Catenate(cpath,db2_name,"",""),"r");
            if (test == NULL)
              { fprintf(stderr,"%s: Could not find GDB %s\n",Prog_Name,db2_name);
                exit (1);
              }
            pwd = Strdup(Catenate(cpath,db2_name,"",""),"Allocating expanded name");
            free(db2_name);
            db2_name = pwd;
          }
        fclose(test);
      }

    free(cpath);
  }

  //  Open DB or DB pair and set up scaffold->contig maps

  { int r, s, hdrs;
    struct stat state;

    ISTWO  = 0;
    if (Open_DB(db1_name,db1) < 0)
      exit (1);

    if (db2_name != NULL)
      { if (Open_DB(db2_name,db2) < 0)
          exit (1);
        ISTWO = 1;
      }
    else
      db2 = db1;

    nacontig = db1->nreads;
    nbcontig = db2->nreads;
    AREADS   = db1->reads;
    BREADS   = db2->reads;

    ascf = (int *) Malloc(sizeof(int)*(4*nacontig+1),"Allocating scaffold map");
    if (ascf == NULL)
      exit (1);
    alen = ascf + nacontig;
    actg = alen + nacontig;
    amap = actg + nacontig;

    s = -1;
    amaxlen = 0;
    for (r = 0; r < nacontig; r++)
      { if (db1->reads[r].origin == 0)
          { s += 1;
            amap[s] = r;
          }
        alen[s] = db1->reads[r].fpulse + db1->reads[r].rlen;
        ascf[r] = s;
        if (alen[s] > amaxlen)
          amaxlen = alen[s];
      }
    nascaff = s+1;
    amap[nascaff] = nacontig;
    actgmax = 0;
    for (r = nacontig-1; r >= 0; r--)
      { alen[r] = alen[ascf[r]];
        actg[r] = r - amap[ascf[r]];
        if (actg[r] > actgmax)
          actgmax = actg[r];
      }

    hdrs = open(Catenate(db1->path,".hdr","",""),O_RDONLY);
    if (hdrs < 0)
      { fprintf(stderr,"%s: Cannot open header file of %s\n",Prog_Name,db1_name);
        exit (1);
      }
    if (fstat(hdrs,&state) < 0)
      { fprintf(stderr,"%s: Cannot fetch size of %s's header file\n",Prog_Name,db1_name);
        exit (1);
      }

    AMAP     = amap;
    HPERMA   = Malloc(sizeof(int)*nascaff,"Allocating header table");
    AHEADERS = Malloc(state.st_size,"Allocating header table");
    if (HPERMA == NULL || AHEADERS == NULL)
      exit (1);

    if (read(hdrs,AHEADERS,state.st_size) < 0)
      { fprintf(stderr,"%s: Cannot read header file of %s\n",Prog_Name,db1_name);
        exit (1);
      }
    AHEADERS[state.st_size-1] = '\0';
    close(hdrs);

    HPERMA[0] = 0;
    for (s = 1; s < nascaff; s++)
      { AHEADERS[AREADS[AMAP[s]].coff-1] = '\0';
        HPERMA[s] = s;
      }

    qsort(HPERMA,nascaff,sizeof(int),HSORTA);

    if (db2_name != NULL)
      { bscf = (int *) Malloc(sizeof(int)*(4*nbcontig+1),"Allocating scaffold map");
        if (bscf == NULL)
          exit (1);
        blen = bscf + nbcontig;
        bctg = blen + nbcontig;
        bmap = bctg + nbcontig;

        s = -1;
        bmaxlen = 0;
        for (r = 0; r < nbcontig; r++)
          { if (db2->reads[r].origin == 0)
              { s += 1;
                bmap[s] = r;
              }
            blen[s] = db2->reads[r].fpulse + db2->reads[r].rlen;
            bscf[r] = s;
            if (blen[s] > bmaxlen)
              bmaxlen = blen[s];
          }
        nbscaff = s+1;
        bmap[nbscaff] = nbcontig;
        bctgmax = 0;
        for (r = nbcontig-1; r >= 0; r--)
          { blen[r] = blen[bscf[r]];
            bctg[r] = r - bmap[bscf[r]];
            if (bctg[r] > bctgmax)
              bctgmax = bctg[r];
          }

        hdrs = open(Catenate(db2->path,".hdr","",""),O_RDONLY);
        if (hdrs < 0)
          { fprintf(stderr,"%s: Cannot open header file of %s\n",Prog_Name,db2_name);
            exit (1);
          }
        if (fstat(hdrs,&state) < 0)
          { fprintf(stderr,"%s: Cannot fetch size of %s's header file\n",Prog_Name,db2_name);
            exit (1);
          }

        BMAP     = bmap;
        HPERMB   = Malloc(sizeof(int)*nbscaff,"Allocating header table");
        BHEADERS = Malloc(state.st_size,"Allocating header table");
        if (HPERMB == NULL || BHEADERS == NULL)
          exit (1);

        if (read(hdrs,BHEADERS,state.st_size) < 0)
          { fprintf(stderr,"%s: Cannot read header file of %s\n",Prog_Name,db1_name);
            exit (1);
          }
        BHEADERS[state.st_size-1] = '\0';
        close(hdrs);

        HPERMB[0] = 0;
        for (s = 1; s < nbscaff; s++)
          { BHEADERS[BREADS[BMAP[s]].coff-1] = '\0';
            HPERMB[s] = s;
          }

        qsort(HPERMB,nbscaff,sizeof(int),HSORTB);
      }
    else
      { bmap = amap;
        blen = alen;
        bctg = actg;
        bscf = ascf;
        nbscaff = nascaff;
        bctgmax = actgmax;
        bmaxlen = amaxlen;

        BMAP     = bmap;
        HPERMB   = HPERMA;
        BHEADERS = AHEADERS;
      }
  }

  //  Setup up contig reporting ranges

  if (argc > 2)
    { contig_range(argv[2],db1,nacontig,nascaff,amap,AHEAD,HPERMA,&abeg,&aend,&afst,&alst);
      if (argc > 3)
        contig_range(argv[3],db2,nbcontig,nbscaff,bmap,BHEAD,HPERMB,&bbeg,&bend,&bfst,&blst);
      else
        { bbeg = 0;
          bend = db2->nreads-1;
          bfst = 0;
          blst = db2->reads[bend].rlen;
        }
    }
  else
    { abeg = bbeg = 0;
      aend = db1->nreads-1;
      bend = db2->nreads-1;
      afst = bfst = 0;
      alst = db1->reads[aend].rlen;
      blst = db2->reads[bend].rlen;
    }
  if (afst < 0)
    { abst = AHEAD(abeg);
      aest = AHEAD(aend);
    }
  if (bfst < 0)
    { bbst = BHEAD(bbeg);
      best = BHEAD(bend);
    }

  //  Read the file and display selected records
  
  { int        j;
    uint16    *trace;
    Work_Data *work;
    int        tmax;
    int64      tps;

    int        aread, bread;
    int        aoffs, boffs;
    int        alens, blens;
    char      *abuffer, *bbuffer, *root;
    int        ar_wide, br_wide;
    int        ai_wide, bi_wide;
    int        ac_wide, bc_wide;
    int        mn_wide, mx_wide;
    int        tp_wide;

    aln->path = &(ovl->path);
    if (ALIGN || REFERENCE)
      { work = New_Work_Data();
        abuffer = New_Read_Buffer(db1);
        bbuffer = New_Read_Buffer(db2);
        if (abuffer == NULL || bbuffer == NULL)
          exit (1);
      }
    else
      { abuffer = NULL;
        bbuffer = NULL;
        work = NULL;
      }

    tmax = input->info['T']->given.max;
    trace = (uint16 *) Malloc(2*sizeof(uint16)*tmax,"Allocating trace vector");
    if (trace == NULL)
      exit (1);
    ovl->path.trace = (void *) trace;

    ar_wide = Number_Digits((int64) nascaff);
    ai_wide = Number_Digits((int64) amaxlen);
    ac_wide = Number_Digits((int64) actgmax+1);

    br_wide = Number_Digits((int64) nbscaff);
    bi_wide = Number_Digits((int64) bmaxlen);
    bc_wide = Number_Digits((int64) bctgmax+1);

    if (db1->maxlen < db2->maxlen)
      { mn_wide = Number_Digits((int64) db1->maxlen);
        mx_wide = bi_wide;
        if (tspace > 0)
          tp_wide = Number_Digits((int64) db1->maxlen/tspace+2);
        else
          tp_wide = 0;
      }
    else
      { mn_wide = Number_Digits((int64) db2->maxlen);
        mx_wide = ai_wide;
        if (tspace > 0)
          tp_wide = Number_Digits((int64) db2->maxlen/tspace+2);
        else
          tp_wide = 0;
      }
    ar_wide += (ar_wide-1)/3;
    br_wide += (br_wide-1)/3;
    ai_wide += (ai_wide-1)/3;
    bi_wide += (bi_wide-1)/3;
    mn_wide += (mn_wide-1)/3;
    tp_wide += (tp_wide-1)/3;

    root  = Root(argv[1],".las");
    printf("\n%s: ",root);
    Print_Number(novl,0,stdout);
    printf(" records\n");
    free(root);

    //  For each record do

    for (j = 0; j < novl; j++)

       //  Read it in

      { Read_Aln_Overlap(input,ovl);
        ovl->path.tlen = Read_Aln_Trace(input,(uint8 *) trace);

        aread = ovl->aread;
        bread = ovl->bread;

        if (aread >= db1->nreads)
          { fprintf(stderr,"%s: A-read is out-of-range of DB %s\n",Prog_Name,argv[1]);
            exit (1);
          }
        if (bread >= db2->nreads)
          { fprintf(stderr,"%s: B-read is out-of-range of DB %s\n",Prog_Name,argv[1+ISTWO]);
            exit (1);
          }

        //  Determine if it should be displayed

        if (afst < 0)
          { char *ahead = AHEADERS+AREADS[aread].coff;

            if (strcmp(ahead,abst) < 0)
              continue;
            if (strcmp(aest,ahead) < 0)
              continue;
          }
        else
          { if (aread < abeg || (aread == abeg && ovl->path.aepos <= afst))
              continue;
            if (aread > aend)
              break;
            if (aread == aend && ovl->path.abpos >= alst)
              continue;
          }
        if (bfst < 0)
          { char *bhead = BHEADERS+BREADS[bread].coff;

            if (strcmp(bhead,bbst) < 0)
              continue;
            if (strcmp(best,bhead) < 0)
              continue;
          }
        else
          { if (bread < bbeg || (bread == bbeg && ovl->path.bepos <= bfst))
              continue;
            if (bread > bend || (bread == bend && ovl->path.bbpos >= blst))
              continue;
          }

        aln->alen  = db1->reads[aread].rlen;
        aln->blen  = db2->reads[bread].rlen;
        aoffs = db1->reads[aread].fpulse;
        alens = alen[aread];
        boffs = db2->reads[bread].fpulse;
        blens = blen[bread];
        aln->flags = ovl->flags;
        tps        = ovl->path.tlen/2;

        //  Display it
            
        if (ALIGN || REFERENCE)
          printf("\n");

        Print_Number((int64) ascf[aread]+1,ar_wide+1,stdout);
        printf(".%0*d",ac_wide,actg[aread]+1);
        printf("  ");
        Print_Number((int64) bscf[bread]+1,br_wide+1,stdout);
        printf(".%0*d",bc_wide,bctg[bread]+1);
        if (COMP(ovl->flags))
          printf(" c");
        else
          printf(" n");
        if (ovl->path.abpos+aoffs == 0)
          printf("   <");
        else
          printf("   [");
        Print_Number((int64) ovl->path.abpos+aoffs,ai_wide,stdout);
        printf("..");
        Print_Number((int64) ovl->path.aepos+aoffs,ai_wide,stdout);
        if (ovl->path.aepos+aoffs == alens)
          printf("> x ");
        else
          printf("] x ");
        if (ovl->path.bbpos+boffs == 0)
          printf("<");
        else
          printf("[");
        if (COMP(ovl->flags))
          { Print_Number((int64) (blens - (ovl->path.bbpos+boffs)),bi_wide,stdout);
            printf("..");
            Print_Number((int64) (blens - (ovl->path.bepos+boffs)),bi_wide,stdout);
          }
        else
          { Print_Number((int64) ovl->path.bbpos+boffs,bi_wide,stdout);
            printf("..");
            Print_Number((int64) ovl->path.bepos+boffs,bi_wide,stdout);
          }
        if (ovl->path.bepos+boffs == blens)
          printf(">");
        else
          printf("]");

        printf("  ~  %5.2f%% ",(200.*ovl->path.diffs) /
               ((ovl->path.aepos - ovl->path.abpos) + (ovl->path.bepos - ovl->path.bbpos)) );
        printf("  (");
        Print_Number(alens,ai_wide,stdout);
        printf(" x ");
        Print_Number(blens,bi_wide,stdout);
        printf(" bps,");
        Print_Number((int64) ovl->path.diffs,mn_wide,stdout);
        printf(" diffs, ");
        Print_Number(tps,tp_wide,stdout);
        printf(" trace pts)\n");

        if (ALIGN || REFERENCE)
          { char *aseq, *bseq;
            int   amin,  amax;
            int   bmin,  bmax;
            int   self;

            Decompress_TraceTo16(ovl);

            self = (ISTWO == 0) && (aread == bread) && !COMP(ovl->flags);

            amin = ovl->path.abpos - BORDER;
            if (amin < 0) amin = 0;
            amax = ovl->path.aepos + BORDER;
            if (amax > aln->alen) amax = aln->alen;
            if (COMP(aln->flags))
              { bmin = (aln->blen-ovl->path.bepos) - BORDER;
                if (bmin < 0) bmin = 0;
                bmax = (aln->blen-ovl->path.bbpos) + BORDER;
                if (bmax > aln->blen) bmax = aln->blen;
              }
            else
              { bmin = ovl->path.bbpos - BORDER;
                if (bmin < 0) bmin = 0;
                bmax = ovl->path.bepos + BORDER;
                if (bmax > aln->blen) bmax = aln->blen;
                if (self)
                  { if (bmin < amin)
                      amin = bmin;
                    if (bmax > amax)
                      amax = bmax;
                  }
              }

            aseq = Load_Subread(db1,aread,amin,amax,abuffer,0);
            if (!self)
              bseq = Load_Subread(db2,bread,bmin,bmax,bbuffer,0);
            else
              bseq = aseq;

            aln->aseq = aseq - amin;
            if (COMP(aln->flags))
              { Complement_Seq(bseq,bmax-bmin);
                aln->bseq = bseq - (aln->blen - bmax);
              }
            else if (self)
              aln->bseq = aln->aseq;
            else
              aln->bseq = bseq - bmin;

            Compute_Trace_PTS(aln,work,tspace,GREEDIEST);

            if (Gap_Improver(aln,work))
              exit (1);

            { int *trace = aln->path->trace;
              int  tlen  = aln->path->tlen;
              int  i;

              aln->path->abpos += aoffs;
              aln->path->aepos += aoffs;
              aln->alen = alens;
              aln->path->bbpos += boffs;
              aln->path->bepos += boffs;
              aln->blen = blens;

              aln->aseq -= aoffs;
              aln->bseq -= boffs;
              for (i = 0; i < tlen; i++)
                if (trace[i] < 0)
                  trace[i] -= aoffs;
                else
                  trace[i] += boffs;
            }

            if (REFERENCE)
              Print_Reference(stdout,aln,work,INDENT,WIDTH,BORDER,UPPERCASE,mx_wide);
            if (ALIGN)
              Print_Alignment(stdout,aln,work,INDENT,WIDTH,BORDER,UPPERCASE,mx_wide);
          }
      }

    free(trace);
    if (ALIGN)
      { free(bbuffer-1);
        free(abuffer-1);
        Free_Work_Data(work);
      }
  }

  if (ISTWO)
    free(bscf);
  free(ascf);
  Close_DB(db1);
  if (ISTWO)
    Close_DB(db2);

  oneFileClose(input);

  exit (0);
}
