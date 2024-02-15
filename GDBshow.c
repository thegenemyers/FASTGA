/*******************************************************************************************
 *
 *  Display a specified set of reads of a database in fasta format.
 *
 *  Author:  Gene Myers
 *  Origin:  Reworking of DAZZ_DB program DBshow.c from the FASTGA module specialized to "gdb" ("dam")
 *  Date  :  Jan 2024
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "DB.h"

#undef DEBUG_RANGE

static char *Usage =
            "[-hU] [-w<int(80)>] <source:path[.gdb] [ <ranges:FILE> | <range> ... ]";

#define LAST_SYMBOL  '#'
#define SCAF_SYMBOL  '@'
#define MAX_BUFFER   10001

typedef struct
  { FILE  *input;
    int    lineno;
  } File_Iter;

static void init_file_iterator(FILE *input, File_Iter *it)
{ it->input  = input;
  it->lineno = 1;
  rewind(input);
}

static char *next_entry(File_Iter *it)
{ static char nbuffer[MAX_BUFFER];
  char *eol;

  do
    { if (fgets(nbuffer,MAX_BUFFER,it->input) == NULL)
        { if ( ! feof(it->input))
            { fprintf(stderr,"%s: Encountered IO error on line %d of read list file.\n",
                             Prog_Name,it->lineno);
              exit (1);
            }
          return (NULL);
        }
      if ((eol = index(nbuffer,'\n')) == NULL)
        { fprintf(stderr,"%s: Line %d in read list file is longer than %d chars!\n",
                         Prog_Name,it->lineno,MAX_BUFFER-1);
          exit (1);
        }
      *eol = '\0';
      it->lineno += 1;
      for (eol = nbuffer; isspace(*eol); eol++) ;
    }
  while (*eol == '\0');
  return (nbuffer);
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
//      indicates val[0] or val[2], respectively, is a scaffold index.

static int  ncontig;  // # of contigs (across all scaffolds)
static int  nscaff;   // # of scaffolds
static int *smap;     // smap[r] = idx of 1st contig in scaffold r in [0,nscaff]

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

static int range(char *x, int *v)
{ int t, sep;

  v[2] = -1;
  x = white(x);
  if (*x == SCAF_SYMBOL)
    { x += 1;
      sep = ' ';
    }
  else
    sep = '.';
  x = address(x,v,sep);
  if (x == NULL)
    return (0);
  t = 1;
  if (*x == '-')
    { x = address(x+1,v+2,sep);
      t = 2;
    }
  else if (*x == '_')
    { x = address(x+1,v+2,'-');
      t = 3;
    }
  if (x == NULL || *x != '\0')
    return (0);
  return (t);
}

 //  Interpret value pair s.c into an absolute contig range (1-based)

static int interpret_address(int s, int c)
{ if (c == -1)
    { if (s == 0)
        return (ncontig);
      else if (s > ncontig)
        { fprintf(stderr,"%s: Absolute contig index '%d' is out of bounds\n",Prog_Name,s);
          exit (1);
        }
      else
        return (s);
    }
  else
    { if (s == 0)
        s = nscaff;
      else if (s > nscaff)
        { fprintf(stderr,"%s: Scaffold index '%d' is out of bounds\n",Prog_Name,s);
          exit (1);
        }
      if (c == 0)
        return (smap[s]);
      else if (c == -2)
        return (s);
      else if (c > smap[s]-smap[s-1])
        { fprintf(stderr,"%s: Relative contig index '%d' is out of bounds\n",Prog_Name,c);
          exit (1);
        }
      else
        return (smap[s-1] + c);
    }
}

int main(int argc, char *argv[])
{ DAZZ_DB    _db, *db = &_db;
  FILE       *hdrs;
  FILE       *input;
  File_Iter  _iter, *iter = &_iter;
  int         rvals[4];

  int         UPPER;
  int         DOSEQ;
  int         WIDTH;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];
    char *eptr;

    ARG_INIT("GDBshow")

    WIDTH = 80;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("hU")
            break;
          case 'w':
            ARG_NON_NEGATIVE(WIDTH,"Line width")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    UPPER = 1+flags['U'];
    DOSEQ = 1-flags['h'];

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"     <range> =    <contig>[-<contig>]     |  <contig>_<int>-(<int>|#)\n");
        fprintf(stderr,"             | @[<scaffold>[-<scaffold>]] | @<scaffold>_<int>-(<int>|#)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"        <contig>   = (<int>|#)[.(<int>|#)]\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"        <scaffold> =  <int>|#\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -h: Show only the header lines.\n");
        fprintf(stderr,"      -w: Print -w bp per line (default is 80).\n");
        fprintf(stderr,"      -U: Use upper case for DNA (default is lower case).\n");
        exit (1);
      }
  }

  //  Open GDB and it .hdr file

  { int status;

    status = Open_DB(argv[1],db);
    if (status < 0)
      exit (1);
    if (status == 0)
      { fprintf(stderr,"%s: Cannot open %s as a .gdb\n",Prog_Name,argv[1]);
        exit (1);
      }
    ncontig = db->nreads;

    hdrs = Fopen(Catenate(db->path,".hdr","",""),"r");
    if (hdrs == NULL)
      exit (1);
  }

  //  Establish map from scaffold to absolute index of 1st contig

  { int r, s;

    smap = (int *) Malloc(sizeof(int)*(ncontig+1),"Allocating scaffold map");
    if (smap == NULL)
      exit (1);

    s = -1;
    for (r = 0; r < ncontig; r++)
      if (db->reads[r].origin == 0)
        { s += 1;
          smap[s] = r;
        }
    nscaff = s+1;
    smap[nscaff] = r;
  }

  //  Determine if extra arg is a range or a file, and if a file then initialize it

  input = NULL;
  if (argc == 3 && range(argv[2],rvals) <= 0)
    { input = fopen(argv[2],"r");
      if (input == NULL)
        { fprintf(stderr,"%s: '%s' is neither a valid read range nor an openable file\n",
                         Prog_Name,argv[2]);
          exit (1);
        }
      init_file_iterator(input,iter);
    }

  //  Display each read (and/or QV streams) in the active DB according to the
  //    range pairs in pts[0..reps) and according to the display options.

  { DAZZ_READ  *reads, *r;
    char       *read, *x;
    int         c, b, e, i;
    int         substr, sunits;
    int         len, fst, lst;
    char        header[MAX_NAME];
    char       *nstring;

    nstring = Malloc(WIDTH+1,"Allocating write buffer\n");
    if (nstring == NULL)
      exit (1);

    if (UPPER == 2)
      for (c = 0; c < WIDTH; c++)
        nstring[c] = 'N';
    else
      for (c = 0; c < WIDTH; c++)
        nstring[c] = 'n';
    nstring[WIDTH] = '\0';

    read   = New_Read_Buffer(db);
    if (read == NULL)
      exit (1);
    reads  = db->reads;
    substr = 0;

    c = 2;
    while (1)
      { substr = 0;
        fst = lst = 0;
        sunits = 0;
        if (argc == 2)      //  No args: show everything
          { if (c > 2)
              break;
            b = 1;
            e = ncontig;
            c = 3;
          }
        else
          { if (input != NULL)        //  Get next entry from either command line or read input file
              { x = next_entry(iter);
                if (x == NULL)
                  break;
              }
            else
              { if (c >= argc)
                  break;
                x = argv[c];
                c += 1;
              }

            i = range(x,rvals);       // Parse the entry

#ifdef DEBUG_RANGE
            printf(" %d: '%s' %d %d %d %d\n",i,x,rvals[0],rvals[1],rvals[2],rvals[3]);
#endif

            if (i == 0)
              { if (input == NULL)
                  fprintf(stderr,"%s: Command line argument %s not a valid read range\n",
                                 Prog_Name,argv[c]);
                else
                  fprintf(stderr,"%s: Line %d of read list file %s is not a valid read range\n",
                                 Prog_Name,iter->lineno,argv[2]);
                exit (1);
              }

            sunits = (rvals[1] == -2);

            //  Convert entry into an absolute contig range or an absolute contig and substring

            b = e = interpret_address(rvals[0],rvals[1]);
            if (i == 2)
              { e = interpret_address(rvals[2],rvals[3]);
                if (b > e)
                  { fprintf(stderr,"%s: Read range in '%s' is empty\n",Prog_Name,x);
                    exit (1);
                  }
              }
            else if (i == 3)
              { e = b;        //  Type 3: intepret rval[2] and rval[3] as the substring interval
                substr = 1;
                fst = rvals[2];
                lst = rvals[3];
                if (sunits)
                  len = reads[smap[b]-1].fpulse + reads[smap[b]-1].rlen;
                else
                  len = reads[b-1].rlen;
                if (lst == 0) 
                  lst = len;
                else if (lst > len)
                  { fprintf(stderr,"%s: Substring interval in '%s' is out of bounds (max = %d)\n",
                                   Prog_Name,x,len);
                    exit (1);
                  }
                if (fst >= lst)
                  { fprintf(stderr,"%s: Substring interval in '%s' is empty\n",Prog_Name,x);
                    exit (1);
                  }
              }
          }

        if (sunits)

          { for (i = b-1; i < e; i++)
              { int cbeg, cend, wpos;
                int j, u, w, f, l;

                r   = reads + smap[i];
                len = reads[smap[i+1]-1].fpulse + reads[smap[i+1]-1].rlen;
                if (!substr)
                  { fst = 0;
                    lst = len;
                  }

                fseeko(hdrs,r->coff,SEEK_SET);
                fgets(header,MAX_NAME,hdrs);
                header[strlen(header)-1] = '\0';
                printf("> %s",header);
                if (substr)
                  printf(" :: [%d,%d]",fst,lst);
                printf("\n");

                if (DOSEQ)
                  { cend = 0;
                    wpos = 0;
                    for (u = smap[i]; u < smap[i+1]; u++, r++)
                      { cbeg = r->fpulse;
                        if (cend < fst)
                          cend = fst;
                        if (cend < lst && cbeg > fst)
                          { if (cbeg <= lst)
                              len = cbeg-cend;
                            else
                              len = lst-cend;

                            for (j = 0; j+(w = WIDTH-wpos) <= len; j += w)
                              { printf("%.*s\n",w,nstring);
                                wpos = 0;
                              }
                            if (j < len)
                              { printf("%.*s",len-j,nstring);
                                wpos += len-j;
                              }
                          }

                        cend = cbeg + r->rlen;
                        if (cbeg < lst && cend > fst)
                          { if (Load_Read(db,u,read,UPPER))
                              exit (1);
                            f = fst-cbeg;
                            l = lst-cbeg;
                            if (f < 0) f = 0; 
                            if (l > r->rlen) l = r->rlen;

                            for (j = f; j+(w = WIDTH-wpos) <= l; j += w)
                              { printf("%.*s\n",w,read+j);
                                wpos = 0;
                              }
                            if (j < l)
                              { printf("%.*s",l-j,read+j);
                                wpos += l-j;
                              }
                          }
                      }
                    printf("\n");
                  }
              }
          }

        else

          { for (i = b-1; i < e; i++)
              { r   = reads + i;
                len = r->rlen;
                if (!substr)
                  { fst = 0;
                    lst = len;
                  }

                if (len > 0)
                  { fseeko(hdrs,r->coff,SEEK_SET);
                    fgets(header,MAX_NAME,hdrs);
                    header[strlen(header)-1] = '\0';
                    printf("> %s :: Contig %d[%d,%d]\n",
                           header,r->origin+1,r->fpulse+fst,r->fpulse+lst);
                  }

                if (DOSEQ)
                  { int j;
    
                    if (Load_Read(db,i,read,UPPER))
                      exit (1);
                    for (j = fst; j+WIDTH < lst; j += WIDTH)
                      printf("%.*s\n",WIDTH,read+j);
                    if (j < lst)
                      printf("%.*s\n",lst-j,read+j);
                  }
                fflush(stdout);
              }
         }
      }
  }

  if (input != NULL)
    fclose(input);

  fclose(hdrs);
  Close_DB(db);

  exit (0);
}
