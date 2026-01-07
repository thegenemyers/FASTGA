/*******************************************************************************************
 *
 *  Module for parsing & interpreting what part of the genomes to display.
 *    
 *    Contig_Range *interpret_range(char *selection, GDB *gdb, Hash_Table *hash)
 *
 *  Output is array xSEQ, xBEG, xEND (x = A or B) that indicate if a contig
 *    is to be displayed and the range thereof.
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <zlib.h>

#include "GDB.h"
#include "hash.h"
#include "select.h"

#undef DEBUG_PARSE_SEQ

  //  Syntax symbols

#define SCAF_SYMBOL '@'
#define CONT_SYMBOL '.'
#define POST_SYMBOL ':'
#define LAST_SYMBOL '#'
#define RANG_SYMBOL '-'
#define DLIM_SYMBOL ','

//  Read next line into a buffer and return a pointer to the buffer
//    the length of the line.  NB: replaces '\n' with '\0'.

static char *RE = "RE";

static char *read_line(void *input, int gzipd, int nline, char *spath)
{ static char *buffer;
  static int   bmax = 0;
  int len;

  if (bmax == 0)
    { bmax = 500;
      buffer = (char *) malloc(bmax);
      if (buffer == NULL)
        { EPRINTF("Out of memory reading %s",spath);
          EXIT(RE);
        }
    }

  if (gzipd)
    { if (gzgets(input,buffer,bmax) == NULL)
        { if (gzeof(input))
            { free(buffer);
              bmax = 0;
              return (NULL);
            }
          EPRINTF("Could not read next line, %d, of file %s",nline,spath);
          EXIT(RE);
        }
    }
  else
    { if (fgets(buffer,bmax,input) == NULL)
        { if (feof(input))
            { free(buffer);
              bmax = 0;
              return (NULL);
            }
          EPRINTF("Could not read next line, %d, of file %s",nline,spath);
          EXIT(RE);
        }
    }

  len = strlen(buffer);
  while (buffer[len-1] != '\n')
    { bmax = ((int) (1.4*bmax)) + 100;
      buffer = (char *) realloc(buffer,bmax);
      if (buffer == NULL)
        { EPRINTF("Out of memory reading %s",spath);
          EXIT(RE);
        }
      if (gzipd)
        { if (gzgets(input,buffer+len,bmax-len) == NULL)
            { if (gzeof(input))
                EPRINTF("Last line %d of file %s does not end with new-line",
                        nline,spath);
              else
                EPRINTF("Could not read next line, %d, of file %s",
                        nline,spath);
              EXIT(RE);
            }
        }
      else
        { if (fgets(buffer+len,bmax-len,input) == NULL)
            { if (feof(input))
                EPRINTF("Last line %d of file %s does not end with new-line",
                        nline,spath);
              else
                EPRINTF("Could not read next line, %d, of file %s",
                        nline,spath);
              EXIT(RE);
            }
        }
      len += strlen(buffer+len);
    }
  buffer[--len] = '\0';

  return (buffer);
}

//  Parse a range expression.
//    val[0] = val[4] = scaffold
//    val[1] = val[5] = contig
//    val[2] = val[6] = position
//    val[3] = val[7] = position precision
//    val[8] = sign
//  val[0-3] is the 1st location of a range, val[4-7] the 2nd location
//  If any item is -1 then denotes last (#), and if -2 then not present

static char *src;

static int follow[128] =      //  isspace !isprint + - . :
  { 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 0, 0, 0, 0, 0, 0, 0,

    0, 0, 0, 1, 0, 1, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0,

    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,

    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1
  };

static inline char *white(char *x)
{ while (isspace(*x))
    x += 1;
  return (x);
}

static inline char *get_int(char *x, int64 *v, int *n)
{ *v = 0;
  *n = 1;
  while (isdigit(*x))
    { *v  = 10*(*v) + (*x++ - '0');
      *n *= 10;
    }
  return (x);
}

static char *get_bps(char *x, int64 *v)
{ int   a, n;
  int64 r, p, m;

  x = get_int(x,&r,&n);
  x = white(x);
  if (*x == '.')
    { x = white(x+1);
      if (isdigit(*x))
        x = get_int(x,&p,&n);
      else
        { EPRINTF("Location . not followed by integer\n\t%.*s ^ %s",(int) (x-src),src,x);
          EXIT(NULL);
        }
      x = white(x);
    }
  else
    { p = 0;
      n = 1;
    }
  a = *x++;
  if (a == 'G')
    m = 1000000000;
  else if (a == 'M')
    m = 1000000;
  else if (a == 'k')
    m = 1000;
  else
    { m = 1;
      x -= 1;
    }
  if (p >= m)
    { EPRINTF("Location precision has more digits than multiplier %c\n\t%.*s ^ %s",
              (int) a,(int) (x-src),src,x);
      EXIT(NULL);
    }
  m /= n;
  *v = (r*n + p) * m;
  v[1] = m;
  return (x);
}

static char *get_location(char *x, int64 *v, Hash_Table *hash)
{ int   a, i;
  char *y;

  v[0] = v[1] = v[2] = v[3] = -2;
  x = white(x);

  if (*x == SCAF_SYMBOL)
    { x = white(x+1);
      if (*x == LAST_SYMBOL)
        { v[0] = -1;
          x += 1;
        }
      else if (isdigit(*x))
        { x = get_int(x,v,&a);
          if (v[0] == 0)
            { EPRINTF("Scaffold index cannot be 0\n\t%.*s ^ %s",(int) (x-src),src,x);
              EXIT(NULL);
            }
        }
      else
        { y = x;
          while ( ! follow[(int) (*x)] )
            x += 1;
          a = *x;
          *x = '\0';
          i = Hash_Lookup(hash,y);
          *x = a;
          if (i < 0)
            { EPRINTF("Could not parse scaffold item\n\t%.*s ^ %s",(int) (y-src),src,y);
              EXIT(NULL);
            }
          v[0] = i+1;
        }
      x = white(x);
    }
  if (*x == CONT_SYMBOL)
    { x = white(x+1);
      if (*x == LAST_SYMBOL)
        { v[1] = -1;
          x += 1;
        }
      else if (isdigit(*x))
        { x = get_int(x,v+1,&a);
          if (v[1] == 0)
            { EPRINTF("Contig index cannot be 0\n\t%.*s ^ %s",(int) (x-src),src,x);
              EXIT(NULL);
            }
        }
      else
        { EPRINTF("Contig is not an integer or #-sign\n\t%*s ^ %.s",(int) (x-src),src,x);
          EXIT(NULL);
        }
      x = white(x);
    }
  if (v[0] >= -1 || v[1] >= -1)
    { if (*x == POST_SYMBOL)
        { x = white(x+1);
          if (*x == LAST_SYMBOL)
            { v[2] = -1;
              x += 1;
            }
          else if (isdigit(*x))
            x = get_bps(x,v+2);
          else
            { EPRINTF("Position is not an integer or #-sign\n\t%.*s ^ %s",
                              (int) (x-src),src,x);
              EXIT(NULL);
            }
        }
    }
  else if (*x == LAST_SYMBOL)
    { v[2] = -1;
      x += 1;
    }
  else if (isdigit(*x))
    x = get_bps(x,v+2);
  else
    { EPRINTF("Empty location\n\t%.*s ^ %s",(int) (x-src),src,x);
      EXIT(NULL);
    }

  return (x);
}

static int get_focus(char *x, int64 *v, Hash_Table *hash1, Hash_Table *hash2)
{ src = x;
 
  v[8] = 0;
  x = get_location(x,v,hash1);
  if (x == NULL)
    return (1);
  if (v[2] < -1)
    { EPRINTF("Focus location must give a position\n\t%.*s ^ %s",(int) (x-src),src,x);
      EXIT(1);
    }
  x = white(x);
  if (*x != DLIM_SYMBOL)
    { EPRINTF("Focus locations must be separatd by a '%c'\n\t%.*s ^ %s",
                     DLIM_SYMBOL,(int) (x-src),src,x);
      EXIT(1);
    }
  x = get_location(white(x+1),v+4,hash2);
  if (x == NULL)
    return (1);
  if (v[6] < -1)
    { EPRINTF("Focus location must give a position\n\t%.*s ^ %s",(int) (x-src),src,x);
      EXIT(1);
    }
  x = white(x);
  if (*x != '\0')
    { EPRINTF("Focus syntax is not complete\n\t%.*s ^ %s",(int) (x-src),src,x);
      EXIT(1);
    }

  return (0);
}

static int get_range(char *x, int64 *v, Hash_Table *hash)
{ char *y;

  src = x;

  y = x + strlen(x) - 1;        //  Clip tailing +/- from string
  while (isspace(*y))
    y -= 1;
  if (*y == '+')
    { v[8] = +1;
      *y = '\0';
    }
  else if (*y == '-')
    { v[8] = -1;
      *y = '\0';
    }
  else
    v[8] = 0;

  x = get_location(x,v,hash);
  if (x == NULL)
    return (1);

  x = white(x);
  if (*x == RANG_SYMBOL)
    { x = get_location(x+1,v+4,hash);
      if (x == NULL)
        return (1);
      x = white(x);
    }
  else
    v[4] = v[5] = v[6] = -2;

  if (*x != '\0')
    { EPRINTF("Range syntax is not complete\n\t%.*s ^ %s",(int) (x-src),src,x);
      EXIT(1);
    }

  return (0);
}

  //  Determine any missing fields so that v[0] is the scaffold, v[1] is the absolute
  //    contig index, and v[2] is the contig relative position.

static int complete_address(int64 *v, GDB *gdb, int first)
{ int64 s, c, p, q;
  int   fc, ec;
  int64 cl;

  int nscaff  = gdb->nscaff;
  int ncontig = gdb->ncontig;

  GDB_CONTIG *contig = gdb->contigs;
  GDB_SCAFFOLD *scaff = gdb->scaffolds;

  s = v[0];
  c = v[1];
  p = q = v[2];

  if (s < -1)
    { if (c < -1)
        { if (p == -1)
            { s = nscaff-1;
              c = ncontig-1;
              p = contig[c].clen;
            }
          else
            { for (s = 0; s < nscaff; s++)
                if (p > scaff[s].slen)
                  p -= scaff[s].slen;
                else
                  break;
              if (s >= nscaff && p > v[3])
                { EPRINTF("Position %lld is larger than genome",q);
                  EXIT(1);
                }
              fc = scaff[s].fctg;
              ec = scaff[s].ectg;
              for (c = fc; c < ec; c++)
                if (p > contig[c].clen)
                  p -= contig[c].clen;
                else
                  break;
            }
        }
      else
        { if (c == -1)
            { s = nscaff-1;
              c = ncontig-1;
            }
          else
            { if (c > ncontig)
                { EPRINTF("Contig %lld is > %d, the # of contigs",c,ncontig);
                  EXIT(1);
                }
              c = c-1;
              for (s = 0; s < nscaff; s++)
                if (c < scaff[s].ectg)
                  break;
            }
          cl = contig[c].clen;
          if (p < -1)
            { if (first)
                p = 0;
              else
                p = cl;
            }
          else if (p == -1)
            p = cl;
          else
            { if (p > cl+v[3])
                { EPRINTF("Position %lld beyond contig %lld of length %lld",p,c,cl);
                  EXIT(1);
                }
            }
        }
    }
  else
    { if (s == -1)
        s = gdb->nscaff-1;
      else
        { if (s >= nscaff)
            { EPRINTF("Scaffold %lld does not exist, only %d scaffolds",s,nscaff);
              EXIT(1);
            }
          s = s-1;
        }
      fc = scaff[s].fctg;
      ec = scaff[s].ectg;
      if (c < -1)
        { if (p < -1)
            { if (first)
                { c = fc;
                  p = 0;
                }
              else
                { c = ec-1;
                  p = contig[c].clen;
                }
            }
          else if (p == -1)
            { c = ec-1;
              p = contig[c].clen;
            }
          else
            { for (c = fc; c < ec; c++)
                if (p < contig[c].sbeg)
                  break;
              c -= 1;
              p -= contig[c].sbeg;
              if (c == ec-1 && p > contig[c].clen + v[3])
                { EPRINTF("Position %lld is beyond scaffold %lld of length %lld",
                                 q,s,scaff[s].slen);
                  EXIT(1);
                }
            }
        }
      else
        { if (c == -1)
            c = ec-1;
          else
            { if (c > ec-fc)
                { EPRINTF("Contig %lld is > %d, the # of contigs in scaffold %lld",
                                 c,ec-fc,s);
                  EXIT(1);
                }
              c += fc - 1;
            }
          cl = contig[c].clen;
          if (p < -1)
            { if (first)
                p = 0;
              else
                p = cl;
            }
          else if (p == -1)
            p = cl;
          else
            { if (p > cl+v[3])
                { EPRINTF("Position %lld beyond contig %lld of length %lld",p,c,cl);
                  EXIT(1);
                }
            }
        }
    }

  v[0] = s;
  v[1] = c;
  v[2] = p;
  return (0);
}

int interpret_point(Selection *s, char *x, GDB *gdb1, Hash_Table *hash1,
                                           GDB *gdb2, Hash_Table *hash2)
{ int64 v[9];


  if (get_focus(x,v,hash1,hash2))
    return (1);

#ifdef DEBUG_PARSE_SEQ
  fprintf(stderr,"RAW: ");
  for (int i = 0; i < 9; i++)
    fprintf(stderr," %lld",v[i]);
  fprintf(stderr,"\n");
#endif

  if (complete_address(v,gdb1,1))
    return (1);
  if (complete_address(v+4,gdb2,1))
    return (1);

  s->type   = POINT_SELECTION;
  s->orient = 0;
  s->s1 = v[0];
  s->c1 = v[1];
  s->p1 = v[2];
  s->s2 = v[4];
  s->c2 = v[5];
  s->p2 = v[6];

#ifdef DEBUG_PARSE_SEQ
  printf("SEL: %d %d,%d,%lld  %d,%d,%lld (%d)\n",
         s->type,s->s1,s->c1,s->p1,s->s2,s->c2,s->p2,s->orient);
#endif

  return (0);
}


  //  Convert the raw selection v[0..6] into a selection record of 2 locations

int interpret_range(Selection *s, char *x, GDB *gdb, Hash_Table *hash)
{ int64 v[9];
  char *y;
  int   a, special;

  special = 10;
  y = white(x);
  a = *y;
  if (a == '@' || a == '.')
    { y = white(y+1);
      if (*y == '\0')
        special = 0;
      else
        { if (*y == '-')
            special = -1;
          else if (*y == '+')
            special = 1;
          y = white(y+1);
          if (*y != '\0')
            special = 10;
        }
    }

  if (special < 10)
    { if (a == '@')
        s->type = SCAFF_SELECTION;
      else
        s->type = CONTG_SELECTION;
      s->s1 = 0;
      s->c1 = 0;
      s->p1 = 0;
      s->s2 = gdb->nscaff-1;
      s->c2 = gdb->ncontig-1;
      s->p2 = gdb->contigs[s->c2].clen;
      s->orient = special;
    }

  else
    { if (get_range(x,v,hash))
        return (1);

#ifdef DEBUG_PARSE_SEQ
      printf("'%s'\n",x);
      fprintf(stderr,"RAW: ");
      for (int i = 0; i < 9; i++)
        fprintf(stderr," %lld",v[i]);
      fprintf(stderr,"\n");
#endif

      if (v[0] < -1)
        s->type = CONTG_SELECTION;
      else
        s->type = SCAFF_SELECTION;
      s->orient = v[8];

      if (v[4] < -1 && v[5] < -1 && v[6] < -1)
        { if (v[2] >= -1)
            { EPRINTF("Must specify a range, not a point");
              EXIT(1);
            }
          v[4] = v[0];
          v[5] = v[1];
        }

      else
        { if (v[4] < -1)
            { v[4] = v[0];
              if (v[5] < -1)
                v[5] = v[1];
            }
        }

      if (complete_address(v,gdb,1))
        return (1);
      if (complete_address(v+4,gdb,0))
        return (1);

      s->s1 = v[0];
      s->c1 = v[1];
      s->p1 = v[2];
      s->s2 = v[4];
      s->c2 = v[5];
      s->p2 = v[6];
    }
#ifdef DEBUG_PARSE_SEQ
  printf("SEL: %d %d,%d,%lld  %d,%d,%lld (%d)\n",
         s->type,s->s1,s->c1,s->p1,s->s2,s->c2,s->p2,s->orient);
#endif
  return (0);
}

static int get_selection_contigs_from_file(FILE *fp, GDB *gdb, Hash_Table *hash,
                                           Contig_Range *chord, char *filename, int ordered)
{ char       c, *p, *q, *line;
  int        i, pbeg, pend, pfst, plst, ori;
  int        nline, order;
  Selection _s, *s = &_s;

  order = 1;
  nline = 1;
  while ((line = read_line(fp,0,nline++,filename)) != NULL)
    { if (line == RE)
        EXIT(1);
      p = q = white(line);
      while (*q != '\0' && !isspace(*q))
        q += 1;
      if (p == q)      // empty line
        continue;

      c = *q;
      *q = '\0';
      if (interpret_range(s,p,gdb,hash))
        EXIT(1);
      *q = c;
      pbeg = s->c1;
      pend = s->c2;
      pfst = s->p1;
      plst = s->p2;
      ori  = s->orient;
#ifdef DEBUG_PARSE_SEQ
      fprintf(stderr,"%s: CTG_RANGE pbeg = %-12d pend = %-12d pfst = %-12d plst = %-12d ori = %d\n",
                     Prog_Name,pbeg,pend,pfst,plst,ori);
#endif
      if (ordered)
        { for (i = pbeg; i <= pend; i++)
            if (chord[i].order)
              { EPRINTF("Overlapping contigs in selection ranges");
                EXIT(1);
              }
        }
      else if (ori != 0)
        { for (i = pbeg; i <= pend; i++)
            if (chord[i].order && ori*chord[i].orient < 0)
              { EPRINTF("Conflicting sign for contig in selection expression");
                EXIT(1);
              }
        }

      for (i = pbeg+1; i < pend; i++)
        { chord[i].order  = order;
          chord[i].beg    = 0;
          chord[i].end    = gdb->contigs[i].clen;
          chord[i].orient = ori;
        }
      if (pbeg != pend)
        { if (chord[pend].order)
            { if (chord[pend].end < plst)
                chord[pend].end = plst;
            }
          else
            { chord[pend].order = order;
              chord[pend].end   = plst;
            }
          chord[pend].beg = 0;
          if (chord[pbeg].order)
            { if (chord[pbeg].beg > pfst)
                chord[pbeg].beg = pfst;
            }
          else
            { chord[pbeg].order = order;
              chord[pbeg].beg   = pfst;
            }
          chord[pbeg].end = gdb->contigs[pbeg].clen;
          chord[pbeg].orient = ori;
          chord[pend].orient = ori;
        }
      else
        { if (chord[pend].order)
            { if (chord[pend].end < plst)
                chord[pend].end = plst;
              if (chord[pbeg].beg > pfst)
                chord[pbeg].beg = pfst;
            }
          else
            { chord[pend].order = order;
              chord[pend].end   = plst;
              chord[pbeg].beg   = pfst;
            }
          chord[pbeg].orient = ori;
        }

      order += 1;
    }

  fclose(fp);
  return (0);
}

Contig_Range *get_selection_contigs(char *x, GDB *gdb, Hash_Table *hash, int ordered)
{ char         *e;
  int           i, pbeg, pend, pfst, plst, ori;
  int           order;
  Contig_Range *chord;
  FILE         *fp;
  Selection    _s, *s = &_s;

  chord = (Contig_Range *) Malloc(sizeof(Contig_Range)*gdb->ncontig,"Allocating sequence array");

  if (x == NULL)
    { for (i = 0; i < gdb->ncontig; i++)
        { chord[i].order  = 1;
          chord[i].beg    = 0;
          chord[i].end    = gdb->contigs[i].clen;
          chord[i].orient = 0;
        }
    }
  else
    { x = white(x);
      if (*x == '\0')
        { EPRINTF("Empty range");
          EXIT(NULL);
        }

      for (i = 0; i < gdb->ncontig; i++)
        chord[i].order = 0;

      fp = fopen(x,"r");
      if (fp != NULL)
        get_selection_contigs_from_file(fp,gdb,hash,chord,x,ordered);
      else
        { order = 1;
          while (x != NULL)
            { e = index(x,DLIM_SYMBOL);
              if (e != NULL)
                *e = '\0';
              
              if (interpret_range(s,x,gdb,hash))
                EXIT(NULL);
              pbeg = s->c1;
              pend = s->c2;
              pfst = s->p1;
              plst = s->p2;
              ori  = s->orient;
#ifdef DEBUG_PARSE_SEQ
              fprintf(stderr,
                     "%s: CTG_RANGE pbeg = %-12d pend = %-12d pfst = %-12d plst = %-12d ori = %d\n",
                     Prog_Name,pbeg,pend,pfst,plst,ori);
#endif
              if (ordered)
                { for (i = pbeg; i < pend; i++)
                    if (chord[i].order)
                      { EPRINTF("Overlapping contigs in selection ranges");
                        EXIT(NULL);
                      }
                }
              else if (ori != 0)
                { for (i = pbeg; i <= pend; i++)
                    if (chord[i].order && ori*chord[i].orient < 0)
                      { EPRINTF("Conflicting sign for contig in selection expression");
                        EXIT(NULL);
                      }
                }

              for (i = pbeg+1; i < pend; i++)
                { chord[i].order  = order;
                  chord[i].beg    = 0;
                  chord[i].end    = gdb->contigs[i].clen;
                  chord[i].orient = ori;
                }
              if (pbeg != pend)
                { if (chord[pend].order)
                    { if (chord[pend].end < plst)
                        chord[pend].end = plst;
                    }
                  else
                    { chord[pend].order = order;
                      chord[pend].end   = plst;
                    }
                  chord[pend].beg = 0;
                  if (chord[pbeg].order)
                    { if (chord[pbeg].beg > pfst)
                        chord[pbeg].beg = pfst;
                    }
                  else
                    { chord[pbeg].order = order;
                      chord[pbeg].beg   = pfst;
                    }
                  chord[pbeg].end = gdb->contigs[pbeg].clen;
                  chord[pbeg].orient = ori;
                  chord[pend].orient = ori;
                }
              else
                { if (chord[pend].order)
                    { if (chord[pend].end < plst)
                        chord[pend].end = plst;
                      if (chord[pbeg].beg > pfst)
                        chord[pbeg].beg = pfst;
                    }
                  else
                    { chord[pend].order = order;
                      chord[pend].end   = plst;
                      chord[pbeg].beg   = pfst;
                    }
                  chord[pbeg].orient = ori;
                }

              order += 1;
   
              if (e != NULL)
                *e++ = DLIM_SYMBOL;
              x = e;
            }
        }
    }

#ifdef DEBUG_PARSE_SEQ
  fprintf(stderr,"%s: Sequence range to plot\n",Prog_Name);
  fprintf(stderr,"%s: CONTIG ORDER  BEGIN        END\n",Prog_Name);
  for (i = 0; i < gdb->ncontig; i++)
    { if (chord[i].order)
        fprintf(stderr,"%s: %-8d %-6d %-12d %-12d\n",
                       Prog_Name,i,chord[i].order,chord[i].beg,chord[i].end);
    }
#endif

  return (chord);
}

static Selection *get_selection_list_from_file(FILE *fp, GDB *gdb, Hash_Table *hash,
                                               int *nlen, char *filename)
{ char       c, *p, *q, *line;
  int        len, nline;
  Selection *list;

  nline = 1;
  while ((line = read_line(fp,0,nline,filename)) != NULL)
    { if (line == RE)
        EXIT(NULL);
      nline += 1;
    }

  list = (Selection *) Malloc(sizeof(Selection)*(nline-1),"Allocating sequence array");

  rewind(fp);

  len = 0;
  while ((line = read_line(fp,0,nline,filename)) != NULL)
    { p = q = white(line);
      while (*q != '\0' && !isspace(*q))
        q += 1;
      if (p == q)      // empty line
        continue;

      c = *q;
      *q = '\0';
      if (interpret_range(list+len,p,gdb,hash))
        EXIT(NULL);
      *q = c;
#ifdef DEBUG_PARSE_SEQ
      fprintf(stderr, "%s: %c %d,%d,%lld  %d,%d,%lld (%d)\n",
              Prog_Name,list[len].type ? '@' : '.',list[len].s1,list[len].c1,list[len].p1,
              list[len].s2,list[len].c2,list[len].p2,list[len].orient);
#endif
      len += 1;
    }
  fclose(fp);

  *nlen = len;
  return (list);
}

Selection *get_selection_list(char *x, GDB *gdb, Hash_Table *hash, int *nlen)
{ char      *e, *y;
  int        len;
  Selection *list;
  FILE      *fp;

  if (x == NULL)
    { list = (Selection *) Malloc(sizeof(Selection),"Allocating sequence array");
      *nlen = 1;
      list[0].type = CONTG_SELECTION;
      list[0].s1 = 0;
      list[0].c1 = 0;
      list[0].p1 = 0;
      list[0].s2 = gdb->nscaff-1;
      list[0].c2 = gdb->ncontig-1;
      list[0].p2 = gdb->contigs[gdb->ncontig-1].clen;
      list[0].orient = 0;
      return (list);
    }

  x = white(x);
  if (*x == '\0')
    { EPRINTF("Empty range");
      EXIT(NULL);
    }

  fp = fopen(x,"r");
  if (fp != NULL)
    return (get_selection_list_from_file(fp,gdb,hash,nlen,x));

   len = 1;
   y = x-1;
   while ((y = index(y+1,DLIM_SYMBOL)) != NULL)
     len += 1;

   list = (Selection *) Malloc(sizeof(Selection)*len,"Allocating sequence array");

   len = 0;
   while (x != NULL)
     { e = index(x,DLIM_SYMBOL);
       if (e != NULL)
         *e = '\0';
       if (interpret_range(list+len,x,gdb,hash))
         EXIT(NULL);
#ifdef DEBUG_PARSE_SEQ
      fprintf(stderr, "%s: %c %d,%d,%lld  %d,%d,%lld (%d)\n",
              Prog_Name,list[len].type ? '@' : '.',list[len].s1,list[len].c1,list[len].p1,
              list[len].s2,list[len].c2,list[len].p2,list[len].orient);
#endif
       len += 1;

       if (e != NULL)
         *e++ = DLIM_SYMBOL;
       x = e;
     }

  *nlen = len;
  return (list);
}
