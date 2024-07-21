/*******************************************************************************************
 *
 *  Module for parsing & interpreting what part of the genomes to display.
 *    
 *    Contig_Range *interpret_selection(char *selection, GDB *gdb, Hash_Table *hash)
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
#define LAST_SYMBOL '#'
#define FILE_SYMBOL ':'
#define DLIM_SYMBOL ','

//  Read next line into a buffer and return a pointer to the buffer
//    the length of the line.  NB: replaces '\n' with '\0'.

static char *read_line(void *input, int gzipd, int nline, char *spath)
{ static char *buffer;
  static int   bmax = 0;
  int len;

  if (bmax == 0)
    { bmax = 500;
      buffer = (char *) malloc(bmax);
      if (buffer == NULL)
        { fprintf(stderr,"%s: Out of memory reading %s\n",Prog_Name,spath);
          exit (1);
        }
    }

  if (gzipd)
    { if (gzgets(input,buffer,bmax) == NULL)
        { if (gzeof(input))
            { free(buffer);
              bmax = 0;
              return (NULL);
            }
          fprintf(stderr,"%s: Could not read next line, %d, of file %s\n",Prog_Name,nline,spath);
          exit (1);
        }
    }
  else
    { if (fgets(buffer,bmax,input) == NULL)
        { if (feof(input))
            { free(buffer);
              bmax = 0;
              return (NULL);
            }
          fprintf(stderr,"%s: Could not read next line, %d, of file %s\n",Prog_Name,nline,spath);
          exit (1);
        }
    }

  len = strlen(buffer);
  while (buffer[len-1] != '\n')
    { bmax = ((int) (1.4*bmax)) + 100;
      buffer = (char *) realloc(buffer,bmax);
      if (buffer == NULL)
        { fprintf(stderr,"%s: Out of memory reading %s\n",Prog_Name,spath);
          exit (1);
        }
      if (gzipd)
        { if (gzgets(input,buffer+len,bmax-len) == NULL)
            { if (gzeof(input))
                fprintf(stderr,"%s: Last line %d of file %s does not end with new-line\n",
                               Prog_Name,nline,spath);
              else
                fprintf(stderr,"%s: Could not read next line, %d, of file %s\n",
                               Prog_Name,nline,spath);
              exit (1);
            }
        }
      else
        { if (fgets(buffer+len,bmax-len,input) == NULL)
            { if (feof(input))
                fprintf(stderr,"%s: Last line %d of file %s does not end with new-line\n",
                               Prog_Name,nline,spath);
              else
                fprintf(stderr,"%s: Could not read next line, %d, of file %s\n",
                               Prog_Name,nline,spath);
              exit (1);
            }
        }
      len += strlen(buffer+len);
    }
  buffer[--len] = '\0';

  return (buffer);
}

//  Parse read range grabbing up to 4 values from it as follows:
//    type 1
//          val[0][.val[1]]
//    type 2
//          val[0][.val[1]] - val[2][.val[3]]
//    type 3
//          val[0][.val[1]] _ val[2] - val[3]
//  Return 0 if there is a syntactic error, otherwise return the type.
//  0 values indicate $, -1 indicates not present, -2 (only in val[1] and possibly val[3])
//      indicates val[0] and possibly val[2], respectively, are scaffold indices

static inline char *white(char *x)
{ while (isspace(*x))
    x += 1;
  return (x);
}

static inline char *getint(char *x, int *v)
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

static int range(char *src, int *v, Hash_Table *hash)
{ int   t, sep;
  char *x, *y, *z;
  int   i, j;

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
  i = Hash_Lookup(hash,x);
  if (i >= 0)
    { v[1] = -2;
      v[0] = i+1;
      return (1);
    }
  y = rindex(x,'_');
  if (y != NULL)
    { *y = '\0';
      i = Hash_Lookup(hash,x);
      if (i >= 0)
        { z = address(y+1,v+2,'-');
          if (z != NULL && *z == '\0')
            { v[1] = -2;
              v[0] = i+1;
              *y = '_';
              return (2);
            }
        }
      *y = '_';
    }
  for (y = index(x,'-'); y != NULL; y = index(y+1,'-'))
    { *y = '\0';
      i = Hash_Lookup(hash,x);
      if (i >= 0)
        { j = Hash_Lookup(hash,y+1);
          if (j >= 0)
            { v[1] = v[3] = -2;
              v[0] = i+1;
              v[2] = j+1;
              *y = '-';
              return (2);
            }
          else if (address(y+1,v+2,sep))
            { v[1] = -2;
              v[0] = i+1;
              *y = '-';
              return (2);
            }
        }
      else if (address(x,v,sep) == y)
        { j = Hash_Lookup(hash,y+1);
          if (j >= 0)
            { v[3] = -2;
              v[2] = j+1;
              *y = '-';
              return (2);
            }
        }
      *y = '-';
    }
  return (0);
}

 //  Interpret value pair s.c into a 1-based contig or scaffold index, accordingly

static int interpret_address(int s, int c, GDB *gdb)
{ if (c == -1)
    { if (s == 0)
        return (gdb->ncontig);
      if (s > gdb->ncontig)
        { fprintf(stderr,"%s: Absolute contig index '%d' is out of bounds\n",Prog_Name,s);
          exit (1);
        }
      return (s);
    }
  if (s == 0)
    s = gdb->nscaff;
  else if (s > gdb->nscaff)
    { fprintf(stderr,"%s: Scaffold index '%d' is out of bounds\n",Prog_Name,s);
      exit (1);
    }
  if (c == 0)
    return (gdb->scaffolds[s-1].ectg);
  if (c == -2)
    return (s);
  s -= 1;
  if (c > gdb->scaffolds[s].ectg-gdb->scaffolds[s].fctg)
    { fprintf(stderr,"%s: Relative contig index '%d' is out of bounds\n",Prog_Name,c);
      exit (1);
    }
  return (gdb->scaffolds[s].fctg + c);
}

  //  Interpret argument x as a range and map to contig level ranges

static void interpret_as_contigs(char *x, GDB *gdb, Hash_Table *hash,
                                 int *pbeg, int *pend, int *pfst, int *plst)
{ int t, vals[4], len, sunit;
  int beg, end, fst, lst;

  t = range(x,vals,hash);
  if (t == 0)
    { fprintf(stderr,"%s: Expression %s not a valid selection\n",Prog_Name,x);
      exit (1);
    }

  sunit = (vals[1] == -2);

  beg = end = interpret_address(vals[0],vals[1],gdb) - 1;
  fst = 0;
  if (t == 2)
    { end = interpret_address(vals[2],vals[3],gdb) - 1;
      if (beg > end)
        { fprintf(stderr,"%s: Object range in '%s' is empty\n",Prog_Name,x);
          exit (1);
        }
    }

  if (sunit)
    { beg = gdb->scaffolds[beg].fctg;
      end = gdb->scaffolds[end].ectg-1;
    }
  if (t <= 2)
    lst = gdb->contigs[end].clen;
  else // t == 3
    { fst = vals[2]; //  Type 3: intepret vals[2] and vals[3] as the substring interval
      lst = vals[3];
      if (sunit)
        len = gdb->scaffolds[gdb->contigs[end].scaf].slen;
      else
        len = gdb->contigs[beg].clen;
      if (lst == 0)
        lst = len;
      if (fst >= lst)
        { fprintf(stderr,"%s: Substring interval in '%s' is empty\n",Prog_Name,x);
          exit (1);
        }
      if (sunit)
        { for ( ; beg <= end; beg++)
            if (fst < gdb->contigs[beg].sbeg + gdb->contigs[beg].clen)
              break;
          fst -= gdb->contigs[beg].sbeg;
          if (fst < 0)
            fst = 0;
          for (; end >= beg; end--)
            if (lst > gdb->contigs[end].sbeg)
              break;
          lst -= gdb->contigs[end].sbeg;
          if (lst > gdb->contigs[end].clen)
            lst = gdb->contigs[end].clen;
        }
    }
  *pbeg = beg;
  *pend = end;
  *pfst = fst;
  *plst = lst;
}

static void interpret_as_selection(char *x, GDB *gdb, Hash_Table *hash, Selection *s)
{ int t, vals[4], len;

  t = range(x,vals,hash);
  if (t == 0)
    { fprintf(stderr,"%s: Expression %s not a valid selection\n",Prog_Name,x);
      exit (1);
    }

  s->beg = s->end = interpret_address(vals[0],vals[1],gdb) - 1; 
  if (vals[1] == -2)
    { len = gdb->scaffolds[s->end].slen;
      s->type = SCAFF_RANGE;
    }
  else
    { len = gdb->contigs[s->end].clen;
      s->type = CONTG_RANGE;
    }
  if (t < 3)
    { if (t == 2)
        { s->end = interpret_address(vals[2],vals[3],gdb)-1;
          if (s->beg > s->end)
             { fprintf(stderr,"%s: Object range in '%s' is empty\n",Prog_Name,x);
               exit (1);
             }
        }
      return;
    }
  s->src   = s->beg;  // == s->end
  s->beg   = vals[2];
  s->end   = vals[3];
  s->type += 1;
  if (s->end == 0)
    s->end = len;
  if (s->beg >= s->end)
    { fprintf(stderr,"%s: Substring interval in '%s' is empty\n",Prog_Name,x);
      exit (1);
    }
  return;
}

static void get_selection_contigs_from_file(FILE *fp, GDB *gdb, Hash_Table *hash,
                                            Contig_Range *chord, char *filename, int ordered)
{ char c, *p, *q, *line;
  int  i, pbeg, pend, pfst, plst;
  int  nline, order;

  order = 1;
  nline = 1;
  while ((line = read_line(fp,0,nline++,filename)) != NULL)
    { p = q = white(line);
      while (*q != '\0' && !isspace(*q))
        q += 1;
      if (p == q)      // empty line
        continue;

      c = *q;
      *q = '\0';
      interpret_as_contigs(p, gdb, hash, &pbeg, &pend, &pfst, &plst);
      *q = c;
#ifdef DEBUG_PARSE_SEQ
      fprintf(stderr, "%s: CTG_RANGE pbeg = %-12d pend = %-12d pfst = %-12d plst = %-12d\n",
              Prog_Name,pbeg,pend,pfst,plst);
#endif
      if (ordered)
        for (i = pbeg; i < pend; i++)
          if (chord[i].order)
            { fprintf(stderr,"%s: Overlapping ranges in selection expression\n",Prog_Name);
              exit (1);
            }

      for (i = pbeg+1; i < pend; i++)
        { chord[i].order = order;
          chord[i].beg   = 0;
          chord[i].end   = gdb->contigs[i].clen;
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
        }

      if (ordered)
        order += 1;
    }

  fclose(fp);
}

Contig_Range *get_selection_contigs(char *x, GDB *gdb, Hash_Table *hash, int ordered)
{ char         *e;
  int           i, pbeg, pend, pfst, plst;
  int           order;
  Contig_Range *chord;
  FILE         *fp;

  chord = (Contig_Range *) Malloc(sizeof(Contig_Range)*gdb->ncontig,"Allocating sequence array");

  if (x == NULL || strcmp(x,"@") == 0)
    { for (i = 0; i < gdb->ncontig; i++)
        { chord[i].order = 1;
          chord[i].beg   = 0;
          chord[i].end   = gdb->contigs[i].clen;
        }
    }
  else
    { if (*x == '\0')
        { fprintf(stderr,"%s: Empty range\n",Prog_Name);
          exit (1);
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

              interpret_as_contigs(x, gdb, hash, &pbeg, &pend, &pfst, &plst);
#ifdef DEBUG_PARSE_SEQ
              fprintf(stderr, "%s: CTG_RANGE pbeg = %-12d pend = %-12d pfst = %-12d plst = %-12d\n",
                      Prog_Name,pbeg,pend,pfst,plst);
#endif
              if (ordered)
                for (i = pbeg; i < pend; i++)
                  if (chord[i].order)
                    { fprintf(stderr,"%s: Overlapping ranges in selection expression\n",Prog_Name);
                      exit (1);
                    }

              for (i = pbeg+1; i < pend; i++)
                { chord[i].order = order;
                  chord[i].beg   = 0;
                  chord[i].end   = gdb->contigs[i].clen;
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
                }

              if (ordered)
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
    nline += 1;

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
      interpret_as_selection(line,gdb,hash,list+len);
      *q = c;
#ifdef DEBUG_PARSE_SEQ
      if (list[len].type % 2)
        fprintf(stderr, "%s: %c %d[%d..%d]\n",
              Prog_Name,list[len].type >= 2 ? '@' : ' ',list[len].src,list[len].beg,list[len].end);
      else
        fprintf(stderr, "%s: %c %d - %d\n",
              Prog_Name,list[len].type >= 2 ? '@' : ' ',list[len].beg,list[len].end);
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

  if (x == NULL || strcmp(x,"@") == 0)
    { list = (Selection *) Malloc(sizeof(Selection),"Allocating sequence array");
      *nlen = 1;
      list[0].beg = 0;
      if (x == NULL)
        { list[0].type = CONTG_RANGE;
          list[0].end  = gdb->ncontig-1;
        }
      else
        { list[0].type = SCAFF_RANGE;
          list[0].end  = gdb->nscaff-1;
        }
      return (list);
    }
  if (*x == '\0')
    { fprintf(stderr,"%s: Empty range\n",Prog_Name);
      exit (1);
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

       interpret_as_selection(x,gdb,hash,list+len);
#ifdef DEBUG_PARSE_SEQ
       if (list[len].type % 2)
         fprintf(stderr, "%s: %c %d[%d..%d]\n",
               Prog_Name,list[len].type >= 2 ? '@' : ' ',list[len].src,list[len].beg,list[len].end);
       else
         fprintf(stderr, "%s: %c %d - %d\n",
               Prog_Name,list[len].type >= 2 ? '@' : ' ',list[len].beg,list[len].end);
#endif
       len += 1;

       if (e != NULL)
         *e++ = DLIM_SYMBOL;
       x = e;
     }

  *nlen = len;
  return (list);
}
