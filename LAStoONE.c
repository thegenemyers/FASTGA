/*******************************************************************************************
 *
 *  Utility for displaying the information in the overlaps of a .las file in a very
 *    simple to parse format.
 *
 *  Author:    Gene Myers
 *  Creation:  July 2013
 *  Last Mod:  Jan 2015
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

#include "DB.h"
#include "align.h"
#include "ONElib.h"

static char *Usage =
    "[-cto] <src1:db|dam> [<src2:db|dam>] <align:las> [<reads:FILE> | <reads:range> ...]";

static char *One_Schema =
  "P 3 dal                  This is a 1-code las file from daligner\n"
  "D X 1 3 INT              Data prolog: trace spacing\n"
  "O P 2 3 INT 8 INT_LIST   A-read and B-read list\n"
  ".                        All per B-read/l.a.s:\n"
  "D O 1 6 STRING               Orientation [+-]\n"
  "D C 1 6 STRING               Chain directive [>+-.]\n"
  "D A 1 8 INT_LIST             A-read alignment intervals: (ab,ae)\n"
  "D B 1 8 INT_LIST             B-read alignment intervals: (bb,be)\n"
  "D D 1 8 INT_LIST             Diffs in alignment\n"
  "D L 2 3 INT 8 INT_LIST    A-read length and then each B-read length\n"
  ".                         One line per B-read:\n"
  "D T 1 8 INT_LIST              Trace segment length\n"
  "D Q 1 8 INT_LIST              Trace segment diffs\n";

static Overlap   *ovls;
static uint16    *trace;
static int64     *list;
static char      *string;
static OneFile   *file1;
static DAZZ_READ *read1, *read2;

static int     OVERLAP;
static int     DOCOORDS;
static int     DOTRACE;

static void output_pile(Overlap *optr)
{ int i, k;
  Overlap *o;

  i = 0;
  for (o = ovls; o < optr; o++)
    list[i++] = o->bread+1;
  oneInt(file1,0) = ovls->aread+1;
  oneWriteLine(file1,'P',i,list);

  i = 0;
  for (o = ovls; o < optr; o++)
    if (COMP(o->flags))
      string[i++] = 'c';
    else
      string[i++] = 'n';
  oneWriteLine(file1,'O',i,string);

  i = 0;
  for (o = ovls; o < optr; o++)
    if (CHAIN_NEXT(o->flags))
      string[i++] = '-';
    else if (BEST_CHAIN(o->flags))
      string[i++] = '>';
    else if (CHAIN_START(o->flags))
      string[i++] = '+';
    else
      string[i++] = '.';
  oneWriteLine(file1,'C',i,string);
    
  if (DOCOORDS)
    { i = 0;
      for (o = ovls; o < optr; o++)
        { list[i++] = o->path.abpos;
          list[i++] = o->path.aepos;
        }
      oneWriteLine(file1,'A',i,list);

      i = 0;
      for (o = ovls; o < optr; o++)
        { list[i++] = o->path.bbpos;
          list[i++] = o->path.bepos;
        }
      oneWriteLine(file1,'B',i,list);

      oneInt(file1,0) = read1[ovls->aread].rlen;
      i = 0;
      for (o = ovls; o < optr; o++)
        list[i++] = read2[o->bread].rlen;
      oneWriteLine(file1,'L',i,list);

      i = 0;
      for (o = ovls; o < optr; o++)
        list[i++] = o->path.diffs;
      oneWriteLine(file1,'D',i,list);
    }
    
  if (DOTRACE)
    { uint16 *trace;
      int     tlen;

      for (o = ovls; o < optr; o++)
        { trace = (uint16 *) o->path.trace;
          tlen  = o->path.tlen;

          i = 0;
          for (k = 1; k < tlen; k += 2)
            list[i++] = trace[k];
          oneWriteLine(file1,'T',i,list);

          i = 0;
          for (k = 0; k < tlen; k += 2)
            list[i++] = trace[k];
          oneWriteLine(file1,'Q',i,list);
        }
    }
}

static int ORDER(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);
  return (x-y);
}

int main(int argc, char *argv[])
{ DAZZ_DB   _db1, *db1 = &_db1; 
  DAZZ_DB   _db2, *db2 = &_db2; 
  OneSchema *schema;
  char      *command;

  FILE   *input;
  int64   novl, omax, tmax;
  int     tspace, tbytes;
  int     reps, *pts;
  int     input_pts;

  int     ISTWO;

  //  Process options and capture command line for provenance

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("LA2ONE")

    { int   n, t;
      char *c;

      n = 0;
      for (t = 1; t < argc; t++)
        n += strlen(argv[t])+1;

      command = Malloc(n+1,"Allocating command string");
      if (command == NULL)
        exit (1);

      c = command;
      if (argc >= 1)
        { c += sprintf(c,"%s",argv[1]);
          for (t = 2; t < argc; t++)
            c += sprintf(c," %s",argv[t]);
        }
      *c = '\0';
    }

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("cto")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    OVERLAP   = flags['o'];
    DOCOORDS  = flags['c'];
    DOTRACE   = flags['t'];

    if (DOTRACE)
      DOCOORDS = 1;

    if (argc <= 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      Output pile reads, orientation, and chains by default");
        fprintf(stderr," (P, O, C lines)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -c: Ootput also aligned intervals, read lengths, and diffs");
          fprintf(stderr," (B, E, L, and D lines)\n");
        fprintf(stderr,"      -t: Output also traces (T and Q lines)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -o: Output proper overlaps only\n");

        exit (1);
      }
  }

  //  Open trimmed DB or DB pair

  { int   status;
    char *pwd, *root;
    FILE *input;

    ISTWO  = 0;
    status = Open_DB(argv[1],db1);
    if (status < 0)
      exit (1);
    if (db1->part > 0)
      { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[1]);
        exit (1);
      }

    if (argc > 3)
      { pwd   = PathTo(argv[3]);
        root  = Root(argv[3],".las");
        if ((input = fopen(Catenate(pwd,"/",root,".las"),"r")) != NULL)
          { ISTWO = 1;
            fclose(input);
            status = Open_DB(argv[2],db2);
            if (status < 0)
              exit (1);
            if (db2->part > 0)
              { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[2]);
                exit (1);
              }
            Trim_DB(db2);
          }
        else
          db2 = db1;
        free(root);
        free(pwd);
      }
    else
      db2 = db1;
    Trim_DB(db1);
  }

  //  Process read index arguments into a sorted list of read ranges

  input_pts = 0;
  if (argc == ISTWO+4)
    { if (argv[ISTWO+3][0] != LAST_READ_SYMBOL || argv[ISTWO+3][1] != '\0')
        { char *eptr, *fptr;
          int   b, e;

          b = strtol(argv[ISTWO+3],&eptr,10);
          if (eptr > argv[ISTWO+3] && b > 0)
            { if (*eptr == '-')
                { if (eptr[1] != LAST_READ_SYMBOL || eptr[2] != '\0')
                    { e = strtol(eptr+1,&fptr,10);
                      input_pts = (fptr <= eptr+1 || *fptr != '\0' || e <= 0);
                    }
                }
              else
                input_pts = (*eptr != '\0');
            }
          else
            input_pts = 1;
        }
    }

  if (input_pts)
    { int v, x;
      FILE *input;

      input = Fopen(argv[ISTWO+3],"r");
      if (input == NULL)
        exit (1);

      reps = 0;
      while ((v = fscanf(input," %d",&x)) != EOF)
        if (v == 0)
          { fprintf(stderr,"%s: %d'th item of input file %s is not an integer\n",
                           Prog_Name,reps+1,argv[2]);
            exit (1);
          }
        else
          reps += 1;

      reps *= 2;
      pts   = (int *) Malloc(sizeof(int)*reps,"Allocating read parameters");
      if (pts == NULL)
        exit (1);

      rewind(input);
      for (v = 0; v < reps; v += 2)
        { fscanf(input," %d",&x);
          pts[v] = pts[v+1] = x;
        }

      fclose(input);
    }

  else
    { pts  = (int *) Malloc(sizeof(int)*2*argc,"Allocating read parameters");
      if (pts == NULL)
        exit (1);

      reps = 0;
      if (argc > 3+ISTWO)
        { int   c, b, e;
          char *eptr, *fptr;

          for (c = 3+ISTWO; c < argc; c++)
            { if (argv[c][0] == LAST_READ_SYMBOL)
                { b = db1->nreads;
                  eptr = argv[c]+1;
                }
              else
                b = strtol(argv[c],&eptr,10);
              if (eptr > argv[c])
                { if (b <= 0)
                    { fprintf(stderr,"%s: %d is not a valid index\n",Prog_Name,b);
                      exit (1);
                    }
                  if (*eptr == '\0')
                    { pts[reps++] = b;
                      pts[reps++] = b;
                      continue;
                    }
                  else if (*eptr == '-')
                    { if (eptr[1] == LAST_READ_SYMBOL)
                        { e = INT32_MAX;
                          fptr = eptr+2;
                        }
                      else
                        e = strtol(eptr+1,&fptr,10);
                      if (fptr > eptr+1 && *fptr == 0 && e > 0)
                        { pts[reps++] = b;
                          pts[reps++] = e;
                          if (b > e)
                            { fprintf(stderr,"%s: Empty range '%s'\n",Prog_Name,argv[c]);
                              exit (1);
                            }
                          continue;
                        }
                    }
                }
              fprintf(stderr,"%s: argument '%s' is not an integer range\n",Prog_Name,argv[c]);
              exit (1);
            }

          qsort(pts,reps/2,sizeof(int64),ORDER);

          b = 0;
          for (c = 0; c < reps; c += 2)
            if (b > 0 && pts[b-1] >= pts[c]-1) 
              { if (pts[c+1] > pts[b-1])
                  pts[b-1] = pts[c+1];
              }
            else
              { pts[b++] = pts[c];
                pts[b++] = pts[c+1];
              }
          pts[b++] = INT32_MAX;
          reps = b;
        }
      else
        { pts[reps++] = 1;
          pts[reps++] = INT32_MAX;
        }
    }

  //  Initiate file reading and read header
  
  { char  *over, *pwd, *root;

    pwd   = PathTo(argv[2+ISTWO]);
    root  = Root(argv[2+ISTWO],".las");
    over  = Catenate(pwd,"/",root,".las");
    input = Fopen(over,"r");
    if (input == NULL)
      exit (1);

    if (fread(&novl,sizeof(int64),1,input) != 1)
      SYSTEM_READ_ERROR
    if (fread(&tspace,sizeof(int),1,input) != 1)
      SYSTEM_READ_ERROR

    if (tspace <= TRACE_XOVR && tspace != 0)
      tbytes = sizeof(uint8);
    else
      tbytes = sizeof(uint16);

    free(pwd);
    free(root);
  }

  schema = oneSchemaCreateFromText(One_Schema);
  file1  = oneFileOpenWriteNew("-",schema,"dal",true,1);
  oneAddProvenance(file1,Prog_Name,"1.0","%s >?.dal",command);

  //  Scan to determine max trace length and max pile size

  { int     in, npt, idx;
    int     j, ar, al;
    int     tlen;
    int64   odeg;
    Overlap _ovl, *ovl = &_ovl;

    in  = 0;
    npt = pts[0];
    idx = 1;

    //  For each record do

    omax = tmax = 0;
    odeg = 0;

    al = 0;
    for (j = 0; j < novl; j++)

       //  Read it in

      { Read_Overlap(input,ovl);
        tlen = ovl->path.tlen;
        fseeko(input,tlen*tbytes,SEEK_CUR);

        //  Determine if it should be displayed

        ar = ovl->aread+1;
        if (in)
          { while (ar > npt)
              { npt = pts[idx++];
                if (ar < npt)
                  { in = 0;
                    break;
                  }
                npt = pts[idx++];
              }
          }
        else
          { while (ar >= npt)
              { npt = pts[idx++];
                if (ar <= npt)
                  { in = 1;
                    break;
                  }
                npt = pts[idx++];
              }
          }
        if (!in)
          continue;

        //  If -o check display only overlaps

        if (OVERLAP)
          { if (ovl->path.abpos != 0 && ovl->path.bbpos != 0)
              continue;
            if (ovl->path.aepos != db1->reads[ovl->aread].rlen &&
                ovl->path.bepos != db2->reads[ovl->bread].rlen)
              continue;
          }

        if (ar != al)
          { if (odeg > omax)
              omax = odeg;
            al = ar;
          }
        odeg += 1;
        if (tlen > tmax)
          tmax = tlen;
      }
    if (odeg > omax)
      omax = odeg;
  }

  //  Read the file and display selected records
  
  { int        j;
    Overlap   *optr;
    uint16    *tptr;
    int        in, npt, idx, ar, last;

    rewind(input);
    fread(&novl,sizeof(int64),1,input);
    fread(&tspace,sizeof(int),1,input);

    ovls   = Malloc(sizeof(Overlap)*omax,"Allocating alignment array");
    trace  = Malloc(sizeof(uint16)*omax*tmax,"Allocating trace buffer");
    string = Malloc(sizeof(int64)*omax,"Allocating 1-string");
    if (tmax > 2*omax)
      list = Malloc(sizeof(int64)*tmax,"Allocating 1-list");
    else
      list = Malloc(sizeof(int64)*omax*2,"Allocating 1-list");
    if (ovls == NULL || trace == NULL || string == NULL || list == NULL)
      exit (1);

    read1 = db1->reads;
    read2 = db2->reads;

    if (DOTRACE)
      { oneInt(file1,0) = tspace;
        oneWriteLine(file1,'X',0,NULL);
      }

    //  For each record do

    in  = 0;
    npt = pts[0];
    idx = 1;

    optr = ovls;
    tptr = trace; 
    last = -1;
    for (j = 0; j < novl; j++)

       //  Read it in

      { Read_Overlap(input,optr);
        ar = optr->aread+1;

        if (in)

          { if (ar == last)
              { optr->path.trace = (void *) tptr;
                Read_Trace(input,optr,tbytes);
                if (tbytes == 1)
                  Decompress_TraceTo16(optr);
                tptr += sizeof(uint16)*optr->path.tlen;
                optr += 1;
              }

            else
              { output_pile(optr);

                while (ar > npt)
                  { npt = pts[idx++];
                    if (ar < npt)
                      { in = 0;
                        break;
                      }
                    npt = pts[idx++];
                  }

                if (in)
                  { ovls[0] = *optr++;
                    tptr    = trace;
                    optr    = ovls;
                    last    = ar;

                    optr->path.trace = (void *) tptr;
                    Read_Trace(input,optr,tbytes);
                    if (tbytes == 1)
                      Decompress_TraceTo16(optr);
                    tptr += sizeof(uint16)*optr->path.tlen;
                    optr += 1;
                  }
                else
                  { fseeko(input,optr->path.tlen*tbytes,SEEK_CUR);
                    optr = ovls;
                    tptr = trace;
                  }
              }
          }

        else
          { while (ar >= npt)
              { npt = pts[idx++];
                if (ar <= npt)
                  { in = 1;
                    break;
                  }
                npt = pts[idx++];
              }

            if (in)
              { last = ar;

                optr->path.trace = (void *) tptr;
                Read_Trace(input,optr,tbytes);
                if (tbytes == 1)
                  Decompress_TraceTo16(optr);
                tptr += sizeof(uint16)*optr->path.tlen;
                optr += 1;
              }
            else
              fseeko(input,optr->path.tlen*tbytes,SEEK_CUR);
          }
      }

    if (in)
      output_pile(optr);

    free(string);
    free(list);
    free(trace);
    free(ovls);
  }

  oneFileClose(file1);
  oneSchemaDestroy(schema);

  Close_DB(db1);
  if (ISTWO)
    Close_DB(db2);

  free(command);

  exit (0);
}
