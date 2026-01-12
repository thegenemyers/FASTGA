/*******************************************************************************************
 *
 *  Utility to convert a .1aln alignment file into a PSL formated file.
 *
 *  Author:    Gene Myers
 *  Creation:  Oct 2023
 *  Last Mod:  Feb 2023 - RD converted to .1aln (from .las)
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

#include "GDB.h"
#include "align.h"
#include "alncode.h"

#undef DEBUG_THREADS

static char *Usage = " [-T<int(8)>} <alignment:path>[.1aln]";

static int  NTHREADS;  // -T
static int  ISTWO;     // one gdb or two?

static int TSPACE;   // Trace spacing

//  THREAD ROUTINE TO GENERATE A SECTION OF THE DESIRED .PAF FILE


static inline void itoa(int x, char *buf, FILE *out)
{ char *s;

  s = buf;
  while (x >= 10)
    { *s++ = '0' + (x % 10);
      x /= 10;
    }
  fputc('0' + x, out);
  while (s > buf)
    fputc(*--s, out);
}

static inline void ltoa(int64 x, char *buf, FILE *out)
{ char *s;

  s = buf;
  while (x >= 10)
    { *s++ = '0' + (x % 10);
      x /= 10;
    }
  fputc('0' + x, out);
  while (s > buf)
    fputc(*--s, out);
}

static inline void stoa(char *s, FILE *out)
{ while (*s != '\0')
    fputc(*s++,out);
}

typedef struct
  { int64    beg;
    int64    end;
    GDB      gdb1;
    GDB      gdb2;
    OneFile *in;
    FILE    *out;
  } Packet;

void *gen_psl(void *args)
{ Packet *parm = (Packet *) args;

  int64    beg  = parm->beg;
  int64    end  = parm->end;
  GDB     *gdb1 = &(parm->gdb1);
  GDB     *gdb2 = &(parm->gdb2);
  OneFile *in   = parm->in;
  FILE    *out  = parm->out;

  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;

  uint16       *trace;
  int           tmax;
  GDB_CONTIG   *contig1, *contig2;
  GDB_SCAFFOLD *scaff1, *scaff2;
  char         *ahead, *bhead;
  int           acontig, bcontig;
  int           ascaff, bscaff;
  int           bmin, bmax;
  int           plen;
  int          *parr, *aarr, *barr;
  char         *aseq, *bseq, *bact;
  char          buf[32];
  int64         aoff, boff;
  Path         *path;
  Work_Data    *work;

  if (beg >= end) return (NULL);

  work = New_Work_Data();
  aseq = New_Contig_Buffer(gdb1);
  bseq = New_Contig_Buffer(gdb2);
  if (aseq == NULL || bseq == NULL)
    exit (1);

  contig1 = gdb1->contigs;
  contig2 = gdb2->contigs;
  scaff1  = gdb1->scaffolds;
  scaff2  = gdb2->scaffolds;
  ahead   = gdb1->headers;
  bhead   = gdb2->headers;

  plen = 10000;
  parr = (int *) Malloc(sizeof(int)*plen*3,"Allocating block position array");
  if (parr == NULL)
    exit (1);
  aarr = parr + plen;
  barr = aarr + plen;

  aln->path = path = &(ovl->path);
  aln->aseq = aseq;

  tmax  = in->info['T']->given.max;
  trace = (uint16 *) Malloc(2*sizeof(uint16)*tmax,"Allocating trace vector");
  if (trace == NULL)
    exit (1);
  path->trace = (void *) trace;

  //  For each alignment do

  if (!oneGoto(parm->in,'A',beg+1))
    { fprintf(stderr,"%s: Can't locate to object %lld in aln file\n",Prog_Name,beg+1);
      exit (1);
    }
  oneReadLine(parm->in);

  aoff = 0;
  for (acontig = -1; beg < end; beg++)
    { Read_Aln_Overlap(in,ovl);
      path->tlen  = Read_Aln_Trace(in,(uint8 *) trace, NULL);
      path->trace = trace;

      Decompress_TraceTo16(ovl);

      if (acontig != ovl->aread)
        { acontig = ovl->aread;
          Get_Contig(gdb1,acontig,NUMERIC,aseq);
          aln->alen = contig1[acontig].clen;
          aoff      = contig1[acontig].sbeg;
        }
      bcontig = ovl->bread;
      aln->blen  = contig2[bcontig].clen;
      aln->flags = ovl->flags;

      ascaff = contig1[acontig].scaf;
      bscaff = contig2[bcontig].scaf;

      if (COMP(aln->flags))
        boff = contig2[bcontig].sbeg + contig2[bcontig].clen;
      else
        boff = contig2[bcontig].sbeg;

      if (COMP(aln->flags))
        { bmin = (aln->blen-path->bepos);
          if (bmin < 0) bmin = 0;
          bmax = (aln->blen-path->bbpos);
          if (bmax > aln->blen) bmax = aln->blen;
        }
      else
        { bmin = path->bbpos;
          if (bmin < 0) bmin = 0;
          bmax = path->bepos;
          if (bmax > aln->blen) bmax = aln->blen;
        }

      bact = Get_Contig_Piece(gdb2,bcontig,bmin,bmax,NUMERIC,bseq);
      if (COMP(aln->flags))
        { Complement_Seq(bact,bmax-bmin);
          aln->bseq = bact - (aln->blen-bmax);
        }
      else
        aln->bseq = bact - bmin; 

      Compute_Trace_PTS(aln,work,TSPACE,GREEDIEST,1,-1);
      Gap_Improver(aln,work);

      { int     i, j, x, p, q;
        int    *t, T;
        int     M;
        int     I, D, S, X;
        int     IB, DB;
        int     bcnt, bmat;

        t = (int *) path->trace;
        T = path->tlen;

        { // trim trailing INDEL
          // better fix needed
          int trim;
          
          trim = 0;
          while (t[T-1] == -path->aepos-1)
            { trim++;
              T--;
            }
          if (trim > 0)
            { // deletion / target insertion
              // trim target
              path->bepos -= trim;
              path->diffs -= trim;
            }
          
          trim = 0;
          while (t[T-1] == path->bepos+1)
            { trim++;
              T--;
            }
          if (trim > 0)
            { // insertion / target deletion
              // trim query
              path->aepos -= trim;
              path->diffs -= trim;
            }
        }

        M = path->aepos - path->abpos;
        //N = path->bepos - path->bbpos;
        I = D = 0;
        IB = DB = 0;
        p = 0; 
        for (x = 0; x < T; x++)
          { q = p;
            if ((p = t[x]) < 0)
              { I += 1;
                if (p != q)
                  IB += 1;
              }
            else
              { D += 1;
                if (p != q)
                  DB += 1;
              }
          }
        S = path->diffs - (I+D);
        //X = (M+N - (I+D+2*S))/2;
        X = M - D - S; // N - I - S

        itoa(X,buf,out);
        fputc('\t',out);
        itoa(S,buf,out);
        stoa("\t0\t0\t",out);
        itoa(DB,buf,out);
        fputc('\t',out);
        itoa(D,buf,out);
        fputc('\t',out);
        itoa(IB,buf,out);
        fputc('\t',out);
        itoa(I,buf,out);
        fputc('\t',out);
        fputc(COMP(ovl->flags)?'-':'+',out);

        fputc('\t',out);
        stoa(ahead+scaff1[ascaff].hoff,out);
        fputc('\t',out);
        ltoa(scaff1[ascaff].slen,buf,out);
        fputc('\t',out);
        ltoa(aoff+path->abpos,buf,out);
        fputc('\t',out);
        ltoa(aoff+path->aepos,buf,out);

        fputc('\t',out);
        stoa(bhead+scaff2[bscaff].hoff,out);
        fputc('\t',out);
        ltoa(scaff2[bscaff].slen,buf,out);
        fputc('\t',out);
        if (COMP(aln->flags))
          { boff = contig2[bcontig].sbeg + contig2[bcontig].clen;
            ltoa(boff-path->bepos,buf,out);
            fputc('\t',out);
            ltoa(boff-path->bbpos,buf,out);
          }
        else
          { boff = contig2[bcontig].sbeg;
            ltoa(boff+path->bbpos,buf,out);
            fputc('\t',out);
            ltoa(boff+path->bepos,buf,out);
          }
        bcnt = 0;
        i = path->abpos+1;
        j = path->bbpos+1;
        for (x = 0; x < T; x++)
          { if ((p = t[x]) < 0)
              { bmat = -(p+i);
                i += bmat;
                j += bmat+1;
              }
            else
              { bmat = p-j;
                i += bmat+1;
                j += bmat;
              }
            if (bmat > 0)
              bcnt += 1;
          }
        bmat = (path->aepos - i)+1;
        if (bmat > 0)
          bcnt += 1;
        fprintf(out,"\t%d\t",bcnt);

        if (bcnt > plen)
          { plen = ((int) (1.2*bcnt)) + 1000;
            parr = (int *) Realloc(parr, sizeof(int)*plen*3, "Reallocating block position array");
            if (parr == NULL)
              exit (1);
            aarr = parr + plen;
            barr = aarr + plen;
        }

        bcnt = 0;
        i = path->abpos+1;
        j = path->bbpos+1;
        for (x = 0; x < T; x++)
          { if ((p = t[x]) < 0)
              { bmat = -(p+i);
                if (bmat > 0)
                  { aarr[bcnt] = i-1;
                    barr[bcnt] = j-1;
                  }
                i += bmat;
                j += bmat+1;
              }
            else
              { bmat = p-j;
                if (bmat > 0)
                  { aarr[bcnt] = i-1;
                    barr[bcnt] = j-1;
                  }
                i += bmat+1;
                j += bmat;
              }
            if (bmat > 0)
              parr[bcnt++] = bmat;
          }
        bmat = (path->aepos - i)+1;
        if (bmat > 0)
          { aarr[bcnt] = i-1;
            barr[bcnt] = j-1;
            parr[bcnt++] = bmat;
          }
        
        if (COMP(aln->flags))
          for (i = bcnt-1; i >= 0; i--)
            { itoa(parr[i],buf,out);
              fputc(',',out);
            }
        else
          for (i = 0; i < bcnt; i++)
            { itoa(parr[i],buf,out);
              fputc(',',out);
            }
        fputc('\t',out);

        if (COMP(aln->flags))
          for (i = bcnt-1; i >= 0; i--)
            { ltoa(scaff1[ascaff].slen-(aoff+aarr[i]+parr[i]),buf,out);
              fputc(',',out);
            }
        else
          for (i = 0; i < bcnt; i++)
            { ltoa(aoff+aarr[i],buf,out);
              fputc(',',out);
            }
        fputc('\t',out);
        
        if (COMP(aln->flags))
          { boff = contig2[bcontig].sbeg + contig2[bcontig].clen;
            for (i = bcnt-1; i >= 0; i--)
              { ltoa(boff-(barr[i]+parr[i]),buf,out);
                fputc(',',out);
              }
          }
        else
          { boff = contig2[bcontig].sbeg;
            for (i = 0; i < bcnt; i++)
              { ltoa(boff+barr[i],buf,out);
                fputc(',',out);
              }
          }
        fputc('\n',out);
      }
    }

  free(parr);
  free(trace);
  free(bseq-1);
  free(aseq-1);
  Free_Work_Data(work);

  return (NULL);
}

int main(int argc, char *argv[])
{ Packet    *parm;
  GDB       _gdb1, *gdb1 = &_gdb1;
  GDB       _gdb2, *gdb2 = &_gdb2;
  FILE     **units1;
  FILE     **units2;
  OneFile   *input;
  int64      novl; // RD need this out of scope below so I can use it when setting up threads

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("ALNtoPSL")

    (void) flags;

    NTHREADS = 8;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
        exit (1);
      }

    parm = Malloc(sizeof(Packet)*NTHREADS,"Allocating thread records");
    if (parm == NULL)
      exit (1);
  }

  //  Initiate .1aln file reading and read header information

  { char       *pwd, *root, *cpath;
    char       *src1_name, *src2_name;
    char       *head, *sptr, *eptr;
    int         s;

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".1aln");
    input = open_Aln_Read(Catenate(pwd,"/",root,".1aln"),NTHREADS,&novl,&TSPACE,
                          &src1_name,&src2_name,&cpath);
    if (input == NULL)
      exit (1);
    free(root);
    free(pwd);


    ISTWO = (src2_name != NULL);

    Skip_Skeleton(input);
    units1 = Get_GDB(gdb1,src1_name,cpath,NTHREADS,NULL);
    if (ISTWO)
      { Skip_Skeleton(input);
        units2 = Get_GDB(gdb2,src2_name,cpath,NTHREADS,NULL);
      }
    else
      { gdb2   = gdb1;
        units2 = units1;
      }

    free(src1_name);
    free(src2_name);
    free(cpath);

    head  = gdb1->headers;
    for (s = 0; s < gdb1->nscaff; s++)
      { sptr = head + gdb1->scaffolds[s].hoff;
        for (eptr = sptr; *eptr != '\0'; eptr++)
          if (isspace(*eptr))
            break;
        *eptr = '\0';
      }

    if (ISTWO)
      { head  = gdb2->headers;
        for (s = 0; s < gdb2->nscaff; s++)
          { sptr = head + gdb2->scaffolds[s].hoff;
            for (eptr = sptr; *eptr != '\0'; eptr++)
              if (isspace(*eptr))
                break;
            *eptr = '\0';
          }
      }
  }

  //  Divide .1aln into NTHREADS parts

  { int p;

    for (p = 0; p < NTHREADS; p++)
      { parm[p].beg = (p * novl) / NTHREADS;
        if (p > 0)
          parm[p-1].end = parm[p].beg;
      }
    parm[NTHREADS-1].end = novl;
  }

  //   Use NTHREADS to produce .paf for each part and then cat the parts to stdout

  { int   p, x;
    char *oprefix;
    char *buffer;
#ifndef DEBUG_THREADS
    pthread_t threads[NTHREADS];
#endif

    //  Setup thread packets

    oprefix = Strdup(Numbered_Suffix("_psl.",getpid(),"."),"Allocating temp file prefix");

    for (p = 0; p < NTHREADS; p++)
      { parm[p].gdb1 = *gdb1;
        parm[p].gdb2 = *gdb2;
        parm[p].gdb1.seqs = units1[p];
        parm[p].gdb2.seqs = units2[p];
        parm[p].in  = input + p ;
        parm[p].out = fopen(Numbered_Suffix(oprefix,p,".psl"),"w+");
        if (parm[p].out == NULL)
          { fprintf(stderr,"%s: Cannot open %s.%d.paf for reading & writing\n",
                           Prog_Name,oprefix,p);
            exit (1);
          }
        unlink(Numbered_Suffix(oprefix,p,".psl"));
      }

      Numbered_Suffix(NULL,0,NULL);
      free(oprefix);

    //  Launch and then gather threads

#ifdef DEBUG_THREADS
    for (p = 0; p < NTHREADS; p++)
      gen_psl(parm+p);
#else
    for (p = 1; p < NTHREADS; p++)
      pthread_create(threads+p,NULL,gen_psl,parm+p);
    gen_psl(parm);
    for (p = 1; p < NTHREADS; p++)
      pthread_join(threads[p],NULL);
#endif

    if (NTHREADS > 1)
      { for (p = 1; p < NTHREADS; p++)
          { fclose(units1[p]);
            if (ISTWO)
              fclose(units2[p]);
          }
        free(units1);
        if (ISTWO)
	  free(units2);
      }

    //  Concatenate thread generated paf parts to stdout

#define PUSH_BLOCK 0x100000

    buffer = Malloc(PUSH_BLOCK,"Allocating seek block");

    for (p = 0; p < NTHREADS; p++)
      { FILE *o = parm[p].out;

        rewind(o);
        x = PUSH_BLOCK;
        while (x >= PUSH_BLOCK)
          { x = fread(buffer,1,PUSH_BLOCK,o);
            if (x > 0)
              fwrite(buffer,x,1,stdout);
          }
        fclose(o);
      }

    free(buffer);
  }

  Close_GDB(gdb1);
  if (ISTWO)
    Close_GDB(gdb2);

  oneFileClose(input);
  
  free(Command_Line);
  free(Prog_Name);
  Catenate(NULL,NULL,NULL,NULL);

  exit (0);
}
