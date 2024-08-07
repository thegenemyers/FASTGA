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
  char         *aseq, *bseq, *bact;
  int64         aoff, boff;
  Path         *path;
  Work_Data    *work;

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
      path->tlen  = Read_Aln_Trace(in,(uint8 *) trace);
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

      Compute_Trace_PTS(aln,work,TSPACE,GREEDIEST);
      Gap_Improver(aln,work);

      { int     i, j, x, p, q;
        int    *t, T;
        int     M, N;
        int     I, D, S, X;
        int     IB, DB;
        int     bcnt, bmat;

        t = (int *) path->trace;
        T = path->tlen;

        M = path->aepos - path->abpos;
        N = path->bepos - path->bbpos;
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
        X = (M+N - (I+D+2*S))/2;

        fprintf(out,"%d\t%d\t0\t0\t%d\t%d\t%d\t%d\t%c",X,S,IB,I,DB,D,COMP(ovl->flags)?'-':'+');
        fprintf(out,"\t%s\t%lld\t%lld\t%lld",
                  ahead+scaff1[ascaff].hoff,scaff1[ascaff].slen,aoff+path->abpos,aoff+path->aepos);
        if (COMP(aln->flags))
          fprintf(out,"\t%s\t%lld\t%lld\t%lld",
                   bhead+scaff2[bscaff].hoff,scaff2[bscaff].slen,boff-path->bepos,boff-path->bbpos);
        else
          fprintf(out,"\t%s\t%lld\t%lld\t%lld",
                   bhead+scaff2[bscaff].hoff,scaff2[bscaff].slen,boff+path->bbpos,boff+path->bepos);

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
              fprintf(out,"%d,",bmat);
          }
        bmat = (path->aepos - i)+1;
        if (bmat > 0)
          fprintf(out,"%d,",bmat);
        fprintf(out,"\t");

        i = path->abpos+1;
        j = path->bbpos+1;
        for (x = 0; x < T; x++)
          { if ((p = t[x]) < 0)
              { bmat = -(p+i);
                if (bmat > 0)
                  fprintf(out,"%d,",i);
                i += bmat;
                j += bmat+1;
              }
            else
              { bmat = p-j;
                if (bmat > 0)
                  fprintf(out,"%d,",i);
                i += bmat+1;
                j += bmat;
              }
          }
        bmat = (path->aepos - i)+1;
        if (bmat > 0)
          fprintf(out,"%d,",i);
        fprintf(out,"\t");

        i = path->abpos+1;
        j = path->bbpos+1;
        for (x = 0; x < T; x++)
          { if ((p = t[x]) < 0)
              { bmat = -(p+i);
                if (bmat > 0)
                  fprintf(out,"%d,",j);
                i += bmat;
                j += bmat+1;
              }
            else
              { bmat = p-j;
                if (bmat > 0)
                  fprintf(out,"%d,",j);
                i += bmat+1;
                j += bmat;
              }
          }
        bmat = (path->aepos - i)+1;
        if (bmat > 0)
          fprintf(out,"%d,",j);
        fprintf(out,"\n");
      }
    }

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
    char       *spath, *tpath;
    int         type;
    FILE       *test;

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".1aln");
    input = open_Aln_Read(Catenate(pwd,"/",root,".1aln"),NTHREADS,
                          &novl,&TSPACE,&src1_name,&src2_name,&cpath);
    if (input == NULL)
      exit (1);
    free(root);
    free(pwd);

    test = fopen(src1_name,"r");
    if (test == NULL)
      { if (*src1_name != '/')
          test = fopen(Catenate(cpath,"/",src1_name,""),"r");
        if (test == NULL)
          { fprintf(stderr,"%s: Could not find GDB %s\n",Prog_Name,src1_name);
            exit (1);
          }
        pwd = Strdup(Catenate(cpath,"/",src1_name,""),"Allocating expanded name");
        free(src1_name);
        src1_name = pwd;
      }
    fclose(test);

    if (src2_name != NULL)
      { test = fopen(src2_name,"r");
        if (test == NULL)
          { if (*src2_name != '/')
              test = fopen(Catenate(cpath,"/",src2_name,""),"r");
            if (test == NULL)
              { fprintf(stderr,"%s: Could not find GDB %s\n",Prog_Name,src2_name);
                exit (1);
              }
            pwd = Strdup(Catenate(cpath,"/",src2_name,""),"Allocating expanded name");
            free(src2_name);
            src2_name = pwd;
          }
        fclose(test);
      }

    free(cpath);

    //  Prepare GDBs from sources if necessary

    units1 = NULL;
    units2 = NULL;
    ISTWO = 0;
    type  = Get_GDB_Paths(src1_name,NULL,&spath,&tpath,0);
    if (type != IS_GDB)
      units1 = Create_GDB(gdb1,spath,type,NTHREADS,NULL);
    else
      { Read_GDB(gdb1,tpath);
        if (gdb1->seqs == NULL)
          { fprintf(stderr,"%s: GDB %s must have sequence data\n",Prog_Name,tpath);
            exit (1);
          }
      }
    free(spath);
    free(tpath);

    if (src2_name != NULL)
      { type = Get_GDB_Paths(src2_name,NULL,&spath,&tpath,0);
        if (type != IS_GDB)
          units2 = Create_GDB(gdb2,spath,type,NTHREADS,NULL);
        else
          { Read_GDB(gdb2,tpath);
            if (gdb2->seqs == NULL)
              { fprintf(stderr,"%s: GDB %s must have sequence data\n",Prog_Name,tpath);
                exit (1);
              }
          }
        free(spath);
        free(tpath);
        ISTWO = 1;
      }
    else
      gdb2 = gdb1;

    free(src1_name);
    free(src2_name);
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
        if (p > 0)
          { if (units1 != NULL)
              parm[p].gdb1.seqs = units1[p];
            else
              { parm[p].gdb1.seqs = fopen(gdb1->seqpath,"r");
                if (parm[p].gdb1.seqs == NULL)
                  { fprintf(stderr,"%s: Cannot open another copy of GDB %s\n",
                                   Prog_Name,gdb1->seqpath);
                    exit (1);
                  }
              }
            if (ISTWO)
              { if (units2 != NULL)
                  parm[p].gdb2.seqs = units2[p];
                else
                  { parm[p].gdb2.seqs = fopen(gdb2->seqpath,"r");
                    if (parm[p].gdb2.seqs == NULL)
                      { fprintf(stderr,"%s: Cannot open another copy of GDB %s\n",
                                       Prog_Name,gdb2->seqpath);
                        exit (1);
                      }
                  }
              }
            else
              { if (units1 != NULL)
                  parm[p].gdb2.seqs = units1[p];
                else
                  parm[p].gdb2.seqs = parm[p].gdb1.seqs;
              }
          }
        parm[p].in  = input + p ;
        parm[p].out = fopen(Numbered_Suffix(oprefix,p,".psl"),"w+");
        if (parm[p].out == NULL)
          { fprintf(stderr,"%s: Cannot open %s.%d.paf for reading & writing\n",
                           Prog_Name,oprefix,p);
            exit (1);
          }
        unlink(Numbered_Suffix(oprefix,p,".psl"));
      }

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

    for (p = 1; p < NTHREADS; p++)
      { fclose(parm[p].gdb1.seqs);
        if (ISTWO)
          fclose(parm[p].gdb2.seqs);
      }
    if (units1 != NULL)
      free(units1);
    if (units2 != NULL)
      free(units2);

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

  exit (0);
}
