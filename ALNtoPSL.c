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

#include "DB.h"
#include "align.h"
#include "alncode.h"

#undef DEBUG_THREADS

static char *Usage = " [-T<int(8)>} <alignment:path>[.1aln]";

static int  NTHREADS;  // -T
static int  ISTWO;     // one gdb or two?

static int TSPACE;   // Trace spacing

static int  AMAX,  BMAX;   //  Max sequence length in A and B
static int *AMAP, *BMAP;   //  Map contig to its scaffold
static int *ALEN, *BLEN;   //  Map contig to its scaffold length

static char *AHEADER, *BHEADER;  //  Tables of all headers

//  THREAD ROUTINE TO GENERATE A SECTION OF THE DESIRED .PAF FILE

typedef struct
  { int64    beg;
    int64    end;
    DAZZ_DB  db1;
    DAZZ_DB  db2;
    OneFile *in;
    FILE    *out;
  } Packet;

void *gen_psl(void *args)
{ Packet *parm = (Packet *) args;

  int64    beg = parm->beg;
  int64    end = parm->end;
  DAZZ_DB *db1 = &(parm->db1);
  DAZZ_DB *db2 = &(parm->db2);
  OneFile *in  = parm->in;
  FILE    *out = parm->out;

  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;

  uint16    *trace;
  int        tmax;
  DAZZ_READ *reads1, *reads2;
  int        aread, bread;
  int        bmin, bmax;
  char      *aseq, *bseq, *bact;
  int        aoff, boff;
  Path      *path;
  Work_Data *work;

  work = New_Work_Data();
  aseq = New_Read_Buffer(db1);
  bseq = New_Read_Buffer(db2);
  if (aseq == NULL || bseq == NULL)
    exit (1);

  reads1 = db1->reads;
  reads2 = db2->reads;

  aln->path = path = &(ovl->path);
  aln->aseq = aseq;

  tmax  = in->info['T']->given.max;
  trace = (uint16 *) Malloc(2*sizeof(uint16)*tmax,"Allocating trace vector");
  if (trace == NULL)
    exit (1);
  path->trace = (void *) trace;

  //  For each alignment do

  if (!oneGotoObject(parm->in,beg))
    { fprintf(stderr,"%s: Can't locate to object %lld in aln file\n",Prog_Name,beg);
      exit (1);
    }
  oneReadLine(parm->in);

  aoff = 0;
  for (aread = -1; beg < end; beg++)
    { Read_Aln_Overlap(in,ovl);
      path->tlen = Read_Aln_Trace(in,(uint8 *) trace);

      Decompress_TraceTo16(ovl);

      if (aread != ovl->aread)
        { aread = ovl->aread;
          if (Load_Read(db1,aread,aseq,0))
            exit (1);
          aln->alen = reads1[aread].rlen;
          aoff      = reads1[aread].fpulse;
        }
      bread = ovl->bread;
      aln->blen  = reads2[bread].rlen;
      boff       = reads2[bread].fpulse;
      aln->flags = ovl->flags;

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

      bact = Load_Subread(db2,bread,bmin,bmax,bseq,0);
      if (COMP(aln->flags))
        { Complement_Seq(bact,bmax-bmin);
          aln->bseq = bact - (aln->blen-bmax);
        }
      else
        aln->bseq = bact - bmin; 

      Compute_Trace_PTS(aln,work,TSPACE,GREEDIEST);
      if (Gap_Improver(aln,work))
        exit (1);

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
        fprintf(out,"\t%s\t%d\t%d\t%d",
                    AHEADER+reads1[aread].coff,ALEN[aread],aoff+path->abpos,aoff+path->aepos);
        fprintf(out,"\t%s\t%d\t%d\t%d",
                    BHEADER+reads2[bread].coff,BLEN[bread],boff+path->bbpos,boff+path->bepos);

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
  DAZZ_DB   _db1, *db1 = &_db1;
  DAZZ_DB   _db2, *db2 = &_db2;
  char      *db1_name;
  char      *db2_name;
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
    FILE       *test;

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".1aln");
    input = open_Aln_Read(Catenate(pwd,"/",root,".1aln"), NTHREADS,
                          &novl, &TSPACE, &db1_name, &db2_name, &cpath) ;
    if (input == NULL)
      exit (1);
    free(root);
    free(pwd);

    test = fopen(db1_name,"r");
    if (test == NULL)
      { if (*db1_name != '/')
          test = fopen(Catenate(cpath,db1_name,"",""),"r");
        if (test == NULL)
          { fprintf(stderr,"%s: Could not find .gdb %s\n",Prog_Name,db1_name);
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
          { if (*db2_name != '/')
              test = fopen(Catenate(cpath,db2_name,"",""),"r");
            if (test == NULL)
              { fprintf(stderr,"%s: Could not find .gdb %s\n",Prog_Name,db2_name);
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

  //  Open DB or DB pair

  { int   s, r, gdb;

    ISTWO  = 0;
    gdb    = Open_DB(db1_name,db1);
    if (gdb < 0)
      exit (1);

    if (db2_name != NULL)
      { gdb = Open_DB(db2_name,db2);
        if (gdb < 0)
          exit (1);
        ISTWO = 1;
      }
    else
      db2  = db1;

    //  Build contig to scaffold maps in global vars

    AMAX = db1->maxlen;
    BMAX = db2->maxlen;

    if (ISTWO)
      AMAP = (int *) Malloc(sizeof(int)*2*(db1->treads+db2->treads),"Allocating scaffold map");
    else
      AMAP = (int *) Malloc(sizeof(int)*2*db1->treads,"Allocating scaffold map");
    if (AMAP == NULL)
      exit (1);
    ALEN = AMAP + db1->treads;

    s = -1;
    for (r = 0; r < db1->treads; r++)
      { if (db1->reads[r].origin == 0)
          s += 1;
        ALEN[s] = db1->reads[r].fpulse + db1->reads[r].rlen;
        AMAP[r] = s;
      }
    for (r = db1->treads-1; r >= 0; r--)
      ALEN[r] = ALEN[AMAP[r]];

    if (ISTWO)
      { BMAP = ALEN + db1->treads;
        BLEN = BMAP + db2->treads;

        s = -1;
        for (r = 0; r < db2->treads; r++)
          { if (db2->reads[r].origin == 0)
              s += 1;
            BLEN[s] = db2->reads[r].fpulse + db2->reads[r].rlen;
            BMAP[r] = s;
          }
        for (r = db2->treads-1; r >= 0; r--)
          BLEN[r] = BLEN[BMAP[r]];
      }
    else
      { BMAP = AMAP;
        BLEN = ALEN;
      }
  }

  //  Preload all scaffold headers

  { int r, hdrs;
    struct stat state;

    hdrs = open(Catenate(db1->path,".hdr","",""),O_RDONLY);
    if (hdrs < 0)
      { fprintf(stderr,"%s: Cannot open header file of %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    if (fstat(hdrs,&state) < 0)
      { fprintf(stderr,"%s: Cannot fetch size of %s's header file\n",Prog_Name,argv[1]);
        exit (1);
      }

    AHEADER = Malloc(state.st_size,"Allocating header table");
    if (AHEADER == NULL)
      exit (1);

    if (read(hdrs,AHEADER,state.st_size) < 0)
      { fprintf(stderr,"%s: Cannot read header file of %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    AHEADER[state.st_size-1] = '\0';
    close(hdrs);

    for (r = 1; r < db1->nreads; r++)
      if (db1->reads[r].origin == 0)
        AHEADER[db1->reads[r].coff-1] = '\0';

    if (ISTWO)
      { hdrs = open(Catenate(db2->path,".hdr","",""),O_RDONLY);
        if (hdrs < 0)
          { fprintf(stderr,"%s: Cannot open header file of %s\n",Prog_Name,argv[1]);
            exit (1);
          }
        if (fstat(hdrs,&state) < 0)
          { fprintf(stderr,"%s: Cannot fetch size of %s's header file\n",Prog_Name,argv[1]);
            exit (1);
          }

        BHEADER = Malloc(state.st_size,"Allocating header table");
        if (BHEADER == NULL)
          exit (1);

        if (read(hdrs,BHEADER,state.st_size) < 0)
          { fprintf(stderr,"%s: Cannot read header file of %s\n",Prog_Name,argv[1]);
            exit (1);
          }
        BHEADER[state.st_size-1] = '\0';
        close(hdrs);

        for (r = 1; r < db2->nreads; r++)
          if (db2->reads[r].origin == 0)
            BHEADER[db2->reads[r].coff-1] = '\0';
      }
    else
      BHEADER = AHEADER;
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
      { parm[p].db1   = *db1;
        parm[p].db2   = *db2;
        if (p > 0)
          { parm[p].db1.bases = fopen(Catenate(db1->path,"","",".bps"),"r");
            if (parm[p].db1.bases == NULL)
              { fprintf(stderr,"%s: Cannot open another copy of DB\n",Prog_Name);
                exit (1);
              }
            parm[p].db2.bases = fopen(Catenate(db2->path,"","",".bps"),"r");
            if (parm[p].db2.bases == NULL)
              { fprintf(stderr,"%s: Cannot open another copy of DB\n",Prog_Name);
                exit (1);
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
      { fclose(parm[p].db1.bases);
        fclose(parm[p].db2.bases);
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

  free(AMAP);
  Close_DB(db1);
  if (ISTWO)
    Close_DB(db2);

  oneFileClose(input);

  exit (0);
}
