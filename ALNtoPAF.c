/*******************************************************************************************
 *
 *  Utility to convert a .1aln alignment file into a PAF formated file.
 *
 *  Author:    Gene Myers
 *  Creation:  Nov 2023
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
#include "DNAsource.h"

#undef DEBUG_THREADS

static char *Usage = " [-mx] [-T<int(8)>] <alignment:path>[.1aln]";

static int  CIGAR_M;   // -m
static int  CIGAR_X;   // -x
static int  CIGAR;     // -m or -x
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

void *gen_paf(void *args)
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
  int        alast;
  int        aoff, boff;
  char      *aseq, *bseq;
  Path      *path;
  Work_Data *work;
  int        blocksum, iid;

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

  if (!oneGotoObject (parm->in, beg))
    { fprintf(stderr,"%s: Can't locate to object %lld in aln file\n",Prog_Name,beg) ;
      exit (1);
    }
  oneReadLine(parm->in);

  for (alast = -1; beg < end; beg++)
    { Read_Aln_Overlap(in,ovl);
      path->tlen  = Read_Aln_Trace(in,(uint8 *) trace);
      path->trace = trace;

      aread = ovl->aread;
      aln->alen = reads1[aread].rlen;
      bread = ovl->bread;
      aln->blen  = reads2[bread].rlen;
      aln->flags = ovl->flags;

      aoff = reads1[aread].fpulse;
      fprintf(out,"%s",AHEADER + reads1[aread].coff);
      fprintf(out,"\t%d",ALEN[aread]);
      fprintf(out,"\t%d",aoff + path->abpos);
      fprintf(out,"\t%d",aoff + path->aepos);

      fprintf(out,"\t%c",COMP(aln->flags)?'-':'+');

      fprintf(out,"\t%s",BHEADER + reads2[bread].coff);
      fprintf(out,"\t%d",BLEN[bread]);
      if (COMP(aln->flags))
        { boff = reads2[bread].fpulse + reads2[bread].rlen;
          fprintf(out,"\t%d",boff - path->bepos);
          fprintf(out,"\t%d",boff - path->bbpos);
        }
      else
        { boff = reads2[bread].fpulse;
          fprintf(out,"\t%d",boff + path->bbpos);
          fprintf(out,"\t%d",boff + path->bepos);
        }

      blocksum = (path->aepos-path->abpos) + (path->bepos-path->bbpos);
      iid      = (blocksum - path->diffs)/2;
     
      fprintf(out,"\t%d",iid);
      fprintf(out,"\t%d",blocksum/2);
      fprintf(out,"\t255");

      fprintf(out,"\tdv:f:%.04f",1.*((path->aepos-path->abpos)-iid)/(path->aepos-path->abpos));
      fprintf(out,"\tdf:i:%d",path->diffs);

      if (CIGAR)
        { int  bmin, bmax;
          char *bact;

          Decompress_TraceTo16(ovl);

          if (aread != alast)
            { if (Load_Read(db1,aread,aseq,0))
                exit (1);
            }

          if (COMP(aln->flags))
            { bmin = (aln->blen-path->bepos);
              bmax = (aln->blen-path->bbpos);
            }
          else
            { bmin = path->bbpos;
              bmax = path->bepos;
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

          if (CIGAR_M)
            { int    k, h, p, x, blen;
              int32 *t = (int32 *) path->trace;
              int    T = path->tlen;
              int    ilen, dlen;

              ilen = dlen = 0;
              fprintf(out,"\tcg:Z:");
              k = path->abpos+1;
              h = path->bbpos+1;
              for (x = 0; x < T; x++)
                { if ((p = t[x]) < 0)
                    { blen = -(p+k);
                      k += blen;
                      h += blen+1;
                      if (dlen > 0)
                        fprintf(out,"%dI",dlen);
                      dlen = 0;
                      if (blen == 0)
                        ilen += 1;
                      else
                        { if (ilen > 0)
                            fprintf(out,"%dD",ilen);
                          fprintf(out,"%dM",blen);
                          ilen = 1;
                        }
                    }
                  else
                    { blen = p-h;
                      k += blen+1;
                      h += blen;
                      if (ilen > 0)
                        fprintf(out,"%dD",ilen);
                      ilen = 0;
                      if (blen == 0)
                        dlen += 1;
                      else
                        { if (dlen > 0)
                            fprintf(out,"%dI",dlen);
                          fprintf(out,"%dM",blen);
                          dlen = 1;
                        }
                    }
                }
              if (dlen > 0)
                fprintf(out,"%dI",dlen);
              if (ilen > 0)
                fprintf(out,"%dD",ilen);
              blen = (path->aepos - k)+1;
              if (blen > 0)
                fprintf(out,"%dM",blen);
            }

          else  //  CIGAR_X
            { int    k, h, p, x, b, blen;
              int32 *t = (int32 *) path->trace;
              int    T = path->tlen;
              int    ilen, dlen;
              int    xlen, elen;
              char  *A, *B;

              A = aln->aseq-1;
              B = aln->bseq-1;
              ilen = dlen = 0;
              fprintf(out,"\tcg:Z:");
              k = path->abpos+1;
              h = path->bbpos+1;
              for (x = 0; x < T; x++)
                { if ((p = t[x]) < 0)
                    { blen = -(p+k);
                      if (dlen > 0)
                        fprintf(out,"%dI",dlen);
                      dlen = 0;
                      if (blen == 0)
                        ilen += 1;
                      else
                        { if (ilen > 0)
                            fprintf(out,"%dD",ilen);
                          elen = xlen = 0;
                          for (b = 0; b < blen; b++, k++, h++)
                            if (A[k] == B[h])
                              { if (xlen > 0)
                                  fprintf(out,"%dX",xlen);
                                xlen = 0;
                                elen += 1;
                              }
                            else
                              { if (elen > 0)
                                  fprintf(out,"%d=",elen);
                                elen = 0;
                                xlen += 1;
                              }
                          if (xlen > 0)
                            fprintf(out,"%dX",xlen);
                          if (elen > 0)
                            fprintf(out,"%d=",elen);
                          ilen = 1;
                        }
                      h += 1;
                    }
                  else
                    { blen = p-h;
                      if (ilen > 0)
                        fprintf(out,"%dD",ilen);
                      ilen = 0;
                      if (blen == 0)
                        dlen += 1;
                      else
                        { if (dlen > 0)
                            fprintf(out,"%dI",dlen);
                          elen = xlen = 0;
                          for (b = 0; b < blen; b++, k++, h++)
                            if (A[k] == B[h])
                              { if (xlen > 0)
                                  fprintf(out,"%dX",xlen);
                                xlen = 0;
                                elen += 1;
                              }
                            else
                              { if (elen > 0)
                                  fprintf(out,"%d=",elen);
                                elen = 0;
                                xlen += 1;
                              }
                          if (xlen > 0)
                            fprintf(out,"%dX",xlen);
                          if (elen > 0)
                            fprintf(out,"%d=",elen);
                          dlen = 1;
                        }
                      k += 1;
                    }
                }
              if (dlen > 0)
                fprintf(out,"%dI",dlen);
              if (ilen > 0)
                fprintf(out,"%dD",ilen);
              blen = (path->aepos - k)+1;
              if (blen > 0)
                { elen = xlen = 0;
                  for (b = 0; b < blen; b++, k++, h++)
                    if (A[k] == B[h])
                      { if (xlen > 0)
                          fprintf(out,"%dX",xlen);
                        xlen = 0;
                        elen += 1;
                      }
                    else
                      { if (elen > 0)
                          fprintf(out,"%d=",elen);
                        elen = 0;
                        xlen += 1;
                      }
                  if (xlen > 0)
                    fprintf(out,"%dX",xlen);
                  if (elen > 0)
                    fprintf(out,"%d=",elen);
                }
            }
        }

      fprintf(out,"\n");
      alast = aread;
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
  char      *db1_name, *db2_name;
  int        hdrs1, hdrs2;
  FILE     **bps1, **bps2;
  OneFile   *input;
  int64      novl;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("ALNtoPAF")

    NTHREADS = 8;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("mx")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    CIGAR_X = flags['x'];
    CIGAR_M = flags['m'];
    CIGAR   = CIGAR_X || CIGAR_M;

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -m: produce Cigar string tag with M's\n");
        fprintf(stderr,"      -x: produce Cigar string tag with X's and ='s\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
        exit (1);
      }

    if (CIGAR_X + CIGAR_M > 1)
      { fprintf(stderr,"%s: Only one of -m, -x, or -t can be set\n",Prog_Name);
        exit (1);
      }

    parm = Malloc(sizeof(Packet)*NTHREADS,"Allocating thread records");
    if (parm == NULL)
      exit (1);
  }

  //  Initiate .1aln file reading and read header information

  { char       *pwd, *root, *cpath, *spath;
    char       *src1_name, *src2_name;
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

    ISTWO = 0;
    type  = get_dna_paths(src1_name,src1_name,&spath,&db1_name,0);
    if (type != IS_GDB)
      bps1 = make_temporary_gdb(spath,db1_name,db1,&hdrs1,NTHREADS);
    else
      { if (Open_DB(db1_name,db1) < 0)
          exit (1);
      }

    free(spath);
    if (src2_name != NULL)
      { type = get_dna_paths(src2_name,src2_name,&spath,&db2_name,0);
        if (type != IS_GDB)
          bps2 = make_temporary_gdb(spath,db2_name,db2,&hdrs2,NTHREADS);
        else
          { if (Open_DB(db2_name,db2) < 0)
              exit (1);
          }
        free(spath);
        ISTWO = 1;
      }
    else
      { db2 = db1;
        bps2 = bps1;
      }
  }

  //  Setup scaffold->contig maps & scaffold name dictionary

  { int   s, r;

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

  { int r;
    char *eptr;
    struct stat state;

    if (fstat(hdrs1,&state) < 0)
      { fprintf(stderr,"%s: Cannot fetch size of %s's header file\n",Prog_Name,argv[1]);
        exit (1);
      }

    AHEADER = Malloc(state.st_size,"Allocating header table");
    if (AHEADER == NULL)
      exit (1);

    if (read(hdrs1,AHEADER,state.st_size) < 0)
      { fprintf(stderr,"%s: Cannot read header file of %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    close(hdrs1);

    for (r = 0; r < db1->nreads; r++)
      if (db1->reads[r].origin == 0)
        { for (eptr = AHEADER + db1->reads[r].coff; *eptr != '\n'; eptr++)
            if (isspace(*eptr))
              break;
          *eptr = '\0';
        }

    if (ISTWO)
      { if (fstat(hdrs2,&state) < 0)
          { fprintf(stderr,"%s: Cannot fetch size of %s's header file\n",Prog_Name,argv[2]);
            exit (1);
          }

        BHEADER = Malloc(state.st_size,"Allocating header table");
        if (BHEADER == NULL)
          exit (1);

        if (read(hdrs2,BHEADER,state.st_size) < 0)
          { fprintf(stderr,"%s: Cannot read header file of %s\n",Prog_Name,argv[2]);
            exit (1);
          }
        close(hdrs2);

        for (r = 0; r < db2->nreads; r++)
          if (db2->reads[r].origin == 0)
            { for (eptr = BHEADER + db2->reads[r].coff; *eptr != '\n'; eptr++)
                if (isspace(*eptr))
                  break;
              *eptr = '\0';
            }
      }
    else
      BHEADER = AHEADER;
  }

  //  Divide .1aln into NTHREADS parts

  { int p;

    for (p = 0; p < NTHREADS ; p++)
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

    oprefix = Strdup(Numbered_Suffix("_paf.",getpid(),"."),"Allocating temp file prefix");

    for (p = 0; p < NTHREADS; p++)
      { parm[p].db1   = *db1;
        parm[p].db2   = *db2;
        if (p > 0)
          { parm[p].db1.bases = bps1[p];
            parm[p].db2.bases = bps2[p];
          }
        parm[p].in  = input + p ;
        parm[p].out = fopen(Numbered_Suffix(oprefix,p,".paf"),"w+");
        if (parm[p].out == NULL)
          { fprintf(stderr,"%s: Cannot open %s.%d.paf for reading & writing\n",
                           Prog_Name,oprefix,p);
            exit (1);
          }
        unlink(Numbered_Suffix(oprefix,p,".paf"));
      }

    //  Launch and then gather threads

#ifdef DEBUG_THREADS
    for (p = 0; p < NTHREADS; p++)
      gen_paf(parm+p);
#else
    for (p = 1; p < NTHREADS; p++)
      pthread_create(threads+p,NULL,gen_paf,parm+p);
    gen_paf(parm);
    for (p = 1; p < NTHREADS; p++)
      pthread_join(threads[p],NULL);
#endif

    for (p = 1; p < NTHREADS; p++)
      { fclose(bps1[p]);
        fclose(bps2[p]);
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
