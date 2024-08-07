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

#include "GDB.h"
#include "align.h"
#include "alncode.h"

#undef DEBUG_THREADS

static char *Usage = " [-mx] [-T<int(8)>] <alignment:path>[.1aln]";

static int  CIGAR_M;   // -m
static int  CIGAR_X;   // -x
static int  CIGAR;     // -m or -x
static int  NTHREADS;  // -T
static int  ISTWO;     // one gdb or two?

static int  TSPACE;   // Trace spacing

//  THREAD ROUTINE TO GENERATE A SECTION OF THE DESIRED .PAF FILE

typedef struct
  { int64    beg;
    int64    end;
    GDB      gdb1;
    GDB      gdb2;
    OneFile *in;
    FILE    *out;
  } Packet;

void *gen_paf(void *args)
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
  GDB_CONTIG   *contigs1, *contigs2;
  GDB_SCAFFOLD *scaff1, *scaff2;
  char         *ahead, *bhead;
  int           acontig, bcontig;
  int           ascaff, bscaff;
  int           alast;
  int64         aoff, boff;
  char         *aseq, *bseq;
  Path         *path;
  Work_Data    *work;
  int           blocksum, iid;

  work = New_Work_Data();
  aseq = New_Contig_Buffer(gdb1);
  bseq = New_Contig_Buffer(gdb2);
  if (aseq == NULL || bseq == NULL)
    exit (1);

  contigs1 = gdb1->contigs;
  contigs2 = gdb2->contigs;
  scaff1   = gdb1->scaffolds;
  scaff2   = gdb2->scaffolds;
  ahead    = gdb1->headers;
  bhead    = gdb2->headers;

  aln->path = path = &(ovl->path);
  aln->aseq = aseq;

  tmax  = in->info['T']->given.max;
  trace = (uint16 *) Malloc(2*sizeof(uint16)*tmax,"Allocating trace vector");
  if (trace == NULL)
    exit (1);
  path->trace = (void *) trace;

  //  For each alignment do

  if (!oneGoto(in,'A',beg+1))
    { fprintf(stderr,"%s: Can't locate to object %lld in aln file\n",Prog_Name,beg+1);
      exit (1);
    }
  oneReadLine(in);

  for (alast = -1; beg < end; beg++)
    { Read_Aln_Overlap(in,ovl);
      path->tlen  = Read_Aln_Trace(in,(uint8 *) trace);
      path->trace = trace;

      acontig = ovl->aread;
      aln->alen = contigs1[acontig].clen;
      bcontig = ovl->bread;
      aln->blen  = contigs2[bcontig].clen;
      aln->flags = ovl->flags;
      ascaff = contigs1[acontig].scaf;
      bscaff = contigs2[bcontig].scaf;

      aoff = contigs1[acontig].sbeg;
      fprintf(out,"%s",ahead + scaff1[ascaff].hoff);
      fprintf(out,"\t%lld",scaff1[ascaff].slen);
      fprintf(out,"\t%lld",aoff + path->abpos);
      fprintf(out,"\t%lld",aoff + path->aepos);

      fprintf(out,"\t%c",COMP(aln->flags)?'-':'+');

      fprintf(out,"\t%s",bhead + scaff2[bscaff].hoff);
      fprintf(out,"\t%lld",scaff2[bscaff].slen);

      if (COMP(aln->flags))
        { boff = contigs2[bcontig].sbeg + contigs2[bcontig].clen;
          fprintf(out,"\t%lld",boff - path->bepos);
          fprintf(out,"\t%lld",boff - path->bbpos);
        }
      else
        { boff = contigs2[bcontig].sbeg;
          fprintf(out,"\t%lld",boff + path->bbpos);
          fprintf(out,"\t%lld",boff + path->bepos);
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

          if (acontig != alast)
            Get_Contig(gdb1,acontig,NUMERIC,aseq);

          if (COMP(aln->flags))
            { bmin = (aln->blen-path->bepos);
              bmax = (aln->blen-path->bbpos);
            }
          else
            { bmin = path->bbpos;
              bmax = path->bepos;
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
      alast = acontig;
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

  { char       *pwd, *root, *cpath;
    char       *src1_name, *src2_name;
    char       *spath, *tpath;
    char       *head, *sptr, *eptr;
    int         type, s;
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
      if (CIGAR)
        units1 = Create_GDB(gdb1,spath,type,NTHREADS,NULL);
      else
        Create_GDB(gdb1,spath,type,0,NULL);
    else
      { Read_GDB(gdb1,tpath);
        if (CIGAR && gdb1->seqs == NULL)
          { fprintf(stderr,"%s: GDB %s must have sequence data\n",Prog_Name,tpath);
            exit (1);
          }
      }
    head  = gdb1->headers;
    for (s = 0; s < gdb1->nscaff; s++)
      { sptr = head + gdb1->scaffolds[s].hoff;
        for (eptr = sptr; *eptr != '\0'; eptr++)
          if (isspace(*eptr))
            break;
        *eptr = '\0';
      }
    free(spath);
    free(tpath);

    if (src2_name != NULL)
      { type = Get_GDB_Paths(src2_name,NULL,&spath,&tpath,0);
        if (type != IS_GDB)
          if (CIGAR)
            units2 = Create_GDB(gdb2,spath,type,NTHREADS,NULL);
          else
            Create_GDB(gdb2,spath,type,0,NULL);
        else
          { Read_GDB(gdb2,tpath);
            if (CIGAR && gdb2->seqs == NULL)
              { fprintf(stderr,"%s: GDB %s must have sequence data\n",Prog_Name,tpath);
                exit (1);
              }
          }
        head  = gdb2->headers;
        for (s = 0; s < gdb2->nscaff; s++)
          { sptr = head + gdb2->scaffolds[s].hoff;
            for (eptr = sptr; *eptr != '\0'; eptr++)
              if (isspace(*eptr))
                break;
            *eptr = '\0';
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
      { parm[p].gdb1 = *gdb1;
        parm[p].gdb2 = *gdb2;
        if (p > 0 && CIGAR)
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
        parm[p].in  = input + p;
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

    if (CIGAR)
      { for (p = 1; p < NTHREADS; p++)
          { fclose(parm[p].gdb1.seqs);
            if (ISTWO)
            fclose(parm[p].gdb2.seqs);
          }
        if (units1 != NULL)
          free(units1);
        if (units2 != NULL)
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

  exit (0);
}
