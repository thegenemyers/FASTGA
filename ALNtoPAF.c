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

char OpComp[256] =
  { ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',

    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ',  '=', ' ', ' ',

    ' ', ' ', ' ', ' ', 'I', ' ', ' ', ' ',
    ' ', 'D', ' ', ' ', ' ',  'M', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    'X', ' ', ' ', ' ', ' ', ' ', ' ', ' ',

    ' ', ' ', ' ', ' ', 'i', ' ', ' ', ' ',
    ' ', 'd', ' ', ' ', ' ',  'm', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
     'x', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
  };

typedef struct
  { int64 n, m;
    char *op;
    int  *ln;
  } CigarList;

static void cigar_push(CigarList *cigar, char op, int ln)
{ if (cigar->n >= cigar->m)
    { cigar->m = cigar->m * 1.2 + 100;
      cigar->op = Realloc(cigar->op,sizeof(char)*cigar->m,"Reallocating CIGAR char buffer");
      cigar->ln = Realloc(cigar->ln,sizeof(int)*cigar->m,"Reallocating CIGAR lens buffer");
      if (cigar->op == NULL || cigar->ln == NULL)
        exit (1);
    }
  cigar->op[cigar->n]   = op;
  cigar->ln[cigar->n++] = ln;
}

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
  CigarList _cig, *cig = &_cig;

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
  
  cig->m  = 1<<16;
  cig->op = Malloc(sizeof(char)*cig->m,"Allocating CIGAR char buffer");
  cig->ln = Malloc(sizeof(int)*cig->m,"Allocating CIGAR lens buffer");
  if (cig->op == NULL || cig->ln == NULL)
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

      if (CIGAR)
        { int  bmin, bmax;
          int  del;
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

          cig->n = 0;
          del = 0;

          if (CIGAR_M)
            { int    k, h, p, x, blen;
              int32 *t = (int32 *) path->trace;
              int    T = path->tlen;
              int    ilen, dlen;

              ilen = dlen = 0;
              k = path->abpos+1;
              h = path->bbpos+1;
              for (x = 0; x < T; x++)
                { if ((p = t[x]) < 0)
                    { blen = -(p+k);
                      k += blen;
                      h += blen+1;
                      if (dlen > 0)
                        //fprintf(out,"%dI",dlen);
                        cigar_push(cig,'I',dlen);    
                      dlen = 0;
                      if (blen == 0)
                        ilen += 1;
                      else
                        { if (ilen > 0)
                            { //fprintf(out,"%dD",ilen);
                              cigar_push(cig,'D',ilen);
                              del += ilen;
                            }
                          //fprintf(out,"%dM",blen);
                          cigar_push(cig,'M',blen);
                          ilen = 1;
                        }
                    }
                  else
                    { blen = p-h;
                      k += blen+1;
                      h += blen;
                      if (ilen > 0)
                        { //fprintf(out,"%dD",ilen);
                          cigar_push(cig,'D',ilen);
                          del += ilen;
                        }
                      ilen = 0;
                      if (blen == 0)
                        dlen += 1;
                      else
                        { if (dlen > 0)
                            //fprintf(out,"%dI",dlen);
                            cigar_push(cig,'I',dlen);
                          //fprintf(out,"%dM",blen);
                          cigar_push(cig,'M',blen);
                          dlen = 1;
                        }
                    }
                }
              if (dlen > 0)
                //fprintf(out,"%dI",dlen);
                cigar_push(cig,'I',dlen);
              if (ilen > 0)
                { //fprintf(out,"%dD",ilen);
                  cigar_push(cig,'D',ilen);
                  del += ilen;
                }
              blen = (path->aepos - k)+1;
              if (blen > 0)
                { //fprintf(out,"%dM",blen);
                  cigar_push(cig,'M',blen);
                }
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
              k = path->abpos+1;
              h = path->bbpos+1;
              for (x = 0; x < T; x++)
                { if ((p = t[x]) < 0)
                    { blen = -(p+k);
                      if (dlen > 0)
                        cigar_push(cig,'I',dlen);
                      dlen = 0;
                      if (blen == 0)
                        ilen += 1;
                      else
                        { if (ilen > 0)
                            { cigar_push(cig,'D',ilen);
                              del += ilen;
                            }
                          elen = xlen = 0;
                          for (b = 0; b < blen; b++, k++, h++)
                            if (A[k] == B[h])
                              { if (xlen > 0)
                                  cigar_push(cig,'X',xlen);
                                xlen = 0;
                                elen += 1;
                              }
                            else
                              { if (elen > 0)
                                  cigar_push(cig,'=',elen);
                                elen = 0;
                                xlen += 1;
                              }
                          if (xlen > 0)
                            cigar_push(cig,'X',xlen);
                          if (elen > 0)
                            cigar_push(cig,'=',elen);
                          ilen = 1;
                        }
                      h += 1;
                    }
                  else
                    { blen = p-h;
                      if (ilen > 0)
                        { cigar_push(cig,'D',ilen);
                          del += ilen;
                        }
                      ilen = 0;
                      if (blen == 0)
                        dlen += 1;
                      else
                        { if (dlen > 0)
                            cigar_push(cig,'I',dlen);
                          elen = xlen = 0;
                          for (b = 0; b < blen; b++, k++, h++)
                            if (A[k] == B[h])
                              { if (xlen > 0)
                                  cigar_push(cig,'X',xlen);
                                xlen = 0;
                                elen += 1;
                              }
                            else
                              { if (elen > 0)
                                  cigar_push(cig,'=',elen);
                                elen = 0;
                                xlen += 1;
                              }
                          if (xlen > 0)
                            cigar_push(cig,'X',xlen);
                          if (elen > 0)
                            cigar_push(cig,'=',elen);
                          dlen = 1;
                        }
                      k += 1;
                    }
                }
              if (dlen > 0)
                cigar_push(cig,'I',dlen);
              if (ilen > 0)
                { cigar_push(cig,'D',ilen);
                  del += ilen;
                }
              blen = (path->aepos - k)+1;
              if (blen > 0)
                { elen = xlen = 0;
                  for (b = 0; b < blen; b++, k++, h++)
                    if (A[k] == B[h])
                      { if (xlen > 0)
                          cigar_push(cig,'X',xlen);
                        xlen = 0;
                        elen += 1;
                      }
                    else
                      { if (elen > 0)
                          cigar_push(cig,'=',elen);
                        elen = 0;
                        xlen += 1;
                      }
                  if (xlen > 0)
                    cigar_push(cig,'X',xlen);
                  if (elen > 0)
                    cigar_push(cig,'=',elen);
                }
            }
          
          blocksum = (path->aepos-path->abpos) + del;
          iid      = blocksum - path->diffs;

          fprintf(out,"\t%d",iid);
          fprintf(out,"\t%d",blocksum);
          fprintf(out,"\t255");

          fprintf(out,"\tdv:f:%.04f",1.*((path->aepos-path->abpos)-iid)/(path->aepos-path->abpos));
          fprintf(out,"\tdf:i:%d",path->diffs);

          { int i;

            fprintf(out,"\tcg:Z:");
            if (COMP(aln->flags))
              for (i = cig->n-1; i >= 0; i--)
                fprintf(out,"%d%c",cig->ln[i],OpComp[(int) cig->op[i]]);
            else
              for (i = 0; i < cig->n; i++)
                fprintf(out,"%d%c",cig->ln[i],cig->op[i]);
          }
        }
      
      else
        { blocksum = (path->aepos-path->abpos) + (path->bepos-path->bbpos);
          iid      = (blocksum - path->diffs)/2;

          // the residue matches and alignment block length are approximate
          fprintf(out,"\t%d",iid);
          fprintf(out,"\t%d",blocksum/2);
          fprintf(out,"\t255");

          fprintf(out,"\tdv:f:%.04f",1.*((path->aepos-path->abpos)-iid)/(path->aepos-path->abpos));
          fprintf(out,"\tdf:i:%d",path->diffs);
        }

      fprintf(out,"\n");
      alast = acontig;
    }

  free(cig->op);
  free(cig->ln);
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
    char       *head, *sptr, *eptr;
    int         s;

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".1aln");
    input = open_Aln_Read(Catenate(pwd,"/",root,".1aln"),NTHREADS,
                          &novl,&TSPACE,&src1_name,&src2_name,&cpath);
    if (input == NULL)
      exit (1);
    free(root);
    free(pwd);

    ISTWO = (src2_name != NULL);

    if (CIGAR)
      units1 = Get_GDB(gdb1,src1_name,cpath,NTHREADS);
    else
      Get_GDB(gdb1,src1_name,cpath,0);

    if (ISTWO)
      { if (CIGAR)
          units2 = Get_GDB(gdb2,src2_name,cpath,NTHREADS);
        else
          Get_GDB(gdb2,src2_name,cpath,0);
      }
    else
      { gdb2   = gdb1;
        units2 = units1;
      }

    free(src1_name);
    free(src2_name);
    free(cpath);

    //  Truncate GDB headers to first white-space

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
        if (CIGAR)
          { parm[p].gdb1.seqs = units1[p];
            parm[p].gdb2.seqs = units2[p];
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

    if (CIGAR && NTHREADS > 1)
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

  exit (0);
}
