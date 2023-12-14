/*******************************************************************************************
 *
 *  Utility for displaying the overlaps in a .las file in a variety of ways including
 *    a minimal listing of intervals, a cartoon, and a full out alignment.
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

#include "DB.h"
#include "align.h"

static int BORDER = 0;

static char *Usage = " [-mxt] [-T<int(8)>] <src1:db|dam> [ <src2:db|dam> ] <align:las>";

int main(int argc, char *argv[])
{ DAZZ_DB   _db1, *db1 = &_db1;
  DAZZ_DB   _db2, *db2 = &_db2;
  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;

  int  CIGAR_M;
  int  CIGAR_X;
  int  CIGAR;
  int  TRACE;

  FILE   *input;
  FILE   *ahdr, *bhdr;
  int     ISTWO;
  int64   novl;
  int     tspace, tbytes, small;
  int    *amap, *bmap;
  int    *alen, *blen;

  //  Process options

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("LAStoPAF")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("mxt")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    CIGAR_X = flags['x'];
    CIGAR_M = flags['m'];
    TRACE   = flags['t'];
    CIGAR   = CIGAR_X || CIGAR_M;

    if (argc != 3 && argc != 4)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -m: produce Cigar string tag with M's\n");
        fprintf(stderr,"      -x: produce Cigar string tag with X's and ='s\n");
        fprintf(stderr,"      -t: produce LAS trace and diff list tags\n");
        exit (1);
      }

    if (CIGAR_X + CIGAR_M + TRACE > 1)
      { fprintf(stderr,"%s: Only one of -m, -x, or -t can be set\n",Prog_Name);
        exit (1);
      }
  }

  //  Open trimmed DB or DB pair

  { int  r, s;
    int  status;

    ISTWO  = 0;
    status = Open_DB(argv[1],db1);
    if (status < 0)
      exit (1);
    if (db1->part > 0)
      { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    if (status == 0)
      { fprintf(stderr,"%s: Cannot be called on a .db: %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    if (argc > 3)
      { status = Open_DB(argv[2],db2);
        if (status < 0)
          exit (1);
        if (db2->part > 0)
          { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[2]);
            exit (1);
          }
        if (status == 0)
          { fprintf(stderr,"%s: Cannot be called on a .db: %s\n",Prog_Name,argv[2]);
            exit (1);
          }
        Trim_DB(db2);
      }
    else
      db2 = db1;
    Trim_DB(db1);

    ahdr = fopen(Catenate(db1->path,".hdr","",""),"r");
    if (ahdr == NULL)
      { fprintf(stderr,"%s: Cannot open .hdr file for %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    if (ISTWO)
      { bhdr = fopen(Catenate(db2->path,".hdr",NULL,NULL),"r");
        if (bhdr == NULL)
          { fprintf(stderr,"%s: Cannot open .hdr file for %s\n",Prog_Name,argv[2]);
            exit (1);
          }
      }
    else
      bhdr = ahdr;

    if (ISTWO)
      amap = (int *) Malloc(sizeof(int)*2*(db1->treads+db2->treads),"Allocating scaffold map");
    else
      amap = (int *) Malloc(sizeof(int)*2*db1->treads,"Allocating scaffold map");
    if (amap == NULL)
      exit (1);
    alen = amap + db1->treads;

    s = -1;
    for (r = 0; r < db1->treads; r++)
      { if (db1->reads[r].fpulse == 0)
          s += 1;
        alen[s] = db1->reads[r].fpulse + db1->reads[r].rlen;
        amap[r] = s;
      }
    for (r = 0; r < db1->treads; r++)
      alen[r] = alen[amap[r]];

    if (ISTWO)
      { bmap = alen + db1->treads;
        blen = bmap + db2->treads;

        s = -1;
        for (r = 0; r < db2->treads; r++)
          { if (db2->reads[r].fpulse == 0)
              s += 1;
            blen[s] = db2->reads[r].fpulse + db2->reads[r].rlen;
            bmap[r] = s;
          }
        for (r = 0; r < db1->treads; r++)
          blen[r] = blen[bmap[r]];
      }
    else
      { bmap = amap;
        blen = alen;
      }
  }

  //  Initiate file reading and read (novl, tspace) header
  
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
    if (tspace < 0)
      { fprintf(stderr,"%s: Garbage .las file, trace spacing < 0 !\n",Prog_Name);
        exit (1);
      }
    if (tspace <= TRACE_XOVR && tspace != 0)
      { small  = 1;
        tbytes = sizeof(uint8);
      }
    else
      { small  = 0;
        tbytes = sizeof(uint16);
      }

    free(pwd);
    free(root);
  }

  //  Read the file and display selected records
  
  { int        j;
    uint16    *trace;
    int        tmax;
    DAZZ_READ *reads1, *reads2;
    int        aread, bread;
    int        aoff, boff;
    int        alast;
    char      *aseq, *bseq;
    Path      *path;
    Work_Data *work;
    char       aheader[MAX_NAME], bheader[MAX_NAME];
    int        blocksum, iid;

    work = New_Work_Data();
    aseq = New_Read_Buffer(db1);
    bseq = New_Read_Buffer(db2);

    reads1 = db1->reads;
    reads2 = db2->reads;

    tmax  = 1000;
    trace = (uint16 *) Malloc(sizeof(uint16)*tmax,"Allocating trace vector");
    if (trace == NULL)
      exit (1);

    aln->path = path = &(ovl->path);
    aln->aseq = aseq;
    path->trace = (void *) trace;

    //  For each alignment do

    alast = -1;
    for (j = 0; j < novl; j++)
      { Read_Overlap(input,ovl);
        if (ovl->path.tlen > tmax)
          { tmax = ((int) 1.2*ovl->path.tlen) + 100;
            trace = (uint16 *) Realloc(trace,sizeof(uint16)*tmax,"Allocating trace vector");
            if (trace == NULL)
              exit (1);
          }
        path->trace = (void *) trace;
        Read_Trace(input,ovl,tbytes);

        aread = ovl->aread;
        aln->alen = reads1[aread].rlen;
        bread = ovl->bread;
        aln->blen  = reads2[bread].rlen;
        aln->flags = ovl->flags;

        aoff = reads1[aread].fpulse;
        boff = reads2[bread].fpulse;

        if (aread != alast)
          { fseeko(ahdr,reads1[aread].coff,SEEK_SET);
            fgets(aheader,MAX_NAME,ahdr);
            aheader[strlen(aheader)-1] = '\0';
          }
        printf("%s",aheader+1);

        printf("\t%d",alen[aread]);
        printf("\t%d",aoff + path->abpos);
        printf("\t%d",aoff + path->aepos);

        printf("\t%c",COMP(aln->flags)?'-':'+');

        fseeko(bhdr,reads2[bread].coff,SEEK_SET);
        fgets(bheader,MAX_NAME,bhdr);
        bheader[strlen(bheader)-1] = '\0';
        printf("\t%s",bheader+1);

        printf("\t%d",blen[bread]);
        if (COMP(aln->flags))
          { printf("\t%d",blen[bread] - (boff + path->bepos));
            printf("\t%d",blen[bread] - (boff + path->bbpos));
          }
        else
          { printf("\t%d",boff + path->bbpos);
            printf("\t%d",boff + path->bepos);
          }

        blocksum = (path->aepos-path->abpos) + (path->bepos-path->bbpos);
        iid      = (blocksum - path->diffs)/2;
       
        printf("\t%d",iid);
        printf("\t%d",blocksum/2);
        printf("\t255");

        printf("\tdv:F%.04f",1.*((path->aepos-path->abpos)-iid)/(path->aepos-path->abpos));
        printf("\tdf:I%d",path->diffs);

        if (TRACE)
          { int i;

            if (small)
              Decompress_TraceTo16(ovl);

            printf("\ttz:Z%d",trace[1]);
            for (i = 3; i < path->tlen; i+= 2)
              printf(",%d",trace[i]);
            printf("\ttd:Z%d",trace[0]);
            for (i = 2; i < path->tlen; i+= 2)
              printf(",%d",trace[i]);
          }

        if (CIGAR)
          { int  bmin, bmax;
            char *bact;

            if (small && !TRACE)
              Decompress_TraceTo16(ovl);

            if (aread != alast)
              Load_Read(db1,aread,aseq,0);

            if (COMP(aln->flags))
              { bmin = (aln->blen-path->bepos) - BORDER;
                if (bmin < 0) bmin = 0;
                bmax = (aln->blen-path->bbpos) + BORDER;
                if (bmax > aln->blen) bmax = aln->blen;
              }
            else
              { bmin = path->bbpos - BORDER;
                if (bmin < 0) bmin = 0;
                bmax = path->bepos + BORDER;
                if (bmax > aln->blen) bmax = aln->blen;
              }

            bact  = Load_Subread(db2,bread,bmin,bmax,bseq,0);
            if (COMP(aln->flags))
              { Complement_Seq(bact,bmax-bmin);
                aln->bseq = bact - (aln->blen-bmax);
              }
            else
              aln->bseq = bact - bmin; 

            Compute_Trace_PTS(aln,work,tspace,GREEDIEST);

            if (CIGAR_M)
              { int    k, h, p, x, blen;
                int32 *t = (int32 *) path->trace;
                int    T = path->tlen;
                int    ilen, dlen;

                ilen = dlen = 0;
                printf("\tcg:Z");
                k = path->abpos+1;
                h = path->bbpos+1;
                for (x = 0; x < T; x++)
                  { if ((p = t[x]) < 0)
                      { blen = -(p+k);
                        k += blen;
                        h += blen+1;
                        if (dlen > 0)
                          printf("%dD",dlen);
                        dlen = 0;
                        if (blen == 0)
                          ilen += 1;
                        else
                          { if (ilen > 0)
                              printf("%dI",ilen);
                            printf("%dM",blen);
                            ilen = 1;
                          }
                      }
                    else
                      { blen = p-h;
                        k += blen+1;
                        h += blen;
                        if (ilen > 0)
                          printf("%dI",ilen);
                        ilen = 0;
                        if (blen == 0)
                          dlen += 1;
                        else
                          { if (dlen > 0)
                              printf("%dD",dlen);
                            printf("%dM",blen);
                            dlen = 1;
                          }
                      }
                  }
                if (dlen > 0)
                  printf("%dD",dlen);
                if (ilen > 0)
                  printf("%dI",ilen);
                blen = (path->aepos - k)+1;
                if (blen > 0)
                  printf("%dM",blen);
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
                printf("\tcg:Z");
                k = path->abpos+1;
                h = path->bbpos+1;
                for (x = 0; x < T; x++)
                  { if ((p = t[x]) < 0)
                      { blen = -(p+k);
                        if (dlen > 0)
                          printf("%dD",dlen);
                        dlen = 0;
                        if (blen == 0)
                          ilen += 1;
                        else
                          { if (ilen > 0)
                              printf("%dI",ilen);
                            elen = xlen = 0;
                            for (b = 0; b < blen; b++, k++, h++)
                              if (A[k] == B[h])
                                { if (xlen > 0)
                                    printf("%dX",xlen);
                                  xlen = 0;
                                  elen += 1;
                                }
                              else
                                { if (elen > 0)
                                    printf("%d=",elen);
                                  elen = 0;
                                  xlen += 1;
                                }
                            if (xlen > 0)
                              printf("%dX",xlen);
                            if (elen > 0)
                              printf("%d=",elen);
                            ilen = 1;
                          }
                        h += 1;
                      }
                    else
                      { blen = p-h;
                        if (ilen > 0)
                          printf("%dI",ilen);
                        ilen = 0;
                        if (blen == 0)
                          dlen += 1;
                        else
                          { if (dlen > 0)
                              printf("%dD",dlen);
                            elen = xlen = 0;
                            for (b = 0; b < blen; b++, k++, h++)
                              if (A[k] == B[h])
                                { if (xlen > 0)
                                    printf("%dX",xlen);
                                  xlen = 0;
                                  elen += 1;
                                }
                              else
                                { if (elen > 0)
                                    printf("%d=",elen);
                                  elen = 0;
                                  xlen += 1;
                                }
                            if (xlen > 0)
                              printf("%dX",xlen);
                            if (elen > 0)
                              printf("%d=",elen);
                            dlen = 1;
                          }
                        k += 1;
                      }
                  }
                if (dlen > 0)
                  printf("%dD",dlen);
                if (ilen > 0)
                  printf("%dI",ilen);
                blen = (path->aepos - k)+1;
                if (blen > 0)
                  printf("%dM",blen);
              }
          }

        printf("\n");
        alast = aread;
      }

    free(trace);
    free(bseq-1);
    free(aseq-1);
  }

  if (ISTWO)
    fclose(bhdr);
  fclose(ahdr);

  free(amap);
  Close_DB(db1);
  if (ISTWO)
    Close_DB(db2);

  exit (0);
}
