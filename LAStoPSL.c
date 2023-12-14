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
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "align.h"

static int BORDER = 0;

static char *Usage = " <src1:db|dam> [ <src2:db|dam> ] <align:las>";

int main(int argc, char *argv[])
{ DAZZ_DB   _db1, *db1 = &_db1;
  DAZZ_DB   _db2, *db2 = &_db2;
  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;

  FILE   *input;
  int     ISTWO;
  int64   novl;
  int     tspace, tbytes, small;
  FILE   *ahdr, *bhdr;
  int    *amap, *bmap;
  int    *alen, *blen;

  //  Process options

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("LAStoPSL")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 3 && argc != 4)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  //  Open trimmed DB or DB pair

  { int status;
    int s, r;

    ISTWO  = (argc > 3);
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
    if (ISTWO)
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
      { bhdr = fopen(Catenate(db2->path,".hdr","",""),"r");
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
    int        bmin, bmax;
    char      *aseq, *bseq, *bact;
    int        aoff, boff;
    Path      *path;
    Work_Data *work;
    char       aheader[MAX_NAME], bheader[MAX_NAME];

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

    aread = -1;
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

        if (small)
          Decompress_TraceTo16(ovl);

        if (aread != ovl->aread)
          { aread = ovl->aread;
            Load_Read(db1,aread,aseq,0);
            aln->alen = reads1[aread].rlen;
            aoff      = reads1[aread].fpulse;
            fseeko(ahdr,reads1[aread].coff,SEEK_SET);
            fgets(aheader,MAX_NAME,ahdr);
            aheader[strlen(aheader)-1] = '\0';
          }
        bread = ovl->bread;
        aln->blen  = reads2[bread].rlen;
        boff       = reads2[bread].fpulse;
        fseeko(bhdr,reads2[bread].coff,SEEK_SET);
        fgets(bheader,MAX_NAME,bhdr);
        bheader[strlen(bheader)-1] = '\0';
        aln->flags = ovl->flags;

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

        bact = Load_Subread(db2,bread,bmin,bmax,bseq,0);
        if (COMP(aln->flags))
          { Complement_Seq(bact,bmax-bmin);
            aln->bseq = bact - (aln->blen-bmax);
          }
        else
          aln->bseq = bact - bmin; 

        Compute_Trace_PTS(aln,work,tspace,GREEDIEST);

        // Print_Alignment(stdout,aln,work,4,100,BORDER,0,9);

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

          printf("%d %d 0 0 %d %d %d %d +%c",X,S,IB,I,DB,D,COMP(ovl->flags)?'-':'+');
          printf(" %s %d %d %d",aheader+1,alen[aread],aoff+path->abpos,aoff+path->aepos);
          printf(" %s %d %d %d",bheader+1,blen[bread],boff+path->bbpos,boff+path->bepos);

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
          printf(" %d ",bcnt);

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
                printf("%d,",bmat);
            }
          bmat = (path->aepos - i)+1;
          if (bmat > 0)
            printf("%d,",bmat);
          printf(" ");

          i = path->abpos+1;
          j = path->bbpos+1;
          for (x = 0; x < T; x++)
            { if ((p = t[x]) < 0)
                { bmat = -(p+i);
                  if (bmat > 0)
                    printf("%d,",i);
                  i += bmat;
                  j += bmat+1;
                }
              else
                { bmat = p-j;
                  if (bmat > 0)
                    printf("%d,",i);
                  i += bmat+1;
                  j += bmat;
                }
            }
          bmat = (path->aepos - i)+1;
          if (bmat > 0)
            printf("%d,",i);
          printf(" ");

          i = path->abpos+1;
          j = path->bbpos+1;
          for (x = 0; x < T; x++)
            { if ((p = t[x]) < 0)
                { bmat = -(p+i);
                  if (bmat > 0)
                    printf("%d,",j);
                  i += bmat;
                  j += bmat+1;
                }
              else
                { bmat = p-j;
                  if (bmat > 0)
                    printf("%d,",j);
                  i += bmat+1;
                  j += bmat;
                }
            }
          bmat = (path->aepos - i)+1;
          if (bmat > 0)
            printf("%d,",j);
          printf("\n");
        }
      }

    free(trace);
    free(bseq-1);
    free(aseq-1);
  }

  free(amap);
  Close_DB(db1);
  if (ISTWO)
    Close_DB(db2);

  exit (0);
}
