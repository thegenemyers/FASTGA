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

static char *Usage = " <src1:db|dam> [ <src2:db|dam> ] <align:las>";

int main(int argc, char *argv[])
{ DAZZ_DB   _db1, *db1 = &_db1;
  DAZZ_DB   _db2, *db2 = &_db2;
  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;

  FILE   *input;
  FILE   *ahdr, *bhdr;
  int     ISTWO;
  int64   novl;
  int     tspace, tbytes, small;

  //  Process options

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("LAStoPAF")

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
    int        aread, bread;
    int        bmin, bmax;
    char      *aseq, *bseq, *bact;
    Path      *path;
    Work_Data *work;

    work = New_Work_Data();
    aseq = New_Read_Buffer(db1);
    bseq = New_Read_Buffer(db2);

    tmax  = 1000;
    trace = (uint16 *) Malloc(sizeof(uint16)*tmax,"Allocating trace vector");
    if (trace == NULL)
      exit (1);

    aln->path = path = &(ovl->path);
    aln->aseq = aseq;
    path->trace = (void *) trace;

    //  For each alignment do

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
            aln->alen = db1->reads[aread].rlen;
          }
        aln->blen  = db2->reads[ovl->bread].rlen;
        aln->flags = ovl->flags;

        bread = ovl->bread;
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

        { char header[MAX_NAME];
          int  blocksum;

          fseeko(ahdr,db1->reads[aread].coff,SEEK_SET);
          fgets(header,MAX_NAME,ahdr);
          header[strlen(header)-1] = '\0';
          printf("%s",header+1);

          printf("\t%d",aln->alen);
          printf("\t%d",path->abpos);
          printf("\t%d",path->aepos);

          printf("\t%c",COMP(aln->flags)?'-':'+');

          fseeko(bhdr,db2->reads[bread].coff,SEEK_SET);
          fgets(header,MAX_NAME,bhdr);
          header[strlen(header)-1] = '\0';
          printf("\t%s",header+1);

          printf("\t%d",aln->blen);
          if (COMP(aln->flags))
            { printf("\t%d",aln->blen-path->bepos);
              printf("\t%d",aln->blen-path->bbpos);
            }
          else
            { printf("\t%d",path->bbpos);
              printf("\t%d",path->bepos);
            }

          blocksum = (path->aepos-path->abpos) + (path->bepos-path->bbpos);

          printf("\t%d",(blocksum - (2*path->diffs-path->tlen))/2);
          printf("\t%d",blocksum/2);
          printf("\t255\n");
        }
      }

    free(trace);
    free(bseq-1);
    free(aseq-1);
  }

  if (ISTWO)
    fclose(bhdr);
  fclose(ahdr);

  Close_DB(db1);
  if (ISTWO)
    Close_DB(db2);

  exit (0);
}
