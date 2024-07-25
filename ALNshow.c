/*******************************************************************************************
 *
 *  Utility for displaying the overlaps in a .1aln file in a variety of ways including
 *    a minimal listing of intervals, and a full out alignment.
 *
 *  Author:    Gene Myers
 *             Modified by Richard Durbin to work on .1aln files
 *  Creation:  July 2013
 *  Last Mod:  March 2023
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
#include <sys/uio.h>
#include <fcntl.h>

#include "GDB.h"
#include "hash.h"
#include "align.h"
#include "select.h"
#include "alncode.h"

static char *Usage[] =
    { "[-arU] [-i<int(4)>] [-w<int(100)>] [-b<int(10)>] ",
      "    <alignments:path>[.1aln] [<selection>|<FILE> [<selection>|<FILE>]]"
    };

int main(int argc, char *argv[])
{ GDB       _gdb1, *gdb1 = &_gdb1; 
  GDB       _gdb2, *gdb2 = &_gdb2; 
  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;

  OneFile *input;
  int64    novl;
  int      tspace;

  int     ALIGN, REFERENCE;
  int     INDENT, WIDTH, BORDER, UPPERCASE;
  int     ISTWO;

  GDB_SCAFFOLD *ascaffs;
  GDB_SCAFFOLD *bscaffs;

  GDB_CONTIG   *acontigs;
  GDB_CONTIG   *bcontigs;

  Hash_Table   *ahash;
  Hash_Table   *bhash;

  Contig_Range *ACHORD;
  Contig_Range *BCHORD;

  int     nascaff, nacontig, amaxlen, actgmax;
  int     nbscaff, nbcontig, bmaxlen, bctgmax;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("ALNshow")

    INDENT    = 4;
    WIDTH     = 100;
    BORDER    = 10;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("arU")
            break;
          case 'i':
            ARG_NON_NEGATIVE(INDENT,"Indent")
            break;
          case 'w':
            ARG_POSITIVE(WIDTH,"Alignment width")
            break;
          case 'b':
            ARG_NON_NEGATIVE(BORDER,"Alignment border")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    ALIGN     = flags['a'];
    REFERENCE = flags['r'];
    UPPERCASE = flags['U'];

    if (argc < 2 || argc > 4)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"     <selection> = <range> [ , <range> ]*\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"     <range> =     <contig>[-<contig>]    |  <contig>_<int>-(<int>|#)\n");
        fprintf(stderr,"             | @[<scaffold>[-<scaffold>]] | @<scaffold>_<int>-(<int>|#)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"        <contig>   = (<int>|#)[.(<int>|#)]\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"        <scaffold> =  <int>|<string>|#\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -a: Show the alignment of each LA with -w columns in each row.\n");
        fprintf(stderr,"      -r: Show the alignment of each LA with -w bp's of A in each row.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -U: Show alignments in upper case.\n");
        fprintf(stderr,"      -i: Indent alignments by -i spaces.\n");
        fprintf(stderr,"      -w: Width of each row of alignment in symbols (-a) or bps (-r).\n");
        fprintf(stderr,"      -b: # of bordering bp.s to show on each side of LA.\n");
        exit (1);
      }
  }

  //  Initiate .las file reading and read header information
  
  { char  *pwd, *root, *cpath;
    char  *src1_name, *src2_name;
    char  *spath, *tpath;
    int    type;
    FILE  *test;

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".1aln");
    input = open_Aln_Read(Catenate(pwd,"/",root,".1aln"),1,
			   &novl,&tspace,&src1_name,&src2_name,&cpath) ;
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
    type  = Get_GDB_Paths(src1_name,NULL,&spath,&tpath,0);
    if (type != IS_GDB)
      if (ALIGN || REFERENCE)
        Create_GDB(gdb1,spath,type,1,NULL);
      else
        Create_GDB(gdb1,spath,type,0,NULL);
    else
      { Read_GDB(gdb1,tpath);
        if ((ALIGN || REFERENCE) && gdb1->seqs == NULL)
          { fprintf(stderr,"%s: GDB %s must have sequence data\n",Prog_Name,tpath);
            exit (1);
          }
      }
    free(spath);
    free(tpath);

    if (src2_name != NULL)
      { type = Get_GDB_Paths(src2_name,NULL,&spath,&tpath,0);
        if (type != IS_GDB)
          if (ALIGN || REFERENCE)
            Create_GDB(gdb2,spath,type,1,NULL);
          else
            Create_GDB(gdb2,spath,type,0,NULL);
        else
          { Read_GDB(gdb2,tpath);
            if ((ALIGN || REFERENCE) && gdb2->seqs == NULL)
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

  //  Set up scaffold->contig maps & scaffold name dictionary

  { int   s;
    char *head, *sptr, *eptr;

    nacontig = gdb1->ncontig;
    nbcontig = gdb2->ncontig;
    nascaff  = gdb1->nscaff;
    nbscaff  = gdb2->nscaff;
    ascaffs  = gdb1->scaffolds;
    bscaffs  = gdb2->scaffolds;
    acontigs = gdb1->contigs;
    bcontigs = gdb2->contigs;

    ahash = New_Hash_Table(nascaff,0);
    head  = gdb1->headers;
    amaxlen = 0;
    actgmax = 0;
    for (s = 0; s < nascaff; s++)
      { sptr = head + ascaffs[s].hoff;
        for (eptr = sptr; *eptr != '\0'; eptr++)
          if (isspace(*eptr))
            break;
        *eptr = '\0';
        if (Hash_Lookup(ahash,sptr) < 0)
          Hash_Add(ahash,sptr);
        else
          { fprintf(stderr,"%s: Duplicate scaffold name: %s\n",Prog_Name,sptr);
            exit (1);
          }
        if (ascaffs[s].slen > amaxlen)
          amaxlen = ascaffs[s].slen;
        if (ascaffs[s].ectg - ascaffs[s].fctg > actgmax)
          actgmax = ascaffs[s].ectg - ascaffs[s].fctg;
      }

    if (ISTWO)
      { bhash = New_Hash_Table(nbscaff,0);
        head  = gdb2->headers;
        bmaxlen = 0;
        bctgmax = 0;
        for (s = 0; s < nbscaff; s++)
          { sptr = head + bscaffs[s].hoff;
            for (eptr = sptr; *eptr != '\0'; eptr++)
              if (isspace(*eptr))
                break;
            *eptr = '\0';
            if (Hash_Lookup(bhash,sptr) < 0)
              Hash_Add(bhash,sptr);
            else
              { fprintf(stderr,"%s: Duplicate scaffold name: %s\n",Prog_Name,sptr);
                exit (1);
              }
            if (bscaffs[s].slen > bmaxlen)
              bmaxlen = bscaffs[s].slen;
            if (bscaffs[s].ectg - bscaffs[s].fctg > bctgmax)
              bctgmax = bscaffs[s].ectg - bscaffs[s].fctg;
          }
      }
    else
      { bhash   = ahash;
        bmaxlen = amaxlen;
        bctgmax = actgmax;
      }
  }

  //  Setup up contig reporting ranges

  { char *aseq, *bseq;

    aseq = NULL;
    bseq = NULL;
    if (argc == 3)
      aseq = argv[2];
    else if (argc == 4)
      { aseq = argv[2];
        bseq = argv[3];
      }

    ACHORD = get_selection_contigs(aseq,gdb1,ahash,0);
    BCHORD = get_selection_contigs(bseq,gdb2,bhash,0);
  }

  //  Read the file and display selected records
  
  { int           j;
    uint16       *trace;
    Work_Data    *work;
    Contig_Range *aptr;
    Contig_Range *bptr;
    int           tmax;
    int64         tps;

    int        aread, bread;
    int        ascaf, bscaf;
    int        aoffs, boffs;
    int        alens, blens, bclen;
    char      *abuffer, *bbuffer, *root;
    int        ar_wide, br_wide;
    int        ai_wide, bi_wide;
    int        ac_wide, bc_wide;
    int        mn_wide, mx_wide;
    int        tp_wide;

    aln->path = &(ovl->path);
    if (ALIGN || REFERENCE)
      { work = New_Work_Data();
        abuffer = New_Contig_Buffer(gdb1);
        bbuffer = New_Contig_Buffer(gdb2);
      }
    else
      { abuffer = NULL;
        bbuffer = NULL;
        work = NULL;
      }

    tmax = input->info['T']->given.max;
    trace = (uint16 *) Malloc(2*sizeof(uint16)*tmax,"Allocating trace vector");
    if (trace == NULL)
      exit (1);
    ovl->path.trace = (void *) trace;

    ar_wide = Number_Digits((int64) nascaff);
    ai_wide = Number_Digits((int64) amaxlen);
    ac_wide = Number_Digits((int64) actgmax+1);

    br_wide = Number_Digits((int64) nbscaff);
    bi_wide = Number_Digits((int64) bmaxlen);
    bc_wide = Number_Digits((int64) bctgmax+1);

    if (gdb1->maxctg < gdb2->maxctg)
      { mn_wide = Number_Digits((int64) gdb1->maxctg);
        mx_wide = bi_wide;
        if (tspace > 0)
          tp_wide = Number_Digits((int64) gdb1->maxctg/tspace+2);
        else
          tp_wide = 0;
      }
    else
      { mn_wide = Number_Digits((int64) gdb2->maxctg);
        mx_wide = ai_wide;
        if (tspace > 0)
          tp_wide = Number_Digits((int64) gdb2->maxctg/tspace+2);
        else
          tp_wide = 0;
      }
    ar_wide += (ar_wide-1)/3;
    br_wide += (br_wide-1)/3;
    ai_wide += (ai_wide-1)/3;
    bi_wide += (bi_wide-1)/3;
    mn_wide += (mn_wide-1)/3;
    tp_wide += (tp_wide-1)/3;

    root  = Root(argv[1],".las");
    printf("\n%s: ",root);
    Print_Number(novl,0,stdout);
    printf(" records\n");
    free(root);

    //  For each record do

    for (j = 0; j < novl; j++)

       //  Read it in

      { Read_Aln_Overlap(input,ovl);
        ovl->path.tlen  = Read_Aln_Trace(input,(uint8 *) trace);
        ovl->path.trace = trace;

        //  Determine if it should be displayed

        aread = ovl->aread;
        if (aread >= nacontig)
          { fprintf(stderr,"%s: A-read is out-of-range of GDB %s\n",Prog_Name,argv[1]);
            exit (1);
          }

        aptr = ACHORD+aread;
        if ( ! aptr->order)
          continue;

        bread = ovl->bread;
        if (bread >= nbcontig)
          { fprintf(stderr,"%s: B-read is out-of-range of GDB %s\n",Prog_Name,argv[1+ISTWO]);
            exit (1);
          }

        bptr = BCHORD+bread;
        if ( ! bptr->order)
          continue;
        if (ovl->path.aepos <= aptr->beg || ovl->path.abpos >= aptr->end)
          continue;
        if (ovl->path.bepos <= bptr->beg || ovl->path.bbpos >= bptr->end)
          continue;

        //  Display it

        ascaf = acontigs[aread].scaf;
        bscaf = bcontigs[bread].scaf;

        aln->alen  = acontigs[aread].clen;
        aln->blen  = bcontigs[bread].clen;
        aoffs      = acontigs[aread].sbeg;
	boffs      = bcontigs[bread].sbeg;
	bclen      = bcontigs[bread].clen;
        aln->flags = ovl->flags;
        tps        = ovl->path.tlen/2;

        alens = ascaffs[ascaf].slen;
        blens = bscaffs[bscaf].slen;

        if (ALIGN || REFERENCE)
          printf("\n");

        Print_Number((int64) ascaf+1,ar_wide+1,stdout);
        printf(".%0*d",ac_wide,(aread - ascaffs[ascaf].fctg)+1);
        printf("  ");
        Print_Number((int64) bscaf+1,br_wide+1,stdout);
        printf(".%0*d",bc_wide,(bread - bscaffs[bscaf].fctg)+1);
        if (COMP(ovl->flags))
          printf(" c");
        else
          printf(" n");
        if (ovl->path.abpos+aoffs == 0)
          printf("   <");
        else
          printf("   [");
        Print_Number((int64) ovl->path.abpos+aoffs,ai_wide,stdout);
        printf("..");
        Print_Number((int64) ovl->path.aepos+aoffs,ai_wide,stdout);
        if (ovl->path.aepos+aoffs == alens)
          printf("> x ");
        else
          printf("] x ");
        if (COMP(ovl->flags))
          { if ((bclen-ovl->path.bbpos)+boffs == blens)
              printf("<");
            else
              printf("[");
            Print_Number((int64) boffs+(bclen-ovl->path.bbpos),bi_wide,stdout);
            printf("..");
            Print_Number((int64) boffs+(bclen-ovl->path.bepos),bi_wide,stdout);
            if ((bclen-ovl->path.bepos)+boffs == 0)
              printf(">");
            else
              printf("]");
          }
        else
          { if (ovl->path.bbpos+boffs == 0)
              printf("<");
            else
              printf("[");
            Print_Number((int64) ovl->path.bbpos+boffs,bi_wide,stdout);
            printf("..");
            Print_Number((int64) ovl->path.bepos+boffs,bi_wide,stdout);
            if (ovl->path.bepos+boffs == blens)
              printf(">");
            else
              printf("]");
          }

        printf("  ~  %5.2f%% ",(200.*ovl->path.diffs) /
               ((ovl->path.aepos - ovl->path.abpos) + (ovl->path.bepos - ovl->path.bbpos)) );
        printf("  (");
        Print_Number(alens,ai_wide,stdout);
        printf(" x ");
        Print_Number(blens,bi_wide,stdout);
        printf(" bps,");
        Print_Number((int64) ovl->path.diffs,mn_wide,stdout);
        printf(" diffs, ");
        Print_Number(tps,tp_wide,stdout);
        printf(" trace pts)\n");

        if (ALIGN || REFERENCE)
          { char *aseq, *bseq;
            int   amin,  amax;
            int   bmin,  bmax;
            int   self;

            Decompress_TraceTo16(ovl);

            self = (ISTWO == 0) && (aread == bread) && !COMP(ovl->flags);

            amin = ovl->path.abpos - BORDER;
            if (amin < 0) amin = 0;
            amax = ovl->path.aepos + BORDER;
            if (amax > aln->alen) amax = aln->alen;
            if (COMP(aln->flags))
              { bmin = (aln->blen-ovl->path.bepos) - BORDER;
                if (bmin < 0) bmin = 0;
                bmax = (aln->blen-ovl->path.bbpos) + BORDER;
                if (bmax > aln->blen) bmax = aln->blen;
              }
            else
              { bmin = ovl->path.bbpos - BORDER;
                if (bmin < 0) bmin = 0;
                bmax = ovl->path.bepos + BORDER;
                if (bmax > aln->blen) bmax = aln->blen;
                if (self)
                  { if (bmin < amin)
                      amin = bmin;
                    if (bmax > amax)
                      amax = bmax;
                  }
              }

            aseq = Get_Contig_Piece(gdb1,aread,amin,amax,NUMERIC,abuffer);
            if (!self)
              bseq = Get_Contig_Piece(gdb2,bread,bmin,bmax,NUMERIC,bbuffer);
            else
              bseq = aseq;

            aln->aseq = aseq - amin;
            if (COMP(aln->flags))
              { Complement_Seq(bseq,bmax-bmin);
                aln->bseq = bseq - (aln->blen - bmax);
              }
            else if (self)
              aln->bseq = aln->aseq;
            else
              aln->bseq = bseq - bmin;

            Compute_Trace_PTS(aln,work,tspace,GREEDIEST);

            Gap_Improver(aln,work);

            { int *trace = aln->path->trace;
              int  tlen  = aln->path->tlen;
              int  i;

              aln->path->abpos += aoffs;
              aln->path->aepos += aoffs;
              aln->alen = alens;
              aln->path->bbpos += boffs;
              aln->path->bepos += boffs;
              aln->blen = blens;

              aln->aseq -= aoffs;
              aln->bseq -= boffs;
              for (i = 0; i < tlen; i++)
                if (trace[i] < 0)
                  trace[i] -= aoffs;
                else
                  trace[i] += boffs;
            }

            if (REFERENCE)
              Print_Reference(stdout,aln,work,INDENT,WIDTH,BORDER,UPPERCASE,mx_wide);
            if (ALIGN)
              Print_Alignment(stdout,aln,work,INDENT,WIDTH,BORDER,UPPERCASE,mx_wide);
          }
      }

    free(trace);
    if (ALIGN)
      { free(bbuffer-1);
        free(abuffer-1);
        Free_Work_Data(work);
      }
  }

  oneFileClose(input);

  free(BCHORD);
  free(ACHORD);


  if (ISTWO)
    Free_Hash_Table(bhash);
  Free_Hash_Table(ahash);

  if (ISTWO)
    Close_GDB(gdb2);
  Close_GDB(gdb1);

  free(Prog_Name);
  free(Command_Line);
  Catenate(NULL,NULL,NULL,NULL);

  exit (0);
}
