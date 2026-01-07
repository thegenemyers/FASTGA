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
    { "[-anrU] [-i<int(4)>] [-w<int(100)>] [-b<int(10)>] ",
      "    <alignments:path>[.1aln] [<selection>|<FILE> [<selection>|<FILE>]]"
    };

int main(int argc, char *argv[])
{ GDB       _gdb1, *gdb1 = &_gdb1; 
  GDB       _gdb2, *gdb2 = &_gdb2; 
  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;
  Path      *path = &(_ovl.path);

  char    *src1_name, *src2_name;

  OneFile *input;
  int64    novl;
  int      tspace;

  int     ALIGN, REFERENCE, USE_NAMES;
  int     INDENT, WIDTH, BORDER, UPPERCASE;
  int     ISTWO;

  GDB_SCAFFOLD *ascaffs;
  GDB_SCAFFOLD *bscaffs;

  GDB_CONTIG   *acontigs;
  GDB_CONTIG   *bcontigs;

  char         *aheaders;
  char         *bheaders;

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
            ARG_FLAGS("anrU")
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
    USE_NAMES = flags['n'];

    if (argc < 2 || argc > 4)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"  <selection> = <range>[+-] [ , <range>[+-] ]*\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"     <range> = <object/position> [ - <object/position> ]  | @ | .\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"        <object/position> = @ <scaffold> [ . <contig>] [ : <position> ]\n");
        fprintf(stderr,"                          |                . <contig>  [ : <position> ]\n");
        fprintf(stderr,"                          |                                <position>\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"           <scaffold> = # | <int> | <identifier>\n");
        fprintf(stderr,"           <contig>   = # | <int>\n");
        fprintf(stderr,"           <position> = # | <int> [ . <int> ] [kMG]\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -a: Show the alignment of each LA with -w columns in each row.\n");
        fprintf(stderr,"      -r: Show the alignment of each LA with -w bp's of A in each row.\n");
        fprintf(stderr,"      -n: Use scaffold names (-a or -r only)\n");
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

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".1aln");
    input = open_Aln_Read(Catenate(pwd,"/",root,".1aln"),1,&novl,&tspace,
                          &src1_name,&src2_name,&cpath) ;
    if (input == NULL)
      exit (1);
    free(root);
    free(pwd);

    if (ALIGN || REFERENCE)
      { Skip_Skeleton(input);
        Get_GDB(gdb1,src1_name,cpath,1,NULL);
      }
    else if (input->lineType == 'g')
      Read_Skeleton(input,src1_name,gdb1);
    else
      Get_GDB(gdb1,src1_name,cpath,0,NULL);

    ISTWO = 0;
    if (src2_name != NULL)
      { if (ALIGN || REFERENCE)
          { Skip_Skeleton(input);
            Get_GDB(gdb2,src2_name,cpath,1,NULL);
          }
        else if (input->lineType == 'g')
          Read_Skeleton(input,src2_name,gdb2);
        else
          Get_GDB(gdb2,src2_name,cpath,0,NULL);
        ISTWO = 1;
      }
    else
      gdb2 = gdb1;

    free(cpath);
  }

  //  Set up scaffold name dictionaries

  { int   s;
    char *sptr, *eptr;

    nacontig = gdb1->ncontig;
    nbcontig = gdb2->ncontig;
    nascaff  = gdb1->nscaff;
    nbscaff  = gdb2->nscaff;
    ascaffs  = gdb1->scaffolds;
    bscaffs  = gdb2->scaffolds;
    acontigs = gdb1->contigs;
    bcontigs = gdb2->contigs;
    aheaders = gdb1->headers;
    bheaders = gdb2->headers;

    ahash = New_Hash_Table(nascaff,0);
    amaxlen = 0;
    actgmax = 0;
    for (s = 0; s < nascaff; s++)
      { sptr = aheaders + ascaffs[s].hoff;
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
        bmaxlen = 0;
        bctgmax = 0;
        for (s = 0; s < nbscaff; s++)
          { sptr = bheaders + bscaffs[s].hoff;
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
    int        aslen, bslen;
    int        aclen, bclen;
    int        ab, ae;
    int        bb, be;
    int        reverse;
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

    if (ALIGN || REFERENCE)
      { ar_wide = br_wide = ai_wide = bi_wide = ac_wide = bc_wide = mn_wide = tp_wide = 0;

        if (amaxlen > bmaxlen)
          mx_wide = Number_Digits((int64) amaxlen);
        else
          mx_wide = Number_Digits((int64) bmaxlen);
      }
    else
      { mx_wide = 0;

        ar_wide = Number_Digits((int64) nascaff);
        ai_wide = Number_Digits((int64) amaxlen);
        ac_wide = Number_Digits((int64) actgmax+1);

        br_wide = Number_Digits((int64) nbscaff);
        bi_wide = Number_Digits((int64) bmaxlen);
        bc_wide = Number_Digits((int64) bctgmax+1);

        if (gdb1->maxctg < gdb2->maxctg)
          { mn_wide = Number_Digits((int64) gdb1->maxctg);
            if (tspace > 0)
              tp_wide = Number_Digits((int64) gdb1->maxctg/tspace+2);
            else
              tp_wide = 0;
          }
        else
          { mn_wide = Number_Digits((int64) gdb2->maxctg);
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
      }

    root  = Root(argv[1],".1aln");
    printf("\n%s: ",root);
    Print_Number(novl,0,stdout);
    printf(" records\n");
    free(root);

    //  For each record do

    for (j = 0; j < novl; j++)

       //  Read it in

      { Read_Aln_Overlap(input,ovl);
        path->tlen  = Read_Aln_Trace(input,(uint8 *) trace,NULL);
        path->trace = trace;

        //  Determine if it should be displayed

        aread = ovl->aread;
        if (aread >= nacontig)
          { fprintf(stderr,"%s: A-read is out-of-range of GDB %s\n",Prog_Name,src1_name);
            exit (1);
          }

        aptr = ACHORD+aread;
        if ( ! aptr->order)
          continue;

        bread = ovl->bread;
        if (bread >= nbcontig)
          { fprintf(stderr,"%s: B-read is out-of-range of GDB %s\n",Prog_Name,src2_name);
            exit (1);
          }

        bptr = BCHORD+bread;
        if ( ! bptr->order)
          continue;
        if (path->aepos <= aptr->beg || path->abpos >= aptr->end)
          continue;
        if (path->bepos <= bptr->beg || path->bbpos >= bptr->end)
          continue;

        if (bptr->orient != 0)
          { if ((aptr->orient >= 0 && bptr->orient < 0) || (aptr->orient < 0 && bptr->orient > 0))
              { if ( ! COMP(ovl->flags))
                  continue;
              }
           else
              { if (COMP(ovl->flags))
                  continue;
              }
          }
           
        //  Display it

        ascaf = acontigs[aread].scaf;
        bscaf = bcontigs[bread].scaf;

        aln->alen  = acontigs[aread].clen;
        aln->blen  = bcontigs[bread].clen;
        aoffs      = acontigs[aread].sbeg;
	aclen      = acontigs[aread].clen;
	boffs      = bcontigs[bread].sbeg;
	bclen      = bcontigs[bread].clen;
        aln->flags = ovl->flags;
        tps        = ovl->path.tlen/2;

        aslen = ascaffs[ascaf].slen;
        bslen = bscaffs[bscaf].slen;

        reverse = (aptr->orient < 0);

        if (ALIGN || REFERENCE)
          printf("\n");

        if (USE_NAMES)
          printf("%s",aheaders+ascaffs[ascaf].hoff);
        else
          Print_Number((int64) ascaf+1,ar_wide+1,stdout);
        printf(".%0*d%c",ac_wide,(aread - ascaffs[ascaf].fctg)+1,reverse?'c':'n');
        printf("  ");
        if (USE_NAMES)
          printf("%s",bheaders+bscaffs[bscaf].hoff);
        else
          Print_Number((int64) bscaf+1,br_wide+1,stdout);
        printf(".%0*d%c",bc_wide,(bread - bscaffs[bscaf].fctg)+1,
                         ((COMP(ovl->flags) == 0) == reverse)?'c':'n');

        if (reverse)
          { ab = aoffs + path->aepos;
            ae = aoffs + path->abpos;
          }
        else
          { ab = aoffs + path->abpos;
            ae = aoffs + path->aepos;
          }

        if (ab == 0 || ab == aslen)
          printf("   <");
        else
          printf("   [");
        Print_Number((int64) ab,ai_wide,stdout);
        printf("..");
        Print_Number((int64) ae,ai_wide,stdout);
        if (ae == 0 || ae == aslen)
          printf("> x ");
        else
          printf("] x ");

        if (COMP(ovl->flags))
          { bb = boffs+(bclen-path->bbpos);
            be = boffs+(bclen-path->bepos);
          }
        else
          { bb = boffs+path->bbpos;
            be = boffs+path->bepos;
          }
        if (reverse)
          { int x = bb; bb = be; be = x; }

        if (bb == 0 || bb == bslen)
          printf("<");
        else
          printf("[");
        Print_Number((int64) bb,bi_wide,stdout);
        printf("..");
        Print_Number((int64) be,bi_wide,stdout);
        if (be == 0 || be == bslen)
          printf(">");
        else
          printf("]");

        if (ALIGN || REFERENCE)
          { char *aseq, *bseq;
            int   amin,  amax;
            int   bmin,  bmax;
            int   self;

            Decompress_TraceTo16(ovl);

            self = (ISTWO == 0) && (aread == bread) && !COMP(ovl->flags);

            amin = path->abpos - BORDER;
            if (amin < 0) amin = 0;
            amax = path->aepos + BORDER;
            if (amax > aclen) amax = aclen;
            if (COMP(aln->flags))
              { bmin = (bclen-path->bepos) - BORDER;
                if (bmin < 0) bmin = 0;
                bmax = (bclen-path->bbpos) + BORDER;
                if (bmax > bclen) bmax = bclen;
              }
            else
              { bmin = path->bbpos - BORDER;
                if (bmin < 0) bmin = 0;
                bmax = path->bepos + BORDER;
                if (bmax > bclen) bmax = bclen;
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
                aln->bseq = bseq - (bclen - bmax);
              }
            else if (self)
              aln->bseq = aln->aseq;
            else
              aln->bseq = bseq - bmin;

            Compute_Trace_PTS(aln,work,tspace,GREEDIEST,1,-1);

            Gap_Improver(aln,work);

            printf("  ~  %5.2f%% ",(200.*ovl->path.diffs) /
                   ((path->aepos - path->abpos) + (path->bepos - path->bbpos)) );
            printf("  (");
            Print_Number(aslen,ai_wide,stdout);
            printf(" x ");
            Print_Number(bslen,bi_wide,stdout);
            printf(" bps, ");
            Print_Number((int64) path->diffs,0,stdout);
            printf(" diffs, ");
            Print_Number(tps,0,stdout);
            printf(" trace pts)\n");

            if (reverse)
              { int  *trace = path->trace;
                int   tlen  = path->tlen;
                int x, j, h;

                Complement_Seq(aseq,amax-amin);
                if (!self)
                  Complement_Seq(bseq,bmax-bmin);
                x = path->abpos;
                path->abpos = aclen - path->aepos;
                path->aepos = aclen - x;
                x = path->bbpos;
                path->bbpos = bclen - path->bepos;
                path->bepos = bclen - x;
                aln->aseq = aseq - (aclen-amax); 
                if (COMP(aln->flags))
                  aln->bseq = bseq - bmin;
                else
                  aln->bseq = bseq - (bclen-bmax); 
                aclen += 2;
                bclen += 2;
                for (j = 0; j < tlen; j++)
                  { x = trace[j]; 
                    if (x < 0)
                      trace[j] = -(aclen+x);
                    else
                      trace[j] = bclen-x;
                  }
                aclen -= 2;
                bclen -= 2;
                for (j = 0, h = tlen-1; j < h; j++, h--)
                  { x = trace[j];
                    trace[j] = trace[h];
                    trace[h] = x;
                  }
              }

            { int *trace = path->trace;
              int  tlen  = path->tlen;
              int  i;

              path->abpos += aoffs;
              path->aepos += aoffs;
              if (reverse)
                aln->alen = 2*aoffs + aclen;
              else
                aln->alen = 0;
              path->bbpos += boffs;
              path->bepos += boffs;
              if ((COMP(ovl->flags) == 0) == reverse)
                aln->blen = 2*boffs + bclen;
              else
                aln->blen = 0;

              aln->aseq -= aoffs;
              aln->bseq -= boffs;
              for (i = 0; i < tlen; i++)
                if (trace[i] < 0)
                  trace[i] -= aoffs;
                else
                  trace[i] += boffs;
            }

            if (REFERENCE)
              Print_Reference(stdout,aln,work,INDENT,WIDTH,BORDER,UPPERCASE,mx_wide,reverse);
            if (ALIGN)
              Print_Alignment(stdout,aln,work,INDENT,WIDTH,BORDER,UPPERCASE,mx_wide,reverse);
          }
        else
          { printf("  ~  %5.2f%% ",(200.*ovl->path.diffs) /
                   ((ovl->path.aepos - ovl->path.abpos) + (ovl->path.bepos - ovl->path.bbpos)) );
            printf("  (");
            Print_Number(aslen,ai_wide,stdout);
            printf(" x ");
            Print_Number(bslen,bi_wide,stdout);
            printf(" bps,");
            Print_Number((int64) ovl->path.diffs,mn_wide,stdout);
            printf(" diffs, ");
            Print_Number(tps,tp_wide,stdout);
            printf(" trace pts)\n");
         }
      }

    free(trace);
    if (ALIGN || REFERENCE)
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

  free(src1_name);
  free(src2_name);

  free(Prog_Name);
  free(Command_Line);
  Catenate(NULL,NULL,NULL,NULL);

  exit (0);
}
