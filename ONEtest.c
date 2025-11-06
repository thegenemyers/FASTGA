#include <stdlib.h>
#include <stdio.h>
#include <strings.h>

#include "ONEaln.h"

int main(int argc, char *argv[])
{ AlnReader *reader;
  AlnGDB    *gdb1, *gdb2;
  AlnRecord *ar;
  int        self;
  int        i, j, x, flags[128];
  int        nscaf, nctg;
  char       snippet[20];
  char      *cigar, *cstag;
  int       *indel;

  ARG_INIT("ONEtest");

  if (argc != 2)
    exit (1);

  reader = alnOpenReader(argv[1],3,1);
  if (reader == NULL)
    exit (1);

  for (i = 0; i < 3; i++)
    { printf("Tspace = %d\n",alnTraceSpacing(reader+i));
      printf("A-cnt  = %d\n",alnCount(reader+i));
      printf("T-max  = %d\n",alnTraceMax(reader+i));
      printf("T-cnt  = %d\n",alnTraceCount(reader+i));
    }

  for (i = 0; i < 3; i++)
    { printf("Search %d\n",i);
      if (!alnGoto(reader+i,alnCount(reader)-3))
        { printf("Goto fail\n");
          exit (1);
        }
      while ( ! alnNext(reader+i))
        printf("Another %c\n",'A'+i); 
    }

  gdb1 = alnGDB1(reader);
  gdb2 = alnGDB2(reader);

  self = 0;
  if (gdb1 == gdb2)
    { self = 1;
      printf("SELF\n");
    }

  printf("Scafs = %d\n",gdbScaffoldCount(gdb1));
  printf("Conts = %d\n",gdbContigCount(gdb1));
  printf("Gaps  = %d\n",gdbGapCount(gdb1));
  printf("Cmax  = %d\n",gdbContigMax(gdb1));
  printf("Gmax  = %d\n",gdbGapMax(gdb1));

  if (!self)
    { printf("Scafs = %d\n",gdbScaffoldCount(gdb1));
      printf("Conts = %d\n",gdbContigCount(gdb1));
      printf("Gaps  = %d\n",gdbGapCount(gdb1));
      printf("Cmax  = %d\n",gdbContigMax(gdb1));
      printf("Gmax  = %d\n",gdbGapMax(gdb1));
    }

  nscaf = gdbScaffoldCount(gdb1);
  printf("Showing first 5 of %d\n",nscaf);
  for (i = 1; i <= 5; i++)
    { printf("Scaffold %d:%s (%d)\n",i,gdbScaffoldName(gdb1,i),gdbScaffoldLen(gdb1,i));
      if (gdbGapLen(gdb1,i,0) > 0)
        printf("  Gap    0: %d\n",gdbGapLen(gdb1,i,0));
      nctg = gdbScaffoldContigs(gdb1,i);
      for (j = 1; j <= nctg; j++)
        { printf("  Contig %d: %d (@%d)\n",j,gdbContigLen(gdb1,i,j),gdbContigStart(gdb1,i,j));
          x = gdbContigStart(gdb1,i,j);
          gdbScaffoldSeq(gdb1,i,x,x+11,snippet);
          printf("  %s\n",snippet);
          if (j < nctg || gdbGapLen(gdb1,i,j) > 0)
            printf("  Gap    %d: %d\n",j,gdbGapLen(gdb1,i,j));
        }
    }

  alnGoto(reader,1);
  while (!alnEOF(reader))
    { ar = alnAlignment(reader,true);
      printf("\n%d[%d..%d] vs %d[%d..%d]\n",ar->seq1,ar->bpos1,ar->epos1,
                                            ar->seq2,ar->bpos2,ar->epos2);
      printf("\nT-points:");
      for (j = 0; j < ar->tlen; j++)
        printf(" %d",ar->tpoints[j]);
      printf("\n\nT-diffs: ");
      for (j = 0; j < ar->tlen; j++)
        printf(" %d",ar->tdiffs[j]);
      printf("\n");

      cigar = alnCreateCigar(ar,true,false);
      printf("\nCigar A vs B = %s\n",cigar);
      cigar = alnCreateCigar(ar,true,true);
      printf("\nCigar B vs A = %s\n",cigar);
      free(cigar);

      cstag = alnCreateCStag(ar,false,false);
      printf("\nCS-tag A vs B = %s\n",cstag);
      cstag = alnCreateCStag(ar,false,true);
      printf("\nCS-tag B vs A = %s\n",cstag);
      free(cstag);

      indel = alnCreateIndelArray(ar,false);
       printf("\nIndels A vs B =");
      for (j = 0; indel[j] != 0; j++)
        printf(" %d",indel[j]);
      printf("\n");
      indel = alnCreateIndelArray(ar,true);
       printf("\nIndels B vs A =");
      for (j = 0; indel[j] != 0; j++)
        printf(" %d",indel[j]);
      printf("\n");

      printf("\nAlignment A vs B\n");
      alnShowAlignment(ar,stdout,8,100,10,9,false,false);

      printf("\nAlignment B vs A\n");
      alnShowAlignment(ar,stdout,8,100,10,9,false,true);

      alnNext(reader);
    }
    

  alnCloseReader(reader);

  exit (0);
}
