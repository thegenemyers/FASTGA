/*****************************************************************************************\
*                                                                                         *
*  Genome selection parser & interpreter                                                  *
*                                                                                         *
*  Author:  Gene Myers (with assist from Chenxi Zhou                                      *
*  Date  :  June 2024                                                                     *
*                                                                                         *
\*****************************************************************************************/

#ifndef _GDB_SELECTION

#define _GDB_SELECTION

  //  get_selection_contigs/list take a selection expression 'expr', the GDB to which
  //    the selection is to be applied, and a hash table of all the scaffold names
  //    in the GDB.  The GDB does not need to have its base pair component, only the
  //    scaffold skeleton is required.
  //
  //  get_selection_contigs returns an array of records, one per contig, that indicate whether
  //    the contig is in the selection (order > 0) and if so, the range [beg,end]
  //    that is part of the selection.  If ordered is non-zero, then each range must be
  //    disjoint and order gives the ordinal position of the range in the expression.
  //    The array is allocated by the routine and given to the user to dispose of.
  //
  //  get_selection_list returns an array of records whose length is set in *nlen.
  //    Each triple (s1,c1,p1) or (s2,c2,p2) gives a location in the genome where
  //    the contig and scaffold coordinates are always absolute and the position is
  //    always relative to the contig.  If the type is CONTG_ ro SCAFF_SELECTION then
  //    the selection is the range between the location 1 and 2, where a CONTG selection
  //    is of the contig sequencess in the range, and a SCAFF selection is of the scaffold
  //    sequences in the range (a scaffold sequence is a single sequence of contig sequences
  //    with N's in the gaps).
  //    The type POINT_SELECTION is only returned by interpret_point which is designed for
  //    the ALNview viewer.  If it is set then point 1 and point 2 are a coordinate pair of
  //    locations in two genomes.

typedef struct
  { int order;     //  0 if out, order in selection expression if > 0
    int beg;       //  if -1 then not in selection, otherwise [beg..end] of contig is in.
    int end;
    int orient;    //  +1 if forward, -1 if reverse, 0 if none
  } Contig_Range;

Contig_Range *get_selection_contigs(char *expr, GDB *gdb, Hash_Table *hash, int ordered);

#define CONTG_SELECTION    0x0
#define SCAFF_SELECTION    0x1
#define POINT_SELECTION    0x2

typedef struct
  { int   type;
    int   orient;        //  +1 if forward, -1 if reverse, 0 if none
    int   s1, s2;
    int   c1, c2;
    int64 p1, p2;
  } Selection;

Selection *get_selection_list(char *expr, GDB *gdb, Hash_Table *hash, int *nlen);

  //  The routines below are for the interactive app ALNview to interpret individual
  //    ranges and the focal point within a dot plot:

    //  interpret_range sets the record pointed at by selection according to the range
    //    expression assumed to be in 'expr'.  

int interpret_range(Selection *select, char *expr, GDB *gdb, Hash_Table *hash);

    //  interpret_point sets the record pointed at by selection according to the pair/matrix
    //    expression assumed to be in 'expr'.  

int interpret_point(Selection *select, char *expr, GDB *gdb1, Hash_Table *hash1,
                                                   GDB *gdb2, Hash_Table *hash2);

#endif  // _GDB_SELECTION
