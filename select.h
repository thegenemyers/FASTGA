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

  //  get_selection_contigs/list take a selection expression 'selection', the GDB to which
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
  //    Each record gives a selection either as a src[beg..end] if a SUBSTR type
  //    or objects beg-end if a RANGE type.  The type also indicates if the objects
  //    are scaffolds or contigs.  The array is allocated by the routine and
  //    given to the user to dispose of.

typedef struct
  { int order;     //  0 if out, order in selection expression if > 0
    int beg;       //  if -1 then not in selection, otherwise [beg..end] of contig is in.
    int end;
  } Contig_Range;

Contig_Range *get_selection_contigs(char *selection, GDB *gdb, Hash_Table *hash, int ordered);

#define CONTG_RANGE    0x0
#define CONTG_SUBSTR   0x1
#define SCAFF_RANGE    0x2
#define SCAFF_SUBSTR   0x3

typedef struct
  { int type;     //  One of 4 values above
    int src;      //  if SUBSTR then src[beg..end] where src is a CONTG or SCAFF
    int beg;      //  if RANGE then beg-end of CONTG or SCAFF
    int end;
  } Selection;

Selection *get_selection_list(char *selection, GDB *gdb, Hash_Table *hash, int *nlen);

#endif  // _GDB_SELECTION
