/*******************************************************************************************
 *
 *  Interface to read and process the contents of a .1aln file.
 *
 *  Author :  Gene Myers
 *  Date   :  October 2025
 *
 ********************************************************************************************/

#ifndef _1ALN_DEFS

#define _1ALN_DEFS

#include <stdio.h>
#include "gene_core.h"
#include "ONElib.h"
#include "GDB.h"
#include "align.h"
#include "alncode.h"

  //  OVERVIEW:

  //  This package is intended to give one a simple interface to read and manipulate the
  //    information in a .1aln file of alignment records produced by FastGA or FasTAN.
  //    A simple use pattern for reading all the alignment records in a file is as follows:

  /*        AlnReader *reader = alnOpenReader(filename,1,true);
            if (reader == NULL)
              printf("Error: %s\n",alnError());
            while (!alnEOF(reader))
              { AlnRecord *align = alnAlignment(reader,true);
                if (align == NULL)
                  printf("Error: %s\n",alnError());

                // do stuff with current alignment record align, e.g.

                alnShowAlignment(align,stdout,8,100,10,5,false,false);

                alnNext(reader);
              }
            alnCloseReader(reader);
  */


  //  ERROR HANDLING:

  //  Two macro variables CHECK_ARGS and HALT_ON_ERROR found at the top of the ONEaln.c file
  //    determine which errors are detected and how they are handled according to whether
  //    the variable is defined or undefined.  In the copy you downloaded they are both
  //    defined, change them before compilation to get the desired behavior.

  //  When CHECK_ARGS is defined, all arguments are checked, e.g. the range of an index, or
  //    the state of a reader.  This creates overhead for simple routines like GetContigLen
  //    so if your application is well debugged and ready for release, we suggest you undef
  //    the variable, so that only errors not under you control, e.g. when opening or reading a
  //    file, are caught.  In the error listing at the end of this file, those routines whose
  //    return value must be checked regardless of the setting of CHECK_ARGS are noted.

  //  When HALT_ON_ERROR is undefined, the detection of an error results in an error message
  //    being placed in a special buffer you can access with alnError(), and the routine in
  //    question returns with a documented error value.  But if your program is not interactive
  //    you can for convenience define HALT_ON_ERROR in which case the detection of an error
  //    results in the error message being written to stderr and your program halting with
  //    exit value 1.  Thus you do not need to check any of the return values of the routines in
  //    the package.

  //  When routines return with an error value, you can access an error message by calling
  //    alnError.  The returned pointer is to a globally shared message string, so its
  //    contents are determined by the last routine to detect an error and write into it.
  //    A list of all the error messages emited by this library are found at the end of this
  //    header file.

char *alnError();

  //  .1ALN FILE READER

typedef void *AlnReader;

  //  Open the .1aln file 'name' for reading with nthreads threads.  NULL is returned
  //    if there is an error, otherwise the return pointer points to the first element
  //    in an array of nthreads AlnReader's.  This first element is called the master.
  //    If you want to have CIGAR strings or print Alignments or otherwise need access
  //    to the sequence of the source genomes, you must set see_seq to true.  Otherwise
  //    set it to false for better memory efficiency.

AlnReader *alnOpenReader(char *name, int nthreads, bool see_seq);

  //  Each reader conceptually has a "cursor" pointing at an alignment record that can
  //    be advanced one record at a time with alnNext or jumped to the idx'th record
  //    in the file with alnGoto.  To be crystal clear, alnGoto(..,1) takes you to the
  //    first alignment record, i.e. indexing begins at 1.  alnNext returns true if and
  //    only if the reader has been advanced to the end of the file.  AlnGoto returns
  //    true if and only if it was successful, otherwise an error has occurred.

bool alnNext(AlnReader *reader);

bool alnGoto(AlnReader *reader, int idx);

  //  Return true iff at the end of file

bool alnEOF(AlnReader *reader);

  //  Close all readers associated with the supplied master reader.

void alnCloseReader(AlnReader *reader);

  //  Return (a) the total number of alignments in the file, (b) the maximum number of
  //    trace intervals in any alignment record, (c) the total number of trace intervals
  //    in the file, and (d) the trace point spacing, respectively.

int  alnCount       (AlnReader *reader);
int  alnTraceMax    (AlnReader *reader);
int  alnTraceCount  (AlnReader *reader);
int  alnTraceSpacing(AlnReader *reader);

  //  GENOME DATA BASE: GDB

typedef void AlnGDB;

  //  From the .1aln reader get the first and second genome data bases over which
  //    the alignments were found.  The two pointers are equal if there was only one genome,
  //    i.e. a self-comparison.  A genome data base gives one the structure of the genome
  //    as a collections of scaffolded contigs separated by gaps, as well as access to
  //    the contig's sequence if 'see_seq' was true in the call that created the reader.

AlnGDB *alnGDB1(AlnReader *reader);
AlnGDB *alnGDB2(AlnReader *reader);

  // The Count... routines return the # of (a) scaffolds, (b) contigs, and (c) gaps in a GDB.
  //   The Max... routines return the maximum # of contigs/gaps in any scaffold of a GDB.

int gdbScaffoldCount(AlnGDB *g);
int gdbContigCount  (AlnGDB *g);
int gdbGapCount     (AlnGDB *g);
int gdbContigMax    (AlnGDB *g);
int gdbGapMax       (AlnGDB *g);

  // Scaffolds are numbered starting at 1, and the first contig of a scaffold is indexed by 1.
  //   In rare, somewhat abnormal cases there can be a run of N's at the start or end of
  //   of a scaffold, the length of these are obtained by index 0 and ContigsInScaffold(g,s).
  //   The routines return -1 if an error occurs.

int gdbScaffoldLen    (AlnGDB *g, int s);        // Length of the s'th scaffold
int gdbScaffoldContigs(AlnGDB *g, int s);        // Number of contigs in the s'th scaffold
int gdbGapLen         (AlnGDB *g, int s, int p); // Length of the p'th gap of the s'th scaffold
int gdbContigLen      (AlnGDB *g, int s, int c); // Length of the c'th contig of the s'th scaffold
int gdbContigStart    (AlnGDB *g, int s, int c); // Start position of the c'th contig of the
                                                 //   s'th scaffold (see getScaffoldSeq)

  // Get the scaffold name of the s'th scaffold.  The string returned is in a buffer internal
  //   to the GDB and is overwritten each time this routine is called.  You should copy it to
  //   memory you control if you wish it to persist beyond the last call.  NULL is returned
  //   if an error occurs.

char *gdbScaffoldName(AlnGDB *g, int s);

  // Get the sequence spanning interval [beg,end] of the s'th scaffold of g.  If buffer = NULL
  //   then the routine allocates an array of (end-beg)+1 bytes, places the requested sequence
  //   there, and returns a pointer to it.  In this case the user must subsequently free this
  //   string.  Otherwise, if the user supplies a non-NULL buffer pointer, then it must be of
  //   length not less than (end-beg)+1, the requested sequence is placed there, and the buffer
  //   pointer is returned.  NULL is returned on an error which includes not asking for sequence
  //   access when opening the reader and requesting an interval not wholly within a contig.

  // Sequences positions are *between* base pairs, e.g. position 0 is before the first bp and
  //   position 1 is between the first and second bp.  In this way the subsequence between
  //   an interval [a,b] is of length b-a and consists of the a+1'st through b'th bp of the
  //   selected sequence.

char *gdbScaffoldSeq(AlnGDB *g, int s, int beg, int end, char *buffer);

  //  ALIGNMENT RECORDS

typedef struct
  { int  seq1,  seq2;    //  seq1[bpos1..epos1] aligns to seq2[bpos2..epos2]
    int  bpos1, epos1;   //  with diffs differences (indels & substitutions).
    int  bpos2, epos2;   //  
    int  diffs;          //  

    int  tlen;           //  tpoints[0..tlen) contains the trace point intervals in bseq
    int *tpoints;        //  tdiffs[0..tlen) is the # of diffs in each tp interval
    int *tdiffs;
  } AlnRecord;

  // Load the alignment record at the reader's current cursor/position.  The return pointer
  //   points at a pre-allocated buffer internal to the reader that is reused/overwritten
  //   each time alnAlignment is called.  This includes the memory for the tpoints and tdiffs
  //   arrays of the trace.  You should copy it and the two arrays just mentioned if you
  //   wish them to persist beyond the last call.  If see_trace is true then the trace point
  //   arrays are sought and loaded so one can create cigar strings, etc.  NULL is returned
  //   on error.

AlnRecord *alnAlignment(AlnReader *reader, bool see_trace);

  // Create a CIGAR string for the the given ALN record.  This is only possible if see_seq
  //   was true when the alignment reader was opened, and see_trace was true when the alignment
  //   record was fetched.  For this routine, the returned string must be freed by the user.
  //   NULL is returned on an error.

  // If show_x is true then 'X' and '=' are used to model stretches without indels, otherwise
  //   'M' is used.  When reversed is false, the CIGAR string is one that transforms the first
  //   sequence into the second where the first string is in the forward orientation.  If
  //   reversed is true, then the CIGAR string is one that transforms the second sequence
  //   into the first where the second sequence is in the forward orientation, i.e. the roles
  //   of the first and second sequence are reversed.

char *alnCreateCigar(AlnRecord *align, bool show_x, bool reversed);

  // Create a CStag string for the the given ALN record.  This is only possible if see_seq
  //   was true when the alignment reader was opened, and see_trace was true when the alignment
  //   record was fetched.  For this routine, the returned string must be freed by the user.
  //   NULL is returned on an error.

  // If short_form is false then the normal CStag is created where the sequence of a matching
  //   segment is given, otherwise just the length of a matching segment is given.
  //   When reversed is false, the CS tag is one that transforms the first sequence
  //   into the second where the first string is in the forward orientation.  If reversed
  //   is false, then the CS tag is one that transforms the second sequence into the
  //   first where the second sequence is in the forward orientation, i.e. the roles of
  //   the first and second sequence are reversed.

char *alnCreateCStag(AlnRecord *align, bool short_form, bool reversed);

  // Create an "indel array" for the the given ALN record.  This is only possible if see_seq
  //   was true when the alignment reader was opened, and see_trace was true when the alignment
  //   record was fetched.  For this routine, the returned integer is in a preallocated buffer
  //   of the reader, so you must make a copy if you wish it to persist beyond the next call to
  //   the package.  NULL is returned on an error.

  // An indel array gives the locations at which a dash (-) should be inserted into either the
  //   first or second sequence in order to expose the alignment between the two sequences
  //   implied by the trace point intervals.  A postive number x indicates one should place a
  //   dash before the x'th character of the first sequence, and a negative number -x indicates
  //   placing a dash before the x'th character of the second sequence.  The absolute values
  //   of the dash locations is increasing along the array and the array is terminated by a 0 value.

  // The indices in the indel array are relative to the (sub)sequences that are aligned.
  //   This includes complementing the B/second sequence if necessary, so that the position 1
  //   refers to the first character in the alignment for both sequences.

  // When reversed is false, the indel array is one that transforms the first sequence
  //   into the second where the first string is in the forward orientation.  If reversed
  //   is true, then the indel array is one that transforms the second sequence into the
  //   first where the second sequence is in the forward orientation, i.e. the roles of
  //   the first and second sequence are reversed.

int  *alnCreateIndelArray(AlnRecord *align, bool reversed);

  // Write a BLAST-like text display of an alignment to FILE where.  For this routine to work
  //   it must be that see_seq was true when the alignment reader was opened, and see_trace
  //   was true when the alignment record was fetched.  True is returned on error.

  // The display is indented by indent spaces, each displayed segment shows width columns of
  //   the alignment, border bp's before and after the aligned portion are also displayed, and
  //   the display width for coordinates is given by coord.  Upper_case controls the case of
  //   the base pairs.

  // When reversed is false, the 1st sequence is displayed above the second with the first
  //   sequence in the forward direction.  When reversed is true, the 2nd sequence is displayed
  //   above the first with the second sequence in the forward direction, i.e. the roles of
  //   the first and second sequence are reversed.

bool  alnShowAlignment(AlnRecord *align, FILE *where,
                       int indent, int width, int border, int coord,
                       bool upper_case, bool reversed);

/* 
   ERROR MESSAGES:

   Find below a list of all error message followed by a list of the routines that can
     potentially encounter the error.  Starred error message are controlled by CHECK_ARGS
     and so do not occur if this macro variable is set to undefined.

    Could not open <path>/<root>.1aln
    Out of memory (Opening .1aln file)
        alnOpenReader

  * Index out of range
  * Could not seek to <idx>th alignment
        alnGoto

  * Scaffold index out of range
        gdbScaffoldLen
        gdbContigLen
        gdbGapLen
        gdbScaffoldContigs
        gdbScaffoldStart
        gdbScaffoldName
        gdbScaffoldSeq

  * Contig index out of range
        gdbContigLen
        gdbScaffoldStart

  * Gap index out of range
        gdbGapLen

  * Sequence access not requested on open
    Could not seek GDB sequence file
    Could not read GDB sequence file
        gdbScaffoldSeq
        alnCreateCigar
        alnCreateCStag
        alnCreateIndelArray
        alnShowAlignment

  * Scaffold interval intersects a gap
    Out of memory (Allocating sequence)
        gdbScaffoldSeq

    T-line not present in alignment record
    T- and X-lines do not have the same length
    D-line not present in alignment record
        alnAlignment

    Bad alignment between trace points
    Trace point out of bounds
    Alignment end points not in band
    Self comparison can cross main diagonal
    Out of memory (Enlarging trace vector)
    Out of memory (Enlarging DP vector)
  * Trace information was not requested when record fetched
        alnCreateCigar
        alnCreateCStag
        alnCreateIndelArray
        alnShowAlignment

    Out of memory (Allocating CIGAR array)
    Out of memory (Expanding CIGAR array)
        alnCreateCigar
        alnCreateCStag

    Out of memory (Allocating CIGAR string)
        alnCreateCigar

    Out of memory (Allocating CS-tag string)
        alnCreateCStag

  In summary the following routines can return an error even if CHECK_ARGS is undefined, as they
    can encounter a read/seek failure, a memory allocation failure, or missing fields in a
    .1aln file

     alnOpenReader
     alnAlignment
     gdbScaffoldSeq
     alnCreateCigar
     alnCreateCStag
     alnCreateIndelArray
     alnShowAlignment

*/

#endif  //  _1ALN_DEFS
