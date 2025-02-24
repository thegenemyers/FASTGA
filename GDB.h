/*******************************************************************************************
 *
 *  Genome data base module.  Auxiliary routines to open and manipulate a data base for
 *    which the sequence and genome information are separated into two separate files, and the
 *    sequence is compressed into 2-bits for each base.  Derived from the daligner codes.
 *
 *  Author :  Gene Myers
 *  Date   :  May 2024
 *
 ********************************************************************************************/

#ifndef _GDB_DEFS

#define _GDB_DEFS

#include <stdio.h>
#include "gene_core.h"
#include "ONElib.h"

OneSchema *make_Seq_Schema();  //  Make a .1seq schema

/*******************************************************************************************
 *
 *  GDB IN-CORE DATA STRUCTURES
 *
 ********************************************************************************************/

typedef struct
  { int64   clen;   //  Length of the contig's sequence
    int64   sbeg;   //  Left index of the contig's sequence in its' scaffold
    int64   boff;   //  Offset (in bytes) of the contig sequence either in memory
                    //     or in the .bps file (see seqstate of GB)
    int     scaf;   //  Index/# of scaffold contig is in
  } GDB_CONTIG;

typedef struct
  { int64   slen;   //  Length of the scaffold (including gaps)
    int     fctg;   //  Index/# of first contig in scaffold
    int     ectg;   //  Index/# of last contig + 1
    int64   hoff;   //  Offset (in bytes) of scaffold header string in headers block
  } GDB_SCAFFOLD;

//  The GDB record holds all information about the current state of an active GDB.

typedef struct
  { int            nprov;
    OneProvenance *prov;

    int           nscaff;     //  # of scaffolds
    GDB_SCAFFOLD *scaffolds;  //  array [0..nscaff) of scaffold records

    int           ncontig;    //  # of contigs
    int64         maxctg;     //  length of maximum contig
    GDB_CONTIG   *contigs;    //  array [0..ncontig) of contig records

    int           hdrtot;     //  total bytes in header block
    char         *headers;    //  memory block of all headers, '\0'-terminated.

    char         *srcpath;    //  Absolute path to origin of GDB (a FASTA or 1-file)
    char         *seqpath;    //  filename of .bps file

    int64         seqtot;     //  total # of bases
    int           seqstate;   //  One of the 5 format constants below
    int           seqsrc;     //  One of the 4 file types below
    void         *seqs;       //  file pointer if EXTERNAL, mem pointer if not
                              //     NULL => not present

    float         freq[4];    //  frequency of A, C, G, T, respectively
  } GDB; 

  //  Sequence format types:

#define EXTERNAL  -1    //  seqs are on .bps file (or not present if seqs == NULL)
#define COMPRESSED 0    //  seqs are in memory, 2-bit compressed
#define NUMERIC    1    //  seqs are in memory, numeric (ACGT = 0123)
#define LOWER_CASE 2    //  seqs are in memory, lower case acgt
#define UPPER_CASE 3    //  seqs are in memory, upper case ACGT

  //  Source file types

#define IS_FA    0   // Fasta file, uncompressed
#define IS_FA_GZ 1   // Fasta file, compressed
#define IS_ONE   2   // 1-code sequence file
#define IS_GDB   3   // Only valid for Get_GDB_Paths


/*******************************************************************************************
 *
 *  GDB ROUTINES
 *
 ********************************************************************************************/

// When in interactive mode (see gene_code.h) all routines return an error value on error
//   and place an error message at EPLACE.
// In bacth mode the error message is output and the routine exits, shutting down the calling
//    program.

  //  Read or create a temporary GDB for the source file found at either <source> or
  //     <cpath>/<source>.  If num_bps > 0 then return an array of num_bps FILE pointers
  //     to the .bps file for the database (necessary as the .bps is already unlinked
  //     at when creating a temporary GDB).  If num_bps > 1 then it is the responsibility
  //     of the caller to free this array after closing all the units save for the first.
  //     The first FILE pointer is shared with the GDB and so will be closed when the GDB
  //     is closed.

FILE **Get_GDB(GDB *gdb, char *source, char *cpath, int num_bps);

  //  Interpret source & target arguments returning complete path + extension names of
  //    the source and target in spath & tpath.  Returns the type of source.
  //  If target == NULL then tpath has the same path & root as spath
  //  Else If target is a directory then tpath has the same path as target
  //  Else tpath has the same path and root as target.
  //  If no_gdb is non-zero then the source cannot be a gdb.
  //  In interactive mode returns a negative number if an error occurs.
  //  The user is responsible for freeing both tpath and spath.

int Get_GDB_Paths(char *source, char *target, char **spath, char **tpath, int no_gdb);

  // Create a GDB from the source file 'spath' which is of type 'ftype'.  The GDB is created
  //   in the record 'gdb' supplied by the user.  The GDB is in seqstate EXTERNAL.
  // num_bps = 0:
  //   Do not create a .bps sequence file
  // num_bps > 0:
  //   Create a .bps file with its name consistent with a GDB with file name 'tpath'.
  //   But if tpath == NULL then create a temporary uniquely named .bps that is
  //     unlinked (i.e. disappears on program exit).
  //   Return an array of num_bps open FILE pointers to the .bps file.
  // If num_bps > 1 then it is the responsibility of the caller to free this array after
  //   closing all the units save for the first.  The first FILE pointer is shared with
  //   the GDB and so will be closed when the GDB is closed.
  // In interactive mode a NULL value is returned if there is an error.

FILE **Create_GDB(GDB *gdb, char *spath, int ftype, int num_bps, char *tpath, int nthresh);

  // Open the given database "path" into the supplied GDB record "gdb".
  //   Initially the sequence data, if any, stays in the .bps file with a FILE pointer to it.
  // Interactive return values:
  //     0: Open of GDB proceeded without mishap
  //     1: The GDB could not be opened, a message why was placed in EPLACE.

int Read_GDB(GDB *gdb, char *spath);

  // The GDB moves its' sequence data from the .bps file to an in-memory block and also converts
  //   it to the requested format.  The boff field of contig records are adjusted, if necessary,
  //   so that they give the offset in the memory block of the associated string.
  // In interactive mode, 1 is returned on error, 0 otherwise.

int Load_Sequences(GDB *gdb, int stype);

  // Write the given gdb to the file 'tpath'.  The GDB must have seqstate EXTERNAL and tpath
  //   must be consistent with the name of the .bps file.

int Write_GDB(GDB *gdb, char *tpath);

  // Shut down an open 'gdb' by freeing all associated space and closing any open file pointers

void Close_GDB(GDB *gdb);


/*******************************************************************************************
 *
 *  CONTIG READING ROUTINES
 *
 ********************************************************************************************/

  // Allocate and return a buffer big enough for the largest contig in 'gdb'.
  //   If it cannot allocate memory then returns NULL (in interactive mode).
  //   **NB** free(x-1) if x is the value returned as *prefix* and suffix sentinal bytes
  //   are added for the convenience of the alignment algorithms.

char *New_Contig_Buffer(GDB *gdb);

  // Return a pointer to the i'th contig in gdb in the format specified by stype.
  // If buffer is NULL:
  //     The GDB's sequences should be in memory in the same format and a pointer directly
  //     to the sequence in this memory block is returned.
  // If buffer != NULL:
  //     The selected contig sequence is constructed in the memory block pointed at by
  //     buffer in the format requested.  This includes reading the compressed sequence
  //     from disk if this is where the GDB's sequences currently reside.
  // The sentinel bytes appropriate to the format are guaranteed to be present and correctly
  //      set.
  // If an error occurs then NULL is returned in interactive mode.

char *Get_Contig(GDB *gdb, int i, int stype, char *buffer);

  // Like Get_Contig except that only the interval from beg to end is returned.
  // The UNCOMPRESSED format is not possible.  If buffer is NULL then one must
  // note carefully that the piece will not be surrounded by sentinnel values.

char *Get_Contig_Piece(GDB *gdb, int i, int beg, int end, int stype, char *buffer);

#endif // _GDB_DEFS
