/******************************************************************************************
 *
 *  File: ONElib.h
 *    Header for ONE file reading and writing
 *
 *  Authors: Richard Durbin (rd109@cam.ac.uk), Gene Myers (myers@mpi-cbg.de)
 *  Copyright (C) Richard Durbin, Gene Myers, 2019-
 *
 * HISTORY:
 * Last edited: Dec  3 06:08 2022 (rd109)
 * * Dec  3 06:01 2022 (rd109): remove oneWriteHeader(), switch to stdarg for oneWriteComment etc.
 *   * Dec 27 09:46 2019 (gene): style edits
 *   * Created: Sat Feb 23 10:12:43 2019 (rd109)
 *
 *****************************************************************************************/

#ifndef ONE_DEFINED
#define ONE_DEFINED

#include <stdio.h>    // for FILE etc.
#include <stdarg.h>   // for formatted writing in oneWriteComment(), oneAddProvenance()
#include <inttypes.h> // for standard size int types and their PRI print macros
#include <stdbool.h>  // for standard bool types
#include <limits.h>   // for INT_MAX etc.
#include <pthread.h>

/***********************************************************************************
 *
 *    DATA TYPES
 *
 **********************************************************************************/

// Basic Types
#ifndef U8_DEFINED
#define U8_DEFINED

typedef int8_t        I8;
typedef int16_t       I16;
typedef int32_t       I32;
typedef int64_t       I64;
typedef unsigned char U8;

#endif // U8_DEFINED

typedef enum { oneINT = 1, oneREAL, oneCHAR, oneSTRING,
	       oneINT_LIST, oneREAL_LIST, oneSTRING_LIST, oneDNA } OneType;
extern char* oneTypeString[] ; 
// = { 0, "INT", "REAL", "CHAR", "STRING", "INT_LIST", "REAL_LIST", "STRING_LIST", "DNA" } ;

typedef union
  { I64    i;
    double r;
    char   c;
    I64    len; // For lists : top 8 bits encode excess bytes, low 56 length
  } OneField;

typedef struct
  { char *program;
    char *version;
    char *command;
    char *date;
  } OneProvenance;

typedef struct
  { char *filename; 
    I64   count;
  } OneReference;

typedef struct
  { I64 count;
    I64 max;
    I64 total;
    I64 groupCount;
    I64 groupTotal;
  } OneCounts;

  // OneCodecs are a private package for binary one file compression

typedef void OneCodec; // forward declaration of opaque type for compression codecs

  // DNAcodec is a special pre-existing compressor one should use for DNA.
  // It compresses every base to 2-bits, where any non-ACGT letter is
  // effectively converted to an A.  Compression is case insensitive,
  // but decompression always delivers lower-case.

extern  OneCodec *DNAcodec;

  // Record for a particular line type.  There is at most one list element.

typedef struct
  { OneCounts accum;            // counts read or written to this moment
    OneCounts given;            // counts read from header
    I64       gCount;           // used internally to calculate groupCount and groupTotal
    I64       gTotal;
    I64       oCount;           // # of objects in prefix before first group (if any)
    I64       oTotal;           // + of objects in prefix (these 2 are for thread parallel apps)

    int       nField;           // number of fields
    OneType  *fieldType;        // type of each field
    int       listEltSize;      // size of list field elements (if present, else 0)
    int       listField;        // field index of list
    char     *comment;          // the comment on the definition line in the schema
    
    bool      isUserBuf;        // flag for whether buffer is owned by user
    I64       bufSize;          // system buffer and size if not user supplied
    void     *buffer;

    OneCodec *listCodec;       // compression codec and flags
    bool      isUseListCodec;  // on once enough data collected to train associated codec
    char      binaryTypePack;   // binary code for line type, bit 8 set.
                                //     bit 0: list compressed
    I64       listTack;         // accumulated training data for this threads codeCodec (master)
  } OneInfo;

  // the schema type - the first record is the header spec, then a linked list of primary classes

typedef struct OneSchema
  {
    char      *primary ;
    int        nSecondary ;
    char     **secondary ;
    OneInfo   *info[128] ;
    int        nFieldMax ;
    char       objectType ;
    char       groupType ;
    struct OneSchema *nxt ;
  } OneSchema ;

typedef struct OneHeaderText
  { char *text ;
    struct OneHeaderText *nxt ;
  } OneHeaderText ;

  // The main OneFile type - this is the primary handle used by the end user

typedef struct
  {
    // this field may be set by the user

    bool           isCheckString;      // set if want to validate string char by char

    // these fields may be read by user - but don't change them!

    char          *fileType;
    char          *subType;
    char           lineType;           // current lineType
    char           objectType;         // line designation character for primary objects
    char           groupType;          // line designation character for groups (optional)
    I64            line;               // current line number
    I64            byte;               // current byte position when writing binary
    I64            object;             // current object - incremented when object line read
    I64            group;              // current group - incremented when group line read
    OneProvenance *provenance;         // if non-zero then count['!'] entries
    OneReference  *reference;          // if non-zero then count['<'] entries
    OneReference  *deferred;           // if non-zero then count['>'] entries
    OneField      *field;              // used to hold the current line - accessed by macros
    OneInfo       *info[128];          // all the per-linetype information
    I64            codecTrainingSize;  // amount of data to see before building codec

    // fields below here are private to the package

    FILE  *f;

    bool   isWrite;                // true if open for writing
    bool   isHeaderOut;            // true if header already written
    bool   isBinary;               // true if writing a binary file
    bool   inGroup;                // set once inside a group
    bool   isLastLineBinary;       // needed to deal with newlines on ascii files
    bool   isIndexIn;              // index read in
    bool   isBig;                  // are we on a big-endian machine?
    bool   isNoAsciiHeader;        // backdoor for ONEview to avoid writing header in ascii

    char   lineBuf[128];           // working buffers
    char   numberBuf[32];
    int    nFieldMax;
    I64    codecBufSize;
    char  *codecBuf;
    I64    nBits;                  // number of bits of list currently in codecBuf
    I64    intListBytes;           // number of bytes per integer in the compacted INT_LIST
    I64    linePos;                // current line position
    OneHeaderText *headerText;     // arbitrary descriptive text that goes with the header

    char   binaryTypeUnpack[256];  // invert binary line code to ASCII line character.
    int    share;                  // index if slave of threaded write, +nthreads > 0 if master
    int    isFinal;                // oneFinalizeCounts has been called on file
    pthread_mutex_t fieldLock;     // Mutexs to protect training accumumulation stats when threadded
    pthread_mutex_t listLock;
  } OneFile;                      //   the footer will be in the concatenated result.


/***********************************************************************************
 *
 *    ROUTINES FOR READING & WRITING ONE FILES IN BOTH ASCII & BINARY (TRANSPARENTLY)
 *
 **********************************************************************************/

//  CREATING AND DESTROYING SCHEMAS

OneSchema *oneSchemaCreateFromFile (char *path) ;
OneSchema *oneSchemaCreateFromText (char *text) ;

  // These functions create a schema handle that can be used to open One-code data files 
  //   for reading and writing.  A schema file is itself a One-code file, consisting of
  //   a set of objects, one per primary file type.  Valid lines in this file are:
  //      P <primary file type>   // a short string
  //      S <secondary file type> // a short string - any number of these
  //      O <char> <field_list>   // definition of object type
  //      G <char> <field_list>   // definition of group type - first field must be an int
  //      D <char> <field_list>   // definition of line
  //   <char> must be a lower or upper case letter.
  //   <field_list> is a list of field types from:
  //      CHAR, INT, REAL, STRING, INT_LIST, REAL_LIST, STRING_LIST, DNA
  //      Only one list type (STRING, *_LIST or DNA) is allowed per line type.
  //   All the D lines following an O line apply to that object type.
  //   By convention comments on each line explain the definition.
  //   Example, with lists and strings preceded by their length in OneCode style
  //      P 3 seq                            this is a sequence file
  //      O S 1 3 DNA                        the DNA sequence - each S line starts an object
  //      D Q 1 6 STRING                     the phred encoded quality score + ASCII 33
  //      D N 4 4 REAL 4 REAL 4 REAL 4 REAL  signal to noise ratio in A, C, G, T channels
  //      G g 2 3 INT 6 STRING               group designator: number of objects, name
  // The ...FromText() alternative writes the text to a temp file and reads it with 
  //   oneSchemaCreateFromFile(). This allows code to set the schema.
  // Internally a schema is a linked list of OneSchema objects, with the first holding
  //   the (hard-coded) schema for the header and footer, and the remainder each 
  //   corresponding to one primary file type.

void oneSchemaDestroy (OneSchema *schema) ;

//  READING ONE FILES:

OneFile *oneFileOpenRead (const char *path, OneSchema *schema, char *type, int nthreads) ;

  // Open ONE file 'path', either binary or ascii encoded, for reading.
  //   If the file doesn't have a header, then 'type' must be specified,
  //   otherwise, if 'type' is non-zero it must match the header type.
  //   All header information (if present) is read.
  // 'schema' is also optional.  If it is NULL then the file must contain its own schema.  
  //   If 'schema' is present then it must support 'type', and if the file contains its 
  //   own schema, then that must be a subset of the one for this type in 'schema'.
  // If nthreads > 1 then nthreadds OneFiles are generated as an array and the pointer
  //   to the first, called the master, is returned.  The other nthreads-1 files are
  //   called slaves.  The package routines are aware of when a OneFile argument is a
  //   slave or master in a parallel group.  The master recieves provenance, counts, etc.
  //   The slaves only read data and have the virtue of sharing indices and codecs with
  //   the master if relevant.

bool oneFileCheckSchema (OneFile *vf, char *textSchema) ;

  // Checks if file schema is consistent with text schema.  Mismatches are reported to stderr.
  // Filetype and all linetypes in text must match.  File schema can contain additional linetypes.
  // e.g. if (! oneFileCheckSchema (vf, "P 3 seq\nD S 1 3 DNA\nD Q 1 6 STRING\nD P 0\n")) die () ;
  // This is provided to enable a program to ensure that its assumptions about data layout
  // are satisfied.

char oneReadLine (OneFile *vf) ;

  // Read the next ONE formatted line returning the line type of the line, or 0
  //   if at the end of the data section.  The content macros immediately below are
  //   used to access the information of the line most recently read.

void   *_oneList (OneFile *vf) ;                // lazy codec decompression if required
void   *_oneCompressedList (OneFile *vf) ;      // lazy codec compression if required

#define oneInt(vf,x)        ((vf)->field[x].i)
#define oneReal(vf,x)       ((vf)->field[x].r)
#define oneChar(vf,x)       ((vf)->field[x].c)
#define _LF(vf)             ((vf)->info[(int)(vf)->lineType]->listField)
#define oneLen(vf)          ((vf)->field[_LF(vf)].len & 0xffffffffffffffll)
#define oneString(vf)       (char *) _oneList(vf)
#define oneDNAchar(vf)      (char *) _oneList(vf)
#define oneDNA2bit(vf)      (U8 *) _oneCompressedList(vf)
#define oneIntList(vf)      (I64 *) _oneList(vf)
#define oneRealList(vf)     (double *) _oneList(vf)
#define oneNextString(vf,s) (s + strlen(s) + 1)

  // Access field information.  The index x of a list object is not required as there is
  //   only one list per line, stored in ->buffer.
  //   A "string list" is implicitly supported, get the first string with oneString, and
  //   subsequent strings sequentially with oneNextString, e.g.:
  //
  //       char *s = oneString(vf);
  //       for (i = 0; i < oneLen(vf); i++)
  //         { // do something with i'th string
  //           s = oneNextString(vf,s);
  //         }

char *oneReadComment (OneFile *vf);

  // Can be called after oneReadLine() to read any optional comment text after the fixed fields.
  // Returns NULL if there is no comment.

//  WRITING ONE FILES:

OneFile *oneFileOpenWriteNew (const char *path, OneSchema *schema, char *type,
			      bool isBinary, int nthreads);
OneFile *oneFileOpenWriteFrom (const char *path, OneFile *vfIn,
			       bool isBinary, int nthreads);

  // Create a new oneFile that will be written to 'path'.  For the 'New' variant supply
  //   the file type, subtype (if non-zero), and whether it should be binary or ASCII.
  //   For the 'From' variant, specify binary or ASCII, schema and all other header 
  //   information is inherited from 'vfIn', where the count stats are from vfIn's 
  //   accumulation (assumes vfIn has been fully read or written) if 'useAccum is true, 
  //   and from vfIn's header otherwise.
  // If nthreads > 1 then nthreads OneFiles are generated as an array and the pointer
  //   to the first, called the master, is returned.  The other nthreads-1 files are
  //   called slaves.  The package routines are aware of when a OneFile argument is a
  //   slave or master in a parallel group.  The slaves are expected to only write data
  //   lines, with the master adding provenance, producing the header, and then some
  //   segment of the initial data lines.  Upon close the final result is effectively
  //   the concatenation of the master, followed by the output of each slave in sequence.

bool oneInheritProvenance (OneFile *vf, OneFile *source);
bool oneInheritReference  (OneFile *vf, OneFile *source);
bool oneInheritDeferred   (OneFile *vf, OneFile *source);

  // Add all provenance/reference/deferred entries in source to header of vf.  Must be
  //   called before first call to oneWriteLine.

bool oneAddProvenance (OneFile *vf, char *prog, char *version, char *format, ...);
bool oneAddReference  (OneFile *vf, char *filename, I64 count);
bool oneAddDeferred   (OneFile *vf, char *filename);

  // Append provenance/reference/deferred to header information.  Must be called before
  //   first call to oneWriteLine.

  // For ASCII output, if you want the header to contain count information then you must
  //   create and fill the relevant OneCounts objects before the first call to oneWriteLine.
  //   For BINARY output, the OneCounts information is accumulated and written automatically.

void oneWriteLine (OneFile *vf, char lineType, I64 listLen, void *listBuf);

  // Set up a line for output just as it would be returned by oneReadLine and then call
  //   this routine to output the line (ASCII or binary).
  // Use the macros above on the l.h.s. of assignments to fill fields (e.g. oneInt(vf,2) = 3).
  // For lists, give the length in the listLen argument, and either place the list data in your
  //   own buffer and give it as listBuf, or put in the line's buffer and set listBuf == NULL.

void oneWriteLineFrom (OneFile *vf, OneFile *source) ; // copies a line from source into vf
void oneWriteLineDNA2bit (OneFile *vf, char lineType, I64 listLen, U8 *dnaBuf);

// Minor variants of oneWriteLine().
// Use oneWriteLineDNA2bit for DNA lists if your DNA is already 2-bit compressed.

void oneWriteComment (OneFile *vf, char *format, ...); // can not include newline \n chars

  // Adds a comment to the current line. Extends line in ascii, adds special line type in binary.

// CLOSING FILES (FOR BOTH READ & WRITE)

void oneFileClose (OneFile *vf);

  // Close vf (opened either for reading or writing). Finalizes counts, merges theaded files,
  // and writes footer if binary. Frees all non-user memory associated with vf.

//  GOTO & BUFFER MANAGEMENT

void oneUserBuffer (OneFile *vf, char lineType, void *buffer);

  // A buffer is used to capture the list element of each line type that has one.
  //   This routine allows you to reassign the buffer to one you've allocated, or
  //   to revert to a default system buffer if 'buffer' = NULL.  The previous buffer
  //   (if any) is freed.  The user must ensure that a buffer they supply is large
  //   enough. BTW, this buffer is overwritten with each new line read of the given type.

bool oneGotoObject (OneFile *vf, I64 i);

  // Goto i'th object in the file. This only works on binary files, which have an index.

I64  oneGotoGroup  (OneFile *vf, I64 i);

  // Goto the first object in group i. Return the size (in objects) of the group, or 0
  //   if an error (i out of range or vf has not group type). Only works for binary files.

/***********************************************************************************
 *
 *    A BIT ABOUT THE FORMAT OF BINARY FILES
 *
 **********************************************************************************/

 //   <bin file> <- <ASCII Prolog> <$-line> <binary data> <footer> <^-line> <footer-size:int64>
 //
 // '$'-line flags file is binary and gives endian
 // The data block ends with a blank line consisting of '\n'
 //
 // EWM: Removed '-' line, simply write off_t to footer start
 //
 //   <ASCII Prolog> <- <'1'-line> [<'2'-line>] ( <'!'-line> | <'<'-line> | <'>'-line> )*
 //
 // The ASCII prolog contains the type, subtype, provenance, reference, and deferred lines
 //   in the ASCII format.  The ONE count statistic lines for each data line type are found
 //   in the footer along with binary ';' and ':' lines that encode their compressors as
 //   needed.  The footer also contains binary '&' and '*' lines that encode the object index
 //   and group indices, respectively.
 //
 //   <Binary line> <- <Binary line code + tags> <fields> [<list data>]
 //
 // Line codes are >= 128 for binary encoded lines.  The low two order bits of these are flags,
 //   so each binary-encoded line type has 4 codes and a table maps these to the ASCII code.
 //   Bit 0 indicates if the fields of the line type are compressed, and Bit 1 indicates if
 //   the list data (if present) is compressed.
 //
 // If a field is a list, then the field array element for that field is the list's length
 //   where the low 56 bits encode length, and the high 8 bits encode the # of high-order
 //   0-bytes in every list element if an INT_LIST (0 otherwise).

#endif  // ONE_DEFINED

/******************* end of file **************/
