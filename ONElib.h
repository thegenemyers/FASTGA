/******************************************************************************************
 *
 *  File: ONElib.h
 *    Header for ONE file reading and writing
 *
 *  Authors: Richard Durbin (rd109@cam.ac.uk), Gene Myers (gene.myers@gmail.com)
 *  Copyright (C) Richard Durbin, Gene Myers, 2019-
 *
 * HISTORY:
 * Last edited: Dec  1 00:28 2024 (rd109)
 * * Dec  3 06:01 2022 (rd109): remove oneWriteHeader(), switch to stdarg for oneWriteComment etc.
 *   * Dec 27 09:46 2019 (gene): style edits
 *   * Created: Sat Feb 23 10:12:43 2019 (rd109)
 *
 *****************************************************************************************/

#ifndef ONE_DEFINED
#define ONE_DEFINED

#include <stdio.h>    // for FILE etc.
#include <stdarg.h>   // for formatted writing in oneWriteComment(), oneAddProvenance()
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

typedef signed char        I8;
typedef signed short       I16;
typedef signed int         I32;
typedef signed long long   I64;
typedef unsigned char      U8;

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
  } OneCounts;

typedef struct
  { I64  count, count0, maxCount ; // used for all contained types
    I64  total, total0, maxTotal ; // used for all contained list types
    char type ;
    bool isList ;
  } OneStat ;

  // OneCodecs are a private package for binary one file compression

typedef void OneCodec; // forward declaration of opaque type for compression codecs

  // DNAcodec is a special pre-existing compressor one should use for DNA.
  // It compresses every base to 2-bits, where any non-ACGT letter is
  // effectively converted to an A.  Compression is case insensitive,
  // but decompression always delivers lower-case.

extern  OneCodec *DNAcodec;

  // Record for a particular line type.  There is at most one list element.

typedef struct
  { bool      isObject;         // set if this is an object type (O in schema)
    I64      *index;            // index for objects
    I64       indexSize;        // size of the index, if present
    bool      contains[128];    // contains[k] is true if linetype k contained in this object
    OneStat  *stats;            // 0-terminated list of stats for all contained types within the object
    bool      isFirst;          // if set then set count0 for any objects closed by this linetype
    bool      isClosed;         // set if has been closed in this file
    OneCounts accum;            // counts read or written to this moment
    OneCounts given;            // counts read from header

    int       nField;           // number of fields
    OneType  *fieldType;        // type of each field
    int       listEltSize;      // size of list field elements (if present, else 0)
    int       listField;        // field index of list
    
    bool      isUserBuf;        // flag for whether buffer is owned by user
    I64       bufSize;          // system buffer and size if not user supplied
    void     *buffer;

    OneCodec *listCodec;        // compression codec and flags
    bool      isUseListCodec;   // on once enough data collected to train associated codec
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
    int        nFieldMax ;
    OneInfo   *info[128] ;
    OneInfo   *currentObject ;    // needed for parsing G lines
    int        nDefn ;            // number of O,D,G definition lines			    
    int        defnOrder[128] ;   // so can write out O,D,G lines in same order they were given
    char      *defnComment[128] ; // comment on the definition line
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

    char          *fileName;           // name of file
    char          *fileType;           // primary file type
    char          *subType;            // secondary file type
    char           lineType;           // current lineType
    I64            line;               // current line number
    I64            byte;               // current byte position when writing binary
    OneProvenance *provenance;         // if non-zero then count['!'] entries
    OneReference  *reference;          // if non-zero then count['<'] entries
    OneReference  *deferred;           // if non-zero then count['>'] entries
    OneField      *field;              // used to hold the current line - accessed by macros
    OneInfo       *info[128];          // all the per-linetype information
    int            nDefn;              // number of O,D,G definition lines			    
    int            defnOrder[128];     // so can write out O,D,G lines in same order they were given
    char          *defnComment[128] ;  // comment on the definition line
    I64            codecTrainingSize;  // amount of data to see before building codec

    // fields below here are private to the package

    FILE  *f;

    bool   isWrite;                // true if open for writing
    bool   isHeaderOut;            // true if header already written
    bool   isBinary;               // true if writing a binary file
    bool   inGroup;                // set once inside a group
    bool   isLastLineBinary;       // needed to deal with newlines on ascii files
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
    OneInfo *openObjects[128];     // stack of infos for open objects
    int    objectFrame;            // index into openObjects, pointing to current object info

    char   binaryTypeUnpack[256];  // invert binary line code to ASCII line character.
    int    share;                  // index if slave of threaded write, +nthreads > 0 if master
    bool   isFinal;                // oneFinalizeCounts has been called on file
    pthread_mutex_t fieldLock;     // Mutexs to protect training accumumulation stats when threaded
    pthread_mutex_t listLock;
    FILE* *tempReadFiles;          // array of file pointers to be used by oneFileReopen()
  } OneFile;                       // the footer will be in the concatenated result.


/***********************************************************************************
 *
 *    ROUTINES FOR READING & WRITING ONE FILES IN BOTH ASCII & BINARY (TRANSPARENTLY)
 *
 **********************************************************************************/

//  CREATING AND DESTROYING SCHEMAS

OneSchema *oneSchemaCreateFromFile (const char *path) ;
OneSchema *oneSchemaCreateFromText (const char *text) ;

  // These functions create a schema handle that can be used to open ONEcode data files 
  //   for reading and writing.  A schema file is itself a ONEcode file, consisting of
  //   a set of objects, one per primary file type.  Valid lines in this file are:
  //      P <primary file type>   // a short string
  //      S <secondary file type> // a short string - any number of these
  //      O <char> <field_list>   // definition of object type - these are indexed
  //      G <char> <field_list>   // definition of group type - indexed and accumulate stats
  //      D <char> <field_list>   // definition of normal line type
  //   <char> must be a lower or upper case letter. By convention upper case letters are used
  //      for objects and records within objects, and lower case letters for groups and records not
  //      assigned to objects, including global and group information.
  //   <field_list> is a list of field types from:
  //      CHAR, INT, REAL, STRING, INT_LIST, REAL_LIST, STRING_LIST, DNA
  //      Only one list type (STRING, *_LIST or DNA) is allowed per line type.
  //   All the D lines following an O line apply to that object.
  //   By convention comments on each schema definition line explain the definition.
  //   Example, with lists and strings preceded by their length as required in ONEcode
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

bool oneFileWriteSchema (OneFile *of, char *filename) ;

  // Utility to write the schema of an open oneFile in a form that can be read by
  //   oneSchemaCreateFromFile().

//  READING ONE FILES:

char* oneErrorString  (void) ; // gives information on errors for routines that fail
                               // e.g. if oneFileOpenRead() or oneFileOpenWrite() return NULL
                               // or oneFileCheckSchema*() returns false.

OneFile *oneFileOpenRead (const char *path, OneSchema *schema, const char *type, int nthreads) ;

  // Open ONE file 'path', either binary or ascii encoded, for reading.
  //   If the file doesn't have a header, then 'type' must be specified,
  //   otherwise, if 'type' is non-zero it must match a header ptype (primary or secondary).
  //   All header information (if present) is read.
  // 'schema' is also optional.  If it is NULL then the file must contain its own schema.  
  //   If 'schema' is present then it must support 'type', and if the file contains its 
  //   own schema, then that must match 'schema' where they share line types (there can be
  //   additional record types in both the file's schema, and in the argument 'schema').
  // If nthreads > 1 then nthreads OneFiles are generated as an array and the pointer
  //   to the first, called the master, is returned.  The other nthreads-1 files are
  //   called slaves.  The package routines are aware of when a OneFile argument is a
  //   slave or master in a parallel group.  The master recieves provenance, counts, etc.
  //   The slaves only read data and have the virtue of sharing indices and codecs with
  //   the master if relevant.

bool oneFileCheckSchema (OneFile *of, OneSchema *schema, bool isRequired) ;
bool oneFileCheckSchemaText (OneFile *of, const char *textSchema) ;

  // Checks if file schema is consistent with provided schema.  Mismatches are reported to stderr.
  // Filetype and all linetypes must match.  File schema can contain additional linetypes.
  // If isRequired is true then file schema must have all line types in supplied schema.
  // e.g. if (!oneFileCheckSchemaText (of, "P 3 seq\nD S 1 3 DNA\nD Q 1 6 STRING\nD P 0\n")) die () ;
  // This is provided to enable a program to ensure that its assumptions about data layout
  // are satisfied.
  // It is also used by oneFileOpenRead() with isRequired false to check consistency.

// ACCESSING GENERAL INFORMATION ABOUT THE CONTENTS OF A FILE:
 
bool  oneStats (OneFile *of, char lineType, I64 *count, I64 *max, I64 *total) ;

  // Report number of lines of specified lineType, maximum list length, total list length

bool  oneStatsContains (OneFile *of, char objectType, char lineType, I64 *maxCount, I64 *maxTotal) ;

  // Report the largest count of lineType within an objectType, and the highest total list length

#define oneFileName(of)         ((of)->fileName)
#define oneReferenceCount(of)   ((of)->info['<'] ? (of)->info['<']->accum.count : 0)

// READING DATA:

char oneReadLine (OneFile *of) ;

  // Read the next ONE formatted line returning the line type of the line, or 0
  //   if at the end of the data section.  The content macros immediately below are
  //   used to access the information of the line most recently read.

void   *_oneList (OneFile *of) ;                // lazy codec decompression if required
void   *_oneCompressedList (OneFile *of) ;      // lazy codec compression if required

#define oneInt(of,x)        ((of)->field[x].i)
#define oneReal(of,x)       ((of)->field[x].r)
#define oneChar(of,x)       ((of)->field[x].c)
#define _LF(of)             ((of)->info[(int)(of)->lineType]->listField)
#define oneLen(of)          ((of)->field[_LF(of)].len & 0xffffffffffffffll)
#define oneString(of)       (char *) _oneList(of)
#define oneDNAchar(of)      (char *) _oneList(of)
#define oneDNA2bit(of)      (U8 *) _oneCompressedList(of)
#define oneIntList(of)      (I64 *) _oneList(of)
#define oneRealList(of)     (double *) _oneList(of)
#define oneNextString(of,s) (s + strlen(s) + 1)

  // Access field information.  The index x of a list object is not required as there is
  //   only one list per line, stored in ->buffer.
  //   A "string list" is implicitly supported, get the first string with oneString, and
  //   subsequent strings sequentially with oneNextString, e.g.:
  //
  //       char *s = oneString(of);
  //       for (i = 0; i < oneLen(of); i++)
  //         { // do something with i'th string
  //           s = oneNextString(of,s);
  //         }

char *oneReadComment (OneFile *of);

  // Can be called after oneReadLine() to read any optional comment text after the fixed fields.
  // Returns NULL if there is no comment.

//  WRITING ONE FILES:

OneFile *oneFileOpenWriteNew (const char *path, OneSchema *schema, const char *type,
			      bool isBinary, int nthreads);
OneFile *oneFileOpenWriteFrom (const char *path, OneFile *ofIn,
			       bool isBinary, int nthreads);

  // Create a new oneFile that will be written to 'path'.  For the 'New' variant supply
  //   the file type, subtype (if non-zero), and whether it should be binary or ASCII.
  //   For the 'From' variant, specify binary or ASCII, schema and all other header 
  //   information is inherited from 'ofIn', where the count stats are from ofIn's 
  //   accumulation (assumes ofIn has been fully read or written) if 'useAccum is true, 
  //   and from ofIn's header otherwise.
  // If nthreads > 1 then nthreads OneFiles are generated as an array and the pointer
  //   to the first, called the master, is returned.  The other nthreads-1 files are
  //   called slaves.  The package routines are aware of when a OneFile argument is a
  //   slave or master in a parallel group.  The slaves are expected to only write data
  //   lines, with the master adding provenance, producing the header, and then some
  //   segment of the initial data lines.  Upon close the final result is effectively
  //   the concatenation of the master, followed by the output of each slave in sequence.
  //   Although slave files are created, they are removed from the file system by unlinking
  //   them immediately after creation - this means they are cleaned up automatically
  //   when the OneFile is closed, or if the program crashes.
  // If the path is a directory, then this creates a temporary file in the directory and
  //   immediately unlinks it. Access to this file must be via oneFileReopenRead() after
  //   completion of writing it; this returns a OneFile (array) like oneFileOpenRead(), 
  //   with the size of the array determined by value of nthreads in the original call
  //   to oneFileOpenWriteXXX().

OneFile *oneFileReopenRead (OneFile *of);  // see end of preceding paragraph

bool oneInheritProvenance (OneFile *of, OneFile *source);
bool oneInheritReference  (OneFile *of, OneFile *source);
bool oneInheritDeferred   (OneFile *of, OneFile *source);

  // Add all provenance/reference/deferred entries in source to header of of.  Must be
  //   called before first call to oneWriteLine.

bool oneAddProvenance (OneFile *of, const char *prog, const char *version, char *format, ...);
bool oneAddReference  (OneFile *of, const char *filename, I64 count);
bool oneAddDeferred   (OneFile *of, const char *filename);

  // Append provenance/reference/deferred to header information.  Must be called before
  //   first call to oneWriteLine.

  // For ASCII output, if you want the header to contain count information then you must
  //   create and fill the relevant OneCounts objects before the first call to oneWriteLine.
  //   For BINARY output, the OneCounts information is accumulated and written automatically.

void oneWriteLine (OneFile *of, char lineType, I64 listLen, void *listBuf);

  // Set up a line for output just as it would be returned by oneReadLine and then call
  //   this routine to output the line (ASCII or binary).
  // Use the macros above on the l.h.s. of assignments to fill fields (e.g. oneInt(of,2) = 3).
  // For lists, give the length in the listLen argument, and either place the list data in your
  //   own buffer and give it as listBuf, or put in the line's buffer and set listBuf == NULL.

void oneWriteLineFrom (OneFile *of, OneFile *source) ; // copies a line from source into of
void oneWriteLineDNA2bit (OneFile *of, char lineType, I64 listLen, U8 *dnaBuf);

// Minor variants of oneWriteLine().
// Use oneWriteLineDNA2bit for DNA lists if your DNA is already 2-bit compressed.

void oneWriteComment (OneFile *of, char *format, ...); // can not include newline \n chars

  // Adds a comment to the current line. Extends line in ascii, adds special line type in binary.

// CLOSING FILES (FOR BOTH READ & WRITE):

void oneFileClose (OneFile *of);

  // Close of (opened either for reading or writing). Finalizes counts, merges theaded files,
  // and writes footer if binary. Frees all non-user memory associated with of.

//  GOTO & BUFFER MANAGEMENT:

void oneUserBuffer (OneFile *of, char lineType, void *buffer);

  // A buffer is used to capture the list element of each line type that has one.
  //   This routine allows you to reassign the buffer to one you've allocated, or
  //   to revert to a default system buffer if 'buffer' = NULL.  The previous buffer
  //   (if any) is freed.  The user must ensure that a buffer they supply is large
  //   enough. BTW, this buffer is overwritten with each new line read of the given type.

#define oneObject(of,i)  ((of) && (of)->info[i] ? (of)->info[i]->accum.count : -1)

  // Returns the number of the object of type lineType currently in.  Works in read
  // or write mode.  0 if no objects of this type read/written yet. -1 if lineType illegal.

bool oneGoto (OneFile *of, char lineType, I64 i);

  // Goto i'th object in the file. This only works on binary files, which have an index.
  // The first object is numbered 1. Setting i == 0 goes to the first data line of the file
  // after the header.

I64 oneCountUntilNext (OneFile *of, char countType, char nextType) ;

  // returns the number of countType object lines before the next nextType object line
  // returns -1 on error, e.g. not reading a binary file, types are not object types

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
 //   <ASCII Prolog> <- <'1'-line> [<'2'-line>] ( <'!'-line> | <'<'-line> | <'>'-line> )*
 //
 // The ASCII prolog contains the type, subtype, provenance, reference, and deferred lines
 //   in the ASCII format.  The ONE count statistic lines for each data line type are found
 //   in the footer along with binary ';' lines that encode their compressors as needed.
 //   The footer also contains binary '&' lines that encode the byte index for object types.
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
