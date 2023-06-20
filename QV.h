/*******************************************************************************************
 *
 *  Compressor/decompressor for .quiv files: customized Huffman codes for each stream based on
 *    the histogram of values occuring in a given file.  The two low complexity streams
 *    (deletionQV and substitutionQV) use a Huffman coding of the run length of the prevelant
 *    character.
 *
 *  Author:   Gene Myers
 *  Date:     Jan 18, 2014
 *  Modified: July 25, 2014
 *
 ********************************************************************************************/

#ifndef _QV_COMPRESSOR

#include <stdio.h>

#define _QV_COMPRESSOR

  //  The defined constant INTERACTIVE (set in DB.h) determines whether an interactive or
  //    batch version of the routines in this library are compiled.  In batch mode, routines
  //    print an error message and exit.  In interactive mode, the routines place the error
  //    message in EPLACE (also defined in DB.h) and return an error value, typically NULL
  //    if the routine returns a pointer, and an unusual integer value if the routine returns
  //    an integer.
  //  Below when an error return is described, one should understand that this value is returned
  //    only if the routine was compiled in INTERACTIVE mode.

  //  A PacBio compression scheme

typedef struct
  { void    *delScheme;   //  Huffman scheme for deletion QVs
    void    *insScheme;   //  Huffman scheme for insertion QVs
    void    *mrgScheme;   //  Huffman scheme for merge QVs
    void    *subScheme;   //  Huffman scheme for substitution QVs
    void    *dRunScheme;  //  Huffman scheme for deletion run lengths (if delChar > 0)
    void    *sRunScheme;  //  Huffman scheme for substitution run lengths (if subChar > 0)
    int      delChar;     //  If > 0, run-encoded deletion value
    int      subChar;     //  If > 0, run-encoded substitution value
    int      flip;        //  Need to flip multi-byte integers
    char    *prefix;      //  Header line prefix
  } QVcoding;

  // Read the next nlines of input, and QVentry returns a pointer to the first line if needed.
  //   If end-of-input is encountered before any further input, -1 is returned.  If there is
  //   an error than -2 is returned.  Otherwise the length of the line(s) read is returned.

int       Read_Lines(FILE *input, int nlines);
char     *QVentry();

  // Get and set the line counter for error reporting

void      Set_QV_Line(int line);
int       Get_QV_Line();

  // Read up to the next num entries or until eof from the .quiva file on input and record
  //   frequency statistics.  Copy these entries to the temporary file temp if != NULL.
  //   If there is an error then -1 is returned, otherwise the number of entries read.

int       QVcoding_Scan(FILE *input, int num, FILE *temp);
void      QVcoding_Scan1(int rlen, char *del, char *tag, char *ins, char *mrg, char *sub);

  // Given QVcoding_Scan has been called at least once, create an encoding scheme based on
  //   the accumulated statistics and return a pointer to it.  The returned encoding object
  //   is *statically allocated within the routine.  If lossy is set then use a lossy scaling
  //   for the insertion and merge streams.  If there is an error, then NULL is returned.

QVcoding *Create_QVcoding(int lossy);

  //  Read/write a coding scheme to input/output.  The encoding object returned by the reader
  //    is *statically* allocated within the routine.  If an error occurs while reading then
  //    NULL is returned.

QVcoding *Read_QVcoding(FILE *input);
void      Write_QVcoding(FILE *output, QVcoding *coding);

  //  Free all the auxiliary storage associated with coding (but not the object itself!)

void      Free_QVcoding(QVcoding *coding);

  //  Assuming the file pointer is positioned just beyond an entry header line, read the
  //    next set of 5 QV lines, compress them according to 'coding', and output.  If lossy
  //    is set then the scheme is a lossy one.  A negative value is returned if an error
  //    occurred, and the sequence length otherwise.

int      Compress_Next_QVentry(FILE *input, FILE *output, QVcoding *coding, int lossy);
void     Compress_Next_QVentry1(int rlen, char *del, char *tag, char *ins, char *mrg, char *sub,
                                FILE *output, QVcoding *coding, int lossy);

  //  Assuming the input is position just beyond the compressed encoding of an entry header,
  //    read the set of compressed encodings for the ensuing 5 QV vectors, decompress them,
  //    and place their decompressed values into entry which is a 5 element array of character
  //    pointers.  The parameter rlen computed from the preceeding header line, critically
  //    provides the length of each of the 5 vectors.  A non-zero value is return only if an
  //    error occured.

int      Uncompress_Next_QVentry(FILE *input, char **entry, QVcoding *coding, int rlen);

#endif // _QV_COMPRESSOR
