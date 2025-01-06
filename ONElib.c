/*****************************************************************************************
 *
 *  File: ONElib.c
 *    implementation for ONElib.h
 *
 *  Author: Richard Durbin (rd109@cam.ac.uk), Gene Myers (gene.myers@gmail.com)
 *  Copyright (C) Richard Durbin, Cambridge University and Eugene Myers 2019-
 *
 * HISTORY:
 * Last edited: Dec  8 16:30 2024 (rd109)
 * * May  1 00:23 2024 (rd109): moved to OneInfo->index and multiple objects/groups
 * * Apr 16 18:59 2024 (rd109): major change to object and group indexing: 0 is start of data
 * * Mar 11 02:49 2024 (rd109): fixed group bug found by Gene
 * * Mar 11 02:48 2024 (rd109): added oneFileWriteSchema() to write schema files for bare text parsing
 * * Mar 10 07:16 2024 (rd109): changed oneOpenFileRead semantics to prioritize file schema
 * * Dec 20 21:29 2022 (rd109): changed DNA compression to little-endian: natural on Intel, Apple
 * * Apr 23 00:31 2020 (rd109): global rename of VGP to ONE, Vgp to One, vgp to one
 * * Apr 20 11:27 2020 (rd109): added VgpSchema to make schema dynamic
 * * Dec 27 09:46 2019 (gene):  style edits + compactify code
 * * Jul  8 04:28 2019 (rd109): refactored to use info[]
 * * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *
 ****************************************************************************************/

#ifdef __linux__
#define _GNU_SOURCE  // needed for vasprintf() on Linux
#endif

#include <sys/errno.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <stdarg.h>
#include <time.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/uio.h>
#include <sys/stat.h>
#include <math.h>

#define DEBUG
#ifdef DEBUG
#include <assert.h>
#else
#define assert(x) 0
#endif

#include "ONElib.h"

// set major and minor code versions

#define MAJOR 2
#define MINOR 1

//  utilities with implementation at the end of the file

static void  die(char *format, ...);                  //  print message to stderr and exit -1
static void *myalloc(size_t size);                    //  allocate block, die if malloc fails
static void *mycalloc(size_t number, size_t size);    //  allocate & zero # objects of size
static void *mydup(size_t n, void *x, size_t size);
#define new(n,Type)  (Type *) myalloc((n)*sizeof(Type))    // use these not myalloc
#define new0(n,Type)  (Type *) mycalloc((n),sizeof(Type))
#define dup(n,x,Type) (Type *) mydup(n,x,sizeof(Type))
#define resize(x,nOld,nNew,Type) { Type* z = new((nNew),Type) ; if (nOld < nNew) memcpy(z,x,(nOld)*sizeof(Type)) ; else memcpy(z,x,(nNew)*sizeof(Type)) ; free(x) ; x = z ; }

// global required for parallelisation

static pthread_mutex_t mutexInit = PTHREAD_MUTEX_INITIALIZER;

// forward declarations of serialisation functions lower in the file
// RD 220818: I think that many of int below should be I64, e.g. for len, ilen etc.

OneCodec *vcCreate();
void      vcAddToTable(OneCodec *vc, int len, char *bytes);
void      vcAddHistogram(OneCodec *vc, OneCodec *vh);
void      vcCreateCodec(OneCodec *vc, int partial);
void      vcDestroy(OneCodec *vc);
int       vcMaxSerialSize();
int       vcSerialize(OneCodec *vc, void *out);
OneCodec *vcDeserialize(void *in);
int       vcEncode(OneCodec *vc, int ilen, char *ibytes, char *obytes);
int       vcDecode(OneCodec *vc, int ilen, char *ibytes, char *obytes);

// forward declarations of 64-bit integer encoding/decoding

static inline int ltfWrite (I64 x, FILE *f) ;
static inline I64 ltfRead (FILE *f) ;

// error handling

static char errorString[1024] = "" ;

char *oneErrorString (void) { return errorString ; }

/***********************************************************************************
 *
 *    ONE_FILE CREATION & DESTRUCTION
 *
 **********************************************************************************/

char* oneTypeString[] = { 0, "INT", "REAL", "CHAR", "STRING", 
			  "INT_LIST", "REAL_LIST", "STRING_LIST", "DNA" } ;

/******************* OneInfo ********************/

static OneInfo *infoCreate (int nField)
{ OneInfo *vi = new0 (1, OneInfo) ;
  vi->nField = nField ;
  if (nField) vi->fieldType = new (nField, OneType) ;
  return vi;
}

static OneInfo *infoDeepCopy (OneInfo *vi0)
{ OneInfo *vi = new (1, OneInfo) ;
  *vi = *vi0 ;
  if (vi0->nField) vi->fieldType = dup (vi->nField, vi0->fieldType, OneType) ;
  if (vi0->listCodec && vi->listCodec != DNAcodec) vi->listCodec = vcCreate() ;
  if (vi0->index) vi->index = dup (vi->indexSize, vi0->index, I64) ;
  if (vi0->stats)
    { int n = 1 ; OneStat *s ; for (s = vi->stats ; s->type ; ++s) ++n ;
      vi->stats = dup (n, vi0->stats, OneStat) ;
    }
  return vi ;
}

static void infoDestroy (OneInfo *vi)
{ if (vi->buffer && ! vi->isUserBuf) free (vi->buffer) ;
  if (vi->listCodec) vcDestroy (vi->listCodec) ;
  if (vi->fieldType) free (vi->fieldType) ;
  if (vi->index) free (vi->index) ;
  if (vi->stats) free (vi->stats) ;
  free (vi);
}

/******************* OneSchema ********************/

static void schemaAddGroup (OneSchema *vs, char t)
{
  if (!vs->currentObject) die ("G %c line when no object defined", t) ;
  vs->currentObject->contains[(int)t] = true ;
  vs->defnOrder[vs->nDefn++] = t | 0x80 ; // definition order
  return ;
}

// a utility to set the OneInfo list information

static int listEltSize[9] = { 0, 0, 0, 0, 1, sizeof(I64), sizeof(double), 1, 1 } ;

static void schemaAddInfoFromArray (OneSchema *vs, int n, OneType *a, char t, char type)
{
  // use during the bootstrap, while parsing schema files, and while parsing ~ lines in other files
  
  if (vs->info[(int) t])
    die ("duplicate schema specification for linetype %c in filetype %s", t, vs->primary) ;
  
  OneInfo *vi = infoCreate (n) ;
  vs->info[(int)t] = vi ;
  if ((t >= 'A' && t <= 'Z') || (t >= 'a' && t <= 'z'))
    vs->defnOrder[vs->nDefn++] = t ; // definition order

  if (isalpha(t) && type == 'O')
    { vi->isObject = true ;
      vs->currentObject = vi ;
    }
  else
    { if (type != 'D')
	die ("type %c not 'O' or 'D' in schemaAddInfo", type) ;
      if (vs->primary && !isalpha(t)) // allow non-alphabetic lines in header
	die ("non-alphabetic linetype %c (ascii %d) in schema for filetype %s",t,t,vs->primary) ;
      if (vs->currentObject)
	vs->currentObject->contains[(int)t] = true ;
    }
    
  if (n > vs->nFieldMax) vs->nFieldMax = n ;
  
  memcpy (vi->fieldType, a, n*sizeof(OneType)) ;
  int i ;
  for (i = 0 ; i < n ; ++i)
    if (a[i] >= oneSTRING)
      { if (vi->listEltSize)
	  die ("OneFile schema error; multiple list types for linetype definition %c", t) ;
	vi->listEltSize = listEltSize[vi->fieldType[i]] ;
	vi->listField = i ;
	if (a[i] == oneDNA)
	  { vi->listCodec = DNAcodec ; vi->isUseListCodec = true ; }
	else if (t != '/') // make a listCodec for any list type except for commments
	  vi->listCodec = vcCreate () ; 
      }

  // need a binary packing code for any linetype that might appear in binary context
  if (t >= 'A' && t <= 'Z') vi->binaryTypePack = ((t-'A') << 1) | (char) 0x80 ;
  else if (t >= 'a' && t <= 'z') vi->binaryTypePack = ((26+t-'a') << 1) | (char) 0x80 ;
  else if (t == ';') vi->binaryTypePack = (52 << 1) | (char) 0x80 ; // list codec
  else if (t == '&') vi->binaryTypePack = (53 << 1) | (char) 0x80 ; // byte index
  else if (t == '/') vi->binaryTypePack = (54 << 1) | (char) 0x80 ; // comment - binary only
  else if (t == '.') vi->binaryTypePack = (55 << 1) | (char) 0x80 ; // blank line
  // don't need for #, +, @, % because these lines are always written in ASCII
}

static void schemaAddInfoFromLine (OneSchema *vs, OneFile *vf, char t, char type)
{ // assumes field specification is in the STRING_LIST of the current vf line
  // need to set vi->comment separately
  
  static OneType a[32] ;
  int            i ;
  OneType        j ;
  char          *s = oneString(vf) ;
  int            n = oneLen(vf) ;
  
  if (n > 32)
    die ("line specification %d fields too long - need to recompile", n) ;

  for (i = 0 ; i < n ; ++i, s = oneNextString(vf,s))
    { a[i] = 0 ;
      for (j = oneINT ; j <= oneDNA ; ++j)
	if (!strcmp (s, oneTypeString[j])) a[i] = j ;
      if (!a[i])
	die ("ONE schema error: bad field %d of %d type %s in line %d type %c",
	     i, n, s, vf->line, t) ;
    }

  schemaAddInfoFromArray (vs, n, a, t, type) ;  
}

static OneSchema *schemaLoadRecord (OneSchema *vs, OneFile *vf)
{ char *s;

  // parse a schema specfication line from vf and add into vs
  // return value is vs unless a new primary type is declared, in which case vs->nxt

  switch (vf->lineType)
    {
    case '.':  // ignore - blank or comment line in schema file
      break ;
    case 'P':
      if (oneLen(vf) == 0) die ("schema: primary name must have at least one letter") ;
      OneSchema *vsNxt = new0 (1, OneSchema) ;
      vs->nxt = vsNxt ;
      vs = vsNxt ;
      s = oneString(vf);
      vs->primary = new0 (oneLen(vf)+1, char) ;
      strcpy (vs->primary, s) ;
      vs->nFieldMax = 4 ; // needed for header
      break ;
    case 'S':
      if (oneLen(vf) == 0) die ("schema: secondary name must have at least one letter") ;
      if (vs->nSecondary)
	{ char **temp = vs->secondary ;
	  vs->secondary = new (vs->nSecondary+1, char*) ;
	  memcpy (vs->secondary, temp, vs->nSecondary*sizeof(char*)) ;
	  free (temp) ;
	}
      else
	vs->secondary = new (1, char*) ;
      s = oneString(vf);
      vs->secondary[vs->nSecondary] = new0 (oneLen(vf)+1, char) ;
      strcpy (vs->secondary[vs->nSecondary++], s) ;
      break ;
    case 'G': // group another object type
      schemaAddGroup (vs, oneChar(vf,0)) ;
      if (oneReadComment (vf)) vs->defnComment[vs->nDefn-1] = strdup (oneReadComment(vf)) ;
      break ;
    case 'O': // object type
    case 'D': // standard record type
      schemaAddInfoFromLine (vs, vf, oneChar(vf,0), vf->lineType) ;
      if (oneReadComment (vf)) vs->defnComment[vs->nDefn-1] = strdup (oneReadComment(vf)) ;
      break ;
    default:
      die ("unrecognized schema line %d starting with %c", vf->line, vf->lineType) ;
    }

  return vs ;
}

static void oneFileDestroy (OneFile *vf) ; // need a forward declaration here

static bool isBootStrap = false ;

OneSchema *oneSchemaCreateFromFile (const char *filename)
{
  FILE *fs = fopen (filename, "r") ;
  if (!fs) return 0 ;
  fclose(fs);

  OneSchema *vs = new0 (1, OneSchema) ;

  isBootStrap = true ;
  OneFile *vf = new0 (1, OneFile) ;      // shell object to support bootstrap
  // bootstrap specification of linetypes to read schemas
  { OneInfo *vi ;
    vi = vf->info['P'] = infoCreate (1) ;  // to define the schema for parsing a .schema file
    vi->fieldType[0] = oneSTRING ; vi->listEltSize = 1 ; vi->listField = 0 ;
    vi = vf->info['O'] = infoCreate (2) ;  // object type specification
    vi->fieldType[0] = oneCHAR ;
    vi->fieldType[1] = oneSTRING_LIST ; vi->listEltSize = 1 ; vi->listField = 1 ;
    vi = vf->info['D'] = infoCreate (2) ;  // line type specification
    vi->fieldType[0] = oneCHAR ;
    vi->fieldType[1] = oneSTRING_LIST ; vi->listEltSize = 1 ; vi->listField = 1 ;
    vf->info['/'] = infoCreate (0) ;       // to store comments
    vf->field = new (2, OneField) ;
  }

  // first load the universal header and footer (non-alphabetic) line types 
  // do this by writing their schema into a temporary file and parsing it into the base schema
  { errno = 0 ;
    static char template[64] ;
#define VALGRIND_MACOS
#ifdef VALGRIND_MACOS // MacOS valgrind is missing functions to make temp files it seems
    sprintf (template, "/tmp/OneSchema.%d", getpid()) ;
    vf->f = fopen (template, "w+") ;
    if (errno) die ("failed to open temporary file %s errno %d\n", template, errno) ;
#else
    strcpy (template, "/tmp/OneSchema.XXXXXX") ;
    int fd = mkstemp (template) ;
    if (errno) die ("failed to open temporary file %s errno %d\n", template, errno) ;
    vf->f = fdopen (fd, "w+") ;
    if (errno) die ("failed to assign temporary file to stream: errno %d\n", errno) ;
#endif
    unlink (template) ;                  // this ensures that the file is removed on closure
    if (errno) die ("failed to remove temporary file %s errno %d\n", template, errno) ;
  }

  // NB if you change the header spec and add a record with more than 4 fields,
  //    change the assignment of ->nFieldMax in the 'P' section of schemaLoadRecord() above
  
  fprintf (vf->f, "D 1 3 6 STRING 3 INT 3 INT         line 1: primary type, major, minor version\n") ;
  fprintf (vf->f, "D 2 1 6 STRING                     optional subtype: subtype\n") ;
  fprintf (vf->f, "D # 2 4 CHAR 3 INT                 count: linetype, count\n") ;
  fprintf (vf->f, "D @ 2 4 CHAR 3 INT                 max: linetype, list max\n") ;
  fprintf (vf->f, "D + 2 4 CHAR 3 INT                 total: linetype, list total\n") ;
  fprintf (vf->f, "D %% 4 4 CHAR 4 CHAR 4 CHAR 3 INT  group maxes: group, #/+, linetype, value\n") ;
  fprintf (vf->f, "D ! 1 11 STRING_LIST               provenance: program, version, command, date\n") ;
  fprintf (vf->f, "D < 2 6 STRING 3 INT               reference: filename, object count\n") ;
  fprintf (vf->f, "D > 1 6 STRING                     deferred: filename\n") ;
  fprintf (vf->f, "D ~ 3 4 CHAR 4 CHAR 11 STRING_LIST embedded schema linetype definition\n") ;
  fprintf (vf->f, "D . 0                              blank line, anywhere in file\n") ;
  fprintf (vf->f, "D $ 1 3 INT                        binary file - goto footer: isBigEndian\n") ;
  fprintf (vf->f, "D ^ 0                              binary file: end of footer designation\n") ;
  fprintf (vf->f, "D - 1 3 INT                        binary file: offset of start of footer\n") ;
  fprintf (vf->f, "D & 2 4 CHAR 8 INT_LIST            binary file: li->index\n") ;
  fprintf (vf->f, "D ; 2 4 CHAR 6 STRING              binary file: list codec\n") ;
  fprintf (vf->f, "D / 1 6 STRING                     binary file: comment\n") ;
  if (fseek (vf->f, 0, SEEK_SET)) die ("ONE schema failure: cannot rewind tmp file") ;
  while (oneReadLine (vf))
    schemaLoadRecord (vs, vf) ;

  // next reuse the temp file to load the schema for reading schemas
  if (fseek (vf->f, 0, SEEK_SET)) die ("ONE schema failure: cannot rewind tmp file") ;
  fprintf (vf->f, "P 3 def                      this is the primary file type for schemas\n") ;
  fprintf (vf->f, "O P 1 6 STRING               primary type name\n") ;
  fprintf (vf->f, "D S 1 6 STRING               secondary type name\n") ;
  fprintf (vf->f, "D O 2 4 CHAR 11 STRING_LIST  define linetype for object type (indexed)\n") ;
  fprintf (vf->f, "D G 1 4 CHAR                 define linetype for grouping another object\n") ;
  fprintf (vf->f, "D D 2 4 CHAR 11 STRING_LIST  define linetype for other records\n") ;
  fprintf (vf->f, "\n") ; // terminator
  if (fseek (vf->f, 0, SEEK_SET)) die ("ONE schema failure: cannot rewind tmp file") ;
  OneSchema *vs0 = vs ;  // need this because loadInfo() updates vs on reading P lines
  vf->line = 0 ;
  while (oneReadLine (vf))
    vs = schemaLoadRecord (vs, vf) ;
  OneSchema *vsDef = vs ; // will need this to destroy it once the true schema is read
  oneFileDestroy (vf) ;   // this will also effectively remove the temp file on closing

  // finally read the schema itself
  if (!(vf = oneFileOpenRead (filename, vs0, "def", 1)))
    return 0 ;
  vs = vs0 ; // set back to vs0, so next filetype spec will replace vsDef
  vs->nxt = 0 ;
  oneSchemaDestroy (vsDef) ; // no longer need this, and can destroy because unlinked from vs0
  while (oneReadLine (vf))
    vs = schemaLoadRecord (vs, vf) ;
  oneFileDestroy (vf) ;

  isBootStrap = false ;
  
  return vs0 ;
}

static char *schemaFixNewlines (const char *text)
{ // replace literal "\n" by '\n' chars in text
  char *newText = strdup (text) ;
  char *s = newText, *t = s ;
  while (*s)
    if (*s == '\\' && s[1] == 'n')
      { *t++ = '\n' ; s += 2 ; }
    else
      *t++ = *s++ ;
  *t = 0 ;
  return newText ;
}
  
OneSchema *oneSchemaCreateFromText (const char *text) // write to temp file and call CreateFromFile()
{
  static char template[64] ;
  sprintf (template, "/tmp/OneTextSchema-%d.schema", getpid()) ;

  errno = 0 ;
  FILE *f = fopen (template, "w") ;
  if (!f) die ("failed to open temporary file %s for writing schema to", template) ;
  char *fixedText = schemaFixNewlines (text) ;
  char *s = fixedText ;
  while (*s && *s != 'P')
    { while (*s && *s != '\n') ++s ;
      if (*s == '\n') ++s ;
    }
  if (!*s) die ("no P line in schema text") ;
  fprintf (f, "%s\n", s) ;
  free (fixedText) ;
  fclose (f) ;
  if (errno) die ("failed to write temporary file %s errno %d\n", template, errno) ;

  OneSchema *vs = oneSchemaCreateFromFile (template) ;

  errno = 0 ;
  unlink (template) ;  // delete temporary file - not ideal: will remain if schemaCreate crashes
  if (errno) die ("failed to remove temporary file %s errno %d\n", template, errno) ;

  return vs ;
}

static OneSchema *oneSchemaCreateDynamic (char *fileType, char *subType)
{ // this is clean, but it seems a bit wasteful to create a temp file
  char *text ;
  assert (fileType && strlen(fileType) > 0) ;
  assert (!subType || strlen(subType) > 0) ;
  text = new (32 + strlen(fileType) + (subType?strlen(subType):0), char) ;
  if (subType)
    sprintf (text, "P %ld %s\nS %ld %s\n", strlen(fileType),fileType, strlen(subType), subType) ;
  else
    sprintf (text, "P %ld %s\n", strlen(fileType), fileType) ;
  OneSchema *vs = oneSchemaCreateFromText (text) ;
  free (text) ;
  return vs ;
}

void oneSchemaDestroy (OneSchema *vs)
{ int i ;
  while (vs)
    { for (i = 0 ; i < 128 ; ++i) if (vs->info[i]) infoDestroy (vs->info[i]) ;
      if (vs->nSecondary)
	{ for (i = 0 ; i < vs->nSecondary ; ++i) free (vs->secondary[i]) ;
	  free (vs->secondary) ;
	}
      free(vs->primary);
      for (i = 0 ; i < vs->nDefn ; ++i)
	if (vs->defnComment[i]) free (vs->defnComment[i]) ;
      OneSchema *t = vs->nxt ;
      free (vs) ;
      vs = t ;
    }
}

static void writeInfoSpec (FILE *f, OneFile *vf, char ci, char *comment) // also used in writeHeader()
{
  if (f == vf->f) fprintf (f, "\n~ ") ; // writing the schema into the file header
  else fprintf (f, "\n") ;              // just writing a schema file

  if (ci & 0x80)
    fprintf (f, "G %c 0", ci & 0x7f) ;
  else
    { OneInfo *vi = vf->info[(int) ci] ;
      if (vi->isObject)
	fprintf (f, "O %c %d", ci, vi->nField) ;
      else
	fprintf (f, "D %c %d", ci, vi->nField) ;
      int i ;
      for (i = 0 ; i < vi->nField ; ++i)
	fprintf (f, " %d %s",
		 (int)strlen(oneTypeString[vi->fieldType[i]]), oneTypeString[vi->fieldType[i]]) ;
    }
  if (comment)
    fprintf (f, " %s", comment) ;
}

bool oneFileWriteSchema (OneFile *vf, char *filename)
{
  int i ;
  FILE *f ;

  if (!vf) die ("oneFileWriteSchema() passed a NULL OneFile object") ;
  if (!strcmp (filename, "-"))
    f = stdout ;
  else if (!(f= fopen (filename, "w")))
    { snprintf (errorString, 1024, "failed to open %s to write schema into\n", filename) ;
      return false ;
    }

  fprintf (f, "P %d %s", (int)strlen(vf->fileType), vf->fileType) ;
  if (vf->subType) fprintf (f, "\nS %d %s", (int)strlen(vf->subType), vf->subType) ;

  for (i = 0 ; i < vf->nDefn ; ++i)
    writeInfoSpec (f, vf, vf->defnOrder[i], vf->defnComment[i]) ;
  
  fprintf (f, "\n") ;
  fclose (f) ;
  return true ;
}

/*************************************/

static void initialiseStats (OneFile *vf)
{
  int      i, j, k ;
  OneInfo *li, *lj ;

  // first ensure the contains[] arrays follow the defn lines
  if (!vf->info[0]) { vf->info[0] = infoCreate(0) ; vf->info[0]->isObject = true ; }
  OneInfo *currentInfo = vf->info[0] ;
  for (i = 0 ; i < vf->nDefn ; ++i)
    { k = vf->defnOrder[i] ;
      if (k & 0x80) currentInfo->contains[k & 0x7f] = true ;
      else if (vf->info[k]->isObject) currentInfo = vf->info[k] ;
      else currentInfo->contains[k] = true ;
    }

  // next ensure the contains[] arrays are complete by recursion
  bool isDone = false ;
  while (!isDone)
    { isDone = true ; // set to false if we have to change anything
      for (i = 'A' ; i <= 'z' ; ++i)
	if (vf->info[i] && vf->info[i]->isObject)
	  { li = vf->info[i] ;
	    for (j = 'A' ; j <= 'z' ; ++j)
	      if (li->contains[j] && vf->info[j] && vf->info[j]->isObject)
		{ lj = vf->info[j] ;
		  for (k = 'A' ; k <= 'z' ; ++k)
		    if (lj->contains[k] && !li->contains[k])
		      { isDone = false ;
			li->contains[k] = true ;
		      }
		}
	  }
    }

  // next create the stats arrays
  for (i = 'A' ; i <= 'z' ; ++i)
    if (vf->info[i] && vf->info[i]->isObject)
      { li = vf->info[i] ;
	if (li->stats) free (li->stats) ;
	int n = 0 ; for (j = 'A' ; j <= 'z' ; ++j) if (li->contains[j]) ++n ;
	OneStat *s = li->stats = new0 (n+1, OneStat) ;
	for (j = 'A' ; j <= 'z' ; ++j)
	  if (li->contains[j])
	    { s->type = j ;
	      if (vf->info[j]->listEltSize) s->isList = true ;
	      ++s ;
	    }
      }

  // finally set isFirst
  for (i = 'A' ; i <= 'z' ; ++i)
    if (vf->info[i])
      vf->info[i]->isFirst = true ;
}

static inline void setCodecBuffer (OneInfo *vi)
{
  vi->bufSize = vcMaxSerialSize() + 1; // +1 for added but unused 0-terminator
  vi->buffer  = new (vi->bufSize, void);
}

static OneFile *oneFileCreate (OneSchema **vsp, const char *type)
{ // searches through the linked list of vs to find type, either as primary or a secondary
  // if found fills and returns vf, else returns 0
  
  int         i, j ;
  OneFile    *vf = new0 (1, OneFile) ;
  char       *secondary = 0 ;
  OneSchema  *vs = *vsp ;

  { static bool isFirst = true ;
    if (isFirst)
      { if (sizeof(I64) != 8) die ("ONElib compile error: sizeof(I64) = %d != 8", sizeof(I64)) ;
        if (sizeof(I32) != 4) die ("ONElib compile error: sizeof(I32) = %d != 4", sizeof(I32)) ;
        if (sizeof(I16) != 2) die ("ONElib compile error: sizeof(I16) = %d != 2", sizeof(I16)) ;
        if (sizeof(I8) != 1) die ("ONElib compile error: sizeof(I8) = %d != 1", sizeof(I8)) ;
        isFirst = false ;
      }
  }
  
  // transfer header info
  for (i = 0 ; i < 128 ; ++i)
    if (vs->info[i]) vf->info[i] = infoDeepCopy (vs->info[i]) ;
  
  // find type in schema 
  while ((vs = vs->nxt))
    if (!strcmp (type, vs->primary))
      break ;
    else if (vs->nSecondary)
      { for (j = 0 ; j < vs->nSecondary ; ++j)
	  if (!strcmp (type, vs->secondary[j])) break ;
	if (j < vs->nSecondary) { secondary = vs->secondary[j] ; break ; }
      }
  if (!vs)
    { oneFileDestroy (vf) ;
      return 0 ; // failed to find a match
    }
  
  // transfer info from matched schema
  for (i = 0 ; i < 128 ; ++i)
    if (vs->info[i]) vf->info[i] = infoDeepCopy (vs->info[i]) ;

  initialiseStats (vf) ; // builds info->stats records once info is complete
  
  // build binaryTypeUnpack[]
  for (i = 0 ; i < 128 ; ++i)
    if (vf->info[i] && vf->info[i]->binaryTypePack)
      { U8 x = vf->info[i]->binaryTypePack ;
	vf->binaryTypeUnpack[x] = i ;
	vf->binaryTypeUnpack[x+1] = i ;
      }
  
  // set other information
  vf->fileType  = new (strlen(vs->primary)+1, char);
  strcpy (vf->fileType, vs->primary) ;
  if (secondary)
    { vf->subType  = new (strlen(secondary)+1, char);
      strcpy (vf->subType, secondary) ;
    }
  vf->nFieldMax = vs->nFieldMax ;
  vf->field = new (vf->nFieldMax, OneField) ;
  vf->nDefn = vs->nDefn ;
  memcpy (vf->defnOrder, vs->defnOrder, 128*sizeof(int)) ;
  for (i = 0 ; i < 128 ; ++i)
    if (vs->defnComment[i]) vf->defnComment[i] = strdup (vs->defnComment[i]) ;

  // setup for compression

  vf->codecTrainingSize = 100000;
  setCodecBuffer (vf->info[';']) ;
  
  // determine endian of machine
  { int   t = 1;
    char *b = (char *) (&t);
    vf->isBig = (*b == 0);
  }

  *vsp = vs ;
  return vf ;
}

static void provRefDefCleanup (OneFile *vf)
{ int n ;

  if (vf->provenance)
    { OneProvenance *p = vf->provenance ;
      for (n = vf->info['!']->accum.count ; n-- ; p++)
	{ free (p->program) ;
	  free (p->version) ;
	  free (p->command) ;
	  free (p->date) ;
	}
      free (vf->provenance) ;
    }
  if (vf->reference)
    { OneReference *r = vf->reference ;
      for (n = vf->info['<']->accum.count ; n-- ; r++) 
	free (r->filename) ;
      free (vf->reference) ;
    }
  if (vf->deferred)
    { OneReference *r = vf->deferred ;
      for (n = vf->info['>']->accum.count ; n-- ; r++) 
	free (r->filename) ;
      free (vf->deferred) ;
    }
}

static void oneFileCleanupSlaves (OneFile *vf)
{ int      i, j;
  OneInfo *li, *lx;

  for (i = 0; i < 128 ; i++)
    { lx = vf->info[i];
      if (lx != NULL)
	{ for (j = 1; j < vf->share; j++)
	    { li = vf[j].info[i];
	      if (li != lx) // the index OneInfos are shared
		{ if (li->listCodec == lx->listCodec) li->listCodec  = NULL;
		  if (li->index == lx->index) li->index = NULL ;
		  infoDestroy(li);
		}
	    }
	}
    }

  for (j = 1; j < vf->share; j++)
    { provRefDefCleanup (&vf[j]) ;
      if (vf[j].codecBuf   != NULL) free (vf[j].codecBuf);
      if (vf[j].f          != NULL) fclose (vf[j].f);
    }
}

static void oneFileDestroy (OneFile *vf)
{ int      i, j;

  if (vf->share)
    oneFileCleanupSlaves (vf) ;

  provRefDefCleanup (vf) ;
  if (vf->codecBuf != NULL) free (vf->codecBuf);
  if (vf->f != NULL && vf->f != stdout) fclose (vf->f);

  for (i = 0; i < 128 ; i++)
    if (vf->info[i] != NULL)
      infoDestroy (vf->info[i]);

  if (vf->field) free (vf->field) ;

  if (vf->headerText)
    { OneHeaderText *t = vf->headerText ;
      while (t)
	{ free (t->text) ;
	  { OneHeaderText *nxt = t->nxt ; free (t) ; t = nxt ; }
	}
    }

  for (j = 0 ; j < vf->nDefn ; ++j)
    if (vf->defnComment[j]) free (vf->defnComment[j]) ;

  free(vf->fileName);
  free(vf->fileType);
  free(vf->subType);

  free (vf) ;
}

/***********************************************************************************
 *
 *    ASCII PARSING UTILITIES: error reporting, lexical level
 *
 **********************************************************************************/

static void parseDie (OneFile *vf, char *format, ...)
{ va_list args;

  fprintf (stderr, "OneFile parse error: ");

  va_start (args, format);
  vfprintf (stderr, format, args);
  va_end (args);

  vf->lineBuf[vf->linePos] = '\0';
  fprintf (stderr, ", line %lld: %s\n", vf->line, vf->lineBuf);

  exit (1);
}

static inline char vfGetc(OneFile *vf)
{ char c = getc(vf->f);
  if (vf->linePos < 127)
    vf->lineBuf[vf->linePos++] = c;
  return c;
}

static inline void eatWhite (OneFile *vf)
{ char x = vfGetc(vf);
  if (x == ' ') // 200414: removed option to have tab instead of space
    return;
  parseDie (vf, "failed to find expected space separation character lineType %c", vf->lineType);
}

static inline char readChar(OneFile *vf)
{ eatWhite(vf);
  return vfGetc(vf);
}

static inline char *readBuf(OneFile *vf)
{ char x, *cp, *endBuf;

  eatWhite (vf);
  endBuf = vf->numberBuf + 32;
  for (cp = vf->numberBuf; cp < endBuf ; cp++)
    { x = vfGetc(vf);
      if (isspace(x) || x == '\0' || x == EOF)
        break;
      *cp = x;
    }
  if (cp >= endBuf)
    { cp[-1] = 0;
      parseDie (vf, "overlong item %s", vf->numberBuf);
    }
  else
    { ungetc (x, vf->f);
      vf->linePos -= 1;
      *cp = 0;
    }
  return vf->numberBuf;
}

static inline I64 readInt(OneFile *vf)
{ char *e, *b;
  I64   x;

  b = readBuf(vf);
  x = strtoll(b, &e, 10);
  if (e == b)
    parseDie (vf, "empty int field");
  if (*e != '\0')
    parseDie (vf, "bad int");
  return x;
}

static inline double readReal(OneFile *vf)
{ char  *e, *b;
  double x;

  b = readBuf(vf);
  x = strtod (b, &e);
  if (e == b)
    parseDie (vf, "empty real field");
  if (*e != '\0')
    parseDie (vf, "bad real");
  return (x);
}

static inline void readString(OneFile *vf, char *buf, I64 n)
{ eatWhite (vf);
  if (vf->isCheckString)
    { char *cp = buf;
      --cp;
      while (n-- && (*++cp = vfGetc (vf)))
        if (*cp == '\n' || *cp == EOF)
      break;
      if (++n)
        parseDie (vf, "line too short %d", buf);
      *++cp = 0;
    }
  else
    { if ((I64) fread (buf, 1, n, vf->f) != n)
	die ("ONE parse error: failed to read %d byte string", n);
      buf[n] = 0 ;
    }
}

static inline void readFlush (OneFile *vf) // reads to the end of the line and stores as comment
{ char       x;
  int        n = 0;
  OneInfo   *li = vf->info['/'] ;

  // check the first character - if it is newline then done
  x = getc (vf->f) ; 
  if (x == '\n')
    return ;
  else if (x != ' ')
    parseDie (vf, "comment not separated by a space") ;

  // else the remainder of the line is a comment
  if (!li->bufSize)
    { li->bufSize = 1024 ;
      li->buffer = new (li->bufSize, char) ;
    }
  while ((x = getc (vf->f)) && x != '\n')
    if (x == EOF)
      parseDie (vf, "premature end of file");
    else
      { if ((n+1) >= li->bufSize)
	  { char *s = new (2*li->bufSize, char) ;
	    memcpy (s, li->buffer, li->bufSize) ;
	    free (li->buffer) ;
	    li->buffer = s ;
	    li->bufSize *= 2 ;
	  }
	((char*)li->buffer)[n] = x ;
	++n ;
      }
  ((char*)li->buffer)[n] = 0 ; // string terminator
}


/***********************************************************************************
 *
 *    LIST BUFFER & COUNT MANAGEMENT: error reporting, lexical level
 *
 **********************************************************************************/

  //  Ensure line type t buffer can handles size+nStrings, and accumulate total and max

static inline void updateTotalAndBuffer (OneFile *vf, char t, I64 size, I64 nStrings)
{ OneInfo *li;

  li = vf->info[(int) t];
  li->accum.total += size;
  if (size > li->accum.max)
    li->accum.max = size;
  size += nStrings;             // need to allocate space for terminal 0s
  if ( ! li->isUserBuf && size > li->bufSize)   // expand buffer
    { if (li->buffer != NULL) free (li->buffer);
      li->bufSize = size + 0x10000 ;
      li->buffer  = new (li->bufSize*li->listEltSize, void);
    }
}

/***********************************************************************************
 *
 *   BINARY INT LIST COMPACTION & UNCOMPACTION
 *
 **********************************************************************************/

static char *compactIntList (OneFile *vf, OneInfo *li, I64 len, char *buf, int *usedBytes)
{ char *y;
  int   d, k;
  I64   z, i, mask, *ibuf;
  
  if (buf != li->buffer && !li->isUserBuf) // copy into li->buffer so can corrupt it
    { if ((I64) (li->bufSize) < len)
	{ if (li->buffer != NULL)
	    free (li->buffer);
	  li->bufSize = len + 1;
	  li->buffer = new (li->bufSize * sizeof(I64), void);
	}
      memcpy (li->buffer, buf, len*sizeof(I64)) ;
      buf = li->buffer ;
    }

  ibuf = (I64 *) buf;

  for (i = len-1; i > 0; i--)  // convert to differences - often a big win, else harmless
    ibuf[i] -= ibuf[i-1];
  
  mask = 0;                    // find how many top bytes can be skipped
  for (i = 1; i < len; i++)
    if (ibuf[i] >= 0) 
      mask |= ibuf[i];
    else
      mask |= -(ibuf[i]+1);

  k = sizeof(I64) ;
  mask >>= 7;
  for (d = 1; d < k; d++)
    { if (mask == 0)
        break;
      mask >>= 8;
    }
  *usedBytes = d ;

  z = k - d;   // number of 0 bytes
  if (z == 0) return (char*)&ibuf[1] ;

  y = li->buffer ;
  buf += sizeof(I64) ; --len ; // don't record the first element of buf, which is not a diff
  if (vf->isBig)     // copy d bytes per I64, ignoring z before or after depending on isBig
    while (len--)
      { buf += z;
        for (k = 0; k < d; k++)
          *y++ = *buf++;
      }
  else
    while (len--)
      { for (k = 0; k < d; k++)
          *y++ = *buf++;
        buf += z;
      }
 
  return li->buffer ;
}

static void decompactIntList (OneFile *vf, I64 len, char *buf, int usedBytes)
{ int   d, z, k;
  char *s, *t;

  z = sizeof(I64) - usedBytes ;
  
  if (z > 0)                          // decompacts in place
    { buf += sizeof(I64) ; --len ;    // don't decompact 0th element
      d = usedBytes;
      s = buf + d*len;
      t = s + z*len; 
      if (vf->isBig)
	while (s > buf)
          { for (k = 0; k < d; k++)
              *--t = *--s;
            if (*s & 0x80)
              for (k = 0; k < z; k++)
                *--t = 0xff;
            else
              for (k = 0; k < z; k++)
                *--t = 0x0;
          }
      else
	while (s > buf)
          { if (s[-1] & 0x80)
              for (k = 0; k < z; k++)
                *--t = 0xff;
            else
              for (k = 0; k < z; k++)
                *--t = 0;
            for (k = 0; k < d; k++)
              *--t = *--s;
          }
      buf -= sizeof(I64) ; ++len ;
    }
  
  { I64 i, *x = (I64 *) buf;             // revert differencing
    for (i = 1; i < len; i++)
      x[i] += x[i-1];
  }
}

// read and write compressed fields

static inline int writeCompressedFields (FILE *f, OneField *field, OneInfo *li)
{
  int i, n = 0 ;
  
  for (i = 0 ; i < li->nField ; ++i)
    switch (li->fieldType[i])
      {
      case oneREAL: fwrite (&field[i].r, 8, 1, f) ; n += 8 ; break ;
      case oneCHAR: putc (field[i].c, f) ; ++n ; break ;
      default: // includes INT and all the LISTs, which store their length in field as an INT
	n += ltfWrite (field[i].i, f) ;
      }

  return n ;
}

static inline void readCompressedFields (FILE *f, OneField *field, OneInfo *li)
{
  int i ;
  
  for (i = 0 ; i < li->nField ; ++i)
    switch (li->fieldType[i])
      {
      case oneREAL:
	if (fread (&field[i].r, 8, 1, f) != 1) die ("failed to read a REAL") ;
	break ;
      case oneCHAR:
	field[i].c = fgetc (f) ;
	break ;
      default: // includes INT and all the LISTs, which store their length in field as an INT
	field[i].i = ltfRead (f) ;
      }
}

/***********************************************************************************
 *
 *  ONE_READ_LINE:
 *      Reads the next line and returns false at end of file or on error. The line is
 *      parsed according to its linetype and contents accessed by macros that follow.
 *      The top bit of the first character determines whether the line is binary or ascii
 *
 **********************************************************************************/

  //  Read a string list, first into new allocs, then into sized line buffer.
  //    Annoyingly inefficient, but we don't use it very much.

static void readStringList(OneFile *vf, char t, I64 len)
{ int    j;
  I64    totLen, sLen;
  char **string, *buf;

  totLen = 0;
  string = new (len, char *);
  for (j = 0; j < len ; ++j)
    { sLen = readInt (vf);
      totLen += sLen;
      string[j] = new (sLen+1, char);
      readString (vf, string[j], sLen);
    }

  updateTotalAndBuffer (vf, t, totLen, len);

  buf = (char *) vf->info[(int) t]->buffer;
  for (j = 0; j < len ; ++j)
    { strcpy (buf, string[j]);
      buf += strlen(buf) + 1;
      free (string[j]);
    }
  free (string);
}

bool addProvenance(OneFile *vf, OneProvenance *from, int n) ; // need forward declaration

char oneReadLine (OneFile *vf)
{ bool      isAscii;
  U8        x;
  char      t;
  OneInfo  *li;

  assert (!vf->isWrite) ;
  assert (!vf->isFinal) ;

  vf->linePos = 0;                 // must come before first vfGetc()
  x = vfGetc (vf);                 // read first char
  if (feof (vf->f) || x == '\n')   // blank line (x=='\n') is end of records marker before footer
    { vf->lineType = 0 ;           // additional marker of end of file
      return 0;
    }

  vf->line += 1;      // otherwise assume this is a good line, and die if not
  if (x & 0x80)
    { isAscii = false;
      t = vf->binaryTypeUnpack[x];
    }
  else
    { isAscii = true;
      t = x;
    }
  vf->lineType = t;

  // if (!isBootStrap) fprintf (stderr, "reading line %d type %c\n", (int)vf->line, t) ;

  li = vf->info[(int) t];
  if (li == NULL)
    parseDie (vf, "unknown line type %c (%d was %d) line %d", t, t, x, (int)vf->line);
  if (li->accum.count >= 0) // after goto set to -1 for unindexed linetypes - can't know the count
    li->accum.count += 1 ;  // includes update of indexed type counts

  if (vf->info['/']->bufSize) // clear the comment buffer
    *(char*)(vf->info['/']->buffer) = 0 ;

  vf->nBits = 0 ;        // will use for any compressed data read in
  
  if (isAscii)           // read field by field according to ascii spec
    { int     i, j;
      I64    *ilst, len;
      double *rlst;

      for (i = 0; i < li->nField; i++)
        switch (li->fieldType[i])
	  {
	  case oneINT:
            vf->field[i].i = readInt (vf);
	    //	    printf ("  field %d int %d\n", i, (int)oneInt(vf,i)) ; 
	    break;
          case oneREAL:
            vf->field[i].r = readReal (vf);
            break;
          case oneCHAR:
            vf->field[i].c = readChar (vf);
	    //	    printf ("  field %d char %c\n", i, (int)oneChar(vf,i)) ;
            break;
	  case oneSTRING:
	  case oneDNA:
            len = readInt (vf);
            vf->field[i].len = len;
            updateTotalAndBuffer (vf, t, len, 1);
            readString (vf, (char*) li->buffer, len);
            break;
          case oneINT_LIST:
            len = readInt (vf);
            vf->field[i].len = len;
            updateTotalAndBuffer (vf, t, len, 0);
            ilst = (I64 *) li->buffer;
            for (j = 0; j < len; ++j)
              ilst[j] = readInt(vf);
            break;
          case oneREAL_LIST:
            len = readInt (vf);
            vf->field[i].len = len;
            updateTotalAndBuffer (vf, t, len, 0);
            rlst = (double *) li->buffer;
            for (j = 0; j < len; ++j)
              rlst[j] = readReal (vf);
            break;
          case oneSTRING_LIST: // STRING_LIST - inefficient for now - also used for binary
            len = readInt (vf);
            vf->field[i].len = len;
	    //	    printf ("  field %d string list len %d\n", i, (int)oneLen(vf)) ;
            readStringList (vf, t, len);
            break;
	  }
      readFlush (vf);
    }

  else        // binary - block read fields and list, potentially compressed
    { 
      // read the fields

      if (li->nField > 0)
	readCompressedFields (vf->f, vf->field, li) ;

      // read the list if there is one

      if (li->listEltSize > 0)
        { I64 listLen = oneLen(vf);

          if (listLen > 0)
            { li->accum.total += listLen;
	      if (listLen > li->accum.max)
		li->accum.max = listLen;

	      if (li->fieldType[li->listField] == oneINT_LIST)
		{ *(I64*)li->buffer = ltfRead (vf->f) ;
		  if (listLen == 1) goto doneLine ;
		  vf->intListBytes = getc(vf->f) ;
		}

	      if (li->fieldType[li->listField] == oneSTRING_LIST) // handle as ASCII
                readStringList (vf, t, listLen);
              else if (x & 0x1)    				  // list is compressed
                { vf->nBits = ltfRead (vf->f) ;
		  size_t bytes = (vf->nBits+7) >> 3 ;
		  if (bytes > (size_t) vf->codecBufSize)
		    { if (vf->codecBuf) free (vf->codecBuf) ;
		      vf->codecBufSize = bytes + 1 ;
		      vf->codecBuf = new (vf->codecBufSize, void) ;
		    }
                  if (fread (vf->codecBuf, bytes, 1, vf->f) != 1)
                    die ("ONE read error: fail to read compressed list");
                }
              else if (li->fieldType[li->listField] == oneINT_LIST)
                { I64 listSize  = (listLen-1) * vf->intListBytes ;
                  if ((I64) fread (&(((I64*)li->buffer)[1]), 1, listSize, vf->f) != listSize)
                    die ("ONE read error: failed to read list size %lld", listSize);
		  decompactIntList (vf, listLen, li->buffer, vf->intListBytes);
                }
	      else
                { I64 listSize  = listLen * li->listEltSize ;
                  if ((I64) fread (li->buffer, 1, listSize, vf->f) != listSize)
                    die ("ONE read error: failed to read list size %lld", listSize);
                }
            }

          if (li->fieldType[li->listField] == oneSTRING)
            ((char *) li->buffer)[listLen] = '\0'; // 0 terminate
        }

    doneLine:

      { U8 peek = getc(vf->f) ; // check if next line is a comment - if so then read it
	ungetc(peek, vf->f) ;
	if (peek & 0x80)
	  peek = vf->binaryTypeUnpack[peek];
	if (peek == '/') // a comment
	  { OneField keepField0 = vf->field[0] ;
	    I64 keepNbits = vf->nBits ; // will be reset in readLine
	    oneReadLine (vf) ; // read comment line into vf->info['/']->buffer
	    vf->lineType = t ;
	    vf->field[0] = keepField0 ;
	    vf->nBits = keepNbits ;
	  }
      }
    }

  return t;
}

char *oneReadComment (OneFile *vf)
{ char *comment = (char*)(vf->info['/']->buffer) ;

  if (comment && *comment != 0)
    return comment ;
  else
    return 0 ;
}

void *_oneList (OneFile *vf)
{
  OneInfo *li = vf->info[(int) vf->lineType] ;

  if (vf->nBits)
    { if (li->fieldType[li->listField] == oneINT_LIST) // first elt is already in buffer
	{ vcDecode (li->listCodec, vf->nBits, vf->codecBuf, (char*)&(((I64*)li->buffer)[1])) ;
	  decompactIntList (vf, oneLen(vf), li->buffer, vf->intListBytes) ;
	}
      else
	vcDecode (li->listCodec, vf->nBits, vf->codecBuf, li->buffer) ;
      vf->nBits = 0 ; // so we don't do it again
    }
  
  return li->buffer ;
}

void *_oneCompressedList (OneFile *vf)
{
  OneInfo *li = vf->info[(int) vf->lineType] ;

  if (!vf->nBits && oneLen(vf) > 0)      // need to compress
    vf->nBits = vcEncode (li->listCodec, oneLen(vf),
			  vf->info[(int) vf->lineType]->buffer, vf->codecBuf);

  return (void*) vf->codecBuf ;
}

OneFile *readThreadMake (OneFile *vfOld, OneSchema *vs0, FILE **files)
{
  int   i, j ;
  int   nthreads = vfOld->share ;
  off_t startOff ;

  OneFile *vf = new0 (nthreads, OneFile) ;
  vf[0] = *vfOld ;
  free (vfOld) ; // NB free() not oneFileDestroy because don't want deep destroy

  startOff = ftello (vf->f) ;
  for (i = 1; i < nthreads; i++)
    { OneSchema *vs = vs0 ; // needed because vs will have changed to map to the relevant page
      OneFile   *v = oneFileCreate(&vs, vf->fileType); // need to do this after header is read
      vf[i] = *v ;
      free (v) ;
      v = vf+i;
      
      v->share = -i ; // so this slave knows its own identity

      v->f = files[i] ;
      if (fseeko (v->f, startOff, SEEK_SET) != 0)
	die ("ONE file error: can't seek to start of data in thread file");
      
      for (j = 0; j < 128; j++)
	{ OneInfo *li = v->info[j];
	  if (li != NULL)
	    { OneInfo *l0 = vf->info[j];
	      li->given = l0->given; // copy the given data
	      li->accum = l0->accum; // copy the accum data (needed for reopen)
	      if (li->listCodec) vcDestroy (li->listCodec) ; // share the codec
	      li->listCodec  = l0->listCodec;
	      if (li->listEltSize > 0) // need a private buffer
		{ li->bufSize = l0->bufSize;
		  if (li->buffer) free (li->buffer) ;
		  li->buffer  = new (l0->bufSize*l0->listEltSize, void);
		}
	      if (l0->isObject) li->isObject = true ;
	      if (l0->index) // share the index
		{ if (li->index) free(li->index) ;
		  li->index = l0->index ;
		}
	      if (l0->stats) // copy the group data
		{ OneStat *s0 = l0->stats, *s = li->stats ;
		  for ( ; s0->type ; ++s0, ++s) *s = *s0 ;
		}
	    }
	}

      v->codecBufSize = vf->codecBufSize;
      if (v->codecBuf) free (v->codecBuf) ;
      v->codecBuf = new (v->codecBufSize, void); // need a private codec buffer

      if (vf->subType != NULL)
	v->subType = strdup (vf->subType) ;
      else
	v->subType = NULL;
    }

  return vf ;
}

/***********************************************************************************
 *
 *   ONE_FILE_OPEN_READ:
 *     Opens file for reading and reads all lines of header if it exists.
 *     If there is no header then fileType must be given (i.e. non-zero),
 *     otherwise if fileType is non-zero then it must match the type in the header.
 *     If the header contains a $-line then it is binary and the routine reads
 *     the footer that includes decompressors and indexes.
 *
 **********************************************************************************/

OneFile *oneFileOpenRead (const char *path, OneSchema *vsArg, const char *fileType, int nthreads)
{
  OneFile   *vf ;
  off_t      startOff = 0, footOff;
  OneSchema *vsFile ;                // will be used to build schema from file
  OneSchema *vs0 ;                   // needed when making slave thread entries
  bool       isBareFile = false ;

  assert (fileType == NULL || strlen(fileType) > 0) ;

  // first open the file, read first header line if it exists, and create the OneFile object
  
  { FILE *f ;
    int   curLine = 0 ;
    U8    c ;

    if (strcmp (path, "-") == 0)
      f = stdin;
    else
      { f = fopen (path, "r");
	if (f == NULL)
	  return NULL;
      }
    
#define OPEN_ERROR1(x) \
    { snprintf (errorString, 1024, "ONEcode file open error %s: %s\n", path, x) ; \
      fclose(f) ; return NULL; }
#define OPEN_ERROR3(x,y,z) \
    { int nChar = snprintf (errorString, 1024, "ONEcode file open error %s: ", path) ; \
    nChar += snprintf (errorString+nChar, 1024-nChar, x,y,z) ; \
    snprintf (errorString+nChar, 1024-nChar, "\n") ; \
    fclose(f) ; return NULL ; }
    
    c = getc(f);
    if (feof(f))
      OPEN_ERROR1("file is empty") ;

    if (c == '1')
      { int  major, minor, slen;
	char *primaryName ;

	if (fscanf (f, " %d", &slen) != 1)
	  OPEN_ERROR1("line 1: failed to read type name length") ;
	if (slen == 0)
	  OPEN_ERROR1("line 1: type name is empty string") ;
        primaryName = new0 (slen+1, char);
	if (fscanf (f, " %s %d %d", primaryName, &major, &minor) != 3)
	  OPEN_ERROR1("line 1: failed to read remainder of line") ;
	while (getc(f) != '\n')
	  if (feof(f))
	    OPEN_ERROR1("end of file before end of line 1") ;
	++curLine ;
	if (major != MAJOR)
	  OPEN_ERROR3("major version file %d != code %d", major, MAJOR) ;
	if (minor > MINOR)
	  OPEN_ERROR3("minor version file %d > code %d", minor, MINOR) ;
	vs0 = vsFile = oneSchemaCreateDynamic (primaryName, 0) ; // create a shell schema
	vf = oneFileCreate (&vsFile, primaryName)  ;
	free (primaryName) ;
      }
    else
      { ungetc (c, f) ;
	isBareFile = true ;
	if (!fileType || !vsArg)
	  OPEN_ERROR1("attempting to open a bare oneFile without giving the type and/or schema") ;
	vs0 = vsArg ;
	vf = oneFileCreate (&vsArg, fileType) ;
	if (!vf)
	  OPEN_ERROR1("failed to find given type in given schema") ;
      }
    
    vf->f = f;
    vf->line = curLine;
    vf->fileName = strdup(path) ;
  }

  // read header and (optionally) footer
  // recognise end of header by peeking at the first char to check if alphabetic 
 
  vf->isCheckString = true;   // always check strings while reading header
  I64 maxIndexSize = 0 ;      // needed to make buffer space for reading in indices
  while (true)
    { U8 peek = getc(vf->f);

      if (feof(vf->f))       // loop exit at end of file
        break;
      ungetc(peek, vf->f);

      if (peek & 0x80)
        peek = vf->binaryTypeUnpack[peek];

      if (peek == '&')
	{ vf->info['&']->bufSize = maxIndexSize ; // make the buffer to read in the indexes
	  vf->info['&']->buffer  = new (maxIndexSize, I64) ;
	}

      if (isalpha(peek) || peek == '\n')  // '\n' to check for end of binary file, i.e. empty file
        break;    // loop exit at standard data line

      if (isBareFile) // can't have any special header lines
	{ snprintf (errorString, 1024,
		   "ONEcode file open error %s: if header exists it must begin with '1' line\n",
		   path) ;
	  oneFileDestroy (vf) ;
	  return 0 ;
	}

      oneReadLine(vf);  // can't fail because we checked file eof already

      switch (vf->lineType)
	{
	case '1':
          parseDie(vf, "1 should be first line in header");
          break;

        case '2':
	  vf->subType = strdup (oneString(vf)) ;
	  break;

	case '.': // blank line for spacing and header text
	  { char *text = oneReadComment (vf) ;
	    if (text)
	      {	OneHeaderText *t ;
		if (vf->headerText)
		  { t = vf->headerText ; while (t->nxt) t = t->nxt ;
		    t->nxt = new0 (1, OneHeaderText) ; t = t->nxt ;
		  }
		else
		  t = vf->headerText = new0 (1, OneHeaderText) ;
		t->text = strdup (text) ;
	      }
	    break ;
	  }

	case '~': // schema definition line
	  { char t = oneChar(vf,1) ;
	    if (!(t >= 'A' && t <= 'Z') && !(t >= 'a' && t <= 'z'))
	      die ("type symbol %c in schema definition line %d is not a letter", t, vf->line) ;
	    if (oneChar(vf,0) == 'G')
	      { schemaAddGroup (vsFile, t) ;
		if (oneReadComment (vf)) vf->defnComment[vf->nDefn] = strdup (oneReadComment (vf)) ;
		vf->defnOrder[vf->nDefn++] = t | 0x80 ; // definition order
	      }
	    else
	      { int oldMax = vf->nFieldMax ;
		schemaAddInfoFromLine (vsFile, vf, t, oneChar(vf,0)) ;
		if (oneReadComment (vf)) vf->defnComment[vf->nDefn] = strdup (oneReadComment (vf)) ;
		OneInfo *vi = vsFile->info[(int)t] ;
		vf->info[(int)t] = infoDeepCopy (vi) ;
		vf->defnOrder[vf->nDefn++] = t ; // definition order
		if (vi->binaryTypePack)
		  { U8 x = vi->binaryTypePack ;
		    vf->binaryTypeUnpack[x] = t ;
		    vf->binaryTypeUnpack[x+1] = t ;
		  }
		if (vsFile->nFieldMax > oldMax)
		  { free (vf->field) ;
		    vf->nFieldMax = vsFile->nFieldMax ;
		    vf->field = new (vf->nFieldMax, OneField) ;
		  }
	      }
	  }
	  break ;

        case '#': // count information
        case '@':
        case '+':
        case '%':
          { char      c = oneChar(vf,0);
            OneInfo *li = vf->info[(int) c] ;
            if (li == NULL) parseDie (vf, "unknown line type %c", c);
            switch (vf->lineType)
            { case '#':
                li->given.count = oneInt(vf,1);
		if (vf->isBinary && li && li->isObject)  // allocate space for indices
		  { li->indexSize = li->given.count + 1 ; // +1 because 1..n
		    li->index = new (li->indexSize, I64) ;
		    if (li->indexSize > maxIndexSize) maxIndexSize = li->indexSize ;
		  }
                break;
              case '@':
                li->given.max = oneInt(vf,1);
                li->bufSize = li->given.max + 1; // allow for string terminators
                li->buffer = new (li->bufSize*li->listEltSize, void);
		break;
              case '+':
                li->given.total = oneInt(vf,1);
                break;
              case '%':
                { if (!li->isObject) parseDie (vf, "% on a non-object type %c", c) ;
		  if (!li->stats) initialiseStats (vf) ;
		  int j = oneChar(vf,2);
		  OneStat *s ;
		  for (s = li->stats ; s->type && s->type != j ; ++s)
		  if (!s->type) parseDie (vf, "unknown line type %c", j);
		  c = oneChar(vf,1);
		  if (c == '#')
		    s->maxCount = oneInt(vf,3);
		  else if (c == '+')
		    s->maxTotal = oneInt(vf,3);
		  else
		    parseDie (vf, "unrecognised symbol %c", c);
		}
		break;
            } /*  */
	  }
	  break;

        case '!':     // NB need to copy the strings
          { OneProvenance p ;
	    p.program = oneString(vf) ;
	    p.version = p.program + strlen(p.program) + 1 ;
	    p.command = p.version + strlen(p.version) + 1 ;
	    p.date    = p.command + strlen(p.command) + 1 ;
            vf->info['!']->accum.count -= 1; // to avoid double counting
            addProvenance (vf, &p, 1) ;
          }
	  break;

        case '<':
          vf->info['<']->accum.count -= 1; // to avoid double counting
          oneAddReference (vf, oneString(vf), oneInt(vf,1));
          break;

        case '>':
          vf->info['>']->accum.count -= 1; // to avoid double counting
          oneAddDeferred (vf, oneString(vf));
          break;

        // Below here are binary file header types - requires given.count/given.max first

        case '$':  // read footer - goto end, find offset to start of footer and go there
          if (oneInt(vf,0) != vf->isBig)
            die ("ONE file error: endian mismatch - convert file to ascii");
          vf->isBinary = true;

          startOff = ftello (vf->f);
          if (fseek (vf->f, -sizeof(off_t), SEEK_END) != 0)
            die ("ONE file error: can't seek to final line");

          if (fread (&footOff, sizeof(off_t), 1, vf->f) != 1)
            die ("ONE file error: can't read footer offset");

          if (fseeko (vf->f, footOff, SEEK_SET) != 0)
            die ("ONE file error: can't seek to start of footer");

          break;

        case '^':    // end of footer - return to where we jumped from header
          if (fseeko (vf->f, startOff, SEEK_SET) != 0)
            die ("ONE file error: can't seek back");
          break;

        case '&': // read index
	  { char c = oneChar(vf,0) ;
	    OneInfo *li = vf->info[(int)c] ;
	    assert (li->indexSize == oneLen(vf)) ;
	    assert (li->index) ;
	    memcpy (li->index, oneIntList(vf), oneLen(vf)*sizeof(I64)) ; // space allocated above
	  }
          break;

        case ';':
          vf->info[(int) oneChar(vf,0)]->listCodec = vcDeserialize (oneString(vf));
          break;

        default:
          parseDie (vf, "unknown header line type %c", vf->lineType);
          break;
      }
    }
  vf->isCheckString = false;   // user can set this back to true if they wish

  if (!isBareFile && vsArg && !oneFileCheckSchema (vf, vsArg, false)) // check schema intersection
    { snprintf (errorString, 1024, "ONEcode file open error %s: schema mismatch to code requirement\n", path) ;
      oneFileDestroy (vf) ;
      return NULL ;
    }

  if (!vf->info[0] || !vf->info[0]->isObject) // will be true if initialiseStats() has been called
    initialiseStats (vf) ; // call here in case not called above for a % line
  
  // allocate codec buffer - always allocate enough to handle fields of all line types

  { I64 size = vf->nFieldMax * sizeof(OneField) ;
    int i ;
    
    for (i = 0; i < 128; ++i)
      if (vf->info[i])
	{ OneInfo *li = vf->info[i];
	  if (li->listCodec && size < li->given.max * li->listEltSize)
	    size = li->given.max * li->listEltSize;
	}
    if (size >= vf->codecBufSize)
      { if (vf->codecBuf) free (vf->codecBuf) ;
	vf->codecBufSize = size+1;
	vf->codecBuf     = new (vf->codecBufSize, void);  // add one for worst case codec usage
      }
  }

  // if parallel, allocate a OneFile array for parallel thread objects, switch vf to head of array

  if (nthreads > 1) // should we allow multiple threads for a bare file, which has no index?
    { int i ;
      FILE **files = new (nthreads, FILE*) ;

      if (strcmp (path, "-") == 0) die ("ONE error: parallel input incompatible with stdin as input");

      for (i = 1 ; i < nthreads ; ++i) files[i] = fopen (path, "r") ;
vf->share = nthreads ;
      vf = readThreadMake (vf, vs0, files) ;
      free (files) ;
    }

  if (!isBareFile)
    oneSchemaDestroy (vs0) ;
    
  return vf;
}

static void oneFinalize (OneFile *vf) ; // forward declaration

static OneSchema *oneSchema (OneFile *vf)
{
  OneSchema *vs = oneSchemaCreateDynamic (vf->fileType, vf->subType) ;
  vs->nFieldMax = vf->nFieldMax ;
  int i ;
  for (i = 0 ; i < vf->nDefn ; ++i)
    { int t = vf->defnOrder[i] ;
      if (t & 0x80) // a group
	schemaAddGroup (vs, t & 0x7f) ;
      else
	vs->info[t] = infoDeepCopy (vf->info[t]) ;
      if (vf->defnComment[i]) vs->defnComment[i] = strdup (vf->defnComment[i]) ;
      vs->defnOrder[vs->nDefn++] = t ;
    }
    
  return vs ;
}

OneFile *oneFileReopenRead (OneFile *vf)
{
  if (!vf || !vf->isWrite) return 0 ;
  oneFinalize (vf) ;    // merges in data from any slaves, completes indices etc.
  oneFileCleanupSlaves (vf) ;
  vf->isFinal = false ; // so we can now read it again
  vf->isWrite = false ; // now it will be readonly
  oneGoto (vf, 0, 0) ; // go to start of data
  if (vf->share == 1)
    return vf ;
  else
    { OneSchema *vs0 = oneSchema (vf) ; // need this because of how oneFileCreate() works
      return readThreadMake (vf, vs0, vf->tempReadFiles) ; // use the cached file handles
      oneSchemaDestroy (vs0) ;
    }
}

/***********************************************************************************
 *
 *   ONE_USER_BUFFER / GOTO
 *
 **********************************************************************************/

  // This lets the user reassign the buffer that lists in a particular line type are read into.
  //   If this is not set, a default buffer is provided.  If buffer == NULL then the package
  //   reverts to the default buffer.  This routine can be called repeatedly.
  // NB the package doesn't check the size of a user supplied buffer - the user must allocate
  //   enough memory for all forthcoming list data.  For a binary file li->given.max+1 gives this.

void oneUserBuffer (OneFile *vf, char lineType, void *buffer)
{ OneInfo *li;

  li = vf->info[(int) lineType];
  if (buffer != NULL)
    { if ( ! li->isUserBuf && li->buffer != NULL)
        { free (li->buffer);
          li->bufSize = 0;
        }
      li->buffer    = buffer;
      li->isUserBuf = true;
    }
  else
    { if (li->isUserBuf)
        { li->bufSize = li->given.max + 1;
          li->buffer  = new (li->given.max*li->listEltSize, void);
        }
      li->isUserBuf = false;
    }
}

bool oneGoto (OneFile *of, char lineType, I64 i)
{
  OneInfo *li = of->info[(int)lineType] ;
  if (!li || !li->index || i < 0 || i > li->given.count) return false ;

  I64 byte = li->index[i] ;
  if (fseek (of->f, byte, SEEK_SET) != 0) return false ;

  li->accum.count = i ? i-1 : 0 ;

  int j, k ;
  for (k = 0 ; k < of->nDefn ; ++k)
    { j = of->defnOrder[k] ;
      if (!(j & 0x80) && j != lineType) // must set vj->accum.count
	{ OneInfo *lj = of->info[j] ;
	  if (i == 0) // start of data for all linetypes
	    lj->accum.count = 0 ;
	  else if (!lj->index) // we can't establish the location - disable count
	    lj->accum.count = -1 ;
	  else if (lj->index[lj->given.count] < byte) // after the start of the last object
	    lj->accum.count = lj->given.count ;
	  else // binary search
	    { int i0 = 0, i1 = lj->given.count ;
	      while (i1 > i0+1)
		{ i = (i1+i0)/2 ;
		  if (lj->index[i] < byte) i0 = i ;
		  else i1 = i ;
		}
	      lj->accum.count = i ;
	    }
	}
    }
  
  return true ;
}

I64 oneCountUntilNext (OneFile *of, char countType, char nextType)
// returns the number of countType object lines before the next nextType object line
// returns -1 on error, e.g. not reading a binary file, types are not object types
{
  if (of->isWrite || !of->isBinary) return -1 ;
  OneInfo *ci = of->info[(int)countType], *ni = of->info[(int)nextType] ;
  if (!ci || !ci->index || !ni || !ni->index) return -1 ;
  if (ni->accum.count == ni->given.count) return ci->given.count - ci->accum.count ;
  I64 nb = ni->index[ni->accum.count + 1] ;
  I64 ix = ci->accum.count + 1 ;
  while (ix <= ci->given.count && ci->index[ix] < nb) ++ix ;
  return ix - 1 ;
}

/***********************************************************************************
 *
 *   ONE_OPEN_WRITE_(NEW | FROM)
 *
 **********************************************************************************/

static inline void allocateIndices (OneFile *vf, int i, I64 size)
{
  OneInfo *li = vf->info[i] ;
  assert (li != NULL) ;

  if (li->isObject)
    { li->indexSize = size ;
      if (li->index) free(li->index) ;
      li->index = new (size, I64) ;
    }
}

OneFile *oneFileOpenWriteNew (const char *path, OneSchema *vs, const char *fileType,
                              bool isBinary, int nthreads)
{ OneFile   *vf ;
  FILE      *f, **tempReadFiles = 0 ;
  OneSchema *vs0 = vs ; // needed here because call to oneFileCreate changes vs
  char      *tempPath, *template ; // used for temporary files (thread files and if path is a dir)

  tempPath = new(strlen(path)+12, char) ;
  strcpy (tempPath, path) ;
  template = tempPath + strlen(tempPath) ;
	  
  if (strcmp (path, "-") == 0)
    f = stdout;
  else
    { struct stat status ;
      if (stat (path, &status) >= 0 && (status.st_mode & S_IFMT) == S_IFDIR) // a directory
	{ *template++ = '/' ;
	  strcpy (template, "oneXXXXXX") ;
	  int fd = mkstemp (tempPath) ;
	  if (fd == -1) return NULL ;
	  f = fdopen (fd, "w+") ;
	  if (nthreads > 1)
	    { tempReadFiles = new (nthreads, FILE*) ;
	      int i ;
	      for (i = 1 ; i < nthreads ; ++i)
		tempReadFiles[i] = fopen (tempPath, "r") ;
	    }
	  if (unlink(tempPath) < 0)
	    die ("ONEfile error: failed to unlink temporary file %s for parallel write", tempPath) ;
	}
      else
	{ f = fopen (path, "w");
	  if (f == NULL) return NULL ;
	}
    }

  vf = oneFileCreate (&vs, fileType) ;
  if (!vf) return NULL ;

  initialiseStats (vf) ;
  
  vf->f = f;
  vf->fileName = strdup (path) ;
  vf->isWrite  = true;
  vf->isBinary = isBinary;
  vf->isLastLineBinary = true; // we don't want to add a newline before the first true line
  
  vf->codecBufSize = vf->nFieldMax*sizeof(OneField) + 1;
  vf->codecBuf     = new (vf->codecBufSize, void);
  int ii ;
  for (ii = 0 ; ii < vf->nDefn ; ++ii)
    if (!(vf->defnOrder[ii] & 0x80))
      allocateIndices (vf, vf->defnOrder[ii], 0x10000) ;

  if (nthreads > 1)
    { OneFile *v, *vf0 = vf ;
      int      i ;
      char     name[100] ;

      vf->share = nthreads ;
      if (tempReadFiles) vf->tempReadFiles = tempReadFiles ;
      vf->fieldLock = mutexInit;
      vf->listLock  = mutexInit;
      vf = new (nthreads, OneFile);
      vf[0] = *vf0 ;
      free (vf0) ; // NB free() not oneFileDestroy because don't want deep destroy
      
      for (i = 1; i < nthreads; i++)
	{ vs = vs0 ; // needed because vs will have changed in prevous oneFileCreate call
	  v = oneFileCreate (&vs, fileType);

	  v->isWrite  = true;
	  v->isBinary = isBinary;
          v->isLastLineBinary = isBinary;
	  
	  v->codecBufSize = vf->codecBufSize;
	  v->codecBuf     = new (v->codecBufSize, void);
	  v->codecTrainingSize /= 3*nthreads;
	  for (ii = 0 ; ii < v->nDefn ; ++ii)
	    if (!(vf->defnOrder[ii] & 0x80))
	      allocateIndices (v, v->defnOrder[ii], 0x10000) ; // make separate indices for each thread

          v->share = -i; // this is the key mark for the i'th slave

	  strcpy (template, "oneXXXXXX") ;
	  int fd = mkstemp (tempPath) ;
	  if (fd == -1) return NULL ;
	  f = fdopen (fd, "w+") ;
          if (f == NULL)
            die ("ONEfile error: cannot create temporary file %s for parallel write", name) ;
	  if (unlink(tempPath) < 0)
	    die ("ONEfile error: failed to unlink temporary file %s for parallel write", name) ;
	  v->f = f ;

	  vf[i] = *v ;
	  free (v) ;
	}
    }

  return vf;
}

OneFile *oneFileOpenWriteFrom (const char *path, OneFile *vfIn, bool isBinary, int nthreads)
{
  // first build a schema from vfIn
  OneSchema *vs0 = oneSchemaCreateDynamic (vfIn->fileType, vfIn->subType) ;
  OneSchema *vs = vs0->nxt ; // this is the actual schema - vs0 is for the header

  int i, k;
  for (k = 0 ; k < vfIn->nDefn ; ++k)
    { i = vfIn->defnOrder[k] ;
      if (i & 0x80) schemaAddGroup (vs, (char)(i & 0x7f)) ;
      else
	{ OneInfo *li = vfIn->info[i] ;
	  if (li->isObject) schemaAddInfoFromArray (vs, li->nField, li->fieldType, (char)i, 'O') ;
	  else schemaAddInfoFromArray (vs, li->nField, li->fieldType, (char)i, 'D') ;
	}
      if (vfIn->defnComment[k]) vs->defnComment[k] = strdup (vfIn->defnComment[k]) ;
    }

  // use it to open the file
  OneFile *vf = oneFileOpenWriteNew (path, vs0, vfIn->subType ? vfIn->subType : vfIn->fileType,
				     isBinary, nthreads);
  oneSchemaDestroy (vs0) ;
  if (!vf) return NULL ;

  oneInheritProvenance (vf, vfIn);
  oneInheritReference  (vf, vfIn);
  oneInheritDeferred   (vf, vfIn);

  if (vfIn->headerText)
    { OneHeaderText *tin = vfIn->headerText ;
      OneHeaderText *t = new0 (1, OneHeaderText) ;
      vf->headerText = t ;
      while (tin)
	{ t->text = strdup(tin->text) ; tin = tin->nxt ;
	  if (tin) t = t->nxt = new0 (1, OneHeaderText) ;
	}
    }     
  
  // set info[]->given, and resize codecBuf accordingly
  I64 size = vf->codecBufSize;
  for (i = 0; i < 128 ; ++i)
    if (vf->info[i] && vfIn->info[i]->given.count)
      { OneInfo *vi = vf->info[i];
	vi->given = vfIn->info[i]->given ;
	if (vi->listCodec)
	  { I64 sz = vi->given.max * vi->listEltSize;
	    if (sz >= size)
	      size = sz+1;
	  }
	for (k = 0 ; k < nthreads ; ++k)
	  allocateIndices (&(vf[k]), i, vi->given.count+1) ; // resize the indices
      }
  if (size > vf->codecBufSize)
    for (i = 0 ; i < nthreads ; ++i)
      { OneFile *v = &(vf[i]) ;
	if (v->codecBuf) free (v->codecBuf) ;
	v->codecBufSize = size;
	v->codecBuf     = new (size, void);
      }

  return vf ;
}

bool oneFileCheckSchema (OneFile *vf, OneSchema *vs, bool isRequired)
{
  bool isMatch = true ;
  int  i, j ;

  if (!vf || !vs) return false ;

  if (vs->nxt) // the textSchema contained at least one 'P' line to define a file type
    { while (vs && (!vs->primary || strcmp (vs->primary, vf->fileType))) vs = vs->nxt ;
      if (!vs)
	{ snprintf (errorString, 1024, "OneSchema mismatch: file type %s not found in schema\n",
		   vf->fileType) ;
	  return false ;
	}
    }

  // at this point vs->primary matches vf->fileType

  for (i = 'A' ; i <= 'z' ; ++i)
    { OneInfo *vis = vs->info[i] ;
      OneInfo *vif = vf->info[i] ;
      if (isRequired && vis && !vif)
	{ snprintf (errorString, 1024, "OneSchema mismatch: record type %c missing in file schema\n", i) ;
	  isMatch = false ;
	}
      else if (vis && vif)
	{ if (vif->isObject != vis->isObject)
	    { snprintf (errorString, 1024, "OneSchema mismatch: object type %c file %d != schema %d\n",
		       i, vif->isObject, vis->isObject) ;
	      isMatch = false ;
	    }
	  if (vif->nField != vis->nField)
	    { snprintf (errorString, 1024, "OneSchema mismatch: number of fields for type %c file %d != schema %d\n",
		       i, vif->nField, vis->nField) ;
	      isMatch = false ;
	    }
	  else
	    for (j = 0 ; j < vif->nField ; ++j)
	      if (vif->fieldType[j] != vis->fieldType[j])
		{ snprintf (errorString, 1024, "OneSchema mismatch: field %d for type %c file %s != schema %s\n",
			   j,i,oneTypeString[vif->fieldType[j]],oneTypeString[vis->fieldType[j]]);
		  isMatch = false ;
		}
	}
    }

  return isMatch ;
}

bool oneFileCheckSchemaText (OneFile *vf, const char *textSchema)
{
  char      *fixedText = schemaFixNewlines (textSchema) ;
  OneSchema *vs = oneSchemaCreateFromText (fixedText) ;
  bool       isCheck = oneFileCheckSchema (vf, vs, true) ;

  free (fixedText) ;
  oneSchemaDestroy (vs) ;
  return isCheck ;
}

/***********************************************************************************
 *
 *    SETTING UP PROVENANCE, REFERENCES, & DEFERRALS
 *
 **********************************************************************************/

bool addProvenance(OneFile *vf, OneProvenance *from, int n)
{ I64 i ;
  OneInfo   *l = vf->info['!'];
  I64         o = l->accum.count;
  OneProvenance *p;

  if (n == 0)
    return (false);
  assert (!vf->isHeaderOut) ;

  l->accum.count += n;

  p = new(o+n, OneProvenance);
  if (o > 0)
    memcpy (p, vf->provenance, o*sizeof(OneProvenance));
  memcpy (p+o, from, n*sizeof(OneProvenance));
  free (vf->provenance);
  vf->provenance = p;

  // finally create self-owned copy of all fields

  p = p+o ;
  for (i = 0 ; i < n ; ++i, ++p)
    { p->program = strdup(p->program) ;
      p->version = strdup(p->version) ;
      p->command = strdup(p->command) ;
      p->date = strdup(p->date) ;
    }

  return (true);
}

bool oneInheritProvenance(OneFile *vf, OneFile *source)
{ return (addProvenance(vf, source->provenance, source->info['!']->accum.count)); }

bool oneAddProvenance(OneFile *vf, const char *prog, const char *version, char *format, ...)
{ va_list args ;
  OneProvenance p;
  time_t t = time(NULL);

  p.program = (char*) prog; // cast to keep compiler happy - this is safe!
  p.version = (char*) version; // cast to keep compiler happy - this is safe!
  va_start (args, format) ;
  if (vasprintf (&p.command, format, args) == -1) die ("vasprintf failure") ;
  va_end (args) ;
  p.date = new (20, char);
  strftime(p.date, 20, "%F_%T", localtime(&t));
  addProvenance (vf, &p, 1);
  free (p.command) ;
  free (p.date) ;
  return true ; // always added something
}

static bool addReference(OneFile *vf, OneReference *from, int n, bool isDeferred)
{ I64        o;
  OneInfo  *l;
  OneReference *r, **t;
  I64 i ;

  if (n == 0)
    return false;
  assert (!vf->isHeaderOut) ;

  if (isDeferred)
    { l = vf->info['>'];
      t = &(vf->deferred);
    }
  else
    { l = vf->info['<'];
      t = &(vf->reference);
    }
  o = l->accum.count;
  l->accum.count += n;

  r = new (o+n, OneReference);
  if (o > 0)
    memcpy (r, *t, o*sizeof(OneReference));
  memcpy (r+o, from, n*sizeof(OneReference));
  free (*t);
  *t = r;

  r += o ; // make self-owned copy of filename strings
  for (i = 0 ; i < n ; ++i, ++r)
    r->filename = strdup (r->filename) ;

  return true;
}

bool oneInheritReference(OneFile *vf, OneFile *source)
{ return (addReference(vf, source->reference, source->info['<']->accum.count, false)); }

bool oneAddReference(OneFile *vf, const char *filename, I64 count)
{ OneReference ref;
  ref.filename = (char*) filename; // cast to keep compiler happy - this is safe!
  ref.count    = count;
  return (addReference(vf, &ref, 1, false));
}

bool oneInheritDeferred (OneFile *vf, OneFile *source)
{ return (addReference (vf, source->deferred, source->info['>']->accum.count, true)); }

bool oneAddDeferred (OneFile *vf, const char *filename)
{ OneReference ref;
  ref.filename = (char *) filename; // cast to keep compiler happy - this is safe!
  return (addReference (vf, &ref, 1, true));
}

bool oneStats (OneFile *of, char lineType, I64 *count, I64 *max, I64 *total)
{
  OneInfo    *info = of->info[(int)lineType] ;
  if (!info)  return false ;
  
  OneCounts   counts = of->isWrite ? info->accum : info->given ;
  if (count) *count = counts.count ;
  if (max)   *max   = counts.max ;
  if (total) *total = counts.total ;
  return true ;
}

bool  oneStatsContains (OneFile *of, char objectType, char lineType, I64 *maxCount, I64 *maxTotal)
{
  OneInfo   *info = of->info[(int)objectType] ;
  if (!info || !info->stats)  return false ; 
  OneStat   *s ;
  for (s = info->stats ; s->type && s->type != lineType ; ++s) ;
  if (!s->type) return false ;
  if (maxCount) *maxCount = s->maxCount ;
  if (maxTotal) *maxTotal = s->maxTotal ;
  return true ;
}

/***********************************************************************************
 *
 *   ONE_WRITE_HEADER / FOOTER
 *
 **********************************************************************************/

static bool writeCounts (OneFile *vf, int i) // always write counts in ascii
{
  OneInfo *li = vf->info[i] ;

  if (li != NULL && li->given.count > 0)
    { fprintf (vf->f, "# %c %lld\n", i, li->given.count);
      if (li->given.max > 0)
	fprintf (vf->f, "@ %c %lld\n", i, li->given.max);
      if (li->given.total > 0)
	fprintf (vf->f, "+ %c %lld\n", i, li->given.total);
      if (li->isObject)
	{ OneStat *s ;
	  for (s = li->stats ; s->type ; ++s)
	    { if (s->maxCount)
		fprintf (vf->f, "%% %c # %c %lld\n", i, s->type, s->maxCount);
	      if (s->maxTotal)
		fprintf (vf->f, "%% %c + %c %lld\n", i, s->type, s->maxTotal);
	    }
	}
      return true ;
    }
  else
    return false ; 
}

static void writeHeader (OneFile *vf)
{ int         i,n;

  assert (vf->isWrite) ;
  assert (vf->share >= 0) ;

  vf->isLastLineBinary = false; // header is in ASCII

  fprintf (vf->f, "1 %lu %s %d %d", strlen(vf->fileType), vf->fileType, MAJOR, MINOR);
  if (vf->subType)
    fprintf (vf->f, "\n2 %lu %s", strlen(vf->subType), vf->subType);

  // any header text on '.' lines
  if (vf->headerText)
    { OneHeaderText *t = vf->headerText ;
      while (t)
	{ fprintf (vf->f, "\n. %s", t->text) ;
	  t = t->nxt ;
	}
    }

  // provenance
  if (vf->info['!']->accum.count)
    { OneProvenance *p = vf->provenance; 
      n = vf->info['!']->accum.count;
      for (i = 0; i < n; i++, p++)
	fprintf (vf->f, "\n! 4 %lu %s %lu %s %lu %s %lu %s",
		 strlen(p->program), p->program, strlen(p->version), p->version,
		 strlen(p->command), p->command, strlen(p->date), p->date);
    }

  fprintf (vf->f, "\n.") ; // always have a spacer after this

  // reference and deferred
  if (vf->info['<']->accum.count || vf->info['>']->accum.count)
    { OneReference *r = vf->reference;
      n = vf->info['<']->accum.count;
      for (i = 0; i < n; i++, r++)
	fprintf (vf->f, "\n< %lu %s %lld", strlen(r->filename), r->filename, r->count);
      
      r = vf->deferred;
      n = vf->info['>']->accum.count;
      for (i = 0; i < n; i++, r++)
	fprintf (vf->f, "\n> %lu %s", strlen(r->filename), r->filename);
      fprintf (vf->f, "\n.") ;
    }

  // write the schema into the header - no need for file type, version etc. since already given
  for (i = 0 ; i < vf->nDefn ; ++i)
    writeInfoSpec (vf->f, vf, vf->defnOrder[i], vf->defnComment[i]) ;

  if (vf->isBinary)         // defer writing rest of header
    fprintf (vf->f, "\n$ %d", vf->isBig);
  else                      // write counts based on those supplied in info[i].given
    { fprintf (vf->f, "\n.\n") ;
      for (i = 0 ; i < vf->nDefn ; ++i)
	if (!(vf->defnOrder[i] & 0x80))
	  writeCounts (vf, vf->defnOrder[i]) ;
      fprintf (vf->f, ".") ; // need to set up an incomplete line
    }
  fflush (vf->f);

  vf->isHeaderOut = true;
}

/***********************************************************************************
 *
 *   ONE_WRITE_LINE
 *
 **********************************************************************************/

static int writeStringList (OneFile *vf, char t, int len, char *buf)
{ OneInfo *li;
  int       j, nByteWritten = 0;
  I64       sLen, totLen;

  totLen = 0;
  for (j = 0; j < len; j++)
    { sLen = strlen (buf);
      totLen += sLen;
      nByteWritten += fprintf (vf->f, " %lld %s", sLen, buf);
      buf += sLen + 1;
    }

  li = vf->info[(int) t];
  li->accum.total += totLen;
  if (li->accum.max < totLen)
    li->accum.max = totLen;

  return nByteWritten ;
}

// code to track counts for objects

static inline void startObject (OneFile *vf, OneInfo *li)
{
  OneStat *s ;
  for (s = li->stats ; s->type ; ++s)
    { s->count = vf->info[(int)s->type]->accum.count ;
      if (s->isList) s->total = vf->info[(int)s->type]->accum.total ;
    }
  if (li->accum.count == 1) // must record count/total before first object in file
    for (s = li->stats ; s->type ; ++s)
      { s->count0 = vf->info[(int)s->type]->accum.count ;
	if (s->isList) s->total0 = vf->info[(int)s->type]->accum.total ;
      }
  vf->openObjects[++vf->objectFrame] = li ;
}

static inline void endObject (OneFile *vf, OneInfo *li)
{
  OneStat *s ;
  for (s = li->stats ; s->type ; ++s)
    { if (vf->info[(int)s->type]->accum.count - s->count > s->maxCount)
	s->maxCount = vf->info[(int)s->type]->accum.count - s->count ;
      if (s->isList && vf->info[(int)s->type]->accum.total - s->total > s->maxTotal)
	s->maxTotal = vf->info[(int)s->type]->accum.total - s->total ;
    }
  --vf->objectFrame ;
}

static inline void closeObjects (OneFile *vf, char t) // set count0 for any objects terminated by t
{
  int i, *ik = vf->defnOrder ;
  for (i = 0 ; i < vf->nDefn ; ++i, ++ik)
    if (!(*ik & 0x80))
      { OneInfo *li = vf->info[*ik] ;
	if (li->isObject && !li->isClosed && !li->contains[(int)t])
	  { OneStat *s = li->stats ;
	    while (s->type)
	      { s->count0 = vf->info[(int)s->type]->accum.count ;
		if (s->isList) s->total0 = vf->info[(int)s->type]->accum.total ;
		++s ;
	      }
	    li->isClosed = true ;
	  }
      }
  vf->info[(int) t]->isFirst = false ;
}
 
// process is to fill fields by assigning to macros, then call - list contents are in buf
// NB in ASCII mode adds '\n' before writing line not after, so oneWriteComment() can add to line
// first call will write initial header

void oneWriteLine (OneFile *vf, char t, I64 listLen, void *listBuf)
{ I64      i, j;
  OneInfo *li;

  // fprintf (stderr, "write type %c char %c listLen %d\n", t, oneChar(vf,0), (int) listLen) ;
  
  assert (vf->isWrite) ;
  assert (!vf->isFinal || !isalpha(t)) ;
  
  li = vf->info[(int) t];
  if (!li) die ("oneWriteLine() attempting to write unkown linetype %c", t) ;

  if (li->isFirst) closeObjects (vf, t) ;
  while (vf->objectFrame && !(vf->openObjects[vf->objectFrame]->contains[(int)t]))
    endObject (vf, vf->openObjects[vf->objectFrame]) ;
  li->accum.count += 1;
  if (li->isObject) startObject (vf, li) ;

  if (li->listEltSize > 0)  // need to write the list
    { assert (listLen >= 0) ;
      vf->field[li->listField].len = listLen ;
      if (listBuf == NULL) listBuf = li->buffer;
    }

  // BINARY - block write and optionally compress

  if (vf->isBinary)
    { U8  x;

      if (!vf->isHeaderOut && vf->share >= 0) // no header on slaves
	{ writeHeader (vf) ;
	  if (!vf->isLastLineBinary) // copied from below because need to set vf->byte before writing the index for object 0
	    { fputc ('\n', vf->f) ;
	      vf->byte = ftello (vf->f) ;
	    }
	  for (i = 'A' ; i <= 'z' ; i++) // write index[0] to be here at start of data
	    if (vf->info[i] && vf->info[i]->index)
	      vf->info[i]->index[0] = vf->byte ;  // OK to only do this for master - slaves start at 0
	}
      else if (!vf->isLastLineBinary)
	{ fputc ('\n', vf->f) ;
	  vf->byte = ftello (vf->f) ;
	}

      if (li->isObject) // update index
	{ if (li->accum.count >= li->indexSize)
	    { I64 oldSize = li->indexSize ;
	      li->indexSize = (oldSize << 2) + 0x10000 ;
	      resize (li->index, oldSize, li->indexSize, I64) ;
	    }
	  li->index[li->accum.count] = vf->byte ;
          // assert (ftello (vf->f) == vf->byte) ; // beware - very costly
	}

      // write the line character
      
      x = li->binaryTypePack;   //  Binary line code + compression flags
      if (li->isUseListCodec)
        x |= 0x01;
      fputc (x, vf->f);
      ++vf->byte ;

      // write the fields

      if (li->nField > 0)
	vf->byte += writeCompressedFields (vf->f, vf->field, li) ;

      // write the list if there is one

      if (li->listEltSize && listLen > 0)
        { I64 nBits, listSize;
	  int listBytes ;

	  li->accum.total += listLen;
          if (listLen > li->accum.max)
            li->accum.max = listLen;
	  
	  if (li->fieldType[li->listField] == oneINT_LIST)
	    { vf->byte += ltfWrite (*(I64*)listBuf, vf->f) ;
	      if (listLen == 1) goto doneLine ; // finish writing this line here
	      listBuf = compactIntList (vf, li, listLen, listBuf, &listBytes) ;
	      --listLen ;
	      fputc ((char)listBytes, vf->f) ;
	      vf->byte++ ;
	    }
	  else
	    listBytes = li->listEltSize ;
	  listSize  = listLen * listBytes;
	  
	  if (li->fieldType[li->listField] == oneSTRING_LIST) // handle as ASCII
	    vf->byte += writeStringList (vf, t, listLen, listBuf);
	  else if (x & 0x1)
	    { if (listSize >= vf->codecBufSize)
		{ free (vf->codecBuf);
		  vf->codecBufSize = listSize+1;
		  vf->codecBuf     = new (vf->codecBufSize, void);
		}
	      nBits = vcEncode (li->listCodec, listSize, listBuf, vf->codecBuf);
	      vf->byte += ltfWrite (nBits, vf->f) ;
	      if (fwrite (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1)
		die ("ONE write error: failed to write compressed list nBits %lld", nBits);
	      vf->byte += ((nBits+7) >> 3) ;
	    }
	  else
	    { if (fwrite (listBuf, listSize, 1, vf->f) != 1)
		die ("ONE write error: failed to write list field %d listLen %lld listSize %lld listBuf %lx",
		     li->listField, listLen, listSize, listBuf);
	      vf->byte += listSize;
	      if (li->listCodec != NULL)
		{ vcAddToTable (li->listCodec, listSize, listBuf);
		  li->listTack += listSize;
		  
		  if (li->listTack > vf->codecTrainingSize)
		    { if (vf->share == 0)
			{ vcCreateCodec (li->listCodec, 1);
			  li->isUseListCodec = true;
			}
		      else
			{ OneFile  *ms;
			  OneInfo *lx;
			  
			  if (vf->share < 0)
			    { ms = vf + vf->share;
			      lx = ms->info[(int) t]; 
			    }
			  else
			    { ms = vf;
			      lx = li;
			    }
			  
			  pthread_mutex_lock(&ms->listLock);
			  
			  if ( ! li->isUseListCodec)
			    
			    { if (vf->share < 0)
				{ lx->listTack += li->listTack;
				  li->listTack = 0;
				}
			      if (lx->listTack > ms->codecTrainingSize)
				{ for (i = 1; i < ms->share; i++)
				    vcAddHistogram (lx->listCodec,
						    ms[i].info[(int) t]->listCodec);
				  vcCreateCodec (lx->listCodec, 1);
				  for (i = 1; i < ms->share; i++)
				    { OneCodec *m = ms[i].info[(int) t]->listCodec;
				      ms[i].info[(int) t]->listCodec = lx->listCodec;
				      vcDestroy (m);
				    }
				  lx->isUseListCodec = true;
				  for (i = 1; i < ms->share; i++)
				    ms[i].info[(int) t]->isUseListCodec = true;
				}
			    }
			  
			  pthread_mutex_unlock(&ms->listLock);
			}
		    }
		}
	    }
	}

    doneLine:

      vf->isLastLineBinary = true;
    }

  // ASCII - write field by field

  else
    { if (!vf->isHeaderOut && !vf->isNoAsciiHeader) writeHeader (vf) ;

      if (!vf->isLastLineBinary)      // terminate previous ascii line
	fputc ('\n', vf->f);

      ++vf->line ; // only really needed when closing the file to see if we need to terminate it
      
      fputc (t, vf->f);

      for (i = 0; i < li->nField; i++)
        switch (li->fieldType[i])
	  {
	  case oneINT:
            fprintf (vf->f, " %lld", vf->field[i].i);
            break;
          case oneREAL:
            fprintf (vf->f, " %f", vf->field[i].r);
            break;
          case oneCHAR:
            fprintf (vf->f, " %c", vf->field[i].c);
            break;
          case oneSTRING:
	  case oneDNA:
          case oneINT_LIST:
          case oneREAL_LIST:
          case oneSTRING_LIST:
            li->accum.total += listLen;
            if (listLen > li->accum.max)
              li->accum.max = listLen;

	    fprintf (vf->f, " %lld", listLen);
            if (li->fieldType[i] == oneSTRING || li->fieldType[i] == oneDNA)
              { if (listLen > INT_MAX)
                  die ("ONE write error: string length %lld > current max %d", listLen, INT_MAX);
                fprintf (vf->f, " %.*s", (int) listLen, (char *) listBuf);
              }
            else if (li->fieldType[i] == oneINT_LIST)
              { I64 *b = (I64 *) listBuf;
                for (j = 0; j < listLen ; ++j)
                  fprintf (vf->f, " %lld", b[j]);
              }
            else if (li->fieldType[i] == oneREAL_LIST)
              { double *b = (double *) listBuf;
                for (j = 0; j < listLen ; ++j)
                  fprintf (vf->f, " %f", b[j]);
              }
            else // vSTRING_LIST
              writeStringList (vf, t, listLen, listBuf);
            break;
        }
      vf->isLastLineBinary = false;
    }
}

int Uncompress_DNA(char *s, int len, char *t) ; // forward declaration for temp solution below

void oneWriteLineDNA2bit (OneFile *vf, char lineType, I64 len, U8 *dnaBuf) // NB len in bp
{ // temporary solution
  char *s = new(len, char) ;
  Uncompress_DNA ((char*)dnaBuf, len, s) ;
  oneWriteLine (vf, lineType, len, s) ;
  free (s) ;
}

void oneWriteComment (OneFile *vf, char *format, ...)
{
  char *comment ;
  
  va_list args ;
  va_start (args, format) ; 
  if (vasprintf (&comment, format, args) == -1) die ("vasprintf failure") ;
  va_end (args) ;

  if (vf->isCheckString) // then check no newlines in format
    { char *s = format ;
      while (*s) if (*s++ == '\n') die ("newline in comment string: %s", comment) ;
    }

  if (vf->isLastLineBinary) // write a comment line
    oneWriteLine (vf, '/', strlen(comment), comment) ;
  else // write on same line after space
    { fputc (' ', vf->f) ;
      fprintf (vf->f, "%s", comment) ;
    }
  free (comment) ;
}

/***********************************************************************************
 *
 *    MERGING, FOOTER HANDLING, AND CLOSE
 *
 **********************************************************************************/

static void oneWriteFooter (OneFile *vf)
{ int      i,k,n;
  off_t    footOff;
  OneInfo *li;
  char    *codecBuf ;

  footOff = ftello (vf->f);
  if (footOff < 0)
    die ("ONE write error: failed footer ftell");

  //  first the per-linetype information
  codecBuf = new (vcMaxSerialSize()+1, char) ; // +1 for added up unused 0-terminator
  bool isWrittenIndexCodec = false ;
  for (k = 0; k < vf->nDefn ; ++k)
    { i  = vf->defnOrder[k] ;
      if (i & 0x80) continue ; // skip the 'G' lines
      li = vf->info[i];
      if (li->accum.count > 0)
        { li->given = li->accum ;
	  writeCounts (vf, i) ;
	  if (li->index)
	    { oneChar(vf,0) = (char) i ;
	      oneWriteLine (vf, '&', li->accum.count+1, li->index) ;
	    }
	  if (vf->info['&']->isUseListCodec && !isWrittenIndexCodec)
	    { oneChar(vf,0) = '&' ;
              n = vcSerialize (vf->info['&']->listCodec, codecBuf);
              oneWriteLine (vf, ';', n, codecBuf);
	      isWrittenIndexCodec = true ;
	    }
          if (li->isUseListCodec && li->listCodec != DNAcodec)
            { oneChar(vf,0) = i;
              n = vcSerialize (li->listCodec, codecBuf);
              oneWriteLine (vf, ';', n, codecBuf);
            }
        }
    }

  li = vf->info['/'] ;		// may need to write list codec for comments
  if (li->isUseListCodec)
    { oneChar(vf,0) = '/' ;
      n = vcSerialize (li->listCodec, codecBuf);
      oneWriteLine (vf, ';', n, codecBuf);
    }

  // NB we don't consider here the possibility of writing a codec for codecs.
  // This should be OK. Each codec is ~300 bytes, and we can have at most 56 of them.
  // Not enough to trigger codec building.
  
  free (codecBuf) ;

  fprintf (vf->f, "^\n"); // end of footer marker

  if (fwrite (&footOff, sizeof(off_t), 1, vf->f) != 1)
    die ("ONE write error: failed writing footer offset");
}

  // After all input has been read, or all data has been written, this routine will finish
  //   accumulating counts/statistics for the file and merge thread stats into those for
  //   the master file (if a parallel OneFile).

void oneFinalizeCounts(OneFile *vf)
{ int      ii, j, k ;
  OneInfo *li, *lk;

  if (vf->share < 0)
    die ("ONE write error: cannot call oneFileClose on a slave OneFile");

  vf->isFinal = true; // needed to prevent infinite recursion

  if (vf->share == 0)
    { while (vf->objectFrame)
	endObject (vf, vf->openObjects[vf->objectFrame]) ; // terminate open objects
      return; 
    }

  int nthreads = vf->share; // if we get here then nthreads > 1
  
  // first we need to complete any objects left open at the end of files and update max count/total
  OneFile *vk, *vk1 ;
  OneStat *s, *s1 ;
  for (k = 1 ; k < nthreads ; ++k)
    { vk = &vf[k] ; vk1 = &vf[k-1] ;
      while (vk1->objectFrame) // try to complete objects open at end of preceding file k-1
	{ OneInfo *li1 = vk1->openObjects[vk1->objectFrame] ; // the object to close in file k-1
	  // We must find the corresponding object in li. The following is a bit ugly.
	  { int i ;
	    for (i = 'A' ; i <= 'z' ; ++i) if (vk1->info[i] == li1) break ; // found it
	    if (i <= 'z') li = vk->info[i] ; else die ("failed to find li") ;
	  }
	  if (li->isClosed)                                     // yes we can close it
	    for (s = li->stats, s1 = li1->stats ; s->type ; ++s, ++s1)
	      { if (vk1->info[(int)s->type]->accum.count - s1->count + s->count0 > s->maxCount)
		  s->maxCount = vk1->info[(int)s->type]->accum.count - s1->count + s->count0 ;
		if (s->isList && vk1->info[(int)s->type]->accum.total - s1->total + s->total0 > s->maxTotal)
		  s->maxTotal = vk1->info[(int)s->type]->accum.total - s1->total + s->total0 ;
	      }
	  else // add to this file k's open list, adjusting count/total to include previous counts
	    { for (s = li->stats, s1 = li1->stats ; s->type ; ++s, ++s1)
		{ s->count -= vk1->info[(int)s->type]->accum.count - s1->count ;
		  if (s->isList) s->total -= vk1->info[(int)s->type]->accum.total - s1->total ;
		}
	      vk->openObjects[++vk->objectFrame] = li ;
	    }
	  --vk1->objectFrame ;
	}
      if (k == nthreads-1) // close any remaining open objects at the end of the final file
	while (vk->objectFrame)
	  endObject (vk, vk->openObjects[vk->objectFrame]) ;
      
      // now see if we need to update max count/total in vf for anything
      for (ii = 0 ; ii < vf->nDefn ; ++ii)
	{ int i = vf->defnOrder[ii] ;
	  if (i & 0x80) continue ; // skip 'G' lines
	  li = vf->info[i] ;
	  if (li->isObject)
	    for (s = li->stats, s1 = vk->info[i]->stats ; s->type ; ++s, ++s1)
	      { if (s1->maxCount > s->maxCount) s->maxCount = s1->maxCount ;
		if (s->isList && s1->maxTotal > s->maxTotal) s->maxTotal = s1->maxTotal ;
	      }
	}
    }
  
  // next update the li->accum - must have fixed up the max count/total first since they use accum
  for (ii = 0 ; ii < vf->nDefn ; ++ii)
    { int i = vf->defnOrder[ii] ;
      if (i & 0x80) continue ; // skip 'G' lines
      li = vf->info[i] ;
      I64 n0 = li->accum.count ;
      for (k = 1 ; k < nthreads ; ++k)
	{ lk = vf[k].info[i] ;
	  li->accum.count += lk->accum.count ;
	  li->accum.total += lk->accum.total ;
	  if (lk->accum.max > li->accum.max) li->accum.max = lk->accum.max ;
	}

      // finally stitch together the index - need to have fixed li->accum.count first
      if (vf->isBinary && li->isObject)
	{ I64 oldIndexSize = li->indexSize ;
	  li->indexSize = li->accum.count+1 ;
	  resize (li->index, oldIndexSize, li->indexSize, I64) ;
	  I64 off = ftello(vf->f) ;
	  I64 n = n0 ;
	  for (k = 1 ; k < nthreads ; ++k)
	    { I64  nk = vf[k].info[i]->accum.count ;
	      I64 *kIndex = vf[k].info[i]->index ;
	      for (j = 1 ; j <= nk ; ++j)
		li->index[++n] = kIndex[j] + off;
	      off += ftello(vf[k].f);
	    }
	}
    }
}

//

static void oneFinalize (OneFile *vf)
{
  if (!vf->isFinal)
    oneFinalizeCounts (vf);

  if (!vf->isHeaderOut && (vf->isBinary || !vf->isNoAsciiHeader)) writeHeader (vf) ;
      
  if (vf->share > 0)
    { int  i, nread ;
      char *buf;
      buf = new (10000000, char);
      for (i = 1; i < vf->share; i++)
	{ if (fseek (vf[i].f, 0L, SEEK_SET) != 0)
	    die ("ONEfile error: failed to rewind parallel file %d", i) ;
	  while (!feof(vf[i].f) && (nread = fread (buf,1,10000000,vf[i].f)) > 0)
	    if ((int) fwrite(buf,1,nread,vf->f) != nread)
	      die ("ONE write error: while cat'ing thread bits (oneFileClose)");
	}
      free(buf);
    }

  if (vf->isBinary || vf->line)
    fputc ('\n', vf->f) ; // terminate last line - end of data marker if binary
  if (vf->isBinary) // write the footer
    { if (!vf->isLastLineBinary)
	fputc ('\n', vf->f);  // need an extra '\n' to ensure end of data marker
      oneWriteFooter (vf);
    }
}

void oneFileClose (OneFile *vf)
{
  assert (vf->share >= 0) ;

  if (vf->isWrite)
    oneFinalize (vf) ;
  
  oneFileDestroy (vf);
}

/***********************************************************************************
 *
 *  Length limited Huffman Compressor/decompressor with special 2-bit compressor for DNA
 *  Author:          Gene Myers
 *  Creation date:   June 27, 2019
 *
 *  inline both compression.h and compression.c here
 *
 **********************************************************************************/

#undef  DEBUG
#undef  TEST

  //  To create a compressor, get an initially empty object with vcCreate, then
  //    add a significant corpus of the byte data to be compressed with vcAddToTable,
  //    and finally create a Huffman codec based on this corpus by calling
  //    vcCreateCodec.  The parameter "partial" should be set if not all the data
  //    to be compressed has been scanned.  At this point you have a compressor ready
  //    to operate.  You can destroy/free it with vcDestroy.

OneCodec *vcCreate();
void      vcAddToTable(OneCodec *vc, int len, char *bytes);
void      vcCreateCodec(OneCodec *vc, int partial);
void      vcDestroy(OneCodec *vc);

  //  In the instance of accumulating data over multiple threads, vcAddHistogram, will
  //    add the counts in the table for vh, to the table for vc.

void      vcAddHistogram(OneCodec *vc, OneCodec *vh);

  //  A diagnostic routine: shows you the compression scheme and if the distribution
  //    of the scanned corpus is available, it shows you that too.  Output to file 'to'.

void      vcPrint(OneCodec *vc, FILE *to);

  //  You can encode and decode where ibytes/ilen are the input and the output
  //    is placed at obytes and the length of the compressed/decompressed result
  //    is returned as the value of the function.  For vcEncode, ilen is the size
  //    of the uncompressed input in bytes, and the return value is the size of
  //    the compressed output in **bits**.  The converse is true for vcDecode, i.e
  //    ilen is the number of bits in the compressed input, and the return value
  //    is the number of bytes in the uncompressed output.  The routines are endian safe.

int       vcEncode(OneCodec *vc, int ilen, char *ibytes, char *obytes);
int       vcDecode(OneCodec *vc, int ilen, char *ibytes, char *obytes);

  //  Rather than directly reading or writing an encoding of a compressor, the routines
  //    below serialize or deserialize the compressor into/outof a user-supplied buffer.
  //    vcMaxSerialSize gives the maximum size of a serialized compressor so the user
  //    can arrange a buffer of the appropriate size.  vcSerialize serializes vc into
  //    buffer 'out' and returns the # of bytes in the encoding.  vcDeserialize will reverse
  //    the process given a serialization.  The routines are endian-safe.

int       vcMaxSerialSize();
int       vcSerialize(OneCodec *vc, void *out);
OneCodec *vcDeserialize(void *in);

typedef unsigned long long  uint64;
typedef unsigned int        uint32;
typedef unsigned short      uint16;
typedef unsigned char       uint8;

#define HUFF_CUTOFF  12     //  This cannot be larger than 16 !

  //  Endian flipping macros

#define FLIP64(p)	\
{ uint8 x = p[0];	\
  p[0] = p[7];		\
  p[7] = x;		\
  x     = p[1];		\
  p[1] = p[6];		\
  p[6] = x;		\
  x     = p[2];		\
  p[2] = p[5];		\
  p[5] = x;		\
  x     = p[3];		\
  p[3] = p[4];		\
  p[4] = x;		\
}

#define FLIP32(p)	\
{ uint8 x = p[0];	\
  p[0] = p[3];		\
  p[3] = x;		\
  x     = p[1];		\
  p[1] = p[2];		\
  p[2] = x;		\
}

#define FLIP16(p)	\
{ uint8 x = p[0];	\
  p[0] = p[1];		\
  p[1] = x;		\
}

/*******************************************************************************************
 *
 *  Routines for computing a length-limited Huffman Encoding Scheme
 *
 ********************************************************************************************/

#define EMPTY        0      //  Compressor just created, histogram zero'd
#define FILLED       1      //  Compressor histogram being filled, no codec
#define CODED_WITH   2      //  Compressor has a codec (can no longer accumulate histogram)
#define CODED_READ   3      //  Compressor has codec but no histogram as was created by read

typedef struct
  { int    state;            //  1 of the 4 states immediately above
    int    isbig;            //  endian of the current machine
    uint16 codebits[256];    //  Code esc_code is the special code for
    uint8  codelens[256];    //    non-Huffman exceptions
    char   lookup[0x10000];  //  Lookup table (just for decoding)
    int    esc_code;         //  The special escape code (-1 if not partial)
    int    esc_len;          //  The length in bits of the special code (if present)
    uint64 hist[256];        //  Byte distribution for codec
  } _OneCodec;

  //  The special "predefined" DNA compressor

static _OneCodec _DNAcodec = { .state = CODED_READ };
OneCodec  *DNAcodec = (OneCodec *) &_DNAcodec;

  //  Create an EMPTY compressor object with zero'd histogram and determine machine endian

OneCodec *vcCreate()
{ _OneCodec *v;
  int i;

  v = (_OneCodec *) malloc(sizeof(_OneCodec));
  if (v == NULL) die ("vcCreate: Could not allocate compressor") ;

  v->state = EMPTY;
  for (i = 0; i < 256; i++)
    v->hist[i] = 0;

  { uint32 t;
    uint8 *b;

    t = 1;
    b = (uint8 *) (&t);
    v->isbig = (b[0] == 0);
  }

  return ((OneCodec *) v);
}

  //  Free a compressor object

void vcDestroy(OneCodec *vc)
{ _OneCodec *v = (_OneCodec *) vc;
  if (vc != DNAcodec)
    free(v);
}

  //  Add the frequencies of bytes in bytes[0..len) to vc's histogram
  //    State becomes FILLED

void vcAddToTable(OneCodec *vc, int len, char *bytes)
{ _OneCodec *v = (_OneCodec *) vc;
  uint8 *data = (uint8 *) bytes;
  int i;

  for (i = 0; i < len; i++)
    v->hist[(int) data[i]] += 1;
  if (v->state < FILLED)
    v->state = FILLED;
}

  //  Add the frequencies of bytes in bytes[0..len) to vc's histogram
  //    State becomes FILLED

void vcAddHistogram(OneCodec *vc, OneCodec *vh)
{ _OneCodec *v = (_OneCodec *) vc;
  _OneCodec *h = (_OneCodec *) vh;
  int i;

  if (v->state >= CODED_WITH) die("vcAddHistogram: Compressor already has a codec");
  if (h->state == CODED_READ) die("vcAddHistogram: Source compressor doesn't have a histogram");

  for (i = 0; i < 256; i++)
    v->hist[i] += h->hist[i];
  v->state = FILLED;
}

  //  Check vc has a non-empty distribution histogram and if so then build
  //    length-limited Huffman tables for the bytes that occur in the histogram,
  //    plus a special escape code if partial is set and there is at least one byte
  //    with a zero count in the histogram.  The algorithm is by Larmore & Hirschberg,
  //    JACM 73, 3 (1990).

uint64 *HIST;

int HSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);
  return (HIST[x] - HIST[y]);
}

void vcCreateCodec(OneCodec *vc, int partial)
{ _OneCodec *v = (_OneCodec *) vc;

  uint64  *hist;
  char    *look;
  uint8   *lens;
  uint16  *bitv;

  int      code[256];
  int      leng[256];
  uint16   bits[256];
  int      ncode, dcode, ecode;

  int      i;

  if (v->state >= CODED_WITH) die("vcCreateCoder: Compressor already has a codec");
  if (v->state == EMPTY) die("vcCreateCoder: Compressor has no byte distribution data");

  hist  = v->hist;
  look  = v->lookup;
  lens  = v->codelens;
  bitv  = v->codebits;

  ecode = -partial;
  ncode = 0;
  for (i = 0; i < 256; i++)
    if (hist[i] > 0)
      code[ncode++] = i;
    else if (ecode < 0)
      { ecode = i;
        code[ncode++] = i;
      }
  dcode = 2*ncode;

  if (ecode < 0)
    partial = 0;

  HIST = hist;
  qsort(code,ncode,sizeof(int),HSORT);

#ifdef DEBUG
  fprintf(stderr,"\nSorted Codes %d:\n",ncode);
  for (i = 0; i < ncode; i++)
    fprintf(stderr," %3d: %3d %10llu\n",i,code[i],hist[code[i]]);
#endif

  { uint8   matrix[HUFF_CUTOFF][dcode];
    uint64  count1[dcode], count2[dcode], countb[ncode];
    uint64 *lcnt, *ccnt, *swp;
    int     llen, span;
    int     j, k, n, L;

    for (n = 0; n < ncode; n++)
      { count1[n] = countb[n] = hist[code[n]];
        leng[n] = 0;
      }

#ifdef DEBUG
    fprintf(stderr,"\nCoin Filter:\n");
    fprintf(stderr,"  Row %2d:",HUFF_CUTOFF);
    for (n = 0; n < ncode; n++)
      fprintf(stderr," %lld*",countb[n]);
    fprintf(stderr,"\n");
#endif

    lcnt = count1;
    ccnt = count2;
    llen = ncode-1;
    for (L = HUFF_CUTOFF-1; L > 0; L--)
      { j = 0;
        k = 0;
        for (n = 0; j < ncode || k < llen; n++)
          { if (k >= llen || (j < ncode && countb[j] <= lcnt[k] + lcnt[k+1]))
              { ccnt[n] = countb[j];
                matrix[L][n] = 1;
                j += 1;
              }
            else
              { ccnt[n] = lcnt[k] + lcnt[k+1];
                matrix[L][n] = 0;
                k += 2;
              }
          }
        llen = n-1;
        swp  = lcnt;
        lcnt = ccnt;
        ccnt = swp;

#ifdef DEBUG
        fprintf(stderr,"  Row %2d:",L);
        for (n = 0; n <= llen; n++)
          fprintf(stderr," %lld%c",lcnt[n],matrix[L][n]?'*':'+');
        fprintf(stderr,"\n");
#endif
      }

    span = 2*(ncode-1);
    for (L = 1; L < HUFF_CUTOFF; L++)
      { j = 0;
        for (n = 0; n < span; n++)
          { if (matrix[L][n])
              leng[j++] += 1;
          }
        span = 2*(span-j);
      }
    for (n = 0; n < span; n++)
      leng[n] += 1;

#ifdef DEBUG
    fprintf(stderr,"\nBack Trace:\n");
    span = 2*(ncode-1);
    for (L = 1; L < HUFF_CUTOFF; L++)
      { j = 0;
        fprintf(stderr,"  Row %2d:",L);
        for (n = 0; n < span; n++)
          { if (matrix[L][n])
              j += 1;
            fprintf(stderr," %c",matrix[L][n]?'*':'+');
          }
        fprintf(stderr,"\n");
        span = 2*(span-j);
      }
    fprintf(stderr,"  Length:");
    for (n = 0; n < ncode; n++)
      fprintf(stderr," %d",leng[n]);
    fprintf(stderr,"\n");
#endif
  }

  { int    n, llen;
    uint16 lbits;

    llen  = leng[0];
    lbits = bits[0] = (1 << llen) - 1;
    for (n = 1; n < ncode; n++)
      { while ((lbits & 0x1) == 0)
          { lbits >>= 1;
            llen -= 1;
          }
        lbits -= 1;
        while (llen < leng[n])
          { lbits = (lbits << 1) | 0x1;
            llen += 1;
          }
        bits[n] = lbits;
      }

#ifdef DEBUG
    { int j;

      fprintf(stderr,"\nCodes:\n");
      for (n = 0; n < ncode; n++)
        { fprintf(stderr,"   %3d: %2d ",code[n],leng[n]);
          for (j = leng[n]-1; j >= 0; j--)
            fprintf(stderr,"%x",(bits[n]>>j)&0x1);
          fprintf(stderr,"\n");
        }
    }
#endif
  }

  for (i = 0; i < 256; i++)
    { lens[i] = 0;
      bitv[i] = 0;
    }

  for (i = 0; i < ncode; i++)
    { lens[code[i]] = leng[i];
      bitv[code[i]] = bits[i];
    }

  { int    j, powr;    //  Fill in a decoder table giving the next Huffman code
    uint16 base;       //    that is a prefix of the next 16 bits

    for (i = 0; i < 256; i++)
      { if (lens[i] > 0)
          { base = (bitv[i] << (16-lens[i]));
            powr = (1 << (16-lens[i]));
            for (j = 0; j < powr; j++)
              look[base+j] = i;
          }
      }
  }

  if (partial)
    { v->esc_code = ecode;
      v->esc_len  = lens[ecode];
      lens[ecode] = 0;
    }
  else
    v->esc_code = -1;
  v->state = CODED_WITH;
}

  //  For debug, give a nice print out of the distribution histogram (if present)
  //     and the Huffman codec

void vcPrint(OneCodec *vc, FILE *to)
{ _OneCodec *v = (_OneCodec *) vc;

  uint64  total_bits, ucomp_bits, count;
  uint16  mask, code, *bits;
  uint64 *hist;
  uint8  *lens;
  int     clen;
  int     hashist;
  int     i, k;

  if (vc == DNAcodec)
    { fprintf(to,"    DNAcompressor\n");
      return;
    }

  if (v->state < CODED_WITH) die("vcPrint: Compressor has no codec");
  hashist = (v->state == CODED_WITH);

  bits = v->codebits;
  lens = v->codelens;
  hist = v->hist; // only needed if hashist, but compiler warning if assignment is conditional
      
  if (hashist)
    { total_bits = 0;
      ucomp_bits = 0;

      count = 0;
      for (i = 0; i < 256; i++)
        count += hist[i];

      fprintf(to,"\nHistogram:\n");
      for (i = 0; i < 256; i++)
        if (hist[i] > 0)
          { if (isprint(i))
              fprintf(to,"      %c: %12llu %5.1f%%\n",i,hist[i],(hist[i]*100.)/count);
            else
              fprintf(to,"    %3d: %12llu %5.1f%%\n",i,hist[i],(hist[i]*100.)/count);
          }
    }

  fprintf(to,"\nCode Table:\n");
  for (i = 0; i < 256; i++)
    { clen = lens[i];
      if (i == v->esc_code)
        clen = v->esc_len;
      if (clen > 0)
        { mask = (1 << clen);
          code = bits[i];
          if (isprint(i))
            fprintf(to,"   %c: %2d ",i,clen);
          else
            fprintf(to," %3d: %2d ",i,clen);
          for (k = 0; k < clen; k++)
            { mask >>= 1;
              if (code & mask)
                fprintf(to,"1");
              else
                fprintf(to,"0");
            }
          if (i == v->esc_code)
            fprintf(to," ***\n");
          else
            { fprintf(to,"\n");
              if (hashist)
                { total_bits += clen*hist[i];
                  ucomp_bits += (hist[i]<<3);
                }
            }
        }
    }
  if (hashist)
    fprintf(to,"\nTotal Bytes = %llu (%.2f%%)\n",(total_bits-1)/8+1,(100.*total_bits)/ucomp_bits);
}


/*******************************************************************************************
 *
 *  Read and Write Huffman Schemes (actually just (de)serialize)
 *
 ********************************************************************************************/

  //  Maximum # of bytes in a serialized compressor code

int vcMaxSerialSize()
{ return (257 + 2*sizeof(int) + 256*sizeof(uint16)); }

  //  Code the compressor into blob 'out' and return number of bytes in the code

int vcSerialize(OneCodec *vc, void *out)
{ _OneCodec *v = (_OneCodec *) vc;
  
  int     i;
  uint16 *bits;
  uint8  *lens, *o;

  if (vc == DNAcodec)
    return (0);

  if (v->state < CODED_WITH) die("vcWrite: Compressor does not have a codec");

  lens = v->codelens;
  bits = v->codebits;
  o    = (uint8 *) out;

  //  Only need to record endian, escape code, code lengths, and codes for those
  //    with non-zero length

  *o++ = v->isbig;
  memcpy(o,&(v->esc_code),sizeof(int));
  o += sizeof(int);
  memcpy(o,&(v->esc_len),sizeof(int));
  o += sizeof(int);
  for (i = 0; i < 256; i++)
    { *o++ = lens[i];
      if (lens[i] > 0 || i == v->esc_code)
        { memcpy(o,bits+i,sizeof(uint16));
          o += sizeof(uint16);
        }
    }
  return (o - (uint8 *) out);
}

  //  Create a compressor object from the serialized code in blob 'in'.
  //    The compressor does not have the original histogram from which
  //    its codec was created.  If the endian of the current machine and
  //    the one that serialized the compressor don't match, then all relevant
  //    items are byte-flipped.

OneCodec *vcDeserialize(void *in)
{ _OneCodec *v;

  char    *look;
  uint8   *lens, *ip;
  uint16  *bits, base;
  int      i, j, powr;

  v = (_OneCodec *) malloc(sizeof(_OneCodec));
  if (v == NULL) die("vcRead: Could not allocate compressor");

  v->state = CODED_READ;
  lens = v->codelens;
  bits = v->codebits;
  look = v->lookup;
  ip   = (uint8 *) in;

  { uint32 t;
    uint8 *b;

    t = 1;
    b = (uint8 *) (&t);
    v->isbig = (b[0] == 0);
  }

  if (v->isbig != *ip++)  // If endians out and in don't match then flip item bytes as needed
    { FLIP32(ip)
      memcpy(&(v->esc_code),ip,sizeof(int));
      ip += sizeof(int);
      FLIP32(ip)
      memcpy(&(v->esc_len),ip,sizeof(int));
      ip += sizeof(int);
      for (i = 0; i < 256; i++)
        { lens[i] = *ip++;
          if (lens[i] > 0 || i == v->esc_code)
            { FLIP16(ip)
              memcpy(bits+i,ip,sizeof(uint16));
              ip += sizeof(uint16);
            }
          else
            bits[i] = 0;
        }
    }
  else
    { memcpy(&(v->esc_code),ip,sizeof(int));
      ip += sizeof(int);
      memcpy(&(v->esc_len),ip,sizeof(int));
      ip += sizeof(int);
      for (i = 0; i < 256; i++)
        { lens[i] = *ip++;
          if (lens[i] > 0 || i == v->esc_code)
            { memcpy(bits+i,ip,sizeof(uint16));
              ip += sizeof(uint16);
            }
          else
            bits[i] = 0;
        }
    }

  if (v->esc_code >= 0)
    lens[v->esc_code] = v->esc_len;
  for (i = 0; i < 256; i++)
    { if (lens[i] > 0)
        { base = (bits[i] << (16-lens[i]));
          powr = (1 << (16-lens[i]));
          for (j = 0; j < powr; j++)
            look[base+j] = i;
        }
    }
  if (v->esc_code >= 0)
    lens[v->esc_code] = 0;

  return ((OneCodec *) v);
}


/*******************************************************************************************
 *
 *  Encoders and Decoders
 *
 ********************************************************************************************/

static uint8 Number[128] =
    { 0, 1, 2, 3, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
    };

  //  Compress DNA into 2-bits per base
  //  Richard switched to little-endian December 2022 - big-endian remains in comments
  //  should detect endianness and check and/or switch

int Compress_DNA(int len, char *s, char *t)
{ int    i, j;
  uint8 *s0, *s1, *s2, *s3;

  s0 = (uint8 *) s;
  s1 = s0+1;
  s2 = s1+1;
  s3 = s2+1;

  len -= 3;
  for (i = j = 0; i < len; i += 4)
    t[j++] = Number[s0[i]] | (Number[s1[i]] << 2) | (Number[s2[i]] << 4) | (Number[s3[i]] << 6) ;
      // (Number[s0[i]] << 6) | (Number[s1[i]] << 4) | (Number[s2[i]] << 2) | Number[s3[i]];
  switch (i-len)
  { case 0:
      t[j++] = Number[s0[i]] | (Number[s1[i]] << 2) | (Number[s2[i]] << 4) ;
	// (Number[s0[i]] << 6) | (Number[s1[i]] << 4) | (Number[s2[i]] << 2);
      break;
    case 1:
      t[j++] = Number[s0[i]] | (Number[s1[i]] << 2) ;
        // (Number[s0[i]] << 6) | (Number[s1[i]] << 4);
      break;
    case 2:
      t[j++] = Number[s0[i]] ;
        // (Number[s0[i]] << 6);
      break;
    default:
      break;
  }

  return ((len+3)<<1);
}

  //  Encode ibytes[0..ilen) according to compressor vc and place in obytes
  //  Return the # of bits used.

int vcEncode(OneCodec *vc, int ilen, char *ibytes, char *obytes)
{ _OneCodec *v = (_OneCodec *) vc;

  uint64  c, ocode, *ob;
  int     n, k, rem, tbits, ibits, esc, elen;
  uint8  *clens, x, *bcode, *bb;
  uint16 *cbits;

  if (vc == DNAcodec)
    return (Compress_DNA(ilen,ibytes,obytes));

  if (v->state < CODED_WITH) die("vcEncode: Compressor does not have a codec");

  esc   = v->esc_code;
  elen  = v->esc_len;
  clens = v->codelens;
  cbits = v->codebits;
  ibits = (ilen << 3);
  bcode = (uint8 *) &ocode;

#define OCODE(L,C)				\
{ rem -= L;					\
  if (rem <= 0)					\
    { ocode |= (C >> (-rem));			\
      *ob++ = ocode;				\
      if (rem < 0)				\
        { rem   += 64;				\
          ocode = (C << rem);			\
        }					\
      else					\
        { rem   = 64;				\
          ocode = 0;				\
        }					\
    } 						\
  else						\
    ocode |= (C << rem);			\
}

  ob    = (uint64 *) obytes;
  tbits = 2;
  rem   = 62;
  if (v->isbig)
    ocode = 0x4000000000000000llu;
  else
    ocode = 0;
  for (k = 0; k < ilen; k++)
    { x = ibytes[k];
      n = clens[x];
      if (n == 0)
        { if (esc < 0) die("Compression lib: No code for %c(%x) and no escape code",x,x);
          c = cbits[esc];
          tbits += 8+elen;
          if (tbits > ibits)
            break;
          OCODE(elen,c);
          c = x;
          OCODE(8,c);
        }
      else
        { tbits += n;
          if (tbits > ibits)
            break;
          c = cbits[x];
          OCODE(n,c);
        }
    }
  
  if (k < ilen)
    { *obytes = 0xff;
      memcpy(obytes+1,ibytes,ilen);
      return (ibits+8);
    }

  bb = (uint8 *) ob;
  if (v->isbig)
    { rem = ((71-rem)>>3);
      for (k = 0; k < rem; k++)
        *bb++ = bcode[k];
    }
  else
    { rem = 7 - ((63-rem)>>3);
      for (k = 7; k >= rem; k--)
        *bb++ = bcode[k];
    }

  if (tbits >= 64 && !v->isbig)
    { x = obytes[7];
      obytes[7] = obytes[0];
      obytes[0] = x;
    }

  return (tbits);
}

  //  Uncompress read from 2-bits per base into [0-3] per byte representation

static char Base[4] = { 'a', 'c', 'g', 't' };

int Uncompress_DNA(char *s, int len, char *t)
{ int   i, tlen, byte;
  char *t0, *t1, *t2, *t3;

  t0 = t;
  t1 = t0+1;
  t2 = t1+1;
  t3 = t2+1;

  tlen = len-3;
  for (i = 0; i < tlen; i += 4)
    { byte = *s++;
      t0[i] = Base[byte & 0x3];        // Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 2) & 0x3]; // Base[(byte >> 4) & 0x3];
      t2[i] = Base[(byte >> 4) & 0x3]; // Base[(byte >> 2) & 0x3];
      t3[i] = Base[(byte >> 6) & 0x3]; // Base[byte & 0x3];
    }

  switch (i-tlen)
  { case 0:
      byte = *s++;
      t0[i] = Base[byte & 0x3];        // Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 2) & 0x3]; // Base[(byte >> 4) & 0x3];
      t2[i] = Base[(byte >> 4) & 0x3]; // Base[(byte >> 2) & 0x3];
      break;
    case 1:
      byte = *s++;
      t0[i] = Base[byte & 0x3];        // Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 2) & 0x3]; // Base[(byte >> 4) & 0x3];
      break;
    case 2:
      byte = *s++;
      t0[i] = Base[byte & 0x3];        // Base[(byte >> 6) & 0x3];
      break;
    default:
      break;
  }

  return (len);
}

  //  Decode ilen bits in ibytes, into obytes according to vc's codec
  //  Return the number of bytes decoded.

int vcDecode(OneCodec *vc, int ilen, char *ibytes, char *obytes)
{ _OneCodec *v = (_OneCodec *) vc;

  char   *look;
  uint8  *lens, *q;
  uint64  icode, ncode, *p;
  int     rem, nem;
  uint8   c, *o;
  int     n, k, elen, inbig, esc;

  if (vc == DNAcodec)
    return (Uncompress_DNA(ibytes,ilen>>1,obytes));

  if (v->state < CODED_WITH) die("vcDecode: Compressor does not have a codec");

  if (*((uint8 *) ibytes) == 0xff)
    { int olen = (ilen>>3)-1;
      memcpy(obytes,ibytes+1,olen);
      return (olen);
    }

  p = (uint64 *) ibytes;

  inbig = (*ibytes & 0x40);
  if (!inbig && ilen >= 64)
    { uint8 x = ibytes[7];
      ibytes[7] = ibytes[0];
      ibytes[0] = x;
    }

  if (inbig != v->isbig)
    { q = (uint8 *) ibytes;
      for (k = 64; k <= ilen; k += 64)
        { FLIP64(q)
          q += 8;
        }
    }

  lens = v->codelens;
  look = v->lookup;
  esc  = v->esc_code;
  elen = v->esc_len;

#define GET(n)						\
  ilen  -= n;						\
  icode <<= n;						\
  rem   -= n;						\
  while (rem < 16)					\
    { int z = 64-rem;					\
      icode |= (ncode >> rem);				\
      if (nem > z)					\
        { nem -= z;					\
          ncode <<= z;					\
          rem = 64;					\
          break;					\
        }						\
      else						\
        { rem += nem; 					\
          if (rem >= ilen)				\
            break;					\
          else if (ilen-rem < 64)			\
            { nem = ilen-rem;				\
              q = (uint8 *) p;				\
              ncode = 0;				\
              for (k = 0; k < nem; k += 8)		\
                ncode |= (((uint64) (*q++)) << (56-k));	\
            }						\
          else						\
            { ncode = *p++;				\
              nem   = 64;				\
            }						\
	}						\
    }
 
  if (ilen < 64)
    { q = (uint8 *) ibytes;
      icode = 0;
      for (k = 0; k < ilen; k += 8)
        icode |= (((uint64) (*q++)) << (56-k));
    }
  else
    icode = *p++;
  o = (uint8 *) obytes;
  icode <<= 2;
  ilen -= 2;
  rem   = 62;
  if (rem > ilen)
    rem = ilen;
  ncode = 0;
  nem   = 0;
  while (ilen > 0)
    { c = look[icode >> 48];
      if (c == esc)
        { GET(elen)
          c = (icode >> 56);
          GET(8);
        }
      else
        { n = lens[(int) c];
          GET(n)
        }
      *o++ = c;
    }

  return (o - (uint8 *) obytes);
}

//////////////////////////////////////////////////////////////////////////////////////
//
// integer compression for write/read of fields
//
// top bit of first byte: number is negative
// second bit: one-byte: next six bits give number (make negative if top bit set)
// third bit: two-byte: next 13 bits give number (make negative if top bit set)
// if second and third bits are not set, remaining 5 bits give number of bytes to read

static inline int intGet (unsigned char *u, I64 *pval)
{
  switch (u[0] >> 5)
    {
    case 2: case 3: // single byte positive
      *pval = (I64) (u[0] & 0x3f) ; return 1 ;
    case 6: case 7: // single byte negative
      *pval =  (I64) u[0] | 0xffffffffffffff00 ; return 1 ;
    case 1: // two bytes positive
      *pval = (I64) (u[0] & 0x1f) << 8 | (I64)u[1] ; return 2 ;
      *pval = - ((I64) (u[0] & 0x1f) << 8 | (I64)u[1]) ; return 2 ;
    case 0:
      switch (u[0] & 0x07)
	{
	case 0: die ("int packing error") ; break ;
	case 1: *pval = *(I64*)(u+1) & 0x0000000000ffff ; return 3 ;
	case 2: *pval = *(I64*)(u+1) & 0x00000000ffffff ; return 4 ;
	case 3: *pval = *(I64*)(u+1) & 0x000000ffffffff ; return 5 ;
	case 4: *pval = *(I64*)(u+1) & 0x0000ffffffffff ; return 6 ;
	case 5: *pval = *(I64*)(u+1) & 0x00ffffffffffff ; return 7 ;
	case 6: *pval = *(I64*)(u+1) & 0xffffffffffffff ; return 8 ;
	case 7: *pval = *(I64*)(u+1) ; return 9 ;
	}
      break ;
    case 4:
      switch (u[0] & 0x07)
	{
	case 0: die ("int packing error") ; break ;
	case 1: *pval = *(I64*)(u+1) | 0xffffffffffff0000 ; return 3 ;
	case 2: *pval = *(I64*)(u+1) | 0xffffffffff000000 ; return 4 ;
	case 3: *pval = *(I64*)(u+1) | 0xffffffff00000000 ; return 5 ;
	case 4: *pval = *(I64*)(u+1) | 0xffffff0000000000 ; return 6 ;
	case 5: *pval = *(I64*)(u+1) | 0xffff000000000000 ; return 7 ;
	case 6: *pval = *(I64*)(u+1) | 0xff00000000000000 ; return 8 ;
	case 7: *pval = *(I64*)(u+1) ; return 9 ;
	}
      break ;
    }
  return 0 ; // shouldn't get here, but needed for compiler happiness
}

static inline int intPut (unsigned char *u, I64 val)
{
  if (val >= 0)
    { if (     !(val & 0xffffffffffffffc0)) { *u = val | 0x40 ;  return 1 ; }
      else if (!(val & 0xffffffffffffe000)) { *u++ = (val >> 8) | 0x20 ; *u = val & 0xff ; return 2 ; }
      else if (!(val & 0xffffffffffff0000)) { *u++ = 1 ; *(I64*)u = val ; return 3 ; }
      else if (!(val & 0xffffffffff000000)) { *u++ = 2 ; *(I64*)u = val ; return 4 ; }
      else if (!(val & 0xffffffff00000000)) { *u++ = 3 ; *(I64*)u = val ; return 5 ; }
      else if (!(val & 0xffffff0000000000)) { *u++ = 4 ; *(I64*)u = val ; return 6 ; }
      else if (!(val & 0xffff000000000000)) { *u++ = 5 ; *(I64*)u = val ; return 7 ; }
      else if (!(val & 0xff00000000000000)) { *u++ = 6 ; *(I64*)u = val ; return 8 ; }
      else                                  { *u++ = 7 ; *(I64*)u = val ; return 9 ; }
    }
  else
    { if (     !(~val & 0xffffffffffffffc0)) { *u = val | 0x40 ;  return 1 ; }
      //     else if (!(~val & 0xffffffffffffe000)) { *u++ = (val >> 8) | 0x20 ; *u = val & 0xff ; return 2 ; }
      else if (!(~val & 0xffffffffffff0000)) { *u++ = 0x81 ; *(I64*)u = val ; return 3 ; }
      else if (!(~val & 0xffffffffff000000)) { *u++ = 0x82 ; *(I64*)u = val ; return 4 ; }
      else if (!(~val & 0xffffffff00000000)) { *u++ = 0x83 ; *(I64*)u = val ; return 5 ; }
      else if (!(~val & 0xffffff0000000000)) { *u++ = 0x84 ; *(I64*)u = val ; return 6 ; }
      else if (!(~val & 0xffff000000000000)) { *u++ = 0x85 ; *(I64*)u = val ; return 7 ; }
      else if (!(~val & 0xff00000000000000)) { *u++ = 0x86 ; *(I64*)u = val ; return 8 ; }
      else                                   { *u++ = 0x87 ; *(I64*)u = val ; return 9 ; }
    }
}

static inline I64 ltfRead (FILE *f)
{
  unsigned char u[16] ;
  I64 val = 0 ;

  u[0] = getc (f) ;
  if (u[0] & 0x40)
    { if (u[0] & 0x80) // negative
	val = (char) u[0] ;
      else
	val = (I64) (u[0] & 0x3f) ;
    }
    // { intGet (u, &val) ;
    //   printf ("read %d n 1 u %02x\n", (int)val, u[0]) ;
    // }
  else if (u[0] & 0x20)
    { u[1] = getc (f) ; intGet (u, &val) ;
      //      printf ("read %d n 2 u %02x %02x\n", (int)val, u[0], u[1]) ;
    }
  else
    { int n = 1 + (u[0] & 0x0f) ;
      unsigned char *v = &u[1] ;
      while (n--) *v++ = getc(f) ;
      n = intGet (u, &val) ;
      //      printf ("read %d n %d u", (int)val, n) ;
      //      { int i ; for (i = 0 ; i< n ; ++i) printf (" %02x", u[i]) ; putchar ('\n') ; }
    }

  return val ;
}

static inline int ltfWrite (I64 x, FILE *f)
{
  unsigned char u[16] ;
  int n = intPut (u, x) ;

  //  printf ("write %d n %d u", (int)x, n) ;
  //  { int i ; for (i = 0 ; i< n ; ++i) printf (" %02x", u[i]) ; putchar ('\n') ; }

  fwrite (u, 1, n, f) ;
  return n ;
}

#if defined(TEST_LTF) || defined(TEST_INT)

// these are the original routines from James Bonfield on which this is based
// incorporated here for credit and for comparisons in the TEST functionality

/***********************************************************************************
 *
 *    LTF encoding for integers
 *    adapted from htslib/cram/cram_io.h with copyright statement:

Copyright (c) 2012-2019 Genome Research Ltd.
Author: James Bonfield <jkb@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 *
 **********************************************************************************/

/* 64-bit itf8 variant */

static inline int ltf8_put(char *cp, int64_t val) {
    unsigned char *up = (unsigned char *)cp;
    if        (!(val & ~((1LL<<7)-1))) {
        *up = val;
        return 1;
    } else if (!(val & ~((1LL<<(6+8))-1))) {
        *up++ = (val >> 8 ) | 0x80;
        *up   = val & 0xff;
        return 2;
    } else if (!(val & ~((1LL<<(5+2*8))-1))) {
        *up++ = (val >> 16) | 0xc0;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 3;
    } else if (!(val & ~((1LL<<(4+3*8))-1))) {
        *up++ = (val >> 24) | 0xe0;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 4;
    } else if (!(val & ~((1LL<<(3+4*8))-1))) {
        *up++ = (val >> 32) | 0xf0;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 5;
    } else if (!(val & ~((1LL<<(2+5*8))-1))) {
        *up++ = (val >> 40) | 0xf8;
        *up++ = (val >> 32) & 0xff;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 6;
    } else if (!(val & ~((1LL<<(1+6*8))-1))) {
        *up++ = (val >> 48) | 0xfc;
        *up++ = (val >> 40) & 0xff;
        *up++ = (val >> 32) & 0xff;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 7;
    } else if (!(val & ~((1LL<<(7*8))-1))) {
        *up++ = (val >> 56) | 0xfe;
        *up++ = (val >> 48) & 0xff;
        *up++ = (val >> 40) & 0xff;
        *up++ = (val >> 32) & 0xff;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 8;
    } else {
        *up++ = 0xff;
        *up++ = (val >> 56) & 0xff;
        *up++ = (val >> 48) & 0xff;
        *up++ = (val >> 40) & 0xff;
        *up++ = (val >> 32) & 0xff;
        *up++ = (val >> 24) & 0xff;
        *up++ = (val >> 16) & 0xff;
        *up++ = (val >> 8 ) & 0xff;
        *up   = val & 0xff;
        return 9;
    }
}

static inline int ltf8_get(char *cp, int64_t *val_p) {
    unsigned char *up = (unsigned char *)cp;

    if (up[0] < 0x80) {
        *val_p =   up[0];
        return 1;
    } else if (up[0] < 0xc0) {
        *val_p = (((uint64_t)up[0]<< 8) |
                   (uint64_t)up[1]) & (((1LL<<(6+8)))-1);
        return 2;
    } else if (up[0] < 0xe0) {
        *val_p = (((uint64_t)up[0]<<16) |
                  ((uint64_t)up[1]<< 8) |
                   (uint64_t)up[2]) & ((1LL<<(5+2*8))-1);
        return 3;
    } else if (up[0] < 0xf0) {
        *val_p = (((uint64_t)up[0]<<24) |
                  ((uint64_t)up[1]<<16) |
                  ((uint64_t)up[2]<< 8) |
                   (uint64_t)up[3]) & ((1LL<<(4+3*8))-1);
        return 4;
    } else if (up[0] < 0xf8) {
        *val_p = (((uint64_t)up[0]<<32) |
                  ((uint64_t)up[1]<<24) |
                  ((uint64_t)up[2]<<16) |
                  ((uint64_t)up[3]<< 8) |
                   (uint64_t)up[4]) & ((1LL<<(3+4*8))-1);
        return 5;
    } else if (up[0] < 0xfc) {
        *val_p = (((uint64_t)up[0]<<40) |
                  ((uint64_t)up[1]<<32) |
                  ((uint64_t)up[2]<<24) |
                  ((uint64_t)up[3]<<16) |
                  ((uint64_t)up[4]<< 8) |
                   (uint64_t)up[5]) & ((1LL<<(2+5*8))-1);
        return 6;
    } else if (up[0] < 0xfe) {
        *val_p = (((uint64_t)up[0]<<48) |
                  ((uint64_t)up[1]<<40) |
                  ((uint64_t)up[2]<<32) |
                  ((uint64_t)up[3]<<24) |
                  ((uint64_t)up[4]<<16) |
                  ((uint64_t)up[5]<< 8) |
                   (uint64_t)up[6]) & ((1LL<<(1+6*8))-1);
        return 7;
    } else if (up[0] < 0xff) {
        *val_p = (((uint64_t)up[1]<<48) |
                  ((uint64_t)up[2]<<40) |
                  ((uint64_t)up[3]<<32) |
                  ((uint64_t)up[4]<<24) |
                  ((uint64_t)up[5]<<16) |
                  ((uint64_t)up[6]<< 8) |
                   (uint64_t)up[7]) & ((1LL<<(7*8))-1);
        return 8;
    } else {
        *val_p = (((uint64_t)up[1]<<56) |
                  ((uint64_t)up[2]<<48) |
                  ((uint64_t)up[3]<<40) |
                  ((uint64_t)up[4]<<32) |
                  ((uint64_t)up[5]<<24) |
                  ((uint64_t)up[6]<<16) |
                  ((uint64_t)up[7]<< 8) |
                   (uint64_t)up[8]);
        return 9;
    }
}

#include <sys/resource.h>
#ifndef RUSAGE_SELF     /* to prevent "RUSAGE_SELF redefined" gcc warning, fixme if this is more intricate */
#define RUSAGE_SELF 0
#endif

void timeUpdate (FILE *f)
{
  static bool isFirst = 1 ;
  static struct rusage rOld, rFirst ;
  struct rusage rNew ;
  int secs, usecs ;

  getrusage (RUSAGE_SELF, &rNew) ;
  if (!isFirst)
    { secs = rNew.ru_utime.tv_sec - rOld.ru_utime.tv_sec ;
      usecs =  rNew.ru_utime.tv_usec - rOld.ru_utime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "user\t%d.%06d", secs, usecs) ;
      secs = rNew.ru_stime.tv_sec - rOld.ru_stime.tv_sec ;
      usecs =  rNew.ru_stime.tv_usec - rOld.ru_stime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "\tsystem\t%d.%06d", secs, usecs) ;
      fprintf (f, "\tmax_RSS\t%ld", rNew.ru_maxrss - rOld.ru_maxrss) ;
      fputc ('\n', f) ;
    }
  else
    { rFirst = rNew ;
      isFirst = false ;
    }

  rOld = rNew ;
}

int main (int argc, char *argv[])
{
  I64   i, j, x, n, tot, mod ;
  FILE *f ;
  static unsigned char buffer[9*(1<<20)] ;

  if (argc < 3) die ("usage: ./test <start> <n> [mod]") ;
  
//  { int   t = 1;
//    char *b = (char *) (&t);
//    if (*b == 0) printf ("bigEndian\n") ; else printf ("smallEndian\n") ;
//  }

  timeUpdate (0) ;
  
  x = atoi(argv[1]) ;
  n = atoi(argv[2]) ;
  if (argc == 4) mod = atoi(argv[3]) ; else mod = 0 ;
  tot = 0 ;
#ifdef TEST_INT
  f = fopen ("int.test", "w") ;
#endif
#ifdef TEST_LTF
  f = fopen ("ltf.test", "w") ;
#endif
  if (argc > 4)
    for (i = 0 ; i < n ; ++i) { tot += ltfWrite (x++, f) ; if (x == mod) x = 0 ; }
  else
    { while (n)
	{ int m = (n > 1<<20) ? 1<<20 : n ;
	  unsigned char *u = buffer ;
	  for (i = 0 ; i < m ; ++i)
	    {
#ifdef TEST_INT
	      u += intPut (u, x++) ;
#endif
#ifdef TEST_LTF
	      u += ltf8_put ((char*)u, x++) ;
#endif
	      if (x == mod) x = 0 ;
	    }
	  tot += (u-buffer) ;
	  fwrite (buffer, 1, (u-buffer), f) ;
	  n -= m ;
	}
    }
  fclose (f) ;
  printf ("wrote %lld bytes: ", tot) ;
  timeUpdate (stdout) ;

  x = atoi(argv[1]) ;
  n = atoi(argv[2]) ;
  tot = 0 ;
#ifdef TEST_INT
  f = fopen ("int.test", "r") ;
#endif
#ifdef TEST_LTF
  f = fopen ("ltf.test", "r") ;
#endif
  if (argc > 4)
    {  for (j = 0 ; j < n ; ++j)
	if ((i = ltfRead (f)) != j)
	  die ("ltf wrote %d read %d", (int)j, (int)i) ;
    }
  else
    { unsigned char *u = buffer, *v = buffer ; // u is current pos, v is end of section read in
      while (n)
	{ int m = (n > 1<<20) ? 1<<20 : n ;
//	  printf ("attempting to read %d chars: n %d (v-u) %d\n", (int)(9*m - (v-u)), m, (int)(v-u)) ;
	  unsigned char *v0 = v, *u0 = u ;
	  v += fread (v, 1, 9*m - (v-u), f) ;
//	  printf ("  read %d chars\n", (int)(v-v0)) ;
	  for (i = 0 ; i < m ; ++i)
#ifdef TEST_INT
	    u += intGet (u, &j) ;
#endif
#ifdef TEST_LTF
	    u += ltf8_get ((char*)u, &j) ;
#endif
//	  printf ("  intGet m %d ints used %d chars\n", m, (int)(u-u0)) ;
	  n -= m ;
	  m = v - u ;
//	  printf ("  memmove %d\n", m) ;
	  memmove (buffer, u, m) ;
	  tot += (u-buffer) ;
	  v -= (u-buffer) ; u = buffer ;
	}
    }
  fclose (f) ;

  printf ("read  %lld bytes: ", tot) ;
  timeUpdate (stdout) ;
}
#endif // TEST_LTF

/***********************************************************************************
 *
 *    UTILITIES: die(), memory allocation
 *    adapted from Richard's utilities
 *
 **********************************************************************************/

static void die(char *format, ...)
{ va_list args;

  va_start (args, format);
  fprintf (stderr, "FATAL ERROR: ");
  vfprintf (stderr, format, args);
  fprintf (stderr, "\n");
  va_end (args);
  exit (-1);
}

static I64 nAlloc = 0;
static I64 totalAlloc = 0;

static void *myalloc(size_t size)
{ void *p;

  p = malloc(size);
  if (p == NULL && size != 0 )
    die("ONElib myalloc failure requesting %d bytes - totalAlloc %lld", size, totalAlloc);
  nAlloc     += 1;
  totalAlloc += size;
  return (p);
}

static void *mycalloc(size_t number, size_t size)
{ void *p;

  p = calloc(number,size);
  if (p == NULL && size > 0) die("mycalloc failure requesting %d objects of size %d", number, size);
  nAlloc     += 1;
  totalAlloc += size*number;
  return p;
}

static void *mydup(size_t n, void *x, size_t size)
{
  void *p = myalloc (n*size);
  memcpy (p, x, n*size);
  return p;
}

/********************* end of file ***********************/
