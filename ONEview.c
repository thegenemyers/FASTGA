/*  File: ONEview.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 30 21:54 2024 (rd109)
 * * May 15 02:26 2024 (rd109): incorporate rd utilities so stand alone
 * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "ONElib.h"

#include <assert.h>

#include <string.h>		/* strcmp etc. */
#include <stdlib.h>		/* for exit() */
#include <stdarg.h>             /* for variable length argument lists */

// forward declarations of utilities at end of file from RD's utils.[ch]

void die (char *format, ...) ;
char *commandLine (int argc, char **argv) ;

void *myalloc (size_t size) ;
void *mycalloc (size_t number, size_t size) ;
#define	new(n,type)	(type*)myalloc((n)*sizeof(type))
#define	new0(n,type)	(type*)mycalloc((n),sizeof(type))
#define resize(x,nOld,nNew,T) { T* z = new((nNew),T) ; if (nOld < nNew) memcpy(z,x,(nOld)*sizeof(T)) ; else memcpy(z,x,(nNew)*sizeof(T)) ; free(x) ; x = z ; }

void timeUpdate (FILE *f) ;	/* print time usage since last call to file */
void timeTotal (FILE *f) ;	/* print full time usage since first call to timeUpdate */

// end of utils declarations

typedef struct IndexListStruct {
  I64 i0, iN ;
  struct IndexListStruct *next ;
} IndexList ;

static IndexList *parseIndexList (char *s)
{
  IndexList *ol, *ol0 = ol = new0 (1, IndexList) ;
  while (*s)
    { while (*s >= '0' && *s <= '9') ol->i0 = ol->i0*10 + (*s++ - '0') ;
      if (*s == '-')
	{ ++s ; while (*s >= '0' && *s <= '9') ol->iN = ol->iN*10 + (*s++ - '0') ;
	  if (ol->iN <= ol->i0) die ("end index %lld <= start index %lld", ol->iN, ol->i0) ;
	}
      else
	ol->iN = ol->i0 + 1 ;
      if (*s == ',') { ol->next = new0 (1, IndexList) ; ol = ol->next ; ++s ; }
      else if (*s) die ("unrecognised character %c at %s in object list\n", *s, s) ;
    }
  return ol0 ; 
}

static void transferLine (OneFile *vfIn, OneFile *vfOut, size_t *fieldSize)
{ memcpy (vfOut->field, vfIn->field, fieldSize[(int)vfIn->lineType]) ;
  oneWriteLine (vfOut, vfIn->lineType, oneLen(vfIn), oneString(vfIn)) ;
  char *s = oneReadComment (vfIn) ; if (s) oneWriteComment (vfOut, "%s", s) ;
}

int main (int argc, char **argv)
{
  I64 i ;
  char *fileType = 0 ;
  char *outFileName = "-" ;
  char *schemaFileName = 0 ;
  bool  isNoHeader = false, isHeaderOnly = false, isWriteSchema = false, 
    isBinary = false, isVerbose = false ;
  char  indexType = 0 ;
  IndexList *objList = 0 ;
  
  timeUpdate (0) ;

  char *command = commandLine (argc, argv) ;
  --argc ; ++argv ;		/* drop the program name */

  if (!argc)
    { fprintf (stderr, "ONEview [options] onefile\n") ;
      fprintf (stderr, "  -t --type <abc>           file type, e.g. seq, aln - required if no header\n") ;
      fprintf (stderr, "  -S --schema <schemafile>      schema file name for reading file\n") ;
      fprintf (stderr, "  -h --noHeader                 skip the header in ascii output\n") ;
      fprintf (stderr, "  -H --headerOnly               only write the header (in ascii)\n") ;
      fprintf (stderr, "  -s --writeSchema              write a schema file based on this file\n") ;
      fprintf (stderr, "  -b --binary                   write in binary (default is ascii)\n") ;
      fprintf (stderr, "  -o --output <filename>        output file name (default stdout)\n") ;
      fprintf (stderr, "  -i --index T x[-y](,x[-y])*   write specified objects/groups of type T\n") ;
      fprintf (stderr, "  -v --verbose                  write commentary including timing\n") ;
      fprintf (stderr, "index only works for binary files; '-i A 0-10' outputs first 10 objects of type A\n") ;
      exit (0) ;
    }
  
  while (argc && **argv == '-')
    if (!strcmp (*argv, "-t") || !strcmp (*argv, "--type"))
      { fileType = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-S") || !strcmp (*argv, "--schema"))
      { schemaFileName = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-h") || !strcmp (*argv, "--noHeader"))
      { isNoHeader = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-H") || !strcmp (*argv, "--headerOnly"))
      { isHeaderOnly = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-s") || !strcmp (*argv, "--writeSchema"))
      { isWriteSchema = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-b") || !strcmp (*argv, "--binary"))
      { isBinary = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-v") || !strcmp (*argv, "--verbose"))
      { isVerbose = true ; --argc ; ++argv ; }
    else if ((!strcmp (*argv, "-o") || !strcmp (*argv, "--output")) && argc >= 2)
      { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if ((!strcmp (*argv, "-i") || !strcmp (*argv, "--index")) && argc >= 3)
      { indexType = *argv[1] ; objList = parseIndexList (argv[2]) ; argc -= 3 ; argv += 3 ; }
    else die ("unknown option %s - run without arguments to see options", *argv) ;

  if (isBinary) isNoHeader = false ;
  if (isHeaderOnly) isBinary = false ;
    
  if (argc != 1)
    die ("need a single data one-code file as argument") ;

  OneSchema *vs = 0 ;
  if (schemaFileName && !(vs = oneSchemaCreateFromFile (schemaFileName)))
    die ("failed to read schema file %s", schemaFileName) ;
  OneFile *vfIn = oneFileOpenRead (argv[0], vs, fileType, 1) ; /* reads the header */
  if (!vfIn) die ("failed to open one file %s", argv[0]) ;

  if (objList)
    { if (!vfIn->isBinary)
	die ("%s is ascii - you can only access objects and groups by index in binary files", argv[0]) ;
      if (!vfIn->info[(int)indexType])
	die ("requested index type %c is not present in the schema", indexType) ;
      if (!vfIn->info[(int)indexType]->index)
	die ("no index for line type %c", indexType) ;
    }

  if (isWriteSchema)
    { oneFileWriteSchema (vfIn, outFileName) ; }
  else
    { OneFile *vfOut = oneFileOpenWriteFrom (outFileName, vfIn, isBinary, 1) ;
      if (!vfOut) die ("failed to open output file %s", outFileName) ;
      if (!isBinary) // need to copy across the object stats, so they write out into the header
	for (i = 0 ; i < vfIn->nDefn ; ++i)
	  { int k = vfIn->defnOrder[i] ;
	    if (!(k & 0x80) && vfIn->info[k]->stats)
	      { int n = 1 ; OneStat *s ;
		for (s = vfIn->info[k]->stats ; s->type ; ++s) ++n ;
		vfOut->info[k]->stats = new0 (n, OneStat) ;
		memcpy (vfOut->info[k]->stats, vfIn->info[k]->stats, n*sizeof(OneStat)) ;
	      }
	}
	
      if (isNoHeader) vfOut->isNoAsciiHeader = true ; // will have no effect if binary

      if (!isHeaderOnly)
	{ oneAddProvenance (vfOut, "ONEview", "0.0", command) ;
      
	  static size_t fieldSize[128] ;
	  for (i = 0 ; i < 128 ; ++i)
	    if (vfIn->info[i]) fieldSize[i] = vfIn->info[i]->nField*sizeof(OneField) ;
      
	  if (objList)
	    while (objList)
	      { if (!oneGoto (vfIn, indexType, objList->i0))
		  die ("can't locate to object %c %lld", indexType, objList->i0 ) ;
		if (!oneReadLine (vfIn))
		  die ("can't read object %c %lld", indexType, objList->i0) ;
		if (objList->i0 == 0) // write up until 1st object of indexType
		  { while (vfIn->lineType && vfIn->lineType != indexType)
		      { transferLine (vfIn, vfOut, fieldSize) ;
			oneReadLine (vfIn) ;
		      }
		    ++objList->i0 ;
		  }
		bool isInside = true ;
		while (vfIn->lineType && objList->i0 < objList->iN) // lineType 0 is end of file
		  { if (isInside) transferLine (vfIn, vfOut, fieldSize) ;
		    oneReadLine (vfIn) ;
		    if (!vfIn->info[(int)indexType]->contains[(int) vfIn->lineType])
		      isInside = false ;
		    if (vfIn->lineType == indexType)
		      { ++objList->i0 ; isInside = true ; }
		  }
		objList = objList->next ;
	      }
	  else
	    while (oneReadLine (vfIn))
	      transferLine (vfIn, vfOut, fieldSize) ;
	}
      oneFileClose (vfOut) ;
    }
      
  oneFileClose (vfIn) ;
  if (vs) oneSchemaDestroy (vs) ;
  
  free (command) ;
  if (isVerbose) timeTotal (stderr) ;

  exit (0) ;
}


/*********** utilities from RD's utils.[ch] ***************/

void die (char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "FATAL ERROR: ") ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;

  exit (-1) ;
}

char *commandLine (int argc, char **argv)
{
  int i, totLen = 0 ;
  for (i = 0 ; i < argc ; ++i) totLen += 1 + strlen(argv[i]) ;
  char *buf = new (totLen, char) ;
  strcpy (buf, argv[0]) ;
  for (i = 1 ; i < argc ; ++i) { strcat (buf, " ") ; strcat (buf, argv[i]) ; }
  return buf ;
}

long totalAllocated = 0 ;

void *myalloc (size_t size)
{
  void *p = (void*) malloc (size) ;
  if (!p) die ("myalloc failure requesting %d bytes - totalAllocated %ld", size, totalAllocated) ;
  totalAllocated += size ;
  return p ;
}

void *mycalloc (size_t number, size_t size)
{
  void *p = (void*) calloc (number, size) ;
  if (!p) die ("mycalloc failure requesting %d objects of size %d - totalAllocated %ld", number, size, totalAllocated) ;
  totalAllocated += size*number ;
  return p ;
}

/***************** rusage for timing information ******************/

#include <sys/resource.h>
#include <sys/time.h>
#ifndef RUSAGE_SELF     /* to prevent "RUSAGE_SELF redefined" gcc warning, fixme if this is more intricate */
#define RUSAGE_SELF 0
#endif

#ifdef RUSAGE_STRUCTURE_DEFINITIONS
struct rusage {
  struct timeval ru_utime; /* user time used */
  struct timeval ru_stime; /* system time used */
  long ru_maxrss;          /* integral max resident set size */
  long ru_ixrss;           /* integral shared text memory size */
  long ru_idrss;           /* integral unshared data size */
  long ru_isrss;           /* integral unshared stack size */
  long ru_minflt;          /* page reclaims */
  long ru_majflt;          /* page faults */
  long ru_nswap;           /* swaps */
  long ru_inblock;         /* block input operations */
  long ru_oublock;         /* block output operations */
  long ru_msgsnd;          /* messages sent */
  long ru_msgrcv;          /* messages received */
  long ru_nsignals;        /* signals received */
  long ru_nvcsw;           /* voluntary context switches */
  long ru_nivcsw;          /* involuntary context switches */
};

struct timeval {
  time_t       tv_sec;   /* seconds since Jan. 1, 1970 */
  suseconds_t  tv_usec;  /* and microseconds */
} ;
#endif /* RUSAGE STRUCTURE_DEFINITIONS */

static struct rusage rOld, rFirst ;
static struct timeval tOld, tFirst ;

void timeUpdate (FILE *f)
{
  static bool isFirst = 1 ;
  struct rusage rNew ;
  struct timeval tNew ;
  int secs, usecs ;

  getrusage (RUSAGE_SELF, &rNew) ;
  gettimeofday(&tNew, 0) ;
  if (!isFirst)
    { secs = rNew.ru_utime.tv_sec - rOld.ru_utime.tv_sec ;
      usecs =  rNew.ru_utime.tv_usec - rOld.ru_utime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "user\t%d.%06d", secs, usecs) ;
      secs = rNew.ru_stime.tv_sec - rOld.ru_stime.tv_sec ;
      usecs =  rNew.ru_stime.tv_usec - rOld.ru_stime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "\tsystem\t%d.%06d", secs, usecs) ;
      secs = tNew.tv_sec - tOld.tv_sec ;
      usecs =  tNew.tv_usec - tOld.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "\telapsed\t%d.%06d", secs, usecs) ;
      fprintf (f, "\tallocated\t%.2f", totalAllocated/1000000000.0) ;   
      fprintf (f, "\tmax_RSS\t%ld", rNew.ru_maxrss - rOld.ru_maxrss) ;
      fputc ('\n', f) ;
    }
  else
    { rFirst = rNew ;
      tFirst = tNew ;
      isFirst = false ;
    }

  rOld = rNew ;
  tOld = tNew ;
}

void timeTotal (FILE *f) { rOld = rFirst ; tOld = tFirst ; timeUpdate (f) ; }

/********************* end of file ***********************/
