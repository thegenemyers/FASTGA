/*******************************************************************************************
 *
 *  Interface to read and process the contents of a .1aln file.
 *
 *  Author :  Gene Myers
 *  Date   :  October 2025
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <strings.h>
#include "ONEaln.h"
#include "GDB.h"
#include "alncode.h"

#define CHECK_ARGS
#define HALT_ON_ERROR

#define BORDER_MAX 100

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

#ifdef HALT_ON_ERROR

static char *alnMessage = NULL;

#else

static char alnMessage[1000];

#endif

char *alnError()
{ return (alnMessage); }

 //  .1ALN FILE READER

  //  Open the .1aln file 'name' for reading with nthreads threads.  NULL is returned
  //    if there is an error, otherwise the return pointer points to the first element
  //    in an array of nthreads AlnReader's.  This first element is called the master.
  //    If you want to have CIGAR strings or print Alignments or otherwise need access
  //    to the sequence of the source genomes, you must set see_seq to true.  Otherwise
  //    set it to false for better memory efficiency.

typedef struct
  { int  op;
    int  len;
  } Cigar;

typedef struct
  { AlnRecord align;
    int       nthreads;
    GDB      *gpt1, *gpt2;
    GDB       gdb1, gdb2;   //  gdb1 headers and pointers to said, gpt1 == gpt2 if self
    int       tspace;
    int       contig1;      //  the contigs for the current AlnRecord
    int       contig2;
    int       acnt;         //  # of alignments in the .1aln file

    OneFile  *file;         //  OneFile reader
    bool      read;         //  The current 'A' record has been read (advanced through)
    bool      eof;          //  The end of the file has been reached.

    Work_Data *work;
    char      *aseq;
    char      *bseq;

    int        cigar_max;
    Cigar     *cigar;
  } _AlnReader;

AlnReader *alnOpenReader(char *name, int nthreads, bool see_seq)
{ char    *pwd, *root;
  OneFile *input;
  char    *src1, *src2, *cpath;
  char    *last_error;
  int      tspace;
  FILE   **units1, **units2;
  GDB     _gdb1, *gdb1 = &_gdb1;
  GDB     _gdb2, *gdb2 = &_gdb2;
  char    *head, *eptr;
  int64    acnt, tmax;
  int      i;

  AlnReader  *a;
  _AlnReader *b;
  int        *c;
  char       *d;

  last_error   = Error_Buffer;
  Error_Buffer = alnMessage;

  pwd   = PathTo(name);
  root  = Root(name,".1aln");
  input = open_Aln_Read(Catenate(pwd,"/",root,".1aln"),nthreads,&acnt,&tspace,
                        &src1,&src2,&cpath);
  if (input == NULL)
    { EPRINTF("Could not open %s/%s.1aln",pwd,root);
      free(root);
      free(pwd);
      Error_Buffer = last_error;
      EXIT(NULL);
    }
  free(root);
  free(pwd);

  oneStats(input,'T',NULL,&tmax,NULL);

  units1 = units2 = NULL;
  if (see_seq)
    { Skip_Skeleton(input);
      units1 = Get_GDB(gdb1,src1,cpath,nthreads,NULL);
      if (units1 == NULL)
        goto closeout;
    }
  else
    { if (input->lineType == 'g')
        { if (Read_Skeleton(input,src1,gdb1))
            goto closeout;
        }
      else
        { units1 = Get_GDB(gdb1,src1,cpath,0,NULL);
          if (units1 == NULL)
            goto closeout;
        }
    }

  if (src2 != NULL)
    { if (see_seq)
        { Skip_Skeleton(input);
          units2 = Get_GDB(gdb2,src2,cpath,nthreads,NULL);
          if (units2 == NULL)
            { Close_GDB(gdb1);
              goto closeout;
            }
        }
      else
        { if (input->lineType == 'g')
            { if (Read_Skeleton(input,src2,gdb2))
                { Close_GDB(gdb1);
                  goto closeout;
                }
            }
          else
            { units2 = Get_GDB(gdb2,src2,cpath,0,NULL);
              if (units2 == NULL)
                { Close_GDB(gdb1);
                  goto closeout;
                }
            }
        }
    }
  else
    { gdb2   = gdb1;
      units2 = units1;
    }

  head = gdb1->headers;
  for (i = 0; i < gdb1->nscaff; i++)
    { for (eptr = head + gdb1->scaffolds[i].hoff; *eptr != '\0'; eptr++)
        if (isspace(*eptr))
          break;
      *eptr = '\0';
    }

  if (src2 != NULL)
    { head  = gdb2->headers;
      for (i = 0; i < gdb2->nscaff; i++)
        { for (eptr = head + gdb2->scaffolds[i].hoff; *eptr != '\0'; eptr++)
            if (isspace(*eptr))
              break;
          *eptr = '\0';
        }
      free(src2);
    }
  free(src1);
  free(cpath);

  a = (AlnReader *)  Malloc(sizeof(AlnReader)*nthreads,"Opening .1aln file");
  b = (_AlnReader *) Malloc(sizeof(_AlnReader)*nthreads,"Opening .1aln file");
  c = (int *)        Malloc(sizeof(int)*3*nthreads*tmax,"Opening .1aln file");
  if (see_seq)
    d = (char *) Malloc(2*nthreads*(tmax*tspace+2*BORDER_MAX+2),"Opening .1aln file");
  else
    d = NULL;
  if (a == NULL || b == NULL || c == NULL || (see_seq && d == NULL))
    { free(d);
      free(c);
      free(b);
      free(a);
      if (units1 != units2)
        free(units2);
      free(units1);
      goto closeout;
    }

  for (i = 0; i < nthreads; i++)
    { a[i] = b+i;

      b[i].nthreads = nthreads;
      b[i].file     = input+i;
      b[i].read     = false;
      b[i].eof      = false;
      b[i].tspace   = tspace;
      b[i].acnt     = acnt;

      b[i].gpt1 = &(b[i].gdb1);
      b[i].gdb1 = *gdb1;
      if (see_seq)
        b[i].gdb1.seqs = units1[i];
      else
        b[i].gdb1.seqs = NULL;

      if (gdb1 == gdb2)
        b[i].gpt2 = &(b[i].gdb1);
      else
        { b[i].gpt2 = &(b[i].gdb2);
          b[i].gdb2 = *gdb2;
          if (see_seq)
            b[i].gdb2.seqs = units2[i];
          else
            b[i].gdb2.seqs = NULL;
        }

      b[i].align.tpoints = c + 3*i*tmax;
      b[i].align.tdiffs  = c + 3*i*tmax + tmax;

      if (see_seq)
        { b[i].work = New_Work_Data();
          b[i].aseq = d + 2*i*(tmax*tspace+2*BORDER_MAX+2);
          b[i].bseq = b[i].aseq + (tmax*tspace+2*BORDER_MAX+2);
        }
      else
        { b[i].work = NULL;
          b[i].aseq = NULL;
          b[i].bseq = NULL;
        }

      b[i].cigar_max = 0;
      b[i].cigar     = NULL;
    }

  if (units1 != units2)
    free(units2);
  free(units1);

  for (i = 0; i < nthreads; i++)
    { oneGoto(input+i,'A',1);
      if ( ! oneReadLine(input+i))
        b[i].eof = b[i].read = true;
    }

  Error_Buffer = last_error;

  return ((AlnReader *) a);

closeout:
  oneFileClose(input);
  Error_Buffer = last_error;
  EXIT(NULL);
}

  //  Each reader conceptually has a "cursor" pointing at an alignment record that can
  //    be advanced one record at a time with alnNext or jumped to the idx'th record
  //    in the file with alnGoto.  To be crystal clear, alnGoto(..,1) takes you to the
  //    first alignment record, i.e. indexing begins at 1.  alnNext returns true if and
  //    only if the reader has been advnaced to the end of the file.  AlnGoto returns
  //    true if and only if it was successful.

bool alnNext(AlnReader *reader)
{ _AlnReader *r = (_AlnReader *) *reader;

  if (r->read)
    { if (r->eof)
        r->read = true;
    }
  else
    { while (1)
        { if (oneReadLine(r->file))
            { if (r->file->lineType == 'A')
                break;
            }
          else
            { r->eof  = true;
              r->read = true;
              break;
            }
        }
    }
  return (r->eof);
}

bool alnGoto(AlnReader *reader, int idx)
{ _AlnReader *r = (_AlnReader *) *reader;
  bool valid;
#ifdef CHECK_ARGS
  char *last_error;
#else
  (void) valid;
#endif

#ifdef CHECK_ARGS
  last_error = Error_Buffer;
  Error_Buffer = alnMessage;
  if (idx <= 0 || idx > r->acnt)
    { EPRINTF("Index out of range (alnGoto)");
      Error_Buffer = last_error;
      EXIT(false);
    }
#endif
  valid = oneGoto(r->file,'A',idx);
#ifdef CHECK_ARGS
  if (!valid)
    { EPRINTF("Could not seek to %d'th alignment",idx);
      Error_Buffer = last_error;
      EXIT(false);
    }
  Error_Buffer = last_error;
#endif
  oneReadLine(r->file);
  r->read = false;
  r->eof  = false;
  return (true);
}

  //  Return true iff at the end of file

bool alnEOF(AlnReader *reader)
{ _AlnReader *r = (_AlnReader *) *reader;
  return (r->eof);
}

  //  Close all readers associated with a supplied master reader.

void alnCloseReader(AlnReader *reader)
{ _AlnReader *r = (_AlnReader *) *reader;
  _AlnReader *q;
  int i;

  for (i = 0; i < r->nthreads; i++)
    { q = (_AlnReader *) reader[i];
      if (q->work != NULL)
        Free_Work_Data(q->work);
      if (q->gdb1.seqs != NULL)
        fclose(q->gdb1.seqs);
      if (q->gpt1 != q->gpt2)
        if (q->gdb2.seqs != NULL)
          fclose(q->gdb2.seqs);
    }
  Close_GDB(r->gpt1);
  if (r->gpt1 != r->gpt2)
    Close_GDB(r->gpt2);
  oneFileClose(r->file);
  free(r->aseq);
  free(r->align.tpoints);
  free(r);
  free(reader);
}
 

  //  Return (a) the total number of alignments in the file, (b) the maximum number of
  //    trace intervals in any alignment record, (c) the total number of trace intervals
  //    in the file, and (d) the trace point spacing, respectively.

int  alnCount       (AlnReader *reader);
int  alnTraceMax    (AlnReader *reader);
int  alnTraceCount  (AlnReader *reader);
int  alnTraceSpacing(AlnReader *reader);

int alnCount(AlnReader *reader)
{ _AlnReader *r = (_AlnReader *) *reader;

  return (r->acnt);
}

int alnTraceMax(AlnReader *reader)
{ _AlnReader *r = (_AlnReader *) *reader;
  int64 val;

  oneStats(r->file,'T',NULL,&val,NULL);

  return (val);
}

int  alnTraceCount(AlnReader *reader)
{ _AlnReader *r = (_AlnReader *) *reader;
  int64 val;

  oneStats(r->file,'T',NULL,NULL,&val);

  return (val);
}

int alnTraceSpacing(AlnReader *reader)
{ _AlnReader *r = (_AlnReader *) *reader;

  return (r->tspace);
}

  //  GENOME DATA BASE: GDB

  //  From the .1aln reader get the first and second genome data bases over which
  //    the alignments were found.  The two pointers are equal if there was only one genome,
  //    i.e. a self-comparison.  A genome data base gives one the structure of the genome
  //    as a collections of scaffolded contigs separated by gaps, as well as access to
  //    the contig's sequence if 'see_seq' was true in the call that created the reader.

AlnGDB *alnGDB1(AlnReader *reader)
{ _AlnReader *r = (_AlnReader *) *reader;
  return (r->gpt1);
}

AlnGDB *alnGDB2(AlnReader *reader)
{ _AlnReader *r = (_AlnReader *) *reader;
  return (r->gpt2);
}

  // The Count... routines return the # of (a) scaffolds, (b) contigs, and (c) gaps in a AlnGDB.
  //   The Max... routines return the maximum # of contigs/gaps in any scaffold of a AlnGDB.

int gdbScaffoldCount(AlnGDB *g)
{ return (((GDB *) g)->nscaff); }

int gdbContigCount(AlnGDB *g)
{ return (((GDB *) g)->ncontig); }

int gdbGapCount(AlnGDB *g)
{ GDB          *d = (GDB *) g;
  GDB_SCAFFOLD *s = d->scaffolds;
  GDB_CONTIG   *c = d->contigs;
  int           i, x;

  x = d->ncontig - d->nscaff;
  for (i = 0; i < d->nscaff; i++)
    { if (c->sbeg > 0)
        x += 1;
      c += s->ectg-s->fctg; 
      if (c[-1].sbeg + c[-1].clen < s->slen)
        x += 1;
      s += 1;
    }
    
  return (x);
}

int gdbContigMax(AlnGDB *g)
{ GDB          *d = (GDB *) g;
  GDB_SCAFFOLD *s = d->scaffolds;
  int           i, x, cmax;

  cmax = 0;
  for (i = 0; i < d->nscaff; i++)
    { x = s->ectg - s->fctg;
      if (x > cmax)
        cmax = x;
      s += 1;
    }
      
  return (cmax);
}

int gdbGapMax(AlnGDB *g)
{ GDB          *d = (GDB *) g;
  GDB_SCAFFOLD *s = d->scaffolds;
  int           i, x, gmax;

  gmax = 0;
  for (i = 0; i < d->nscaff; i++)
    { x = (s->ectg - s->fctg) - 1;
      if (d->contigs[s->fctg].sbeg > 0)
        x += 1;
      if (d->contigs[s->ectg-1].sbeg + d->contigs[s->ectg-1].clen < s->slen)
        x += 1;
      if (x > gmax)
        gmax = x;
      s += 1;
    }
      
  return (gmax);
}

  // Scaffolds are numbered starting at 1, and the first contig of a scaffold is indexed by 1.
  //   In rare, somewhat abnormal cases there can be a run of N's at the start or end of
  //   of a scaffold, the length of these are obtained by index 0 and ContigsInScaffold(g,s).
  //   If checking is on these routines return -1 if an index is out of bounds.

int gdbScaffoldLen(AlnGDB *g, int s)          // Length of the s'th scaffold
{ GDB *d = (GDB *) g;

#ifdef CHECK_ARGS
  char *last_error;

  last_error = Error_Buffer;
  Error_Buffer = alnMessage;
  if (s < 1 || s > d->nscaff)
    { EPRINTF("Scaffold index out of range (gdbScaffoldLen)");
      Error_Buffer = last_error;
      EXIT(-1);
    }
  Error_Buffer = last_error;
#endif

  return (d->scaffolds[s-1].slen);
}

int gdbContigLen(AlnGDB *g, int s, int c)     // Length of the c'th contig of the s'th scaffold
{ GDB *d = (GDB *) g;
  GDB_SCAFFOLD *scaf = d->scaffolds + (s-1);

#ifdef CHECK_ARGS
  char *last_error;

  last_error = Error_Buffer;
  Error_Buffer = alnMessage;
  if (s < 1 || s > d->nscaff)
    { EPRINTF("Scaffold index out of range (gdbContigLen)");
      Error_Buffer = last_error;
      EXIT(-1);
    }
  if (c < 1 || c > scaf->ectg - scaf->fctg)
    { EPRINTF("Contig index out of range (gdbContigLen)");
      Error_Buffer = last_error;
      EXIT(-1);
    }
  Error_Buffer = last_error;
#endif

  return (d->contigs[scaf->fctg+(c-1)].clen);
}

int gdbGapLen(AlnGDB *g, int s, int p)       // Length of the p'th gap of the s'th scaffold
{ GDB          *d = (GDB *) g;
  GDB_SCAFFOLD *scaf = d->scaffolds + (s-1);
  GDB_CONTIG   *ctg  = d->contigs;

#ifdef CHECK_ARGS
  char *last_error;

  last_error = Error_Buffer;
  Error_Buffer = alnMessage;
  if (s < 1 || s > d->nscaff)
    { EPRINTF("Scaffold index out of range (gdbGapLen)");
      Error_Buffer = last_error;
      EXIT(-1);
    }
  if (p < 0 || p > scaf->ectg - scaf->fctg)
    { EPRINTF("Gap index out of range (gdbGapLen)");
      Error_Buffer = last_error;
      EXIT(-1);
    }
  Error_Buffer = last_error;
#endif

  ctg += scaf->fctg + (p-1);
  if (p == 0)
    return (ctg[1].sbeg);
  if (p >= scaf->ectg-scaf->fctg)
    return (scaf->slen - (ctg->sbeg + ctg->clen));
  return (ctg[1].sbeg - (ctg->sbeg + ctg->clen));
}

int gdbScaffoldContigs(AlnGDB *g, int s)      // Number of contigs in the s'th scaffold
{ GDB *d = (GDB *) g;
  GDB_SCAFFOLD *scaf = d->scaffolds + (s-1);

#ifdef CHECK_ARGS
  char *last_error;

  last_error = Error_Buffer;
  Error_Buffer = alnMessage;
  if (s < 1 || s > d->nscaff)
    { EPRINTF("Scaffold index out of range (gdbScaffoldLen)");
      Error_Buffer = last_error;
      EXIT(-1);
    }
  Error_Buffer = last_error;
#endif

  return (scaf->ectg - scaf->fctg);
}

int gdbContigStart(AlnGDB *g, int s, int c)  // Start position of the c'th contig of the s'th
{ GDB *d = (GDB *) g;                        //   scaffold (see getScaffoldSeq)
  GDB_SCAFFOLD *scaf = d->scaffolds + (s-1);

#ifdef CHECK_ARGS
  char *last_error;

  last_error = Error_Buffer;
  Error_Buffer = alnMessage;
  if (s < 1 || s > d->nscaff)
    { EPRINTF("Scaffold index out of range (gdbContigStart)");
      Error_Buffer = last_error;
      EXIT(-1);
    }
  if (c < 1 || c > scaf->ectg - scaf->fctg)
    { EPRINTF("Contig index out of range (gdbContigStart)");
      Error_Buffer = last_error;
      EXIT(-1);
    }
  Error_Buffer = last_error;
#endif

  return (d->contigs[scaf->fctg+(c-1)].sbeg);
}

  // Get the scaffold name of the s'th scaffold.  The string returned is in a buffer internal
  //   to the AlnGDB and is overwritten each time this routine is called.  You should copy it to
  //   memory you control if you wish it to persist beyond the last call.

char *gdbScaffoldName(AlnGDB *g, int s)
{ GDB *d = (GDB *) g;

#ifdef CHECK_ARGS
  char *last_error;

  last_error = Error_Buffer;
  Error_Buffer = alnMessage;
  if (s < 1 || s > d->nscaff)
    { EPRINTF("Scaffold index out of range (gdbScaffoldName)");
      Error_Buffer = last_error;
      EXIT(NULL);
    }
  Error_Buffer = last_error;
#endif

  return (d->headers + d->scaffolds[s-1].hoff);
}

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

static bool getContigSeq(GDB *d, int c, int beg, int end, uint8 *seq)
{ int cbeg, cend, clen;
  int len, byte;
  int x;

  len  = end-beg;
  cbeg = (beg>>2);
  cend = ((end+3)>>2);
  clen = cend-cbeg;

  x = fseeko(d->seqs,cbeg+d->contigs[c].boff,SEEK_SET);
  if (x < -1)
    { EPRINTF("Could not seek sequence file");
      EXIT(true);
    }
  x = fread(seq,clen,1,d->seqs);
  if (x != 1)
    { EPRINTF("Could not read sequence file");
      EXIT(true);
    }

  byte = seq[--clen];
  switch (end&0x3)
  { case 0:
      seq[--len] = (byte>>6) & 0x3;
    case 3:
      seq[--len] = (byte>>4) & 0x3;
    case 2:
      seq[--len] = (byte>>2) & 0x3;
    case 1:
      seq[--len] = byte & 0x3;
      break;
  }
  for ( ; len >= 4; )
    { byte = seq[--clen];
      seq[--len] = (byte>>6) & 0x3;
      seq[--len] = (byte>>4) & 0x3;
      seq[--len] = (byte>>2) & 0x3;
      seq[--len] = byte & 0x3;
    }
  if (len > 0)
    { byte = seq[--clen];
      for (x = 6; len > 0; x -= 2)
        seq[--len] = (byte>>x) & 0x3;
    }
  return (false);
}

char *gdbScaffoldSeq(AlnGDB *g, int s, int beg, int end, char *buffer)
{ static char letter[4] = { 'a', 'c', 'g', 't' };

  GDB  *d = (GDB *) g;
  int   c, x, len, flip;
  char *seq;
  char *last_error;

  last_error = Error_Buffer;
  Error_Buffer = alnMessage;

#ifdef CHECK_ARGS
  if (d->seqs == NULL)
    { EPRINTF("Sequence access not requested on open (gdbScaffoldSeq)");
      Error_Buffer = last_error;
      EXIT(NULL);
    }
  if (s < 1 || s > d->nscaff)
    { EPRINTF("Scaffold index out of range (gdbScaffoldSeq)");
      Error_Buffer = last_error;
      EXIT(NULL);
    }
#endif

  flip = 0;
  if (beg > end)
    { flip = beg;
      beg  = end;
      end  = flip;
      flip = 1;
    }

  for (c = d->scaffolds[s-1].fctg; c < d->scaffolds[s-1].ectg; c++)
    if (d->contigs[c].sbeg > beg)
      break;
  c -= 1;

  x = d->contigs[c].sbeg;
#ifdef CHECK_ARGS
  if (c < 0 || end > x + d->contigs[c].clen)
    { EPRINTF("Scaffold interval intersects a gap (gdbScaffoldSeq)");
      Error_Buffer = last_error;
      EXIT(NULL);
    }
#endif

  beg -= x;
  end -= x;
  len  = end-beg; 

  if (buffer == NULL)
    { seq = Malloc(len+1,"Allocating sequence (gdbScaffoldSeq)");
      if (seq == NULL)
        { Error_Buffer = last_error;
          EXIT(NULL);
        }
    }
  else
    seq = buffer;

  if (getContigSeq(d,c,beg,end,(uint8 *) seq))
    { if (buffer == NULL)
        free(seq);
      Error_Buffer = last_error;
      EXIT(NULL);
    }

  if (flip)
    { int y, u;

      for (x = 0, y = len-1; x <= y; x++, y--)
        { u = letter[3-seq[x]];
          seq[x] = letter[3-seq[y]];
          seq[y] = u;
        }
    }
  else
    for (x = 0; x < len; x++)
      seq[x] = letter[(int) seq[x]];

  seq[len] = '\0';

  Error_Buffer = last_error;
  return (seq);
}

  //  ALIGNMENT RECORDS

  // Load the alignment record at the reader's current cursor/position.  The return pointer
  //   points at a pre-allocated buffer internal to the reader that is reused/overwritten
  //   each time alnAlignment is called.  This includes the memory for the tpoints and tdiffs
  //   arrays of the trace.  You should copy it and the two arrays just mentioned if you
  //   wish them to persist beyond the last call.

AlnRecord *alnAlignment(AlnReader *reader, bool see_trace)
{ _AlnReader *r  = (_AlnReader *) *reader;
  OneFile    *of = r->file;
  GDB_CONTIG *c1, *c2;
  bool        comp;
  int         tlen, xlen;
  int64      *trace64;
  int         j;
  char       *last_error;

  if (r->eof)
    return (NULL);

  last_error = Error_Buffer;
  Error_Buffer = alnMessage;

  r->contig1     = oneInt(of,0);
  r->align.bpos1 = oneInt(of,1);
  r->align.epos1 = oneInt(of,2);
  r->contig2     = oneInt(of,3);
  r->align.bpos2 = oneInt(of,4);
  r->align.epos2 = oneInt(of,5);
  r->align.diffs = -1;
  r->align.tlen  = 0;
  
  tlen = 0;
  xlen = 0;
  comp = false;
  while (1)
    if ( ! oneReadLine(of))
      { r->eof = true;
        break;
      }
    else if (of->lineType == 'T')
      { if (see_trace)
          { tlen = oneLen(of);
            trace64 = oneIntList(of);
            for (j = 0; j < tlen; j++)
              r->align.tpoints[j] = trace64[j];
          }
      }
    else if (of->lineType == 'X')
      { if (see_trace)
          { xlen = oneLen(of);
            trace64 = oneIntList(of);
            for (j = 0; j < xlen; j++)
              r->align.tdiffs[j] = trace64[j];
          }
      }
    else if (of->lineType == 'R')
      comp = true;
    else if (of->lineType == 'D')
      r->align.diffs = oneInt(of,0);
    else if (of->lineType == 'A')
      break;

  r->read = true;
  if (see_trace)
    { if (tlen == 0)
        { EPRINTF("T-line not present in alignment record");
          Error_Buffer = last_error;
          EXIT(NULL);
        }
      if (xlen != tlen)
        { EPRINTF("T- and X-lines do not have the same length");
          Error_Buffer = last_error;
          EXIT(NULL);
        }
    }
  else
    tlen = 0;
  if (r->align.diffs < 0)
    { EPRINTF("D-line not present in alignment record");
      Error_Buffer = last_error;
      EXIT(NULL);
    }

  r->align.tlen = tlen;

  c1 = r->gpt1->contigs+r->contig1;
  r->align.bpos1 += c1->sbeg;
  r->align.epos1 += c1->sbeg;
  r->align.seq1   = c1->scaf+1;

  c2 = r->gpt2->contigs+r->contig2;
  if (comp)
    { r->align.bpos2 = (c2->clen + c2->sbeg) - r->align.bpos2;
      r->align.epos2 = (c2->clen + c2->sbeg) - r->align.epos2;
    }
  else
    { r->align.bpos2 += c2->sbeg;
      r->align.epos2 += c2->sbeg;
    }
  r->align.seq2 = c2->scaf+1;

  Error_Buffer = last_error;
  return (&r->align);
}

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

static char dna[4] = { 'a', 'c', 'g', 't' };

static void print_bit_seq(uint8 *seq, int len)
{ int i;

  for (i = 0; i < len; i++)
    { if (i%100 == 0 && i > 0)
        putchar('\n');
      putchar(dna[seq[i]]);
    }
  putchar('\n');
}

static Cigar *cigar_core(AlnRecord *align, Alignment *aln, bool show_x, int *clen)
{ _AlnReader *reader;
  Path       *path;
  int         amin, amax;
  int         bmin, bmax;
  GDB        *g1, *g2;
  char       *aseq, *bseq;
  int         cnum, cmax;
  Cigar      *cigar;

  (void) print_bit_seq;

  reader = (_AlnReader *) align;

  g1 = reader->gpt1;
  g2 = reader->gpt2;

  path = aln->path;

  amin = align->bpos1 - g1->contigs[reader->contig1].sbeg;
  amax = align->epos1 - g1->contigs[reader->contig1].sbeg;
  bmin = align->bpos2 - g2->contigs[reader->contig2].sbeg;
  bmax = align->epos2 - g2->contigs[reader->contig2].sbeg;

  aseq = reader->aseq+1;
  bseq = reader->bseq+1;
  path->trace = (uint16 *) (align->tdiffs + (align->tdiffs-align->tpoints));

  { int j, x;
    uint16 *trace = (uint16 *) path->trace;

    path->tlen = 2*align->tlen;
    j = 0;
    for (x = 1; x < path->tlen; x += 2)
      trace[x] = align->tpoints[j++];
    j = 0;
    for (x = 0; x < path->tlen; x += 2)
      trace[x] = align->tdiffs[j++];
  }

  if (getContigSeq(g1, reader->contig1, amin, amax, (uint8 *) aseq))
    EXIT(NULL);
  amax -= amin;
  aseq[-1] = 4;
  aseq[amax] = 4;

  amin %= reader->tspace;
  amax += amin;
  path->abpos = amin;
  path->aepos = amax;
  aln->aseq   = aseq-amin;
  aln->alen   = amax;

  if (bmin > bmax)
    { if (getContigSeq(g2, reader->contig2, bmax, bmin, (uint8 *) bseq))
        EXIT(NULL);
      Complement_Seq(bseq,bmin-bmax);
      bmax = bmin-bmax;
      aln->flags = COMP_FLAG;
    }
  else
    { if (getContigSeq(g2, reader->contig2, bmin, bmax, (uint8 *) bseq))
        EXIT(NULL);
      bmax -= bmin;
      aln->flags = 0;
    }
  bseq[-1] = 4;
  bseq[bmax] = 4;

  path->bbpos = 0;
  path->bepos = bmax;
  aln->bseq   = bseq;
  aln->blen   = bmax;

  if (Compute_Trace_PTS(aln,reader->work,reader->tspace,GREEDIEST,1,-1))
    EXIT(NULL);

  if (Gap_Improver(aln,reader->work))
    EXIT(NULL);

  cmax  = reader->cigar_max;
  cigar = reader->cigar;
  if (cigar == NULL)
    { cmax  = 3*path->tlen + 500;
      reader->cigar = cigar = Malloc(sizeof(Cigar)*cmax,"Allocating CIGAR array");
      if (cigar == NULL)
        EXIT(NULL);
    }
  cnum  = 0;

#define ADD2CIGAR(LEN,OP)								\
{ if (cnum >= cmax)									\
    { cmax = cmax * 1.2 + 100;								\
      cigar = (Cigar *) Realloc(cigar,sizeof(Cigar)*cmax,"Expanding CIGAR array");	\
      reader->cigar = cigar;								\
      if (cigar == NULL)								\
        EXIT(NULL);									\
    }											\
  cigar[cnum].op =  OP;									\
  cigar[cnum].len = LEN;								\
  cnum += 1;										\
}

  { int    k, h, p, x, b, blen;
    int32 *t = (int32 *) path->trace;
    int    T = path->tlen;
    int    ilen, dlen;
    int    xlen, elen;
    char  *A, *B;

    A = aln->aseq-1;
    B = aln->bseq-1;
    ilen = dlen = 0;
    k = path->abpos+1;
    h = 1;
    for (x = 0; x < T; x++)
      { if ((p = t[x]) < 0)
          { blen = -(p+k);
            if (dlen > 0)
              ADD2CIGAR(dlen,'D');
            dlen = 0;
            if (blen == 0)
              ilen += 1;
            else
              { if (ilen > 0)
                  ADD2CIGAR(ilen,'I');
                if (show_x)
                  { elen = xlen = 0;
                    for (b = 0; b < blen; b++, k++, h++)
                      if (A[k] == B[h])
                        { if (xlen > 0)
                            ADD2CIGAR(xlen,'X');
                          xlen = 0;
                          elen += 1;
                        }
                      else
                        { if (elen > 0)
                            ADD2CIGAR(elen,'=');
                          elen = 0;
                          xlen += 1;
                        }
                    if (xlen > 0)
                      ADD2CIGAR(xlen,'X');
                    if (elen > 0)
                      ADD2CIGAR(elen,'=');
                  }
                else
                  ADD2CIGAR(blen,'M');
                ilen = 1;
              }
            h += 1;
          }
        else
          { blen = p-h;
            if (ilen > 0)
              ADD2CIGAR(ilen,'I');
            ilen = 0;
            if (blen == 0)
              dlen += 1;
            else
              { if (dlen > 0)
                  ADD2CIGAR(dlen,'D');
                if (show_x)
                  { elen = xlen = 0;
                    for (b = 0; b < blen; b++, k++, h++)
                      if (A[k] == B[h])
                        { if (xlen > 0)
                            ADD2CIGAR(xlen,'X');
                          xlen = 0;
                          elen += 1;
                        }
                      else
                        { if (elen > 0)
                            ADD2CIGAR(elen,'=');
                          elen = 0;
                          xlen += 1;
                        }
                    if (xlen > 0)
                      ADD2CIGAR(xlen,'X');
                    if (elen > 0)
                      ADD2CIGAR(elen,'=');
                  }
                else
                  ADD2CIGAR(blen,'M');
                dlen = 1;
              }
            k += 1;
          }
      }
    if (dlen > 0)
      ADD2CIGAR(dlen,'D');
    if (ilen > 0)
      ADD2CIGAR(ilen,'I');
    blen = (path->aepos - k)+1;
    if (blen > 0)
      { if (show_x)
          { elen = xlen = 0;
            for (b = 0; b < blen; b++, k++, h++)
              if (A[k] == B[h])
                { if (xlen > 0)
                    ADD2CIGAR(xlen,'X');
                  xlen = 0;
                  elen += 1;
                }
              else
                { if (elen > 0)
                    ADD2CIGAR(elen,'=');
                  elen = 0;
                  xlen += 1;
                }
            if (xlen > 0)
              ADD2CIGAR(xlen,'X');
            if (elen > 0)
              ADD2CIGAR(elen,'=');
          }
        else
          ADD2CIGAR(blen,'M');
      }
  }

  reader->cigar_max = cmax;

  *clen = cnum;
  return (cigar);
}

char *alnCreateCigar(AlnRecord *align, bool show_x, bool reverse)
{ char        buf[12];
  char       *s;
  int         n;
  Alignment  _aln, *aln = &_aln;
  Path       _path, *path = &_path;
  Cigar      *cigar;
  int         clen, slen;
  char       *cstring, *cs;
  int         beg, end, step;
  int         i;
  char       *last_error;

  last_error = Error_Buffer;
  Error_Buffer = alnMessage;

#ifdef CHECK_ARGS
  if (align->tlen == 0)
    { EPRINTF("Trace information was not requested when record fetched");
      Error_Buffer = last_error;
      EXIT(NULL);
    }
#endif

  aln->path = path;

  cigar = cigar_core(align,aln,show_x,&clen);
  if (cigar == NULL)
    { Error_Buffer = last_error;
      EXIT(NULL);
    }

  slen = 2*clen+1;
  for (i = 0; i < clen; i++)
    { n = cigar[i].len;
      while (n >= 10)
        { slen += 1;
           n /= 10;
        }
    }

  cstring = cs = Malloc(slen,"Allocating CIGAR string");
  if (cstring == NULL)
    { Error_Buffer = last_error;
      free(cigar);
      EXIT(NULL);
    }
  Error_Buffer = last_error;

  beg  = 0;
  end  = clen;
  step = 1;
  if (reverse)
    { if (COMP(aln->flags))
        { beg  = clen-1;
          end  = -1;
          step = -1;
        }
      for (i = 0; i < clen; i++)
        if (cigar[i].op == 'I')
          cigar[i].op = 'D';
        else if (cigar[i].op == 'D')
          cigar[i].op = 'I';
    }

  for (i = beg; i != end; i += step)
    { s = buf;
      n = cigar[i].len;
      while (n >= 10)
        { *s++ = '0' + (n%10);
           n /= 10;
        }
      *cs++ = '0'+n;
      while (s > buf)
        *cs++ = *--s;
      *cs++ = cigar[i].op;
    }
  *cs = '\0';

  return (cstring);
}

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

char *alnCreateCStag(AlnRecord *align, bool short_form, bool reverse)
{ char        buf[12];
  char       *s;
  int         n;
  Alignment  _aln, *aln = &_aln;
  Path       _path, *path = &_path;
  Cigar      *cigar;
  int         clen, slen;
  char       *cstring, *cs;
  int         i, j, w;
  int         beg, end, step;
  uint8      *A, *B;
  char       *last_error;

  last_error = Error_Buffer;
  Error_Buffer = alnMessage;

#ifdef CHECK_ARGS
  if (align->tlen == 0)
    { EPRINTF("Trace information was not requested when record fetched");
      Error_Buffer = last_error;
      EXIT(NULL);
    }
#endif

  aln->path = path;

  cigar = cigar_core(align,aln,1,&clen);
  if (cigar == NULL)
    { Error_Buffer = last_error;
      EXIT(NULL);
    }

  if (short_form)
    { for (i = 0; i < clen; i++)
        if (cigar[i].op == '=')
          cigar[i].op = 'M';
    }

  slen = clen+1;
  for (i = 0; i < clen; i++)
    { if (cigar[i].op == 'M')
        { n = cigar[i].len;
          while (n >= 10)
            { slen += 1;
               n /= 10;
            }
          slen += 1;
        }
      else if (cigar[i].op == 'X')
        slen += 2*cigar[i].len;
      else
        slen += cigar[i].len;
    }

  cstring = cs = Malloc(slen,"Allocating CS-tag string");
  if (cstring == NULL)
    { Error_Buffer = last_error;
      EXIT(NULL);
    }
  Error_Buffer = last_error;

  A    = (uint8 *) aln->aseq+aln->path->abpos;
  B    = (uint8 *) aln->bseq;
  beg  = 0;
  end  = clen;
  step = 1;
  if (reverse)
    { if (COMP(aln->flags))
        { Complement_Seq((char *) A,path->aepos-path->abpos);
          Complement_Seq((char *) B,path->bepos);
          beg  = clen-1;
          end  = -1;
          step = -1;
        }
      B = A;
      A = (uint8 *) aln->bseq;
      for (i = 0; i < clen; i++)
        if (cigar[i].op == 'I')
          cigar[i].op = 'D';
        else if (cigar[i].op == 'D')
          cigar[i].op = 'I';
    }

  cs = cstring;
  for (i = beg; i != end; i += step)
    { w = cigar[i].len;
      switch (cigar[i].op)
      { case '=':
          *cs++ = '=';
          for (j = 0; j < w; j++)
            *cs++ = dna[*A++];
          B += w;
          break;
        case 'M':
          *cs++ = ':';
          s = buf;
          n = cigar[i].len;
          while (n >= 10)
            { *s++ = '0' + (n%10);
               n /= 10;
            }
          *cs++ = '0'+n;
          while (s > buf)
            *cs++ = *--s;
          A += w;
          B += w;
          break;
        case 'X':
          *cs++ = '*';
          for (j = 0; j < w; j++)
            { *cs++ = dna[*A++];
              *cs++ = dna[*B++];
            }
          break;
        case 'I':
          *cs++ = '+';
          for (j = 0; j < w; j++)
            *cs++ = dna[*B++];
          break;
        case 'D':
          *cs++ = '-';
          for (j = 0; j < w; j++)
            *cs++ = dna[*A++];
          break;
        default:
          break;
      }
    }
  *cs++ = '\0';

  return (cstring);
}

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
  //   is false, then the indel array is one that transforms the second sequence into the
  //   first where the second sequence is in the forward orientation, i.e. the roles of the
  //   first and sequence sequence are reversed.

int  *alnCreateIndelArray(AlnRecord *align, bool reversed)
{ _AlnReader *reader;
  Alignment  _aln, *aln = &_aln;
  Path       _path, *path = &_path;
  int         amin, amax;
  int         bmin, bmax;
  GDB        *g1, *g2;
  char       *aseq, *bseq;
  int        *trace;
  char       *last_error;

  reader = (_AlnReader *) align;

  last_error = Error_Buffer;
  Error_Buffer = alnMessage;

#ifdef CHECK_ARGS
  if (align->tlen == 0)
    { EPRINTF("Trace information was not requested when record fetched");
      Error_Buffer = last_error;
      EXIT(NULL);
    }
#endif

  g1 = reader->gpt1;
  g2 = reader->gpt2;

  aln->path = path;

  amin = align->bpos1 - g1->contigs[reader->contig1].sbeg;
  amax = align->epos1 - g1->contigs[reader->contig1].sbeg;
  bmin = align->bpos2 - g2->contigs[reader->contig2].sbeg;
  bmax = align->epos2 - g2->contigs[reader->contig2].sbeg;

  aseq = reader->aseq+1;
  bseq = reader->bseq+1;
  path->trace = (uint16 *) (align->tdiffs + (align->tdiffs-align->tpoints));

  { int j, x;
    uint16 *trace = (uint16 *) path->trace;

    path->tlen = 2*align->tlen;
    j = 0;
    for (x = 1; x < path->tlen; x += 2)
      trace[x] = align->tpoints[j++];
    j = 0;
    for (x = 0; x < path->tlen; x += 2)
      trace[x] = align->tdiffs[j++];
  }

  if (getContigSeq(g1, reader->contig1, amin, amax, (uint8 *) aseq))
    { Error_Buffer = last_error;
      EXIT(NULL);
    }
  amax -= amin;
  aseq[-1] = 4;
  aseq[amax] = 4;

  amin %= reader->tspace;   
  amax += amin;
  path->abpos = amin;
  path->aepos = amax;
  aln->aseq   = aseq-amin;
  aln->alen   = amax;

  if (bmin > bmax)
    { if (getContigSeq(g2, reader->contig2, bmax, bmin, (uint8 *) bseq))
        { Error_Buffer = last_error;
          EXIT(NULL);
        }
      Complement_Seq(bseq,bmin-bmax);
      bmax = bmin-bmax;
      aln->flags = COMP_FLAG;
    }
  else
    { if (getContigSeq(g2, reader->contig2, bmin, bmax, (uint8 *) bseq))
        { Error_Buffer = last_error;
          EXIT(NULL);
        }
      bmax -= bmin;
      aln->flags = 0;
    }
  bseq[-1] = 4;
  bseq[bmax] = 4;

  path->bbpos = 0;
  path->bepos = bmax;
  aln->bseq   = bseq;
  aln->blen   = bmax;

  if (Compute_Trace_PTS(aln,reader->work,reader->tspace,GREEDIEST,1,-1))
    { Error_Buffer = last_error;
      EXIT(NULL);
    }

  if (Gap_Improver(aln,reader->work))
    { Error_Buffer = last_error;
      EXIT(NULL);
    }

  trace = (int *) path->trace;

  { int i;

    for (i = 0; i < path->tlen; i++)
      if (trace[i] < 0)
        trace[i] += amin;
    trace[path->tlen] = 0;
  }

  if (reversed)
    { int i, j, x;

      if (COMP(aln->flags))
        { bmax = bmax+2;
          amax = (amax-amin)+2;
          for (i = 0; i < path->tlen; i++)
            if (trace[i] < 0)
              trace[i] = amax+trace[i];
            else
              trace[i] = trace[i]-bmax;
          for (i = 0, j = path->tlen-1; i < j; i++, j--)
            { x = trace[i];
              trace[i] = trace[j];
              trace[j] = x;
            }
        }
      else
        for (i = 0; i < path->tlen; i++)
          trace[i] = -trace[i];
    }

  Error_Buffer = last_error;
  return (trace);
}

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
                       bool upper_case, bool reversed)
{ _AlnReader *reader;
  Alignment  _aln, *aln = &_aln;
  Path       _path, *path = &_path;
  int         amin, amax;
  int         bmin, bmax;
  int         aoff, boff;
  GDB        *g1, *g2;
  char       *aseq, *bseq;
  char       *last_error;
  int        *trace, tlen;
  int         abc, aec, bbc, bec;

  reader = (_AlnReader *) align;

  last_error = Error_Buffer;
  Error_Buffer = alnMessage;

#ifdef CHECK_ARGS
  if (align->tlen == 0)
    { EPRINTF("Trace information was not requested when record fetched");
      Error_Buffer = last_error;
      EXIT(NULL);
    }
#endif

  g1 = reader->gpt1;
  g2 = reader->gpt2;

  aln->path = path;

  aoff = g1->contigs[reader->contig1].sbeg;
  boff = g2->contigs[reader->contig2].sbeg;

  amin = path->abpos = align->bpos1 - aoff;
  amax = path->aepos = align->epos1 - aoff;
  bmin = path->bbpos = align->bpos2 - boff;
  bmax = path->bepos = align->epos2 - boff;

  aln->alen = g1->contigs[reader->contig1].clen;
  aln->blen = g2->contigs[reader->contig2].clen;

  aseq = reader->aseq+1;
  bseq = reader->bseq+1;
  path->trace = (uint16 *) (align->tdiffs + (align->tdiffs-align->tpoints));

  { int     j, x;
    uint16 *trace = (uint16 *) path->trace;

    path->tlen = 2*align->tlen;
    j = 0;
    for (x = 1; x < path->tlen; x += 2)
      trace[x] = align->tpoints[j++];
    j = 0;
    for (x = 0; x < path->tlen; x += 2)
      trace[x] = align->tdiffs[j++];
  }

  amin -= border;
  amax += border;
  if (amin < 0)
    amin = 0;
  if (amax > aln->alen)
    amax = aln->alen;

  if (getContigSeq(g1, reader->contig1, amin, amax, (uint8 *) aseq))
    { Error_Buffer = last_error;
      EXIT(true);
    }
  aln->aseq = aseq - amin;

  if (bmin > bmax)
    { int x;

      bmax -= border;
      bmin += border;
      if (bmax < 0)
        bmax = 0;
      if (bmin > aln->blen)
        bmin = aln->blen;
      if (getContigSeq(g2, reader->contig2, bmax, bmin, (uint8 *) bseq))
        EXIT(NULL);
      Complement_Seq(bseq,bmin-bmax);

      x = path->bbpos;
      path->bbpos = path->bepos;
      path->bepos = x;

      x = bmin-bmax;
      bmin = path->bbpos - (bmin-path->bepos);
      bmax = bmin + x;

      aln->flags  = COMP_FLAG;
    }
  else
    { bmin -= border;
      bmax += border;
      if (bmin < 0)
        bmin = 0;
      if (bmax > aln->blen)
        bmax = aln->blen;
      if (getContigSeq(g2, reader->contig2, bmin, bmax, (uint8 *) bseq))
        EXIT(NULL);
      aln->flags = 0;
    }
  aln->bseq = bseq - bmin;

  abc = aln->aseq[path->abpos-1];
  aec = aln->aseq[path->aepos];
  aln->aseq[path->abpos-1] = 4;
  aln->aseq[path->aepos] = 4;

  bbc = aln->bseq[path->bbpos-1];
  bec = aln->bseq[path->bepos];
  aln->bseq[path->bbpos-1] = 4;
  aln->bseq[path->bepos] = 4;

  if (Compute_Trace_PTS(aln,reader->work,reader->tspace,GREEDIEST,1,-1))
    { Error_Buffer = last_error;
      EXIT(true);
    }

  if (Gap_Improver(aln,reader->work))
    { Error_Buffer = last_error;
      EXIT(true);
    }

  aln->aseq[path->abpos-1] = abc;
  aln->aseq[path->aepos] = aec;
  aln->bseq[path->bbpos-1] = bbc;
  aln->bseq[path->bepos] = bec;

  aseq[-1] = 4;
  aseq[amax-amin] = 4;
  bseq[-1] = 4;
  bseq[bmax-bmin] = 4;

  trace = (int *) path->trace;
  tlen  = path->tlen;

  if (reversed)
    { int   i, j, x;
      char *s;

      if (COMP(aln->flags))
        { Complement_Seq(aseq,amax-amin);
          Complement_Seq(bseq,bmax-bmin);

          printf("shift = %d/%d\n",(path->abpos - amin) - (amax-path->aepos),
                                   (path->bbpos - bmin) - (bmax-path->bepos));

          aln->aseq -= (path->abpos - amin) - (amax-path->aepos);
          aln->bseq -= (path->bbpos - bmin) - (bmax-path->bepos);

          bmax = path->bbpos + path->bepos + 2;
          amax = path->abpos + path->aepos + 2;
          for (i = 0; i < path->tlen; i++)
            if (trace[i] < 0)
              trace[i] = amax+trace[i];
            else
              trace[i] = trace[i]-bmax;

          for (i = 0, j = path->tlen-1; i < j; i++, j--)
            { x = trace[i];
              trace[i] = trace[j];
              trace[j] = x;
            }

        }
      else
        for (i = 0; i < path->tlen; i++)
          trace[i] = -trace[i];

      s = aln->aseq;
      aln->aseq = aln->bseq;
      aln->bseq = s;

      x = aln->alen;
      aln->alen = aln->blen;
      aln->blen = x;

      x = path->abpos;
      path->abpos = path->bbpos;
      path->bbpos = x;

      x = path->aepos;
      path->aepos = path->bepos;
      path->bepos = x;

      x = aoff;
      aoff = boff;
      boff = x;
    }

  { int  i;

    if (COMP(aln->flags))
      aln->blen = 2*boff + (path->bbpos+path->bepos);
    path->abpos += aoff;
    path->aepos += aoff;
    path->bbpos += boff;
    path->bepos += boff;

    aln->aseq -= aoff;
    aln->bseq -= boff;
    for (i = 0; i < tlen; i++)
      if (trace[i] < 0)
        trace[i] -= aoff;
      else
        trace[i] += boff;
  }

  if (Print_Alignment(where,aln,reader->work,indent,width,border,upper_case,coord,0))
    { Error_Buffer = last_error;
      EXIT(true);
    }
  
  return (false);
}


#ifdef TEST

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

  ARG_INIT("ONEalnTEST");

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
    { printf("Scafs = %d\n",gdbScaffoldCount(gdb2));
      printf("Conts = %d\n",gdbContigCount(gdb2));
      printf("Gaps  = %d\n",gdbGapCount(gdb2));
      printf("Cmax  = %d\n",gdbContigMax(gdb2));
      printf("Gmax  = %d\n",gdbGapMax(gdb2));
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

#endif
