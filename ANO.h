/*******************************************************************************************
 *
 *  .1bed module.  Auxiliary routines to open and manipulate a ONEcode bed-style file.
 *
 *  Author :  Gene Myers
 *  Date   :  December 2025
 *
 ********************************************************************************************/
  
#ifndef _ANO_DEFS

#define _ANO_DEFS

#include <stdio.h>
#include "gene_core.h"
#include "ONElib.h"
#include "GDB.h"

/*******************************************************************************************
 *
 *  ANO IN-CORE DATA STRUCTURES
 *
 ********************************************************************************************/

typedef struct
  { int64  beg;     //  [beg,end] interval, in contig coords (externally coords are for scaffolds)
    int64  end;
    int    score;   //  score of anotation
    int    orient;  //  orientation of interval: + = 0, - = 1
    char  *label;   //  NULL if none
  } ANO_PAIR;

typedef struct _ANO
  { int            nprov;     //  provenance of file
    OneProvenance *prov;

    int            shared;     //  the gdb is not owned by this bed
    GDB           *gdb;        //  the GDB to which this ANO applies

    int           nints;      //  # of bed intervals
    ANO_PAIR     *masks;      //  array [0..nints) of mask records

    int64        *moff;       //  [moff[i],moff[i+1]) are the ANO intervals for contig i

    int           maxlab;     //  length of longest label
    char         *labels;     //  pointer to memory block of all labels
    char         *headers;    //  pointer to memory block of all headers
  } ANO;


/*******************************************************************************************
 *  
 *  ANO ROUTINES
 *  
 ********************************************************************************************/

OneSchema *make_ANO_Schema();    //  Make a .1bed schema

  //  Read the ANO file spath and build a ANO data structure around the record bed.
  //    If gdb is not NULL then the skeleton in the bed file must be equal (up to
  //    scaffold permuation) to the skeleton of gdb and the order of scaffolds will
  //    be that of gdb.  The ANO_PAIRs for each contig will be sorted on beg. 

int Read_ANO(ANO *bed, char *spath, GDB *gdb);

  //  Create and write a .1bed file tpath for the ANO data structure bed.  Each M-line
  //    will contain at most width pairs

int Write_ANO(ANO *bed, char *tpath, int width);

  //  Create an unlabelled bed that is the union of the beds sbeds[i] for i = 0 to nway-1.
  //    The input beds are checked for equality where they must all have the same
  //    permutation of scaffolds (this can be achieved by supplying a "master" gdb to
  //    each Read_ANO call).

int ANO_Union(ANO *tbed, int nway, ANO **sbeds);

  //  Free the storage for bed.

void Free_ANO(ANO *bed);

  //   For debug purposes, list the contents of a bed.

void Show_ANO(ANO *bed);

#endif   //  _ANO_DEFS
