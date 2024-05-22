#include <stdio.h>
#include <stdlib.h>
#include "DB.h"

#define IS_FA    0  //  Return types
#define IS_FA_GZ 1
#define IS_ONE   2
#define IS_GDB   3

//  Interpret source & target arguments returning complete path + extension names of
//    the source and target in spath & tpath.  Returns the type of source.
//  If target == source then tpath is a name for a temporary file
//  If target == NULL then tpath has the same path & root as spath
//  Else If target is a path then tpath has the same root as spath
//  If no_gdb then the source cannot be a gdb, otherwise if the source is a gdb
//    then tpath = NULL on return.

int get_dna_paths(char *source, char *target, char **spath, char **tpath, int no_gdb);

FILE **make_temporary_gdb(char *spath, char *tpath, DAZZ_DB *db, int *hdrs, int NTHREADS);
