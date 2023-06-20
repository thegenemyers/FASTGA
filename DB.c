/*******************************************************************************************
 *
 *  Compressed data base module.  Auxiliary routines to open and manipulate a data base for
 *    which the sequence and read information are separated into two separate files, and the
 *    sequence is compressed into 2-bits for each base.  Support for tracks of additional
 *    information, and trimming according to the current partition.
 *
 *  Author :  Gene Myers
 *  Date   :  July 2013
 *  Revised:  April 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <dirent.h>
#include <limits.h>
#include <sys/stat.h>

#include "DB.h"
#include "QV.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif


/*******************************************************************************************
 *
 *  GENERAL UTILITIES
 *
 ********************************************************************************************/

char *Prog_Name;

#ifdef INTERACTIVE

char Ebuffer[1000];

#endif

int Count_Args(char *var)
{ int   cnt, lev;
  char *s;

  cnt = 1;
  lev = 0;
  for (s = var; *s != '\0'; s++)
    if (*s == ',')
      { if (lev == 0)
          cnt += 1;
      }
    else if (*s == '(')
      lev += 1;
    else if (*s == ')')
      lev -= 1;
  return (cnt);
}

static char *MyCatenate(char *path, char *sep, char *root, char *suffix)
{ static char *cat = NULL;
  static int   max = -1;
  int   len;

  if (path == NULL || root == NULL || sep == NULL || suffix == NULL)
    return (NULL);
  len =  strlen(path);
  len += strlen(sep);
  len += strlen(root);
  len += strlen(suffix);
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      cat = (char *) realloc(cat,max+1);
      if (cat == NULL)
        { EPRINTF(EPLACE,"%s: Out of memory (Making path name for %s)\n",Prog_Name,root);
          return (NULL);
        }
    }
  sprintf(cat,"%s%s%s%s",path,sep,root,suffix);
  return (cat);
}

static char *MyNumbered_Suffix(char *left, int num, char *right)
{ static char *sfx = NULL;
  static int   max = -1;
  int   len;

  if (left == NULL || right == NULL)
    return (NULL);
  len =  strlen(left);
  len += strlen(right) + 40;
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      sfx = (char *) realloc(sfx,max+1);
      if (sfx == NULL)
        { EPRINTF(EPLACE,"%s: Out of memory (Making number suffix for %d)\n",Prog_Name,num);
          return (NULL);
        }
    }
  sprintf(sfx,"%s%d%s",left,num,right);
  return (sfx);
}


/*******************************************************************************************
 *
 *  DB STUB HANDLING ROUTINES
 *
 ********************************************************************************************/

  // Read the contents of the DB stub file at "path" and return it encoded in a DAZZ_STUB
  //   structure.  This is allocated by the routine.  "path" is assumed to be the complete
  //   name of the file.

DAZZ_STUB *Read_DB_Stub(char *path, int what)
{ FILE      *dbfile;
  DAZZ_STUB *stub;

  char  buf1[MAX_NAME+100];
  char  buf2[MAX_NAME+100];
  int   nread;

  int   i;
  int   nfiles;
  int   nblocks;
  int64 size;
  int   all, cutoff;

  dbfile = Fopen(path,"r");
  if (dbfile == NULL)
    EXIT(NULL);

  stub = Malloc(sizeof(DAZZ_STUB),"Allocating DB stub record");
  if (stub == NULL)
    { fclose(dbfile);
      EXIT(NULL);
    }

  stub->nreads  = NULL;
  stub->fname   = NULL;
  stub->prolog  = NULL;
  stub->ublocks = NULL;
  stub->tblocks = NULL;

  if (fscanf(dbfile,DB_NFILE,&nfiles) != 1)
    goto stub_trash;

  if (what & DB_STUB_NREADS)
    { stub->nreads = (int *) Malloc(sizeof(int)*(nfiles+1),"Allocating DB stub record");
      if (stub->nreads == NULL)
        goto stub_error;
      stub->nreads += 1;
    }

  if (what & DB_STUB_FILES)
    { stub->fname = (char **) Malloc(sizeof(char *)*(nfiles+1),"Allocating DB stub record");
      if (stub->fname == NULL)
        goto stub_error;
      stub->fname += 1;

      stub->nfiles  = nfiles;
      for (i = 0; i < nfiles; i++)
        stub->fname[i] = NULL;
    }

  if (what & DB_STUB_PROLOGS)
    { stub->prolog = (char **) Malloc(sizeof(char *)*(nfiles+1),"Allocating DB stub record");
      if (stub->prolog == NULL)
        goto stub_error;
      stub->prolog += 1;

      for (i = 0; i < nfiles; i++)
        stub->prolog[i] = NULL;
    }

  for (i = 0; i < nfiles; i++)
    { if (fscanf(dbfile,DB_FDATA,&nread,buf1,buf2) != 3)
        goto stub_trash;
      if (what & DB_STUB_NREADS)
        stub->nreads[i] = nread;
      if (what & DB_STUB_FILES)
        { stub->fname[i]   = Strdup(buf1,"Alloacting DB stub record");
          if (stub->fname[i] == NULL)
            goto stub_error;
        }
      if (what & DB_STUB_PROLOGS)
        { stub->prolog[i] = Strdup(buf2,"Alloacting DB stub record");
          if (stub->prolog[i] == NULL)
            goto stub_error;
        }
    }

  if (fscanf(dbfile,DB_NBLOCK,&nblocks) != 1)
    goto stub_trash;

  if (fscanf(dbfile,DB_PARAMS,&size,&cutoff,&all) != 3)
    goto stub_trash;

  if (what & DB_STUB_BLOCKS)
    { stub->ublocks  = (int *) Malloc(sizeof(int)*(nblocks+1),"Allocating DB stub record");
      stub->tblocks  = (int *) Malloc(sizeof(int)*(nblocks+1),"Allocating DB stub record");
      if (stub->ublocks == NULL || stub->tblocks == NULL)
        goto stub_error;

      for (i = 0; i <= nblocks; i++)
        if (fscanf(dbfile,DB_BDATA,stub->ublocks+i,stub->tblocks+i) != 2)
          goto stub_trash;
    }

  fclose(dbfile);

  stub->nfiles  = nfiles;
  stub->all     = all;
  stub->cutoff  = cutoff;
  stub->bsize   = size;
  stub->nblocks = nblocks;
  return (stub);

stub_trash:
  EPRINTF(EPLACE,"%s: Stub file %s is junk\n",Prog_Name,path);
stub_error:
  Free_DB_Stub(stub);
  fclose(dbfile);
  EXIT(NULL);
}

  // Read the DB stub file "path" and extract the read index range [*first,*last)
  //   for block n, for the trimmed DB if trim is set, the untrimmed DB otherwise.

int Fetch_Block_Range(char *path, int trim, int n, int *first, int *last)
{ FILE *dbfile;
  char  buffer[2*MAX_NAME+100];
  int   nfiles;
  int   nblocks;
  int64 size;
  int   all, cutoff;
  int   tfirst, tlast;
  int   ufirst, ulast;
  int   i;

  dbfile = Fopen(path,"r");
  if (dbfile == NULL)
    EXIT(1);
  if (fscanf(dbfile,DB_NFILE,&nfiles) != 1)
    goto stub_error;
  for (i = 0; i < nfiles; i++)
    if (fgets(buffer,2*MAX_NAME+100,dbfile) == NULL)
      goto stub_error;
  if (fscanf(dbfile,DB_NBLOCK,&nblocks) != 1)
    goto stub_error;

  if (n < 0 || n >= nblocks)
    { *first = *last = -1;
      return (0);
    }

  if (fscanf(dbfile,DB_PARAMS,&size,&cutoff,&all) != 3)
    goto stub_error;
  for (i = 1; i <= n; i++)
    if (fscanf(dbfile,DB_BDATA,&ufirst,&tfirst) != 2)
      goto stub_error;
  if (fscanf(dbfile,DB_BDATA,&ulast,&tlast) != 2)
    goto stub_error;

  fclose(dbfile);

  if (trim)
    { *first = tfirst;
      *last  = tlast;
    }
  else
    { *first = ufirst;
      *last  = ulast;
    }

  return (0);

stub_error:
  fclose(dbfile);
  EPRINTF(EPLACE,"%s: Stub file %s is junk\n",Prog_Name,path);
  EXIT(1);
}

  // Free a DAZZ_STUB data structure returned by Read_DB_Stub

void Free_DB_Stub(DAZZ_STUB *stub)
{ int i;

  if (stub == NULL)
    return;
  if (stub->fname != NULL)
    { for (i = 0; i < stub->nfiles; i++)
        free(stub->fname[i]);
      free(stub->fname-1);
    }
  if (stub->prolog != NULL)
    { for (i = 0; i < stub->nfiles; i++)
        free(stub->prolog[i]);
      free(stub->prolog-1);
    }
  if (stub->nreads != NULL)
    free(stub->nreads-1);
  free(stub->ublocks);
  free(stub->tblocks);
  free(stub);
}


/*******************************************************************************************
 *
 *  DB OPEN, TRIM, SIZE_OF, LIST_FILES & CLOSE ROUTINES
 *
 ********************************************************************************************/


// Open the given database or dam, "path" into the supplied DAZZ_DB record "db". If the name has
//   a part # in it then just the part is opened.  The index array is allocated (for all or
//   just the part) and read in.
// Return status of routine:
//    -1: The DB could not be opened for a reason reported by the routine to EPLACE
//     0: Open of DB proceeded without mishap
//     1: Open of DAM proceeded without mishap

static char *atrack_name = ".@arw";
static char *qtrack_name = ".@qvs";

int Open_DB(char* path, DAZZ_DB *db)
{ DAZZ_DB dbcopy;
  char   *root, *pwd, *bptr, *fptr, *cat;
  int     nreads;
  FILE   *index, *dbvis, *bases;
  int     status, plen, isdam;
  int     part, cutoff, all;
  int     ufirst, tfirst, ulast, tlast;

  status = -1;
  dbcopy = *db;

  plen = strlen(path);
  if (strcmp(path+(plen-4),".dam") == 0)
    { root = Root(path,".dam");
      isdam = 1;
    }
  else
    { if (strcmp(path+(plen-3),".db") == 0)
        isdam = -1;
      else
        isdam = 0;
      root = Root(path,".db");
    }
  pwd = PathTo(path);

  bptr = rindex(root,'.');
  if (bptr != NULL && bptr[1] != '\0' && bptr[1] != '-')
    { part = strtol(bptr+1,&fptr,10);
      if (*fptr != '\0' || part == 0)
        part = 0;
      else
        *bptr = '\0';
    }
  else
    part = 0;

  if (isdam > 0)
    cat = MyCatenate(pwd,"/",root,".dam");
  else
    cat = MyCatenate(pwd,"/",root,".db");
  if (cat == NULL)
    return (-1);
  if ((dbvis = fopen(cat,"r")) == NULL)
    { if (isdam < 0)
        { EPRINTF(EPLACE,"%s: Could not open DB %s\n",Prog_Name,path);
          goto error;
        }
      if (isdam > 0)
        { EPRINTF(EPLACE,"%s: Could not open DAM %s\n",Prog_Name,path);
          goto error;
        }
      cat = MyCatenate(pwd,"/",root,".dam");
      if (cat == NULL)
        return (-1);
      if ((dbvis = fopen(cat,"r")) == NULL)
        { EPRINTF(EPLACE,"%s: Could not open %s as a DB or a DAM\n",Prog_Name,path);
          goto error;
        }
      isdam = 1;
    }
  if (isdam < 0)
    isdam = 0;

  if ((index = Fopen(MyCatenate(pwd,PATHSEP,root,".idx"),"r")) == NULL)
    goto error1;
  if (fread(db,sizeof(DAZZ_DB),1,index) != 1)
    { EPRINTF(EPLACE,"%s: Index file (.idx) of %s is junk\n",Prog_Name,root);
      goto error2;
    }

  { int   p, nblocks, nfiles;
    int64 size;
    char  fname[MAX_NAME], prolog[MAX_NAME];

    nblocks = 0;
    if (fscanf(dbvis,DB_NFILE,&nfiles) != 1)
      { EPRINTF(EPLACE,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
        goto error2;
      }
    for (p = 0; p < nfiles; p++)
      if (fscanf(dbvis,DB_FDATA,&tlast,fname,prolog) != 3)
        { EPRINTF(EPLACE,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
          goto error2;
        }
    if (fscanf(dbvis,DB_NBLOCK,&nblocks) != 1)
      if (part == 0)
        { cutoff = 0;
          all    = DB_ALL;
        }
      else
        { EPRINTF(EPLACE,"%s: DB %s has not yet been partitioned, cannot request a block !\n",
                         Prog_Name,root);
          goto error2;
        }
    else
      { if (fscanf(dbvis,DB_PARAMS,&size,&cutoff,&all) != 3)
          { EPRINTF(EPLACE,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
            goto error2;
          }
        if (part > nblocks)
          { EPRINTF(EPLACE,"%s: DB %s has only %d blocks\n",Prog_Name,root,nblocks);
            goto error2;
          }
      }

    if (part > 0)
      { for (p = 1; p <= part; p++)
          if (fscanf(dbvis,DB_BDATA,&ufirst,&tfirst) != 2)
            { EPRINTF(EPLACE,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
              goto error2;
            }
        if (fscanf(dbvis,DB_BDATA,&ulast,&tlast) != 2)
          { EPRINTF(EPLACE,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
            goto error2;
          }
      }
    else
      { ufirst = tfirst = 0;
        ulast  = db->ureads;
        tlast  = db->treads;
      }
  }

  db->trimmed = 0;
  db->tracks  = NULL;
  db->part    = part;
  db->cutoff  = cutoff;
  db->allarr |= all;
  db->ufirst  = ufirst;
  db->tfirst  = tfirst;

  nreads = ulast-ufirst;
  if (part <= 0)
    { db->reads = (DAZZ_READ *) Malloc(sizeof(DAZZ_READ)*(nreads+2),"Allocating Open_DB index");
      if (db->reads == NULL)
        goto error2;
        
      db->reads += 1;
      if (fread(db->reads,sizeof(DAZZ_READ),nreads,index) != (size_t) nreads)
        { EPRINTF(EPLACE,"%s: Index file (.idx) of %s is junk\n",Prog_Name,root);
          free(db->reads-1);
          goto error2;
        }
    }
  else
    { DAZZ_READ *reads;
      int        i, r, maxlen;
      int64      totlen;

      reads = (DAZZ_READ *) Malloc(sizeof(DAZZ_READ)*(nreads+2),"Allocating Open_DB index");
      if (reads == NULL)
        goto error2;
      reads += 1;

      fseeko(index,sizeof(DAZZ_READ)*ufirst,SEEK_CUR);
      if (fread(reads,sizeof(DAZZ_READ),nreads,index) != (size_t) nreads)
        { EPRINTF(EPLACE,"%s: Index file (.idx) of %s is junk\n",Prog_Name,root);
          free(reads-1);
          goto error2;
        }

      totlen = 0;
      maxlen = 0;
      for (i = 0; i < nreads; i++)
        { r = reads[i].rlen;
          totlen += r;
          if (r > maxlen)
            maxlen = r;
        }

      db->maxlen = maxlen;
      db->totlen = totlen;
      db->reads  = reads;
    }

  ((int *) (db->reads))[-1] = ulast - ufirst;   //  Kludge, need these for DB part
  ((int *) (db->reads))[-2] = tlast - tfirst;

  db->nreads = nreads;
  db->path   = Strdup(MyCatenate(pwd,PATHSEP,root,""),"Allocating Open_DB path");
  if (db->path == NULL)
    { free(db->reads-1);
      goto error2;
    }
  bases = Fopen(MyCatenate(db->path,"","",".bps"),"r");
  if (bases == NULL)
     { free(db->path);
       free(db->reads-1);
       goto error2;
     }
  db->bases = (void *) bases;
  db->loaded = 0;

  status = isdam;

error2:
  fclose(index);
error1:
  fclose(dbvis);
error:
  if (bptr != NULL)
    *bptr = '.';

  free(pwd);
  free(root);

  if (status < 0)
    *db = dbcopy;

  return (status);
}


// Trim the DB or part thereof and all opened tracks according to the cuttof and all settings
//   of the current DB partition.  Reallocate smaller memory blocks for the information kept
//   for the retained reads.

void Trim_DB(DAZZ_DB *db)
{ int         i, j, r, f;
  int         allflag, cutoff, css;
  int64       totlen;
  int         maxlen, nreads;
  DAZZ_TRACK *record;
  DAZZ_READ  *reads;

  if (db->trimmed) return;

  if (db->cutoff <= 0 && (db->allarr & DB_ALL) != 0) return;

  { int load_error;

    load_error = db->loaded;
    for (record = db->tracks; record != NULL; record = record->next)
      if (record->name == atrack_name)
        { if (((DAZZ_ARROW *) record)->loaded)
            load_error = 1;
        }
      else if (record->name != qtrack_name)
        { if (record->loaded)
            load_error = 1;
        }
    if (load_error)
      { EPRINTF(EPLACE,"%s: Cannot load anything before trim (Trim_DB)\n",Prog_Name);
        return;
      }
  }

  cutoff = db->cutoff;
  if ((db->allarr & DB_ALL) != 0)
    allflag = 0;
  else
    allflag = DB_BEST;

  reads  = db->reads;
  nreads = db->nreads;

  for (record = db->tracks; record != NULL; record = record->next)
    if (record->name == qtrack_name)
      { uint16 *table = ((DAZZ_QV *) record)->table;

        j = 0;
        for (i = 0; i < db->nreads; i++)
          if ((reads[i].flags & DB_BEST) >= allflag && reads[i].rlen >= cutoff)
            table[j++] = table[i];
      }
    else if (record->name == atrack_name)
      { DAZZ_ARROW *atrack = (DAZZ_ARROW *) record;
        int64      *aoff   = atrack->aoff;

        for (j = i = 0; i < nreads; i++)
          if ((reads[i].flags & DB_BEST) >= allflag && reads[i].rlen >= cutoff)
            aoff[j++] = aoff[i];
        atrack->aoff = Realloc(aoff,sizeof(int64)*j,NULL);
      }
    else
      { int size;

        size = record->size;
        if (record->data == NULL)
          { char *anno = (char *) record->anno;
            j = 0;
            for (i = r = 0; i < db->nreads; i++, r += size)
              if ((reads[i].flags & DB_BEST) >= allflag && reads[i].rlen >= cutoff)
                { memmove(anno+j,anno+r,size);
                  j += size;
                }
            record->anno = Realloc(record->anno,record->size*j,NULL);
          }
        else if (size == 4)
          { int *anno4 = (int *) (record->anno);
            int *alen  = record->alen;

            j = 0;
            for (i = 0; i < db->nreads; i++)
              if ((reads[i].flags & DB_BEST) >= allflag && reads[i].rlen >= cutoff)
                { anno4[j] = anno4[i];
                  alen[j]  = alen[i];
                  j += 1;
                }
            record->alen = Realloc(record->alen,sizeof(int)*j,NULL);
            record->anno = Realloc(record->anno,record->size*(j+1),NULL);
          }
        else // size == 8
          { int64 *anno8 = (int64 *) (record->anno);
            int *alen    = record->alen;

            j = 0;
            for (i = 0; i < db->nreads; i++)
              if ((reads[i].flags & DB_BEST) >= allflag && reads[i].rlen >= cutoff)
                { anno8[j] = anno8[i];
                  alen[j]  = alen[i];
                  j += 1;
                }
            record->alen = Realloc(record->alen,sizeof(int)*j,NULL);
            record->anno = Realloc(record->anno,record->size*(j+1),NULL);
          }
        record->nreads = j;
      }

  css    = 0;
  totlen = maxlen = 0;
  for (j = i = 0; i < nreads; i++)
    { f = reads[i].flags;
      if ((f & DB_CCS) == 0)
        css = 0;
      r = reads[i].rlen;
      if ((f & DB_BEST) >= allflag && r >= cutoff)
        { totlen += r;
          if (r > maxlen)
            maxlen = r;
          reads[j] = reads[i];
          if (css)
            reads[j++].flags |= DB_CCS;
          else
            reads[j++].flags &= ~DB_CCS;
          css = 1;
        }
    }
  
  db->totlen  = totlen;
  db->maxlen  = maxlen;
  db->nreads  = j;
  db->trimmed = 1;

  if (j < nreads)
    { db->reads = Realloc(reads-1,sizeof(DAZZ_READ)*(j+2),NULL);
      db->reads += 1;
    }
}


// Return the size in bytes of the memory occupied by a given DB

int64 sizeof_DB(DAZZ_DB *db)
{ int64       s;
  DAZZ_TRACK *t;

  s = sizeof(DAZZ_DB)
    + sizeof(DAZZ_READ)*(db->nreads+2)
    + strlen(db->path)+1
    + (db->totlen+db->nreads+4);

  t = db->tracks;
  if (t != NULL && strcmp(t->name,".@qvs") == 0)
    { DAZZ_QV *q = (DAZZ_QV *) t;
      s += sizeof(DAZZ_QV)
         + sizeof(uint16) * db->nreads
         + q->ncodes * sizeof(QVcoding)
         + 6;
      t = t->next;
    }

  for (; t != NULL; t = t->next)
    { s += sizeof(DAZZ_TRACK)
         + strlen(t->name)+1
         + t->size * (db->nreads+1);
      if (t->data != NULL)
        { if (t->size == 8)
            s += sizeof(int)*((int64 *) t->anno)[db->nreads];
          else //  t->size == 4
            s += sizeof(int)*((int *) t->anno)[db->nreads];
        }
    }

  return (s);
}


// For the DB or DAM "path" = "prefix/root.[db|dam]", find all the files for that DB, i.e. all
//   those of the form "prefix/[.]root.part" and call actor with the complete path to each file
//   pointed at by path, and the suffix of the path by extension.  The . proceeds the root
//   name if the defined constant HIDE_FILES is set.  Always the first call is with the
//   path "prefix/root.[db|dam]" and extension "db" or "dam".  There will always be calls for
//   "prefix/[.]root.idx" and "prefix/[.]root.bps".  All other calls are for *tracks* and
//   so this routine gives one a way to know all the tracks associated with a given DB.
//   -1 is returned if the path could not be found, and 1 is returned if an error (reported
//   to EPLACE) occured and INTERACTIVE is defined.  Otherwise a 0 is returned.

int List_DB_Files(char *path, void actor(char *path, char *extension))
{ int            status, plen, rlen, dlen;
  char          *root, *pwd, *name;
  int            isdam;
  DIR           *dirp;
  struct dirent *dp;

  status = 0;
  pwd    = PathTo(path);
  plen   = strlen(path);
  if (strcmp(path+(plen-4),".dam") == 0)
    root = Root(path,".dam");
  else
    root = Root(path,".db");
  rlen = strlen(root);

  if (root == NULL || pwd == NULL)
    { free(pwd);
      free(root);
      EXIT(1);
    }

  if ((dirp = opendir(pwd)) == NULL)
    { EPRINTF(EPLACE,"%s: Cannot open directory %s (List_DB_Files)\n",Prog_Name,pwd);
      status = -1;
      goto error;
    }

  isdam = 0;
  while ((dp = readdir(dirp)) != NULL)     //   Get case dependent root name (if necessary)
    { name = dp->d_name;
      if (strcmp(name,MyCatenate("","",root,".db")) == 0)
        break;
      if (strcmp(name,MyCatenate("","",root,".dam")) == 0)
        { isdam = 1;
          break;
        }
    }
  if (dp == NULL)
    { status = -1;
      closedir(dirp);
      goto error;
    }

  if (isdam)
    actor(MyCatenate(pwd,"/",root,".dam"),"dam");
  else
    actor(MyCatenate(pwd,"/",root,".db"),"db");

  rewinddir(dirp);                         //   Report each auxiliary file
  while ((dp = readdir(dirp)) != NULL)
    { name = dp->d_name;
      dlen = strlen(name);
#ifdef HIDE_FILES
      if (name[0] != '.')
        continue;
      dlen -= 1;
      name += 1;
#endif
      if (dlen < rlen+1)
        continue;
      if (name[rlen] != '.')
        continue;
      if (strncmp(name,root,rlen) != 0)
        continue;
      actor(MyCatenate(pwd,PATHSEP,name,""),name+(rlen+1));
    }
  closedir(dirp);

error:
  free(pwd);
  free(root);
  return (status);
}

void Print_Read(char *s, int width)
{ int i;

  if (s[0] < 4)
    { for (i = 0; s[i] != 4; i++)
        { if (i%width == 0 && i != 0)
            printf("\n");
          printf("%d",s[i]);
        }
      printf("\n");
    }
  else
    { for (i = 0; s[i] != '\0'; i++)
        { if (i%width == 0 && i != 0)
            printf("\n");
          printf("%c",s[i]);
        }
      printf("\n");
    }
}


// Shut down an open 'db' by freeing all associated space, including tracks and QV structures, 
//   and any open file pointers.  The record pointed at by db however remains (the user
//   supplied it and so should free it).

void Close_DB(DAZZ_DB *db)
{ if (db->loaded)
    free(((char *) (db->bases)) - 1);
  else if (db->bases != NULL)
    fclose((FILE *) db->bases);
  if (db->reads != NULL)
    free(db->reads-1);
  free(db->path);

  Close_QVs(db);

  Close_Arrow(db);

  while (db->tracks != NULL)
    Close_Track(db,db->tracks);
}


/*******************************************************************************************
 *
 *  READ AND ARROW BUFFER ALLOCATION, LOAD, & LOAD_ALL
 *
 ********************************************************************************************/

// Allocate and return a buffer big enough for the largest read in 'db', leaving room
//   for an initial delimiter character

char *New_Read_Buffer(DAZZ_DB *db)
{ char *read;

  read = (char *) Malloc(db->maxlen+4,"Allocating New Read Buffer");
  if (read == NULL)
    EXIT(NULL);
  return (read+1);
}

// Load into 'read' the i'th read in 'db'.  As an upper case ASCII string if ascii is 2, as a
//   lower-case ASCII string is ascii is 1, and as a numeric string over 0(A), 1(C), 2(G), and
//   3(T) otherwise.
//
// **NB**, the byte before read will be set to a delimiter character!

int Load_Read(DAZZ_DB *db, int i, char *read, int ascii)
{ FILE      *bases  = (FILE *) db->bases;
  int64      off;
  int        len, clen;
  DAZZ_READ *r = db->reads;

  if (i < 0 || i >= db->nreads)
    { EPRINTF(EPLACE,"%s: Index out of bounds (Load_Read)\n",Prog_Name);
printf("Index = %d\n",i);
      EXIT(1);
    }

  if (db->loaded)
    { len = r[i].rlen;
      strncpy(read,(char *) bases + r[i].boff,len);
      if (ascii == 0)
        { if (*read < 4)
            read[-1] = read[len] = 4;
          else
            { read[len] = '\0';
              Number_Read(read);
              read[-1] = 4;
            }
        }
      else
        { if (*read < 4)
            { read[len] = 4;
              if (ascii == 1)
                Lower_Read(read);
              else
                Upper_Read(read);
              read[-1] = '\0';
            }
          else
            { read[len] = '\0';
              if ((ascii == 1) != islower(*read))
                Change_Read(read);
            }
          read[-1] = '\0';
        }
      return (0);
    }

  off = r[i].boff;
  len = r[i].rlen;

  if (ftello(bases) != off)
    fseeko(bases,off,SEEK_SET);
  clen = COMPRESSED_LEN(len);
  if (clen > 0)
    { if (fread(read,clen,1,bases) != 1)
        { EPRINTF(EPLACE,"%s: Failed read of .bps file (Load_Read)\n",Prog_Name);
          EXIT(1);
        }
    }
  Uncompress_Read(len,read);
  if (ascii == 1)
    { Lower_Read(read);
      read[-1] = '\0';
    }
  else if (ascii == 2)
    { Upper_Read(read);
      read[-1] = '\0';
    }
  else
    read[-1] = 4;
  return (0);
}


// Load into 'read' the subread [beg,end] of the i'th read in 'db' and return a pointer to the
//   the start of the subinterval (not necessarily = to read !!! ).  As a lower case ascii
//   string if ascii is 1, an upper case ascii string if ascii is 2, and a numeric string
//   over 0(A), 1(C), 2(G), and 3(T) otherwise.  A '\0' (or 4) is prepended and appended to
//   the string holding the substring so it has a delimeter for traversals in either direction.
//   A NULL pointer is returned if an error occured and INTERACTIVE is defined.

char *Load_Subread(DAZZ_DB *db, int i, int beg, int end, char *read, int ascii)
{ FILE      *bases  = (FILE *) db->bases;
  int64      off;
  int        len, clen;
  int        bbeg, bend;
  DAZZ_READ *r = db->reads;

  if (i < 0 || i >= db->nreads)
    { EPRINTF(EPLACE,"%s: Index out of bounds (Load_Read)\n",Prog_Name);
      EXIT(NULL);
    }
    
  if (db->loaded)
    { len = end-beg;
      strncpy(read,(char *) bases + r[i].boff + beg,len);
      if (ascii == 0)
        { if (*read < 4)
            read[-1] = read[len] = 4;
          else
            { read[len] = '\0';
              Number_Read(read);
              read[-1] = 4;
            }
        }
      else
        { if (*read < 4)
            { read[len] = 4;
              if (ascii == 1)
                Lower_Read(read);
              else
                Upper_Read(read);
              read[-1] = '\0';
            }
          else
            { read[len] = '\0';
              if ((ascii == 1) != islower(*read))
                Change_Read(read);
            }
          read[-1] = '\0';
        }
      return (read);
    }

  bbeg = beg/4;
  bend = (end-1)/4+1;

  off = r[i].boff + bbeg;
  len = end - beg;

  if (ftello(bases) != off)
    fseeko(bases,off,SEEK_SET);
  clen = bend-bbeg;
  if (clen > 0)
    { if (fread(read,clen,1,bases) != 1)
        { EPRINTF(EPLACE,"%s: Failed read of .bps file (Load_Read)\n",Prog_Name);
          EXIT(NULL);
        }
    }
  Uncompress_Read(4*clen,read);
  read += beg%4;
  read[len] = 4;
  if (ascii == 1)
    { Lower_Read(read);
      read[-1] = '\0';
    }
  else if (ascii == 2)
    { Upper_Read(read);
      read[-1] = '\0';
    }
  else
    read[-1] = 4;

  return (read);
}

// Allocate a block big enough for all the uncompressed sequences, read them into it,
//   reset the 'off' in each read record to be its in-memory offset, and set the
//   bases pointer to point at the block after closing the bases file.  If ascii is
//   non-zero then the reads are converted to ACGT ascii, otherwise the reads are left
//   as numeric strings over 0(A), 1(C), 2(G), and 3(T).

int Load_All_Reads(DAZZ_DB *db, int ascii)
{ FILE      *bases = (FILE *) db->bases;
  int        nreads = db->nreads;
  DAZZ_READ *reads = db->reads;
  void     (*translate)(char *s);

  char  *seq;
  int64  o, off;
  int    i, len, clen;

  if (db->loaded)
    return (0);

  seq = (char *) Malloc(db->totlen+nreads+4,"Allocating All Sequence Reads");
  if (seq == NULL)
    EXIT(1);

  *seq++ = 4;

  if (ascii == 1)
    translate = Lower_Read;
  else
    translate = Upper_Read;

  o = 0;
  for (i = 0; i < nreads; i++)
    { len = reads[i].rlen;
      off = reads[i].boff;
      if (ftello(bases) != off)
        fseeko(bases,off,SEEK_SET);
      clen = COMPRESSED_LEN(len);
      if (clen > 0)
        { if (fread(seq+o,clen,1,bases) != 1)
            { EPRINTF(EPLACE,"%s: Read of .bps file failed (Load_All_Sequences)\n",Prog_Name);
              free(seq-1);
              EXIT(1);
            }
        }
      Uncompress_Read(len,seq+o);
      if (ascii)
        translate(seq+o);
      reads[i].boff = o;
      o += (len+1);
    }
  reads[nreads].boff = o;

  fclose(bases);

  db->bases  = (void *) seq;
  db->loaded = 1;

  return (0);
}


/*******************************************************************************************
 *
 *  ARROW OPEN, LOAD, LOAD_ALL, & CLOSE
 *
 ********************************************************************************************/

DAZZ_DB    *Arrow_DB = NULL;         //  Last db/arw used by "Load_Arrow"
DAZZ_ARROW *Arrow_Ptr;               //    Becomes invalid after closing

// If the Arrow pseudo track is not already in db's track list, then load it and set it up.
//   The database reads must not have been loaded with Load_All_Reads yet.
//   -1 is returned if a .arw file is not present, and 1 is returned if an error (reported
//   to EPLACE) occured and INTERACTIVE is defined.  Otherwise a 0 is returned.

int Open_Arrow(DAZZ_DB *db)
{ int64      *avector;
  DAZZ_ARROW *atrack;
  FILE       *afile;
  DAZZ_READ  *reads;
  int         i, nreads;

  if (db->tracks != NULL && db->tracks->name == atrack_name)
    return (0);
 
  if ((db->allarr & DB_ARROW) == 0)
    { EPRINTF(EPLACE,"%s: The DB is not an Arrow database (Open_Arrow)\n",Prog_Name);
      EXIT(1);
    }
  if (db->loaded)
    { EPRINTF(EPLACE,"%s: Cannot open Arrow vectors after loading all reads (Open_Arrow)\n",
                        Prog_Name);
      EXIT(1);
    }

  afile = Fopen(MyCatenate(db->path,"","",".arw"),"r");
  if (afile == NULL)
    return (-1);

  nreads  = db->nreads;
  avector = (int64 *) Malloc(sizeof(int64)*nreads,"Allocating Arrow index");
  atrack  = (DAZZ_ARROW *) Malloc(sizeof(DAZZ_ARROW),"Allocating Arrow track");
  if (avector == NULL || atrack == NULL)
    { fclose(afile);
      if (avector != NULL)
        free(avector);
      EXIT(1);
    }
  db->tracks     = (DAZZ_TRACK *) atrack;
  atrack->next   = NULL;
  atrack->name   = atrack_name;
  atrack->aoff   = avector;
  atrack->arrow  = (void *) afile;
  atrack->loaded = 0;


  reads = db->reads;
  for (i = 0; i < nreads; i++)
    avector[i] = reads[i].boff;
  return (0);
}

// Load into 'read' the i'th arrow in 'db'.  As an ASCII string if ascii is 1, 
//   and as a numeric string otherwise.

int Load_Arrow(DAZZ_DB *db, int i, char *arrow, int ascii)
{ FILE      *afile;
  int64      off;
  int        len, clen;

  if (db != Arrow_DB)
    { if (db->tracks == NULL || db->tracks->name != atrack_name)
        { EPRINTF(EPLACE,"%s: Arrow data is not available (Load_Arrow)\n",Prog_Name);
          EXIT(1);
        }
      Arrow_Ptr = (DAZZ_ARROW *) db->tracks;
      Arrow_DB  = db;
    }

  if (i < 0 || i >= db->nreads)
    { EPRINTF(EPLACE,"%s: Index out of bounds (Load_Arrow)\n",Prog_Name);
      EXIT(1);
    }

  afile = (FILE *) Arrow_Ptr->arrow;
  off   = Arrow_Ptr->aoff[i];
  len   = db->reads[i].rlen;

  if (ftello(afile) != off)
    fseeko(afile,off,SEEK_SET);
  clen = COMPRESSED_LEN(len);
  if (clen > 0)
    { if (fread(arrow,clen,1,afile) != 1)
        { EPRINTF(EPLACE,"%s: Failed read of .arw file (Load_Arrow)\n",Prog_Name);
          EXIT(1);
        }
    }
  Uncompress_Read(len,arrow);
  if (ascii == 1)
    { Letter_Arrow(arrow);
      arrow[-1] = '\0';
    }
  else
    arrow[-1] = 4;
  return (0);
}

// Allocate a block big enough for all the uncompressed Arrow vectors, read them into it,
//   reset the 'off' in each arrow record to be its in-memory offset, and set the
//   arrow pointer to point at the block after closing the arrow file.  If ascii is
//   non-zero then the arrows are converted to 0123 ascii, otherwise the arrows are left
//   as numeric strings over [0-3].

int Load_All_Arrows(DAZZ_DB *db, int ascii)
{ int        nreads = db->nreads;
  DAZZ_READ *reads = db->reads;
  FILE      *afile;
  int64     *aoff;

  char  *seq;
  int64  o, off;
  int    i, len, clen;

  if (db != Arrow_DB)
    { if (db->tracks == NULL || db->tracks->name != atrack_name)
        { EPRINTF(EPLACE,"%s: Arrow data is not available (Load_All_Arrows)\n",Prog_Name);
          EXIT(1);
        }
      Arrow_Ptr = (DAZZ_ARROW *) db->tracks;
      Arrow_DB  = db;
    }

  if (Arrow_Ptr->loaded)
    return (0);

  afile = (FILE *) Arrow_Ptr->arrow;
  aoff  = Arrow_Ptr->aoff;

  seq = (char *) Malloc(db->totlen+nreads+4,"Allocating All Arrows");
  if (seq == NULL)
    EXIT(1);

  *seq++ = 4;
  o = 0;
  for (i = 0; i < nreads; i++)
    { len = reads[i].rlen;
      off = aoff[i];
      if (ftello(afile) != off)
        fseeko(afile,off,SEEK_SET);
      clen = COMPRESSED_LEN(len);
      if (clen > 0)
        { if (fread(seq+o,clen,1,afile) != 1)
            { EPRINTF(EPLACE,"%s: Read of .bps file failed (Load_All_Sequences)\n",Prog_Name);
              free(seq-1);
              EXIT(1);
            }
        }
      Uncompress_Read(len,seq+o);
      if (ascii)
        Letter_Arrow(seq+o);
      aoff[i] = o;
      o += (len+1);
    }
  aoff[nreads] = o;

  fclose(afile);

  Arrow_Ptr->arrow  = (void *) seq;
  Arrow_Ptr->loaded = 1;

  return (0);
}

// Remove the Arrow pseudo track, all space associated with it, and close the .arw file.

void Close_Arrow(DAZZ_DB *db)
{ DAZZ_ARROW *atrack;

  Arrow_DB = NULL;
  if (db->tracks != NULL && db->tracks->name == atrack_name)
    { atrack = (DAZZ_ARROW *) db->tracks;
      if (atrack->loaded)
        free(atrack->arrow);
      else
        fclose((FILE *) atrack->arrow);
      free(atrack->aoff);
      db->tracks = db->tracks->next;
      free(atrack);
    }
}


/*******************************************************************************************
 *
 *  TRACK CHECK, OPEN, BUFFER ALLOCATION, LOAD, LOAD_ALL & CLOSE ROUTINES
 *     TRACK EXTRAS READING & WRITING
 *
 ********************************************************************************************/

//  Return status of track:
//     1: Track is for trimmed DB
//     0: Track is for untrimmed DB
//    -1: Track is not the right size of DB either trimmed or untrimmed
//    -2: Could not find the track 
//    -3: Error return (if INTERACTIVE mode only)

int Check_Track(DAZZ_DB *db, char *track, int *kind)
{ FILE       *afile;
  int         tracklen, size, ispart;
  int         ureads, treads;

  afile = NULL;
  if (db->part > 0)
    { afile  = fopen(MyCatenate(db->path,MyNumbered_Suffix(".",db->part,"."),track,".anno"),"r");
      ispart = 1;
    }
  if (afile == NULL)
    { afile  = fopen(MyCatenate(db->path,".",track,".anno"),"r");
      ispart = 0;
    }
  if (afile == NULL)
    return (-2);

  if (fread(&tracklen,sizeof(int),1,afile) != 1)
    { EPRINTF(EPLACE,"%s: track files for %s are corrupted\n",Prog_Name,track);
      fclose(afile);
      EXIT(-3);
    }
  if (fread(&size,sizeof(int),1,afile) != 1)
    { EPRINTF(EPLACE,"%s: track files for %s are corrupted\n",Prog_Name,track);
      fclose(afile);
      EXIT(-3);
    }

  if (size == 0)
    *kind = MASK_TRACK;
  else if (size > 0)
    *kind = CUSTOM_TRACK;
  else
    { EPRINTF(EPLACE,"%s: track files for %s are corrupted\n",Prog_Name,track);
      fclose(afile);
      EXIT(-3);
    }
  
  fclose(afile);

  if (ispart)
    { ureads = ((int *) (db->reads))[-1];
      treads = ((int *) (db->reads))[-2];
    }
  else
    { ureads = db->ureads;
      treads = db->treads;
    }

  if (tracklen == ureads)
    return (0);
  else if (tracklen == treads)
    return (1);
  else
    return (-1);
}

// The DB has already been trimmed, but a track over the untrimmed DB needs to be opened.
//   Trim the track by rereading the untrimmed DB index from the file system.

static int Late_Track_Trim(DAZZ_DB *db, DAZZ_TRACK *track, int ispart)
{ int         i, j, r;
  int         allflag, cutoff;
  int         ureads;
  char       *root;
  DAZZ_READ   read;
  FILE       *indx;

  if (db->cutoff <= 0 && (db->allarr & DB_ALL) != 0) return (0);

  cutoff = db->cutoff;
  if ((db->allarr & DB_ALL) != 0)
    allflag = 0;
  else
    allflag = DB_BEST;

  root = rindex(db->path,'/') + 2;
  indx = Fopen(MyCatenate(db->path,"","",".idx"),"r");
  fseeko(indx,sizeof(DAZZ_DB) + sizeof(DAZZ_READ)*db->ufirst,SEEK_SET);
  if (ispart)
    ureads = ((int *) (db->reads))[-1];
  else
    ureads = db->ureads;

  { int    size;

    size = track->size;
    if (track->data == NULL)
      { char *anno = (char *) track->anno;
        j = r = 0;
        for (i = r = 0; i < ureads; i++, r += size)
          { if (fread(&read,sizeof(DAZZ_READ),1,indx) != 1)
              { EPRINTF(EPLACE,"%s: Index file (.idx) of %s is junk\n",Prog_Name,root);
                fclose(indx);
                EXIT(1);
              }
            if ((read.flags & DB_BEST) >= allflag && read.rlen >= cutoff)
              { memmove(anno+j,anno+r,size);
                j += size;
              }
            r += size;
          }
        track->anno = Realloc(track->anno,track->size*j,NULL);
      }
    else if (size == 4)
      { int *anno4 = (int *) (track->anno);
        int *alen  = track->alen;

        j = 0;
        for (i = 0; i < ureads; i++)
          { if (fread(&read,sizeof(DAZZ_READ),1,indx) != 1)
              { EPRINTF(EPLACE,"%s: Index file (.idx) of %s is junk\n",Prog_Name,root);
                fclose(indx);
                EXIT(1);
              }
            if ((read.flags & DB_BEST) >= allflag && read.rlen >= cutoff)
              { anno4[j] = anno4[i];
                alen[j]  = alen[i];
                j += 1;
              }
          }
        track->alen = Realloc(track->alen,sizeof(int)*j,NULL);
        track->anno = Realloc(track->anno,track->size*(j+1),NULL);
      }
    else // size == 8
      { int64 *anno8 = (int64 *) (track->anno);
        int   *alen  = track->alen;

        j = 0;
        for (i = 0; i < ureads; i++)
          { if (fread(&read,sizeof(DAZZ_READ),1,indx) != 1)
              { EPRINTF(EPLACE,"%s: Index file (.idx) of %s is junk\n",Prog_Name,root);
                fclose(indx);
                EXIT(1);
              }
            if ((read.flags & DB_BEST) >= allflag && read.rlen >= cutoff)
              { anno8[j] = anno8[i];
                alen[j]  = alen[i];
                j += 1;
              }
          }
        track->alen = Realloc(track->alen,sizeof(int)*j,NULL);
        track->anno = Realloc(track->anno,track->size*(j+1),NULL);
      }
  }

  fclose(indx);
  return (0);
}

// If track is not already in the db's track list, then allocate all the storage for it,
//   read it in from the appropriate file, add it to the track list, and return a pointer
//   to the newly created DAZZ_TRACK record.  If the track does not exist or cannot be
//   opened for some reason, then NULL is returned.

DAZZ_TRACK *Open_Track(DAZZ_DB *db, char *track)
{ FILE       *afile, *dfile;
  int         tracklen, size;
  int         nreads, ispart;
  int         treads, ureads;
  int64       dmax;
  void       *anno;
  int        *alen;
  void       *data;
  char       *name;
  DAZZ_TRACK *record;

  if (track[0] == '.')
    { EPRINTF(EPLACE,"%s: Track name, '%s', cannot begin with a .\n",Prog_Name,track);
      EXIT(NULL);
    }

  for (record = db->tracks; record != NULL; record = record->next)
    if (strcmp(record->name,track) == 0)
      return (record);

  afile = NULL;
  if (db->part)
    { afile  = fopen(MyCatenate(db->path,MyNumbered_Suffix(".",db->part,"."),track,".anno"),"r");
      ispart = 1;
    }
  if (afile == NULL)
    { afile = fopen(MyCatenate(db->path,".",track,".anno"),"r");
      ispart = 0;
    }
  if (afile == NULL)
    { EPRINTF(EPLACE,"%s: Track '%s' does not exist\n",Prog_Name,track);
      return (NULL);
    }

  dfile  = NULL;
  anno   = NULL;
  alen   = NULL;
  data   = NULL;
  record = NULL;

  if (ispart)
    name = MyCatenate(db->path,MyNumbered_Suffix(".",db->part,"."),track,".data");
  else
    name = MyCatenate(db->path,".",track,".data");
  if (name == NULL)
    goto error;
  dfile = fopen(name,"r");

  if (fread(&tracklen,sizeof(int),1,afile) != 1)
    { EPRINTF(EPLACE,"%s: Track '%s' annotation file is junk\n",Prog_Name,track);
      goto error;
    }
  if (fread(&size,sizeof(int),1,afile) != 1)
    { EPRINTF(EPLACE,"%s: Track '%s' annotation file is junk\n",Prog_Name,track);
      goto error;
    }

  if (size < 0)
    { EPRINTF(EPLACE,"%s: Track '%s' annotation file is junk\n",Prog_Name,track);
      goto error;
    }
  if (size == 0)
    size = 8;

  if (ispart)
    { ureads = ((int *) (db->reads))[-1];
      treads = ((int *) (db->reads))[-2];
    }
  else
    { ureads = db->ureads;
      treads = db->treads;
    }

  if (db->trimmed)
    { if (tracklen != treads && tracklen != ureads)
        { EPRINTF(EPLACE,"%s: Track '%s' not same size as database !\n",Prog_Name,track);
          goto error;
        }
      if ( ! ispart && db->part > 0)
        { if (tracklen == treads)
            fseeko(afile,size*db->tfirst,SEEK_CUR);
          else
            fseeko(afile,size*db->ufirst,SEEK_CUR);
        }
    }
  else
    { if (tracklen != ureads)
        { if (tracklen == treads)
            EPRINTF(EPLACE,"%s: Track '%s' is for a trimmed DB !\n",Prog_Name,track);
          else
            EPRINTF(EPLACE,"%s: Track '%s' not same size as database !\n",Prog_Name,track);
          goto error;
        }
      if ( ! ispart && db->part > 0)
        fseeko(afile,size*db->ufirst,SEEK_CUR);
    }
  if (tracklen == treads)
    nreads = ((int *) (db->reads))[-2];
  else
    nreads = ((int *) (db->reads))[-1];

  anno = (void *) Malloc(size*(nreads+1),"Allocating Track Anno Vector");
  if (anno == NULL)
    goto error;

  if (dfile != NULL)
    { int64 *anno8;
      int   *anno4;
      int64  x, y;
      int    i;

      alen = (int *)  Malloc(sizeof(int)*nreads,"Allocating Track Anno Lengths");
      if (alen == NULL)
        goto error;

      if (fread(anno,size,nreads+1,afile) != (size_t) (nreads+1))
        { EPRINTF(EPLACE,"%s: Track '%s' annotation file is junk\n",Prog_Name,track);
          goto error;
        }

      dmax = 0;
      if (size == 4)
        { anno4 = (int *) anno;
          y = anno4[0];
          for (i = 1; i <= nreads; i++)
            { x = anno4[i];
              y = x-y;
              if (y > dmax)
                dmax = y; 
              alen[i-1] = y;
              y = x;
            }
        }
      else
        { anno8 = (int64 *) anno;
          y = anno8[0];
          for (i = 1; i <= nreads; i++)
            { x = anno8[i];
              y = x-y;
              if (y > dmax)
                dmax = y; 
              alen[i-1] = y;
              y = x;
            }
        }
    }
  else
    { dmax = 0;
      if (fread(anno,size,nreads,afile) != (size_t) nreads)
        { EPRINTF(EPLACE,"%s: Track '%s' annotation file is junk\n",Prog_Name,track);
          goto error;
        }
    }

  fclose(afile);

  record = (DAZZ_TRACK *) Malloc(sizeof(DAZZ_TRACK),"Allocating Track Record");
  if (record == NULL)
    goto error;
  record->name = Strdup(track,"Allocating Track Name");
  if (record->name == NULL)
    goto error;
  if (dfile == NULL)
    record->data = NULL;
  else
    record->data = (void *) dfile;
  record->anno   = anno;
  record->alen   = alen;
  record->size   = size;
  record->nreads = nreads;
  record->loaded = 0;
  record->dmax   = dmax;

  if (db->trimmed && tracklen != treads)
    { if (Late_Track_Trim(db,record,ispart))
        goto error;
    }

  if (db->tracks != NULL && (db->tracks->name == qtrack_name || db->tracks->name == atrack_name))
    { record->next     = db->tracks->next;
      db->tracks->next = record;
    }
  else
    { record->next = db->tracks;
      db->tracks   = record;
    }

  return (record);

error:
  if (record != NULL)
    free(record);
  if (data != NULL)
    free(data);
  if (alen != NULL)
    free(alen);
  if (anno != NULL)
    free(anno);
  if (dfile != NULL)
    fclose(dfile);
  fclose(afile);
  EXIT (NULL);
}

// Allocate a data buffer large enough to hold the longest read data block that will occur
//   in the track.  If cannot allocate memory then return NULL if INTERACTIVE is defined,
//   or print error to stderr and exit otherwise.

void *New_Track_Buffer(DAZZ_TRACK *track)
{ void *data;

  data = (void *) Malloc(track->dmax,"Allocating New Track Data Buffer");
  if (data == NULL)
    EXIT(NULL);
  return (data);
}

// Load into 'data' the read data block for read i's "track" data.  Return the length of
//   the data in bytes, unless an error occurs and INTERACTIVE is defined in which case
//   return wtih -1.

int Load_Track_Data(DAZZ_TRACK *track, int i, void *data)
{ FILE      *dfile;
  int64      off;
  int        len;

  if (i < 0 || i >= track->nreads)
    { EPRINTF(EPLACE,"%s: Index out of bounds (Load_Track_Data)\n",Prog_Name);
      EXIT(-1);
    }

  if (track->size == 4)
    off = ((int *) track->anno)[i];
  else
    off = ((int64 *) track->anno)[i];
  len = track->alen[i];

  if (track->loaded)
    { strncpy(data,(void *) track->data + off,len);
      return (len);
    }

  dfile = (FILE *) track->data;
  if (ftello(dfile) != off)
    fseeko(dfile,off,SEEK_SET);
  if (len > 0)
    if (fread(data,len,1,dfile) != 1)
      { EPRINTF(EPLACE,"%s: Failed read of .data file (Load_Track_Data)\n",Prog_Name);
        EXIT(-1);
      }
  return (len);
}

// Allocate a block big enough for all the track data and read the data into it,
//   reset the 'off' in each anno pointer to be its in-memory offset, and set the
//   data pointer to point at the block after closing the data file.  Return with a
//   zero, except when an error occurs and INTERACTIVE is defined in which
//   case return wtih 1.

int Load_All_Track_Data(DAZZ_TRACK *track)
{ FILE  *dfile;
  void  *data;
  int   *alen;
  int64  dlen, off, o;
  int    i, len, nreads;

  if (track->loaded || track->data == NULL)
    return (0);

  nreads = track->nreads;
  dfile  = (FILE *) track->data;
  alen   = track->alen;

  dlen = 0;
  for (i = 0; i < nreads; i++)
    dlen += alen[i];

  data = (void *) Malloc(dlen,"Allocating All Track Data");
  if (data == NULL)
    EXIT(1);

  o = 0;
  if (track->size == 4)
    { int *anno4 = (int *) track->anno;

      for (i = 0; i < nreads; i++)
        { len = alen[i];
          off = anno4[i];
          if (ftello(dfile) != off)
            fseeko(dfile,off,SEEK_SET);
          if (len > 0)
            { if (fread(data+o,len,1,dfile) != 1)
                { EPRINTF(EPLACE,"%s: Read of .data failed (Load_All_Track_Data)\n",Prog_Name);
                  free(data);
                  EXIT(1);
                }
            }
          anno4[i] = o;
          o += len;
        }
      anno4[nreads] = o;
    }
  else
    { int64 *anno8 = (int64 *) track->anno;

      for (i = 0; i < nreads; i++)
        { len = alen[i];
          off = anno8[i];
          if (ftello(dfile) != off)
            fseeko(dfile,off,SEEK_SET);
          if (len > 0)
            { if (fread(data+o,len,1,dfile) != 1)
                { EPRINTF(EPLACE,"%s: Read of .data failed (Load_All_Track_Data)\n",Prog_Name);
                  free(data);
                  EXIT(1);
                }
            }
          anno8[i] = o;
          o += len;
        }
      anno8[nreads] = o;
    }

  fclose(dfile);

  track->data = (void *) data;
  track->loaded = 1;

  return (0);
}


// Assumming file pointer for afile is correctly positioned at the start of a extra item,
//   and aname is the name of the .anno file, decode the value present and places it in
//   extra if extra->nelem == 0, otherwise reduce the value just read into extra according 
//   according the to the directive given by 'accum'.  Leave the read poinrt at the next
//   extra or end-of-file.
//   Returns:
//      1 if at the end of file,
//      0 if item was read and folded correctly,
//     -1 if there was a system IO or allocation error (if interactive), and
//     -2 if the new value could not be reduced into the currenct value of extra (interactive)

int Read_Extra(FILE *afile, char *aname, DAZZ_EXTRA *extra)
{ int   vtype, nelem, accum, slen;
  char *name;
  void *value;

#define EREAD(v,s,n,file,ret)                                                           \
  { if (fread(v,s,n,file) != (size_t) n)                                                \
      { if (ferror(file))                                                               \
          EPRINTF(EPLACE,"%s: System error, read failed!\n",Prog_Name);       		\
        else if (ret)                                                                   \
          return (1);									\
        else										\
          EPRINTF(EPLACE,"%s: The file %s is corrupted\n",Prog_Name,aname);		\
        EXIT(-1);									\
      }                                                                                 \
  }

  EREAD(&vtype,sizeof(int),1,afile,1)
  EREAD(&nelem,sizeof(int),1,afile,0)
  EREAD(&accum,sizeof(int),1,afile,0)
  EREAD(&slen,sizeof(int),1,afile,0)

  if (extra == NULL)
    { if (fseeko(afile,slen+8*nelem,SEEK_CUR) < 0)
        { EPRINTF(EPLACE,"%s: System error, read failed!\n",Prog_Name);
          EXIT(-1);
        }
      return (0);
    }

  name  = (char *) Malloc(slen+1,"Allocating extra name");
  value = Malloc(8*nelem,"Allocating extra value");
  if (value == NULL || name == NULL)
    { free(name);
      free(value);
      EXIT(-1);
    }

  EREAD(name,1,slen,afile,0);
  EREAD(value,8,nelem,afile,0);
  name[slen] = '\0';
  
  if (extra->nelem == 0)
    { extra->vtype = vtype;
      extra->nelem = nelem;
      extra->accum = accum;
      extra->name  = name;
      extra->value = value;
      return (0);
    }

  if (vtype != extra->vtype)
    { EPRINTF(EPLACE,"%s: Type of extra %s does not agree with previous .anno block files\n",
                     Prog_Name,name);
      goto error;
    }
  if (nelem != extra->nelem)
    { EPRINTF(EPLACE,"%s: Length of extra %s does not agree with previous .anno block files\n",
                     Prog_Name,name);
      goto error;
    }
  if (accum != extra->accum)
    { EPRINTF(EPLACE,
           "%s: Reduction indicator of extra %s does not agree with previos .anno block files\n",
           Prog_Name,name);
      goto error;
    }
  if (strcmp(name,extra->name) != 0)
    { EPRINTF(EPLACE,"%s: Expecting extra %s in .anno block file, not %s\n",
                     Prog_Name,extra->name,name);
      goto error;
    }

  if (vtype == DB_INT)
    { int64 *ival = (int64 *) value;
      int64 *eval = (int64 *) (extra->value);
      int    j;

      if (accum == DB_EXACT)
        { for (j = 0; j < nelem; j++)
            if (eval[j] != ival[j])
              { EPRINTF(EPLACE,
                    "%s: Value of extra %s doe not agree with previous .anno block files\n",
                    Prog_Name,name);
                goto error;
              }
        }
      else
        { for (j = 0; j < nelem; j++)
            eval[j] += ival[j];
        }
    }

  else
    { double *ival = (double *) value;
      double *eval = (double *) (extra->value);
      int     j;

      if (accum == DB_EXACT)
        { for (j = 0; j < nelem; j++)
            if (eval[j] != ival[j])
              { EPRINTF(EPLACE,
                    "%s: Value of extra %s doe not agree with previous .anno block files\n",
                    Prog_Name,name);
                goto error;
              }
        }
      else
        { for (j = 0; j < nelem; j++)
            eval[j] += ival[j];
        }
    }

  free(value);
  free(name);
  return (0);

error:
  free(value);
  free(name);
  EXIT(-2);
}

//  Write extra record to end of file afile and advance write pointer
//  If interactive, then return non-zero on error, if bash, then print
//  and halt if an error

int Write_Extra(FILE *afile, DAZZ_EXTRA *extra)
{ int slen; 

  FFWRITE(&(extra->vtype),sizeof(int),1,afile)
  FFWRITE(&(extra->nelem),sizeof(int),1,afile)
  FFWRITE(&(extra->accum),sizeof(int),1,afile)
  slen = strlen(extra->name);
  FFWRITE(&slen,sizeof(int),1,afile)
  FFWRITE(extra->name,1,slen,afile)
  FFWRITE(extra->value,8,extra->nelem,afile)

  return (0);
}

void Close_Track(DAZZ_DB *db, DAZZ_TRACK *track)
{ DAZZ_TRACK *record, *prev;

  prev = NULL;
  for (record = db->tracks; record != NULL; record = record->next)
    { if (track == record)
        { free(record->anno);
          free(record->alen);
          if (record->loaded)
            free(record->data);
          else
            fclose((FILE *) record->data);
          free(record->name);
          if (prev == NULL)
            db->tracks = record->next;
          else
            prev->next = record->next;
          free(record);
          return;
        }
      prev = record;
    }
  return;
}


/*******************************************************************************************
 *
 *  QV OPEN, BUFFER ALLOCATION, LOAD, & CLOSE ROUTINES
 *
 ********************************************************************************************/

DAZZ_DB *Active_DB = NULL;  //  Last db/qv used by "Load_QVentry"
DAZZ_QV *Active_QV;         //    Becomes invalid after closing

int Open_QVs(DAZZ_DB *db)
{ FILE        *quiva, *istub, *indx;
  char        *root;
  uint16      *table;
  DAZZ_QV     *qvtrk;
  QVcoding    *coding, *nx;
  int          ncodes = 0;

  if (db->tracks != NULL && db->tracks->name == qtrack_name)
    return (0);

  if (db->trimmed)
    { EPRINTF(EPLACE,"%s: Cannot load QVs after trimming the DB\n",Prog_Name);
      EXIT(1);
    }

  if (db->reads[db->nreads-1].coff < 0)
    { if (db->part > 0)
        { EPRINTF(EPLACE,"%s: All QVs for this block have not been added to the DB!\n",Prog_Name);
          EXIT(1);
        }
      else
        { EPRINTF(EPLACE,"%s: All QVs for this DB have not been added!\n",Prog_Name);
          EXIT(1);
        }
    }

  //  Open .qvs, .idx, and .db files

  quiva = Fopen(MyCatenate(db->path,"","",".qvs"),"r");
  if (quiva == NULL)
    return (-1);

  istub  = NULL;
  indx   = NULL; 
  table  = NULL;
  coding = NULL;
  qvtrk  = NULL;

  root = rindex(db->path,'/');
  if (root[1] == '.')
    { *root = '\0';
      istub = Fopen(MyCatenate(db->path,"/",root+2,".db"),"r");
      *root = '/';
    }
  else
    istub = Fopen(MyCatenate(db->path,"","",".db"),"r");
  if (istub == NULL)
    goto error;

  { int   first, last, nfiles;
    char  prolog[MAX_NAME], fname[MAX_NAME];
    int   i, j;

    if (fscanf(istub,DB_NFILE,&nfiles) != 1)
      { EPRINTF(EPLACE,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
        goto error;
      }

    if (db->part > 0)
      { int       pfirst, plast;
        int       fbeg, fend;
        int       n, k;
        FILE     *indx;

        //  Determine first how many and which files span the block (fbeg to fend)

        pfirst = db->ufirst;
        plast  = pfirst + db->nreads;

        first = 0;
        for (fbeg = 0; fbeg < nfiles; fbeg++)
          { if (fscanf(istub,DB_FDATA,&last,fname,prolog) != 3)
              { EPRINTF(EPLACE,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
                goto error;
              }
            if (last > pfirst)
              break;
            first = last;
          }
        for (fend = fbeg+1; fend <= nfiles; fend++)
          { if (last >= plast)
              break;
            if (fscanf(istub,DB_FDATA,&last,fname,prolog) != 3)
              { EPRINTF(EPLACE,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
                goto error;
              }
            first = last;
          }

        indx   = Fopen(MyCatenate(db->path,"","",".idx"),"r");
        ncodes = fend-fbeg;
        coding = (QVcoding *) Malloc(sizeof(QVcoding)*ncodes,"Allocating coding schemes");
        table  = (uint16 *) Malloc(sizeof(uint16)*db->nreads,"Allocating QV table indices");
        if (indx == NULL || coding == NULL || table == NULL)
          { ncodes = 0;
            goto error;
          }

        //  Carefully get the first coding scheme (its offset is most likely in a DAZZ_RECORD
        //    in .idx that is *not* in memory).  Get all the other coding schemes normally and
        //    assign the tables # for each read in the block in "tables".

        rewind(istub);
        (void) fscanf(istub,DB_NFILE,&nfiles);

        first = 0;
        for (n = 0; n < fbeg; n++)
          { (void) fscanf(istub,DB_FDATA,&last,fname,prolog);
            first = last;
          }

        for (n = fbeg; n < fend; n++)
          { (void) fscanf(istub,DB_FDATA,&last,fname,prolog);

            i = n-fbeg;
            if (first < pfirst)
              { DAZZ_READ read;

                fseeko(indx,sizeof(DAZZ_DB) + sizeof(DAZZ_READ)*first,SEEK_SET);
                if (fread(&read,sizeof(DAZZ_READ),1,indx) != 1)
                  { EPRINTF(EPLACE,"%s: Index file (.idx) of %s is junk\n",Prog_Name,root);
                    ncodes = i;
                    goto error;
                  }
                fseeko(quiva,read.coff,SEEK_SET);
                nx = Read_QVcoding(quiva);
                if (nx == NULL)
                  { ncodes = i;
                    goto error;
                  }
                coding[i] = *nx;
              }
            else
              { fseeko(quiva,db->reads[first-pfirst].coff,SEEK_SET);
                nx = Read_QVcoding(quiva);
                if (nx == NULL)
                  { ncodes = i;
                    goto error;
                  }
                coding[i] = *nx;
                db->reads[first-pfirst].coff = ftello(quiva);
              }

            j = first-pfirst;
            if (j < 0)
              j = 0;
            k = last-pfirst;
            if (k > db->nreads)
              k = db->nreads;
            while (j < k)
              table[j++] = (uint16) i;

            first = last;
	  }

        fclose(indx);
        indx = NULL;
      }

    else
      { //  Load in coding scheme for each file, adjust .coff of first read in the file, and
        //    record which table each read uses

        ncodes = nfiles;
        coding = (QVcoding *) Malloc(sizeof(QVcoding)*nfiles,"Allocating coding schemes");
        table  = (uint16 *) Malloc(sizeof(uint16)*db->nreads,"Allocating QV table indices");
        if (coding == NULL || table == NULL)
          goto error;
  
        first = 0;
        for (i = 0; i < nfiles; i++)
          { if (fscanf(istub,DB_FDATA,&last,fname,prolog) != 3)
              { EPRINTF(EPLACE,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
                goto error;
              }
  
            fseeko(quiva,db->reads[first].coff,SEEK_SET);
            nx = Read_QVcoding(quiva);
            if (nx == NULL)
              { ncodes = i;
                goto error;
              }
            coding[i] = *nx;
	    db->reads[first].coff = ftello(quiva);

            for (j = first; j < last; j++)
              table[j] = (uint16) i;

            first = last;
          }
      }

    //  Allocate and fill in the DAZZ_QV record and add it to the front of the
    //    track list

    qvtrk = (DAZZ_QV *) Malloc(sizeof(DAZZ_QV),"Allocating QV pseudo-track");
    if (qvtrk == NULL)
      goto error;
    qvtrk->name   = qtrack_name;
    if (qvtrk->name == NULL)
      goto error;
    qvtrk->next   = db->tracks;
    db->tracks    = (DAZZ_TRACK *) qvtrk;
    qvtrk->ncodes = ncodes;
    qvtrk->table  = table;
    qvtrk->coding = coding;
    qvtrk->quiva  = quiva;
  }

  fclose(istub);
  return (0);

error:
  if (qvtrk != NULL)
    free(qvtrk);
  if (table != NULL)
    free(table);
  if (coding != NULL)
    { int i;
      for (i = 0; i < ncodes; i++)
        Free_QVcoding(coding+i);
      free(coding);
    }
  if (indx != NULL)
    fclose(indx);
  if (istub != NULL)
    fclose(istub);
  fclose(quiva);
  EXIT(1);
}

// Allocate and return a buffer of 5 vectors big enough for the largest read in 'db'

char **New_QV_Buffer(DAZZ_DB *db)
{ char **entry;
  char  *qvs;
  int    i;

  qvs   = (char *) Malloc(db->maxlen*5,"Allocating New QV Buffer");
  entry = (char **) Malloc(sizeof(char *)*5,"Allocating New QV Buffer");
  if (qvs == NULL || entry == NULL)
    EXIT(NULL);
  for (i = 0; i < 5; i++)
    entry[i] = qvs + i*db->maxlen;
  return (entry);
}

// Load into entry the QV streams for the i'th read from db.  The parameter ascii applies to
//  the DELTAG stream as described for Load_Read.

int Load_QVentry(DAZZ_DB *db, int i, char **entry, int ascii)
{ DAZZ_READ *reads;
  FILE      *quiva;
  int        rlen;

  if (db != Active_DB)
    { if (db->tracks == NULL || strcmp(db->tracks->name,".@qvs") != 0)
        { EPRINTF(EPLACE,"%s: QV's have not been opened (Load_QVentry)\n",Prog_Name);
          EXIT(1);
        }
      Active_QV = (DAZZ_QV *) db->tracks;
      Active_DB = db;
    }

  if (i < 0 || i >= db->nreads)
    { EPRINTF(EPLACE,"%s: Index out of bounds (Load_QVentry)\n",Prog_Name);
      EXIT(1);
    }

  reads = db->reads;
  quiva = Active_QV->quiva;
  rlen  = reads[i].rlen;

  fseeko(quiva,reads[i].coff,SEEK_SET);
  if (Uncompress_Next_QVentry(quiva,entry,Active_QV->coding+Active_QV->table[i],rlen))
    EXIT(1);

  if (ascii != 1)
    { char *deltag = entry[1];

      if (ascii != 2)
        { char x = deltag[rlen];
          deltag[rlen] = '\0';
          Number_Read(deltag);
          deltag[rlen] = x;
        }
      else
        { int j;
          int u = 'A'-'a';

          for (j = 0; j < rlen; j++)
            deltag[j] = (char) (deltag[j]+u);
        }
    }

  return (0);
}

// Close the QV stream, free the QV pseudo track and all associated memory

void Close_QVs(DAZZ_DB *db)
{ DAZZ_TRACK *track;
  DAZZ_QV    *qvtrk;
  int         i;

  Active_DB = NULL;

  track = db->tracks;
  if (track != NULL && strcmp(track->name,".@qvs") == 0)
    { qvtrk = (DAZZ_QV *) track;
      for (i = 0; i < qvtrk->ncodes; i++)
        Free_QVcoding(qvtrk->coding+i);
      free(qvtrk->coding);
      free(qvtrk->table);
      fclose(qvtrk->quiva);
      db->tracks = track->next;
      free(track);
    }
  return;
}


/*******************************************************************************************
 *
 *  COMMAND LINE @-EXPANSION PARSER
 *    Take a command line argument and interpret the '@' block number ranges.
 *    Parse_Block_Arg produces an Block_Looper iterator object that can then
 *    be invoked multiple times to iterate through all the files implied by
 *    the @ pattern/range.
 *
 ********************************************************************************************/

typedef struct
  { int first, last, next;
    char *root, *pwd, *ppnt;
    int   isDB;
    char *slice;
  } _Block_Looper;

//  Advance the iterator e_parse to the next file, open it, and return the file pointer
//   to it.  Return NULL if at the end of the list of files.

int Next_Block_Exists(Block_Looper *e_parse)
{ _Block_Looper *parse = (_Block_Looper *) e_parse;

  char       *disp;
  struct stat sts;

  if (parse->isDB)
    { if (parse->next+1 > parse->last)
        return (0);
      else
        return (1);
    }

  if (parse->next+1 > parse->last)
    return (0);

  if (parse->next < 0)
    disp = parse->root;
  else
    disp = MyNumbered_Suffix(parse->root,parse->next+1,parse->ppnt);

  if (stat(MyCatenate(parse->pwd,"/",disp,".las"),&sts))
    return (0);
  else
    return (1);
}


FILE *Next_Block_Arg(Block_Looper *e_parse)
{ _Block_Looper *parse = (_Block_Looper *) e_parse;

  char *disp;
  FILE *input;

  if (parse->isDB)
    { fprintf(stderr,"%s: Cannot open a DB block as a file (Next_Block_Arg)\n",Prog_Name);
      exit (1);  //  exit even in interactive mode as this is a programming bug on
    }            //   the part of the caller

  parse->next += 1;
  if (parse->next > parse->last)
    return (NULL);

  if (parse->next < 0)
    disp  = parse->root;
  else
    disp = MyNumbered_Suffix(parse->root,parse->next,parse->ppnt);

  if ((input = fopen(MyCatenate(parse->pwd,"/",disp,".las"),"r")) == NULL)
    { if (parse->last != INT_MAX)
        { EPRINTF(EPLACE,"%s: %s.las is not present\n",Prog_Name,disp);
          EXIT(NULL);
        }
      return (NULL);
    }
  return (input);
}

//  Reset the iterator e_parse to the first file

void Reset_Block_Arg(Block_Looper *e_parse)
{ _Block_Looper *parse = (_Block_Looper *) e_parse;

  parse->next = parse->first - 1;
}

//  Advance the iterator e_parse to the next file

int Advance_Block_Arg(Block_Looper *e_parse)
{ _Block_Looper *parse = (_Block_Looper *) e_parse;

  if (Next_Block_Exists(e_parse))
    { parse->next += 1;
      return (1);
    }
  else
    return (0);
}

//  Return a pointer to the path for the current file

char *Block_Arg_Path(Block_Looper *e_parse)
{ _Block_Looper *parse = (_Block_Looper *) e_parse;

  return (Strdup(parse->pwd,"Allocating block path"));
}

//  Return a pointer to the root name for the current file

char *Block_Arg_Root(Block_Looper *e_parse)
{ _Block_Looper *parse = (_Block_Looper *) e_parse;
  char *name;

  if (parse->next < 0)
    name = parse->root;
  else
    name = MyNumbered_Suffix(parse->root,parse->next,parse->ppnt);
  return (Strdup(name,"Allocating block root"));
}

//  Free the iterator

void Free_Block_Arg(Block_Looper *e_parse)
{ _Block_Looper *parse = (_Block_Looper *) e_parse;

  free(parse->root);
  free(parse->pwd);
  free(parse->slice);
  free(parse);
}

char *Next_Block_Slice(Block_Looper *e_parse, int slice)
{ _Block_Looper *parse = (_Block_Looper *) e_parse;

  if (parse->slice == NULL)
    { int size = strlen(parse->pwd) + strlen(Block_Arg_Root(parse)) + 30;
      parse->slice =  (char *)  Malloc(size,"Block argument slice");
      if (parse->slice == NULL)
        EXIT(NULL);
    }

  if (parse->next+1 > parse->last)
    return (NULL);
  if (parse->next+slice > parse->last)
    slice = parse->last-parse->next;

  if (parse->first < 0)
    sprintf(parse->slice,"%s/%s",parse->pwd,parse->root);
  else
    sprintf(parse->slice,"%s/%s%c%d-%d%s",parse->pwd,parse->root,BLOCK_SYMBOL,parse->next+1,
                                          parse->next+slice,parse->ppnt);
  parse->next += slice;
  return (parse->slice);
}

//  Parse the command line argument and return an iterator to move through the
//    file names, setting it up to report the first file.

static Block_Looper *parse_block_arg(char *arg, int isDB)
{ _Block_Looper *parse;
  char *pwd, *root;
  char *ppnt, *cpnt;
  int   first, last;

  parse = (_Block_Looper *) Malloc(sizeof(_Block_Looper),"Allocating parse node");
  pwd   = PathTo(arg);
  if (isDB)
    { int len = strlen(arg);
      if (strcmp(arg+(len-4),".dam") == 0)
        { root = Root(arg,".dam");
          isDB = 2;
        }
      else
        root = Root(arg,".db");
    }
  else
    root  = Root(arg,".las");
  if (parse == NULL || pwd == NULL || root == NULL)
    goto error;

  ppnt = index(root,BLOCK_SYMBOL);
  if (ppnt == NULL)
    first = last = -1;
  else
    { if (index(ppnt+1,BLOCK_SYMBOL) != NULL)
        { EPRINTF(EPLACE,"%s: Two or more occurences of %c-sign in source name '%s'\n",
                         Prog_Name,BLOCK_SYMBOL,root);
          goto error;
        }
      *ppnt++ = '\0';
      first = strtol(ppnt,&cpnt,10);
      if (cpnt == ppnt)
        { first = 1;
          last  = INT_MAX;
        }
      else
        { if (first < 1)
            { EPRINTF(EPLACE,
                      "%s: Integer following %c-sigan is less than 1 in source name '%s'\n",
                      Prog_Name,BLOCK_SYMBOL,root);
              goto error;
            }
          if (*cpnt == '-')
            { ppnt = cpnt+1;
              last = strtol(ppnt,&cpnt,10);
              if (cpnt == ppnt)
                { EPRINTF(EPLACE,"%s: Second integer must follow - in source name '%s'\n",
                                 Prog_Name,root);
                  goto error;
                }
              if (last < first)
                { EPRINTF(EPLACE,
                          "%s: 2nd integer is less than 1st integer in source name '%s'\n",
                          Prog_Name,root);
                  goto error;
                }
              ppnt = cpnt;
            }
          else
            { last = INT_MAX;
              ppnt = cpnt;
            }
        }
    }

  parse->pwd   = pwd;
  parse->root  = root;
  parse->ppnt  = ppnt;
  parse->first = first;
  parse->last  = last;
  parse->next  = first-1;
  parse->slice = NULL;
  parse->isDB  = isDB;

  if (isDB && first >= 0 && last == INT_MAX)
    { char  buffer[2*MAX_NAME+100];
      char *dbname;
      FILE *dbfile;
      int   i, nfiles, nblocks;

      dbname = MyCatenate(pwd,"/",root,"db"); 
      dbfile = fopen(dbname,"r");
      if (dbfile == NULL)
        { dbname = MyCatenate(pwd,"/",root,"dam"); 
          dbfile = fopen(dbname,"r");
          if (dbfile == NULL)
            { EPRINTF(EPLACE,"%s: Cannot open database %s[db|dam]\n",Prog_Name,root);
              goto error;
            }
        }

      if (fscanf(dbfile,DB_NFILE,&nfiles) != 1)
        SYSTEM_READ_ERROR
      for (i = 0; i < nfiles; i++)
        if (fgets(buffer,2*MAX_NAME+100,dbfile) == NULL)
          SYSTEM_READ_ERROR
      if (fscanf(dbfile,DB_NBLOCK,&nblocks) != 1)
        SYSTEM_READ_ERROR
      fclose(dbfile);

      parse->last = nblocks;
    }

  return ((Block_Looper *) parse);

error:
  free(parse);
  free(root);
  free(pwd);
  EXIT(NULL);
}

Block_Looper *Parse_Block_LAS_Arg(char *arg)
{ return (parse_block_arg(arg, 0)); }

Block_Looper *Parse_Block_DB_Arg(char *arg)
{ return (parse_block_arg(arg, 1)); }

#include "QV.c"
