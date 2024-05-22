/*******************************************************************************************
 *
 *  Identify the type and full path name of a DNA source file and it GDB target
 *
 *  Author:   Gene Myers
 *            Modified by Richard Durbin to also read .1seq files
 *  Origin:   Complete rework/specialization from fasta2DAM of DAZZ_DB for FASTGA module
 *  Creation: Jan 2024
 *  Last Mod: Mar 2024
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <zlib.h>
#include <dirent.h>

#include "DB.h"
#include "alncode.h"
#include "DNAsource.h"

static OneSchema *Schema;

static char *find_1code(char *path, char *root)
{ DIR           *dirp;
  struct dirent *dp;
  OneFile       *of;
  int            rlen;
  char          *extn, *name;

  dirp = opendir(path);
  if (dirp == NULL)
    { fprintf(stderr,"%s: Cannot open directory %s\n",Prog_Name,path);
      exit (1);
    }

  extn = NULL;
  rlen = strlen(root);
  while ((dp = readdir(dirp)) != NULL)
    { name = dp->d_name;
      if (strncmp(root,name,rlen) != 0)
        continue;
      if (name[rlen] != '.')
        continue;
      of = oneFileOpenRead(Catenate(path,"/",name,""),Schema,"seq",1);
      if (of == NULL)
        continue;
      if (strcmp(of->fileType,"seq") != 0 || of->info['S'] == NULL)
        { oneFileClose(of);
          continue;
        }
      if (extn != NULL)
        { fprintf(stderr,"%s: Two 1-code sequence files with root %s found\n",Prog_Name,root);
          exit (1);
        }
      extn = Strdup(name+rlen,"Allocating extension");
      oneFileClose(of);
    }

  return (extn);
}

//  Interpret source & target arguments returning complete path + extension names of
//    the source and target in spath & tpath.  Returns the type of source.
//  If target == source then tpath is a name for a temporary file
//  If target == NULL then tpath has the same path & root as spath
//  Else If target is a path then tpath has the same root as spath
//  If no_gdb then the source cannot be a gdb, otherwise if the source is a gdb
//    then tpath = spath on return.

int get_dna_paths(char *source, char *target, char **spath, char **tpath, int no_gdb)
{ char   *SROOT, *SPATH;
  char   *TROOT, *TPATH;
  char   *SEXTN;
  char   *suffix[9] = { "", ".fa", ".fna", ".fasta",
                        ".gz", ".fa.gz", ".fna.gz", ".fasta.gz", ".gdb" };
  int     suflen[9] = { 0, 3, 4, 6, 3, 6, 7, 9, 4 };

  char *p, *s, *e;
  int   i, exists, tmax;
  FILE *input;
  struct stat status;

  if (no_gdb)
    tmax = 7;
  else
    tmax = 8;

  SPATH = source;
  p = rindex(SPATH,'/');
  if (p == NULL)
    { SROOT = SPATH;
      SPATH = ".";
    }
  else
    { *p++ = '\0';
      SROOT = p;
    }

  Schema = make_Aln_Schema();
  if (Schema == NULL)
    { fprintf(stderr,"%s: Failed to create onecode schema\n",Prog_Name);
      exit (1);
    }

  for (i = tmax; i >= 0; i--)
    { input = fopen(Catenate(SPATH,"/",SROOT,suffix[i]),"r");
      if (input != NULL)
        break;
    }
  if (i < 0)
    { SEXTN = find_1code(SPATH,SROOT);
      if (SEXTN == NULL)
        { fprintf(stderr,"%s: Could not find %s or a fasta or 1-code extension there of\n",
                         Prog_Name,SROOT);
          exit (1);
        }
    }
  else
    { SEXTN = suffix[i];
      fclose(input);
    }

  if (i == 4)
    { e = SROOT + strlen(SROOT);
      for (i = 1; i < 4; i++)
        if (strcmp(e-suflen[i],suffix[i]) == 0)
          break;
      if (i >= 4)
        { fprintf(stderr,"%s: Could not find valid extension of %s\n",Prog_Name,SROOT);
          exit (1);
        }
      e[-suflen[i]] = '\0';
      i += 4;
      SEXTN = suffix[i];
    }
  else if (i == 0)
    { e = SROOT + strlen(SROOT);
      for (i = 1; i <= tmax; i++)
        if (strcmp(e-suflen[i],suffix[i]) == 0 && i != 4)
          break;
      if (i <= tmax)
        { e[-suflen[i]] = '\0';
          SEXTN = suffix[i];
        }
      else
        { OneFile *of;

          of = oneFileOpenRead(Catenate(SPATH,"/",SROOT,""),Schema,"seq",1);
          if (of == NULL)
            { fprintf(stderr,"%s: Could not find valid extension of %s\n",Prog_Name,SROOT);
              exit (1);
            }
          oneFileClose(of);
          p = rindex(SROOT,'.');
          if (p != NULL)
            { SEXTN = Strdup(p,"Allocating extension");
              *p = '\0';
            }
          else
            SEXTN = "";
          i = -1;
        }
    }

  oneSchemaDestroy(Schema);

  *spath = Strdup(Catenate(SPATH,"/",SROOT,SEXTN),"Allocating source name");

  if (i == 8)
    { *tpath = Strdup(*spath,"Allocating target name");
      return (IS_GDB);
    }

  if (target == source)
    { TPATH = ".";
      TROOT = Strdup(Catenate(".",SROOT,".",Numbered_Suffix("",getpid(),"")),
                     "Allocating temp name");
    }
  else if (target == NULL)
    { TPATH = SPATH;
      TROOT = SROOT;
    }
  else // argc == 3
    { TPATH = target;
      exists = (stat(TPATH,&status) >= 0);
      if ( !exists || (status.st_mode & S_IFMT) != S_IFDIR)
        { p = rindex(TPATH,'/');
          if (p == NULL)
            { TROOT = TPATH;
              TPATH = ".";
            }
          else
            { *p++ = '\0';
              TROOT = p;
            }
          s = rindex(TROOT,'.');
          if (s != NULL)
            { if (strcmp(s+1,"gdb") == 0)
                *s = '\0';
              else if (exists)
                { fprintf(stderr,"%s: Existing target %s must be a .gdb\n",Prog_Name,TROOT);
                  exit (1);
                }
            }
          else if (exists)
            { fprintf(stderr,"%s: Existing target %s must be a .gdb\n",Prog_Name,TROOT);
              exit (1);
            }
        }
      else
        TROOT = SROOT;
    }

  *tpath = Strdup(Catenate(TPATH,"/",TROOT,".gdb"),"Allocating target name");

  if (i < 0)
    return (IS_ONE);
  else if (i < 4)
    return (IS_FA);
  else
    return (IS_FA_GZ);
}

FILE **make_temporary_gdb(char *spath, char *tpath, DAZZ_DB *db, int *hdrs, int NTHREADS)
{ char  *command;
  FILE **bps;
  int    i;

  if (NTHREADS <= 0)
    NTHREADS = 1;

  command = Malloc(strlen(spath)+strlen(tpath)+100,"Allocating command string");
  if (command == NULL)
    exit (1);
  sprintf(command,"FAtoGDB %s %s",spath,tpath);
  if (system(command) != 0)
    { fprintf(stderr,"\n%s: Call to FAtoGDB failed\n",Prog_Name);
      exit (1);
    }
  free(command);

  if (Open_DB(tpath,db) < 0)
    goto clean;

  if (NTHREADS > 1)
    { bps = Malloc(NTHREADS*sizeof(FILE *),"Allocating bps IO");
      if (bps == NULL)
        goto clean;

      for (i = 1; i < NTHREADS; i++)
        { bps[i] = fopen(Catenate(db->path,".bps","",""),"r");
          if (bps[i] == NULL)
            { fprintf(stderr,"%s: Cannot open another copy of DB %s\n",Prog_Name,tpath);
              goto clean;
            }
        }
    }
  else
    bps = NULL;

  *hdrs = open(Catenate(db->path,".hdr","",""),O_RDONLY);
  if (*hdrs < 0)
    goto clean;

  unlink(Catenate(db->path,".hdr","",""));
  unlink(Catenate(db->path,".bps","",""));
  unlink(Catenate(db->path,".idx","",""));
  unlink(tpath);

  return (bps);

clean:
  sprintf(command,"GDBrm %s",tpath);
  system(command);
  exit (1);
}
