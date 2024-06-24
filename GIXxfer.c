/********************************************************************************************
 *
 *  Move or copy a GDB/GIX ensemble
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "gene_core.h"
#include "GDB.h"

static char *Usage = "[-vinfx] <source:path>[.1gdb|.gix] <target:path>[.1gdb|.gix]";

int main(int argc, char *argv[])
{ int   VERBOSE;
  int   ASK;
  int   NO_OVER;
  int   FORCE;
  char *op, *command;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];

#ifdef MOVE
    ARG_INIT("GIXmv")
#else
    ARG_INIT("GIXcp")
#endif

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("vinf") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE  = flags['v'];
    ASK      = flags['i'];
    FORCE    = flags['f'];
    NO_OVER  = flags['n'];

    if (argc != 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, list what is being deleted.\n");
        fprintf(stderr,"      -i: prompt for each deletion\n");
        fprintf(stderr,"      -n: do not overwrite existing files.\n");
        fprintf(stderr,"      -f: force operation quietly\n");
        exit (1);
      }
    if (FORCE)
      ASK = VERBOSE = 0;

#ifdef MOVE
    if (FORCE)
      op = "mv -f";
    else
      op = "mv";
#else
    if (FORCE)
      op = "cp -f";
    else
      op = "cp";
#endif
  }

  //  Determine source and target root names, paths, and extensions

  { char *SPATH, *SROOT, *SEXTN;
    char *TPATH, *TROOT, *TEXTN;
    int   HAS_GDB, HAS_GIX;
    GDB  _gdb, *gdb = &_gdb;
    struct stat status;
    char  com[10];
    int   a, yes;
    char *p;
    FILE *input;

    if (FORCE)
      sprintf(com,"mv -fn?");
    else
      sprintf(com,"mv");

    HAS_GDB = HAS_GIX = 0;

    SPATH = argv[1];
    p = rindex(SPATH,'/');
    if (p == NULL)
      { SROOT = SPATH;
        SPATH = ".";
      }
    else
      { *p++ = '\0';
        SROOT = p;
      }

    input = fopen(Catenate(SPATH,"/",SROOT,".1gdb"),"r");
    if (input != NULL)
      { fclose(input);
        HAS_GDB = 1;
        SEXTN   = ".1gdb";
      }
    else if (strcmp(SROOT+(strlen(SROOT)-5),".1gdb") == 0)
      { HAS_GDB = 1;
        SROOT[strlen(SROOT)-5] = '\0';
        SEXTN   = ".1gdb";
      }
    else
      { input = fopen(Catenate(SPATH,"/",SROOT,".gdb"),"r");
        if (input != NULL)
          { fclose(input);
            HAS_GDB = 1;
            SEXTN   = ".gdb";
          }
        else if (strcmp(SROOT+(strlen(SROOT)-4),".gdb") == 0)
          { HAS_GDB = 1;
            SROOT[strlen(SROOT)-4] = '\0';
            SEXTN   = ".gdb";
          }
      }

    input = fopen(Catenate(SPATH,"/",SROOT,".gix"),"r");
    if (input != NULL)
      { fclose(input);
        HAS_GIX = 1;
      }
    else if (strcmp(SROOT+(strlen(SROOT)-4),".gix") == 0)
      { HAS_GIX = 1;
        SROOT[strlen(SROOT)-4] = '\0';
      }

    TEXTN = ".1gdb";
    TPATH = argv[2];
    if (stat(TPATH,&status) == 0 && (status.st_mode & S_IFMT) == S_IFDIR)
      TROOT = SROOT;
    else
      { p = rindex(TPATH,'/');
        if (p == NULL)
          { TROOT = TPATH;
            TPATH = ".";
          }
        else
          { *p++ = '\0';
            TROOT = p;
          }
        p = rindex(TROOT,'.');
        if (p != NULL)
          { if (strcmp(p+1,"gdb") == 0)
              { *p = '\0';
                TEXTN = ".gdb";
              }
            else if (strcmp(p+1,"1gdb") == 0 || strcmp(p+1,"gix") == 0)
              *p = '\0';
          }
      }

    command = Malloc(strlen(SPATH)+strlen(SROOT)+strlen(TPATH)+strlen(TROOT)+100,
                    "Allocating command buffer");

    if (HAS_GIX + HAS_GDB == 0)
      fprintf(stderr,"  There is no GDB or GIX with root %s/%s\n",SPATH,SROOT);

    if (HAS_GIX)
      { yes = 1;
        if (ASK)
          { yes = 0;
#ifdef MOVE
            printf("Move %s/%s.gix? ",SPATH,SROOT);
#else
            printf("Copy %s/%s.gix? ",SPATH,SROOT);
#endif
            fflush(stdout);
            while ((a = getc(stdin)) != '\n')
              if (a == 'y' || a == 'Y')
                yes = 1;
              else if (a == 'n' || a == 'N')
                yes = 0;
            fflush(stdout);
          }
        if (yes)
          { if (!NO_OVER || stat(Catenate(TPATH,"/",TROOT,".gix"),&status) != 0)
              { if (VERBOSE)
#ifdef MOVE
                  { fprintf(stderr,"  Moving %s/%s.gix to %s/%s.gix\n",SPATH,SROOT,TPATH,TROOT);
#else
                  { fprintf(stderr,"  Copying %s/%s.gix to %s/%s.gix\n",SPATH,SROOT,TPATH,TROOT);
#endif
                    fflush(stderr);
                  }
                sprintf(command,"%s %s/%s.gix %s/%s.gix",op,SPATH,SROOT,TPATH,TROOT);
                if (system(command) != 0) goto sys_error;
                for (a = 1; 1; a++)
                 { if (stat(Catenate(SPATH,"/.",SROOT,Numbered_Suffix(".ktab.",a,"")),&status) != 0)
                      break;
                   sprintf(command,"%s %s/.%s.ktab.%d %s/.%s.ktab.%d",
                                   op,SPATH,SROOT,a,TPATH,TROOT,a);
                   if (system(command) != 0) goto sys_error;
                   sprintf(command,"%s %s/.%s.post.%d %s/.%s.post.%d",
                                   op,SPATH,SROOT,a,TPATH,TROOT,a);
                   if (system(command) != 0) goto sys_error;
                 }
              }
          }
      }

    if (HAS_GDB)
      { yes = 1;
        if (ASK)
          { yes = 0;
#ifdef MOVE
            printf("Move %s/%s%s? ",SPATH,SROOT,SEXTN);
#else
            printf("Copy %s/%s%s? ",SPATH,SROOT,SEXTN);
#endif
            fflush(stdout);
            while ((a = getc(stdin)) != '\n')
              if (a == 'y' || a == 'Y')
                yes = 1;
              else if (a == 'n' || a == 'N')
                yes = 0;
            fflush(stdout);
          }
        if (yes)
          { if (!NO_OVER || stat(Catenate(TPATH,"/",TROOT,SEXTN),&status) != 0)
              { if (VERBOSE)
#ifdef MOVE
                  { fprintf(stderr,"  Moving %s/%s%s to %s/%s%s\n",
                                   SPATH,SROOT,SEXTN,TPATH,TROOT,TEXTN);
#else
                  { fprintf(stderr,"  Copying %s/%s%s to %s/%s%s\n",
                                   SPATH,SROOT,SEXTN,TPATH,TROOT,TEXTN);
#endif
                    fflush(stderr);
                  }
                Read_GDB(gdb,Catenate(SPATH,"/",SROOT,SEXTN));
#ifdef MOVE
                sprintf(command,"rm %s/%s%s",SPATH,SROOT,SEXTN);
                if (system(command) != 0) goto sys_error;
#endif

                sprintf(command,"%s %s/.%s.bps %s/.%s.bps",op,SPATH,SROOT,TPATH,TROOT);
                if (system(command) != 0) goto sys_error;

                Write_GDB(gdb,Catenate(TPATH,"/",TROOT,TEXTN));
                Close_GDB(gdb);
              }
          }
      }

    free(command);
  }

  exit (0);

sys_error:
  fprintf(stderr,"%s: Cannot execute command '%s'\n",Prog_Name,command);
  exit (1);
}
