/********************************************************************************************
 *
 *  Remove a list of .db databases
 *     Delete all the files for the given data bases <name>.db ... (there are a couple
 *     of hidden . files for each DB, and these are removed too.)  Do not use "rm" to
 *     remove a database.
 *
 *  Author:  Gene Myers
 *  Date  :  July 2013
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "gene_core.h"

static char *Usage = "[-vifg] <source:path>[.1gdb|.gix] ... ";

int main(int argc, char *argv[])
{ int   VERBOSE;
  int   ASK;
  int   FORCE;
  int   GDB_TOO;
  char *command;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];

    ARG_INIT("GIXrm")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("vifg") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE  = flags['v'];
    ASK      = flags['i'];
    FORCE    = flags['f'];
    GDB_TOO  = flags['g'];

    if (argc <= 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, list what is being deleted.\n");
        fprintf(stderr,"      -i: prompt for each (stub) deletion\n");
        fprintf(stderr,"      -f: force operation quietly\n");
        fprintf(stderr,"      -g: Also delete the associated GDB.\n");
        exit (1);
      }
    if (FORCE)
      ASK = VERBOSE = 0;
  }

  //  Determine source and target root names, paths, and extensions

  { char *PATH, *ROOT, *GEXTN;
    int   HAS_GDB, HAS_GIX;
    char *p, *com;
    FILE *input;
    int   c, a, yes;

    if (FORCE)
      com = "rm -f";
    else
      com = "rm";

    for (c = 1; c < argc; c++)
      { HAS_GDB = HAS_GIX = 0;

        PATH = argv[c];
        p = rindex(PATH,'/');
        if (p == NULL)
          { ROOT = PATH;
            PATH = ".";
          }
        else
          { *p++ = '\0';
            ROOT = p;
          }

        input = fopen(Catenate(PATH,"/",ROOT,".1gdb"),"r");
        if (input != NULL)
          { fclose(input);
            HAS_GDB = 1;
            GEXTN   = ".1gdb";
          }
        else if (strcmp(ROOT+(strlen(ROOT)-5),".1gdb") == 0)
          { HAS_GDB = 1;
            ROOT[strlen(ROOT)-5] = '\0';
            GEXTN   = ".1gdb";
          }
        else
          { input = fopen(Catenate(PATH,"/",ROOT,".gdb"),"r");
            if (input != NULL)
              { fclose(input);
                HAS_GDB = 1;
                GEXTN   = ".gdb";
              }
            else if (strcmp(ROOT+(strlen(ROOT)-4),".gdb") == 0)
              { HAS_GDB = 1;
                ROOT[strlen(ROOT)-4] = '\0';
                GEXTN   = ".gdb";
              }
          }

        input = fopen(Catenate(PATH,"/",ROOT,".gix"),"r");
        if (input != NULL)
          { fclose(input);
            HAS_GIX = 1;
          }
        else if (strcmp(ROOT+(strlen(ROOT)-4),".gix") == 0)
          { HAS_GIX = 1;
            ROOT[strlen(ROOT)-4] = '\0';
          }

        command = Malloc(4*(strlen(PATH)+strlen(ROOT))+100,"Allocating command buffer");

        if (HAS_GIX + HAS_GDB * GDB_TOO == 0)
          { if (GDB_TOO)
              fprintf(stderr,"  Warning: there is no GDB or GIX with root %s/%s\n",PATH,ROOT);
            else
              fprintf(stderr,"  Warning: there is no GIX with root %s/%s\n",PATH,ROOT);
          }

        if (HAS_GIX)
          { yes = 1;
            if (ASK)
              { printf("Remove %s/%s.gix? ",PATH,ROOT);
                fflush(stdout);
                yes = 0;
                while ((a = getc(stdin)) != '\n')
                  if (a == 'y' || a == 'Y')
                    yes = 1;
                  else if (a == 'n' || a == 'N')
                    yes = 0;
                fflush(stdout);
              }
            if (yes)
              { if (VERBOSE)
                  { fprintf(stderr,"  Removing %s/%s.gix\n",PATH,ROOT);
                    fflush(stderr);
                  }
                sprintf(command,"%s %s/%s.gix %s/.%s.ktab.* %s/.%s.post.*",
                                com,PATH,ROOT,PATH,ROOT,PATH,ROOT);
                if (system(command) != 0) goto sys_error;
              }
          }

        if (HAS_GDB && GDB_TOO)
          { yes = 1;
            if (ASK)
              { printf("Remove %s/%s%s? ",PATH,ROOT,GEXTN);
                fflush(stdout);
                yes = 0;
                while ((a = getc(stdin)) != '\n')
                  if (a == 'y' || a == 'Y')
                    yes = 1;
                  else if (a == 'n' || a == 'N')
                    yes = 0;
                fflush(stdout);
              }
            if (yes)
              { if (VERBOSE)
                  { fprintf(stderr,"  Removing %s/%s%s\n",PATH,ROOT,GEXTN);
                    fflush(stderr);
                  }
                sprintf(command,"%s %s/%s%s %s/.%s.bps",com,PATH,ROOT,GEXTN,PATH,ROOT);
                if (system(command) != 0) goto sys_error;
              }
          }

        free(command);
      }
  }

  exit (0);

sys_error:
  fprintf(stderr,"%s: Cannot execute command '%s'\n",Prog_Name,command);
  exit (1);
}
