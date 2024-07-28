/********************************************************************************************
 *
 *  Recreate the .fasta file that gave rise to a given .gdb
 *
 *  Author:  Gene Myers
 *  Origin:  Reworking of DAZZ_DB program DAM2fasta.c for the FASTGA module
 *  Date  :  Jan 2024
 *
 ********************************************************************************************/
/*  Last edited: Jul 25 19:12 2024 (rd109) */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <zlib.h>

#include "GDB.h"
#include "ONElib.h"

extern bool addProvenance(OneFile *of, OneProvenance *from, int n) ; // backdoor - clean up some day

static char *Usage =
         "[-vU] [-w<int(80)>] <source:path>[.1gdb] [ @ | <target:path>[<fa_extn>|.1seq] ]";

int main(int argc, char *argv[])
{ GDB        _gdb, *gdb = &_gdb;
  char       *TPATH, *TROOT;
  int         TEXTN;
  int         is_one;
  OneSchema  *schema;
  OneFile    *outone;
  int         gzip;
  void       *output;

  int UPPER;   // -U
  int WIDTH;   // -w
  int VERBOSE; // -v

  char *suffix[7] = { ".1seq", ".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz" };
  int   suflen[7] = { 5, 3, 4, 6, 6, 7, 9 };

  //  Process arguments

  { int   i, j, k;
    int   flags[128];
    char *eptr;

    ARG_INIT("GDBtoFA")

    WIDTH = 80;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vU")
            break;
          case 'w':
            ARG_NON_NEGATIVE(WIDTH,"Line width")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (flags['U'])
      UPPER = UPPER_CASE;
    else
      UPPER = LOWER_CASE;
    VERBOSE = flags['v'];

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -U: Use upper case for DNA (default is lower case).\n");
        fprintf(stderr,"      -w: Print -w bp per line (default is 80).\n");
        exit (1);
      }
  }

  //  Open db & determine source and target parts

  { struct stat tdesc;
    char       *p, *e, *root;
    int         i;

    Read_GDB(gdb,argv[1]);

    e = argv[1] + strlen(argv[1]);
    if (strcmp(e-5,".1gdb") == 0)
      root = Root(argv[1],".1gdb");
    else
      root = Root(argv[1],".gdb");

    is_one = 0;
    if (argc == 2)
      { if (gdb->seqsrc == IS_ONE)
          is_one = 1;
        else
          gzip   = 0;
        TPATH  = NULL;
        TROOT  = NULL;
        TEXTN  = 0;
      }
    else
      { if (strcmp(argv[2],"@") == 0)
          { TPATH = gdb->srcpath;
            p = rindex(TPATH,'/');
            *p++ = '\0';
            TROOT = p; 
          }
        else
          { TPATH = argv[2];
            if (stat(TPATH,&tdesc) >= 0 && (tdesc.st_mode & S_IFMT) != S_IFDIR)
              TROOT = root;
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
              }
          }

        e = TROOT + strlen(TROOT);
        for (i = 6; i >= 0; i--)
          if (strcmp(e-suflen[i],suffix[i]) == 0)
            break;

        if (i >= 0)
          { TEXTN = i;
            e[-suflen[i]] = '\0';
          }
        else
          { e = gdb->srcpath + strlen(gdb->srcpath);
            for (i = 6; i >= 0; i--)
              if (strcmp(e-suflen[i],suffix[i]) == 0)
                break;
            TEXTN = i;
          }
        is_one = (TEXTN == 0);
        gzip   = (TEXTN >= 4);
      }

    free(root);

    if (VERBOSE && argc > 2)
      { if (is_one)
          { if (strcmp(TPATH,".") == 0)
              fprintf(stderr,"\n  Creating 1-seq file %s%s at current directory\n",
                             TROOT,suffix[TEXTN]);
            else
              fprintf(stderr,"\n  Creating 1-seq file %s%s at directory %s\n",
                             TROOT,suffix[TEXTN],TPATH);
          }

        else
          { if (strcmp(TPATH,".") == 0)
              fprintf(stderr,"\n  Creating%sfasta file %s%s at current directory\n",
                             TEXTN>=4?" gzip'd ":" ",TROOT,suffix[TEXTN]);
            else
              fprintf(stderr,"\n  Creating%sfasta file %s%s at directory %s\n",
                             TEXTN>=4?" gzip'd ":" ",TROOT,suffix[TEXTN],TPATH);
	  }
        fflush(stderr);
      }
  }

  { GDB_CONTIG   *ctg;
    GDB_SCAFFOLD *scf;
    char         *contig, *header;

    scf    = gdb->scaffolds;
    ctg    = gdb->contigs;
    contig = New_Contig_Buffer(gdb);

    if (is_one)

      { int64 spos;
        int   len;
        int   i, k;

        schema = make_Seq_Schema();
        if (schema == NULL)
          { fprintf(stderr,"%s: Failed to create ONEcode schema\n",Prog_Name);
            exit (1);
          }
        if (argc == 2)
          { outone = oneFileOpenWriteNew("-",schema,"seq",0,0);
            if (outone == NULL)
              { fprintf(stderr,"%s: Cannot open ONEcode for writing to stdout\n",Prog_Name);
                exit (1);
              }
          }
        else
          { outone = oneFileOpenWriteNew(Catenate(TPATH,"/",TROOT,".1seq"),schema,"seq",1,0);
            if (outone == NULL)
              { fprintf(stderr,"%s: Cannot open %s/%s.%s for writing\n",
                               Prog_Name,TPATH,TROOT,suffix[TEXTN]);
                exit (1);
              }
          }

        addProvenance(outone,gdb->prov,gdb->nprov);
        oneAddProvenance(outone,Prog_Name,"0.1",Command_Line);

        if (!outone->isBinary)
          { int64 seqcnt, seqmax, seqtot;
            int64 gapcnt, gapmax, gaptot;
            int64 scfcnt, scfmax, scftot;
            int64 seqgrpcnt, gapgrpcnt, scnt, gcnt;
            int64 seqgrptot, gapgrptot, stot, gtot;
            OneStat *s;

            seqcnt = seqmax = seqtot = 0;
            gapcnt = gapmax = gaptot = 0;
            scfcnt = scfmax = scftot = 0;
            seqgrpcnt = gapgrpcnt = 0;
            seqgrptot = gapgrptot = 0;

            for (i = 0; i < gdb->nscaff; i++)
              { len = strlen(gdb->headers + scf[i].hoff);
                if (len > scfmax)
                  scfmax = len;
                scfcnt += 1;
                scftot += len;

                gtot = gcnt = 0;
                stot = scnt = 0;
                spos = 0;
                for (k = scf[i].fctg; k < scf[i].ectg; k++)
                  { if (ctg[k].sbeg > spos)
                      { len = ctg[k].sbeg-spos;
                        if (len > gapmax)
                          gapmax = len;
                        gtot += len;
                        gcnt += 1;
                      }
                    spos = ctg[k].sbeg;

                    len = ctg[k].clen;
                    if (len > seqmax)
                      seqmax = len;
                    stot += len;
                    scnt += 1;
                    spos += len;
                  }
                if (spos < scf[i].slen)
                  { len = scf[i].slen-spos;
                    if (len > gapmax)
                      gapmax = len;
                    gtot += len;
                    gcnt += 1;
                  }

                if (stot > seqgrptot)
                  seqgrptot = stot;
                if (scnt > seqgrpcnt)
                  seqgrpcnt = scnt;
                seqcnt += scnt;
                seqtot += stot;
                if (gtot > gapgrptot)
                  gapgrptot = gtot;
                if (gcnt > gapgrpcnt)
                  gapgrpcnt = gcnt;
                gapcnt += gcnt;
                gaptot += gtot;
              }

            outone->info['S']->given.count = seqcnt;
            outone->info['S']->given.max   = seqmax;
            outone->info['S']->given.total = seqtot;
            outone->info['n']->given.count = gapcnt;
            outone->info['n']->given.max   = gapmax;
            outone->info['n']->given.total = gaptot;
            outone->info['s']->given.count = scfcnt;
            outone->info['s']->given.max   = scfmax;
            outone->info['s']->given.total = scftot;
            for (s = outone->info['s']->stats; s->type != '\0'; s++)
              { if (s->type == 'S')
                  { s->maxCount = seqgrpcnt;
                    s->maxTotal = seqgrptot;
                  }
                if (s->type == 'n')
                  { s->maxCount = gapgrpcnt;
                    s->maxTotal = gapgrptot;
                  }
              }
          }

        for (i = 0; i < gdb->nscaff; i++)
          { header = gdb->headers + scf[i].hoff;

            oneInt(outone,0) = scf[i].slen;
            oneWriteLine(outone,'s',strlen(header),header);

            spos = 0;
            for (k = scf[i].fctg; k < scf[i].ectg; k++)
              { if (ctg[k].sbeg > spos) 
                  { oneChar(outone,0) = 'n';
                    oneInt(outone,1)  = ctg[k].sbeg-spos;
                    oneWriteLine(outone,'n',0,0);
                  }
                spos = ctg[k].sbeg;
    
                if (Get_Contig(gdb,k,UPPER,contig) == NULL)
                  goto clean_up2;
                oneWriteLine(outone,'S',ctg[k].clen,contig);
                spos += ctg[k].clen;
              }
            if (spos < scf[i].slen)
              { oneChar(outone,0) = 'n';
                oneInt(outone,1)  = scf[i].slen-spos;
                oneWriteLine(outone,'n',0,0);
              }
          }

        oneFileClose(outone);
        oneSchemaDestroy(schema);
      }

    else

      { char       *nstring;
        int         i;

        if (argc == 2)
          output = stdout;
        else
          { if (gzip)
              output = gzopen(Catenate(TPATH,"/",TROOT,suffix[TEXTN]),"w");
            else
              output = fopen(Catenate(TPATH,"/",TROOT,suffix[TEXTN]),"w");
            if (output == NULL)
              { fprintf(stderr,"%s: Cannot open %s/%s.%s for writing\n",
                               Prog_Name,TPATH,TROOT,suffix[TEXTN]);
                exit (1);
              }
          }
    
        nstring = Malloc(WIDTH+1,"Allocating write buffer\n");
        if (UPPER == UPPER_CASE)
          for (i = 0; i < WIDTH; i++)
            nstring[i] = 'N';
        else
          for (i = 0; i < WIDTH; i++)
            nstring[i] = 'n';
        nstring[WIDTH] = '\0';
    
        //   For the relevant range of reads, write each to the file
        //     recreating the original headers with the index meta-data about each read

        for (i = 0; i < gdb->nscaff; i++)
          { int   k, j, w, z;
            int   len;
            int64 spos, wpos;
    
            header = gdb->headers + scf[i].hoff;
            if (gzip)
              z = gzprintf(output,">%s\n",header);
            else
              z = fprintf(output,">%s\n",header);
            if (z < 0)
              goto out_error;

#define WIDTH_WRITE(target)				\
  for (j = 0; j+(w = WIDTH-wpos) <= len; j += w)	\
    { if (gzip)						\
        z = gzprintf(output,"%.*s\n",w,target);		\
      else						\
        z = fprintf(output,"%.*s\n",w,target);		\
      if (z < 0)					\
        goto out_error;					\
      wpos = 0;						\
    }							\
  if (j < len)						\
    { if (gzip)						\
        z = gzprintf(output,"%.*s",len-j,target);	\
      else						\
        z = fprintf(output,"%.*s",len-j,target);	\
      if (z < 0)					\
        goto out_error;					\
      if (j == 0)					\
        wpos += len;					\
      else						\
        wpos = len-j;					\
    }

            wpos = 0;
            spos = 0;
            for (k = scf[i].fctg; k < scf[i].ectg; k++)
              { if (ctg[k].sbeg > spos) 
                  { len = ctg[k].sbeg - spos;
                    WIDTH_WRITE(nstring)
                  }
                spos = ctg[k].sbeg;
    
                Get_Contig(gdb,k,UPPER,contig);
 
                len = ctg[k].clen;
                WIDTH_WRITE(contig+j)
                spos += len;
              }
            if (spos < scf[i].slen)
              { len = scf[i].slen - spos;
                WIDTH_WRITE(nstring)
              }
            if (wpos > 0)
              { if (gzip)
                  z = gzprintf(output,"\n");
                else
                  z = fprintf(output,"\n");
                if (z < 0)
                  goto out_error;
              }
          }
        if (output != stdout)
          { if (gzip)
              gzclose(output);
            else
              fclose(output);
          }
        free(nstring);
      }
  }

  Close_GDB(gdb);

  exit (0);

out_error:
  fprintf(stderr,"%s: Failed to write to fasta file\n",Prog_Name);
  fclose(output);
  unlink(Catenate(TPATH,"/",TROOT,suffix[TEXTN]));
  exit (1);

clean_up2:
  oneFileClose(outone);
  unlink(Catenate(TPATH,"/",TROOT,suffix[TEXTN]));
  exit (1);
}
