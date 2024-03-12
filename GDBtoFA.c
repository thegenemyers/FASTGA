/********************************************************************************************
 *
 *  Recreate the .fasta file that gave rise to a given .gdb
 *
 *  Author:  Gene Myers
 *  Origin:  Reworking of DAZZ_DB program DAM2fasta.c for the FASTGA module
 *  Date  :  Jan 2024
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <zlib.h>

#include "DB.h"
#include "alncode.h"

static char *Usage =
         "[-vU] [-w<int(80)>] <source:path>[.gdb] [ @ | <target:path>[<fa_extn>|<1_extn>] ]";

int main(int argc, char *argv[])
{ DAZZ_DB    _db, *db = &_db;
  char       *SPATH, *SROOT;
  char       *TPATH, *TROOT;
  int         TEXTN;
  int         is_one;
  OneSchema  *schema;
  OneFile    *outone;
  int         gzip;
  void       *output;
  char       *command;

  int UPPER;   // -U
  int WIDTH;   // -w
  int VERBOSE; // -v

  char *suffix[7] = { ".1seq", ".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz" };
  int   suflen[7] = { 5, 3, 4, 6, 6, 7, 9 };

  { int   n, i;
    char *c;

    n = 0;
    for (i = 1; i < argc; i++)
      n += strlen(argv[i])+1;

    command = Malloc(n+1,"Allocating command string");
    if (command == NULL)
      exit (1);

    c = command;
    if (argc >= 1)
      { c += sprintf(c,"%s",argv[1]);
        for (i = 2; i < argc; i++)
          c += sprintf(c," %s",argv[i]);
      }
    *c = '\0';
  }

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

    UPPER   = 1 + flags['U'];
    VERBOSE = flags['v'];

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -U: Use upper case for DNA (default is lower case).\n");
        fprintf(stderr,"      -w: Print -w bp per line (default is 80).\n");
        exit (1);
      }
  }

  //  Open db & determine source and target parts

  { int         status, exists;
    DAZZ_STUB  *stub;
    struct stat tdesc;
    char       *p, *e;
    char       *oname, *opath;
    int         i;

    status = Open_DB(argv[1],db);
    if (status < 0)
      exit (1);
    if (status == 0)
      { fprintf(stderr,"%s: Cannot open %s as a .gdb\n",Prog_Name,argv[1]);
        exit (1);
      }

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
    if (strcmp(SROOT+(strlen(SROOT)-4),".gdb") == 0)
      SROOT[strlen(SROOT)-4] = '\0';

    stub   = Read_DB_Stub(Catenate(SPATH,"/",SROOT,".gdb"),DB_STUB_FILES|DB_STUB_PROLOGS);
    oname  = stub->fname[0];
    opath  = stub->prolog[0];

    if (argc == 2)
      { if (strcmp(oname+(strlen(oname)-5),".1seq") == 0)
          { is_one = 1;
            schema = oneSchemaCreateFromText(alnSchemaText);
            if (schema == NULL)
              { fprintf(stderr,"%s: Failed to create ONEcode schema\n",Prog_Name);
                exit (1);
              }
            outone = oneFileOpenWriteNew("-",schema,"seq",0,0);
            if (outone == NULL)
              { fprintf(stderr,"%s: Cannot open ONEcode for writing\n",Prog_Name);
                exit (1);
              }
          }
        else
          { is_one = 0;
            gzip   = 0;
            output = stdout;
          }
        TPATH  = NULL;
        TROOT  = NULL;
        TEXTN  = 0;
      }
    else
      { if (strcmp(argv[2],"@") == 0)
          { TPATH = oname;
            if (opath[0] != '\0')
              { if (TPATH[0] == '.' && TPATH[1] == '/')
                  TPATH = Strdup(Catenate(opath,"/",TPATH+2,""),"");
                else
                  TPATH = Strdup(Catenate(opath,"/",TPATH,""),"");
              }
            p = rindex(TPATH,'/');
            *p++ = '\0';
            TROOT = p; 
          }
        else
          { TPATH = argv[2];
            exists = (stat(TPATH,&tdesc) >= 0);
            if ( !exists || (tdesc.st_mode & S_IFMT) != S_IFDIR)
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
            else
              TROOT = SROOT;
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
          { e = oname + strlen(oname);
            for (i = 6; i >= 0; i--)
              if (strcmp(e-suflen[i],suffix[i]) == 0)
                break;
            TEXTN = i;
          }
        is_one = (TEXTN == 0);
        gzip   = (TEXTN >= 3);
        if (is_one)
          { schema = oneSchemaCreateFromText(alnSchemaText);
            if (schema == NULL)
              { fprintf(stderr,"%s: Failed to create ONEcode schema\n",Prog_Name);
                exit (1);
              }
            outone = oneFileOpenWriteNew(Catenate(TPATH,"/",TROOT,".1seq"),schema,"seq",1,0);
            if (outone == NULL)
              { fprintf(stderr,"%s: Cannot open %s/%s.%s for writing\n",
                               Prog_Name,TPATH,TROOT,suffix[TEXTN]);
                exit (1);
              }
          }
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
      }

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
              fprintf(stderr,"\n  Creating %sfasta file %s%s at current directory\n",
                             TEXTN>=3?"gzip'd ":" ",TROOT,suffix[TEXTN]);
            else
              fprintf(stderr,"\n  Creating %sfasta file %s%s at directory %s\n",
                             TEXTN>=3?"gzip'd ":" ",TROOT,suffix[TEXTN],TPATH);
	  }
        fflush(stderr);
      }
  }

  //  For each file do:

  { DAZZ_READ  *reads;
    char       *read;
    int         nreads;
    FILE       *hdrs;
    char        header[MAX_NAME];

    hdrs = Fopen(Catenate(db->path,".hdr","",""),"r");
    if (hdrs == NULL)
      exit (1);

    reads = db->reads;
    read  = New_Read_Buffer(db);
    if (read == NULL)
      exit (1);

    nreads = db->nreads;

    if (is_one)

      { int64 last;
        int   r, s;

        oneAddProvenance(outone,Prog_Name,"0.1",command);

        if (!outone->isBinary)
          { int64 seqcnt, seqmax, seqtot;
            int64 gapcnt, gapmax, gaptot;
            int64 scfcnt, scfmax, scftot;
            int64 seqgrpcnt, gapgrpcnt, scnt, gcnt;
            int64 seqgrptot, gapgrptot, stot, gtot;

            seqcnt = seqmax = seqtot = 0;
            gapcnt = gapmax = gaptot = 0;
            scfcnt = scfmax = scftot = 0;
            seqgrpcnt = gapgrpcnt = 0;
            seqgrptot = gapgrptot = 0;

            r = 0;
            while (r < nreads)
              { for (s = r+1; s < nreads; s++)
                  if (reads[s].origin == 0)
                    break;
    
                if (fseeko(hdrs,reads[r].coff,SEEK_SET) < 0)
                  { fprintf(stderr,"%s: Cannot seek GDB header file\n",Prog_Name);
                    goto clean_up;
                  }
                if (fgets(header,MAX_NAME,hdrs) == NULL)
                  { fprintf(stderr,"%s: Failed to read from the GDB header file\n",Prog_Name);
                    goto clean_up;
                  }
                if ((int64) strlen(header) >= scfmax)
                  scfmax = strlen(header)-1;
                scfcnt += 1;
                scftot += strlen(header)-1;

                gtot = gcnt = 0;
                stot = scnt = 0;
                last = 0;
                while (r < s)
                  { if (reads[r].fpulse - last > 0) 
                      { if (reads[r].fpulse-last > gapmax)
                          gapmax = reads[r].fpulse-last;
                        gtot += reads[r].fpulse-last;
                        gcnt += 1;
                      }
                    last = reads[r].fpulse + reads[r].rlen;
                    if (reads[r].rlen > 0)
                      { if (reads[r].rlen > seqmax)
                          seqmax = reads[r].rlen;
                        stot += reads[r].rlen;
                        scnt += 1;
                      }
                    r += 1;
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
            outone->info['S']->given.groupCount = seqgrpcnt;
            outone->info['S']->given.groupTotal = seqgrptot;
            outone->info['N']->given.count = gapcnt;
            outone->info['N']->given.max   = gapmax;
            outone->info['N']->given.total = gaptot;
            outone->info['N']->given.groupCount = gapgrpcnt;
            outone->info['N']->given.groupTotal = gapgrptot;
            outone->info['s']->given.count = scfcnt;
            outone->info['s']->given.max   = scfmax;
            outone->info['s']->given.total = scftot;
          }

        r = 0;
        while (r < nreads)
          { for (s = r+1; s < nreads; s++)
              if (reads[s].origin == 0)
                break;

            if (fseeko(hdrs,reads[r].coff,SEEK_SET) < 0)
              { fprintf(stderr,"%s: Cannot seek GDB header file\n",Prog_Name);
                goto clean_up;
              }
            if (fgets(header,MAX_NAME,hdrs) == NULL)
              { fprintf(stderr,"%s: Failed to read from the GDB header file\n",Prog_Name);
                goto clean_up;
              }
            header[strlen(header)-1] = '\0';

            if (outone->isBinary)
              oneInt(outone,0) = 0;
            else
              oneInt(outone,0) = s-r;
            oneInt(outone,1) = reads[s-1].fpulse + reads[s-1].rlen;
            oneWriteLine(outone,'s',strlen(header),header);

            last = 0;
            while (r < s)
              { if (reads[r].fpulse - last > 0) 
                  { oneChar(outone,0) = 'n';
                    oneInt(outone,1)  = reads[r].fpulse-last; 
                    oneWriteLine(outone,'n',0,0);
                  }
                last = reads[r].fpulse + reads[r].rlen;
                if (reads[r].rlen > 0)
                  { if (Load_Read(db,r,read,1))
                      goto clean_up;
                    oneWriteLine(outone,'S',reads[r].rlen,read);
                  }
                r += 1;
              }
          }

        oneFileClose(outone);
        oneSchemaDestroy(schema);
      }

    else

      { char       *nstring;
        int         i, f, z, wpos;
    
        nstring = Malloc(WIDTH+1,"Allocating write buffer\n");
        if (nstring == NULL)
          exit (1);
    
        if (UPPER == 2)
          for (f = 0; f < WIDTH; f++)
            nstring[f] = 'N';
        else
          for (f = 0; f < WIDTH; f++)
            nstring[f] = 'n';
        nstring[WIDTH] = '\0';
    
        //   For the relevant range of reads, write each to the file
        //     recreating the original headers with the index meta-data about each read
    
        wpos = 0;
        for (i = 0; i < nreads; i++)
          { int        j, len, nlen, w;
            DAZZ_READ *r;
    
            r     = reads + i;
            len   = r->rlen;
    
            if (r->origin == 0)
              { if (i != 0 && wpos != 0)
                  { if (gzip)
                      z = gzprintf(output,"\n");
                    else
                      z = fprintf(output,"\n");
                    if (z < 0)
                      goto out_error;
                    wpos = 0;
                  }
                if (fseeko(hdrs,r->coff,SEEK_SET) < 0)
                  { fprintf(stderr,"%s: Cannot seek GDB header file\n",Prog_Name);
                    goto clean_up;
                  }
                if (fgets(header,MAX_NAME,hdrs) == NULL)
                  { fprintf(stderr,"%s: Failed to read from the GDB header file\n",Prog_Name);
                    goto clean_up;
                  }
                if (gzip)
                  z = gzprintf(output,"> %s",header);
                else
                  z = fprintf(output,"> %s",header);
                if (z < 0)
                  goto out_error;
              }
    
            if (r->fpulse != 0)
              { if (r->origin != 0)
                  nlen = r->fpulse - (reads[i-1].fpulse + reads[i-1].rlen);
                else
                  nlen = r->fpulse;
    
                for (j = 0; j+(w = WIDTH-wpos) <= nlen; j += w)
                  { if (gzip)
                      z = gzprintf(output,"%.*s\n",w,nstring);
                    else
                      z = fprintf(output,"%.*s\n",w,nstring);
                    if (z < 0)
                      goto out_error;
                    wpos = 0;
                  }
                if (j < nlen)
                  { if (gzip)
                      z = gzprintf(output,"%.*s",nlen-j,nstring);
                    else
                      z = fprintf(output,"%.*s",nlen-j,nstring);
                    if (z < 0)
                      goto out_error;
                    if (j == 0)
                      wpos += nlen;
                    else
                      wpos = nlen-j;
                  }
              }
    
            if (Load_Read(db,i,read,UPPER))
              exit (1);
    
            for (j = 0; j+(w = WIDTH-wpos) <= len; j += w)
              { if (gzip)
                  z = gzprintf(output,"%.*s\n",w,read+j);
                else
                  z = fprintf(output,"%.*s\n",w,read+j);
                if (z < 0)
                  goto out_error;
                wpos = 0;
              }
            if (j < len)
              { if (gzip)
                  z = gzprintf(output,"%s",read+j);
                else
                  z = fprintf(output,"%s",read+j);
                if (z < 0)
                  goto out_error;
                if (j == 0)
                  wpos += len;
                else
                  wpos = len-j;
              }
          }
        if (wpos > 0)
          { if (gzip)
              z = gzprintf(output,"\n");
            else
              z = fprintf(output,"\n");
            if (z < 0)
              goto out_error;
          }
        if (output != stdout)
          { if (gzip)
              gzclose(output);
            else
              fclose(output);
          }
        free(nstring);
      }

    fclose(hdrs);
  }

  Close_DB(db);

  exit (0);

out_error:
  fprintf(stderr,"%s: Failed to write to fasta file\n",Prog_Name);

clean_up:
  if (output != stdout)
    { if (is_one)
        fclose(output);
      else
        oneFileClose(outone);
      unlink(Catenate(TPATH,"/",TROOT,suffix[TEXTN]));
    }

  exit (1);
}
