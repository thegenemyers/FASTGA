/*******************************************************************************************
 *
 *  Convert a .fasta files to a .gdb:
 *
 *  Author:  Gene Myers
 *  Origin:  Complete rework/specialization from fasta2DAM of DAZZ_DB for FASTGA module
 *  Date  :  Jan 2024
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include "DB.h"

static char *Usage = "[-v] <source:path>[<fa_extn>] [<target:path>[.gdb]]";

static char number[128] =
    { 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 4, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 4, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
    };


static FILE *fzopen(const char *path, const char *mode)
{ gzFile zfp;

  zfp = gzopen(path,mode);
  if (zfp == NULL)
    return (fopen(path,mode));

  return (funopen(zfp,
                  (int(*)(void*,char*,int))gzread,
                  (int(*)(void*,const char*,int))gzwrite,
                  (fpos_t(*)(void*,fpos_t,int))gzseek,
                  (int(*)(void*))gzclose) );
}

int main(int argc, char *argv[])
{ char   *SPATH, *SROOT;
  char   *TPATH, *TROOT;
  int     SEXTN;
  char   *suffix[8] = { "", ".fa", ".fna", ".fasta", ".gz", ".fa.gz", ".fna.gz", ".fasta.gz" };
  int     suflen[8] = { 0, 3, 4, 6, 3, 6, 7, 9 };
  FILE   *input;
  DAZZ_DB db;
  FILE   *bases, *indx, *hdrs;

  int VERBOSE;

  //   Process command line

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("FAtoGDB")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc < 2 || argc > 3)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"           <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        exit (1);
      }
  }

  //  Determine source and target root names, paths, and extensions
  
  { char *p, *s, *e;
    int   i, exists;
    struct stat status;

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
    for (i = 7; i >= 0; i--)
      { input = fzopen(Catenate(SPATH,"/",SROOT,suffix[i]),"r");
        if (input != NULL)
          break;
      }
    if (i < 0)
      {fprintf(stderr,"%s: Could not find a fasta file with base name %s\n",Prog_Name,SROOT);
        exit (1);
      }
    SEXTN = i;
    if (i%4 == 0)
      { e = SROOT + strlen(SROOT);
        for (i = 1; i < 8; i++)
          if (strcmp(e-suflen[i],suffix[i]) == 0 && i != 4)
            break;
        SEXTN += i;
        if (SEXTN >= 8)
          { fprintf(stderr,"%s: Could not find valid fasta extension of %s\n",Prog_Name,SROOT);
            exit (1);
          }
        e[-suflen[i]] = '\0';
      }

    if (argc == 2)
      { TPATH = SPATH;
        TROOT = SROOT;
      }
    else // argc == 3
      { TPATH = argv[2];
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
                  { fprintf(stderr,"%s: Exisiting target %s must be a .gdb\n",Prog_Name,TROOT);
                    exit (1);
                  }
              }
            else if (exists)
              { fprintf(stderr,"%s: Exisiting target %s must be a .gdb\n",Prog_Name,TROOT);
                exit (1);
              }
          }
        else
          TROOT = SROOT;
      }

    if (VERBOSE)
      { if (strcmp(TPATH,".") == 0)
          fprintf(stderr,"\n  Creating genome data base (GDB) %s.gdb in the current directory\n",
                         TROOT);  
        else
          fprintf(stderr,"\n  Creating genome data base (GDB) %s.gdb in directory %s\n",
                         TROOT,TPATH);  
        fflush(stderr);
      }
  }

  { int64      offset, hdrset;
    DAZZ_READ  prec;
    int64      maxlen, totlen, count[4];
    int        nreads;
    int64      rmax, rlen;
    char      *read;
    int        nline, eof;
    int        c, i, x, n;

    nreads  = 0;
    offset  = 0;
    hdrset  = 0;

    //  Buffer for accumulating .fasta sequence over multiple lines

    rmax  = MAX_NAME + 10000000;
    read  = (char *) Malloc(rmax+1,"Allocating line buffer");
    if (read == NULL)
      exit (1);

    totlen = 0;              //  total # of bases in new .fasta files
    maxlen = 0;              //  longest read in new .fasta files
    for (c = 0; c < 4; c++)  //  count of acgt in new .fasta files
      count[c] = 0;

    bases  = Fopen(Catenate(TPATH,"/.",TROOT,".bps"),"w");
    indx   = Fopen(Catenate(TPATH,"/.",TROOT,".idx"),"w");
    hdrs   = Fopen(Catenate(TPATH,"/.",TROOT,".hdr"),"w");
    if (bases == NULL || indx == NULL || hdrs == NULL)
      goto clean_up;

    if (fwrite(&db,sizeof(DAZZ_DB),1,indx) < 1)
      goto out_error;

    //  Get the header of the first line.  If the file is empty skip.

    rlen  = 0;
    nline = 1;
    eof   = (fgets(read,MAX_NAME,input) == NULL);
    if (eof & !feof(input))
      goto in_error;
    if (eof || strlen(read) < 1)
      { fclose(input);
        fprintf(stderr,"Input is empty, terminating!\n");
        goto clean_up;
      }

    // Check that the first line is a header line

    if (read[0] != '>')
      { fprintf(stderr,"%s: Line 1: First header in fasta file %s%s is missing\n",
                       Prog_Name,SROOT,suffix[SEXTN]);
        goto clean_up;
      }
    if (read[strlen(read)-1] != '\n')
      { fprintf(stderr,"%s: Line 1: Fasta header line in %s%s is too long (> %d chars)\n",
                       Prog_Name,SROOT,suffix[SEXTN],MAX_NAME-2);
        goto clean_up;
      }

    //  Read in all the sequences until end-of-file

    while (!eof)
      { int hlen;

        for (rlen++; read[rlen] != '\n'; rlen++)
          if (!isspace(read[rlen]))
            break;
        hlen = strlen(read+rlen);
        if ((int) fwrite(read+rlen,1,hlen,hdrs) < hlen)
          goto out_error;

        rlen = 0;
        while (1)
          { eof = (fgets(read+rlen,MAX_NAME,input) == NULL);
            if (eof & !feof(input))
              goto in_error;
            x = strlen(read+rlen)-1;
            if (read[rlen] == '>')
              { if (read[rlen+x] != '\n')
                  { fprintf(stderr,"%s: Line %d: Fasta header line",Prog_Name,nline);
                    fprintf(stderr," in file %s%s is too long (> %d chars)\n",
                                   SROOT,suffix[SEXTN],MAX_NAME-2);
                    goto clean_up;
                  }
                nline += 1;
                break;
              }
            if (eof)
              break;
            if (read[rlen+x] == '\n')
              nline += 1;
            else
              x += 1;
            rlen += x;
            if (rlen + MAX_NAME > rmax)
              { rmax = ((int64) (1.4 * rmax)) + 10000000 + MAX_NAME;
                read = (char *) realloc(read,rmax+1);
                if (read == NULL)
                  { fprintf(stderr,"%s: Line %d:",Prog_Name,nline);
                    fprintf(stderr," Out of memory allocating line buffer while reading %s%s\n",
                                   SROOT,suffix[SEXTN]);
                    goto clean_up;
                  }
              }
          }
        read[rlen] = '\0';

        n = 0;
        i = -1;
        while (i < rlen)
          { int pbeg, plen, clen;

            while (i < rlen)
              if (number[(int) read[++i]] < 4)
                break;

            // if (i >= rlen) break;

            pbeg = i;
            prec.fpulse = pbeg;
            prec.origin = n++;
            prec.boff   = offset;
            prec.coff   = hdrset;
            prec.flags  = DB_BEST;
            while (i < rlen)
              { x = number[(int) read[i]];
                if (x >= 4) break;
                count[x] += 1;
                read[i++] = (char) x;
              }
            prec.rlen = plen = i-pbeg;
            nreads += 1;
            totlen += plen;
            if (plen > maxlen)
              maxlen = plen;

            Compress_Read(plen,read+pbeg);
            clen = COMPRESSED_LEN(plen);
            if ((int) fwrite(read+pbeg,1,clen,bases) < clen)
              goto out_error;
            offset += clen;

            if (fwrite(&prec,sizeof(DAZZ_READ),1,indx) < 1)
              goto out_error;
          }
        hdrset += hlen;
      }

    fclose(input);

    //  Update relevant fields in db record

    db.ureads = nreads;
    for (c = 0; c < 4; c++)
      db.freq[c] = (float) ((1.*count[c])/totlen);
    db.totlen = totlen;
    db.maxlen = maxlen;
    db.treads = nreads;     //  Start with split into a single huge block
    db.cutoff = 0;
    db.allarr = 1;

    if (fseek(indx,0,SEEK_SET) < 0)
      goto out_error;
    if (fwrite(&db,sizeof(DAZZ_DB),1,indx) < 1)   //  Write the finalized db record into .idx
      goto out_error;

    fclose(indx);
    fclose(bases);
    fclose(hdrs);
  }

  { char *cpath;
    FILE *ostub;

    ostub = Fopen(Catenate(TPATH,"/",TROOT,".gdb"),"w");
    if (ostub == NULL)
      goto clean_up;
    if (fprintf(ostub,DB_NFILE,1) < 0) goto out_error;
    if (SPATH[0] == '/')
      { if (fprintf(ostub,DB_FDATA,db.ureads,Catenate(SPATH,"/",SROOT,suffix[SEXTN]),"") < 0)
          goto out_error;
      }
    else
      { cpath = getcwd(NULL,0);
        if (fprintf(ostub,DB_FDATA,db.ureads,Catenate(SPATH,"/",SROOT,suffix[SEXTN]),cpath) < 0)
          goto out_error;
        free(cpath);
      }
    if (fprintf(ostub,DB_NBLOCK,1) < 0) goto out_error;
    if (fprintf(ostub,DB_PARAMS,db.totlen,0,1) < 0) goto out_error;
    if (fprintf(ostub," %9d %9d\n",0,0) < 0) goto out_error;
    if (fprintf(ostub," %9d %9d\n",db.ureads,db.ureads) < 0) goto out_error;
    fclose(ostub);
  }
  
  exit (0);

  //  Error exit:  Either truncate or remove the .idx, .bps, and .hdr files as appropriate.
  //               Remove the new image file <pwd>/<root>.dbx

in_error:
  fprintf(stderr,"%s: IO error reading fasta file\n",Prog_Name);
  goto clean_up;

out_error:
  fprintf(stderr,"%s: IO error writing gdb files\n",Prog_Name);

clean_up:
  if (indx != NULL)
    { fclose(indx);
      unlink(Catenate(TPATH,"/.",TROOT,".idx"));
    }
  if (bases != NULL)
    { fclose(bases);
      unlink(Catenate(TPATH,"/.",TROOT,".bps"));
    }
  if (hdrs != NULL)
    { fclose(hdrs);
      unlink(Catenate(TPATH,"/.",TROOT,".hdr"));
    }

  exit (1);
}
