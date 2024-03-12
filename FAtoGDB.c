/*******************************************************************************************
 *
 *  Convert a .fasta or .1seq file to a .gdb:
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
#include <zlib.h>
#include <dirent.h>

#include "DB.h"
#include "alncode.h"

static char *Usage = "[-v] <source:path>[<fa_extn>|<1_extn>] [<target:path>[.gdb]]";

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

int main(int argc, char *argv[])
{ char   *SPATH, *SROOT, *SEXTN;
  char   *TPATH, *TROOT;
  int     is_one, gzipd;

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
        fprintf(stderr,"           <1_extn>  = any valid 1-code sequence file type\n");
        exit (1);
      }
  }

  //  Determine source and target root names, paths, and extensions

  { char   *suffix[8] = { "", ".fa", ".fna", ".fasta", ".gz", ".fa.gz", ".fna.gz", ".fasta.gz" };
    int     suflen[8] = { 0, 3, 4, 6, 3, 6, 7, 9 };
    char   *p, *s, *e;
    int     i, j, exists;
    FILE   *input;
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

    Schema = make_Aln_Schema();
    if (Schema == NULL)
      { fprintf(stderr,"%s: Failed to create onecode schema\n",Prog_Name);
        exit (1);
      }

    for (i = 7; i >= 0; i--)
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
      fclose(input);

    if (i >= 0)
      { SEXTN = suffix[i];
        if (i%4 == 0)
          { e = SROOT + strlen(SROOT);
            for (j = 1; j < 8; j++)
              if (strcmp(e-suflen[j],suffix[j]) == 0 && j != 4)
                break;
            i += j;
            if (i < 8)
              { e[-suflen[j]] = '\0';
                SEXTN = suffix[i];
              }
            else
              { OneFile *of;

                of = oneFileOpenRead(Catenate(SPATH,"/",SPATH,""),Schema,"seq",1);
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
      }
    is_one = (i <  0);
    gzipd  = (i >= 4);

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

  { DAZZ_DB   db;
    int64     maxlen, totlen, count[4];
    int64     offset, hdrset;
    int       nreads;
    DAZZ_READ prec;

    bases  = Fopen(Catenate(TPATH,"/.",TROOT,".bps"),"w");
    indx   = Fopen(Catenate(TPATH,"/.",TROOT,".idx"),"w");
    hdrs   = Fopen(Catenate(TPATH,"/.",TROOT,".hdr"),"w");
    if (bases == NULL || indx == NULL || hdrs == NULL)
      goto clean_up;

    if (fwrite(&db,sizeof(DAZZ_DB),1,indx) < 1)
      goto out_error;

    nreads = 0;
    totlen = 0;
    maxlen = 0;
    bzero(count,sizeof(int64)*4);
    offset = 0;
    hdrset = 0;

    if (is_one)
      { int   i;
        int   hlen, rlen, clen;
        bool  isScaffold;
        int64 scafPos, nContig; 
        int   byte_cnt[256];     // times byte occurs in compressed sequence

        OneFile *of = oneFileOpenRead (Catenate(SPATH,"/",SROOT,SEXTN), Schema, "seq", 1) ;
        if (of == NULL)
          goto clean_up ;

        // line sequence must be either (I!S)+, or (s(n!(Sn)*Sn!)+, e.g.
        //   - scaffold files with s (scaffold id) and n (inter-contig N's) and no I lines, or
        //   - straight sequence files with no s and n lines, and only S lines, allowing I lines

        // Ignore N characters in sequences (non-acgt chars within contigs)

        if (of->info['s']->given.count > 0 && of->info['I']->given.count == 0)
          isScaffold = true;
        else if (of->info['s']->given.count == 0 && of->info['n']->given.count == 0)
          isScaffold = false;
        else
          { fprintf(stderr,"%s: 1seq file has incomplete scaffold properties - can't use it\n",
                           Prog_Name);
            fprintf(stderr,"          I.E. #'s' = %lld  #'I' = %lld  #'n' = %lld\n",
                     of->info['s']->given.count, of->info['I']->given.count,
                     of->info['n']->given.count) ;
            goto one_error;
          }

        for (i = 0; i < 256; i++)
          byte_cnt[i] = 0;
        scafPos   = 0;
        nContig   = 0;
        prec.coff = 0;

        while (oneReadLine (of))
          switch (of->lineType)

          { case 's': // scaffold name
              hlen = oneLen(of);
              if (hlen > MAX_NAME-2)
                { fprintf(stderr,"%s: Line %lld: Scaffold name",Prog_Name,of->line);
                  fprintf(stderr," in file %s%s is too long (> %d chars)\n",
                                 SROOT,SEXTN,MAX_NAME-2);
                  goto one_error;
                }
              prec.coff = hdrset;
              if (fwrite (oneString(of), hlen, 1, hdrs) != 1)
                goto out_error;
              if (putc('\n',hdrs) != '\n')
                goto out_error;
              hdrset += hlen + 1;

              scafPos = 0;
              nContig = 0;
              break ;

            case 'n': // non-acgt chars outside sequences
              scafPos += oneInt(of,1);
              break ;

            case 'S': // sequence object
              rlen = oneLen(of);

              prec.origin = nContig;
              prec.rlen   = rlen;
              prec.fpulse = scafPos;
              prec.boff   = offset;
              prec.flags  = DB_BEST;
              if (fwrite (&prec, sizeof(DAZZ_READ),1,indx) != 1)
                goto clean_up ;

              clen = COMPRESSED_LEN(rlen);
              uint8* bytes = oneDNA2bit(of);
              if (fwrite (bytes, clen, 1, bases) != 1)
                goto out_error;
              offset += clen;

              nreads += 1;
              totlen += rlen;
              if (rlen > maxlen)
                maxlen = rlen;
              for (i = 0 ; i < clen; i++)
                byte_cnt[bytes[i]] += 1;
              if (rlen & 0x3)
                count[0] -= (4-(rlen&0x3));

              if (isScaffold)
                { scafPos += rlen;
                  nContig += 1;
                }
              break ;

            case 'I': // identifier of sequence
              hlen = oneLen(of);
              if (hlen > MAX_NAME-2)
                { fprintf(stderr,"%s: Line %lld: Scaffold name",Prog_Name,of->line);
                  fprintf(stderr," in file %s%s is too long (> %d chars)\n",
                                 SROOT,SEXTN,MAX_NAME-2);
                  goto one_error;
                }
              prec.coff = hdrset;
              if (fwrite (oneString(of), hlen, 1, hdrs) != 1)
                goto out_error;
              if (putc ('\n', hdrs) != '\n')
                goto out_error;
              hdrset += hlen + 1;
              break;
            }

        oneFileClose(of);

        for (i = 0; i < 256; i++)
          { count[i & 0x3] += byte_cnt[i];
            count[(i>>2) & 0x3] += byte_cnt[i];
            count[(i>>4) & 0x3] += byte_cnt[i];
            count[(i>>6) & 0x3] += byte_cnt[i];
          }
      }

    else // a fasta file of some sort

      { int64      rmax, rlen;
        char      *read;
        int        nline, eof;
        int        i, x, n;
        void      *input;

        //  Buffer for accumulating .fasta sequence over multiple lines

        rmax  = MAX_NAME + 10000000;
        read  = (char *) Malloc(rmax+1,"Allocating line buffer");
        if (read == NULL)
          goto clean_up;

        if (gzipd)
          input = gzopen(Catenate(SPATH,"/",SROOT,SEXTN),"r");
        else
          input = fopen(Catenate(SPATH,"/",SROOT,SEXTN),"r");

        //  Get the header of the first line.  If the file is empty skip.

        rlen  = 0;
        nline = 1;
        if (gzipd)
          eof = (gzgets(input,read,MAX_NAME) == NULL);
        else
          eof = (fgets(read,MAX_NAME,input) == NULL);
        if (eof)
          { if (gzipd)
              { if (!gzeof(input))
                  goto fna_error;
              }
            else
              { if (!feof(input))
                  goto fna_error;
              }
          }
        if (eof || strlen(read) < 1)
          { fclose(input);
            fprintf(stderr,"Input is empty, terminating!\n");
            goto clean_up;
          }

        // Check that the first line is a header line

        if (read[0] != '>')
          { fprintf(stderr,"%s: Line 1: First header in fasta file %s%s is missing\n",
                           Prog_Name,SROOT,SEXTN);
            goto clean_up;
          }
        if (read[strlen(read)-1] != '\n')
          { fprintf(stderr,"%s: Line 1: Fasta header line in %s%s is too long (> %d chars)\n",
                           Prog_Name,SROOT,SEXTN,MAX_NAME-2);
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
              { if (gzipd)
                  eof = (gzgets(input,read+rlen,MAX_NAME) == NULL);
                else
                  eof = (fgets(read+rlen,MAX_NAME,input) == NULL);
                if (eof)
                  { if (gzipd)
                      { if (!gzeof(input))
                          goto fna_error;
                      }
                    else
                      { if (!feof(input))
                          goto fna_error;
                      }
                  }
                x = strlen(read+rlen)-1;
                if (read[rlen] == '>')
                  { if (read[rlen+x] != '\n')
                      { fprintf(stderr,"%s: Line %d: Fasta header line",Prog_Name,nline);
                        fprintf(stderr," in file %s%s is too long (> %d chars)\n",
                                       SROOT,SEXTN,MAX_NAME-2);
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
                                       SROOT,SEXTN);
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

        if (gzipd)
          gzclose(input);
        else
          fclose(input);
      }

    oneSchemaDestroy(Schema);

    //  Update relevant fields in db record

    { int c;

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
    }

    fclose(indx);
    fclose(bases);
    fclose(hdrs);

    { char *cpath;
      FILE *ostub;

      ostub = Fopen(Catenate(TPATH,"/",TROOT,".gdb"),"w");
      if (ostub == NULL)
        goto clean_up;
      if (fprintf(ostub,DB_NFILE,1) < 0)
        goto out_error;
      if (SPATH[0] == '/')
        { if (is_one)
            { if (fprintf(ostub,DB_FDATA,db.ureads,Catenate(SPATH,"/",SROOT,".1seq"),"") < 0)
                goto out_error;
            }
          else
            { if (fprintf(ostub,DB_FDATA,db.ureads,Catenate(SPATH,"/",SROOT,SEXTN),"") < 0)
                goto out_error;
            }
        }
      else
        { cpath = getcwd(NULL,0);
          if (is_one)
            { if (fprintf(ostub,DB_FDATA,db.ureads,Catenate(SPATH,"/",SROOT,".1seq"),cpath) < 0)
                goto out_error;
            }
          else
            { if (fprintf(ostub,DB_FDATA,db.ureads,Catenate(SPATH,"/",SROOT,SEXTN),cpath) < 0)
                goto out_error;
            }
          free(cpath);
        }
      if (fprintf(ostub,DB_NBLOCK,1) < 0) goto out_error;
      if (fprintf(ostub,DB_PARAMS,db.totlen,0,1) < 0) goto out_error;
      if (fprintf(ostub," %9d %9d\n",0,0) < 0) goto out_error;
      if (fprintf(ostub," %9d %9d\n",db.ureads,db.ureads) < 0) goto out_error;
      fclose(ostub);
    }
  }

  exit (0);

  //  Error exit:  Either truncate or remove the .idx, .bps, and .hdr files as appropriate.
  //               Remove the new image file <pwd>/<root>.dbx

fna_error:
  fprintf(stderr,"%s: IO error reading fasta file\n",Prog_Name);
  goto clean_up;

one_error:
  fprintf(stderr,"%s: IO error reading 1seq file\n",Prog_Name);
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
