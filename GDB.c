/*******************************************************************************************
 *
 *  Genome data base module.  Auxiliary routines to open and manipulate a data base for
 *    which the sequence and genome information are separated into two separate files, and the
 *    sequence is compressed into 2-bits for each base.  Derived from the daligner codes.
 *
 *  Author :  Gene Myers
 *  Date   :  May 2024
 *
 ********************************************************************************************/
/*  Last edited: Jul 25 19:12 2024 (rd109) */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <unistd.h>
#include <dirent.h>
#include <limits.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <zlib.h>

#include "gene_core.h"
#include "GDB.h"

/*******************************************************************************************
 *
 *  GENERAL UTILITIES
 *
 ********************************************************************************************/

static char *gdbSchemaText =
  "1 3 def 1 0                 schema for genome skeleton\n"
  ".\n"
  "P 3 gdb                     GDB\n"
  "D f 4 4 REAL 4 REAL 4 REAL 4 REAL   global: base frequency vector\n"
  "O S 1 6 STRING              id for a scaffold\n"
  "D G 1 3 INT                 gap of given length\n"
  "D C 1 3 INT                 contig of given length\n"
;

static OneSchema *make_GDB_Schema()
{ return (oneSchemaCreateFromText(gdbSchemaText)); }

static char *seqSchemaText =
  "1 3 def 1 0                 schema for aln and FastGA\n"
  ".\n"
  "P 3 seq                     SEQUENCE\n"
  "O s 2 3 INT 6 STRING        length and id for group of sequences = a scaffold\n"
  "G S                         scaffolds (s) group sequence objects (S)\n"
  "D n 2 4 CHAR 3 INT          non-acgt chars outside (between) sequences within scaffold\n"
  "O S 1 3 DNA                 sequence\n"
  "D I 1 6 STRING              identifier of sequence\n"
;

OneSchema *make_Seq_Schema()
{ return (oneSchemaCreateFromText(seqSchemaText)); }

static char *ERROR = "error";   //  Error return for subroutines returning char *

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
        { EPRINTF(EPLACE,"%s: Out of memory (Cating 4 strings)\n",Prog_Name);
          return (NULL);
        }
    }
  sprintf(cat,"%s%s%s%s",path,sep,root,suffix);
  return (cat);
}

/*******************************************************************************************
 *
 *  FILE NAME COMPLETION ROUTINE
 *
 ********************************************************************************************/

static char *find_1seq(char *path, char *root)
{ DIR           *dirp;
  struct dirent *dp;
  FILE          *f;
  int            rlen, one, len;
  char          *n, *extn, *name, type[4];

  dirp = opendir(path);
  if (dirp == NULL)
    { EPRINTF(EPLACE,"%s: Cannot open directory %s\n",Prog_Name,path);
      return (ERROR);
    }

  extn = NULL;
  rlen = strlen(root);
  while ((dp = readdir(dirp)) != NULL)
    { name = dp->d_name;
      if (strncmp(root,name,rlen) != 0)
        continue;
      if (name[rlen] != '.')
        continue;
      n = MyCatenate(path,"/",name,"");
      if (n == NULL)
        return (ERROR);
      f = fopen(n,"r");
      if (f == NULL)
        continue;
      if (fscanf(f," %d %d",&one,&len) != 2)
        goto escape;
      if (one != 1 || len != 3)
        goto escape;
      if (fscanf(f," %s",type) != 1)
        goto escape;
      if (strcmp(type,"seq") != 0)
        goto escape;
      fclose(f);
      if (extn != NULL)
        { EPRINTF(EPLACE,"%s: Two 1-code sequence files with root %s found\n",Prog_Name,root);
          return (ERROR);
        }
      extn = Strdup(name+rlen,"Allocating extension");
      if (extn == NULL)
        return (ERROR);
      continue;
    escape:
      fclose(f);
    }

  return (extn);
}

//  Interpret source & target arguments returning complete path + extension names of
//    the source and target in spath & tpath.  Returns the type of source.
//  If target == NULL then tpath has the same path & root as spath
//  Else If target is a directory then tpath has the same path as target
//  Else tpath has the same path and root as target.
//  If no_gdb then the source cannot be a gdb.
//  In interactive mode returns a negative number if an error occurs.
//  The user is responsible for freeing both tpath and spath.

int Get_GDB_Paths(char *source, char *target, char **spath, char **tpath, int no_gdb)
{ char   *SROOT, *SPATH;
  char   *TROOT, *TPATH;
  char   *SEXTN, *TEXTN;
  char   *suffix[10] = { "", ".fa", ".fna", ".fasta",
                         ".gz", ".fa.gz", ".fna.gz", ".fasta.gz", ".gdb", ".1gdb" };
  int     suflen[10] = { 0, 3, 4, 6, 3, 6, 7, 9, 4, 5 };
  char   *sdir, *tdir, *sdot, *tdot;

  OneSchema *schema;
  char *p, *e;
  int   i, tmax;
  FILE *input;
  struct stat status;

  sdir = NULL;
  sdot = NULL;
  tdir = NULL;
  tdot = NULL;

  if (no_gdb)
    tmax = 7;
  else
    tmax = 9;

  SPATH = source;
  p = rindex(SPATH,'/');
  if (p == NULL)
    { SROOT = SPATH;
      SPATH = ".";
    }
  else
    { *p = '\0';
      sdir = p;
      SROOT = p+1;
    }

  schema = make_Seq_Schema();
  if (schema == NULL)
    { EPRINTF(EPLACE,"%s: Failed to create onecode schema\n",Prog_Name);
      goto error;
    }

  for (i = tmax; i >= 0; i--)
    { input = fopen(MyCatenate(SPATH,"/",SROOT,suffix[i]),"r");
      if (input != NULL)
        break;
    }
  if (i < 0)
    { SEXTN = find_1seq(SPATH,SROOT);
      if (SEXTN == ERROR)
        goto error;
      else if (SEXTN == NULL)
        { EPRINTF(EPLACE,"%s: Could not find %s or a fasta or 1-code extension there of\n",
                         Prog_Name,SROOT);
          oneSchemaDestroy(schema);
          goto error;
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
        { EPRINTF(EPLACE,"%s: Could not find valid extension of %s\n",Prog_Name,SROOT);
          oneSchemaDestroy(schema);
          EXIT(-1);
        }
      e[-suflen[i]] = '\0';
      sdot = e-suflen[i];
      i += 4;
      SEXTN = suffix[i];
    }
  else if (i == 0)
    { e = SROOT + strlen(SROOT);
      for (i = tmax; i >= 1; i--)
        if (strcmp(e-suflen[i],suffix[i]) == 0 && i != 4)
          break;
      if (i > 0)
        { e[-suflen[i]] = '\0';
          sdot = e-suflen[i];
          SEXTN = suffix[i];
        }
      else
        { OneFile *of;

          of = oneFileOpenRead(MyCatenate(SPATH,"/",SROOT,""),schema,"seq",1);
          if (of == NULL)
            { EPRINTF(EPLACE,"%s: Could not find valid extension of %s\n",Prog_Name,SROOT);
              oneSchemaDestroy(schema);
              EXIT(-1);
            }
          oneFileClose(of);
          SEXTN = "";
          i = -1;
        }
    }

  oneSchemaDestroy(schema);

  *spath = Strdup(MyCatenate(SPATH,"/",SROOT,SEXTN),"Allocating source name");
  if (*spath == NULL)
    goto error;

  if (i >= 8)
    { *tpath = Strdup(*spath,"Allocating target name");
      if (*tpath == NULL)
        { free(*spath);
          goto error;
        }

      if (sdir != NULL) *sdir = '/';
      if (sdot != NULL) *sdot = '.';
      if (tdir != NULL) *tdir = '/';
      if (tdot != NULL) *tdot = '.';
      return (IS_GDB);
    }

  TEXTN = ".1gdb";
  if (target == NULL)
    { TPATH = SPATH;
      TROOT = SROOT;
    }
  else // argc == 3
    { TPATH = target;
      if (stat(TPATH,&status) >= 0 && (status.st_mode & S_IFMT) == S_IFDIR)
        TROOT = SROOT;
      else
        { p = rindex(TPATH,'/');
          if (p == NULL)
            { TROOT = TPATH;
              TPATH = ".";
            }
          else
            { *p = '\0';
              tdir = p;
              TROOT = p+1;
            }
          p = rindex(TROOT,'.');
          if (p != NULL)
            { if (strcmp(p+1,"gdb") == 0)
                { *p = '\0';
                  tdot = p;
                  TEXTN = ".gdb";
                }
              else if (strcmp(p+1,"1gdb") == 0)
                { *p = '\0';
                  tdot = p;
                }
            }
        }
    }

  *tpath = Strdup(MyCatenate(TPATH,"/",TROOT,TEXTN),"Allocating target name");
  if (*tpath == NULL)
    { free(*spath);
      goto error;
    }

  if (sdir != NULL) *sdir = '/';
  if (sdot != NULL) *sdot = '.';
  if (tdir != NULL) *tdir = '/';
  if (tdot != NULL) *tdot = '.';

  if (i < 0)
    return (IS_ONE);
  else if (i < 4)
    return (IS_FA);
  else
    return (IS_FA_GZ);

error:
  if (sdir != NULL) *sdir = '/';
  if (sdot != NULL) *sdot = '.';
  if (tdir != NULL) *tdir = '/';
  if (tdot != NULL) *tdot = '.';
  EXIT (-1);
}


/*******************************************************************************************
 *
 *  GDB CREATION FROM SOURCE
 *
 ********************************************************************************************/

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

//  Read next line into a buffer and return a pointer to the buffer and set *plen
//    the length of the line.  NB: replaces '\n' with '\0'.

static char *read_line(void *input, int gzipd, int *plen, int nline, char *spath)
{ static char *buffer;
  static int   bmax = 0;
  int len;

  if (bmax == 0)
    { bmax = 500;
      buffer = (char *) malloc(bmax);
      if (buffer == NULL)
        { EPRINTF(EPLACE,"%s: Out of memory reading %s\n",Prog_Name,spath);
          return (ERROR);
        }
    }

  if (gzipd)
    { if (gzgets(input,buffer,bmax) == NULL)
        { if (gzeof(input))
            return (NULL);
          EPRINTF(EPLACE,"%s: Could not read next line, %d, of file %s\n",Prog_Name,nline,spath);
          return (ERROR);
        }
    }
  else
    { if (fgets(buffer,bmax,input) == NULL)
        { if (feof(input))
            return (NULL);
          EPRINTF(EPLACE,"%s: Could not read next line, %d, of file %s\n",Prog_Name,nline,spath);
          return (ERROR);
        }
    }

  len = strlen(buffer);
  while (buffer[len-1] != '\n')
    { bmax = ((int) (1.4*bmax)) + 100;
      buffer = (char *) realloc(buffer,bmax);
      if (buffer == NULL)
        { EPRINTF(EPLACE,"%s: Out of memory reading %s\n",Prog_Name,spath);
          return (ERROR);
        }
      if (gzipd)
        { if (gzgets(input,buffer+len,bmax-len) == NULL)
            { if (gzeof(input))
                EPRINTF(EPLACE,"%s: Last line %d of file %s does not end with new-line\n",
                               Prog_Name,nline,spath);
              else
                EPRINTF(EPLACE,"%s: Could not read next line, %d, of file %s\n",
                               Prog_Name,nline,spath);
              return (ERROR);
            }
        }
      else
        { if (fgets(buffer+len,bmax-len,input) == NULL)
            { if (feof(input))
                EPRINTF(EPLACE,"%s: Last line %d of file %s does not end with new-line\n",
                               Prog_Name,nline,spath);
              else
                EPRINTF(EPLACE,"%s: Could not read next line, %d, of file %s\n",
                               Prog_Name,nline,spath);
              return (ERROR);
            }
        }
      len += strlen(buffer+len);
    }
  buffer[--len] = '\0';

  if (plen != NULL)
    *plen = len;
  return (buffer);
}

FILE **Create_GDB(GDB *gdb, char *spath, int ftype, int bps, char *tpath, int nthresh)
{ GDB_SCAFFOLD  *scaffs;
  GDB_CONTIG    *contigs;
  OneProvenance *prov;
  char          *headers, *seqpath;
  int            ctgtop, scftop;
  int64          hdrtop;
  int64          count[4];
  FILE          *bases;
  int64          hdrtot, maxctg, seqtot, boff;
  int            ncontig, nscaff, nprov;
  int            len, clen;
  int64          spos;

  int            gzipd, nline;
  void          *input;
  OneSchema     *schema;
  OneFile       *of;

  prov  = NULL;
  input  = NULL;
  schema = NULL;
  of     = NULL;

  //  Establish .bps file if needed

  bases = NULL;
  if (bps != 0)
    { if (tpath == NULL)
        seqpath = Numbered_Suffix("._gdb.",getpid(),".bps");
      else
        { char *path = PathTo(tpath);
          char *root = Root(tpath,".1gdb");
          if (root == NULL || path == NULL)
            { free(path);
              free(root);
              EXIT (NULL);
            }
          seqpath = MyCatenate(path,"/.",root,".bps");
          free(root);
          free(path);
        }
      if (seqpath == NULL)
        EXIT (NULL);
      seqpath = strdup(seqpath);
      if (seqpath == NULL)
        EXIT(NULL);
      bases = Fopen(seqpath,"w+");
      if (bases == NULL)
        { free(seqpath);
          EXIT(NULL);
        }
    }
  else
    seqpath = NULL;

  //  Setup expanding arrays for headers, scaffolds, & contigs

  hdrtot  = 0;
  maxctg  = 0;
  ncontig = 0;
  nscaff  = 0;
  seqtot  = 0;
  boff    = 0;

  count[0] = 0;
  count[1] = 0;
  count[2] = 0;
  count[3] = 0;

  hdrtop  = 100000;
  ctgtop  = 5000;
  scftop  = 1000;
  contigs = malloc(ctgtop*sizeof(GDB_CONTIG));
  scaffs  = malloc(scftop*sizeof(GDB_CONTIG));
  headers = malloc(hdrtop);
  if (*spath != '/')
    gdb->srcpath = strdup(MyCatenate(getcwd(NULL,0),"/",spath,""));
  else
    gdb->srcpath = strdup(spath);
  if (contigs == NULL || scaffs == NULL || headers == NULL || gdb->srcpath == NULL)
    { EPRINTF(EPLACE,"%s: Out of memory creating GDB for %s\n",Prog_Name,spath);
      goto error;
    }

  //  Read either 1-file or fasta and accumulate skeleton

  if (ftype == IS_ONE)

    { int   i, isScaffold;
      int   byte_cnt[256];     // times byte occurs in compressed sequence
      int   psize;
      char *pstr;

      schema = make_Seq_Schema();
      if (schema == NULL)
        { EPRINTF(EPLACE,"%s: Failed to create onecode schema (Create_GDB)\n",Prog_Name);
          goto error;
        }

      of = oneFileOpenRead (spath, schema, "seq", 1) ;
      if (of == NULL)
        { EPRINTF(EPLACE,"%s: Cannot open %s as a 1-code sequence file (Create_GDB)\n",
                         Prog_Name,spath);
          goto error;
        }
    
      for (i = 0; i < 256; i++)
        byte_cnt[i] = 0;

      nprov   = of->info['!']->accum.count;

      psize = nprov * sizeof(OneProvenance);
      for (i = 0; i < nprov; i++)
        psize += strlen(of->provenance[i].program)
               + strlen(of->provenance[i].version)
               + strlen(of->provenance[i].command)
               + strlen(of->provenance[i].date) + 4;
      if (psize > 0)
        { prov = malloc(psize);
          if (prov == NULL)
            goto error;

          pstr = (char *) (prov + nprov);
          for (i = 0; i < nprov; i++)
            { prov[i].program = pstr;
              pstr = stpcpy(pstr,of->provenance[i].program) + 1;
              prov[i].version = pstr;
              pstr = stpcpy(pstr,of->provenance[i].version) + 1;
              prov[i].command = pstr;
              pstr = stpcpy(pstr,of->provenance[i].command) + 1;
              prov[i].date = pstr;
              pstr = stpcpy(pstr,of->provenance[i].date) + 1;
            }
        }
      else
        prov = NULL;
    
      // line sequence must be either (I!S)+, or (s(n!(Sn)*Sn!)+, e.g.
      //   - scaffold files with s (scaffold id) and n (inter-contig N's) and no I lines, or
      //   - straight sequence files with no s and n lines, and only S lines, allowing I lines
    
      // Ignore N characters in sequences (non-acgt chars within contigs)
    
      if (of->info['s']->given.count > 0 && of->info['I']->given.count == 0)
        isScaffold = true;
      else if (of->info['s']->given.count == 0 && of->info['n']->given.count == 0)
        isScaffold = false;
      else
        { EPRINTF(EPLACE,"%s: 1seq file %s has incomplete scaffold properties - can't use it\n",
                         Prog_Name,spath);
          EPRINTF(EPLACE,"          I.E. #'s' = %lld  #'I' = %lld  #'n' = %lld\n",
                   of->info['s']->given.count, of->info['I']->given.count,
                   of->info['n']->given.count) ;
          goto error;
        }
    
      spos = 0;
      while (oneReadLine(of))
        switch (of->lineType)
    
        { case 's': // scaffold name
            if (boff != 0)
              { scaffs[nscaff].ectg = ncontig;
                scaffs[nscaff].slen = spos;
                nscaff += 1;
              }
            if (nscaff >= scftop)
              { scftop = 1.2*nscaff + 500;
                scaffs = realloc(scaffs,scftop*sizeof(GDB_SCAFFOLD));
                if (scaffs == NULL)
                  { EPRINTF(EPLACE,"%s: Out of memory creating GDB for %s\n",Prog_Name,spath);
                    goto error;
                  }
              }
            scaffs[nscaff].hoff = hdrtot;
            scaffs[nscaff].fctg = ncontig;
            spos = 0;
    
            len = oneLen(of);
            memcpy(headers+hdrtot,oneString(of),len);
            hdrtot += len;
            headers[hdrtot++] = '\0';
            break;
    
          case 'n': // non-acgt chars outside sequences
            spos += oneInt(of,1);
            break;
    
          case 'S': // sequence object
            len = oneLen(of);
    
            if (ncontig >= ctgtop)
              { ctgtop = 1.2*ncontig + 1000;
                contigs = realloc(contigs,ctgtop*sizeof(GDB_CONTIG));
                if (contigs == NULL)
                  { EPRINTF(EPLACE,"%s: Out of memory creating GDB for %s\n",Prog_Name,spath);
                    goto error;
                  }
              }
            contigs[ncontig].boff = boff;
            contigs[ncontig].clen = len;
            contigs[ncontig].sbeg = spos;
            contigs[ncontig].scaf = nscaff;
            ncontig += 1;
    
            clen = COMPRESSED_LEN(len);
            uint8* bytes = oneDNA2bit(of);
    
            seqtot += len;
            if (len > maxctg)
              maxctg = len;
            for (i = 0 ; i < clen; i++)
              byte_cnt[bytes[i]] += 1;
            if (len & 0x3)
              count[0] -= (4-(len&0x3));
    
            if (bps)
              { if (fwrite (bytes, clen, 1, bases) != 1)
                  goto error;
              }
            boff += clen;
    
            if (isScaffold)
              spos += len;
            break;
    
          case 'I': // identifier of sequence
            if (boff != 0)
              { scaffs[nscaff].ectg = ncontig;
                scaffs[nscaff].slen = spos;
                nscaff += 1;
              }
            scaffs[nscaff].hoff = hdrtot;
            scaffs[nscaff].fctg = ncontig;
            spos = 0;
    
            len = oneLen(of);
            memcpy(headers+hdrtot,oneString(of),len);
            hdrtot += len;
            headers[hdrtot++] = '\0';
            break;
          }
      if (boff != 0)
        { scaffs[nscaff].ectg = ncontig;
          scaffs[nscaff].slen = spos;
          nscaff += 1;
        }
    
      oneFileClose(of);
      oneSchemaDestroy(schema);
    
      for (i = 0; i < 256; i++)
        { count[i & 0x3] += byte_cnt[i];
          count[(i>>2) & 0x3] += byte_cnt[i];
          count[(i>>4) & 0x3] += byte_cnt[i];
          count[(i>>6) & 0x3] += byte_cnt[i];
        }
    }

  else  //  Fasta reader

    { char          bpsbuf[1024];
      int           bpscur;
      char         *line;
      int           i, x, m, in;
      int           g, lste;
      uint8         byte, bfl;

      byte    = 0;
      bpscur  = 0;

      gzipd = (ftype == IS_FA_GZ);
      if (gzipd)
        input = gzopen(spath,"r");
      else
        input = fopen(spath,"r");
      if (input == NULL)
        { EPRINTF(EPLACE,"%s: Cannot open %s for reading\n",Prog_Name,spath);
          goto error;
        }

      nprov = 0;
      prov  = NULL;
    
      //  Get the header of the first line.  If the file is empty skip.
    
      nline = 1;
      line  = read_line(input,gzipd,&len,nline,spath);
      if (line == ERROR)
        goto error;
      else if (line == NULL)
        { EPRINTF(EPLACE,"%s: Input %s is empty, terminating!\n",Prog_Name,spath);
          goto error;
        }
    
      if (line[0] != '>')
        { EPRINTF(EPLACE,"%s: First header in fasta file %s is missing\n",
                         Prog_Name,spath);
          goto error;
        }
    
      while (line != NULL)
        { for (i = 1; i < len; i++)
            if (!isspace(line[i]))
              break;
          len -= i;
    
          if (nscaff >= scftop)
            { scftop = 1.2*nscaff + 500;
              scaffs = realloc(scaffs,scftop*sizeof(GDB_SCAFFOLD));
              if (scaffs == NULL)
                { EPRINTF(EPLACE,"%s: Out of memory creating GDB for %s\n",Prog_Name,spath);
                  goto error;
                }
            }
          scaffs[nscaff].fctg = ncontig;
          scaffs[nscaff].hoff = hdrtot; 
    
          if (hdrtot + len + 1 > hdrtop)
            { hdrtop = 1.2*(hdrtot+len+1) + 10000;
              headers = realloc(headers,hdrtop);
              if (headers == NULL)
                { EPRINTF(EPLACE,"%s: Out of memory creating GDB for %s\n",Prog_Name,spath);
                  goto error;
                }
            }
          memcpy(headers+hdrtot,line+i,len);
          hdrtot += len;
          headers[hdrtot++] = '\0';

          bfl  = 0;
          in   = 0;
          spos = 0;
          lste = -nthresh;
          clen = 0;
          while (1)
            { nline += 1;
              line   = read_line(input,gzipd,&len,nline,spath);
              if (line == ERROR)
                goto error;
              else if (line == NULL || line[0] == '>')
                { if (in)
                    { if (bps)
                        { if ((clen & 0x7) != 0)
                            { if (bpscur >= 1024)
                                { fwrite(bpsbuf,1024,1,bases);
                                  bpscur = 0;
                                }
                              bpsbuf[bpscur++] = byte;
                              boff += 1;
                            }
                          clen >>= 1;
                        }
                      spos += clen;
                      contigs[ncontig].clen = clen;
                      seqtot += clen;
                      if (clen > maxctg)
                        maxctg = clen;
                      ncontig += 1;
                      if (ncontig >= ctgtop)
                        { ctgtop = 1.2*ncontig + 1000;
                          contigs = realloc(contigs,ctgtop*sizeof(GDB_CONTIG));
                          if (contigs == NULL)
                            { EPRINTF(EPLACE,"%s: Out of memory creating GDB for %s\n",
                                             Prog_Name,spath);
                              goto error;
                            }
                        }
                    }
                  else
                    { if (bps && bfl)
                        { if (bpscur >= 1024)
                            { fwrite(bpsbuf,1024,1,bases);
                              bpscur = 0;
                            }
                          bpsbuf[bpscur++] = byte;
                          boff += 1;
                        }
                      if (spos - lste < nthresh)
                        spos = lste;
                    }
                  break;
                }
    
              if (bps)
                for (i = 0; i < len; i++)
                  { x = number[(int) line[i]];
                    if (x < 4)
                      { if (!in)

                          { if (spos - lste < nthresh)
                              { if (ncontig == 0)
                                  clen = 0;
                                else
                                  { ncontig -= 1;
                                    clen = contigs[ncontig].clen;
                                    seqtot -= clen;
                                    clen <<= 1;
                                  }
                                for (g = spos-lste; g > 0; g--)
                                  { count[0] += 1;
                                    m = (clen & 0x7);
                                    if (m == 0)
                                      byte = 0;
                                    else if (m == 6)
                                      { if (bpscur >= 1024)
                                          { fwrite(bpsbuf,1024,1,bases);
                                            bpscur = 0;
                                          }
                                        bpsbuf[bpscur++] = byte;
                                        boff += 1;
                                      }
                                    clen += 2;
                                  }
                                spos -= (clen >> 1);
                              }
                            else
                              { if (bfl)
                                  { if (bpscur >= 1024)
                                      { fwrite(bpsbuf,1024,1,bases);
                                        bpscur = 0;
                                      }
                                    bpsbuf[bpscur++] = byte;
                                    boff += 1;
                                  }
                                if (ncontig >= ctgtop)
                                  { ctgtop = 1.2*ncontig + 1000;
                                    contigs = realloc(contigs,ctgtop*sizeof(GDB_CONTIG));
                                    if (contigs == NULL)
                                      { EPRINTF(EPLACE,"%s: Out of memory creating GDB for %s\n",
                                                       Prog_Name,spath);
                                        goto error;
                                      }
                                  }
                                contigs[ncontig].sbeg = spos;
                                contigs[ncontig].boff = boff;
                                contigs[ncontig].scaf = nscaff;
                                clen = 0;
                              }
                            in = 1;
                          }

                        count[x] += 1;
                        m = (clen & 0x7);
                        if (m == 0)
                          byte = x;
                        else
                          { byte |= (x << m);
                            if (m == 6)
                              { if (bpscur >= 1024)
                                  { fwrite(bpsbuf,1024,1,bases);
                                    bpscur = 0;
                                  }
                                bpsbuf[bpscur++] = byte;
                                boff += 1;
                              }
                          }
                        clen += 2;
                      }
                    else
                      { if (in)
                          { bfl = ((clen & 0x7) != 0);
                            clen >>= 1;
                            spos += clen;
                            lste  = spos;
                            contigs[ncontig].clen = clen;
                            seqtot += clen;
                            if (clen > maxctg)
                              maxctg = clen;
                            ncontig += 1;
                            in = 0;
                          }
                        spos += 1;
                      }
                  }
              else
                for (i = 0; i < len; i++)
                  { x = number[(int) line[i]];
                    if (x < 4)
                      { if (!in)
                          { if (spos - lste < nthresh)
                              { if (ncontig == 0)
                                  clen = 0;
                                else
                                  { ncontig -= 1;
                                    clen = contigs[ncontig].clen;
                                    seqtot -= clen;
                                    clen += spos-lste;
                                  }
                                spos -= clen;
                              }
                            else
                              { if (ncontig >= ctgtop)
                                  { ctgtop = 1.2*ncontig + 1000;
                                    contigs = realloc(contigs,ctgtop*sizeof(GDB_CONTIG));
                                    if (contigs == NULL)
                                      { EPRINTF(EPLACE,"%s: Out of memory creating GDB for %s\n",
                                                       Prog_Name,spath);
                                        goto error;
                                      }
                                  }
                                contigs[ncontig].sbeg = spos;
                                contigs[ncontig].boff = boff;
                                contigs[ncontig].scaf = nscaff;
                                clen = 0;
                                in   = 1;
                              }
                          }
                        count[x] += 1;
                        clen += 1;
                      }
                    else
                      { if (in)
                          { spos += clen;
                            contigs[ncontig].clen = clen;
                            seqtot += clen;
                            if (clen > maxctg)
                              maxctg = clen;
                            ncontig += 1;
                            in = 0;
                          }
                        spos += 1;
                      }
                  }
            }
    
          if (spos == 0)
            { EPRINTF(EPLACE,"%s: Missing sequence entry at line %d in file %s\n",
                             Prog_Name,nline,spath);
              goto error;
            }
          scaffs[nscaff].slen = spos;
          scaffs[nscaff].ectg = ncontig;
          nscaff += 1;
        }
    
      if (bpscur > 0)
        fwrite(bpsbuf,bpscur,1,bases);

      if (gzipd)
        gzclose(input);
      else
        fclose(input);
    }

  if (bps > 0)
    rewind(bases);

  gdb->nprov = nprov;
  gdb->prov  = prov;

  gdb->nscaff  = nscaff;
  gdb->ncontig = ncontig;
  gdb->maxctg  = maxctg;
  gdb->hdrtot  = hdrtot;
  gdb->seqtot  = seqtot;

  gdb->scaffolds = scaffs;
  gdb->contigs   = contigs;
  gdb->headers   = headers;
  gdb->seqstate  = EXTERNAL;
  gdb->seqsrc    = ftype;
  gdb->seqpath   = seqpath;
  gdb->seqs      = bases;

  gdb->freq[0] = (1.*count[0])/seqtot;
  gdb->freq[1] = (1.*count[1])/seqtot;
  gdb->freq[2] = (1.*count[2])/seqtot;
  gdb->freq[3] = (1.*count[3])/seqtot;

  if (bps <= 1)
    { if (bps == 1 && tpath == NULL)
        unlink(seqpath);
      return ((FILE **) &gdb->seqs);
    }
  else
    { FILE **units;
      int    i;

      units = malloc(sizeof(FILE *)*bps);
      if (units == NULL)
        goto error;
      units[0] = bases;
      for (i = 1; i < bps; i++)
        units[i] = Fopen(seqpath,"r");
      if (tpath == NULL)
        unlink(seqpath);
      return (units);
    }

error:
  if (of != NULL)
    oneFileClose(of);
  if (schema != NULL)
    oneSchemaDestroy(schema);
  if (input != NULL)
    { if (gzipd)
        gzclose(input);
      else
        fclose(input);
    }
  free(prov);
  free(gdb->srcpath);
  free(headers);
  free(scaffs);
  free(contigs);
  if (bases != NULL)
    fclose(bases);
  if (tpath == NULL)
    unlink(seqpath);
  free(seqpath);
  EXIT (NULL);
}

/*******************************************************************************************
 *
 *  GDB OPEN & CLOSE ROUTINES
 *
 ********************************************************************************************/

// Open the given database "path" into the supplied GDB record "gdb".
//   Initially the sequence data, if any, stays in .bps file with a FILE pointer to it.
// Return values in interactive mode:
//     0: Open of GDB proceeded without mishap
//     1: The GDB could not be opened, a message why is in EPLACE

int Read_GDB(GDB *gdb, char *path)
{ OneSchema     *schema;
  OneFile       *of;
  FILE          *seqs;
  GDB_SCAFFOLD  *scf;
  GDB_CONTIG    *ctg;
  OneProvenance *prov;
  char          *hdr, *srcpath, *seqpath;
  int            nscaff, ncontig, nprov;
  int64          len, seqtot, hdrtot, maxctg, boff, spos, psize;

  { char *e;
    char *root, *pwd;
    char *fname;

    pwd = PathTo(path);
    e = path + strlen(path);
    if (strcmp(e-5,".1gdb") == 0)
      root = Root(path,".1gdb");
    else
      root = Root(path,".gdb");
    fname = MyCatenate(pwd,"/",root,".1gdb");
    if (pwd == NULL || root == NULL || fname == NULL)
      { free(root);
        free(path);
        EXIT(1);
      }
    seqs  = fopen(fname,"r");
    if (seqs == NULL)
      { fname = MyCatenate(pwd,"/",root,".1gdb");
        if (fname == NULL)
          { free(root);
            free(path);
            EXIT(1);
          }
        seqs = fopen(fname,"r");
        if (seqs == NULL)
          { EPRINTF(EPLACE,"%s: Cannot find/open GDB file %s\n",Prog_Name,path);
            free(root);
            free(pwd);
            EXIT(1);
          }
        fclose(seqs);
      }
    else
      fclose(seqs);

    schema = make_GDB_Schema();
    if (schema == NULL)
      { EPRINTF(EPLACE,"%s: Failed to create gdb schema\n",Prog_Name);
        EXIT(1);
      }
  
    of = oneFileOpenRead(fname,schema,"gdb",1);
    if (of == NULL)
      { EPRINTF(EPLACE,"%s: Failed to open .1gdb file %s\n",Prog_Name,path);
        oneSchemaDestroy(schema);
        free(root);
        free(pwd);
        EXIT(1);
      }

    seqpath = MyCatenate(pwd,"/.",root,".bps");
    free(root);
    free(pwd);
    if (seqpath == NULL)
      { oneFileClose(of);
        oneSchemaDestroy(schema);
        EXIT(1);
      }
    seqs = fopen(seqpath,"r");
    if (seqs == NULL)
      { EPRINTF(EPLACE,"%s: Failed to open .bps file for GDB %s\n",Prog_Name,path);
        oneFileClose(of);
        oneSchemaDestroy(schema);
        EXIT(1);
      }
  }

  nprov   = of->info['!']->accum.count;
  nscaff  = of->info['S']->given.count;
  ncontig = of->info['C']->given.count;
  hdrtot  = of->info['S']->given.total + nscaff;

  { int i;

    psize = nprov * sizeof(OneProvenance);
    for (i = 0; i < nprov; i++)
      psize += strlen(of->provenance[i].program)
             + strlen(of->provenance[i].version)
             + strlen(of->provenance[i].command)
             + strlen(of->provenance[i].date) + 4;
  }

  seqpath = strdup(seqpath);
  srcpath = strdup(of->reference[0].filename);
  scf   = malloc(sizeof(GDB_SCAFFOLD)*nscaff);
  ctg   = malloc(sizeof(GDB_CONTIG)*ncontig);
  hdr   = malloc(hdrtot);
  if (psize > 0)
    prov  = malloc(psize);
  else
    prov  = NULL;
  if (seqpath == NULL || srcpath == NULL || scf == NULL || ctg == NULL || hdr == NULL
                      || (psize > 0 && prov == NULL))
    { EPRINTF(EPLACE,"%s: Could not allocate memory for GDB (Read_GDB)\n",Prog_Name);
      goto error;
    }

  { int   i;
    char *pstr;

    pstr = (char *) (prov + nprov);
    for (i = 0; i < nprov; i++)
      { prov[i].program = pstr;
        pstr = stpcpy(pstr,of->provenance[i].program) + 1;
        prov[i].version = pstr;
        pstr = stpcpy(pstr,of->provenance[i].version) + 1;
        prov[i].command = pstr;
        pstr = stpcpy(pstr,of->provenance[i].command) + 1;
        prov[i].date = pstr;
        pstr = stpcpy(pstr,of->provenance[i].date) + 1;
      }
  }

  nscaff  = -1;
  ncontig = 0;
  hdrtot  = 0;
  seqtot  = 0;
  maxctg  = 0;
  boff = 0;
  spos = 0;
  while (oneReadLine(of))
    switch (of->lineType)
    { case 'f':
        gdb->freq[0] = oneReal(of,0);
        gdb->freq[1] = oneReal(of,1);
        gdb->freq[2] = oneReal(of,2);
        gdb->freq[3] = oneReal(of,3);
        break;
      case 'S':
        if (nscaff >= 0)
          { scf[nscaff].ectg = ncontig;
            scf[nscaff].slen = spos;
            spos = 0;
          }
        nscaff += 1;
        scf[nscaff].hoff = hdrtot;
        scf[nscaff].fctg = ncontig;
        len = oneLen(of);
        memcpy(hdr+hdrtot,oneString(of),len);
        hdrtot += len;
        hdr[hdrtot++] = '\0';
        break;
      case 'G':
        spos += oneInt(of,0);
        break;
      case 'C':
        len = oneInt(of,0);
        ctg[ncontig].boff = boff;
        ctg[ncontig].sbeg = spos;
        ctg[ncontig].clen = len;
        ctg[ncontig].scaf = nscaff;
        ncontig += 1;
        if (len > maxctg)
          maxctg = len;
        seqtot  += len;
        boff    += COMPRESSED_LEN(len);
        spos    += len;
        break;
    }
  scf[nscaff].ectg = ncontig;
  scf[nscaff].slen = spos;
  nscaff += 1;

  gdb->nprov = nprov;
  gdb->prov  = prov;

  gdb->nscaff    = nscaff;
  gdb->scaffolds = scf;

  gdb->ncontig = ncontig;
  gdb->maxctg  = maxctg;
  gdb->contigs = ctg;

  gdb->srcpath  = srcpath;
  gdb->seqpath  = seqpath;

  gdb->hdrtot  = hdrtot;
  gdb->headers = hdr;

  gdb->seqtot   = seqtot;
  gdb->seqstate = EXTERNAL;
  gdb->seqs     = seqs;

  srcpath += strlen(srcpath);
  if (strcmp(srcpath-3,".gz") == 0)
    gdb->seqsrc = IS_FA_GZ;
  else if (strcmp(srcpath-3,".fa") == 0)
    gdb->seqsrc = IS_FA;
  else if (strcmp(srcpath-4,".fna") == 0)
    gdb->seqsrc = IS_FA;
  else if (strcmp(srcpath-6,".fasta") == 0)
    gdb->seqsrc = IS_FA;
  else
    gdb->seqsrc = IS_ONE;

  oneFileClose(of);
  oneSchemaDestroy(schema);
  return (0);

error:
  free(hdr);
  free(ctg);
  free(scf);
  free(srcpath);
  free(seqpath);
  if (seqs != NULL)
    fclose(seqs);
  oneFileClose(of);
  oneSchemaDestroy(schema);
  EXIT(1);
}

static void Print_Contig(char *s, int width)
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

// The GDB moves it's sequence data from the .bps file to an in-memory block and also converts
//   it to the requested format.  The boff field of contigs are adjusted, if necessary, so that
//   the give the offset in the memory block of the associated string.
// If in interactive mode, NULL is returned on error.

int Load_Sequences(GDB *gdb, int stype)
{ FILE       *b       = (FILE *) gdb->seqs;
  int         ncontig = gdb->ncontig;
  GDB_CONTIG *c       = gdb->contigs;
  void       (*translate)(char *s);
  struct stat state;

  char  *seq;
  int64  off;
  int    i, len, clen;

  if (b == NULL)
    { EPRINTF(EPLACE,"%s: GDB has no sequence data (Load_Sequences)\n",Prog_Name);
      EXIT(1);
    }
  if (stype < COMPRESSED || stype > UPPER_CASE)
    { EPRINTF(EPLACE,"%s: Invalid load type %d (Load_Sequences)\n",Prog_Name,stype);
      EXIT(1);
    }
  if (gdb->seqstate != EXTERNAL)
    { EPRINTF(EPLACE,"%s: GDB's sequencing info already loaded (Load_Sequences)\n",Prog_Name);
      EXIT(1);
    }

  if (stype == COMPRESSED)
    { if (fstat(fileno(b),&state) < 0)
        { EPRINTF(EPLACE,"%s: Cannot fetch size of GDB's base pair file (Load_Sequences)\n",
                         Prog_Name);
          EXIT(1);
        }
      seq = (char *) malloc(state.st_size);
      if (seq == NULL)
        { EPRINTF(EPLACE,"%s: Cannot allocate in-memory sequence array for GDB (Load_Sequences)\n",
                         Prog_Name);
          EXIT(1);
        }
      if (fread(seq,1,state.st_size,b) != (size_t) state.st_size)
        { EPRINTF(EPLACE,"%s: Cannot read sequence file of GDB (Load_Sequences)\n",Prog_Name);
          free(seq);
          EXIT(1);
        }

      gdb->seqstate = COMPRESSED;
      gdb->seqs     = (void *) seq; 
      return (0);
    }

  seq = (char *) malloc(gdb->seqtot+gdb->ncontig+4);
  if (seq == NULL)
    { EPRINTF(EPLACE,"%s: Cannot allocate in-memory sequence array for GDB (Load_Sequences)\n",
                     Prog_Name);
      EXIT(1);
    }

  *seq++ = 0;
  if (stype == LOWER_CASE)
    translate = Lower_Read;
  else if (stype == UPPER_CASE)
    translate = Upper_Read;
  else
    { seq[-1] = 4;
      translate = Upper_Read;
    }

  rewind(b);
  off = 0;
  for (i = 0; i < ncontig; i++)
    { len  = c[i].clen;
      clen = COMPRESSED_LEN(len);
      if (clen > 0)
        { if (fread(seq+off,clen,1,b) != 1)
            { EPRINTF(EPLACE,"%s: Read of .bps file failed (Load_All_Sequences)\n",Prog_Name);
              free(seq-1);
              EXIT(1);
            }
        }
      Uncompress_Read(len,seq+off);
      if (stype > NUMERIC)
        translate(seq+off);
      c[i].boff = off;
      off += (len+1);
    }

  gdb->seqstate = stype;
  gdb->seqs     = (void *) seq; 

  return (0);
}

// Write the given gdb to the file 'tpath'.  The GDB must have seqstate EXTERNAL and tpath
//   must be consistent with the name of the .bps file.

extern bool addProvenance(OneFile *of, OneProvenance *from, int n) ; // backdoor - clean up some day

int Write_GDB(GDB *gdb, char *tpath)
{ OneSchema *schema;
  OneFile   *of;
  bool       binary;
  int64      spos, len;
  char      *head;
  int        s, c;

  if (gdb->seqstate != EXTERNAL)
    { EPRINTF(EPLACE,"%s: GDB must be in EXTERNAL state (Write_GDB)\n",Prog_Name);
      EXIT(1);
    }

  { char *e;
    char *root, *pwd;
    FILE *dbvis;

    pwd = PathTo(tpath);
    e = tpath + strlen(tpath);
    if (strcmp(e-4,".gdb") == 0)
      { root = Root(tpath,".gdb");
        binary = false;
      }
    else
      { root = Root(tpath,".1gdb");
        binary = true;
      }
    head = MyCatenate(pwd,"/.",root,".bps");
    if (root == NULL || pwd == NULL || head == NULL)
      { free(pwd);
        free(root);
        EXIT(1);
      }
    free(pwd);
    free(root);
    if ((dbvis = fopen(head,"r")) == NULL)
      { EPRINTF(EPLACE,"%s: Could not find .bps %s file for write (Write_GDB)\n",Prog_Name,tpath);
        EXIT(1);
      }
    fclose(dbvis);
  }

  schema = make_GDB_Schema();
  if (schema == NULL)
    { EPRINTF(EPLACE,"%s: Failed to create GDB schema (Write_GDB)\n",Prog_Name);
      EXIT(1);
    }

  of = oneFileOpenWriteNew(tpath,schema,"gdb",binary,1);
  if (of == NULL)
    { EPRINTF(EPLACE,"%s: Failed to open GDB file %s (Write_GDB)\n",Prog_Name,tpath);
      oneSchemaDestroy(schema);
      EXIT(1);
    }

  addProvenance(of,gdb->prov,gdb->nprov);
  oneAddProvenance(of,Prog_Name,"0.1",Command_Line);

  oneAddReference(of,gdb->srcpath,1);

  oneReal(of,0) = gdb->freq[0];
  oneReal(of,1) = gdb->freq[1];
  oneReal(of,2) = gdb->freq[2];
  oneReal(of,3) = gdb->freq[3];
  oneWriteLine(of,'f',0,0);

  for (s = 0; s < gdb->nscaff; s++)
    { head = gdb->headers + gdb->scaffolds[s].hoff;
      oneWriteLine(of,'S',strlen(head),head);

      spos = 0;
      for (c = gdb->scaffolds[s].fctg; c < gdb->scaffolds[s].ectg; c++)
        { if (gdb->contigs[c].sbeg > spos)
            { oneInt(of,0) = gdb->contigs[c].sbeg - spos;
              oneWriteLine(of,'G',0,0);
            }
          len = gdb->contigs[c].clen;
          oneInt(of,0) = len;
          oneWriteLine(of,'C',0,0);
          spos = gdb->contigs[c].sbeg + len;
        }
      if (gdb->scaffolds[s].slen > spos)
        { oneInt(of,0) = gdb->scaffolds[s].slen - spos;
          oneWriteLine(of,'G',0,0);
        }
    }

  oneFileClose(of);
  oneSchemaDestroy(schema);
  return (0);
}

// Shut down an open 'gdb' by freeing all associated space and the file pointers

void Close_GDB(GDB *gdb)
{ if (gdb->seqs != NULL)
    { if (gdb->seqstate == EXTERNAL)
        fclose(gdb->seqs);
      else if (gdb->seqstate == COMPRESSED)
        free(gdb->seqs);
      else
        free(gdb->seqs-1);
    }
  free(gdb->headers);
  free(gdb->contigs);
  free(gdb->scaffolds);
  free(gdb->srcpath);
  free(gdb->seqpath);
  if (gdb->nprov > 0)
    free(gdb->prov);
}


/*******************************************************************************************
 *
 *  CREATE OR READ GDB FOR 1-ALN INTERPRETATION
 *
 ********************************************************************************************/

FILE **Get_GDB(GDB *gdb, char *source, char *cpath, int num_bps)
{ int   used_cpath;
  int   i, type;
  char *spath, *tpath;
  FILE *test, **units;

  used_cpath = 0;
  test = fopen(source,"r");
  if (test == NULL)
    { if (*source != '/')
        test = fopen(MyCatenate(cpath,"/",source,""),"r");
      if (test == NULL)
        { EPRINTF(EPLACE,"%s: Could not find GDB %s\n",Prog_Name,source);
          EXIT (NULL);
        }
      source = Strdup(MyCatenate(cpath,"/",source,""),"Allocating expanded name");
      used_cpath = 1;
    }
  fclose(test);

  type  = Get_GDB_Paths(source,NULL,&spath,&tpath,0);
  if (type != IS_GDB)
    units = Create_GDB(gdb,spath,type,num_bps,NULL,0);
  else
    { Read_GDB(gdb,tpath);
      units = (FILE **) &(gdb->seqs);
      if (num_bps > 0)
        { if (gdb->seqs == NULL)
            { EPRINTF(EPLACE,"%s: GDB %s must have sequence data\n",Prog_Name,tpath);
              EXIT (NULL);
            }
          if (num_bps > 1)
            { units = malloc(sizeof(FILE *)*num_bps);
              if (units == NULL)
                { EPRINTF(EPLACE,"%s: Could not allocate units array for GDB's\n",Prog_Name);
                  EXIT(NULL);
                }
              units[0] = gdb->seqs;
              for (i = 1; i < num_bps; i++)
                { units[i] = fopen(gdb->seqpath,"r");
                  if (units[i] == NULL)
                    { EPRINTF(EPLACE,"%s: Cannot open another copye of GDB bps %s\n",
                                     Prog_Name,gdb->seqpath);
                      free(units);
                      EXIT(NULL);
                    }
                }
            }
        }
    }

  if (used_cpath)
    free(source);
  return (units);
}


/*******************************************************************************************
 *
 *  READ BUFFER ALLOCATION, LOAD, & LOAD_ALL
 *
 ********************************************************************************************/

// Allocate and return a buffer big enough for the largest contig in 'gdb'.
//   If cannot allocate memory then return NULL with an error message at EPLACE.
//   **NB** free(x-1) if x is the value returned as *prefix* and suffix 0(4)-bytes
//   are needed by the alignment algorithms.  

char *New_Contig_Buffer(GDB *gdb)
{ char *contig;

  contig = (char *) malloc(gdb->maxctg+4);
  if (contig == NULL)
    { EPRINTF(EPLACE,"%s: Cannot allocate a contig buffer\n",Prog_Name);
      EXIT(NULL);
    }
  return (contig+1);
}

// Return a pointer to the i'th contig in gdb in the format specified by stype.
// If buffer is NULL:
//     The GDB's sequences should be in memory in the same format and a pointer directly
//     to the sequence in this memory block is returned.
// If buffer != NULL:
//     The selected contig sequence is constructed in the memory block pointed at by
//     buffer in the format requested.  This includes reading the compressed sequence
//     from disk if this is where the GDB's sequences currently reside.
// If an error occurs then NULL is returned in interactive mode.

char *Get_Contig(GDB *gdb, int i, int stype, char *buffer)
{ char       *m = (char *) gdb->seqs;
  FILE       *b = (FILE *) gdb->seqs;
  GDB_CONTIG *c = gdb->contigs;
  int64       off;
  int         len, clen;

  (void) Print_Contig;

  if (m == NULL)
    { EPRINTF(EPLACE,"%s: GDB has no sequence data\n",Prog_Name);
      EXIT(NULL);
    }
  if (stype < COMPRESSED || stype > UPPER_CASE)
    { EPRINTF(EPLACE,"%s: Invalid load type %d (Get_Contig)\n",Prog_Name,stype);
      EXIT(NULL);
    }
  if (i < 0 || i >= gdb->ncontig)
    { EPRINTF(EPLACE,"%s: Index %d out of bounds (Get_Contig)\n",Prog_Name,i);
      EXIT(NULL);
    }

  off = c[i].boff;
  len = c[i].clen;

  if (buffer == NULL)
    { if (stype == gdb->seqstate)
        return (m + off);
      EPRINTF(EPLACE,"%s: GDB seq format does not match requested type (Get_Contig)\n",Prog_Name);
      EXIT(NULL);
    }

  if (gdb->seqstate != EXTERNAL)
    { if (gdb->seqstate == COMPRESSED)
        { memcpy(buffer,m + off,COMPRESSED_LEN(len));
          if (stype >= NUMERIC)
            Uncompress_Read(len,buffer);
          if (stype == LOWER_CASE)
            Lower_Read(buffer);
          else if (stype == UPPER_CASE)
            Upper_Read(buffer);
        }
      else
        { memcpy(buffer,m + off,len+1);
          if (stype == COMPRESSED)
            { if (gdb->seqstate > NUMERIC)
                Number_Read(buffer);
              Compress_Read(len,buffer);
            }
          else if (stype != gdb->seqstate)
            { if (stype == NUMERIC)
                Number_Read(buffer);
              else if (gdb->seqstate == NUMERIC)
                { if (stype == UPPER_CASE)
                    Upper_Read(buffer);
                  else
                    Lower_Read(buffer);
                }
              else
                Change_Read(buffer);
            }
        }

      if (stype == NUMERIC)
        buffer[-1] = 4;
      else if (stype != COMPRESSED)
        buffer[-1] = 0;

      return (buffer);
    }

  if (ftello(b) != off)
    fseeko(b,off,SEEK_SET);
  clen = COMPRESSED_LEN(len);
  if (clen > 0)
    { if (fread(buffer,clen,1,b) != 1)
        { EPRINTF(EPLACE,"%s: Failed read of GDB sequence file (Get_Contig)\n",Prog_Name);
          EXIT(NULL);
        }
    }

  if (stype == COMPRESSED)
    return (buffer);

  Uncompress_Read(len,buffer);
  if (stype == LOWER_CASE)
    { Lower_Read(buffer);
      buffer[-1] = '\0';
    }
  else if (stype == UPPER_CASE)
    { Upper_Read(buffer);
      buffer[-1] = '\0';
    }
  else
    buffer[-1] = 4;

  return (buffer);
}


char *Get_Contig_Piece(GDB *gdb, int i, int beg, int end, int stype, char *buffer)
{ FILE       *b  = (FILE *) gdb->seqs;
  char       *m  = (char *) gdb->seqs;
  GDB_CONTIG *c = gdb->contigs;
  int64      off;
  int        len, clen, bbeg;

  if (m == NULL)
    { EPRINTF(EPLACE,"%s: GDB has no sequence data (Get_Contig_Piece)\n",Prog_Name);
      EXIT(NULL);
    }
  if (stype < NUMERIC || stype > UPPER_CASE)
    { EPRINTF(EPLACE,"%s: Invalid load type %d (Get_Contig_Piece)\n",Prog_Name,stype);
      EXIT(NULL);
    }
  if (i < 0 || i >= gdb->ncontig)
    { EPRINTF(EPLACE,"%s: Index %d out of bounds (Get_Contig_Piece)\n",Prog_Name,i);
      EXIT(NULL);
    }


  if (beg < 0 || end > c[i].clen)
    { EPRINTF(EPLACE,"%s: Subrange %d,%d out of bounds (Get_Contig_Piece)\n",Prog_Name,beg,end);
      EXIT(NULL);
    }

  if (buffer == NULL)
    { if (stype == gdb->seqstate)
        return (m + c[i].boff + beg);
      EPRINTF(EPLACE,"%s: GDB seq format does not match requested type (Get_Contig_Piece)\n",
                     Prog_Name);
      EXIT(NULL);
    }
 
  bbeg = beg/4;
  off  = c[i].boff + bbeg;
  len  = end - beg;
  clen = ((end-1)/4+1) - bbeg;

  if (gdb->seqstate != EXTERNAL)
    { if (gdb->seqstate == COMPRESSED)
        { memcpy(buffer,m + off,clen);
          Uncompress_Read(len,buffer);
          buffer += beg%4;
          buffer[len] = 4;
          if (stype == LOWER_CASE)
            Lower_Read(buffer);
          else if (stype == UPPER_CASE)
            Upper_Read(buffer);
        }
      else
        { memcpy(buffer,m + off + (beg%4),len+1);
          buffer[len] = '\0';
          if (stype != gdb->seqstate)
            { if (stype == NUMERIC)
                Number_Read(buffer);
              else if (gdb->seqstate == NUMERIC)
                { buffer[len] = 4;
                  if (stype == UPPER_CASE)
                    Upper_Read(buffer);
                  else
                    Lower_Read(buffer);
                }
              else
                Change_Read(buffer);
            }
        }

      if (stype == NUMERIC)
        buffer[-1] = 4;
      else if (stype != COMPRESSED)
        buffer[-1] = 0;

      return (buffer);
    }

  if (ftello(b) != off)
    fseeko(b,off,SEEK_SET);
  if (clen > 0)
    { if (fread(buffer,clen,1,b) != 1)
        { EPRINTF(EPLACE,"%s: Failed read of GDB sequence file (Get_Contig_Piece)\n",Prog_Name);
          EXIT(NULL);
        }
    }
  Uncompress_Read(4*clen,buffer);
  buffer += beg%4;
  buffer[len] = 4;
  if (stype == LOWER_CASE)
    { Lower_Read(buffer);
      buffer[-1] = '\0';
    }
  else if (stype == UPPER_CASE)
    { Upper_Read(buffer);
      buffer[-1] = '\0';
    }
  else
    buffer[-1] = 4;

  return (buffer);
}
