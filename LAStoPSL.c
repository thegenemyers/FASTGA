/*******************************************************************************************
 *
 *  Utility for displaying the overlaps in a .las file in a variety of ways including
 *    a minimal listing of intervals, a cartoon, and a full out alignment.
 *
 *  Author:    Gene Myers
 *  Creation:  Oct 2023
 *  Last Mod:  Dec 2023
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "align.h"

#undef DEBUG_THREADS

static char *Usage = " [-T<int(8)>} <align:las>";

static int  NTHREADS;  // -T
static int  ISTWO;     // one gdb or two?

static int TSPACE;   // Trace spacing

static int  AMAX,  BMAX;   //  Max sequence length in A and B
static int *AMAP, *BMAP;   //  Map contig to its scaffold
static int *ALEN, *BLEN;   //  Map contig to its scaffold length

static int PTR_SIZE = sizeof(void *);
static int OVL_SIZE = sizeof(Overlap) - sizeof(void *);


//   ROUTINES TO FIND AN OVERLAP RECORD WITHIN A .LAS FILE

#define SEEK_BLOCK  0x100000

  //  Is o possibly a valid Overlap record?

static int plausible_record(Overlap *o)
{ int l1, l2;

  if (o->flags < 0 || o->flags > 2)
   return (0);
  if (o->path.abpos < 0 || o->path.abpos > o->path.aepos || o->path.aepos > AMAX)
    return (0);
  if (2*(o->path.aepos/TSPACE - o->path.abpos/TSPACE + 1) != o->path.tlen)
    return (0);
  if (o->path.bbpos < 0 || o->path.bbpos > o->path.bepos || o->path.bepos > BMAX)
    return (0);
  l1 = o->path.aepos - o->path.abpos;
  l2 = o->path.bepos - o->path.bbpos;
  if (o->path.diffs < 0 || abs(l1-l2) > o->path.diffs)
    return (0);
  if (o->aread < 0 || o->bread < 0)
    return (0);
  return (o->path.tlen);
}

  //  Find first location after offset in input of an Overlap record

static int64 find_ovl_boundary(int64 offset, FILE *input, void *buffer)
{ void *bend, *off;
  int   y, len;

  bend = buffer + SEEK_BLOCK - (OVL_SIZE+PTR_SIZE);

  fseek(input,offset,SEEK_SET);
  fread(buffer,SEEK_BLOCK,1,input);
  for (y = 0; 1; y += 2)
    { off = (buffer+y) - PTR_SIZE;    //  3 linked plausibles very firm
      if (off >= bend)
        return (0);
      len = plausible_record(off);
      if (len == 0)
        continue;
      off += OVL_SIZE + len;
      if (off >= bend)
        return (0);
      len = plausible_record(off);
      if (len == 0)
        continue;
      off += OVL_SIZE + len;
      if (off >= bend)
        return (0);
      len = plausible_record(off);
      if (len == 0)
        continue;
      break;
    }

  return (offset+y);
}


//  THREAD ROUTINE TO GENERATE A SECTION OF THE DESIRED .PAF FILE

typedef struct
  { int64   beg;
    int64   end;
    DAZZ_DB db1;
    DAZZ_DB db2;
    FILE   *in;
    FILE   *out;
  } Packet;

void *gen_paf(void *args)
{ Packet *parm = (Packet *) args;

  int64    beg = parm->beg;
  int64    tot = parm->end - beg;
  DAZZ_DB *db1 = &(parm->db1);
  DAZZ_DB *db2 = &(parm->db2);
  FILE    *in  = parm->in;
  FILE    *out = parm->out;

  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;

  FILE   *ahdr, *bhdr;
  
  int64      bytes;
  uint16    *trace;
  int        tmax;
  DAZZ_READ *reads1, *reads2;
  int        aread, bread;
  int        bmin, bmax;
  char      *aseq, *bseq, *bact;
  int        aoff, boff;
  Path      *path;
  Work_Data *work;
  char       aheader[MAX_NAME], bheader[MAX_NAME];

  work = New_Work_Data();
  aseq = New_Read_Buffer(db1);
  bseq = New_Read_Buffer(db2);
  if (aseq == NULL || bseq == NULL)
    exit (1);

  reads1 = db1->reads;
  reads2 = db2->reads;

  tmax  = 1000;
  trace = (uint16 *) Malloc(sizeof(uint16)*tmax,"Allocating trace vector");
  if (trace == NULL)
    exit (1);

  aln->path = path = &(ovl->path);
  aln->aseq = aseq;
  path->trace = (void *) trace;

  ahdr = fopen(Catenate(db1->path,".hdr","",""),"r");
  if (ahdr == NULL)
    { fprintf(stderr,"%s: Cannot open %s.hdr file for reading\n",Prog_Name,db1->path);
      exit (1);
    }
  if (ISTWO)
    { bhdr = fopen(Catenate(db2->path,".hdr","",""),"r");
      if (bhdr == NULL)
        { fprintf(stderr,"%s: Cannot open %s.hdr file for reading\n",Prog_Name,db2->path);
          exit (1);
        }
    }
  else
    bhdr = ahdr;

  //  For each alignment do

  fseek(in,beg,SEEK_SET);

  aread = -1;
  bytes = 0;
  while (bytes < tot)
    { Read_Overlap(in,ovl);
      if (ovl->path.tlen > tmax)
        { tmax = ((int) 1.2*ovl->path.tlen) + 100;
          trace = (uint16 *) Realloc(trace,sizeof(uint16)*tmax,"Allocating trace vector");
          if (trace == NULL)
            exit (1);
        }
      path->trace = (void *) trace;
      Read_Trace(in,ovl,1);

      bytes += OVL_SIZE + path->tlen;

      Decompress_TraceTo16(ovl);

      if (aread != ovl->aread)
        { aread = ovl->aread;
          if (Load_Read(db1,aread,aseq,0))
            exit (1);
          aln->alen = reads1[aread].rlen;
          aoff      = reads1[aread].fpulse;
          fseeko(ahdr,reads1[aread].coff,SEEK_SET);
          fgets(aheader,MAX_NAME,ahdr);
          aheader[strlen(aheader)-1] = '\0';
        }
      bread = ovl->bread;
      aln->blen  = reads2[bread].rlen;
      boff       = reads2[bread].fpulse;
      fseeko(bhdr,reads2[bread].coff,SEEK_SET);
      fgets(bheader,MAX_NAME,bhdr);
      bheader[strlen(bheader)-1] = '\0';
      aln->flags = ovl->flags;

      if (COMP(aln->flags))
        { bmin = (aln->blen-path->bepos);
          if (bmin < 0) bmin = 0;
          bmax = (aln->blen-path->bbpos);
          if (bmax > aln->blen) bmax = aln->blen;
        }
      else
        { bmin = path->bbpos;
          if (bmin < 0) bmin = 0;
          bmax = path->bepos;
          if (bmax > aln->blen) bmax = aln->blen;
        }

      bact = Load_Subread(db2,bread,bmin,bmax,bseq,0);
      if (COMP(aln->flags))
        { Complement_Seq(bact,bmax-bmin);
          aln->bseq = bact - (aln->blen-bmax);
        }
      else
        aln->bseq = bact - bmin; 

      Compute_Trace_PTS(aln,work,TSPACE,GREEDIEST);
      if (Gap_Improver(aln,work))
        exit (1);

      { int     i, j, x, p, q;
        int    *t, T;
        int     M, N;
        int     I, D, S, X;
        int     IB, DB;
        int     bcnt, bmat;

        t = (int *) path->trace;
        T = path->tlen;

        M = path->aepos - path->abpos;
        N = path->bepos - path->bbpos;
        I = D = 0;
        IB = DB = 0;
        p = 0; 
        for (x = 0; x < T; x++)
          { q = p;
            if ((p = t[x]) < 0)
              { I += 1;
                if (p != q)
                  IB += 1;
              }
            else
              { D += 1;
                if (p != q)
                  DB += 1;
              }
          }
        S = path->diffs - (I+D);
        X = (M+N - (I+D+2*S))/2;

        fprintf(out,"%d\t%d\t0\t0\t%d\t%d\t%d\t%d\t%c",X,S,IB,I,DB,D,COMP(ovl->flags)?'-':'+');
        fprintf(out,"\t%s\t%d\t%d\t%d",aheader,ALEN[aread],aoff+path->abpos,aoff+path->aepos);
        fprintf(out,"\t%s\t%d\t%d\t%d",bheader,BLEN[bread],boff+path->bbpos,boff+path->bepos);

        bcnt = 0;
        i = path->abpos+1;
        j = path->bbpos+1;
        for (x = 0; x < T; x++)
          { if ((p = t[x]) < 0)
              { bmat = -(p+i);
                i += bmat;
                j += bmat+1;
              }
            else
              { bmat = p-j;
                i += bmat+1;
                j += bmat;
              }
            if (bmat > 0)
              bcnt += 1;
          }
        bmat = (path->aepos - i)+1;
        if (bmat > 0)
          bcnt += 1;
        fprintf(out,"\t%d\t",bcnt);

        i = path->abpos+1;
        j = path->bbpos+1;
        for (x = 0; x < T; x++)
          { if ((p = t[x]) < 0)
              { bmat = -(p+i);
                i += bmat;
                j += bmat+1;
              }
            else
              { bmat = p-j;
                i += bmat+1;
                j += bmat;
              }
            if (bmat > 0)
              fprintf(out,"%d,",bmat);
          }
        bmat = (path->aepos - i)+1;
        if (bmat > 0)
          fprintf(out,"%d,",bmat);
        fprintf(out,"\t");

        i = path->abpos+1;
        j = path->bbpos+1;
        for (x = 0; x < T; x++)
          { if ((p = t[x]) < 0)
              { bmat = -(p+i);
                if (bmat > 0)
                  fprintf(out,"%d,",i);
                i += bmat;
                j += bmat+1;
              }
            else
              { bmat = p-j;
                if (bmat > 0)
                  fprintf(out,"%d,",i);
                i += bmat+1;
                j += bmat;
              }
          }
        bmat = (path->aepos - i)+1;
        if (bmat > 0)
          fprintf(out,"%d,",i);
        fprintf(out,"\t");

        i = path->abpos+1;
        j = path->bbpos+1;
        for (x = 0; x < T; x++)
          { if ((p = t[x]) < 0)
              { bmat = -(p+i);
                if (bmat > 0)
                  fprintf(out,"%d,",j);
                i += bmat;
                j += bmat+1;
              }
            else
              { bmat = p-j;
                if (bmat > 0)
                  fprintf(out,"%d,",j);
                i += bmat+1;
                j += bmat;
              }
          }
        bmat = (path->aepos - i)+1;
        if (bmat > 0)
          fprintf(out,"%d,",j);
        fprintf(out,"\n");
      }
    }

  if (ISTWO)
    fclose(bhdr);
  fclose(ahdr);
  fclose(in);

  free(trace);
  free(bseq-1);
  free(aseq-1);
  Free_Work_Data(work);

  return (NULL);
}


int main(int argc, char *argv[])
{ Packet    *parm;
  DAZZ_DB   _db1, *db1 = &_db1;
  DAZZ_DB   _db2, *db2 = &_db2;
  char      *db1_name;
  char      *db2_name;
  FILE      *input;
  char      *iname;

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("LAStoPSL")

    NTHREADS = 8;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
        exit (1);
      }

    parm = Malloc(sizeof(Packet)*NTHREADS,"Allocating thread records");
    if (parm == NULL)
      exit (1);
  }

  //  Initiate .las file reading and read header information

  { char       *pwd, *root, *cpath;
    FILE       *test;
    int64       novl;
    int         nlen;

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".las");
    iname = Strdup(Catenate(pwd,"/",root,".las"),"Allocating input name");
    input = Fopen(iname,"r");
    if (input == NULL)
      exit (1);

    if (fread(&novl,sizeof(int64),1,input) != 1)
      SYSTEM_READ_ERROR
    if (fread(&TSPACE,sizeof(int),1,input) != 1)
      SYSTEM_READ_ERROR
    if (TSPACE < 0)
      { fprintf(stderr,"%s: Garbage .las file, trace spacing < 0 !\n",Prog_Name);
        exit (1);
      }

    free(pwd);
    free(root);

    if (fread(&nlen,sizeof(int),1,input) != 1)
      SYSTEM_READ_ERROR
    db1_name = Malloc(nlen+1,"Allocating name 1\n");
    if (db1_name == NULL)
      exit (1);
    if (fread(db1_name,nlen,1,input) != 1)
      SYSTEM_READ_ERROR
    db1_name[nlen] = '\0';

    if (fread(&nlen,sizeof(int),1,input) != 1)
      SYSTEM_READ_ERROR
    if (nlen == 0)
      db2_name = NULL;
    else
      { db2_name = Malloc(nlen+1,"Allocating name 2\n");
        if (db2_name == NULL)
          exit (1);
        if (fread(db2_name,nlen,1,input) != 1)
          SYSTEM_READ_ERROR
        db2_name[nlen] = '\0';
      }

    if (fread(&nlen,sizeof(int),1,input) != 1)
      SYSTEM_READ_ERROR
    cpath = Malloc(nlen+1,"Allocating name 1\n");
    if (cpath == NULL)
      exit (1);
    if (fread(cpath,nlen,1,input) != 1)
      SYSTEM_READ_ERROR
    cpath[nlen] = '\0';

    test = fopen(db1_name,"r");
    if (test == NULL)
      { test = fopen(Catenate(cpath,db1_name,"",""),"r");
        if (test == NULL)
          { fprintf(stderr,"%s: Could not find .gdb %s\n",Prog_Name,db1_name);
            exit (1);

          }
        pwd = Strdup(Catenate(cpath,db1_name,"",""),"Allocating expanded name");
        free(db1_name);
        db1_name = pwd;
      }
    fclose(test);

    if (db2_name != NULL)
      { test = fopen(db2_name,"r");
        if (test == NULL)
          { test = fopen(Catenate(cpath,db2_name,"",""),"r");
            if (test == NULL)
              { fprintf(stderr,"%s: Could not find .gdb %s\n",Prog_Name,db2_name);
                exit (1);
              }
            pwd = Strdup(Catenate(cpath,db2_name,"",""),"Allocating expanded name");
            free(db2_name);
            db2_name = pwd;
          }
        fclose(test);
      }

    free(cpath);
  }

  //  Open DB or DB pair

  { int   s, r, gdb;

    ISTWO  = 0;
    gdb    = Open_DB(db1_name,db1);
    if (gdb < 0)
      exit (1);

    if (db2_name != NULL)
      { gdb = Open_DB(db2_name,db2);
        if (gdb < 0)
          exit (1);
        ISTWO = 1;
      }
    else
      db2  = db1;

    //  Build contig to scaffold maps in global vars

    AMAX = db1->maxlen;
    BMAX = db2->maxlen;

    if (ISTWO)
      AMAP = (int *) Malloc(sizeof(int)*2*(db1->treads+db2->treads),"Allocating scaffold map");
    else
      AMAP = (int *) Malloc(sizeof(int)*2*db1->treads,"Allocating scaffold map");
    if (AMAP == NULL)
      exit (1);
    ALEN = AMAP + db1->treads;

    s = -1;
    for (r = 0; r < db1->treads; r++)
      { if (db1->reads[r].origin == 0)
          s += 1;
        ALEN[s] = db1->reads[r].fpulse + db1->reads[r].rlen;
        AMAP[r] = s;
      }
    for (r = db1->treads-1; r >= 0; r--)
      ALEN[r] = ALEN[AMAP[r]];

    if (ISTWO)
      { BMAP = ALEN + db1->treads;
        BLEN = BMAP + db2->treads;

        s = -1;
        for (r = 0; r < db2->treads; r++)
          { if (db2->reads[r].origin == 0)
              s += 1;
            BLEN[s] = db2->reads[r].fpulse + db2->reads[r].rlen;
            BMAP[r] = s;
          }
        for (r = db2->treads-1; r >= 0; r--)
          BLEN[r] = BLEN[BMAP[r]];
      }
    else
      { BMAP = AMAP;
        BLEN = ALEN;
      }
  }

  //  Divide .las into NTHREADS parts

  { int   i;
    int64 p, x, ioff, ipar, size;
    char *buffer;
    struct stat info;

    stat(iname,&info);
    size = info.st_size;
    ioff = ftello(input);
    size -= ioff;
    ipar = (ioff % 2);

    buffer = Malloc(SEEK_BLOCK,"Allocating seek block");
    parm[0].beg = ioff;
    for (i = 1; i < NTHREADS; i++)
      { p = (size*i)/NTHREADS + ioff;
        if (p % 2 != ipar)
          p += 1;
        x = find_ovl_boundary(p,input,buffer);
        if (x == 0)
          { fprintf(stderr,"%s: Could not partition .las file ??\n",Prog_Name);
            exit (1);
          }
        parm[i].beg = parm[i-1].end = x;
      }
    parm[NTHREADS-1].end = size;

    free(buffer);

    fclose(input);
  }

  //   Use NTHREADS to produce .paf for each part and then cat the parts to stdout

  { int   p, x;
    char *oprefix;
    char *buffer;
#ifndef DEBUG_THREADS
    pthread_t threads[NTHREADS];
#endif

    //  Setup thread packets

    oprefix = Strdup(Numbered_Suffix("_psl.",getpid(),"."),"Allocating temp file prefix");

    for (p = 0; p < NTHREADS; p++)
      { parm[p].db1   = *db1;
        parm[p].db2   = *db2;
        if (p > 0)
          { parm[p].db1.bases = fopen(Catenate(db1->path,"","",".bps"),"r");
            if (parm[p].db1.bases == NULL)
              { fprintf(stderr,"%s: Cannot open another copy of DB\n",Prog_Name);
                exit (1);
              }
            parm[p].db2.bases = fopen(Catenate(db2->path,"","",".bps"),"r");
            if (parm[p].db2.bases == NULL)
              { fprintf(stderr,"%s: Cannot open another copy of DB\n",Prog_Name);
                exit (1);
              }
          }
        parm[p].in  = fopen(iname,"r");
        if (parm[p].in == NULL)
          { fprintf(stderr,"%s: Cannot open %s for reading\n",Prog_Name,iname);
            exit (1);
          }
        parm[p].out = fopen(Numbered_Suffix(oprefix,p,".psl"),"w+");
        if (parm[p].out == NULL)
          { fprintf(stderr,"%s: Cannot open %s.%d.paf for reading & writing\n",
                           Prog_Name,oprefix,p);
            exit (1);
          }
        unlink(Numbered_Suffix(oprefix,p,".psl"));
      }

    //  Launch and then gather threads

#ifdef DEBUG_THREADS
    for (p = 0; p < NTHREADS; p++)
      gen_paf(parm+p);
#else
    for (p = 1; p < NTHREADS; p++)
      pthread_create(threads+p,NULL,gen_paf,parm+p);
    gen_paf(parm);
    for (p = 1; p < NTHREADS; p++)
      pthread_join(threads[p],NULL);
#endif

    for (p = 1; p < NTHREADS; p++)
      { fclose(parm[p].db1.bases);
        fclose(parm[p].db2.bases);
      }

    //  Concatenate thread generated paf parts to stdout

#define PUSH_BLOCK 0x100000

    buffer = Malloc(PUSH_BLOCK,"Allocating seek block");

    for (p = 0; p < NTHREADS; p++)
      { FILE *o = parm[p].out;

        rewind(o);
        x = PUSH_BLOCK;
        while (x >= PUSH_BLOCK)
          { x = fread(buffer,1,PUSH_BLOCK,o);
            if (x > 0)
              fwrite(buffer,x,1,stdout);
          }
        fclose(o);
      }

    free(buffer);
  }

  free(AMAP);
  Close_DB(db1);
  if (ISTWO)
    Close_DB(db2);
  free(iname);

  exit (0);
}
