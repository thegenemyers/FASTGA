/*******************************************************************************************
 *
 *  Utility for displaying the overlaps in a .las file in a variety of ways including
 *    a minimal listing of intervals, a cartoon, and a full out alignment.
 *
 *  Author:    Gene Myers
 *  Creation:  July 2013
 *  Last Mod:  Jan 2015
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"
#include "align.h"

static char *Usage = " <src1:db|dam> [ <src2:db|dam> ] <align:las>";


/*******************************************************************************************
 *
 *  Splay tree routines for ordered list: INSERT, DELETE, FIND, NEXT
 *
 ********************************************************************************************/

typedef struct vtx
  { struct vtx *L, *R;
    int64       V;
    int         score;
    int         link;
  } NODE;

#ifdef DEBUG_CHAIN

static void PRINT_LIST(NODE *v)
{ if (v == NULL)
    return;
  PRINT_LIST(v->L);
  printf(" %lld:%d:%d",v->V,v->score,v->link);
  PRINT_LIST(v->R);
}

#endif

#ifdef DEBUG_SPLAY

static void PRINT_TREE(NODE *v, int deep, NODE *space)
{ if (v == NULL)
    return;
  PRINT_TREE(v->R,deep+3,space);
  printf("%*s %lld:%d:%d (%ld)\n",deep,"",v->V,v->score,v->link,v-space);
  PRINT_TREE(v->L,deep+3,space);
}

#endif

static NODE *SPLAY(NODE *v, int64 x)    //  Assumes x is in the tree
{ NODE *u, *n;

  if (v == NULL || x == v->V)
    return (v);
  if (x < v->V)
    { u = v->L;
      if (x == u->V)
        { v->L = u->R;
          u->R = v;
          return (u);
        }
      if (x < u->V)
        { n = SPLAY(u->L,x);
          v->L = u->R;
          u->R = v;
          u->L = n->R;
          n->R = u;
        }
      else
        { n = SPLAY(u->R,x);
          v->L = n->R;
          u->R = n->L;
          n->L = u;
          n->R = v;
        }
    }
  else
    { u = v->R;
      if (x == u->V)
        { v->R = u->L;
          u->L = v;
          return (u);
        }
      if (x > u->V)
        { n = SPLAY(u->R,x);
          v->R = u->L;
          u->L = v;
          u->R = n->L;
          n->L = u;
        }
      else
        { n = SPLAY(u->L,x);
          v->R = n->L;
          u->L = n->R;
          n->R = u;
          n->L = v;
        }
    }
  return (n);
}

static NODE *FIND(NODE *v, int64 x)   //  Find v s.t. v->V <= x && x < v->next->V
{ NODE *u;

  if (v == NULL || v->V == x)
    return (v);
  if (x < v->V)
    return (FIND(v->L,x));
  else
    { u = FIND(v->R,x);
      if (u == NULL)
        return (v);
      else
        return (u);
    }
}

static NODE *NEXT(NODE *v, NODE *t, NODE *w)
{ if (v == NULL || t->V == v->V)
    { if (v->R != NULL)
        { w = v->R;
          while (w->L != NULL)
            w = w->L;
        }
      return (w);
    }
  if (t->V < v->V)
    return (NEXT(v->L,t,v));
  else
    return (NEXT(v->R,t,w));
}

static NODE *JOIN(NODE *v, NODE *w)
{ NODE *p;

  if (v == NULL)
    return (w);
  for (p = v; p->R != NULL; p = p->R)
    ;
  v = SPLAY(v,p->V);
  v->R = w;
  return (v);
}

static NODE *INSERT(NODE *v, NODE *new)
{ NODE *u, *p;

  if (v == NULL)
    return (new);
  u = FIND(v,new->V);
  if (u != NULL && u->R == NULL)
    u->R = new;
  else
    { if (u == NULL)
        p = v;
      else  // u->R == NULL
        p = u->R;
      while (p->L != NULL)
        p = p->L;
      p->L = new;
    }
  return (SPLAY(v,new->V));
}

static NODE *DELETE(NODE *v, NODE *old)
{ NODE *u, *w;

  u = FIND(v,old->V);
  if (u == NULL || u->V != old->V)
    return (NULL);
  v = SPLAY(v,old->V);
  w = JOIN(v->L,v->R);
  return (w);
}


  v = NULL;
  new = Node(0,0);
  for (hits) do
    { w = FIND(v,hit->bbpos);
      lnk = w->lnk;
      scr = w->score;
      if (lnk->aepos > hit->abpos)
        { if (lnk->aepos < hit->aepos)
            scr += (hit->aepos - lnk->aepos);
        }
      else
        scr += (hit->aepos-hit->abpos);
      

int main(int argc, char *argv[])
{ DAZZ_DB   _db1, *db1 = &_db1;
  DAZZ_DB   _db2, *db2 = &_db2;
  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;

  FILE   *input;
  int     sameDB, ISTWO;
  int64   novl;
  int     tspace, tbytes, small;

  //  Process options

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("Hitter")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 3 && argc != 4)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }


  //  Open trimmed DB or DB pair

  { int   status;
    char *pwd, *root;
    FILE *input;
    struct stat stat1, stat2;

    ISTWO  = 0;
    status = Open_DB(argv[1],db1);
    if (status < 0)
      exit (1);
    if (db1->part > 0)
      { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    sameDB = 1;
    if (argc > 3)
      { pwd   = PathTo(argv[3]);
        root  = Root(argv[3],".las");
        if ((input = fopen(Catenate(pwd,"/",root,".las"),"r")) != NULL)
          { ISTWO = 1;
            fclose(input);
            status = Open_DB(argv[2],db2);
            if (status < 0)
              exit (1);
            if (db2->part > 0)
              { fprintf(stderr,"%s: Cannot be called on a block: %s\n",Prog_Name,argv[2]);
                exit (1);
              }
            stat(Catenate(db1->path,"","",".idx"),&stat1);
            stat(Catenate(db2->path,"","",".idx"),&stat2);
            if (stat1.st_ino != stat2.st_ino)
              sameDB = 0;
            Trim_DB(db2);
          }
        else
          db2 = db1;
        free(root);
        free(pwd);
      }
    else
      db2 = db1;
    Trim_DB(db1);
  }

  //  Initiate file reading and read (novl, tspace) header
  
  { char  *over, *pwd, *root;

    pwd   = PathTo(argv[2+ISTWO]);
    root  = Root(argv[2+ISTWO],".las");
    over  = Catenate(pwd,"/",root,".las");
    input = Fopen(over,"r");
    if (input == NULL)
      exit (1);

    if (fread(&novl,sizeof(int64),1,input) != 1)
      SYSTEM_READ_ERROR
    if (fread(&tspace,sizeof(int),1,input) != 1)
      SYSTEM_READ_ERROR
    if (tspace < 0)
      { fprintf(stderr,"%s: Garbage .las file, trace spacing < 0 !\n",Prog_Name);
        exit (1);
      }
    if (tspace <= TRACE_XOVR && tspace != 0)
      { small  = 1;
        tbytes = sizeof(uint8);
      }
    else
      { small  = 0;
        tbytes = sizeof(uint16);
      }

    free(pwd);
    free(root);
  }

  //  Read the file and display selected records
  
  { int        j, i;
    uint16    *trace;
    int        tmax;
    int        aread;
    int64     *wgt;
    int       *cnt;

    wgt = (int64 *) Malloc(sizeof(int64)*db2->treads,"Allocating weight vector");
    cnt = (int *) Malloc(sizeof(int)*db2->treads,"Allocating weight vector");

    tmax  = 1000;
    trace = (uint16 *) Malloc(sizeof(uint16)*tmax,"Allocating trace vector");
    if (trace == NULL)
      exit (1);

    aln->path = &(ovl->path);
    ovl->path.trace = (void *) trace;

    //  For each alignment do

    aread = -1;
    for (j = 0; j < novl; j++)
      { Read_Overlap(input,ovl);
        if (ovl->path.tlen > tmax)
          { tmax = ((int) 1.2*ovl->path.tlen) + 100;
            trace = (uint16 *) Realloc(trace,sizeof(uint16)*tmax,"Allocating trace vector");
            if (trace == NULL)
              exit (1);
            ovl->path.trace = (void *) trace;
          }
        Read_Trace(input,ovl,tbytes);

        if (small)
          Decompress_TraceTo16(ovl);

        if (aread != ovl->aread)
          { if (aread >= 0)
              { printf("\nA-contig %d:\n",aread);
                for (i = 0; i < db2->treads; i++)
                  if (wgt[i] > 0)
                    printf("  %3d: %8lld %5d\n",i,wgt[i],cnt[i]);
              }
            for (i = 0; i < db2->treads; i++)
              { wgt[i] = 0;
                cnt[i] = 0;
              }
            aread = ovl->aread;
          }

printf(" %d: %d   %lld\n",ovl->bread,ovl->path.aepos-ovl->path.abpos,wgt[18]);
        wgt[ovl->bread] += ovl->path.aepos - ovl->path.abpos;
        cnt[ovl->bread] += 1;
      }

    printf("\nA-contig %d:\n",aread);
    for (i = 0; i < db2->treads; i++)
      if (wgt[i] > 0)
        printf("  %3d: %8lld %5d\n",i,wgt[i],cnt[i]);

    free(trace);
  }

  Close_DB(db1);
  if (ISTWO)
    Close_DB(db2);

  exit (0);
}
