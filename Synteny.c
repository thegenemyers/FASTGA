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

#undef DEBUG_CHAIN

#define ALIGN_OVERLAP 20


/*******************************************************************************************
 *
 *  Splay tree routines for ordered list: INSERT, DELETE, FIND, NEXT
 *
 ********************************************************************************************/

typedef struct chain
  { struct chain *next;
    struct chain *link;
    struct chain *L, *R;
    int           bepos;
    int           score;
    int           clen;
    int           mark;
    int           dead;
    Overlap       ovl;
  } CHAIN;

typedef struct order
  { int event;
    int which;
  } ORDER;

static void PRINT_LIST(CHAIN *v)
{ if (v == NULL)
    return;
  PRINT_LIST(v->L);
  printf(" %d:%d",v->bepos,v->score);
  PRINT_LIST(v->R);
}

#ifdef DEBUG_SPLAY

static void PRINT_TREE(CHAIN *v, int deep, CHAIN *space)
{ if (v == NULL)
    return;
  PRINT_TREE(v->R,deep+3,space);
  printf("%*s %lld:%d:%d (%ld)\n",deep,"",v->score,v->link,v-space);
  PRINT_TREE(v->L,deep+3,space);
}

#endif

static CHAIN *SPLAY(CHAIN *v, int64 x)    //  Assumes x is in the tree
{ CHAIN *u, *n;

  if (v == NULL || x == v->bepos)
    return (v);
  if (x < v->bepos)
    { u = v->L;
      if (x == u->bepos)
        { v->L = u->R;
          u->R = v;
          return (u);
        }
      if (x < u->bepos)
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
      if (x == u->bepos)
        { v->R = u->L;
          u->L = v;
          return (u);
        }
      if (x > u->bepos)
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

static CHAIN *FIND(CHAIN *v, int x)   //  Find v s.t. v->bepos <= x && x < v->next->bepos
{ CHAIN *u;

  if (v == NULL || v->bepos == x)
    return (v);
  if (x < v->bepos)
    return (FIND(v->L,x));
  else
    { u = FIND(v->R,x);
      if (u == NULL)
        return (v);
      else
        return (u);
    }
}

static CHAIN *NEXT(CHAIN *v, int x, CHAIN *w)
{ if (v == NULL || x == v->bepos)
    { if (v->R != NULL)
        { w = v->R;
          while (w->L != NULL)
            w = w->L;
        }
      return (w);
    }
  if (x < v->bepos)
    return (NEXT(v->L,x,v));
  else
    return (NEXT(v->R,x,w));
}

static CHAIN *JOIN(CHAIN *v, CHAIN *w)
{ CHAIN *p;

  if (v == NULL)
    return (w);
  for (p = v; p->R != NULL; p = p->R)
    ;
  v = SPLAY(v,p->bepos);
  v->R = w;
  return (v);
}

static CHAIN *INSERT(CHAIN *v, CHAIN *new)
{ CHAIN *u, *p;

  if (v == NULL)
    return (new);
  u = FIND(v,new->bepos);
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
  return (SPLAY(v,new->bepos));
}

static CHAIN *DELETE(CHAIN *v, CHAIN *old)
{ CHAIN *u, *w;

  u = FIND(v,old->bepos);
  if (u == NULL || u->bepos != old->bepos)
    return (NULL);
  v = SPLAY(v,old->bepos);
  w = JOIN(v->L,v->R);
  return (w);
}

int EORDER(const void *l, const void *r)
{ ORDER *x = (ORDER *) l;
  ORDER *y = (ORDER *) r;
  int xm, ym;

  xm = abs(x->event);
  ym = abs(y->event);
  if (xm < ym)
    return (-1);
  else if (xm > ym)
    return (1);
  else
    { if (x->event < y->event) 
        return (-1);
      else if (x->event > y->event)
        return (1);
      else
        return (0);
    }
}

void analyze(int ascaf, CHAIN *chain, CHAIN **blist, int *bstack, int btop, ORDER *order)
{ int    i, b, acnt, bscaf;
  CHAIN *prof, *e, *v, *w;

#ifndef DEBUG_CHAIN
  (void) ascaf;
  (void) PRINT_LIST;
#endif

  for (b = 0; b < btop; b++)
    { bscaf = bstack[b];

      acnt = 0;
      for (w = blist[bscaf]; w != NULL; w = w->next)
        { order[2*acnt].event   = -w->ovl.path.abpos;
          order[2*acnt].which   = w-chain;
          order[2*acnt+1].event = w->ovl.path.aepos - ALIGN_OVERLAP;
          order[2*acnt+1].which = w-chain;
          acnt += 1;
        }

      qsort(order,2*acnt,sizeof(ORDER),EORDER);

#ifdef DEBUG_CHAIN
      printf("\nScaf %d vs %d (%d)\n",ascaf,bscaf,acnt);
#endif

      prof = NULL;
      for (i = 0; i < 2*acnt; i++)
        { e = chain + order[i].which;
          if (order[i].event <= 0)
            { v = FIND(prof,e->ovl.path.bbpos);
              e->link = v;
              if (v == NULL)
                { e->score = (e->ovl.path.aepos - e->ovl.path.abpos);
                  e->clen  = 1;
#ifdef DEBUG_CHAIN
                  v = chain;
#endif
                }
              else
                { e->score = v->score + (e->ovl.path.aepos - e->ovl.path.abpos);
                  e->clen  = v->clen + 1;
                }
              e->bepos = e->ovl.path.bepos - ALIGN_OVERLAP;
              e->L     = NULL;
              e->R     = NULL;
#ifdef DEBUG_CHAIN
              printf("  A %8d: %6d   (sc %6d  [%8d,%8d] -> %ld)\n",
                     -order[i].event,order[i].which,e->score,e->ovl.path.bbpos,e->bepos,v-chain);
              fflush(stdout);
#endif
            }
          else
            { v = FIND(prof,e->bepos);
              if (v == NULL)
                prof = INSERT(prof,e);
              else if (v->score <= e->score)
                { w = NEXT(prof,v->bepos,NULL);
                  if (v->bepos >= e->bepos)
                    prof = DELETE(prof,v);
                  while (w != NULL && w->score <= e->score)
                    { v = w;
                      w = NEXT(prof,w->bepos,NULL);
                      prof = DELETE(prof,v);
                    }
                  prof = INSERT(prof,e);
                }
#ifdef DEBUG_CHAIN
              printf("  D %8d: %6d\n",order[i].event,order[i].which);
              fflush(stdout);
              PRINT_LIST(prof);
              printf("\n");
              fflush(stdout);
#endif
            }
        }
    }
}

void collisions(int acnt, CHAIN *chain)
{ int i, j;
  int abeg, aend;
  int cbeg, cend, cov;
  int rmark, lmark;

  for (i = 0; i < acnt; i++)
    chain[i].mark = chain[i].ovl.path.abpos;

  for (i = 0; i < acnt; i++)
    { aend = chain[i].ovl.path.aepos;
      abeg = chain[i].ovl.path.abpos;
      for (j = i+1; j < acnt; j++)
        { if (chain[j].ovl.path.abpos >= aend)
            break;
          if (chain[j].mark < aend)
            chain[j].mark = aend;
	}
      if (i+1 < acnt && chain[i+1].ovl.path.abpos < aend)
        rmark = chain[i+1].ovl.path.abpos;
      else
        rmark = aend;
      lmark = chain[i].mark;
      if (lmark > aend)
        lmark = aend;
      cov = lmark-abeg;
#ifdef DEBUG_REPEAT
      printf("  %d: %d..%d  <%d>  -> %d\n",i,abeg,aend,lmark,cov);
#endif
      for (j = i+1; j < acnt; j++)
        { cbeg = chain[j].ovl.path.abpos;
          cend = chain[j].ovl.path.aepos;
          if (cbeg >= aend)
            break;
          if (cend <= lmark)
            continue;
          if (cend > aend)
            cend = aend;
          if (cbeg < lmark)
            cov += cend - lmark;
          else
            cov += cend - cbeg;
#ifdef DEBUG_REPEAT
          printf("            %d: %d\n",j,cov);
#endif
          lmark = cend;
        }
      chain[i].mark = (100*cov)/(aend-abeg);
#ifdef DEBUG_REPEAT
      fflush(stdout);
#endif
    }
}

int main(int argc, char *argv[])
{ DAZZ_DB   _db1, *db1 = &_db1;
  DAZZ_DB   _db2, *db2 = &_db2;
  Overlap   _ovl, *ovl = &_ovl;
  int       amax, tmax;
  int       bscaf, *scaffold, *invscaff;

  FILE   *input;
  int     sameDB, ISTWO;
  int64   novl;
  int     tspace, tbytes, small;

  //  Process options

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("Synteny")

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

  //  Read the file to get maximums
  
  { int j;
    int aread, acnt;

    amax  = 0;
    tmax  = 0;
    acnt  = 0;
    aread = -1;
    for (j = 0; j < novl; j++)
      { Read_Overlap(input,ovl);
        fseek(input,tbytes*ovl->path.tlen,SEEK_CUR);
        if (ovl->path.tlen > tmax)
          tmax = ovl->path.tlen;
      
        if (aread != ovl->aread && db1->reads[ovl->aread].fpulse == 0)
          { if (acnt > amax)
              amax = acnt;
            acnt = 0;
          }
        aread = ovl->aread;
        acnt += 1;
      }
    if (acnt > amax)
      amax = acnt;

    printf("There are %d alignments for largest scaffold, longest trace is %d\n",amax,tmax);
    fflush(stdout);
  }

  //  Build scaffold map

  { int r;

    scaffold = (int *) Malloc(sizeof(int)*(2*db2->treads+1),"Allocating scaffold map");
    if (scaffold == NULL)
      exit (1);
    invscaff = scaffold + db2->treads;

    bscaf = -1;
    for (r = 0; r < db2->treads; r++)
      { if (db2->reads[r].fpulse == 0)
          { bscaf += 1;
            invscaff[bscaf] = r;
          }
        scaffold[r] = bscaf;
      }
    bscaf += 1;
    invscaff[bscaf] = db2->treads;

    printf("There are %d scaffolds with %d contigs\n",bscaf,db2->treads);
    fflush(stdout);
  }

  //  Read the file and chain
  
  { int    j, b;
    int    aread, acnt, ascaf;
    int    bread, btop;
    int    apulse, bpulse;
    CHAIN *chain, **blist;
    int   *bstack;
    CHAIN *n, *w;
    ORDER *order;

    rewind(input);
    fread(&novl,sizeof(int64),1,input);
    fread(&tspace,sizeof(int),1,input);

    chain  = (CHAIN *) Malloc(sizeof(CHAIN)*amax,"Allocating alignment space");
    blist  = (CHAIN **) Malloc(sizeof(CHAIN *)*2*bscaf,"Allocating alignment space");
    bstack = (int *) Malloc(sizeof(int)*2*bscaf,"Allocating alignment space");
    order  = (ORDER *) Malloc(2*sizeof(ORDER)*amax,"Allocation of order array");
    if (chain == NULL || blist == NULL || bstack == NULL || order == NULL)
      exit (1);

    acnt  = 0;
    ascaf = 0;
    for (b = 0; b < 2*bscaf; b++)
      blist[b] = NULL;
    btop  = 0;

    aread = -1;
    for (j = 0; j <= novl; j++)
      { if (j == novl)
          { ovl->aread = aread+1;
            apulse = 0;
          }
        else
          { Read_Overlap(input,ovl);
            fseek(input,tbytes*ovl->path.tlen,SEEK_CUR);
            apulse = db1->reads[ovl->aread].fpulse;
          }
      
        if (aread != ovl->aread && apulse == 0 && j > 0)
          { CHAIN *u, *v;

            printf("\nScaffold %d\n",ascaf+1);

            collisions(acnt,chain);

            analyze(ascaf,chain,blist,bstack,btop,order);

            for (b = 0; b < btop; b++)
              { int   bs, len, span;
                int   mnum, snum;
                Path *p;
                Path *q;
                int   dela, delb;

                bs = bstack[b];

                for (w = blist[b]; w != NULL; w = w->next)
                  if (w != NULL)
                    w->L = NULL;

                mnum = snum = 0;
                for (w = blist[b]; w != NULL; w = w->next)
                  { if (w->L == NULL)
                      { w->score = (w->ovl.path.aepos - w->ovl.path.abpos);
                        w->clen  = 1;
                      }
                    if (w->link != NULL)
                      { v = w->link;
                        p = &(v->ovl.path);
                        if (v->L != NULL)
                          { if (v->score < w->score + (p->aepos - p->abpos)) 
                              { u = v->L;
                                v->L->link = NULL;
                                v->score = w->score + (p->aepos - p->abpos);
                                v->clen  = w->clen + 1;
                                v->L = w;
                              }
                            else
                              { u = w;
                                w->link = NULL;
                              }
                            len = 0;
                            for (n = u; n != NULL; n = n->L)
                              { n->dead = 1;
                                len    += 1;
                                span    = n->ovl.path.aepos;
                              }
                            span -= u->ovl.path.abpos;
                            if (len > 3 && u->score > .1*span)
                              { // printf("Keeping spur %d with score %6d, %4d links, span = %7d\n",
                                        // snum+1,u->score,len,span);
                                snum += 1;
                                for (n = u; n != NULL; n = n->L)
                                  { n->dead = 2;
                                    n->bepos = snum;
                                  }
                              }
                          }
                        else
                          { v->score = w->score + (p->aepos - p->abpos);
                            v->clen  = w->clen + 1;
                            v->L = w;
                          }
                      }
                    else
                      { len = 0;
                        for (n = w; n != NULL; n = n->L)
                          { len += 1;
                            span = n->ovl.path.aepos;
                          }
                        span -= w->ovl.path.abpos;
                        if (w->score < 10000)
                          { // printf("Removing main with score %6d, %4d links, span = %7d\n",
                                   // w->score,len,span);
                            fflush(stdout);
                            for (n = w; n != NULL; n = n->L)
                              n->dead = 1;
                          }
                        else
                          { mnum += 1;
                            for (n = w; n != NULL; n = n->L)
                              n->bepos = mnum;
                          }
                      }
                  }

/*
                if (bs >= bscaf)
                  printf("\n  Chains with scaffold %d(c)\n",(bs-bscaf)+1);
                else
                  printf("\n  Chains with scaffold %d(n)\n",bs+1);
                for (w = blist[b]; w != NULL; w = w->next)
                  { if (w->dead == 1)
                      continue;
                    p = &(w->ovl.path);
                    printf("    %6ld:",w-chain);
                    printf(" %5d[%9d,%9d] vs %5d[%9d,%9d]",
                           w->ovl.aread+1,w->ovl.path.abpos,w->ovl.path.aepos,
                           w->ovl.bread+1,w->ovl.path.bbpos,w->ovl.path.bepos);
                    if (w->L == NULL)
                      printf("  <%5d += %8d>",w->score,w->score);
                    else
                      printf("  <%5d += %8d>",w->score - w->L->score,w->score);
                    printf("  [%4d]  R=%3d%%",w->clen,w->mark);
                    if (w->dead == 0)
                      printf("  main %3d",w->bepos);
                    else
                      printf("  spur %3d",w->bepos);
                    if (w->link == NULL)
                      printf("  ___\n");
                    else
                      { for (n = w->next; n != NULL; n = n->next)
                          if (n->dead != 1)
                            break;
                        if (n == w->link)
                          printf("  ***      ");
                        else
                          printf("  -> %6ld",w->link-chain);
                        printf("  %6d / %6d\n",w->ovl.path.abpos-w->link->ovl.path.aepos,
                                               w->ovl.path.bbpos-w->link->ovl.path.bepos);
                      }
                  }
*/

                for (w = blist[b]; w != NULL; w = w->next)
                  { if (w->dead == 1)
                      continue;
                    if (w->dead == 0 && w->L == NULL)
                      { int64 slen;
                        int   last;

                        if (bs >= bscaf)
                          last = invscaff[(b-bscaf)+1]-1;
                        else
                          last = invscaff[b+1]-1;
                        slen = db2->reads[last].fpulse + db2->reads[last].rlen;

                        for (u = w; u->link != NULL; u = u->link)
                          ;
                        if (u->score < .001*slen)
                          continue;
                        if (u->score < .01*(w->ovl.path.aepos - u->ovl.path.abpos))
                          continue;
                        if (bs >= bscaf)
                          printf("  vs %4d(c) = %8lld:",(bs-bscaf)+1,slen);
                        else
                          printf("  vs %4d(n) = %8lld:",bs+1,slen);
                        printf(" len = %6d  span = [%9d - %9d] = %9d  covr = %9d = %.2f%%\n",
                               u->clen,u->ovl.path.abpos,w->ovl.path.aepos,
                               w->ovl.path.aepos - u->ovl.path.abpos,u->score,
                               (100.*u->score)/(w->ovl.path.aepos - u->ovl.path.abpos));
                        fflush(stdout);
#define DETAIL
#ifdef DETAIL
                        for (u = w; u != NULL; u = u->link)
                          { p = &(u->ovl.path);
                            printf("    %6ld:",u-chain);
                            printf(" %5d[%9d,%9d] vs %5d[%9d,%9d]",u->ovl.aread,p->abpos,p->aepos,
                                                                   u->ovl.bread,p->bbpos,p->bepos);
                            if (u->L == NULL)
                              printf("  <%5d += %8d>",u->score,u->score);
                            else
                              printf("  <%5d += %8d>",u->score - u->L->score,u->score);
                            printf("  [%4d]  R=%3d%%",u->clen,u->mark);
                            if (u->link == NULL)
                              printf("  ___  END\n");
                            else
                              { q = &(u->link->ovl.path);
                                dela = p->abpos - q->aepos;
                                delb = p->bbpos - q->bepos;
                                if (dela > 50000 || delb > 50000 || abs(dela-delb) > 20000) 
                                  printf("  ___  BREAK  %6d / %6d\n\n",dela,delb);
                                else
                                  { for (n = u->next; n != NULL; n = n->next)
                                      if (n->dead != 1)
                                        break;
                                    if (n == u->link)
                                      printf("  ***      ");
                                    else
                                      printf("  -> %6ld",u->link-chain);
                                    printf("  %6d / %6d\n",dela,delb);
                                  }
                              }
                          }
#endif
                      }
                  }

/*
                for (w = blist[b]; w != NULL; w = w->next)
                  { if (w->dead == 1)
                      continue;
                    if (w->dead == 2 && w->L == NULL)
                      { for (u = w; u->link != NULL; u = u->link)
                          ;
                        printf("  Spur vs %d: len = %6d  span = %9d  covr = %9d\n",
                               bs+1,u->clen,w->ovl.path.aepos - u->ovl.path.abpos,u->score);
                        printf("\n  Spur for scaffold %d\n",bs+1);
                        for (u = w; u != NULL; u = u->link)
                          { p = &(u->ovl.path);
                            printf("    %6ld:",u-chain);
                            printf(" [%9d,%9d] vs [%9d,%9d]",p->abpos,p->aepos,
                                                             p->bbpos,p->bepos);
                            if (u->L == NULL)
                              printf("  <%5d += %8d>",u->score,u->score);
                            else
                              printf("  <%5d += %8d>",u->score - u->L->score,u->score);
                            printf("  [%4d]  R=%3d%%",u->clen,u->mark);
                            if (u->link == NULL)
                              printf("  ___\n");
                            else
                              { q = &(u->link->ovl.path);
                                dela = p->abpos - q->aepos;
                                delb = p->bbpos - q->bepos;
                                if ((dela > 50000 || delb > 50000) && abs(dela-delb) > 10000) 
                                  printf("  ___  BREAK  %6d / %6d\n\n",dela,delb);
                                else
                                  { for (n = u->next; n != NULL; n = n->next)
                                      if (n->dead != 1)
                                        break;
                                    if (n == u->link)
                                      printf("  ***      ");
                                    else
                                      printf("  -> %6ld",u->link-chain);
                                    printf("  %6d / %6d\n",dela,delb);
                                  }
                              }
                          }
                      }
                  }
*/
              }

            if (j >= novl)
              break;

            acnt = 0;
            ascaf += 1;
            for (b = 0; b < btop; b++)
              blist[bstack[b]] = NULL;
            btop = 0;
          }

        aread = ovl->aread;
        bread = ovl->bread;
        ovl->path.abpos += apulse;
        ovl->path.aepos += apulse;
        bpulse = db2->reads[bread].fpulse;
        ovl->path.bbpos += bpulse;
        ovl->path.bepos += bpulse;

        b = scaffold[bread];
        if (COMP(ovl->flags))
          b += bscaf;
        w = blist[b];
        n = chain + acnt++;
        n->ovl  = *ovl;
        n->next = w;
        n->dead = 0;
        blist[b] = n;
        if (w == NULL)
          bstack[btop++] = b;
      }

    free(order);
    free(bstack);
    free(blist);
    free(chain);
  }

  Close_DB(db1);
  if (ISTWO)
    Close_DB(db2);

  exit (0);
}
