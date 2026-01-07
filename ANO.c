/*******************************************************************************************
 *
 *  ONEcode .1ano module.  Auxiliary routines to open and manipulate a ONEcode BED-style file.
 *
 *  Author :  Gene Myers
 *  Date   :  December 2025
 *
 ********************************************************************************************/

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
#include "ANO.h"
#include "GDB.h"

/*******************************************************************************************
 *
 *  GENERAL UTILITIES
 *
 ********************************************************************************************/


static char *anoSchemaText =
  "1 3 def 1 0                 schema for 1ano files\n"
  ".\n"
  "P 3 ano\n"
  ".                           GDB skeleton (may not be presend)\n"
  "O g 0                       groups scaffolds into a GDB skeleton\n"
  "G S                         collection of scaffolds constituting a GDB\n"
  "O S 1 6 STRING              id for a scaffold\n"
  "D G 1 3 INT                 gap of given length\n"
  "D C 1 3 INT                 contig of given length\n"
  ".\n"
  "O M 2 3 INT 8 INT_LIST      scaffold index, beg,end pair list\n"
  "D L 1 6 STRING              optional tab-separated labels for preceeding M\n"
  "D X 1 8 INT_LIST            optional list of scores for the preceeding M\n"
;

OneSchema *make_ANO_Schema()
{ return (oneSchemaCreateFromText(anoSchemaText)); }

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
        { EPRINTF("Out of memory (Cating 4 strings)");
          return (NULL);
        }
    }
  sprintf(cat,"%s%s%s%s",path,sep,root,suffix);
  return (cat);
}

void Show_ANO(ANO *ano)
{ int i;

  for (i = 0; i < ano->nints; i++)
    { printf("Mask: %lld %lld",ano->masks[i].beg,ano->masks[i].end);
      if (ano->masks[i].label != NULL)
        printf(" %s\n",ano->masks[i].label);
      printf("\n");
    }
}


/*******************************************************************************************
 *
 *  ANO ROUTINES
 *
 ********************************************************************************************/

static int PSORT(const void *l, const void *r)
{ ANO_PAIR *x = (ANO_PAIR *) l;
  ANO_PAIR *y = (ANO_PAIR *) r;

  if (x->beg < y->beg)
    return (-1);
  if (x->beg > y->beg)
    return (1);
  return (0);
}

int Read_ANO(ANO *ano, char *path, GDB *gdb)
{ OneSchema *schema;
  OneFile   *of1, *of2;
  GDB       _skel, *skel = &_skel;
  GDB       *mydb;
  int        shared;

  int64          nints, nprov, nscaff;
  int64          labtot, psize, nlines;
  int            maxlab;
  int64         *moff;
  ANO_PAIR      *mask;
  char          *labs;
  char          *source;
  OneProvenance *prov;

  //  Open ANO OneFile and create schema

  { char *e; 
    char *root, *pwd;
    char *fname;               
    FILE *f;
  
    pwd = PathTo(path);        
    e = path + strlen(path);   
    if (e > path+5 && strcmp(e-5,".1ano") == 0)
      root = Root(path,".1ano");
    else 
      root = Root(path,".ano");
    fname = MyCatenate(pwd,"/",root,".1ano");
    if (pwd == NULL || root == NULL || fname == NULL)
      { free(root);
        free(pwd);
        EXIT(1);
      }
    f = fopen(fname,"r");
    if (f == NULL)
      { f = fopen(MyCatenate(pwd,"/",root,".ano"),"r");
        if (f == NULL)
          { EPRINTF("Cannot find/open ANO file %s",path);
            free(root);
            free(pwd);
            EXIT(1);
          } 
        fclose(f);
      }
    else
      fclose(f);

    schema = make_ANO_Schema();
    if (schema == NULL)
      { EPRINTF("Failed to create ANO schema");
        free(root);
        free(pwd);
        EXIT(1);
      }
 
    of1 = oneFileOpenRead(fname,schema,"ano",2);

    if (of1 == NULL)
      { EPRINTF("Failed to open .1ano file %s",path);
        oneSchemaDestroy(schema);
        free(root);
        free(pwd);
        EXIT(1);
      }
    free(root);
    free(pwd);

    of2 = of1+1;
  }

  //  Determine if anchored and whether it has a skeleton

  if (of1->info['<'] == NULL || of1->reference[0].count != 1)
    { EPRINTF(".1ano does not contain a gdb reference");
      oneFileClose(of1);
      oneSchemaDestroy(schema);
      EXIT(1);
    }
  source = of1->reference[0].filename;

  oneReadLine(of1);
  if (of1->lineType != 'g')
    { EPRINTF(".1ano does not contain prefacing GDB skeleton");
      oneFileClose(of1);
      oneSchemaDestroy(schema);
      EXIT(1);
    }
  if (gdb == NULL)
    { mydb   = skel;    //  temporary, will ultimately point at allocated GDB record.
      shared = 0;
      Read_Skeleton(of1,source,skel);
    }
  else
    { mydb   = gdb;
      shared = 1;
    }

  //  Get sizes of things and allocate memory for ANO object

  nprov = of1->info['!']->accum.count;  //  Sigh, "given" counts don't work for provenance
  oneStats(of1,'M',NULL,NULL,&nints);
  oneStats(of1,'L',&nlines,NULL,&labtot);
  nints  /= 2;
  labtot += nlines;
  nscaff  = mydb->nscaff;

  { int i;

    psize = nprov * sizeof(OneProvenance);
    for (i = 0; i < nprov; i++)
      psize += strlen(of1->provenance[i].program)
             + strlen(of1->provenance[i].version)
             + strlen(of1->provenance[i].command)
             + strlen(of1->provenance[i].date) + 4;
  }

  if (!shared)
    { mydb  = malloc(sizeof(GDB));
      *mydb = *skel;
    }
  moff = malloc(sizeof(int64)*(mydb->ncontig+1));
  mask = malloc(sizeof(ANO_PAIR)*(nints+1));
  if (labtot > 0)
    labs = malloc(labtot);
  else
    labs = NULL;
  if (psize > 0)
    prov = malloc(psize);
  else
    prov = NULL;
  if (mydb == NULL || moff == NULL || mask == NULL || (labtot > 0 && labs == NULL)
                                                   || (psize  > 0 && prov == NULL)  )
    { EPRINTF("Could not allocate memory for ANO (Read_ANO)");
      goto error;
    }

  //  Capture provenance

  { int   i;
    char *pstr;

    pstr = (char *) (prov + nprov);
    for (i = 0; i < nprov; i++)
      { prov[i].program = pstr;
        pstr = stpcpy(pstr,of1->provenance[i].program) + 1;
        prov[i].version = pstr;
        pstr = stpcpy(pstr,of1->provenance[i].version) + 1;
        prov[i].command = pstr;
        pstr = stpcpy(pstr,of1->provenance[i].command) + 1;
        prov[i].date = pstr;
        pstr = stpcpy(pstr,of1->provenance[i].date) + 1;
      }
  }

  //  1st pass to count and check things

  { int   *map;
    int64 *count;
    int    s;

    map   = malloc(sizeof(int64)*2*nscaff);
    count = ((int64 *) map) + nscaff;
    if (map == NULL)
      { EPRINTF("Could not allocate memory for ANO (Read_ANO)");
        free(map);
        goto error;
      }

    if (shared)
      { Read_Skeleton(of1,source,skel);
        if (Are_Skeletons_Equal(gdb,skel,map) == 0)
          { EPRINTF("GDB structures not equivalent");
            Close_GDB(skel);
            goto error;
          }
        Close_GDB(skel);
      }
    else
      { for (s = 0; s < nscaff; s++)
          map[s] = s;
      }

    { int   nextL, nextX, len;
      int   scf;
      char *lab;

      for (s = 0; s < nscaff; s++)
        count[s] = 0;

      nextL = 0;
      nextX = 0;
      do
        switch (of1->lineType)
        { case 'M':
            scf = oneInt(of1,0);
            if (scf >= nscaff)
              { EPRINTF("%d'th scaffold not declared in prolog",scf+1);
                free(map);
                goto error;
              }
            len = oneLen(of1) >> 1;
            count[map[scf]] += len;
            nextL = len;
            nextX = len;
            break;
          case 'L':
            len = 1;
            lab = oneString(of1);
            for (lab = index(lab,'\t'); lab != NULL; lab = index(lab,'\t'))
              { lab += 1;
                len += 1;
              }
            if (nextL != len)
              { EPRINTF("L-line must immediately follow an M-line and have same list length");
                free(map);
                goto error;
              }
	    nextL = 0;
            break;
          case 'X':
            len = oneLen(of1);
            if (nextX != len)
              { EPRINTF("X-line must immediately follow an M-line and have same list length");
                free(map);
                goto error;
              }
	    nextX = 0;
            break;
          default:
            EPRINTF("Do not recognize line of type %c(%d)",of1->lineType,of1->lineType);
            free(map);
            goto error;
        }
      while (oneReadLine(of1));
    }

    //  Set up load indices so all intervals for each scaffold are contiguous in mask

    { int64 sum, x;
      int   i;

      sum = 0;
      for (i = 0; i < nscaff; i++)
        { x = count[i];
          count[i] = sum;
          sum += x;
        }
    }

    //  2nd pass this time loading data

    { int64  ltop, *list;
      char  *str;
      int    len, ldx, idx, scf;
      int    n;

      oneReadLine(of2);
      Skip_Skeleton(of2);

      maxlab = 0;
      ltop   = 0;
      do
        switch (of2->lineType)
        { case 'M':
            scf  = map[oneInt(of2,0)];
            idx  = ldx = count[scf];
            len  = oneLen(of2);
            list = oneIntList(of2);
            for (n = 0; n < len; n += 2)
              { if (list[n] < list[n+1])
                  { mask[idx].beg    = list[n];
                    mask[idx].end    = list[n+1];
                    mask[idx].orient = 0;
                  }
                else
                  { mask[idx].beg    = list[n+1];
                    mask[idx].end    = list[n];
                    mask[idx].orient = 1;
                  }
                mask[idx].label = NULL;
                mask[idx].score = 0;
                idx += 1;
              }
            count[scf] = idx;
            break;
          case 'L':
            idx = ldx;
            str = labs+ltop;
            strcpy(str,oneString(of2));
            ltop += oneLen(of2);
            labs[ltop++] = '\0';
            mask[idx++].label = str; 
            for (str = index(str,'\t'); str != NULL; str = index(str,'\t'))
              { *str++ = '\0';
                mask[idx++].label = str;
              }
            break;
          case 'X':
            idx  = ldx;
            len  = oneLen(of2);
            list = oneIntList(of2);
            for (n = 0; n < len; n++)
              mask[idx++].score = list[n];
            break;
          default:
            printf("Should not happen %c(%d)\n",of2->lineType,of2->lineType);
            break;
        }
      while (oneReadLine(of2));
    }

    //  sort ANO intervals for each scaffold if not already sorted

    { int64 beg, j;
      int   i;

      beg = 0;
      for (i = 0; i < nscaff; i++)
        { for (j = beg+1; j < count[i]; j++)
            if (mask[j-1].beg > mask[j].beg)
              { qsort(mask+beg,count[i]-beg,sizeof(ANO_PAIR),PSORT);
                break;
              }
          moff[i] = beg;
          beg = count[i];
        }
      moff[nscaff] = beg;
    }

    //  Map from scaffold coords to contig coords and set contig partition in moff

    { int64 bot, top;
      int   scf, ctg;
      int   i;

      GDB_SCAFFOLD *scaffs = mydb->scaffolds;
      GDB_CONTIG   *contig = mydb->contigs;

      scf = 0;
      ctg = 0;
      bot = 0;
      if (scaffs[0].ectg <= 1)
        top = 0x7fffffffffffffffll;
      else
        top = contig[1].sbeg;
      for (i = 0; i < nints; i++)
        { while (i == count[scf])
            { scf += 1;
              while (ctg < scaffs[scf].fctg)
                { ctg += 1;
                  moff[ctg] = i;
                }
              bot = 0;
              if (scaffs[scf].ectg <= ctg+1)
                top = 0x7fffffffffffffffll;
              else
                top = contig[ctg+1].sbeg;
            }
          while (mask[i].beg >= top)
            { ctg += 1;
              moff[ctg] = i;
              bot = top;
              if (scaffs[scf].ectg <= ctg+1)
                top = 0x7fffffffffffffffll;
              else
                top = contig[ctg+1].sbeg;
            }
          mask[i].beg -= bot;
          mask[i].end -= bot;
        }
      while (ctg < mydb->ncontig)
        { ctg += 1;
          moff[ctg] = nints;
        }
    }

    free(map);
  }

  ano->nprov  = nprov;
  ano->prov   = prov;

  ano->gdb     = mydb;
  ano->shared  = shared;
  ano->nints   = nints;
  ano->moff    = moff;
  ano->masks   = mask;
  ano->labels  = labs;
  ano->maxlab  = maxlab;

  (void) Show_ANO;

  oneFileClose(of1);
  oneSchemaDestroy(schema);
  return (0);

error:
  free(prov);
  free(labs);
  free(mask);
  free(moff);
  if (!shared)
    { Close_GDB(mydb);
      free(mydb);
    }
  oneFileClose(of1);
  oneSchemaDestroy(schema);
  EXIT(1);
}

  // Write the given ano to the file 'tpath'.
  
extern bool addProvenance(OneFile *of, OneProvenance *from, int n) ; // backdoor - clean up some day
  
int Write_ANO(ANO *ano, char *tpath, int width)   
{ OneSchema *schema;                   
  OneFile   *of;                       
  bool       binary;                   
  GDB       *gdb = ano->gdb;

  { char *e;                   
    char *root, *pwd;          
  
    pwd = PathTo(tpath);       
    e = tpath + strlen(tpath); 
    if (e > tpath+4 && strcmp(e-4,".ano") == 0)
      { root = Root(tpath,".ano");
        binary = false;
      }
    else
      { root = Root(tpath,".1ano");
        binary = true;
      }
    free(pwd);
    free(root);
  } 
      
  schema = make_ANO_Schema();
  if (schema == NULL)
    { EPRINTF("Failed to create ANO schema (Write_ANO)");
      EXIT(1);
    }

  of = oneFileOpenWriteNew(tpath,schema,"ano",binary,1);
  if (of == NULL)
    { EPRINTF("Failed to open ANO file %s (Write_ANO)",tpath);
      oneSchemaDestroy(schema);
      EXIT(1);
    }

  addProvenance(of,ano->prov,ano->nprov);
  oneAddProvenance(of,Prog_Name,"0.1",Command_Line);

  oneAddReference(of,gdb->srcpath,1);

  Write_Skeleton(of,gdb);

  { int        idx, scf, nobj, scored;
    ANO_PAIR  *mask;
    int64     *moff;
    char      *s;
    int64      base;
    int        j, p, k, end, len;
    int64      mbuffer[2*width];
    char       mstring[width*(ano->maxlab+1)];

    nobj = gdb->ncontig;
    mask = ano->masks;
    moff = ano->moff;
    for (idx = 0; idx < nobj; idx++)
      { base = gdb->contigs[idx].sbeg;
        scf  = gdb->contigs[idx].scaf;
        end  = moff[idx+1];
        for (j = moff[idx]; j < end; j += width)
          { oneInt(of,0) = scf;
            if (j + width > end)
              len = 2*(end-j);
            else
              len = 2*width;

            p = j;
            for (k = 0; k < len; k += 2)
              if (mask[p].orient)
                { mbuffer[k]   = mask[p].end + base;
                  mbuffer[k+1] = mask[p++].beg + base;
                }
              else
                { mbuffer[k]   = mask[p].beg + base;
                  mbuffer[k+1] = mask[p++].end + base;
                }
            oneWriteLine(of,'M',len,mbuffer);

            len /= 2;
            s = mstring;
            p = j;
            for (k = 0; k < len; k++)
              { if (mask[p].label == NULL)
                  p += 1;
                else
                  s = stpcpy(s,mask[p++].label);
                *s++ = '\t';
              }
            s[-1] = '\0';
            if (s - mstring > len)
              oneWriteLine(of,'L',s-mstring,mstring);

            p = j;
            scored = 0;
            for (k = 0; k < len; k++)
              { mbuffer[k] = mask[p].score;
                if (mask[p].score > 0)
                  scored = 1;
                p += 1;
              }
            if (scored)
              oneWriteLine(of,'X',len,mbuffer);
          }
      }
  }

  oneFileClose(of);
  return (0);
}

void Free_ANO(ANO *ano)
{ free(ano->labels);
  free(ano->masks);
  free(ano->moff);
  free(ano->prov);
  if (!ano->shared)
    { Close_GDB(ano->gdb);
      free(ano->gdb);
    }
}

static void heapify(int s, ANO_PAIR **heap, int hsize)
{ int       c, l, r;
  ANO_PAIR *hs, *hr, *hl;

  c  = s;
  hs = heap[s];
  while ((l = (c<<1)) <= hsize)
    { r  = l+1;
      hl = heap[l];
      if (r > hsize)
        { if (hl->beg < hs->beg)   
            { heap[c] = hl;
              c = l; 
            }            
          break;
        }
      hr = heap[r];  
      if (hl->beg > hr->beg)
        { if (hl->beg < hs->beg)   
            { heap[c] = hl;
              c = l; 
            }            
          else
            break;
        }
      else    
        { if (hr->beg < hs->beg)
            { heap[c] = hr;
              c = r;
            }
          else
            break;
        }
    }
  heap[c] = hs;
}

int ANO_Union(ANO *tano,int nway, ANO **sanos)
{ int            nints, nscaff, nprov, mlen;
  GDB           *gdb;
  int64         *moff;
  ANO_PAIR      *masks;
  OneProvenance *prov;

  gdb    = sanos[0]->gdb;
  nscaff = gdb->nscaff;
  mlen   = gdb->ncontig;

  { int  i, j;
    int *map;

    for (i = 1; i < nway; i++)
      if (sanos[i]->gdb->nscaff != nscaff)
        { EPRINTF("ANOs do not have the same scaffold structure");
          EXIT (1);
        }

    map = malloc(sizeof(int)*nscaff);

    for (i = 1; i < nway; i++)
      { if (!Are_Skeletons_Equal(gdb,sanos[i]->gdb,map))
          { EPRINTF("ANOs do not have the same scaffold structure");
            free(map);
            EXIT(1);
          }
        for (j = 0; j < nscaff; j++)
          if (map[j] != j)
            { EPRINTF("ANOs do not have the same scaffold structure");
              free(map);
              EXIT(1);
            }
      }

    free(map);
  }

  //  1st pass, determine how big union mask will be

  { int64     beg[nway];
    ANO_PAIR *heap[nway+1], *pair;
    int       hsize;
    int64     end;
    int       i, j;

    nints = 0;
    for (i = 0; i < nway; i++)
      beg[i] = sanos[i]->masks[0].beg;
    for (j = 0; j < mlen; j++)
      { hsize = 0;
        for (i = 0; i < nway; i++)
          { ANO_PAIR *mk;
            int64    *mo;

            mk = sanos[i]->masks;
            mo = sanos[i]->moff;
            mk[mo[j]].beg = beg[i];
            if (mo[j+1] - mo[j] > 0)
              heap[++hsize] = mk + mo[j];
            beg[i] = mk[mo[j+1]].beg;
            mk[mo[j+1]].beg = -1;
          }

        if (hsize > 3)
          for (i = hsize/2; i > 1; i--)
            heapify(i,heap,hsize);

        end = -1;
        while (hsize > 0)
          { heapify(1,heap,hsize);
            pair = heap[1];
  
            if (pair->beg > end)
              { nints += 1;
                end = pair->end;
              }
            else if (pair->end > end)
              end = pair->end;
  
            pair += 1;
            if (pair->beg < 0)
              { heap[1] = heap[hsize];
                hsize -= 1;
              }
            else
              heap[1] = pair;
          }
      }
    for (i = 0; i < nway; i++)
      sanos[i]->masks[sanos[i]->moff[mlen]].beg = beg[i];
  }

  //  Allocate union ANO components

  { char *last;
    int   i;
    int   psize;

    nprov = 0;
    psize = 0;
    for (i = 0; i < nway; i++)
      if (sanos[i]->nprov > 0)
        { nprov += sanos[i]->nprov;
          last = sanos[i]->prov[sanos[i]->nprov-1].date;
          psize += (last + strlen(last) + 1) - sanos[i]->prov[0].program;
        }
    psize += nprov * sizeof(OneProvenance);

    moff    = malloc(sizeof(int64)*(mlen+1));
    masks   = malloc(sizeof(ANO_PAIR)*(nints+1));
    if (psize > 0)
      prov  = malloc(psize);
    else
      prov  = NULL;
    if (moff == NULL || masks == NULL || (psize > 0 && prov == NULL) )
      { EPRINTF("Could not allocate memory for ANO (ANO_Union)");
        goto error;
      }
  }

  //  Provenance is the concatenation of the provenance of each source ANO

  { int   i, j, k;
    char *pstr;

    pstr = (char *) (prov + nprov);
    k = 0;
    for (j = 0; j < nway; j++)
      for (i = 0; i < sanos[j]->nprov; i++)
        { prov[k].program = pstr;
          pstr = stpcpy(pstr,sanos[j]->prov[i].program) + 1;
          prov[k].version = pstr;
          pstr = stpcpy(pstr,sanos[j]->prov[i].version) + 1;
          prov[k].command = pstr;
          pstr = stpcpy(pstr,sanos[j]->prov[i].command) + 1;
          prov[k].date = pstr;
          pstr = stpcpy(pstr,sanos[j]->prov[i].date) + 1;
          k += 1;
        }
  }

  //  2nd pass, form the union in the now allocated moff and masks arrays

  { int64     beg[nway];
    ANO_PAIR *heap[nway+1], *pair;
    int       hsize, mcur;
    int64     end;
    int       i, j;

    mcur = 0;
    for (i = 0; i < nway; i++)
      beg[i] = sanos[i]->masks[0].beg;
    for (j = 0; j < mlen; j++)
      { moff[j] = mcur;
        hsize = 0;
        for (i = 0; i < nway; i++)
          { ANO_PAIR *mk;
            int64    *mo;

            mk = sanos[i]->masks;
            mo = sanos[i]->moff;
            mk[mo[j]].beg = beg[i];
            if (mo[j+1] - mo[j] > 0)
              heap[++hsize] = mk + mo[j];
            beg[i] = mk[mo[j+1]].beg;
            mk[mo[j+1]].beg = -1;
          }

        if (hsize == 0)
          continue;

        for (i = hsize/2; i >= 1; i--)
          heapify(i,heap,hsize);

        masks[mcur].beg   = heap[1]->beg;
        masks[mcur].label = NULL;
        end = heap[1]->end;

        while (hsize > 0)
          { heapify(1,heap,hsize);
            pair = heap[1];
  
            if (pair->beg > end)
              { masks[mcur].end = end;
                mcur += 1;
                masks[mcur].beg    = pair->beg;
                masks[mcur].label  = NULL;
                masks[mcur].orient = 0;
                masks[mcur].score  = 0;
                end = pair->end;
              }
            else if (pair->end > end)
              end = pair->end;
  
            pair += 1;
            if (pair->beg < 0)
              { heap[1] = heap[hsize];
                hsize -= 1;
              }
            else
              heap[1] = pair;
          }
        masks[mcur].end = end;
        mcur += 1;
      }
    moff[mlen] = mcur;
    for (i = 0; i < nway; i++)
      sanos[i]->masks[sanos[i]->moff[mlen]].beg = beg[i];
  }

  tano->nprov   = nprov;
  tano->prov    = prov;

  tano->gdb     = gdb;

  tano->nints   = nints;
  tano->moff    = moff;
  tano->masks   = masks;
  tano->labels  = NULL;

  return (0);

error:
  free(prov);
  free(masks);
  free(moff);
  EXIT(1);
}
