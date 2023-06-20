/*******************************************************************************************
 *
 *  Compressor/decompressor for .quiv files: customized Huffman codes for each stream based on
 *    the histogram of values occuring in a given file.  The two low complexity streams
 *    (deletionQV and substitutionQV) use a Huffman coding of the run length of the prevelant
 *    character.
 *
 *  Author:   Gene Myers
 *  Date:     Jan 18, 2014
 *  Modified: July 25, 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "DB.h"

#undef DEBUG

#define MIN_BUFFER 1000

#define HUFF_CUTOFF  16   //  This cannot be larger than 16 !


/*******************************************************************************************
 *
 *  Endian flipping routines
 *
 ********************************************************************************************/

static int LittleEndian;  //  Little-endian machine ?
                          //     Referred by: Decode & Decode_Run
static int Flip;          //  Flip endian of all coded shorts and ints
                          //     Referred by: Decode & Decode_Run & Read_Scheme

static void Set_Endian(int flip)
{ uint32 x = 3;
  uint8 *b = (uint8 *) (&x);

  Flip         = flip;
  LittleEndian = (b[0] == 3);
}

static void Flip_Long(void *w)
{ uint8 *v = (uint8 *) w;
  uint8  x;

  x    = v[0];
  v[0] = v[3];
  v[3] = x;
  x    = v[1];
  v[1] = v[2];
  v[2] = x;
}

static void Flip_Short(void *w)
{ uint8 *v = (uint8 *) w;
  uint8  x;

  x    = v[0];
  v[0] = v[1];
  v[1] = x;
}


/*******************************************************************************************
 *
 *  Routines for computing a Huffman Encoding Scheme
 *
 ********************************************************************************************/

typedef struct
  { int    type;             //  0 => normal, 1 => normal but has long codes, 2 => truncated
    uint32 codebits[256];    //  If type = 2, then code 255 is the special code for
    int    codelens[256];    //    non-Huffman exceptions
    int    lookup[0x10000];  //  Lookup table (just for decoding)
  } HScheme;

typedef struct _HTree
  { struct _HTree *lft, *rgt; 
    uint64         count;
  } HTree;

  //  Establish heap property from node s down (1 is root, siblings of n are 2n and 2n+1)
  //    assuming s is the only perturbation in the tree.

static void Reheap(int s, HTree **heap, int hsize)
{ int      c, l, r;
  HTree   *hs, *hr, *hl;

  c  = s;
  hs = heap[s];
  while ((l = 2*c) <= hsize)
    { r  = l+1;
      hl = heap[l];
      hr = heap[r];
      if (r > hsize || hr->count > hl->count)
        { if (hs->count > hl->count)
            { heap[c] = hl;
              c = l;
            }
          else
            break;
        }
      else
        { if (hs->count > hr->count)
            { heap[c] = hr;
              c = r;
            }
          else
            break;
        }
    }
  if (c != s)
    heap[c] = hs;
}

  //  Given Huffman tree build a table of codes from it, the low-order codelens[s] bits
  //    of codebits[s] contain the code for symbol s.

static void Build_Table(HTree *node, int code, int len, uint32 *codebits, int *codelens)
{ if (node->rgt == NULL)
    { uint64 symbol = (uint64) (node->lft);
      codebits[symbol] = code;
      codelens[symbol] = len;
    }
  else
    { code <<= 1;
      len   += 1;
      Build_Table(node->lft,code,len,codebits,codelens);
      Build_Table(node->rgt,code+1,len,codebits,codelens);
    }
}

  // For the non-zero symbols in hist, compute a huffman tree over them, and then
  //   build a table of the codes.  If inscheme is not NULL, then place all symbols
  //   with code 255 or with more than HUFF_CUTOFF bits in the encoding by inscheme
  //   as a single united entity, whose code signals that the value of these symbols
  //   occur explicitly in 8 (values) or 16 (run lengths) bits following the code.
  //   All the symbols in this class will have the same entry in the code table and
  //   255 is always in this class.

static HScheme *Huffman(uint64 *hist, HScheme *inscheme)
{ HScheme *scheme;
  HTree   *heap[259];
  HTree    node[512];
  int      hsize;
  HTree   *lft, *rgt;
  int      value, range;
  int     i;

  scheme = (HScheme *) Malloc(sizeof(HScheme),"Allocating Huffman scheme record");
  if (scheme == NULL)
    return (NULL);

  hsize = 0;                        //  Load heap
  value = 0;
  if (inscheme != NULL)
    { node[0].count = 0;
      node[0].lft   = (HTree *) (uint64) 255;
      node[0].rgt   = NULL;
      heap[++hsize] = node+(value++);
    }
  for (i = 0; i < 256; i++)
    if (hist[i] > 0)
      { if (inscheme != NULL && (inscheme->codelens[i] > HUFF_CUTOFF || i == 255))
          node[0].count += hist[i];
        else
          { node[value].count = hist[i];
            node[value].lft   = (HTree *) (uint64) i;
            node[value].rgt   = NULL;
            heap[++hsize] = node+(value++);
          }
      }

  for (i = hsize/2; i >= 1; i--)    //  Establish heap property
    Reheap(i,heap,hsize);

  range = value;                    //   Merge pairs with smallest count until have a tree
  for (i = 1; i < value; i++)
    { lft = heap[1];
      heap[1] = heap[hsize--];
      Reheap(1,heap,hsize);
      rgt = heap[1];
      node[range].lft = lft;
      node[range].rgt = rgt;
      node[range].count = lft->count + rgt->count;
      heap[1] = node+(range++);
      Reheap(1,heap,hsize);
    }

  for (i = 0; i < 256; i++)        //  Build the code table
    { scheme->codebits[i] = 0;
      scheme->codelens[i] = 0;
    }

  Build_Table(node+(range-1),0,0,scheme->codebits,scheme->codelens);

  if (inscheme != NULL)            //  Set scheme type and if truncated (2), map truncated codes
    { scheme->type = 2;            //    to code and length for 255
      for (i = 0; i < 255; i++)
        if (inscheme->codelens[i] > HUFF_CUTOFF || scheme->codelens[i] > HUFF_CUTOFF)
          { scheme->codelens[i] = scheme->codelens[255];
            scheme->codebits[i] = scheme->codebits[255];
          }
    }
  else
    { scheme->type = 0;
      for (i = 0; i < 256; i++)
        { if (scheme->codelens[i] > HUFF_CUTOFF)
            scheme->type = 1;
        }
    }

  return (scheme);
}

#ifdef DEBUG

  //  For debug, show the coding table

static void Print_Table(HScheme *scheme, uint64 *hist, int infosize)
{ uint64 total_bits;
  uint32 specval, mask, code, *bits;
  int    speclen, clen, *lens;
  int    i, k;

  total_bits = 0;
  bits = scheme->codebits;
  lens = scheme->codelens;
  if (scheme->type == 2)
    { specval = bits[255];
      speclen = lens[255];
    }
  else
    specval = speclen = 0x7fffffff;

  printf("\nCode Table:\n");
  for (i = 0; i < 256; i++)
    if (lens[i] > 0)
      { clen = lens[i];
        mask = (1 << clen);
        code = bits[i];
        printf(" %3d: %2d ",i,clen);
        for (k = 0; k < clen; k++)
          { mask >>= 1;
            if (code & mask)
              printf("1");
            else
              printf("0");
          }
        if (code == specval && clen == speclen)
          { printf(" ***");
            if (hist != NULL)
              total_bits += (clen+infosize)*hist[i];
          }
        else if (hist != NULL)
          total_bits += clen*hist[i];
        printf("\n");
      }
  if (hist != NULL)
    printf("\nTotal Bytes = %lld\n",(total_bits-1)/8+1);
}

  //  For debug, show the histogram

static void Print_Histogram(uint64 *hist)
{ int    i, low, hgh;
  uint64 count;

  for (hgh = 255; hgh >= 0; hgh--)
    if (hist[hgh] != 0)
      break;
  for (low = 0; low < 256; low++)
    if (hist[low] != 0)
      break;
  count = 0;
  for (i = low; i <= hgh; i++)
    count += hist[i];

  for (i = hgh; i >= low; i--)
    printf("    %3d: %8llu %5.1f%%\n",i,hist[i],(hist[i]*100.)/count);
}

#endif


/*******************************************************************************************
 *
 *  Read and Write Huffman Schemes
 *
 ********************************************************************************************/

  //  Write the code table to out.

static void Write_Scheme(HScheme *scheme, FILE *out)
{ int     i;
  uint8   x;
  uint32 *bits;
  int    *lens;

  lens = scheme->codelens;
  bits = scheme->codebits;

  x = (uint8) (scheme->type);
  fwrite(&x,1,1,out);

  for (i = 0; i < 256; i++)
    { x = (uint8) (lens[i]);
      fwrite(&x,1,1,out);
      if (x > 0)
        fwrite(bits+i,sizeof(uint32),1,out);
    }
}

  //  Allocate and read a code table from in, and return a pointer to it.

static HScheme *Read_Scheme(FILE *in)
{ HScheme *scheme;
  int     *look, *lens;
  uint32  *bits, base;
  int      i, j, powr;
  uint8    x;

  scheme = (HScheme *) Malloc(sizeof(HScheme),"Allocating Huffman scheme record");
  if (scheme == NULL)
    return (NULL);

  lens = scheme->codelens;
  bits = scheme->codebits;
  look = scheme->lookup;

  if (fread(&x,1,1,in) != 1)
    { EPRINTF(EPLACE,"Could not read scheme type byte (Read_Scheme)\n");
      free(scheme);
      return (NULL);
    }
  scheme->type = x;
  for (i = 0; i < 256; i++)
    { if (fread(&x,1,1,in) != 1)
        { EPRINTF(EPLACE,"Could not read length of %d'th code (Read_Scheme)\n",i);
          return (NULL);
        }
      lens[i] = x;
      if (x > 0)
        { if (fread(bits+i,sizeof(uint32),1,in) != 1)
            { EPRINTF(EPLACE,"Could not read bit encoding of %d'th code (Read_Scheme)\n",i);
              free(scheme);
              return (NULL);
            }
        }
      else
        bits[i] = 0;
    }

  if (Flip)
    { for (i = 0; i < 256; i++)
        Flip_Long(bits+i);
    }

  for (i = 0; i < 256; i++)
    { if (lens[i] > 0)
        { base = (bits[i] << (16-lens[i]));
          powr = (1 << (16-lens[i]));
          for (j = 0; j < powr; j++)
            look[base+j] = i;
        }
    }

  return (scheme);
}


/*******************************************************************************************
 *
 *  Encoders and Decoders
 *
 ********************************************************************************************/

  //  Encode read[0..rlen-1] according to scheme and write to out

static void Encode(HScheme *scheme, FILE *out, uint8 *read, int rlen)
{ uint32  x, c, ocode;
  int     n, k, olen, llen;
  int    *nlens;
  uint32 *nbits;
  uint32  nspec;
  int     nslen;

  nlens = scheme->codelens;
  nbits = scheme->codebits;

  if (scheme->type == 2)
    { nspec = nbits[255];
      nslen = nlens[255];
    }
  else
    nspec = nslen = 0x7fffffff;

#define OCODE(L,C)				\
{ int    len  = olen + (L);			\
  uint32 code = (C);				\
						\
  llen = olen;					\
  if (len >= 32)				\
    { olen   = len-32;				\
      ocode |= (code >> olen);			\
      fwrite(&ocode,sizeof(uint32),1,out);	\
      if (olen > 0)				\
        ocode = (code << (32-olen));		\
      else					\
        ocode = 0;				\
    } 						\
  else						\
    { olen   = len;				\
      ocode |= (code << (32-olen));;		\
    }						\
}

  llen  = 0;
  olen  = 0;
  ocode = 0;
  for (k = 0; k < rlen; k++)
    { x = read[k];
      n = nlens[x];
      c = nbits[x];
      OCODE(n,c);
      if (c == nspec && n == nslen)
        OCODE(8,x);
    }

  if (olen > 0)                              //  Tricky: must pad so decoder does not read past
    { fwrite(&ocode,sizeof(uint32),1,out);   //    last integer int the coded output.
      if (llen > 16 && olen > llen)
        fwrite(&ocode,sizeof(uint32),1,out);
    }
  else if (llen > 16)
    fwrite(&ocode,sizeof(uint32),1,out);
}

  //  Encode read[0..rlen-1] according to non-rchar table neme, and run-length table reme for
  //    runs of rchar characters.  Write to out.

static void Encode_Run(HScheme *neme, HScheme *reme, FILE *out, uint8 *read, int rlen, int rchar)
{ uint32  x, c, ocode;
  int     n, h, k, olen, llen;
  int    *nlens, *rlens;
  uint32 *nbits, *rbits;
  uint32  nspec, rspec;
  int     nslen, rslen;

  nlens = neme->codelens;
  nbits = neme->codebits;
  rlens = reme->codelens;
  rbits = reme->codebits;

  if (neme->type == 2)
    { nspec = nbits[255];
      nslen = nlens[255];
    }
  else
    nspec = nslen = 0x7fffffff;

  rspec = rbits[255];
  rslen = rlens[255];

  llen  = 0;
  olen  = 0;
  ocode = 0;
  k     = 0;
  while (k < rlen)
    { h = k;
      while (k < rlen && read[k] == rchar)
        k += 1;
      if (k-h >= 255)
        x = 255;
      else
        x = k-h;
      n = rlens[x];
      c = rbits[x];
      OCODE(n,c);
      if (c == rspec && n == rslen)
        OCODE(16,k-h);
      if (k < rlen)
        { x = read[k];
          n = nlens[x];
          c = nbits[x];
          OCODE(n,c);
          if (c == nspec && n == nslen)
            OCODE(8,x);
          k += 1;
        }
    }

  if (olen > 0)
    { fwrite(&ocode,sizeof(uint32),1,out);
      if (llen > 16 && olen > llen)
        fwrite(&ocode,sizeof(uint32),1,out);
    }
  else if (llen > 16)
    fwrite(&ocode,sizeof(uint32),1,out);
}

  //  Read and decode from in, the next rlen symbols into read according to scheme

static int Decode(HScheme *scheme, FILE *in, char *read, int rlen)
{ int    *look, *lens;
  int     signal, ilen;
  uint64  icode;
  uint32 *ipart;
  uint16 *xpart;
  uint8  *cpart;
  int     j, n, c;

  if (LittleEndian)
    { ipart = ((uint32 *) (&icode));
      xpart = ((uint16 *) (&icode)) + 2;
      cpart = ((uint8  *) (&icode)) + 5;
    }
  else
    { ipart = ((uint32 *) (&icode)) + 1;
      xpart = ((uint16 *) (&icode)) + 1;
      cpart = ((uint8  *) (&icode)) + 2;
    }

  if (scheme->type == 2)
    signal  = 255;
  else
    signal  = 256;
  lens = scheme->codelens;
  look = scheme->lookup;

#define GET								\
  if (n > ilen)								\
    { icode <<= ilen;							\
      if (fread(ipart,sizeof(uint32),1,in) != 1)			\
        { EPRINTF(EPLACE,"Could not read more bits (Decode)\n");	\
          return (1);							\
        }								\
      ilen    = n-ilen;							\
      icode <<= ilen;							\
      ilen    = 32-ilen;						\
    }									\
  else									\
    { icode <<= n;							\
      ilen   -= n;							\
    }

#define GETFLIP								\
  if (n > ilen)								\
    { icode <<= ilen;							\
      if (fread(ipart,sizeof(uint32),1,in) != 1)			\
        { EPRINTF(EPLACE,"Could not read more bits (Decode)\n");	\
          return (1);							\
        }								\
      Flip_Long(ipart);							\
      ilen    = n-ilen;							\
      icode <<= ilen;							\
      ilen    = 32-ilen;						\
    }									\
  else									\
    { icode <<= n;							\
      ilen   -= n;							\
    }

  n     = 16;
  ilen  = 0;
  icode = 0;
  if (Flip)
    for (j = 0; j < rlen; j++)
      { GETFLIP
        c = look[*xpart];
        n = lens[c];
        if (c == signal)
          { GETFLIP
            c = *cpart;
            n = 8;
          }
        read[j] = (char) c;
      }
  else
    for (j = 0; j < rlen; j++)
      { GET
        c = look[*xpart];
        n = lens[c];
        if (c == signal)
          { GET
            c = *cpart;
            n = 8;
          }
        read[j] = (char) c;
      }

  return (0);
}

  //  Read and decode from in, the next rlen symbols into read according to non-rchar scheme
  //    neme, and the rchar runlength shceme reme

static int Decode_Run(HScheme *neme, HScheme *reme, FILE *in, char *read,
                      int rlen, int rchar)
{ int    *nlook, *nlens;
  int    *rlook, *rlens;
  int     nsignal, ilen;
  uint64  icode;
  uint32 *ipart;
  uint16 *xpart;
  uint8  *cpart;
  int     j, n, c, k;

  if (LittleEndian)
    { ipart = ((uint32 *) (&icode));
      xpart = ((uint16 *) (&icode)) + 2;
      cpart = ((uint8  *) (&icode)) + 5;
    }
  else
    { ipart = ((uint32 *) (&icode)) + 1;
      xpart = ((uint16 *) (&icode)) + 1;
      cpart = ((uint8  *) (&icode)) + 2;
    }

  if (neme->type == 2)
    nsignal = 255;
  else
    nsignal = 256;
  nlens = neme->codelens;
  nlook = neme->lookup;

  rlens = reme->codelens;
  rlook = reme->lookup;

  n     = 16;
  ilen  = 0;
  icode = 0;
  if (Flip)
    for (j = 0; j < rlen; j++)
      { GETFLIP
        c = rlook[*xpart];
        n = rlens[c];
        if (c == 255)
          { GETFLIP
            c = *xpart;
            n = 16;
          }
        for (k = 0; k < c; k++)
          read[j++] = (char) rchar;

        if (j < rlen)
          { GETFLIP
            c = nlook[*xpart];
            n = nlens[c];
            if (c == nsignal)
              { GETFLIP
                c = *cpart;
                n = 8;
              }
            read[j] = (char) c;
          }
      }
  else
    for (j = 0; j < rlen; j++)
      { GET
        c = rlook[*xpart];
        n = rlens[c];
        if (c == 255)
          { GET
            c = *xpart;
            n = 16;
          }
        for (k = 0; k < c; k++)
          read[j++] = (char) rchar;

        if (j < rlen)
          { GET
            c = nlook[*xpart];
            n = nlens[c];
            if (c == nsignal)
              { GET
                c = *cpart;
                n = 8;
              }
            read[j] = (char) c;
          }
      }

  return (0);
}


/*******************************************************************************************
 *
 *  Histogrammers
 *
 ********************************************************************************************/

//  Histogram runlengths of symbol runChar in stream[0..rlen-1] into run.

static void Histogram_Seqs(uint64 *hist, uint8 *stream, int rlen)
{ int k;

  for (k = 0; k < rlen; k++)
    hist[stream[k]] += 1;
}

static void Histogram_Runs(uint64 *run, uint8 *stream, int rlen, int runChar)
{ int k, h;

  k = 0;
  while (k < rlen)
    { h = k;
      while (k < rlen && stream[k] == runChar)
        k += 1;
      if (k-h >= 256)
        run[255] += 1;
      else
        run[k-h] += 1;
      if (k < rlen)
        k += 1;
    }
}


/*******************************************************************************************
 *
 *  Reader
 *
 ********************************************************************************************/

static char  *Read = NULL;   //  Referred by:  QVentry, Read_Lines, QVcoding_Scan,
static int    Rmax = -1;     //                Compress_Next_QVentry

static int    Nline;         //  Referred by:  QVcoding_Scan

char *QVentry()
{ return (Read); }

void Set_QV_Line(int line)
{ Nline = line; }

int Get_QV_Line()
{ return (Nline); }

//  If nlines == 1 trying to read a single header, nlines = 5 trying to read 5 QV/fasta lines
//    for a sequence.  Place line j at Read+j*Rmax and the length of every line is returned
//    unless eof occurs in which case return -1.  If any error occurs return -2.

int Read_Lines(FILE *input, int nlines)
{ int   i, rlen;
  int   tmax;
  char *tread;
  char *other;

  if (Read == NULL)
    { tmax  = MIN_BUFFER;
      tread = (char *) Malloc(5*tmax,"Allocating QV entry read buffer");
      if (tread == NULL)
        EXIT(-2);
      Rmax = tmax;
      Read = tread;
    }

  Nline += 1;
  if (fgets(Read,Rmax,input) == NULL)
    return (-1);

  rlen = strlen(Read);
  while (Read[rlen-1] != '\n')
    { tmax  = ((int) 1.4*Rmax) + MIN_BUFFER;
      tread = (char *) Realloc(Read,5*tmax,"Reallocating QV entry read buffer");
      if (tread == NULL)
        EXIT(-2);
      Rmax = tmax;
      Read = tread;
      if (fgets(Read+rlen,Rmax-rlen,input) == NULL)
        { EPRINTF(EPLACE,"Line %d: Last line does not end with a newline !\n",Nline);
          EXIT(-2);
        }
      rlen += strlen(Read+rlen);
    }
  other = Read;
  for (i = 1; i < nlines; i++)
    { other += Rmax;
      Nline += 1;
      if (fgets(other,Rmax,input) == NULL)
        { EPRINTF(EPLACE,"Line %d: incomplete last entry of .quiv file\n",Nline);
          EXIT(-2);
        }
      if (rlen != (int) strlen(other))
        { EPRINTF(EPLACE,"Line %d: Lines for an entry are not the same length\n",Nline);
          EXIT(-2);
        }
    }
  return (rlen-1);
}


/*******************************************************************************************
 *
 *  Tag compression and decompression routines
 *
 ********************************************************************************************/

//  Keep only the symbols in tags[0..rlen-1] for which qvs[k] != rchar and
//    return the # of symbols kept.

static int Pack_Tag(char *tags, char *qvs, int rlen, int rchar)
{ int j, k;

  j = 0;
  for (k = 0; k < rlen; k++)
    if (qvs[k] != rchar)
      tags[j++] = tags[k];
  tags[j] = '\0';
  return (j);
}

  //  Count the # of non-rchar symbols in qvs[0..rlen-1]

static int Packed_Length(char *qvs, int rlen, int rchar)
{ int k, clen;

  clen = 0;
  for (k = 0; k < rlen; k++)
    if (qvs[k] != rchar)
      clen += 1;
  return (clen);
}

  //  Unpack tags by moving its i'th char to position k where qvs[k] is the i'th non-rchar
  //    symbol in qvs.  All other chars are set to rchar.  rlen is the length of qvs and
  //    the unpacked result, clen is the initial length of tags.

static void Unpack_Tag(char *tags, int clen, char *qvs, int rlen, int rchar)
{ int j, k;

  j = clen-1;
  for (k = rlen-1; k >= 0; k--)
    { if (qvs[k] == rchar)
        tags[k] = 'n';
      else
        tags[k] = tags[j--];
    }
}


/*******************************************************************************************
 *
 *  Statistics Scan and Scheme creation and write
 *
 ********************************************************************************************/

  // Read up to the next num entries or until eof from the .quiva file on input and record
  //   frequency statistics.  Copy these entries to the temporary file temp if != NULL.
  //   If there is an error then -1 is returned, otherwise the number of entries read.

static uint64   delHist[256], insHist[256], mrgHist[256], subHist[256], delRun[256], subRun[256];
static uint64   totChar;
static int      delChar, subChar;

  // Referred by:  QVcoding_Scan, Create_QVcoding

void QVcoding_Scan1(int rlen, char *delQV, char *delTag, char *insQV, char *mergeQV, char *subQV)
{
  if (rlen == 0)   //  Initialization call
    { int   i;

      //  Zero histograms

      bzero(delHist,sizeof(uint64)*256);
      bzero(mrgHist,sizeof(uint64)*256);
      bzero(insHist,sizeof(uint64)*256);
      bzero(subHist,sizeof(uint64)*256);

      for (i = 0; i < 256; i++)
        delRun[i] = subRun[i] = 1;

      totChar    = 0;
      delChar    = -1;
      subChar    = -1;
      return;
    }

  //  Add streams to accumulating histograms and figure out the run chars
  //    for the deletion and substition streams

  Histogram_Seqs(delHist,(uint8 *) delQV,rlen);
  Histogram_Seqs(insHist,(uint8 *) insQV,rlen);
  Histogram_Seqs(mrgHist,(uint8 *) mergeQV,rlen);
  Histogram_Seqs(subHist,(uint8 *) subQV,rlen);

  if (delChar < 0)
    { int   k;

      for (k = 0; k < rlen; k++)
        if (delTag[k] == 'n' || delTag[k] == 'N')
          { delChar = delQV[k];
            break;
          }
    }
  if (delChar >= 0)
    Histogram_Runs( delRun,(uint8 *) delQV,rlen,delChar);
  totChar += rlen;
  if (subChar < 0)
    { if (totChar >= 100000)
        { int k;

          subChar = 0;
          for (k = 1; k < 256; k++)
            if (subHist[k] > subHist[subChar])
              subChar = k;
        }
    }
  if (subChar >= 0)
    Histogram_Runs( subRun,(uint8 *) subQV,rlen,subChar);
  return;
}

int QVcoding_Scan(FILE *input, int num, FILE *temp)
{ char *slash;
  int   rlen;
  int   i, r;

  //  Zero histograms

  bzero(delHist,sizeof(uint64)*256);
  bzero(mrgHist,sizeof(uint64)*256);
  bzero(insHist,sizeof(uint64)*256);
  bzero(subHist,sizeof(uint64)*256);

  for (i = 0; i < 256; i++)
    delRun[i] = subRun[i] = 1;

  totChar    = 0;
  delChar    = -1;
  subChar    = -1;

  //  Make a sweep through the .quiva entries, histogramming the relevant things
  //    and figuring out the run chars for the deletion and substition streams

  r = 0;
  for (i = 0; i < num; i++)
    { int well, beg, end, qv;

      rlen = Read_Lines(input,1);
      if (rlen == -2)
        EXIT(-1);
      if (rlen < 0)
        break;

      if (rlen == 0 || Read[0] != '@')
        { EPRINTF(EPLACE,"Line %d: Header in quiva file is missing\n",Nline);
          EXIT(-1);
        }
      slash = index(Read+1,'/');
      if (slash == NULL)
  	{ EPRINTF(EPLACE,"%s: Line %d: Header line incorrectly formatted ?\n",
                         Prog_Name,Nline);
          EXIT(-1);
        }
      if (sscanf(slash+1,"%d/%d_%d RQ=0.%d\n",&well,&beg,&end,&qv) != 4)
        { EPRINTF(EPLACE,"%s: Line %d: Header line incorrectly formatted ?\n",
                         Prog_Name,Nline);
          EXIT(-1);
        }

      if (temp != NULL)
        fputs(Read,temp);

      rlen = Read_Lines(input,5);
      if (rlen < 0)
        { if (rlen == -1)
            EPRINTF(EPLACE,"Line %d: incomplete last entry of .quiv file\n",Nline);
          EXIT(-1);
        }

      if (temp != NULL)
        { fputs(Read,temp);
          fputs(Read+Rmax,temp);
          fputs(Read+2*Rmax,temp);
          fputs(Read+3*Rmax,temp);
          fputs(Read+4*Rmax,temp);
        }

      Histogram_Seqs(delHist,(uint8 *) (Read),rlen);
      Histogram_Seqs(insHist,(uint8 *) (Read+2*Rmax),rlen);
      Histogram_Seqs(mrgHist,(uint8 *) (Read+3*Rmax),rlen);
      Histogram_Seqs(subHist,(uint8 *) (Read+4*Rmax),rlen);

      if (delChar < 0)
        { int   k;
          char *del = Read+Rmax;

          for (k = 0; k < rlen; k++)
            if (del[k] == 'n' || del[k] == 'N')
              { delChar = Read[k];
                break;
              }
        }
      if (delChar >= 0)
        Histogram_Runs( delRun,(uint8 *) (Read),rlen,delChar);
      totChar += rlen;
      if (subChar < 0)
        { if (totChar >= 100000)
            { int k;

              subChar = 0;
              for (k = 1; k < 256; k++)
                if (subHist[k] > subHist[subChar])
                  subChar = k;
            }
        }
      if (subChar >= 0)
        Histogram_Runs( subRun,(uint8 *) (Read+4*Rmax),rlen,subChar);

      r += 1;
    }

  return (r);
}

  //   Using the statistics in the global stat tables, create the Huffman schemes and write
  //   them to output.  If lossy is set, then create a lossy table for the insertion and merge
  //   QVs.

QVcoding *Create_QVcoding(int lossy)
{ static QVcoding coding;

  HScheme *delScheme, *insScheme, *mrgScheme, *subScheme;
  HScheme *dRunScheme, *sRunScheme;

  delScheme  = NULL;
  dRunScheme = NULL;
  insScheme  = NULL;
  mrgScheme  = NULL;
  subScheme  = NULL;
  sRunScheme = NULL;

  //  Check whether using a subtitution run char is a win

  if (totChar < 200000 || subHist[subChar] < .5*totChar)
    subChar = -1;

  //  If lossy encryption is enabled then scale insertions and merge QVs.

  if (lossy)
    { int k;

      for (k = 0; k < 256; k += 2)
        { insHist[k] += insHist[k+1];
          insHist[k+1] = 0;
        }

      for (k = 0; k < 256; k += 4)
        { mrgHist[k] += mrgHist[k+1];
          mrgHist[k] += mrgHist[k+2];
          mrgHist[k] += mrgHist[k+3];
          mrgHist[k+1] = 0;
          mrgHist[k+2] = 0;
          mrgHist[k+3] = 0;
        }
    }

  //  Build a Huffman scheme for each stream entity from the histograms

#define SCHEME_MACRO(meme,hist,label,bits)	\
  scheme = Huffman( (hist), NULL);		\
  if (scheme == NULL)				\
    goto error;					\
  if (scheme->type)				\
    { (meme) = Huffman( (hist), scheme);	\
      free(scheme);				\
    }						\
  else						\
    (meme) = scheme;

#ifdef DEBUG

#define MAKE_SCHEME(meme,hist,label,bits)	\
  SCHEME_MACRO(meme,hist,label,bits)		\
  printf("\n%s\n", (label) );			\
  Print_Histogram( (hist));			\
  Print_Table( (meme), (hist), (bits));	

#else

#define MAKE_SCHEME(meme,hist,label,bits)	\
  SCHEME_MACRO(meme,hist,label,bits)

#endif

  { HScheme *scheme;

    if (delChar < 0)
      { MAKE_SCHEME(delScheme,delHist, "Hisotgram of Deletion QVs", 8);
        dRunScheme = NULL;
      }
    else
      { delHist[delChar] = 0;
        MAKE_SCHEME(delScheme,delHist, "Hisotgram of Deletion QVs less run char", 8);
        MAKE_SCHEME(dRunScheme,delRun, "Histogram of Deletion Runs QVs", 16);
#ifdef DEBUG
        printf("\nRun char is '%c'\n",delChar);
#endif
      }

#ifdef DEBUG
    { int    k;
      uint64 count;

      count = 0;
      for (k = 0; k < 256; k++)
        count += delHist[k];
      printf("\nDelTag will require %lld bytes\n",count/4);
    }
#endif

    MAKE_SCHEME(insScheme,insHist, "Hisotgram of Insertion QVs", 8);
    MAKE_SCHEME(mrgScheme,mrgHist, "Hisotgram of Merge QVs", 8);

    if (subChar < 0)
      { MAKE_SCHEME(subScheme,subHist, "Hisotgram of Subsitution QVs", 8);
        sRunScheme = NULL;
      }
    else
      { subHist[subChar] = 0;
        MAKE_SCHEME(subScheme,subHist, "Hisotgram of Subsitution QVs less run char", 8);
        MAKE_SCHEME(sRunScheme,subRun, "Histogram of Substitution Run QVs", 16);
#ifdef DEBUG
        printf("\nRun char is '%c'\n",subChar);
#endif
      }
  }

  //  Setup endian handling

  Set_Endian(0);

  coding.delScheme  = delScheme;
  coding.insScheme  = insScheme;
  coding.mrgScheme  = mrgScheme;
  coding.subScheme  = subScheme;
  coding.dRunScheme = dRunScheme;
  coding.sRunScheme = sRunScheme;
  coding.delChar    = delChar;
  coding.subChar    = subChar;
  coding.prefix     = NULL;
  coding.flip       = 0;

  return (&coding);

error:
  if (delScheme != NULL)
    free(delScheme);
  if (dRunScheme != NULL)
    free(dRunScheme);
  if (insScheme != NULL)
    free(insScheme);
  if (mrgScheme != NULL)
    free(mrgScheme);
  if (subScheme != NULL)
    free(subScheme);
  if (sRunScheme != NULL)
    free(sRunScheme);
  EXIT(NULL);
}

  // Write the encoding scheme 'coding' to 'output'

void Write_QVcoding(FILE *output, QVcoding *coding)
{
  //   Write out the endian key, run chars, and prefix (if not NULL)

  { uint16 half;
    int    len;

    half = 0x33cc;
    fwrite(&half,sizeof(uint16),1,output);

    if (coding->delChar < 0)
      half = 256;
    else
      half = (uint16) (coding->delChar);
    fwrite(&half,sizeof(uint16),1,output);

    if (coding->subChar < 0)
      half = 256;
    else
      half = (uint16) (coding->subChar);
    fwrite(&half,sizeof(uint16),1,output);

    len = strlen(coding->prefix);
    fwrite(&len,sizeof(int),1,output);
    fwrite(coding->prefix,1,len,output);
  }

  //   Write out the scheme tables

  Write_Scheme(coding->delScheme,output);
  if (coding->delChar >= 0)
    Write_Scheme(coding->dRunScheme,output);
  Write_Scheme(coding->insScheme,output);
  Write_Scheme(coding->mrgScheme,output);
  Write_Scheme(coding->subScheme,output);
  if (coding->subChar >= 0)
    Write_Scheme(coding->sRunScheme,output);
}

  // Read the encoding scheme 'coding' to 'output'

QVcoding *Read_QVcoding(FILE *input)
{ static QVcoding coding;

  // Read endian key, run chars, and short name common to all headers

  { uint16 half;
    int    len;

    if (fread(&half,sizeof(uint16),1,input) != 1)
      { EPRINTF(EPLACE,"Could not read flip byte (Read_QVcoding)\n");
        EXIT(NULL);
      }
    coding.flip = (half != 0x33cc);

    if (fread(&half,sizeof(uint16),1,input) != 1)
      { EPRINTF(EPLACE,"Could not read deletion char (Read_QVcoding)\n");
        EXIT(NULL);
      }
    if (coding.flip)
      Flip_Short(&half);
    coding.delChar = half;
    if (coding.delChar >= 256)
      coding.delChar = -1;

    if (fread(&half,sizeof(uint16),1,input) != 1)
      { EPRINTF(EPLACE,"Could not read substitution char (Read_QVcoding)\n");
        EXIT(NULL);
      }
    if (coding.flip)
      Flip_Short(&half);
    coding.subChar = half;
    if (coding.subChar >= 256)
      coding.subChar = -1;

    //  Read the short name common to all headers

    if (fread(&len,sizeof(int),1,input) != 1)
      { EPRINTF(EPLACE,"Could not read header name length (Read_QVcoding)\n");
        EXIT(NULL);
      }
    if (coding.flip)
      Flip_Long(&len);
    coding.prefix = (char *) Malloc(len+1,"Allocating header prefix");
    if (coding.prefix == NULL)
      EXIT(NULL);
    if (len > 0)
      { if (fread(coding.prefix,len,1,input) != 1)
          { EPRINTF(EPLACE,"Could not read header name (Read_QVcoding)\n");
            EXIT(NULL);
          }
      }
    coding.prefix[len] = '\0';
  }

  //  Setup endian handling

  Set_Endian(coding.flip);

  //  Read the Huffman schemes used to compress the data

  coding.delScheme  = NULL;
  coding.dRunScheme = NULL;
  coding.insScheme  = NULL;
  coding.mrgScheme  = NULL;
  coding.subScheme  = NULL;
  coding.sRunScheme = NULL;

  coding.delScheme = Read_Scheme(input);
  if (coding.delScheme == NULL)
    goto error;
  if (coding.delChar >= 0)
    { coding.dRunScheme = Read_Scheme(input);
      if (coding.dRunScheme == NULL)
        goto error;
    }
  coding.insScheme = Read_Scheme(input);
  if (coding.insScheme == NULL)
    goto error;
  coding.mrgScheme = Read_Scheme(input);
  if (coding.mrgScheme == NULL)
    goto error;
  coding.subScheme = Read_Scheme(input);
  if (coding.subScheme == NULL)
    goto error;
  if (coding.subChar >= 0)
    { coding.sRunScheme = Read_Scheme(input);
      if (coding.sRunScheme == NULL)
        goto error;
    }

  return (&coding);

error:
  if (coding.delScheme != NULL)
    free(coding.delScheme);
  if (coding.dRunScheme != NULL)
    free(coding.dRunScheme);
  if (coding.insScheme != NULL)
    free(coding.insScheme);
  if (coding.mrgScheme != NULL)
    free(coding.mrgScheme);
  if (coding.subScheme != NULL)
    free(coding.subScheme);
  if (coding.sRunScheme != NULL)
    free(coding.sRunScheme);
  EXIT(NULL);
}

  //  Free all the auxilliary storage associated with the encoding argument

void Free_QVcoding(QVcoding *coding)
{ if (coding->subChar >= 0)
    free(coding->sRunScheme);
  free(coding->subScheme);
  free(coding->mrgScheme);
  free(coding->insScheme);
  if (coding->delChar >= 0)
    free(coding->dRunScheme);
  free(coding->delScheme);
  free(coding->prefix);
}


/*******************************************************************************************
 *
 *  Encode/Decode (w.r.t. coding) next entry from input and write to output
 *
 ********************************************************************************************/

void Compress_Next_QVentry1(int rlen, char *del, char *tag, char *ins, char *mrg, char *sub,
                            FILE *output, QVcoding *coding, int lossy)
{ int clen;

  if (coding->delChar < 0)
    { Encode(coding->delScheme, output, (uint8 *) del, rlen);
      clen = rlen;
    }
  else
    { Encode_Run(coding->delScheme, coding->dRunScheme, output,
                 (uint8 *) del, rlen, coding->delChar);
      clen = Pack_Tag(tag,del,rlen,coding->delChar);
    }
  Number_Read(tag);
  Compress_Read(clen,tag);
  fwrite(tag,1,COMPRESSED_LEN(clen),output);

  if (lossy)
    { uint8 *insert = (uint8 *) ins;
      uint8 *merge  = (uint8 *) mrg;
      int    k;

      for (k = 0; k < rlen; k++)
        { insert[k] = (uint8) ((insert[k] >> 1) << 1);
          merge[k]  = (uint8) (( merge[k] >> 2) << 2);
        }
    }

  Encode(coding->insScheme, output, (uint8 *) ins, rlen);
  Encode(coding->mrgScheme, output, (uint8 *) mrg, rlen);
  if (coding->subChar < 0)
    Encode(coding->subScheme, output, (uint8 *) sub, rlen);
  else
    Encode_Run(coding->subScheme, coding->sRunScheme, output,
               (uint8 *) sub, rlen, coding->subChar);
  return;
}

int Compress_Next_QVentry(FILE *input, FILE *output, QVcoding *coding, int lossy)
{ int rlen, clen;

  //  Get all 5 streams, compress each with its scheme, and output

  rlen = Read_Lines(input,5);
  if (rlen < 0)
    { if (rlen == -1)
        EPRINTF(EPLACE,"Line %d: incomplete last entry of .quiv file\n",Nline);
      EXIT (-1);
    }

  if (coding->delChar < 0)
    { Encode(coding->delScheme, output, (uint8 *) Read, rlen);
      clen = rlen;
    }
  else
    { Encode_Run(coding->delScheme, coding->dRunScheme, output,
                 (uint8 *) Read, rlen, coding->delChar);
      clen = Pack_Tag(Read+Rmax,Read,rlen,coding->delChar);
    }
  Number_Read(Read+Rmax);
  Compress_Read(clen,Read+Rmax);
  fwrite(Read+Rmax,1,COMPRESSED_LEN(clen),output);

  if (lossy)
    { uint8 *insert = (uint8 *) (Read+2*Rmax);
      uint8 *merge  = (uint8 *) (Read+3*Rmax);
      int    k;

      for (k = 0; k < rlen; k++)
        { insert[k] = (uint8) ((insert[k] >> 1) << 1);
          merge[k]  = (uint8) (( merge[k] >> 2) << 2);
        }
    }

  Encode(coding->insScheme, output, (uint8 *) (Read+2*Rmax), rlen);
  Encode(coding->mrgScheme, output, (uint8 *) (Read+3*Rmax), rlen);
  if (coding->subChar < 0)
    Encode(coding->subScheme, output, (uint8 *) (Read+4*Rmax), rlen);
  else
    Encode_Run(coding->subScheme, coding->sRunScheme, output,
               (uint8 *) (Read+4*Rmax), rlen, coding->subChar);

  return (rlen);
}

int Uncompress_Next_QVentry(FILE *input, char **entry, QVcoding *coding, int rlen)
{ int clen, tlen;

  //  Decode each stream and write to output

  if (coding->delChar < 0)
    { if (Decode(coding->delScheme, input, entry[0], rlen))
        EXIT(1);
      clen = rlen;
      tlen = COMPRESSED_LEN(clen);
      if (tlen > 0)
        { if (fread(entry[1],tlen,1,input) != 1)
            { EPRINTF(EPLACE,"Could not read deletions entry (Uncompress_Next_QVentry\n");
              EXIT(1);
            }
        }
      Uncompress_Read(clen,entry[1]);
      Lower_Read(entry[1]);
    }
  else
    { if (Decode_Run(coding->delScheme, coding->dRunScheme, input,
                     entry[0], rlen, coding->delChar))
        EXIT(1);
      clen = Packed_Length(entry[0],rlen,coding->delChar);
      tlen = COMPRESSED_LEN(clen);
      if (tlen > 0)
        { if (fread(entry[1],tlen,1,input) != 1)
            { EPRINTF(EPLACE,"Could not read deletions entry (Uncompress_Next_QVentry\n");
              EXIT(1);
            }
        }
      Uncompress_Read(clen,entry[1]);
      Lower_Read(entry[1]);
      Unpack_Tag(entry[1],clen,entry[0],rlen,coding->delChar);
    }

  if (Decode(coding->insScheme, input, entry[2], rlen))
    EXIT(1);

  if (Decode(coding->mrgScheme, input, entry[3], rlen))
    EXIT(1);

  if (coding->subChar < 0)
    { if (Decode(coding->subScheme, input, entry[4], rlen))
        EXIT(1);
    }
  else
    { if (Decode_Run(coding->subScheme, coding->sRunScheme, input,
                     entry[4], rlen, coding->subChar))
        EXIT(1);
    }

  return (0);
}
