/*****************************************************************************************\
*                                                                                         *
*  Hash Table data abstraction.                                                           *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  March 2006                                                                    *
*                                                                                         *
\*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hash.h"
#include "gene_core.h"

/* Hash Table cell or entry, index is implicitly given by position in table->cells array.  */

typedef struct
{ int      next;  // hash bucket link
  char    *text;  // string entry
} Entry;

typedef struct
{ int         veclen;  // size of hash vector
  int         cntmax;  // current max # of entries
  int         count;   // number of entries in hash table
  int         strmax;  // current max of string array
  int         strtop;  // current top of string array
  int        *vector;  // hash vector
  Entry      *cells;   // array where hash cells are allocated
  char       *strings; // array of entry strings
} Table;

#define CELL_RATIO    .4  // ratio of max # of cells to hash vector length
#define STRING_RATIO   6  // expected average entry length (including terminating 0-byte)

#define T(x) ((Table *) x)

void Free_Hash_Table(Hash_Table *hash_table)
{ Table *table = T(hash_table);

  free(table->cells);
  free(table);
}

/* Hash key for a string is xor of each consecutive 3 bytes. */

static int hash_key(char *entry)
{ int i, key, glob;

  key = 0;
  glob = 0;
  for (i = 0; entry[i] != '\0'; i++)
    { glob = (glob << 8) | entry[i];
      if (i % 3 == 2)
        { key = key ^ glob;
          glob = 0;
        }
    }
  if (i % 3 != 0)
    key = key ^ glob;
  return (key);
}

/*  Find the next prime larger than size.  Variant of Sieve of Arosthenes:  First
      time its called computes all primes between 2 and 0xFFFF using the basic
      sieve algorithm.  With these one builds a sieve for the 0xFFFF numbers from
      size upwards, using the primes to x-out sieve elements as being non-prime.
      This will work up to 0x7FFFFFF, beyond the largest positive integer, because 
      it suffices to sieve against the square root of the largest number in the sieve.   */

static int next_prime(int size)
{ static int           firstime = 1;
  static int           Prime[0x4000], Ptop;
  static unsigned char Sieve[0x10000];
  int p, q, n;

  if (firstime)
    { firstime = 0;

      Ptop = 0;
      for (p = 2; p < 0x10000; p++)
        Sieve[p] = 1;
      for (p = 2; p < 0x10000; p++)
        if (Sieve[p])
          { for (q = 2*p; q < 0x10000; q += p)
              Sieve[q] = 0;
            Prime[Ptop++] = p;
          }
    }

  while (size < 0x7FFF0000)
    { for (q = 0; q < 0x10000; q++)
        Sieve[q] = 1;
      for (p = 0; p < Ptop; p++)
        { n = Prime[p];
          if (n >= size) break;
          for (q = ((size-1)/n+1)*n - size; q < 0x10000; q += n)
            Sieve[q] = 0;
        }
      for (q = 0; q < 0x10000; q++)
        if (Sieve[q])
          return (size+q);
      size += 0x10000;
    }

  return (size);
}


/* Diagnostic output of hash table contents. */

void Print_Hash_Table(FILE *file, Hash_Table *hash_table)
{ Table      *table  = T(hash_table);
  int        *vector = table->vector;
  Entry      *cells  = table->cells;
  int         i, c;

  fprintf(file,"\nHASH TABLE %d/%d %d",table->count,table->cntmax,table->veclen);
  if (table->strmax > 0)
    fprintf(file," %d/%d",table->strtop,table->strmax);
  fprintf(file,"\n");
  for (i = 0; i < table->veclen; i++)
    if ((c = vector[i]) >= 0)
      { fprintf(file,"  Vector %4d:\n",i);
        for (; c >= 0; c = cells[c].next)
          fprintf(file,"    %4d: '%s'\n",c,cells[c].text);
      }
}

Hash_Table *New_Hash_Table(int size, int keep)
{ Table *table;
  void  *room;
  int    vlen, smax;
  int    i;

  if (size <= 0)
    { fprintf(stderr,"%s: Table must have > 0 entries (New_Hash_Table)\n",Prog_Name);
      EXIT (NULL);
    }

  vlen = next_prime((int) (size/CELL_RATIO));
  smax = size*STRING_RATIO;

  table = Malloc(sizeof(Table),"Allocating hash table");
  if (keep)
    room = Malloc(size*sizeof(Entry)+vlen*sizeof(int)+smax,"Allocating hash table");
  else
    room = Malloc(size*sizeof(Entry)+vlen*sizeof(int),"Allocating hash table");
  if (table == NULL || room == NULL)
    { if (table != NULL)
        free(table);
      if (room != NULL)
        free (room);
      EXIT (NULL);
    }

  table->cells = (Entry *) room;
  room += size*sizeof(Entry);
  table->vector = (int *) room;

  table->cntmax = size;
  table->veclen = vlen;

  table->count  = 0;
  for (i = 0; i < vlen; i++)
    table->vector[i] = -1;

  if (keep)
    { room += vlen*sizeof(int);
      table->strings = (char *) room;
      table->strmax  = smax;
      table->strtop  = 0;
    }
  else
    table->strmax = 0;

  return ((Hash_Table *) table);
}

/* Double the size of a hash table
   while preserving its contents.    */

static int double_hash_table(Table *table)
{ int   size, vlen, smax;
  void *room;

  size = 2*table->cntmax;
  vlen = next_prime((int) (size/CELL_RATIO));
  smax = (int) (2.1 * table->strtop + 1000);

  if (table->strmax > 0)
    room = Realloc(table->cells,size*sizeof(Entry)+vlen*sizeof(int)+smax,"Expanding hash table");
  else
    room = Realloc(table->cells,size*sizeof(Entry)+vlen*sizeof(int),"Expanding hash table");
  if (room == NULL)
    return (1);

  table->cells = (Entry *) room;
  room += size*sizeof(Entry);
  table->vector = (int *) room;

  table->cntmax = size;
  table->veclen = vlen;

  if (table->strmax > 0)
    { room += vlen*sizeof(int);
      memmove(room,table->strings,table->strtop);
      table->strings = (char *) room;
      table->strmax  = smax;
    }
 
  { int   *vector = table->vector;
    Entry *cells  = table->cells;
    int    c, key;

    for (c = 0; c < vlen; c++)
      vector[c] = -1;
    for (c = 0; c < table->count; c++)
      { key = hash_key(cells[c].text) % size;
        cells[c].next = vector[key];
        vector[key] = c;
      }
  }

  return (0);
}

/* Lookup string 'entry' in table 'table' and return its
   unique nonnegative id, or -1 if it is not in the table. */

int Hash_Lookup(Hash_Table *hash_table, char *entry)
{ Table *table = T(hash_table);
  int    key, chain;

  key   = hash_key(entry) % table->veclen;
  chain = table->vector[key];
  while (chain >= 0)
    { if (strcmp(table->cells[chain].text,entry) == 0)
        return (chain);
      chain = table->cells[chain].next;
    }
  return (-1);
}

/* Add string 'entry' in table 'table' and return its assigned
   uniqe nonnegative id.  Return -1 if an error occurs in INTERACTIVE mode. */

int Hash_Add(Hash_Table *hash_table, char *entry)
{ Table *table = T(hash_table);
  void  *room;
  int    smax, vlen, size;
  int    key, chain, len;

  key   = hash_key(entry) % table->veclen;
  chain = table->vector[key];
  while (chain >= 0)
    { if (strcmp(table->cells[chain].text,entry) == 0)
        return (chain);
      chain = table->cells[chain].next;
    }

  if (table->count+1 > table->cntmax)
    { if (double_hash_table(table))
        EXIT (-1);
      key = hash_key(entry) % table->veclen;
    }

  chain = table->count;
  table->cells[chain].next = table->vector[key];
  table->vector[key] = chain;

  if (table->strmax == 0)
    { table->cells[chain].text = entry;
      return (table->count++);
    }

  len = (int) (strlen(entry) + 1);
  if (table->strtop + len > table->strmax)
    { smax = ((table->strtop+len)*1.1*table->cntmax) / (table->count+1) + 1000;
      vlen = table->veclen;
      size = table->cntmax;

      room = Realloc(table->cells,size*sizeof(Entry)+vlen*sizeof(int)+smax,"Expanding hash table");
      if (room == NULL)
        EXIT (-1);

      table->cells = (Entry *) room;
      room += size*sizeof(Entry);
      table->vector = (int *) room;
      room += vlen*sizeof(int);
      table->strings = room;
      table->strmax = smax;
    }
  strcpy(table->strings + table->strtop, entry);
  table->cells[chain].text = table->strings + table->strtop;
  table->strtop += len;
  return (table->count++);
}

/* Return the current # of entries in the hash table. */

int Get_Hash_Size(Hash_Table *hash_table)
{ return (T(hash_table)->count); }

/* Return the string with unique id i in table. */

char *Get_Hash_String(Hash_Table *hash_table, int i)
{ return (T(hash_table)->cells[i].text); }

/* Clear the contents of hash table, reseting it to be empty. */

void Clear_Hash_Table(Hash_Table *hash_table)
{ Table *table = T(hash_table);
  int    i;

  table->count  = 0;
  table->strtop = 0;
  for (i = 0; i < table->veclen; i++)
    table->vector[i] = -1;
}
