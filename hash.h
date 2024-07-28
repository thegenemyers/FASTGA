/*****************************************************************************************\
*                                                                                         *
*  Hash Table data abstraction.                                                           *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  March 2006                                                                    *
*                                                                                         *
\*****************************************************************************************/

#ifndef _HASH_TABLE

#define _HASH_TABLE

typedef void Hash_Table;

  //  New_Hash_Table:
  //    Create a hash table to accommodate "size" entries, initially.  If more than
  //    this # of entries is added than the hash table will expand by progressive
  //    doubling to accommodate.  If "keep" is non-zero than the hash table records
  //    an entry's string, otherwise it will simply point at the string supplied by
  //    the user.  NULL is return if there is an error in INTERACTIVE mode.
  //  Clear_Hash_Table:
  //    Empty the contents of the hash table (but do not destroy it)
  //  Free_Hash_Table:
  //    Free all space used by the table.

Hash_Table *New_Hash_Table(int size, int keep);
void        Clear_Hash_Table(Hash_Table *table);
void        Free_Hash_Table(Hash_Table *table);

  //  Hash_Lookup:
  //    Lookup "entry" in the hash table.  If found return its index, otherwise return -1
  //  Hash_Add:
  //    Add "entry" to the hash table and return its assigned index
  //    If there is an error then -1 is return in INTERACTIVE mode.

int Hash_Lookup(Hash_Table *table, char *entry);
int Hash_Add(Hash_Table *table, char *entry);

  //  Get_Hash_Size:
  //    Return the # of entries currently in the hash table
  //  Get_Hash_String:
  //    Return a pointer to the string for the entry with index "i"

int   Get_Hash_Size(Hash_Table *table);
char *Get_Hash_String(Hash_Table *table, int i);

  //  Print_Hash_Table:
  //    Print out the contents and structure of a hash table (for debug)

void Print_Hash_Table(FILE *file, Hash_Table *table);

#endif
