# FastGA: A Fast Genome Aligner
  
<font size ="4">**_Author:  Gene Myers_**<br>
**_First:   May 10, 2023_**<br>

FastGA searches for all local DNA alignments between two high quality genomes.
The core assumption is that the genomes are nearly complete involving at most several
thousand contigs with a sequence quality of Q40 or better.
Based on a novel adaptive seed finding algorithm and the wave-based local aligner developed for
[daligner](https://github.com/thegenemyers/DALIGNER), the tool can for example compare
two 2Gbp genomes finding almost all regions over 100bp that are 70% or more similar
in about 4 minutes wall clock time on my MacPro with 8 cores (about 24 CPU minutes).
Moreover, it uses a trace point concept to record all the found alignments in a very space-efficient manner, e.g. just 50MB for over
640,000 local alignments (la's) in our running example.

FastGA uses several components of the Daligner assembler software suite,
namely its concept of a genome "database" and its tools for encoding
and manipulating local alignments.  So in addition to installing this
repository, you need to download the repositories and install
[DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) and
[DALIGNER](https://github.com/thegenemyers/DALIGNER).
No need to panic, simply download and say "make" and that's it.
You should then place the produced binaries on your execution path so
that you then have fasta2DAM, DBsplit, DAM2fasta, DBstats, and DBshow available for
creating and viewing a genome database made from the standard fasta
encoding, and the tool LAshow for viewing your alignments.  In
addition, there are programs therein, [DB2ONE](https://github.com/thegenemyers/DAZZ_DB)
and [LA2ONE](https://github.com/thegenemyers/DALIGNER) to convert the genome database and the
local alignments into [ONEcode](https://github.com/thegenemyers/ONEcode) formated files which are trivial to read and analyze.

FastGA also uses the data encodings used by the [FastK](https://github.com/thegenemyers/FASTK)
k-mer counter suite
for a k-mer table, so in principle one could use tools from FastK to manipulate
the .ktab portion of a genome index produced by **Gindex**.  But unlike the
Daligner components this is not necessary.

The overall workflow for using **FastGA**, assuming you are starting from
say NCBI fasta files of the target genomes, is to first convert each of them into a Dazzler database with [fasta2DAM](https://github.com/thegenemyers/DAZZ_DB) and
[DBsplit](https://github.com/thegenemyers/DAZZ_DB), then build an index of each of them with **Gindex**, and then you are ready to compare pairs of
genomes with **FastGA**.
While you may at first think it a little inconvenient to have to build
a Dazzler database, its actually a good thing as (a) the database is
25% the size of the fasta file, (b) you can delete the fasta file as it can be reconstructed *exactly* from the database with [DAM2fasta](https://github.com/thegenemyers/DAZZ_DB), and (c) you have a tool [DBstats](https://github.com/thegenemyers/DAZZ_DB) to give you summary statistics and another, [DBshow](https://github.com/thegenemyers/DAZZ_DB), to randomly access
and view the contigs of the genome.  Creating the data base for a genome in file say ```mygenome.fasta``` is as simple as issuing the two commands:

```
   fasta2DAM -v mydb mygenome.fasta
   DBsplit mydb
```

The two command sequence above results in the creation of the visible file ```mydb.dam``` and 3 hidden
files ```.mydb.idx```, ```.mydb.hdr```, and ```.mydb.bps```.  The database can be deleted
with ```DBrm  mydb``` which removes the visible file and the 3 hidden files.

The remainder of this document gives the command line syntax and complete technical description
for the tools **Gindex** that creates an index, **Gshow** that prints a listing of an index, and **FastGA** which compares two genomes and outputs a local alignments (.las) file of all the matches it finds above user-supplied thresholds for length and similarity.

<a name="Gindex"></a>

```
1. Gindex [-v] [-T<int(8)>] -f<int> <source>.[dam]
```

Given a Dazzler database, or .dam, for a genome, **Gindex** produces an index of the genome
consisting of a pair of files \<source>.ktab and \<source>.post.  The .ktab file consists of a
FastK table of all the 40-mers in the genome that occur -f or few times.  The .post file
contains successive lists of the genome positions at which the 40-mers in the first table occur
in order of the 40-mers in the table.  The option -f is designed to remove 40-mers from
repetitive regions of the genome from consideration.  It must be specified, we suggest a good
default value is say 10.  With the -v
option the program prints a progress report to the standard output.

The program runs with -T threads, 8 by default and for which it typically sees a 6X or more
speed up in wall clock time.  **Especially note that every index to be compared must be built
with the same number of threads, and this number of threads will be used by FastGA**.

Both the .ktab and .post files refer to -T<sup>2</sup> hidden files with the names ```.<source>.ktab.#```
and ```.<source>.post.#``` where # varies from 1 to T<sup>2</sup>.  Altogether these files typically occupy 8-9GB per Gbp of genome for
the .ktab and 5-6GB per Gbp for the .post file.  So they are rather large, e.g. ~42GB for
a human genome.

<a name="Gshow"></a>

```
2. Gshow <source>.[dam]
```

Simply print to stdout a representation of the index produced by **Gindex**.  Mostly for
debug and illustrative purposes.

```
3. FastGA [-v] [-P<dir(/tmp)] [-o<out:name>] -f<int>
          [-c<int(100)>] [-s<int(500)>] [-a<int(100)>] [-e<float(.7)>]
          <source1>[.dam] <source2>[.dam]
```

Once indices have been built for all relevant genomes, any pair \<source1> and \<source2> can be compared with **FastGA**.  The indices for the genomes must have been built with the same
number of threads and this is the number of threads that FastGA uses to perform the comparison.
The output of the program is a Daligner .las (local alignments) file whose alignments can then
be viewed with LAshow.  If the -o option is given then the name of the .las file is
```<out>.las```, otherwise the name is ```<source1>.<source2>.las``` by default.  FastGA
produces a number of large intermediate files in the directory ```/tmp``` unless an alternate
directory is given with -P parameter.

The traditional definition of the *adaptive seed* at a given position p of source1, is the longest string beginning at that position that is also somewhere in source2.  If the number of
occurences of this string in source2 is greater than -f, then the adaptamer is deemed repetitive and is not considered.  Otherwise an adpatamer seed hit occurs at p and the positions of source2
where the adapamer also occurs.  The option -f must be supplied and it must be less than or equal to the value of -f specified for the index construction.

We practically limit the length of the longest adaptamer to 40bp, the length of the k-mers in
the indices.  The original adaptamer algorithm built a suffix tree or array of source2 and then
in a scan of source1 determined all the adaptamers and their hits.  This is slow and not cache
coherent.  FastGA finds adaptamer seeds and their hits in a novel and highly efficient sequential merge of the two indices, that spends less time on lookup and is cache coherent.  

FastGA then searches for runs or chains of adaptamer seed hits that (a) all lie within a diagonal band of width 128, (b) the spacing between every two consecutive seeds is less than -s(500), and
(c) the seeds in the chain cover at least -c(100) bases in both genomes.  For these chain
hits, FastGA then runs a wave-based local alignment routine that searches for a local alignment
of length at least -a(100) with a similarity of -e(70%) or better that contains at least one
of the seeds in the chain.  All such found alignments are reported in the output .las file in
lexicgraphical order of source1 contig #, and then source 2 contig #, and then the start coordinate of the alignment in source1.  The options -s, -c -a, and -e can be used to modify the default
thresholds for chaining and alignment.