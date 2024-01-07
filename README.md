# FastGA: A Fast Genome Aligner
  
<font size ="4">**_Author:  Gene Myers_**<br>
**_First:   May 10, 2023_**<br>

FastGA searches for all local DNA alignments between two high quality genomes.
The core assumption is that the genomes are nearly complete involving at most several
thousand contigs with a sequence quality of Q40 or better.
Based on a novel adaptive seed finding algorithm and the wave-based local aligner developed for
[daligner](https://github.com/thegenemyers/DALIGNER), the tool can for example compare
two 2Gbp bat genomes finding almost all regions over 100bp that are 70% or more similar
in about 4 minutes wall clock time on my MacPro with 8 cores (about 26 CPU minutes).
Moreover, it uses a trace point concept to record all the found alignments in a very space-efficient manner, e.g. just 57MB for over
759,000 local alignments (la's) in our running example.
The alignments can also be translated into .psl format, .paf format, or 1-code format
with programs described below.

FastGA uses the Daligner assembler's
"genome database" encoding of a genome or assembly.  So in addition to installing this
repository, you need to download and install the repository
[DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB).
No need to panic, simply download and say "make" and that's it.
You should then place the produced binaries on your execution path so
that you then have fasta2DAM, DAM2fasta, DBstats, and DBshow available for
creating and viewing a genome database made from the standard fasta
encoding of an assembly/genome.  In addition, the program fasta2DAM converts a database back
into the original .fasta file from which it was built

FastGA outputs all the local alignments it finds in the Daligner assembler's .las binary
format which very compactly encode alignments.  If one is comfortable with or wishes to work directly with this format, then one should download and install the repository
[DALIGNER](https://github.com/thegenemyers/DALIGNER)
so that they have the tools LAshow, LAsort, and LAmerge for manipulating and
viewing alignments in .las format.  But most user's will likely want to work with a more
popular, albeit voluminous, format such a .psl or .paf.  To this end we provide here
LAStoPSL that converts our .las files into .psl, and LAStoPAF that converts
them into .paf.

For [ONEcode](https://github.com/thegenemyers/ONEcode) enthusiasts there is a program
[DB2ONE](https://github.com/thegenemyers/DAZZ_DB) in the Dazzler repository that
converts a genome database into a 1-code file, and another
[LA2ONE](https://github.com/thegenemyers/DALIGNER) in the Daligner repository that
converts a .las file into a 1-code file.  ONEcode is a general data encoding system
that supports an easy to read and interpret ASCII format coupled with a production-capable
binary form that includes compression.

The overall workflow for using **FastGA**, assuming you are starting from
say NCBI fasta files of the target genomes, is to first convert each of them into a Dazzler database with [fasta2DAM](https://github.com/thegenemyers/DAZZ_DB), then build an index of each of them with **Gindex**, and then you are ready to compare pairs of
genomes with **FastGA**.
While you may at first think it a little inconvenient to have to build
a Dazzler database, its actually a good thing as (a) the database is
25% the size of the fasta file, (b) you can delete the fasta file as it can be reconstructed *exactly* from the database with [DAM2fasta](https://github.com/thegenemyers/DAZZ_DB), and (c) you have a tool [DBstats](https://github.com/thegenemyers/DAZZ_DB) to give you summary statistics and another, [DBshow](https://github.com/thegenemyers/DAZZ_DB), to randomly access
and view the contigs of the genome.  Creating the data base for a genome in file say ```mygenome.fasta``` is as simple as issuing the command:

```
   fasta2DAM -v mydb mygenome.fasta
```

This command results in the creation of the visible file ```mydb.dam``` and 3 hidden
files ```.mydb.idx```, ```.mydb.hdr```, and ```.mydb.bps```.  The database can be deleted
with ```DBrm  mydb``` which removes the visible file and the 3 hidden files.

In the directory EXAMPLE of this repository you will find a text script file showing a complete
example building indices for two genomes, comparing them, and then output the alignments in
three different forms including .psl format.  The directory also contains all or, if too large
for github, a portion, of the three outputs.

The remainder of this document gives the command line syntax and complete technical description
for the tools **Gindex** that creates an index, **Gshow** that prints a listing of an index, **FastGA** which compares two genomes and outputs a local alignments (.las) file of all the matches it finds above user-supplied thresholds for length and similarity, **LAStoPSL** that converts
a .las file into .psl format, and **LAStoPAF** that converts a .las file into .paf format.

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

The program runs with -T threads, 8 by default, for which it typically sees a 6X or more
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
          <source1>[.dam] [<source2>[.dam]]
```

Once indices have been built for all relevant genomes, any pair \<source1> and \<source2> can be compared with **FastGA**.  **The indices for the genomes must have been built with the same
number of threads and this is the number of threads that FastGA uses to perform the comparison.**
The output of the program is a Daligner .las (local alignments) file whose alignments can then
be viewed with LAshow.  If the -o option is given then the name of the .las file is
```<out>.las```, otherwise the name is ```<source1>.<source2>.las``` by default.  FastGA
produces a number of large intermediate files in the directory ```/tmp``` unless an alternate
directory is given with -P parameter.

You can also call FastGA with a single source in which case it compares the genome/assembly against
itself, avoiding the finding and reporting of identity matches.  We anticipate that this mode may be usefule for (a) identifying repetitive regions in a genome, and (b) resolving and identifying haplotype correspondences between the various contigs of an assembly. 

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


```
4. LAStoPSL [-T<int(8)>] <source1>[.dam] [<source2>[.dam]] <alignments>[.las]
```

In order to convert a .las file into a [.psl](https://www.ensembl.org/info/website/upload/psl.html) file, LAStoPSL also needs as arguments the one or two genomes,
as Dazzler .dam's, that were compared by FastGA to produce the file.  Given these 2 or 3 arguments in the
order shown above, the program outputs .psl encoded alignments, one alignment per line to the
standard output.
Warning, the .psl output is almost 15 times larger than the .las file.
LAStoPSL uses 8 threads by default, but this can be changed with the -T option.



```
5. LAStoPAF [-mxtg] [-T<int(8)>] <source1>[.dam] [<source2>[.dam]] <alignments>[.las]
```

In order to convert a .las file into a [.paf](https://github.com/lh3/miniasm/blob/master/PAF.md) file, LAStoPAF also needs as arguments the one or two genomes,
as Dazzler .dam's, that were compared by FastGA to produce the file.  Given these 2 or 3 arguments in the
order shown above, the program outputs .paf encoded records, one alignment per line to the
standard output.

In addition to the standard .paf fields, LAStoPAF outputs a ```dv:F<fraction>``` SAM-tag that gives the divergences of the query from the target and a ```df:I<diffs>``` SAM-tag that gives the number
of differences in an optimal alignment between the two intervals of the query and target.

If the -t option is set then the program also outputs a ```tz:Z<length-list>``` tag and a
```td:Z<diff-list>``` tag encoding the .las trace information from which an explicit alignment
can be efficiently
constructed.  Perhaps more relevant to most users are the -m and -x options that requests LAStoPAF to produce a CIGAR string tag of the form ```cg:Z<cigar-string>``` that explicitly encodes an optimal alignment.  With the -m option, aligned characters regardless of whether they are equal
or not are
encoded with an 'M'.  With the -x option, aligned *equal* characters are encoded with an '=' and
aligned *unequal* characters with an 'X'.
The -t option doubles the time taken and quadruples the size of the .paf file.  Beware, the -m and -x options
increase the time taken by a factor of 10 and the file size by a factor of almost 100 !  This can
be ameliorated somewhat by running LAStoPAF with more threads, controllable with the -T option.

The -g option is an interim flag that allows one to turn on gap reduction when producing the alignments required to effect the -m or -x options.  Alignment computed under simple Levenstein
tend to smear characters of one sequence across a large insert in the other sequence.  Affine
gap costs avoid this



