# FastGA: A Fast Genome Aligner
  
<font size ="4">**_Author:  Gene Myers_**<br>
**_First:   May 10, 2023_**<br>
**_Last:  Feb. 1, 2024_**<br>

- [FastGA](#FastGA)

- [Sub-Process Routines](#subprocess)
  - [FAtoGDB](#FAtoGDB): Convert a FASTA file into a genome database (GDB)
  - [GIXmake](#GIXmake): Build a genome index (GIX) for a given GDB
  - [LAStoPAF](#LAStoPAF): Stream PAF formatted alignments for a given LAS file
  - [LAStoPSL](#LAStoPSL): Stream PSL formatted alignments for a given LAS file
  - [LAStoONE](#LAStoONE): Stream a binary 1-code file for a given LAS file

- [Viewing Utilities](#viewing)
  - [GDBshow](#GDBshow): Display select contigs or substrings thereof from a GDB
  - [GDBstat](#GDBstat): Display various statistics and histograms of the scaffolds & contigs in a GDB
  - [GIXshow](#GIXshow): Display range of a GIX
  - [LASshow](#LASshow): Display selected alignments in an LAS in a variety of forms

- [Additional Utilities](#addons)
  - [GDBtoFA](#GDBtoFA): Converts a GDB back to exactly the FASTA file it was derived from
  - [GIXrm](#GIXrm): Remove GDBs and GIXs including their hidden parts
  - [GIXcp](#GIXcp): Copy GDBs and GIXs including their hidden parts as an ensemble
  - [GIXmv](#GIXmv): Move GDBs and GIXs including their hidden parts as an ensemble
  - [LASreset](#LASreset): Reset an LAS's internal references to the GDB(s) it was computed from


## Overview

**FastGA** searches for all local DNA alignments between two high quality genomes.
The core assumption is that the genomes are nearly complete involving at most several
thousand contigs with a sequence quality of Q40 or better.
Based on a novel adaptive seed finding algorithm and the wave-based local aligner developed for
[daligner (2012)](https://github.com/thegenemyers/DALIGNER), the tool can for example compare
two 2Gbp bat genomes finding almost all regions over 100bp that are 70% or more similar
in about 4.5 minutes wall clock time on my MacPro with 8 cores (about 26 CPU minutes).
Moreover, it uses a trace point concept to record all the found alignments in a very space-efficient manner, e.g. just 53MB for over
635,000 local alignments (la's) in our running example.
These trace point encodings of the alignments can then be swiftly translated into .psl format, .paf format, or 1-code format on demand with programs provided here.

Using **FastGA** can be as simple as calling it with two FASTA files containing genome
assemblies where each entry is a scaffold with runs of N's separating and potentially giving the
estimated distance between the contigs thereof.  By default a PAF file encoding all the
local alignments found between the two genomes is streamed to the standard output.  In the
subdirectory ```EXAMPLE``` you will find a pair of sample input files, an output file, and a text file, ```sample_session``` capturing a session that serves to illustrate the use of FastGA.
Try it for yourself.

Under the surface, a number of intermediate steps take place.  First, the FASTA files
are converted to **genome databases** with extension **.gdb** that compress and index the ASCII DNA sequences into 2-bit compressed form.  This allows FastGA to randomly access contigs
and do so with four times less IO.  Second, a **genome index**
with extension **.gidx** is then built for each genome that is basically a truncated suffix array.
One of the things that makes FastGA fast is that it compares these two indices against
each other directly rather than looking up sequences of one genome in the index of the
other.
Thirdly, FastGA records all the alignments it finds in a file with extension **.las** that uses a very space efficient trace point encoding of each alignment.  Finally in linear time, this trace point representation is converted into the desired PAF output.
Note carefully, that one has the option to keep the results in the very disk efficient LAS file, and then convert it to any of PAF, PSL, ONEcode, or other desired format on demand.  The diagram immediately below summarizes and details the data flow just described.

![Fig. 1](System.Fig1.png)

While the entire set of blue shadowed processes can be fired off by simply calling FastGA, we
provide routines to perform each step under direct control (labeled in blue along dataflow arrows).
In addition we provide utilities labeled in brown what allow one to examine the intermediate
GDB, GIX, and LAS files.  An invocation of FastGA with the -k option or direct application of the sub-process routines create persistent GDB and GIX entities that can be reused saving time if
a given genome is to be compared repeatedly.  The GDB and GIX items are actually an ensemble of a proxy file and a number of hidden
files.  So we provide the utilities GIXmv, GIXcp, and GIXrm to
manipulate these as an ensemble.  Finally, we provide the utility GDBtoFA that perfectly inverts
the process of converting a FASTA into a GDB, providing the option of removing all your fasta
files, compressed or not, for the more efficient GDB representation.

There are three conventions for all the tools in this package designed for your convenience.
First, suffix extensions need not be given for arguments of a known set of types.  For example,
if an argument is a fasta with root name "foo" without extensions, then
our commands will look for ```foo.fa,``` ```foo.fna,``` ```foo.fasta,```
 ```foo.fa.gz,``` ```foo.fna.gz,``` ```foo.fasta.gz``` if you specify
```foo``` as the argument.  Second, optional arguments (those that begin with a -) can
be in any order and in any position relative to the non-optional primary arguments (which must
be given in the order specified).  We find this pretty convenient when for example you
have typed out an entire FastGA command but forgot that you wanted PSL output instead
of the default PAF output.
All you do is append -psl to what you've already typed and then hit return.  So for example,
```FastGA -v Asm1 -T16 Asm2 -psl``` is an acceptable command line.  Finally, if a -v option
is specified for a command then it always means "verbose mode", i.e. output to the standard
error a running discourse of the command's progress.

<a name="FastGA"></a>

## FastGA Reference

```
FastGA [-vk] [-T<int(8)>] [-P<dir(/tmp)] [<format(-paf)>]
          [-f<int(10)>] [-c<int(100)>] [-s<int(500)>] [-a<int(100)>] [-e<float(.7)>]
          <source1:path>[<precursor] [<source2:path>[<precursor>]]
          
    <format> = -paf[mx] | -psl | -one | -las 
        
    <precursor> = .gix | .gdb | <fa_extn>
    
    <fa_extn>   = (.fa|.fna|.fasta)[.gz]
```

Performing a FastGA comparison can be as simple as issuing the command ```FastGA A B``` where A and B are FASTA or gzip'd FASTA files.  By default 8 threads will be used but this can be changed with the -T
parameter.  By default the myriad temporary files produced by FastGA are located in /tmp but this directory
can be changed with the -P option.  All the alignments found by FastGA are streamed to the standard output
and by default will be in PAF format.  You can change this to PSL, ONEcode, or LAS formatted output with
the -psl, -one, and -las options, respectively.
Be careful though, the ONEcode and LAS streams are binary, and in these cases you are meant to capture 
the output in a file with the appropriate extension.
The -paf option can further be modulated with an 'x' or 'm', e.g. -pafx, which further requests that CIGAR
strings detailing the alignments be output (see [LAStoPAF](#LAStoPAF) below).

You can also call FastGA on a single source, e.g. ```FastGA A```, in which FastGA compares A against
itself, carefully avoiding self matches.  This is useful for detecting repetititve regions of the
genome (and their degree of repetitiveness), and for finding homologous regions between haplotypes in a
haplotype-phased but otherwise unsegregated genome assembly.

The one or two source arguments to FastGA can be either a FASTA file, a precomputed genome database
(GDB), or a precomputed genome index (GIX).  FastGA determines this by looking at the extension of
the argument if it is given explicitly, or if only the "root" name is given then it looks first for
a GIX with extenxion .gix, then a GDB with extension .gdb, and finally, for a fasta file with one of
the extensions .fa, .fna, .fasta, .fa.gz, .fna.gz, or .fasta.gz.  If a GIX is not present then FastGA
makes one, and if in turn a GDB is not present than one is made.  The objects so made are removed
upon completion of the execution of FastGA unless the -k option is set in which case they are kept.
Note carefully, any object already in the file system is not affected.  We recommend that one *replace*,
using [FAtoGDB](#FAtoGDB), every FASTA file with its GDB, as the GDB occupies less disk space and its originating FASTA can be reproduced exactly with [GDBtoFA](#GDBtoFA) on demand.
Similarly, if a group of genomes will be compared against each other, we recommend that one build a
GIX for each beforehand with [GIXmake](#GIXmake).  It takes about 30seconds / gigabase to build a GIX,
so building
them prospectively saves time as FastGA need not do so every time it is called on the same genome.
On the otherhand, GIXs are large, occupying 14GB for every gigabase of a genome, so we recommend that
one should build these as a prelude to a series of FastGA invocations and then remove them (but not
their GDBs) afterwords with [GIXrm](#GIXrm).

All the other options control the alignment discovery process.  Generally the defaults are fine and you
shouldn't bother touching these dials unless you are curious or confident.  For those willing to go further
FastGA uses the indices to find adaptive seed hits, where an
**adaptive seed** at a given position p of source1, is the longest string beginning at that position that is also somewhere in source2.  If the number of
occurences of this string in source2 is greater than the -f option, default value 10, then the adaptamer is deemed *repetitive* and is not considered.  Otherwise **adpatamer seed hits** occur at (p,q) for each q in  source2 where the adapamer at p also occurs.  The option -f, if specified, must be less than or equal to the value of -f specified if GIXmake is invoked separately to make the index in advance of calling FastGA.

FastGA then searches for runs or chains of adaptamer seed hits that (a) all lie within a diagonal band of width 128, (b) the spacing between every two consecutive seeds is less than -s(500), and
(c) the seeds in the chain cover at least -c(100) bases in both genomes.  For these **chain
hits**, FastGA then runs a wave-based local alignment routine that searches for a local alignment
of length at least -a(100)bp with a similarity of -e(70%) or better that contains at least one
of the seeds in the chain.  All such found alignments are recorded as a trace-point encoding in
lexicgraphical order of source1 contig #, and then the source2 contig #, and then the start
coordinate of the alignment in source1.  The options -s, -c -a, and -e can be used to modify the default
thresholds for chaining and alignment just described.

<a name="subprocess"></a>

## Sub-Process Routines

<a name="FAtoGDB"></a>

```
1. FAtoGDB [-v] <source:path>[<fa_extn>] [<target:path>[.gdb]]

       <fa_extn> = (.fa|.fna|.fasta)[.gz]
```

FAtoGDB takes as input a FASTA file with extension .fa, .fna, .fasta, .fa.gz, .fna.gz, or .fasta.gz
and produces a **genome database** or GDB with extension .gdb.
The name and location of the resulting GDB is determined as follows:

* If only a source file is given, then the GDB is built in the same directory and with the same core name as the FASTA file.

* If a target path is given and it is a directory, then the GDB is built in the given directory with the same core name as the FASTA file.

* If the target path is given and it is to a file name (that may not exist) then the directory and core name of the GDB are as for this target path.

A few examples: ```FAtoGDB A.fa``` and ```FAtoGDB PATH/A.fna .``` both produce A.gdb in the current directory, ```FAtoGDB PATH/A.fa.gz BLUE/AG.gdb``` produces AG.gdb in the directory BLUE (which must
exist).


<a name="GIXmake"></a>

```
2. GIXmake [-v] [-T<int(8)>] [-P<dir(/tmp)>] [-k<int(40)>] [-f<int(10)>]
            ( <source:path>[.gdb]  |  <source:path>[<fa_extn>] [<target:path>[.gix]] )
            
       <fa_extn> = (.fa|.fna|.fasta)[.gz]
```

Given a source FASTA or GDB file, GIXmake produces a **genome index** or GIX of the source with extension
.gix, creating the intermediate GDB if needed.
The name of the .gix is determined exactly as for FAtoGDB immediately above, save that
you are not allowed to create a GIX with a different core name and location than the GDB
it is derived from.  So note that in the summary command line syntax above you cannot specify
a target if the source is a GDB, you can only do so if one is starting from a FASTA file in which
both the GDB and GIX are created as per the target directive (if present).

The -T option can be used to specify the number of threads to use, where the default is 8.
The -P option similarly allows one to override the default /tmp, as the directory where the
(quite large and numerous) temporary files produced by GIXmake are held during its execution.
When running on an
HPC cluster node it is very important that this directory be on the disk local to the node
running the command.

The genome index basically consists of two parts: (1) a sorted table of the k-mers (k=40 by default) in the underlying genome that occur -f or fewer times (f=10 by default) along with the number of occurrences,
and (2) a list of all the positions in the genome that have a k-mer in the table, in the order in
which their k-mers occur in the table.  The .gix file is actually just a proxy for an ensemble
of -T hidden files with the extension .ktab.\<int\> that contain the k-mer table, and -T hidden files with
the extension .post.\<int\> that contain the position list.  Altogether these files occupy about 13-14GB
per gigabase of the genome and so a GIX is quite large.  Due to this structure we strongly recommend
that when you want to delete, copy, or move a GIX and its GDB, that rather than doing it piecemeal by
hand, you use the utilities [GIXrm](#GIXrm), [GIXcp](#GIXcp), [GIXmv](#GIXmv) that will handle not only
the proxy .gix file, but also the entire ensemble of hidden files as a single entity.

While you can reset the k-mer size with the -k option we
strongly suggest you use the default value of 40 or at least use a bigger value which will cost
you more compute time.  The option -f is designed to
remove k-mers from repetitive regions of the genome: only k-mers that occur -f or fewer times
are kept in the index.  The default value of -f is 10 and again we strongly suggest you use
this default.  Increasing it will improve sensitivity at the expense of more time and space,
decreasing it, the converse.  The effect is quadratic in -f so take care.

<a name="LAStoPAF"></a>

```
3. LAStoPAF [-mx] [-T<int(8)>] <alignments:path>[.las]
```

LAStoPAF converts an LAS file into a [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md) file, streaming the PAF to the standard output.
LAStoPAF uses 8 threads by default, but this can be changed with the -T option.

The command must have access to the one or two GDB's from which the LAS file was derived.
So the path, both relative and absolute, of these is recorded within the LAS file at the time
it is created.  So one should be careful not to move or rename these GDBs, the one exception being
if you move them so that when you call LAStoPAF they are at the same relative location from the
current directory as was true at the time of creation.
If you do have to rename or otherwise move the GDBs, then you can change the LAS file's internal
references to the new GDB locations with [LASreset](#LASreset).

In addition to the standard PAF fields, LAStoPAF outputs a ```dv:F:<fraction>``` SAM-tag that gives the divergences of the query from the target and a ```df:I:<diffs>``` SAM-tag that gives the number
of differences in an optimal alignment between the two intervals of the query and target.

The -m and -x options request LAStoPAF to produce a CIGAR string tag of the form ```cg:Z:<cigar-string>``` that explicitly encodes an optimal alignment.  With the -m option, aligned characters regardless of whether they are equal or not are
encoded with an 'M'.  With the -x option, aligned *equal* characters are encoded with an '=' and
aligned *unequal* characters with an 'X'.

*Beware*, the -m and -x options increase the time taken by LAStoPAF by a factor of 10 and the file size
by a factor of almost 100 !  The time taken can
be ameliorated somewhat by running LAStoPAF with more threads, controllable with the -T option.

<a name="LAStoPSL"></a>

```
4. LAStoPSL [-T<int(8)>] <alignments:path>[.las]
```

LAStoPSL converts an LAS file into a [PSL](https://www.ensembl.org/info/website/upload/psl.html) file,
streaming the PSL to the standard output.
LAStoPSL uses 8 threads by default, but this can be changed with the -T option.

The command must have access to the one or two GDB's from which the LAS file was derived.
So the path, both relative and absolute, of these is recorded within the LAS file at the time
it is created.  So one should be careful not to move or rename these GDBs, the one exception being
if you move them so that when you call LAStoPAF they are at the same relative location from the
current directory as was true at the time of creation.
If you do have to rename or otherwise move the GDBs, then you can change the LAS file's internal
references to the new GDB locations with [LASreset](#LASreset).

*Warning*, the PSL output is almost 15 times larger than the LAS file.

<a name="LAStoONE"></a>

```
5. LAStoONE [-T<int(8)>] <alignments:path>[.las]
```

Not yet.


<a name="viewing"></a>

## Viewing Utilities

<a name="GDBshow"></a>

```
1. GDBshow [-hU] [-w<int(80)>] <source:path[.gdb] [ <ranges:FILE> | <range> ... ]

       <range> =    <contig>[-<contig>]     |  <contig>_<int>-(<int>|#)
               | @[<scaffold>[-<scaffold>]] | @<scaffold>_<int>-(<int>|#)

         <contig>   = (<int>|#)[.(<int>|#)]
         
         <scaffold> =  <int>|<string>|#
```

GDBshow allows one to view a given set of scaffolds/contigs or portions thereof for the source GDB.
If the -h option is set then only the header lines are shown.  By default DNA sequence is lower-case,
80bp per row.  You can request upper-case with -U, and set the line width with -w.  If no arguments
besides the source DB are given, then the contigs of the entire GDB are output (in order).  If a file name follows, then GDB interprets each line of the file as a display directive.  Otherwise each argument after the source GDB is interpreted as a display directive where the syntax and meaning is as follows:

* A contig (index) is either (1) a single integer, say c, denoting the c'th contig in the genome, or (2) a pair of integers separated by a ., say s.c, denoting the c'th contig of the s'th scaffold in the genome.

* A scaffold (index) is either a single integer, say s, denoting the s'th scaffold in the genome, or
   it can be a string in which case it denotes every scaffold whose header prefix (excluding any leading
   white space) matches the string. 

* The special symbol # which may substitute for an integer denotes "the last", e.g. # adresses the last contig in the genome, #.1 addresses the 1st contig of the last scaffold, and 1_500-# selects the substring from 500 to the end of contig 1.

* A range either begins with an @ sign, in which case all indices are to scaffolds, or it does not
 in which case all indices are to contigs.  If the range is (1) a single index, then the given contig or scaffold is displayed, (2) a pair of indices separated by a hyphen, then the range of contigs or scaffolds are displayed (inclusively), or (3) a single index followed by an under-bar and then two integers separated by a hyphen, then the substring
of the contig or scaffold between the indices given by the two integers is displayed.

* If no directives are given then every contig is output in order, i.e. 1-#, and if a single
 @-sign is given then every scaffold is output in order, i.e. @1-#.

* A request to display a contig, prints the contig's scaffold's header with the relative contig
 number and interval of the scaffold appended followed by the sequence of the contig.  A request to display a scaffold, prints the scaffold's header followed by the sequence of the scaffold which is each of its' individual contigs separated by runs of N's of the same length as they occurred in the originating FASTA file.

While contigs and scaffolds are referenced by ordinal numbers starting at 1, sequence positions are
between characters and so start at 0 (the position immediately to the left of the 1st character).
The convenience of this 0-based numbering is that the length of an interval [b,e] is e-b.  Further
note that when the item for a substring is a contig then the interval is respect to the contig,
while when the item is a scaffold it is respect to the scaffold considered as a single string
with intervening N's between contigs.
 
<a name="GDBstat"></a>

```
2. GDBstat [-h[<ctg:int>,<scaf:int>]] [-hlog] <source:path[.gdb]
```
GDBstat gives you summary statistics about the genome in the source GDB.  It always outputs the number,
cumulative base pairs, and average size of scaffolds, contigs, and gaps.  It also outputs the max, min,
and N<x> sizes for scaffolds and contigs, for <x> a multiple of 10.

If the -h option is specified then GDBstat also displays histograms of contig and scaffold sizes.
If the complete option is -hlog then GDBstat displays a histogram with logarithmically scaled
buckets, i.e. 1, 2, 5, 10, 20, 50, 100 ...
If two integers follow the -h then the first sets the histogram bucket size for the contig histogram
and the second for the scaffold histogram.  If nothing follows the -h then GDBstat picks a round
numbered bucket size so that each histogram has about 20 buckets (defined constant NBINS in GDBstat.c if you'd like to change that).

<a name="GIXshow"></a>

```
3. GIXshow <source:path[.gix] [ <address>[-<address>] ]

       <address> = <int> | <dna:string>
```

GIXshow displays all or a range of k-mers in a genome index along with the positions for each k-mer.
If an argument does not follow the source then the entire GIX is output starting with the smallest k-mer.
Otherwise the extra argument denotes a range of k-mers either as integer positions, e.g. the i'th k-mer
(in alphabetical order), or if a dna string is given it specifies the first k-mer whose prefix matches
the string (or the last if it is the second argument of a range).

<a name="LASshow"></a>

```
4. LASshow [-arU] [-i<int(4).] [-w<int(100)>] [-b<int(10)>>
              <alignments:path>[.las] [ <range> [<range>] ]

       <range> =    <contig>[-<contig>]     |  <contig>_<int>-(<int>|#)
               | @[<scaffold>[-<scaffold>]] | @<scaffold>_<int>-(<int>|#)

         <contig>   = (<int>|#)[.(<int>|#)]
         
         <scaffold> =  <int>|<string>|#
```

LAshow produces a printed listing of a subset of the local alignments contained in the specified
LAS file, where one can optionally view the alignment in a BLAST like format.  If just the LAS is given
as an argument then every alignment is displayed.  If a single range is given in addition, then only
those alignments whose interval in the 1st genome intersects the range are displayed.
If a pair of ranges are given then those alignments where its interval in the 1st genome intersects
the 1st range and its interval in the 2nd genome intersects the 2nd range are displayed.
See the documentation for [GDBshow](#GDBshow) for a detailed explanation of the format and meaning
of ranges.

The command must have access to the one or two GDB's from which the LAS file was derived.
So the path, both relative and absolute, of these is recorded within the LAS file at the time
it is created.  So one should be careful not to move or rename these GDBs, the one exception being
if you move them so that when you call LASshow they are at the same relative location from the
current directory as was true at the time of creation.
If you do have to rename or otherwise move the GDBs, then you can change the LAS file's internal
references to the new GDB locations with [LASreset](#LASreset).

If the -a or -r option is set then an alignment of the local alignment is displayed.
The -a option puts exactly -w columns per segment of the display, whereas the -r option
puts exactly -w a-read symbols in each segment of the display.  The -r display mode is
useful when one wants to visually compare two alignments involving the same a-read.
If both the -a, and -r flags are set, then the -a alignment comes first followed by the
-r alignment.  The -i option sets the indent for the alignment displays,
if they are requested.  The -b option sets the number of symbols on
either side of the aligned segments in an alignment display, and -U specifies that
uppercase should be used for DNA sequence instead of the default lowercase.

When examining LAshow output it is important to keep in mind that the coordinates
describing an interval of a read are referring conceptually to positions between bases
starting at 0 for the position to the left of the first base.  That is, a coordinate c
refers to the position between the c-1'st and c'th base, and the interval [b,e] captures
the e-b bases from the b'th to the e-1'st, inclusive.  We give an example with
part of an alignment for which we will explain several additional features:

```
2.1 6.05c <0..3,110] x [162,039..158,933] ~ 6.18% (181,250,738 x 162,401,845 bps, 192 diffs, 32 trace pts)

             0 ..........ctaaccctaaccctaaccctaaccctaaccctaaccctaa
                         |||||||||||||||||||||||||||||||||**|||||
        162049 aaccctaaccctaaccctaaccctaaccctaaccctaacccta--cctaa   5.0%

            40 ccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaacc
               ||||||||||||||||||||||||||||||||||||||||||||||||||
        162001 ccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaacc   0.0%

            90 ctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccct
               |||||||||||||||*||||||||||||||||||||||||||||||||||
        161951 ctaaccctaacccta-ccctaaccctaaccctaaccctaaccctaaccct   2.0%

           140 aaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaa
               ||||||||||||||||||||||*|||||||||||||||||||||||||||
        161902 aaccctaaccctaaccctaacc-taaccctaaccctaaccctaaccctaa   2.0%

      . . . .
```

The display of an LA always begins with a line giving the A-scaffold & contig, then the B-scaffold & contig, then an indication of orientation (i.e. 'n' for same strand, and 'c' for the opposite
strand) followed by the A-interval and B-interval that are aligned in scaffold string coordinates
and then the % identify and then in parentheses
the lengths of the two scaffolds, the number of differences in the alignment, and the number of tracepoints used to encode the alignment between them.
In particular,
note carefully that when the B-item is in the complement orientation (c), then the
B-interval gives the higher coordinate first, the idea being that one will align from
the highest base down to the lowest base in the descending direction on B, complementing
the characters as you go.  Further note that in the alignment display the coordinates at
the start of each line follow this orientation convention and give the coordinate of the
"tick mark" just left of the first character in each line.  It is useful to know if an
interval reaches the end of a read, and to signal this we use an angle-bracket \<\> instead
of a square bracket [].

<a name="addons"></a>

## Additional Utilities

<a name="GDBtoFA"></a>

```
1. GDBtoFA [-vU] [-w<int(80)> <source:path>[.gdb] [ @ | <target:path>[<fa_sten>] ]

       <fa_extn> = (.fa|.fna|.fasta)[.gz]
```
GDBtoFA will produce exactly the FASTA file contents from which the source GDB was derived by
a call to FAtoGDB, i.e. it is an exact inverse operation.  When a GDB is built it records internally
where the source FASTA file is, it's name, and it's extension.  We call this the "origin" in what
follows.  Where the FASTA file is placed by GDBtoFA and what it is named is as follows:

* If there is no target, then the output is streamed to the standard output (uncompressed).

* If the target is the special symbol @, then the FASTA file is built at it's origin directory and given it's origin name and extension (implying compression if it ends with .gz).

* If the target is a directory, then the FASTA file is built at said directory, it's name and extension are as for the origin.

* If the target is a file (that may not exist), then the FASTA file is build at the directory and named as given by the target.  If the target has an extension then that extension is used, otherwise the extension of the origin is used.


<a name="GIXrm"></a>
<a name="GIXcp"></a>
<a name="GIXmv"></a>

```
2.a GIXrm [-vifg] <source:path>[.gix|.gdb] ...
2.b GIXmv [-vinf] <source:path>[.gix|.gdb] <target:path>[.gix|.gdb]
2.c GIXcp [-vinf] <source:path>[.gix|.gdb] <target:path>[.gix|.gdb]
```

A GDB consists of not only a proxy file with a .gdb extension but 3 hidden files as well with extensions
.hdr, .idx, and .bps.
Likewise a GIX consists of not only a proxy file with a .gix extension but 2T hidden files with
extensions of the form .ktab.\<int\> and .post.\<int\> where T is the number of threads used to
create the index.
As such it is cumbersome to remove, move, or copy a GDB or GIX directly with the UNIX OS as it requires
you to utter 2T+5 commands or possibly only one if wild cards are used, albeit this has the potential for
surprising matches affecting unexpected files.
So we provide the commands GIXrm, GIXmv, and GIXcp to remove, move, or copy GDBs and/or GIXs as a
single entity.

The routines operate exactly as for rm, mv, and cp including the meaning/effect of the flags -v, -i, -n,
and -f.  For the GIXmv and GIXcp commands both the GIX and GDB are moved or copied, where just the
GDB will be affected if there is no associated GIX.  On the otherhand, for GIXrm only the GIX will be
deleted unless the -g flag is explicitly set in which case the GDB will also be removed.  We chose
this requirement as a safeguard because if you have replaced your FASTAs with the GDBs as we recommend,
then deleting the GDB is tantamount to deleting your genome source! 

```
3. LASreset <alignments:path>[.las] <DB1:path>[.gdb] [<DB2:path>[.gdb]]
```

In the unfortunate event that the internal reference of an LAS file to its genome databases have become
"stale", LASreset allows you to reset these paths within the given file.
