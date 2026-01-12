# FastGA: A Fast Genome Aligner
  
<font size ="4">**_Authors:  Gene Myers & Chenxi Zhou_**<br>
**_First:   May 10, 2023_**<br>
**_Last:  December 30, 2025_**<br>

- [FastGA](#FastGA) Compare two genomes or a genome against itself and output a .1aln, .paf, or .psl file of all alignments found.

- [Sub-Process Routines](#subprocess)
  - [FAtoGDB](#FAtoGDB): Convert a FASTA or ONEcode sequence file into a genome database (GDB)
  - [GIXmake](#GIXmake): Build a genome index (GIX) for a given GDB
  - [ALNtoPAF](#ALNtoPAF): Stream PAF formatted alignments for a given .1aln file
  - [ALNtoPSL](#ALNtoPSL): Stream PSL formatted alignments for a given .1aln file

- [Viewing Utilities](#viewing)
  - [GDBshow](#GDBshow): Display select contigs or substrings thereof from a GDB
  - [GDBstat](#GDBstat): Display various statistics and histograms of the scaffolds & contigs in a GDB
  - [ANOshow](#ANOshow): Display annotation intervals of select contigs or subranges thereof from an ANO file
  - [ANOstat](#ANOstat): Display various statistics and histograms of about the intervals in an ANO file
  - [GIXshow](#GIXshow): Display range of a GIX
  - [ALNshow](#ALNshow): Display selected alignments in a .1aln file in a variety of forms
  - [ALNplot](#ALNplot): Display alignments in a .1aln or .paf file in a static collinear plot

- [Additional Utilities](#addons)
  - [GDBtoFA](#GDBtoFA): Converts a GDB back to the FASTA or ONEcode sequence file it was derived from
  - [BEDtoANO](#BEDtoANO): Convert a BED formatted file to a .1ano-file
  - [ANOtoBED](#ANOtoBED): Convert an ANO file to a BED file
  - [PAFtoALN](#PAFtoALN): Convert a PAF formatted file with X-CIGAR strings to a .1aln file
  - [PAFtoPSL](#PAFtoPSL): Convert a PAF formatted file with X-CIGAR strings to a .psl file
  - [GIXrm](#GIXrm): Remove GDBs and GIXs including their hidden parts
  - [GIXcp](#GIXcp): Copy GDBs and GIXs including their hidden parts as an ensemble
  - [GIXmv](#GIXmv): Move GDBs and GIXs including their hidden parts as an ensemble
  - [ALNchain](#ALNchain): Alignment filtering by construction of local chains
  - [ALNreset](#ALNreset): Reset a .1aln file's internal references to the GDB(s) it was computed from

- [C-Library for Accessing .1aln Files](#ONEaln)  
  - [Error Handling](#ehandler)
  - [.1aln File Reader](#areader)
  - [Genome Database (GDB)](#gdb)
  - [Alignment Records](#arecord)
  - [Error Messages](#emessage)


## Version 1.5 (December 30, 2025) ONEcode ANO Files

The addition of soft masking in V1.3 has lead to the development of a ONEcode version of a BED file,
called an ANO- or .1ano-file,
that records a collection of (oriented) intervals on a genome, along with a possible label and/or score for each annotation.
There are new routines ANOshow, ANOstat, ANOtoBED, and BEDtoANO for showing the contents of a .1ano file,
giving summary statistics about an .1ano file, and converting between BED files and .1ano files, respectively.

This change augments the interface to **GIXmake**, **FastGA**, **GDBtoFA**, and **GDBshow**, and has changed
the operation of
**FAtoGDB**.  Previously, if FAtoGDB detected an "implicit" mask in the source FASTA file indicated by masked regions being
lower case, and unmasked regions upper case, then this mask was recorded in the GDB.  Now it is recorded
in a separate .1ano file that has the same location and root name as the target GDB.  Further, GIXmake has been upgraded
so that you can specify multiple .1ano files on the command line, and the GIX will be masked with the
union of these.  Note carefully that the only way to change the mask encoded in a GIX is to rebuild the
GIX.  By default FastGA now does not use the soft mask in the GIX's of the genomes, but the -M flag
instructs it to so so.  You can also now follow each FastGA genome argument
with a list of masks to apply (syntactically a primary argument beginning with #) in which case the
GIX's will be rebuilt with the specified mask(s) and FastGA will soft mask accordingly.  Lastly, GDBtoFA and GDBshow now take an optional
\#-sign mask argument that if present masks the result accordingly.

## Version 1.4 (November 1, 2025) ONEaln C-Library

A C-library of routines designed to make it easy to read and access the contents of a .1aln file
has been added.  The interface is described [here](#ONEaln). This library of routines is in ONEaln.c
with the interface declared in ONEaln.h.
**Caution:** several of the modules used by FastGA must also be compiled in, namely, GDB.[ch],
ONElib.[ch], alncode.[ch], align.[ch], and gene_core.[ch].  See the make command for ONEalnTEST in
the Makefile.

## Version 1.3 (July 23, 2025) Soft Masking and Log Files
Soft masking is now supported and taken advantage of by the FastGA suite.  Soft masking is assumed
to be specified in the input FastA files by denoting masked sequence in lower-case and unmasked sequence in upper-case.  Such masks are recorded in our GDB's and in a suitable form within our
GIX indices.  The later required a slight modification to the GIX data structure.  Old GIX's are
still recognized and supported, but if you want masking you must rebuild any GDB's and GIX's that
were produced previously.

Additionally,

* GIXmake has been substantially improved to use much less memory and make better use of threads.

* The IO performance of ALNtoPAF and ALNtoPSL has also been substantially improved.

* A -L log file option has been added to support HPC cluster usage.

We are seeking similar improvements in FastGA proper, better handling of satellitic repeats, and
higher sensitivity for distant genomes without resorting to using LastZ as a subroutine.

## Overview

**FastGA** searches for all local DNA alignments between two high quality genomes.
The core assumption is that the genomes are nearly complete involving at most several
thousand contigs with a sequence quality of Q40 or better.
Based on a novel adaptive seed finding algorithm and the first wave-based local aligner developed for
[daligner (2012)](https://github.com/thegenemyers/DALIGNER), the tool can for example compare
two 2Gbp bat genomes finding almost all regions over 100bp that are 70% or more similar
in about 5.0 minutes wall clock time on my MacPro with 8 cores (about 28 CPU minutes).
Moreover, it uses a trace point concept to record all the found alignments in a compressed
and indexable [ONEcode](https://github.com/thegenemyers/ONEcode) file in a very space-efficient manner, e.g. just 44.5MB for over
635,000 local alignments in our running example.
These trace point encodings of the alignments can then be swiftly translated into .psl or .paf format on demand with programs provided here.

Using **FastGA** can be as simple as calling it with two FASTA files containing genome
assemblies where each entry is a scaffold with runs of N's separating and potentially giving the
estimated distance between the contigs thereof.  By default a PAF file encoding all the
local alignments found between the two genomes is streamed to the standard output.  In the
subdirectory ```EXAMPLE``` you will find a pair of sample input files, an output file, and a text file, ```sample_session``` capturing a session that serves to illustrate the use of FastGA.
Try it for yourself.

Under the surface, a number of intermediate steps take place.  First, the FASTA files
are converted to **genome databases** with extension **.1gdb** that are a [ONEcode](https://github.com/thegenemyers/ONEcode) binary file and associated hidden file containing the
ASCII DNA sequences in 2-bit compressed form.  This allows FastGA to randomly access contigs
and do so with four times less IO and no text parsing.  Second, a **genome index**
with extension **.gidx** is then built for each genome that is basically a truncated suffix array.
One of the things that makes FastGA fast is that it compares these two indices against
each other directly rather than looking up sequences of one genome in the index of the
other.
Third, FastGA records all the alignments it finds in a [ONEcode](https://github.com/thegenemyers/ONEcode) binary file we refer to here as a ALN-formated file with extension **.1aln** that uses a very space efficient trace point encoding of each alignment.  Finally in linear time, this trace point representation is converted into the desired PAF output.
Note carefully, that one has the option to keep the results in the very disk efficient ALN format, and then convert it to any of PAF, PSL, or other desired alignment format on demand.  The diagram immediately below summarizes and details the data flow just described.

![Fig. 1](System.Fig1.png)

While the entire set of blue shadowed processes can be fired off by simply calling FastGA, we
provide routines to perform each step under direct control (labeled in blue along dataflow arrows).
In addition we provide utilities labeled in brown that allow one to examine the intermediate
GDB, GIX, and ALN files.  An invocation of FastGA with the -k option or direct application of the sub-process routines, create persistent GDB and GIX entities that can be reused saving time if
a given genome is to be compared repeatedly.  The GDB and GIX items are actually an ensemble, consisting of a proxy file and a number of hidden
files.  So we provide the utilities GIXmv, GIXcp, and GIXrm to
manipulate these as an ensemble.  Finally, we provide the utility GDBtoFA that inverts
the process of converting a FASTA file into a GDB, providing the option of removing
all your fasta files, compressed or not, for the space efficient GDB representation.

FastGA features the use of the [ONEcode](https://github.com/thegenemyers/ONEcode) data encoding
framework with both its' GDB and ALN files that encode all the alignments found.  As such FastGA also
supports as input ONEcode sequence files that encode a genome, in addition to the usual
Fasta format.  So both FAtoGDB and GDBtoFA (despite their names) also recognize and support
ONEcode SEQ files as well as FASTA.

There are three conventions for all the tools in this package designed for your convenience.
First, suffix extensions need not be given for arguments of a known set of types.  For example,
if an argument is a fasta with root name "foo" without extensions, then
our commands will look for ```foo.fa,``` ```foo.fna,``` ```foo.fasta,```
 ```foo.fa.gz,``` ```foo.fna.gz,``` ```foo.fasta.gz``` if you specify
```foo``` as the argument.  Second, optional arguments (those that begin with a -) can
be in any order and in any position relative to the non-optional primary arguments (which must
be given in the order specified).  We find this convenient when for example you
have typed out an entire FastGA command but forgot that you wanted PSL output instead
of the default PAF output.
All you do is append -psl to what you've already typed and then hit return.  So for example,
```FastGA -v Asm1 -T16 Asm2 -psl``` is an acceptable command line.  Finally, if a -v option
is specified for a command then it always means "verbose mode", i.e. output to standard
error a running discourse of the command's progress, and if a -L option is available and
specified with a file name, then a log file of the command's performance is appended to said
file.

<a name="FastGA"></a>

## FastGA Reference

```
FastGA [-vkMS] [-L:<log:path>] [-T<int(8)>] [-P<dir($TMPDIR)>] [<format(-paf)>]
               [-f<int(10)>] [-c<int(85)> [-s<int(1000)>] [-l<int(100)>] [-i<float(.7)]
                 <source1:path>[<precursor>] (#[<mask>[.1ano]])*
               [ <source2:path>[<precursor>] (#[<mask>[.1ano]])* ] 

         <format> = -paf[mxsS]* | -psl | -1:<align:path>[.1aln]

         <precursor> = .gix | .1gdb | <fa_extn> | <1_extn>

             <fa_extn> = (.fa|.fna|.fasta)[.gz]
             <1_extn>  = any valid 1-code sequence file type
```

Performing a FastGA comparison can be as simple as issuing the command ```FastGA A B``` where A and B are FASTA, gzip'd FASTA, or ONEcode sequence files.  By default 8 threads will be used but this can be changed with the -T
parameter.  By default the myriad temporary files produced by FastGA are located in the directory given
by the environment varialble ```TMPDIR``` (or ```.``` if it is undefined) but this directory
can be changed with the -P option.  In the default case, all the alignments found by FastGA are streamed to the standard output in PAF format.  You can change this to PSL, or ONEcode ALN formatted output by specifying
the -psl and -1, options, respectively.
Note carefully, that the ONEcode -1 option produces binary output and so requires you also specify a file where the binaray will be stored.
The -paf option can further be followed by any combination of 'x', 'm', 's' or 'S', save that
'x' and 'm' do not both occur, and 's' and 'S' do not both occur.  These modulators specifiy that
CIGAR strings (x or m) or CS strings (s or S) be additionally output for each alignment, e.g.
-pafxS would request a CIGAR string with X and = in aligned regions, and a long form CS string
be output.
See the argument description for [ALNtoPAF](#ALNtoPAF) below for more details.

You can also call FastGA on a single source, e.g. ```FastGA A```, in which case FastGA compares A against
itself, carefully avoiding self matches.  This is useful for detecting repetititve regions of a
genome (and their degree of repetitiveness), and for finding homologous regions between haplotypes in an unphased genome assembly, or one that is phased but not split into separate haplotype files.

FastGA uses the adaptamer seed idea of Martin Frith which means that ```FastGA A B``` does not find
the same set of alignments as ```FastGA B A``` as the adaptamers of A and B are not the same.
Specifying the 'S' option (S for symmetric) requests FastGA to use the adaptamers of both genomes
taking a bit more time but producing roughly the same result independent of the order of A and B
on the command line.  The difference is typically that more repetitive alignments in B are found.
So we advise not to use the option when synteny is the goal, and only use it when understanding the
repeat structure of both genomes is of interest.

The one or two source arguments to FastGA can be either a FASTA file, a ONEcode sequence file (e.g. .1seq), a precomputed genome database
(GDB), or a precomputed genome index (GIX).  FastGA determines this by looking at the extension of
the argument if it is given explicitly, or if only the "root" name is given then it looks first for
a GIX with extenxion .gix, then a GDB with extension .1gdb, then a fasta file with one of
the extensions .fa, .fna, .fasta, .fa.gz, .fna.gz, or .fasta.gz, and if all else fails then FastGA
tries to open the file as a ONEcode (sequence) file.  If a GIX is not present then FastGA
makes one, and if in turn a GDB is not present than one is made.  The objects so made are removed
upon completion of the execution of FastGA unless the -k option is set in which case they are kept.
Note carefully, any object already in the file system is not affected.  We recommend that one *replace*,
using [FAtoGDB](#FAtoGDB), every FASTA file with its GDB, as the GDB occupies less disk space and its originating FASTA can be reproduced exactly with [GDBtoFA](#GDBtoFA) on demand.
Similarly, if a group of genomes will be compared against each other, we recommend that one build a
GIX for each beforehand with [GIXmake](#GIXmake).  It takes about 15-30 seconds per gigabase to build a GIX,
so building
them prospectively saves time as FastGA need not do so every time it is called on the same genome.
On the other hand, GIXs are large, occupying 14GB for every gigabase of a genome, so we recommend that
one should build these as a prelude to a series of FastGA invocations and then remove them (but not
their GDBs) afterwards with [GIXrm](#GIXrm).

Each source argument can be followed by a list of masks where each is a #-sign immediately followed by the
name of a .1ano file specifying a collection of intervals, e.g. ```FastGA A #mask1.1ano #mask2 B #``` will
soft mask ```A``` with the union of the intervals in ```mask1.1ano``` and ```mask2.1ano```, and ```B``` will be
masked with the "implicit" mask ```B.1ano``` that was specified in its' FastA file using upper and lower-case 
to denote unmasked versus masked sequence, respectively, or in the M-lines of its' ONEcode file, depending on the file type of its source.  Specifying # masks arguments requires that the GIX for the argument in question be recomputed
even if the GIX already exists at the time of the call.  If no explicit mask arguments are given then by default FastGA ignores the soft mask encoded in the GIX's, but if such mask arguments occur or the -M
flag is set, then FastGA uses the soft masks in the GIX's.

All the other options control the alignment discovery process.  Generally the defaults are fine and you
shouldn't bother touching these dials unless you are curious or confident.  For those willing to go further
FastGA uses the indices to find adaptive seed hits, where an
**adaptive seed** at a given position p of source1, is the longest string beginning at that position that is also somewhere in source2.  If the number of
occurences of this string in source2 is greater than the -f option, default value 10, then the
adaptamer is deemed *repetitive* and is not considered.  Otherwise **adpatamer seed hits** occur
at (p,q) for each q in  source2 where the adaptamer at p also occurs.

FastGA then searches for runs or chains of adaptamer seed hits that (a) all lie within a diagonal band of width 128, (b) the spacing between every pair of consecutive seeds is less than -s(1000), and
(c) the seeds in the chain cover at least -c(85) bases in both genomes.  For these **chain
hits**, FastGA then runs a wave-based local alignment routine that searches for a local alignment
of length at least -l(100)bp with a similarity of -i(70%) or better that contains at least one
of the seeds in the chain.  All such found alignments are recorded as a trace-point encoding in
lexicographical order of source1 contig #, and then the source2 contig #, and then the start
coordinate of the alignment in source1.  The options -s, -c -l, and -i can be used to modify the default
thresholds for chaining and alignment just described.

<a name="subprocess"></a>

## Sub-Process Routines

<a name="FAtoGDB"></a>

```
1. FAtoGDB [-v] [-L:<log:path>] [-n<int>] <source:path>[<extn>] [<target:path>[.1gdb]]

       <extn> = (.fa|.fna|.fasta)[.gz] or any valid 1-code sequence filel type
```

FAtoGDB takes as input a FASTA file with extension .fa, .fna, .fasta, .fa.gz, .fna.gz, or .fasta.gz,
or a ONEcode sequence file (e.g. with extension .1seq)
and produces a **genome database** or GDB with extension .1gdb.
The name and location of the resulting GDB is determined as follows:

* If only a source file is given, then the GDB is built in the same directory and with the same core name as the FASTA file.

* If a target path is given and it is a directory, then the GDB is built in the given directory with the same core name as the FASTA file.

* If the target path is given and it is to a file name (that may not exist) then the directory and core name of the GDB are as for this target path.

A few examples: ```FAtoGDB A.fa``` and ```FAtoGDB PATH/A.fna .``` both produce A.gdb in the current directory, ```FAtoGDB PATH/A.fa.gz BLUE/AG.gdb``` produces AG.gdb in the directory BLUE (which must
exist).

The GDB actually consists of two files.  The first, *visible* file, is a ONEcode binary file with extension
.1gdb that contains all the information about an assembly except for the base-pair sequences which
are kept in a separate hidden file in 2-bit compressed format.  If the visible file has name say,
```foo.1gdb``` then this hidden file has name ```.foo.bps```.  We split the GDB this way as many application do not actually need the sequence, but simply need the sizes of contigs, gaps, & scaffolds and their names which are kept in the "light-weight" .1gdb portion.

Normally, n's in the FASTA file are interpreted as gaps between contigs of a scaffold.  But some projects
produce reconstructions that use n's as an undetermined base.  Typically these occur in isolation or pairs
but not in long, say 200bp, runs of n's.  In this case, use the -n option to specify a threshold above which
a run of n's is considered a gap, and below which the n's are considered unkown bases (and therefore treated
arbitrarily as a's).

If the input is a FASTA and the convention that upper-case bases are unmasked and lower-case bases are masked,
then FAtoGDB will also produce a .1ano file of this "implicit" mask that has the same root name as the
resulting GDB and is in the same directory.  This implicit mask can be referred to with the option/flag -#,
*not followed by a name*, in the tools, e.g. GIXmake, that takes .1ano's as masking arguments.  An implicit
mask is also produced if the source is a .1seq file and it contains M-lines (mask lines).


<a name="GIXmake"></a>

```
2. GIXmake [-v] [-L:<log:path>] [-T<int(8)>] [-P<dir($TMPDIR>] [-k<int(40)>]
           <source:path>[.1gdb|<extn>]  (#[<mask>[.1ano]])*
or
   GIXmake [-v] [-L:<log:path>] [-T<int(8)>] [-P<dir($TMPDIR>] [-k<int(40)>]
           <source:path>[<extn>] <target:path>[.gix] (#[<mask>[.1ano]])*
            
       <extn> = (.fa|.fna|.fasta)[.gz] or any valid 1-code sequence file type
```

Given a source FASTA, ONEcode, or GDB file, GIXmake produces a **genome index** or GIX of the source with extension
.gix, creating the intermediate GDB if needed.
The name of the .gix is determined exactly as for FAtoGDB immediately above, save that
you are not allowed to create a GIX with a different core name and location than the GDB
it is derived from.  So note that in the summary command line syntax above you cannot specify
a target if the source is a GDB, you can only do so if one is starting from a FASTA file in which
both the GDB and GIX are created as per the target directive (if present).

The -T option can be used to specify the number of threads to use, where the default is 8.
The -P option similarly allows one to override the default directory given by the environmnt variable
```TMPDIR``` (or ```.``` if it is undefined), as the directory where the
(quite large and numerous) temporary files produced by GIXmake are held during its execution.
When running on an
HPC cluster node it is very important that this directory be on the disk local to the node
running the command.

A GIX can encode a softmasking of its underlying genome, and if this is desired then one or more arguments
following the source that begin with a #-sign followed immediately be the name of a .1ano file
can be supplied.  The masking is the union of the masks so supplied on the command line.  If # alone is
given then this refers to the "implicit" masking that was encoded in the case of bases in its' FASTA
file or the M-lines of its' 1SEQ file, depending on the type of its source.

The genome index basically consists of a sorted table of all the k-mers (k=40 by default) and their complements that begin with a (12,8) syncmer in the underlying genome where each k-mer instance is paired with the
position it came from, and the longest prefix thereof that should/could be masked.
The .gix file is actually just a proxy for an ensemble
of -T hidden files with the extension .ktab.\<int\> that contain parts
of the k-mer table.  Altogether these files occupy about 13-14GB
per gigabase of the genome and so a GIX is quite large.  Due to this structure we strongly recommend
that when you want to delete, copy, or move a GIX and its GDB, that rather than doing it piecemeal by
hand, you use the utilities [GIXrm](#GIXrm), [GIXcp](#GIXcp), [GIXmv](#GIXmv) that will handle not only
the proxy .gix file, but also the entire ensemble of hidden files as a single entity.

While you can reset the k-mer size with the -k option we
strongly suggest you use the default value of 40 or at least use a bigger value which will cost you more compute time.

<a name="ALNtoPAF"></a>

```
3. ALNtoPAF [-mxsS] [-T<int(8)>] <alignment:path>[.1aln]
```

ALNtoPAF converts a ALN file into a [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md) file, streaming the PAF to the standard output.
ALNtoPAF uses 8 threads by default, but this can be changed with the -T option.

The command must have access to the one or two GDB's from which the ALN file was derived.
So the path, both relative and absolute, of these is recorded within the ALN file at the time
it is created.  So one should be careful not to move or rename these GDBs, the one exception being
if you move them so that when you call ALNtoPAF they are at the same relative location from the
current directory as was true at the time of creation.
If you do have to rename or otherwise move the GDBs, then you can change the ALN file's internal
references to the new GDB locations with [ALNreset](#ALNreset).

In addition to the standard PAF fields, ALNtoPAF outputs a ```dv:F:<fraction>``` SAM-tag that gives the divergences of the query from the target and a ```df:I:<diffs>``` SAM-tag that gives the number
of differences in an optimal alignment between the two intervals of the query and target.

The -m and -x options request ALNtoPAF to produce a CIGAR string tag of the form ```cg:Z:<cigar-string>``` that explicitly encodes an optimal alignment.  With the -m option, aligned characters regardless of whether they are equal or not are
encoded with an 'M'.  With the -x option, aligned *equal* characters are encoded with an '=' and
aligned *unequal* characters with an 'X'.  Only one of 'x' or 'm can be specified.

The -s and -S options request ALNtoPAF to produce a CS string tag of the form ```cs:Z:<cs-string>``` that encodes the differences between the sequences in the short form (-s) or the entire query and reference sequences in the long form (-S).  Only one of 's' or 'S can be specified.

*Beware*, the -m, -x, -s and -S options increase the time taken by ALNtoPAF by a factor of 10 and the file size
by a factor of almost 100 !  The time taken can
be ameliorated somewhat by running ALNtoPAF with more threads, controllable with the -T option.

<a name="ALNtoPSL"></a>

```
4. ALNtoPSL [-T<int(8)>] <alignments:path>[.1aln]
```

ALNtoPSL converts a ALN file into a [PSL](https://www.ensembl.org/info/website/upload/psl.html) file,
streaming the PSL to the standard output.
ALNtoPSL uses 8 threads by default, but this can be changed with the -T option.

The command must have access to the one or two GDB's from which the ALN file was derived.
So the path, both relative and absolute, of these is recorded within the ALN file at the time
it is created.  So one should be careful not to move or rename these GDBs, the one exception being
if you move them so that when you call ALNtoPAF they are at the same relative location from the
current directory as was true at the time of creation.
If you do have to rename or otherwise move the GDBs, then you can change the ALN file's internal
references to the new GDB locations with [ALNreset](#ALNreset).

*Warning*, the PSL output is almost 15 times larger than the ALN file.



<a name="viewing"></a>

## Viewing Utilities

<a name="GDBshow"></a>

```
1. GDBshow [-hU] [-w<int(80)>] [#[<mask>[.1ano]]] <source:path[.1gdb] [ <selection> | <FILE> ]

        <selection> = <range>[+-] [ , <range>[+-] ]*

           <range> = <object/position> [ - <object/position> ]  | @ | .

              <object/position> = @ <scaffold> [ . <contig>] [ : <position> ]
                                |                . <contig>  [ : <position> ]
                                |                                <position>
 
                 <scaffold> = # | <int> | <identifier>
                 <contig>   = # | <int>         
                 <position> = # | <int> [ . <int> ] [kMG]
```

GDBshow allows one to view a given set of scaffolds/contigs or portions thereof for the source GDB.
If the <nobr>-h</nobr> option is set then only the header lines are shown.
By default DNA sequence is lower-case, 80bp per row.  You can request upper-case with -U, and set the line width with -w.  If a #-sign mask argument is present, then the -U flag is ignored, and the output is
upper case, except for masked regions which are in lower case.

If no other arguments
besides the source GDB are given, then all the contigs of the GDB are output (in order).  If a file name follows, then GDB interprets each line of the file as a \<range> selection.  Otherwise the argument is
directly interpreted as a directive as to which scaffolds or contigs or parts thereof to display.
We explain the syntax of a \<selection> above in a bottom up fashion coupled with examples:

* @s denotes the s'th scaffold in the genome, @# the last scaffold in the genome, and @id the scaffold
whose fasta header is id.  In this latter case the header is the part of the faster header line excluding
the '>' in the first column and any white space directly following it, and excluding any tailing white space.

* .c denotes the c'th scaffold in the genome, and .# the last contig. @s.c denotes the c'th contig of the s'th scaffold, so by extension, @s.# denotes the last contig of the s'th scaffold, @#.1 denotes the 1'st contig of the last scaffold, and so on.

* An integer p denotes an absolute position within the genome viewing it as the concatenation of its scaffolds -- not very useful generally.  But .c:p denotes the position p in the c'th contig, @s:p the
position p in the s'th scaffold, .c:# the last position in the c'th contig, @s.c:p the position p in
the c'th contig of the s'th scaffold, and so on.

* A position is considered to be a location *between* base pairs, so the first position is 0, and in general position i is the spot between the i'th and i+1'st base pairs.  In this way the number of symbols
 between position p and q > p is q-p.  Also, the last position is the length of the sequence.
 
* We support a more elaborate syntax for positions where the suffix symbols k, M, and G denote kilobases,
megabases, and gigabases and a decimal fraction is permited.  So for example, 10.1k is short for
position 10,100, and 2.0324M is position 2,032,400.

* A range can be a reference to a contig or scaffold in which case the entire contig or scaffold is
 being selected.  It can also be a pair of contigs or scaffolds separated by a hyphen (-) in which case all the objects from the first to the last are being selected.  For example, @3-@5 would cause GDBshow to display the 3'rd through 5'th scaffolds inclusive and .3-.5 would display the 3'rd through 5'th contigs
 of the genome.  A range can also be a position or a pair of positions separated by a hyphen.  In the
 first case the range is empty which is not useful for an app like GDBshow but as you'll see later it is usefule for
 an app like ALNshow that displays every alignment containing the position.  When a pair of positions is given then all the sequence between the pairs is selected.  For example, @1.3:10k-@1.5:20k selects all but
 the first 10kbp of the 3'rd contig of the first scaffold, all of it's 4'th contig, and the first 20k of it's 5'th contig.

* A range can be @ in which case every scaffold is selected, and it can be . in which case every contig is selected.  A genome can be considered to be a list of its' contig sequences, or a list of its' scaffold sequences each of which is its' contigs concatenated together with intervening runs of N's denoting the gap (and possibly length) between them.  If a range involves an @-sign, that is, a scaffold index, then the selection is over scaffold sequences, otherwise it is over contig sequences.

* The second component of a range when present, is assumed to have the same scaffold or contig prefix
if it is not given.  For example, the expression @1.3:10k-.5:20k is the same as @1.3:10k-@1.5:20k and
@2:100-200 is the same as @2:100-@2:200.

* A range can be followed by an optional sign, + or -.  This indicates the selection of either the 5' or 3' direction of the selected sequence interval and affects the direction of display.  For example, .3-
asks GDBshow to display reverse complement of the 3'rd contig of the genome.

* Lastly, on the command line, one can give a list of ranges separated by commas.  Generally this represents the union of the range intervals, but for GDBshow indicates that it should output each
range in the order given.  When a file is given, then each line of the file is assumed to be a
range so that the file also specifies a list of ranges.
 
<a name="GDBstat"></a>

```
2. GDBstat [-h[<ctg:int>,<scaf:int>]] [-hlog] <source:path[.1gdb]
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

<a name="ANOshow"></a>

```
3. ANOshow <source:path[.1ano] [ <selection> | <FILE> ]
```

In analogy with GDBshow documented above, ANOshow gives you a listing of the intervals and labels (if present)
in a given region or regions of the source genome for the ANO file.  So the ```<selection>``` argument if present has the same interpretation as for GDBshow, save that now it is the intervals in those regions that are displayed.

<a name="ANOstat"></a>

```
4. ANOstat [-h[<itvls:int>,<covrd:int>]] [-hlog] <source:path[.1ano]
```

ANOstat gives you summary statistics about the anotation intervals in an ANO file in analogy with GDBstat
documented above.  Options are the same, save that the -h integer pair set bucket size for a histogram
of interval lengths and covered segment lengths, respectively.

<a name="GIXshow"></a>

```
5. GIXshow <source:path[.gix] [ <address>[-<address>] ]

       <address> = <int> | <dna:string>
```

GIXshow displays all or a range of k-mers in a genome index along with the positions for each k-mer.
If an argument does not follow the source then the entire GIX is output starting with the smallest k-mer.
Otherwise the extra argument denotes a range of k-mers either as integer positions, e.g. the i'th k-mer
(in alphabetical order), or if a dna string is given it specifies the first k-mer whose prefix matches
the string (or the last if it is the second argument of a range).

<a name="ALNshow"></a>

```
6. ALNshow [-arU] [-i<int(4).] [-w<int(100)>] [-b<int(10)>>
              <alignments:path>[.1aln] [ <selection>|<FILE> [<selection>|<FILE>] ]

        <selection> = <range>[+-] [ , <range>[+-] ]*

           <range> = <object/position> [ - <object/position> ]  | @ | .

              <object/position> = @ <scaffold> [ . <contig>] [ : <position> ]
                                |                . <contig>  [ : <position> ]
                                |                                <position>
 
                 <scaffold> = # | <int> | <identifier>
                 <contig>   = # | <int>         
                 <position> = # | <int> [ . <int> ] [kMG]
```

ALNshow produces a printed listing of a subset of the local alignments contained in the specified
ALN file, where one can optionally view the alignments in a BLAST like format.  If just the ALN file is given
as an argument then every alignment is displayed.  If a single selection is given in addition, then only
those alignments whose interval in the 1st genome intersects the selection are displayed.
If a pair of selections are given then those alignments where its interval in the 1st genome intersects
the 1st selection and its interval in the 2nd genome intersects the 2nd selection are displayed.
If a range in the first selection is signed than that determines if the alignment is displayed with
the substring of the 1st genome in the normal orientation (+, also no sign) or the complement
orientation (-).  If a range in the second selection, if present, is signed than only alignments
in which the relative orientation of the first and second substrings agrees with the sign polarity
of the selection pair are displayed.
See the documentation for [GDBshow](#GDBshow) for a detailed explanation of the format and meaning
of selections.

The command must have access to the one or two source files from which the ALN file was derived.
This can be either a Fasta file, a ONEcode SEQ file, or a GDB depending on how the ALN file was created.
For example, if FastGA is run on two FASTA files without the -k option then the recorded sources
will be the FASTA files.  But with the -k option, then the GDB's (that are kept) are recorded.
The path, both relative and absolute, of these sources is recorded within the ALN file at the time
it is created.  So one should be careful not to move, rename, or remove the sources, the one exception being
if you move them so that when you call ALNshow they are at the same relative location from the
current directory as was true at the time of creation.
If you do have to rename or otherwise move the source files, then you can change the ALN file's internal references to their new locations with [ALNreset](#ALNreset).

If the -a or -r option is set then an alignment of the local alignment is displayed.
The -a option puts exactly -w columns per segment of the display, whereas the -r option
puts exactly -w a-read symbols in each segment of the display.  The -r display mode is
useful when one wants to visually compare two alignments involving the same a-read.
If both the -a, and -r flags are set, then the -a alignment comes first followed by the
-r alignment.  If the -n option is set while either of -a or -r is set, then the scaffold names
are used instead of scaffold numbers in the disply.
The -i option sets the indent for the alignment displays,
if they are requested.  The -b option sets the number of symbols on
either side of the aligned segments in an alignment display, and -U specifies that
uppercase should be used for DNA sequence instead of the default lowercase.

When examining ALNshow output it is important to keep in mind that the coordinates
describing an interval of a read are referring conceptually to positions between bases
starting at 0 for the position to the left of the first base.  That is, a coordinate c
refers to the position between the c'th and c+1'st base, and the interval [b,e] captures
the e-b bases from the b+1'st to the e'th, inclusive.  We give an example with
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

The display of a local alignment always begins with a line giving the A-scaffold & contig, then the B-scaffold & contig, then an indication of orientation (i.e. 'n' for same strand, and 'c' for the opposite
strand) followed by the A-interval and B-interval that are aligned in scaffold string coordinates
and then the % identify.  Then in parentheses follows
the lengths of the two scaffolds, the number of differences in the alignment, and the number of tracepoints used to encode the alignment between them.
In particular,
note carefully that when the B-item is in the complement orientation (c), then the
B-interval gives the higher coordinate first, the idea being that one will align from
the highest base down to the lowest base in the descending direction on B, complementing
the characters as you go.  Further note that in the alignment display the coordinates at
the start of each line follow this orientation convention and give the coordinate of the
"tick mark" just left of the first character in each line.  It is useful to know if an
interval reaches the beginning or end of a read, and to signal this we use angle-brackets \<\> instead
of square brackets [].

<a name="ALNplot"></a>

```
7. ALNplot [-vSL] [-T<int(4)>] [-p[:<output:path>[.pdf]]]
               [-a<int(100)>] [-e<float(0.7)>] [-n<int(100000)>]
               [-H<int(600)>] [-W<int>] [-f<int>] [-t<float>]
               <alignment:path>[.1aln|.paf[.gz]]> [<selection>|<FILE> [<selection>|<FILE>]]

        <selection> = <range>[+-] [ , <range>[+-] ]*

           <range> = <object/position> [ - <object/position> ]  | @ | .
            
              <object/position> = @ <scaffold> [ . <contig>] [ : <position> ]
                                |                . <contig>  [ : <position> ]
                                |                                <position>
 
                 <scaffold> = # | <int> | <identifier>
                 <contig>   = # | <int>         
                 <position> = # | <int> [ . <int> ] [kMG]
```

ALNplot produces a static collinear plot of the local alignments contained in the specified ALN file or PAF 
file in a EPS or PDF format file.  If just the ALN/PAF is given as an argument then every alignment is considered 
for plotting.  If a single selection is given in addition, then only those alignments whose interval in the 1st 
genome intersects the selection are considered. If a pair of selections are given then those alignments where 
its interval in the 1st genome intersects the 1st selection and its interval in the 2nd genome intersects the 2nd 
selection are considered. A selection can be comma-delimited to include multiple, interspersed ranges. These ranges
will be placed in order along the axis for plotting. The selection can also be a FILE, with each line representing 
a range. This is equivalent to concatenating all lines and separating them with commas. See the documentation 
for [GDBshow](#GDBshow) for a detailed explanation of the format and meaning of selections.

When the alignment input is an ALN file, the command must have access to the one or two source files from which 
the ALN file was derived. See [ALNshow](#ALNshow) for a detailed explanation.
When the alignment input is a PAF file, each sequence will be treated as a single contig, even if there are gaps. 
Consequently, selecting specific contigs is not possible in this case.

The default output is an EPS file sent to standard output. With the -p option, the output will be a PDF file, 
requiring a software ('[e]ps[to|2]pdf') to convert the EPS to PDF. You can specify the output PDF file name with 
the -p option; otherwise, the output file name will be determined by the input file name. 

You can use the -l and -i options to filter alignment records based on alignment length and identity. By default, 
only the longest 100,000 alignment records are used for plotting to maintain a manageable file size. This limit 
can be adjusted with the -n option, and setting it to 0 will include all alignments.

The program automatically adjusts the display of the output figure based on the input data. If these automatic 
settings are not suitable, you can use the options -S, -L, -H, -W, -f, and -t to manually configure the display 
parameters.

<a name="addons"></a>


## Additional Utilities

<a name="GDBtoFA"></a>

```
1. GDBtoFA [-vU] [-w<int(80)>  [#[<mask>[.1ano]]]
            <source:path>[.1gdb] [ @ | <target:path>[<fa_extn>|.1seq] ]

       <fa_extn> = (.fa|.fna|.fasta)[.gz]

```
GDBtoFA can produce exactly the FASTA file contents from which the source GDB was derived by
a call to FAtoGDB, or if so directed a ONEcode .1seq file, i.e. it is an inverse operation for
FAtoGDB.    When a GDB is built it records internally where the source FASTA or ONEcode sequence
file is, it's name, and it's extension.  We call this the "origin" in what
follows.  Where the FASTA file is placed by GDBtoFA and what it is named is as follows:

* If there is no target, then the output is streamed to the standard output (uncompressed).

* If the target is the special symbol @, then the FASTA file is built at it's origin directory and given it's origin name and extension (implying compression if it ends with .gz).

* If the target is a directory, then the FASTA file is built at said directory, it's name and extension are as for the origin.

* If the target is a file (that may not exist), then the FASTA file is built at the directory and named as given by the target.  If the target has an extension then that extension is used, otherwise the extension of the origin is used.

The -U, -w, and # arguments have the same interpretation as for GDBshow.

<a name="BEDtoANO"></a>

```
2. BEDtoANO [-T<int(8)>] <bed:path>[.bed] <genome:path>[.1gdb|<fa_extn>|<1_extn>]

           <fa_extn> = (.fa|.fna|.fasta)[.gz]
           <1_extn>  = any valid 1-code sequence file type

      -T: Number of threads to use.
```

Convert a BED file into a ONEcode ANO file.  An ANO file can capture all the information in a BED file
save for the fields denoting display information such a thickness and color and exon clustering.  You
must supply the genome to which the BED information applies.

<a name="ANOtoBED"></a>

```
3. ANOtoBED [-v] <source:path>[.1ano] [ <target:path>[.bed] ]
```

Convert a ONEcode ANO file into a BED file.  If no target is supplied then the output is streamed to
standard output.

<a name="PAFtoALN"></a>

```
4. PAFtoALN [-T<int(8)>] <alignments:path>[.paf]
                         <source1:path>[.gdb|<fa_extn>|<1_extn>] [<source2:path>[.gdb|<fa_extn>|<1_extn>]]
                                     
       <fa_extn> = (.fa|.fna|.fasta)[.gz]
       <1_extn>  = any valid 1-code sequence file type
```

PAFtoALN takes a PAF file as its first argument and the two sources that were compared to produce the
alignment in the PAF file as the second and third arguments.
The PAF file must have CIGAR strings that use X and = to describe the alignment of ungapped segments
as opposed to just M.
The number of threads used is 8 by default, but can be set with the -T option.
A .1aln file with the name \<alignments\>.1aln is produced.  (*NB*: CIGAR strings with M could be
accepted but would require the explicit reconstruction of each alignment, a possible to do if
requested/required.)

<a name="PAFtoPSL"></a>

```
5. PAFtoPSL [-T<int(8)>] [-C<str(cg:Z:)>] <alignments:path>[.paf]
```

PAFtoPSL takes an uncompressed PAF file as input. The input file must include CIGAR strings. Supported CIGAR operators include 'M', '=', 'X', 'I', and 'D'. By default, the CIGAR string tag is 'cg:Z:', but you can 
specify a different tag using the '-C' option. The custom tag must be a five-character string where the third 
and fifth characters are ':'. The number of threads used is 8 by default, but can be set with the -T option.
Output is streamed directly to STDOUT.

<a name="GIXrm"></a>
<a name="GIXcp"></a>
<a name="GIXmv"></a>

```
6.a GIXrm [-vifg] <source:path>[.gix|.1gdb] ...
6.b GIXmv [-vinf] <source:path>[.gix|.1gdb] <target:path>[.gix|.1gdb]
6.c GIXcp [-vinf] <source:path>[.gix|.1gdb] <target:path>[.gix|.1gdb]
```

A GDB consists of not only a "skeleton" file with a .1gdb extension but also a hidden file with extension .bps.
Likewise a GIX consists of not only a proxy file with a .gix extension but 2T hidden files with
extensions of the form .ktab.\<int\> and .post.\<int\> where T is the number of threads used to
create the index.
As such it is cumbersome to remove, move, or copy a GDB or GIX directly with the UNIX OS as it requires
you to utter 2T+3 commands or possibly only one if wild cards are used, albeit this has the potential for
surprising matches affecting unexpected files.
So we provide the commands GIXrm, GIXmv, and GIXcp to remove, move, or copy GDBs and/or GIXs as a
single entity.

The routines operate exactly as for rm, mv, and cp including the meaning/effect of the flags -v, -i, -n,
and -f.  For the GIXmv and GIXcp commands both the GIX and GDB are moved or copied, where just the
GDB will be affected if there is no associated GIX.  On the otherhand, for GIXrm only the GIX will be
deleted unless the -g flag is explicitly set in which case the GDB will also be removed.  We chose
this requirement as a safeguard because if you have replaced your FASTAs with the GDBs as we recommend,
then deleting the GDB is tantamount to deleting your genome source! 

<a name="ALNchain"></a>

```
7. ALNchain [-v] [-g<int(10000)>] [-l<int(10000)>] [-p<float(0.1)>] [-q<float(0.1)>]
                [-z<int(1000)>] [-s<int(10000)>] [-n<int(1)>] [-c<float(0.5)>] [-e<0.0>]
                [-f<int(1000)>] [-o<output:path>[.1aln]] <alignments:path>[.1aln]
```

For each pair of sequences, ALNchain generates a subset of alignments to achieve a one-to-one global alignment 
allowing rearrangements by selecting the best-scored local chains, adhering to user-specified constraints. 
We use a linear gap penalty for chaining, where the cost of a gap or overlap between consecutive alignments in 
the chain is defined by -p or -q. The maximum sizes of gaps and overlaps allowed in the chain can be adjusted 
using the -g and -l options. Chains are scored as *C-G\*p-O\*q*, where *C* represents the total number of unique 
sequence positions covered by the alignments. A chain is terminated if its score drops by more than -z.

Chains are selected based on their scores, from highest to lowest. The -s and -n options specify the minimum 
requirements for a chain to be considered. We track the sequence positions covered by all selected chains. For 
any new chain, the number of additional positions it covers on the sequences is calculated. If this number is 
below certain thresholds, determined by -c as a fraction of the chain size and -e as a fraction of the sequence 
size, the chain is not selected. When calculating the sequence positions covered by chains, the -f option is used 
as the upper limit for closing gaps.

<a name="ALNreset"></a>

```
8. ALNreset [-T<int(8)>] <alignments:path>[.1aln]
                 <source1:path>[.1gdb|<fa_extn>|<1_extn>] [<source2:path>[.1gdb|<fa_extn>|<1_extn>]]
                                     
       <fa_extn> = (.fa|.fna|.fasta)[.gz]
       <1_extn>  = any valid 1-code sequence file type
```

In the unfortunate event that the internal references of a 1-code alignment file (.1aln) to its genome source files have become
"stale", ALNreset allows you to reset these paths within the given file.  Note carefully, that the references can be not only to a GDB but also the source 1-code or FASTA files from which a GDB can be
built.


<a name="ONEaln"></a>

## C-Library for Accessing .1aln Files

This package is intended to give one a simple interface to read and manipulate the
information in a .1aln file of alignment records produced by FastGA or FasTAN.
A simple use pattern for reading all the alignment records in a file is as follows:

```
     AlnReader *reader = alnOpenReader(filename,1,true);
     if (reader == NULL)
       printf("Error: %s\n",alnError());
     while (!alnEOF(reader))
       { AlnRecord *align = alnAlignment(reader,true);
         if (align == NULL)
           printf("Error: %s\n",alnError());

         // do stuff with current alignment record align, e.g.

         alnShowAlignment(align,stdout,8,100,10,5,false,false);

         alnNext(reader);
       }
     alnCloseReader(reader);
```

This library of routines is in ONEaln.c
with the interface declared in ONEaln.h.
Several of the modules used by FastGA must also be compiled in, namely, GDB.[ch],
ONElib.[ch], alncode.[ch], align.[ch], and gene_core.[ch].  See the make command for ONEalnTEST in
the Makefile for an example.


<a name="ehandler"></a>

### ERROR HANDLING:

Two macro variables CHECK\_ARGS and HALT\_ON\_ERROR found at the top of the ONEaln.c file
determine which errors are detected and how they are handled according to whether
the variable is defined or undefined.  In the copy you downloaded they are both
defined, change them before compilation to get the desired behavior.

When CHECK\_ARGS is defined, all arguments are checked, e.g. the range of an index, or
the state of a reader.  This creates overhead for simple routines like GetContigLen
so if your application is well debugged and ready for release, we suggest you undef
the variable, so that only errors not under you control, e.g. when opening or reading a
file, are caught.  In the error listing at the end of this document, those routines whose
return value must be checked regardless of the setting of CHECK\_ARGS are noted.

When HALT\_ON\_ERROR is undefined, the detection of an error results in an error message
being placed in a special buffer you can access with alnError(), and the routine in
question returns with a documented error value.  But if your program is not interactive,
you can, for convenience, define HALT\_ON\_ERROR in which case the detection of an error
results in the error message being written to stderr and your program halting with
exit value 1.  In this case you do not need to check any of the return values of the routines in
the package.

When routines return with an error value, you can access an error message by calling
alnError.  The returned pointer is to a globally shared message string, so its
contents are determined by the last routine to detect an error and write into it.
A list of all the error messages emited by this library are found at the end of this
document.

```
     char *alnError()
```

<a name="areader"></a>

### .1ALN FILE READER:

```
     typedef void *AlnReader   //  A multi-thread .1aln reader
```

Open the .1aln file 'name' for reading with nthreads threads.  NULL is returned
if there is an error, otherwise the return pointer points to the first element
in an array of nthreads AlnReader's.  This first element is called the master.
If you want to have CIGAR strings or print alignments or otherwise need access
to the sequence of the source genomes, you must set see_seq to true.  Otherwise
set it to false for better memory efficiency.

```
     AlnReader *alnOpenReader(char *name, int nthreads, bool see_seq)
```

Each reader conceptually has a "cursor" pointing at an alignment record that can
be advanced one record at a time with alnNext or jumped to the idx'th record
in the file with alnGoto.  To be crystal clear, alnGoto(..,1) takes you to the
first alignment record, i.e. indexing begins at 1.  alnNext returns true if and
only if the reader has been advanced to the end of the file.  AlnGoto returns
true if and only if it was successful, otherwise an error has occurred.

```
     bool alnNext(AlnReader *reader)

     bool alnGoto(AlnReader *reader, int idx)
```

Return true iff at the end of file

```
     bool alnEOF(AlnReader *reader)

```

Close all readers associated with the supplied master reader.

```
     void alnCloseReader(AlnReader *reader)
```

Return (a) the total number of alignments in the file, (b) the maximum number of
trace intervals in any alignment record, (c) the total number of trace intervals
in the file, and (d) the trace point spacing, respectively.

```
     int  alnCount       (AlnReader *reader)
     int  alnTraceMax    (AlnReader *reader)
     int  alnTraceCount  (AlnReader *reader)
     int  alnTraceSpacing(AlnReader *reader)
```

<a name="gdb"></a>

### GENOME DATA BASE (GDB):

```
     typedef void AlnGDB    //  A genome database
```

From the .1aln reader get the first and second genome data bases over which
the alignments were found.  The two pointers are equal if there was only one genome,
i.e. a self-comparison.  A genome data base gives one the structure of the genome
as a collection of scaffolded contigs separated by gaps, as well as access to
the contig's sequence if 'see_seq' was true in the call that created the reader.

```
     AlnGDB *alnGDB1(AlnReader *reader)
     AlnGDB *alnGDB2(AlnReader *reader);

```

The Count... routines return the # of (a) scaffolds, (b) contigs, and (c) gaps in a GDB.
The Max... routines return the maximum # of contigs/gaps in any scaffold of a GDB.

```
     int gdbScaffoldCount(AlnGDB *g)
     int gdbContigCount  (AlnGDB *g)
     int gdbGapCount     (AlnGDB *g)
     int gdbContigMax    (AlnGDB *g)
     int gdbGapMax       (AlnGDB *g)
```

Scaffolds are numbered starting at 1, and the first contig of a scaffold is indexed by 1.
In rare, somewhat abnormal cases there can be a run of N's at the start or end
of a scaffold, the length of these are obtained by gdbGapLen(g,s,0) and gdbGapLen(g,s,gdbScaffoldContigs(g,s)).
 The routines return -1 if an error occurs.


```
     int gdbScaffoldLen    (AlnGDB *g, int s)        // Length of the s'th scaffold
     int gdbScaffoldContigs(AlnGDB *g, int s)        // Number of contigs in the s'th scaffold
     int gdbGapLen         (AlnGDB *g, int s, int p) // Length of the p'th gap of the s'th scaffold
     int gdbContigLen      (AlnGDB *g, int s, int c) // Length of the c'th contig of the s'th scaffold
     int gdbContigStart    (AlnGDB *g, int s, int c) // Start position of the c'th contig of the
                                                     //   s'th scaffold (see getScaffoldSeq)
```

Get the scaffold name of the s'th scaffold.  The string returned is in a buffer internal
to the GDB and is overwritten each time this routine is called.  You should copy it to
memory you control if you wish it to persist beyond the last call.  NULL is returned
if an error occurs.

```
     char *gdbScaffoldName(AlnGDB *g, int s)
```

Get the sequence spanning interval [beg,end] of the s'th scaffold of g.  If buffer = NULL
then the routine allocates an array of (end-beg)+1 bytes, places the requested sequence '\0'-terminated
there, and returns a pointer to it.  In this case the user must subsequently free this
string.  Otherwise, if the user supplies a non-NULL buffer pointer, then it must be of
length not less than (end-beg)+1, the
 requested sequence is placed there '\0'-terminated, and the buffer
pointer is returned.  NULL is returned on an error which includes not
asking for sequence
access when opening the reader and requesting an interval not wholly within a contig.

Sequences positions are *between* base pairs, e.g. position 0 is before the first bp and
position 1 is between the first and second bp.  In this way the subsequence between
an interval [a,b] is of length b-a and consists of the a+1'st through b'th bp of the
selected sequence.

```
char *gdbScaffoldSeq(AlnGDB *g, int s, int beg, int end, char *buffer)
```

<a name="arecord"></a>

### ALIGNMENT RECORDS:

```
     typedef struct
       { int  seq1,  seq2;    //  seq1[bpos1..epos1] aligns to seq2[bpos2..epos2]
         int  bpos1, epos1;   //  with diffs differences (indels & substitutions).
         int  bpos2, epos2;   //
         int  diffs;          //

         int  tlen;           //  tpoints[0..tlen) contains the inter trace point lengths in bseq
         int *tpoints;        //  tdiffs[0..tlen) is the # of diffs in each tp interval
         int *tdiffs;
       } AlnRecord;
```

Load the alignment record at the reader's current cursor/position.  The return pointer
points at a pre-allocated buffer internal to the reader that is reused/overwritten
each time alnAlignment is called.  This includes the memory for the tpoints and tdiffs
arrays of the trace.  You should copy it and the two arrays just mentioned if you
wish them to persist beyond the last call.  If see_trace is true then the trace point
arrays are sought and loaded so one can create cigar strings, etc.  NULL is returned
on error.

```
     AlnRecord *alnAlignment(AlnReader *reader, bool see_trace)
```

Create a CIGAR string for the the given ALN record.  This is only possible if see\_seq
was true when the alignment reader was opened, and see\_trace was true when the alignment
record was fetched.  For this routine, the returned string must be freed by the user.
NULL is returned on an error.

If show_x is true then 'X' and '=' are used to model stretches without indels, otherwise
'M' is used.  When reversed is false, the CIGAR string is one that transforms the first
sequence into the second where the first string is in the forward orientation.  If
reversed is true, then the CIGAR string is one that transforms the second sequence
into the first where the second sequence is in the forward orientation, i.e. the roles of
the first and second sequence are reversed.

```
     char *alnCreateCigar(AlnRecord *align, bool show_x, bool reversed)
```

Create a CStag string for the the given ALN record.  This is only possible if see\_seq
was true when the alignment reader was opened, and see\_trace was true when the alignment
record was fetched.  For this routine, the returned string must be freed by the user.
NULL is returned on an error.

If short\_form is false then the normal CStag is created where the sequence of a matching
segment is given, otherwise just the length of a matching segment is given.
When reversed is false, the CS tag is one that transforms the first sequence
into the second where the first string is in the forward orientation.  If reversed
is false, then the CS tag is one that transforms the second sequence into the
first where the second sequence is in the forward orientation, i.e. the roles of the first and second sequence are reversed.
 

```
     char *alnCreateCStag(AlnRecord *align, bool short_form, bool reversed)
```

Create an "indel array" for the the given ALN record.  This is only possible if see\_seq
was true when the alignment reader was opened, and see\_trace was true when the alignment
record was fetched.  For this routine, the returned integer is in a preallocated buffer
of the reader, so you must make a copy if you wish it to persist beyond the next call to
the package. 
NULL is returned on an error.

An indel array gives the locations at which a dash (-) should be inserted into either the
first or second sequence in order to expose the alignment between the two sequences
implied by the trace point intervals.  A postive number x indicates one should place a
dash before the x'th character of the first sequence, and a negative number -x indicates
placing a dash before the x'th character of the second sequence.  The absolute values
of the dash locations is increasing along the array and the array is terminated by a 0 value.

The indices in the indel array are relative to the (sub)sequences that are aligned.
This includes complementing the B/second sequence if necessary, so that the position 1 refers
to the first character in the alignment for both sequences.

When reversed is false, the indel array is one that transforms the first sequence
into the second where the first string is in the forward orientation.  If reversed
is true, then the indel array is one that transforms the second sequence into the
first where the second sequence is in the forward orientation, i.e. the roles of the first and second sequence are reversed.


```
     int  *alnCreateIndelArray(AlnRecord *align, bool reversed)
```

Write a BLAST-like text display of an alignment to FILE where.  For this routine to work
it must be that see\_seq was true when the alignment reader was opened, and see\_trace
was true when the alignment record was fetched.  True is returned on error.
  
The display is indented by indent spaces, each displayed segment shows width columns of
the alignment, border bp's before and after the aligned portion are also displayed, and
the display width for coordinates is given by coord.  Upper_case controls the case of
the base pairs.

When reversed is false, the 1st sequence is displayed above the second with the first
sequence in the forward direction.  When reversed is true, the 2nd sequence is displayed
above the first with the second sequence in the forward direction, i.e. the roles of
the first and second sequence are reversed.

```
     bool  alnShowAlignment(AlnRecord *align, FILE *where,
                            int indent, int width, int border, int coord,
                            bool upper_case, bool reversed)
```

<a name="emessage"></a>

### ERROR MESSAGES:

Find below a list of all error message followed by a list of the routines that can
potentially encounter the error.  Starred error message are controlled by CHECK_ARGS
and so do not occur if this macro variable is set to undefined.

```
    Could not open <path>/<root>.1aln
    Out of memory (Opening .1aln file)
        alnOpenReader

  * Index out of range
  * Could not seek to <idx>th alignment
        alnGoto
        
  * Scaffold index out of range
        gdbScaffoldLen
        gdbContigLen
        gdbGapLen
        gdbScaffoldContigs
        gdbScaffoldStart
        gdbScaffoldName
        gdbScaffoldSeq

  * Contig index out of range
        gdbContigLen
        gdbScaffoldStart

  * Gap index out of range
        gdbGapLen

  * Sequence access not requested on open
    Could not seek GDB sequence file
    Could not read GDB sequence file
        gdbScaffoldSeq
        alnCreateCigar
        alnCreateCStag
        alnCreateIndelArray
        alnShowAlignment

  * Scaffold interval intersects a gap
    Out of memory (Allocating sequence)
        gdbScaffoldSeq

    T-line not present in alignment record
    T- and X-lines do not have the same length
    D-line not present in alignment record
        alnAlignment

    Bad alignment between trace points
    Trace point out of bounds
    Alignment end points not in band
    Self comparison can cross main diagonal
    Out of memory (Enlarging trace vector)
    Out of memory (Enlarging DP vector)
  * Trace information was not requested when record fetched
        alnCreateCigar
        alnCreateCStag
        alnCreateIndelArray
        alnShowAlignment

    Out of memory (Allocating CIGAR array)
    Out of memory (Expanding CIGAR array)
        alnCreateCigar
        alnCreateCStag

    Out of memory (Allocating CIGAR string)
        alnCreateCigar

    Out of memory (Allocating CS-tag string)
        alnCreateCStag
```

In summary the following routines can return an error even if CHECK_ARGS is undefined, as they
can encounter a read/seek failure, a memory allocation failure, or missing fields in a .1aln file

```
     alnOpenReader
     alnAlignment
     gdbScaffoldSeq (iff buffer == NULL)
     alnCreateCigar
     alnCreateCStag
     alnCreateIndelArray
     alnShowAlignment
```




