/*******************************************************************************************
 *
 *  Approximate match filter of Rasmussen, Stoye, Myers (RECOMB 2005)
 *
 *  Author:  Gene Myers
 *  Date  :  August 2004
 *
 ********************************************************************************************/

#ifndef _H_FILTER

#define _H_FILTER

/* Filter hit tubes:

     Coordinates are in terms of the d.p. matrix that go from (0,0)
     to (|A|,|B|), where A is the A-sequence argument and B is the
     B-sequence argument.  A coordinate is a position between chars
     of the sequence.  For example, the interval (3,6) designates
     characters 4 through 6 of a sequence.

     Diagonal k is the set of coordinates (a,b) s.t. a-b = k.

     A hit "tube" is a parallelogram of the d.p. matrix between rows bstart and
     bfinish (of the second sequence), and diagonals diag - (BinSize-1)
     through diag.  All MinMatch or longer matches involving MaxError or fewer
     differences are in the set of tubes returned by Match_Filter, but some
     tubes may be false positives, i.e. do not contain such a match.  The
     parameters BinSize, MinMatch, and MaxError can be set by the user
     (see Set_Filter_Params below).
*/

typedef struct {
  int diag;      /* Max diagonal of hit tube */
  int bstart;    /* B position of start of hit tube */
  int bfinish;   /* B position of end of hit tube */
} HitTube;

/* Match_Filter:
     Match_Filter runs the filtration algorithm on A (of length Alen) versus
     B (of length Blen), where A is indexed and B is scanned.  If the A sequence
     was not the same argument in the previous call, then an index of it is
     built, otherwise the index built in the previous call is reused.  If asym
     is non-zero then only hits on negative diagonals are considered, i.e.
     only the lower half of the d.p matrix is searched for matches.  This is
     useful when A and B are the same (set of) sequences or their complements
     and one does not wish to find i vs. j and then later j vs. i.
     
     An array of *NumHits hit tubes containing all the filtration matches is
     returned.  Freeing the space for the array is the responsibility of the
     user and since it can be quite large, we recommend doing so at the first
     opportunity within the application.
*/
   
HitTube *Match_Filter(char *A, int Alen, char *B, int Blen, int asym, int *NumHits);

/* Set_Filter_Params:

     The filter is controlled by a number of parameters that have a complex interdependence.
     See filter_params.h or our paper for more details.  For an appropriately selected
     set of parameters the filter is guaranteed to find all Erate-matches of length MinMatch
     or longer.  Use the functions in filter_params.c to select an appropriate set of
     parameters.  The filter itself only needs the match rate, k-mer length, and minimum
     length given that there is a set of the parameters for the combination (the routine
     checks this).
 
      Core Match:
        Erate     (.05) -- error rate of matches
        KmerLen   ( 13) -- index tuple lengths (=> index size of 4^KmerLen)
        MinLen    ( 47) -- minimum length match to be detected

          ==>  HitThresh is 9, MaxDiffs is 4, and TrapMin is 60, i.e. one seeks
                 60x4 parallelograms containing 9 or more 13-mers.
 
      Repeat Suppression:
        Suppress ( 0) -- if 0 then count all KmerLen-tuples matches, otherwise
                         do not count tuples occuring more than Suppress
                         times towards a hit
 
      Tube Size:
        BinShift ( 3) -- Hits are counted in a "tube" of
                         Binsize = 2^BinShift + max_diffs(Erate,KmerLen,HitThresh) (=12) diagonals
 
    Apart from determining what matches are guaranteed to be found, these parameters
    also determine time and space performance as follows:
 
       Space: 4*(4^KmerLen) + 4*Alen bytes for the index of the A sequence,
              12*Alen / BinSize bytes working storage, and
              12*NumHits bytes for the returned array.
 
       Time:  O(Alen * (1 + Blen / 4^KmerLen) + NumHits).  The specificity of the
              filtration is roughly BinSize / 4^(KmerLen+HitThresh), i.e.
              the number of false positives contributions to NumHits is roughly
              this specificity times Alen*Blen.
 
    The parameters must satisfy a number of constraints including:
       a. MinLen >= KmerLen > 0
       b. Erate in [0.,1.)
       c. Kmerlen < ceiling( 1 / Erate )
       d. HitThresh = max_hits(Erate,KmerLen,MinLen) >= 1
       e. 2^BinShift >= max_diffs(Erate,KmerLen,HitThresh)
    If one or more of these conditions is not met then Set_Filter_Params returns
    a non-zero value and does not change the current parameter settings.
 
    Not counting overrepresented Kmerlen-tuples towards a filter hit prevents low
    complexity and ubiquitous interspersed repeats from producing hits in the event
    one's primary interest is in aligning the unique regions of the two sequences.
    The parameter Suppress should be large enough that with very high probability
    all shared tuples between uniquely matching segments are not suppressed.
 
    The size of a tube is specificied by BinShift in order to be able to use
    bit-shifting to replace multiplies and divides -- the wider the tube the less
    specific the filter.  In principle tubes need only be max_diffs(Erate,KmerLen,HitThresh)
    diagonals wide, but we use the larger BinSize in order to reduce memory requirements
    without significantly compromising specificity.  The default value, 3, is a good
    compromise.
*/

int Set_Filter_Params(double Erate, int MinLen, int KmerLen, int BinShift, int Suppress); 

/* Reset_Filter:

     All space for the A sequence index is freed.  Normally this index is saved by
     Match_Filter in case the next call requires the same index.  To force an index
     to be created prior to real scans, call Match_Filter with B = NULL.
*/

void Reset_Filter();

#endif // _H_FILTER
