//
//  dnaPolyTail.h
//  fastqSeqStats
//
//  Created by Georgi Tushev on 23/12/14.
//  Copyright (c) 2014 Georgi Tushev. All rights reserved.
//

#ifndef __fastqSeqStats__dnaPolyTail__
#define __fastqSeqStats__dnaPolyTail__

#include <stdio.h>
#include "utilities.h"

// Poly Tail Algorithm Acquired From KENT UCSC source code

int findTailPolyAMaybeMask(DNA *dna, int size, boolean loose);
/* Identify PolyA at end; mask to 'n' if specified.  This allows a few
 * non-A's as noise to be trimmed too.  Returns number of bases trimmed.
 * Leaves first two bases of PolyA in case there's a taa stop codon. */

int findHeadPolyTMaybeMask(DNA *dna, int size, boolean loose);
/* Return size of PolyT at start (if present); mask to 'n' if specified.
 * This allows a few non-T's as noise to be trimmed too, but skips last
 * two tt for revcomp'd taa stop codon. */

int tailPolyASizeLoose(DNA *dna, int size);
/* Return size of PolyA at end (if present).  This allows a few non-A's as
 * noise to be trimmed too, but skips first two aa for taa stop codon.
 * It is less conservative in extending the polyA region than maskTailPolyA. */

int headPolyTSizeLoose(DNA *dna, int size);
/* Return size of PolyT at start (if present).  This allows a few non-T's as
 * noise to be trimmed too, but skips last two tt for revcomp'd taa stop
 * codon.
 * It is less conservative in extending the polyA region than maskHeadPolyT. */

#endif /* defined(__fastqSeqStats__dnaPolyTail__) */