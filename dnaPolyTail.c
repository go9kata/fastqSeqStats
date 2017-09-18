//
//  dnaPolyTail.c
//  fastqSeqStats
//
//  Created by Georgi Tushev on 23/12/14.
//  Copyright (c) 2014 Georgi Tushev. All rights reserved.
//

#include "dnaPolyTail.h"

int findTailPolyAMaybeMask(DNA *dna, int size, boolean loose)
/* Identify PolyA at end; This allows a few
 * non-A's as noise to be trimmed too.  Returns number of bases trimmed.
 * Leaves first two bases of PolyA in case there's a taa stop codon. */
{
    int i;
    int score = 10;
    int bestScore = 10;
    int bestPos = -1;
    int trimSize = 0;
    
    for (i=size-1; i>=0; --i)
    {
        DNA b = dna[i];
        if (b == 'n' || b == 'N')
            continue;
        if (score > 20) score = 20;
        if (b == 'a' || b == 'A')
        {
            score += 1;
            if (score >= bestScore)
            {
                bestScore = score;
                bestPos = i;
            }
            else if (loose && score >= (bestScore - 8))
            {
                /* If loose, keep extending even if score isn't back up to best. */
                bestPos = i;
            }
        }
        else
        {
            score -= 10;
        }
        if (score < 0)
        {
            break;
        }
    }
    if (bestPos >= 0)
    {
        trimSize = size - bestPos - 2;	// Leave two for aa in taa stop codon
        if (trimSize < 0)
            trimSize = 0;
    }
    return trimSize;
}

int tailPolyASizeLoose(DNA *dna, int size)
/* Return size of PolyA at end (if present).  This allows a few non-A's as
 * noise to be trimmed too, but skips first two aa for taa stop codon.
 * It is less conservative in extending the polyA region than maskTailPolyA. */
{
    return findTailPolyAMaybeMask(dna, size, TRUE);
}



int findHeadPolyTMaybeMask(DNA *dna, int size, boolean loose)
/* Return size of PolyT at start (if present); mask to 'n' if specified.
 * This allows a few non-T's as noise to be trimmed too, but skips last
 * two tt for revcomp'd taa stop codon. */
{
    int i;
    int score = 10;
    int bestScore = 10;
    int bestPos = -1;
    int trimSize = 0;
    
    for (i=0; i<size; ++i)
    {
        DNA b = dna[i];
        if (b == 'n' || b == 'N')
            continue;
        if (score > 20) score = 20;
        if (b == 't' || b == 'T')
        {
            score += 1;
            if (score >= bestScore)
            {
                bestScore = score;
                bestPos = i;
            }
            else if (loose && score >= (bestScore - 8))
            {
                /* If loose, keep extending even if score isn't back up to best. */
                bestPos = i;
            }
        }
        else
        {
            score -= 10;
        }
        if (score < 0)
        {
            break;
        }
    }
    if (bestPos >= 0)
    {
        trimSize = bestPos+1 - 2;	// Leave two for aa in taa stop codon
        if (trimSize < 0)
            trimSize = 0;
    }
    return trimSize;
}

int headPolyTSizeLoose(DNA *dna, int size)
/* Return size of PolyT at start (if present).  This allows a few non-T's as
 * noise to be trimmed too, but skips last two tt for revcomp'd taa stop
 * codon.
 * It is less conservative in extending the polyA region than maskHeadPolyT. */
{
    return findHeadPolyTMaybeMask(dna, size, TRUE);
}

