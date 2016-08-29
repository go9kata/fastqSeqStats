//
//  fastqCounter.h
//  fastqSeqStats
//
//  Created by Georgi Tushev on 23/12/14.
//  Copyright (c) 2014 Georgi Tushev. All rights reserved.
//

#ifndef __fastqSeqStats__fastqCounter__
#define __fastqSeqStats__fastqCounter__

#include <stdio.h>
#include <math.h>
#include "utilities.h"
#include "dnaPolyTail.h"

typedef struct
{
    // incrementing counters
    unsigned long int countLines;
    unsigned long int countRecords;
    
    // histogram of Phred score values
    unsigned long int *countPhred;
    
    // Phred score pre sequence length
    double            *avgPhredScore;
    double            *varPhredScore;
    
    // nucleotides per sequence length
    unsigned long int *countAny;
    unsigned long int *countAs;
    unsigned long int *countCs;
    unsigned long int *countGs;
    unsigned long int *countTs;
    unsigned long int *countNs;
    
    // span of poly(A/T)s per sequence length
    unsigned long int *spanPolyAs;
    unsigned long int *spanHeadTs;
    
} Accumulators;

// functions primitives
int InitializeAccumulators(Accumulators *, Parameters *);
/* dynamically initializes all accumulative arrays in the
 * Accumulators struct to 0. Sequence length and number of
 * Phred scores is predifined utilities.h 
 */

int ParseFastqStream(Accumulators *, Parameters *);
/* parsing a FASTQ file stream line by line. Keeps line
 * itterative flag to distinguish between header,sequence and score.
 * Running average and variance of Phred score per nucleotide are calculated.
 */

int PrintResultStats(Accumulators *, Parameters *);
/* prints a tab-delimited text files with the required statistics.
 */

int DestroyAccumulators(Accumulators *, Parameters *);
/* frees all dynamically allocated arrays
 */

int ParseNucleotides(Accumulators *, DNA *);
/* count each base per DNA sequence */

int ParsePolyTails(Accumulators *, DNA *);
/* count loose size of poly tails */

int ParsePhredQualities(Accumulators *, PHRED*, boolean);
/* count each phred quality score per DNA sequence */


#endif /* defined(__fastqSeqStats__fastqCounter__) */