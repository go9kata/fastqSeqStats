//
//  fastqCounter.c
//  fastqSeqStats
//
//  Created by Georgi Tushev on 23/12/14.
//  Copyright (c) 2014 Georgi Tushev. All rights reserved.
//

#include "fastqCounter.h"

int InitializeAccumulators(Accumulators *data, Parameters *config)
/* dynamically initializes all accumulative arrays in the
 * Accumulators struct to 0. Sequence length and number of
 * Phred scores is predifined utilities.h
 */
{
    // incrementing counters
    data->countLines = 0;
    data->countRecords = 0;
    
    // histogram of Phred score values
    data->countPhred = (unsigned long int *)calloc(PHRED_MAX, sizeof(unsigned long int));
    
    // Phred score per sequence length
    data->avgPhredScore = (double *)calloc(config->seqlength, sizeof(double));
    data->varPhredScore = (double *)calloc(config->seqlength, sizeof(double));
    
    // nucleotides per sequence length
    data->countAny = (unsigned long int *)calloc(config->seqlength, sizeof(unsigned long int));
    data->countAs = (unsigned long int *)calloc(config->seqlength, sizeof(unsigned long int));
    data->countCs = (unsigned long int *)calloc(config->seqlength, sizeof(unsigned long int));
    data->countGs = (unsigned long int *)calloc(config->seqlength, sizeof(unsigned long int));
    data->countTs = (unsigned long int *)calloc(config->seqlength, sizeof(unsigned long int));
    data->countNs = (unsigned long int *)calloc(config->seqlength, sizeof(unsigned long int));
    
    // span of poly(A/T)s per sequence length
    if(config->polytail == TRUE)
    {
        data->spanPolyAs = (unsigned long int *)calloc(config->seqlength, sizeof(unsigned long int));
        data->spanHeadTs = (unsigned long int *)calloc(config->seqlength, sizeof(unsigned long int));
    }
    else
    {
        data->spanPolyAs = NULL;
        data->spanHeadTs = NULL;
    }
    
    return 0;
}

int ParseFastqStream(Accumulators *data, Parameters *config)
/* parsing a FASTQ file stream line by line. Keeps line
 * itterative flag to distinguish between header,sequence and score.
 * Running average and variance of Phred score per nucleotide are calculated.
 */
{
    char    line_buf[MAX_LENGTH];
    int     fastq_flag = 0;
    boolean first_record = FALSE;
    
    
    // loop through input stream
    while(fgets(line_buf, sizeof(line_buf), config->source))
    {
        // update line and flag itterators
        data->countLines++;
        fastq_flag++;
        
        // skip header and 3rd FASTQ line
        if((fastq_flag == 1) || (fastq_flag == 3))
            continue;
        
        // process nucleotide counts from 2nd FASTQ line
        if(fastq_flag == 2)
        {
            chompNewLine(line_buf);
            ParseNucleotides(data, line_buf);
            
            if(config->polytail == TRUE)
            {
                ParsePolyTails(data, line_buf);
            }
        }
        
        // process nucleotide counts from 4th FASTQ line
        if(fastq_flag == 4)
        {
            // update FASTQ record
            data->countRecords++;
            fastq_flag = 0;
            first_record = (data->countRecords == 1) ? TRUE : FALSE;
            
            // parse quality score
            chompNewLine(line_buf);
            ParsePhredQualities(data, line_buf, first_record);
        }
        
    }
    
    return 0;
}

int PrintResultStats(Accumulators *data, Parameters *config)
/* prints a tab-delimited text files with the required statistics.
 */
{
    /* print out counters */
    fprintf(config->destination, "#number of lines\tnumber of fastq records\n");
    fprintf(config->destination, "%lu\t%lu\n",data->countLines, data->countRecords);
    
    /* print out histogram */
    fprintf(config->destination, "#Histogram: phred score per base frequency\n");
    for (int k = 0; k < PHRED_MAX; k++)
    {
        if (data->countPhred[k] > 0)
        {
            fprintf(config->destination, "%d\t%lu\n", k, data->countPhred[k]);
        }
    }
    
    /* print out per length stats */
    fprintf(config->destination, "#Stats per sequence length\n");
    fprintf(config->destination, "#base\tAs\tCs\tGs\tTs\tNs\t");
    if(config->polytail == TRUE)
    {
        fprintf(config->destination, "polyTailAs\tpolyHeadTs\t");
    }
    fprintf(config->destination, "avgQ\tvarQ\n");
    
    for (int k = 0; k < config->seqlength; k++)
    {
        if (data->countAny[k] > 0)
        {
            fprintf(config->destination, "%d\t", k);
            //fprintf(config->destination, "%lu\t", data->countAny[k]);
            //fprintf(config->destination, "%lu\t", data->countAs[k]);
            //fprintf(config->destination, "%lu\t", data->countCs[k]);
            //fprintf(config->destination, "%lu\t", data->countGs[k]);
            //fprintf(config->destination, "%lu\t", data->countTs[k]);
            //fprintf(config->destination, "%lu\t", data->countNs[k]);
            
            fprintf(config->destination, "%.4f\t", (double)data->countAs[k]/data->countAny[k]);
            fprintf(config->destination, "%.4f\t", (double)data->countCs[k]/data->countAny[k]);
            fprintf(config->destination, "%.4f\t", (double)data->countGs[k]/data->countAny[k]);
            fprintf(config->destination, "%.4f\t", (double)data->countTs[k]/data->countAny[k]);
            fprintf(config->destination, "%.4f\t", (double)data->countNs[k]/data->countAny[k]);
            if(config->polytail == TRUE)
            {
                fprintf(config->destination, "%.4f\t", (double)data->spanPolyAs[k]/data->countAny[k]);
                fprintf(config->destination, "%.4f\t", (double)data->spanHeadTs[k]/data->countAny[k]);
            }
            fprintf(config->destination, "%.4f\t", data->avgPhredScore[k]);
            fprintf(config->destination, "%.4f\n", (double)(sqrt(data->varPhredScore[k])/sqrt(data->countAny[k])));
        }
    }
    
    
    return 0;
}


int DestroyAccumulators(Accumulators *data, Parameters *config)
/* frees all dynamically allocated arrays
 */
{
    // histogram of Phred score values
    if(data->countPhred != NULL)
    {
        free(data->countPhred);
        data->countPhred = NULL;
    }
    
    // Phred score per sequence length
    if(data->avgPhredScore != NULL)
    {
        free(data->avgPhredScore);
        data->avgPhredScore = NULL;
    }
    if(data->varPhredScore != NULL)
    {
        free(data->varPhredScore);
        data->varPhredScore = NULL;
    }
    
    // nucleotides per sequence length
    if(data->countAny != NULL)
    {
        free(data->countAny);
        data->countAny = NULL;
    }
    if(data->countAs != NULL)
    {
        free(data->countAs);
        data->countAs = NULL;
    }
    if(data->countCs != NULL)
    {
        free(data->countCs);
        data->countCs = NULL;
    }
    if(data->countGs != NULL)
    {
        free(data->countGs);
        data->countGs = NULL;
    }
    if(data->countTs != NULL)
    {
        free(data->countTs);
        data->countTs = NULL;
    }
    if(data->countNs != NULL)
    {
        free(data->countNs);
        data->countNs = NULL;
    }
    
    // span of poly(A/T)s per sequence length
    if(data->spanPolyAs != NULL)
    {
        free(data->spanPolyAs);
        data->spanPolyAs = NULL;
    }
    if(data->spanHeadTs != NULL)
    {
        free(data->spanHeadTs);
        data->spanHeadTs = NULL;
    }
    
    // close open file pointers
    fclose(config->source);
    fclose(config->destination);
    
    return 0;
}


int ParseNucleotides(Accumulators *data, DNA *dna)
/* count each base per DNA sequence */
{
    char    dna_base;
    size_t  dna_length = strlen(dna);
    
    for(int i = 0; i < dna_length; i++)
    {
        // current base
        dna_base = dna[i];
        
        // update counting arrays
        data->countAny[i]++;
        switch(dna_base)
        {
            case 'A':
                data->countAs[i]++;
                break;
            case 'a':
                data->countAs[i]++;
                break;
            case 'C':
                data->countCs[i]++;
                break;
            case 'c':
                data->countCs[i]++;
                break;
            case 'G':
                data->countGs[i]++;
                break;
            case 'g':
                data->countGs[i]++;
                break;
            case 'T':
                data->countTs[i]++;
                break;
            case 't':
                data->countTs[i]++;
                break;
            default:
                data->countNs[i]++;
                break;
        }
    }
    return 0;
}


int ParsePolyTails(Accumulators *data, DNA *dna)
/* count loose size of poly tails */
{
    int     spanPolyTailA = 0;
    int     spanPolyHeadT = 0;
    int     dna_length = (int)strlen(dna);
    
    // calculate current poly size
    spanPolyTailA = tailPolyASizeLoose(dna, dna_length);
    spanPolyHeadT = headPolyTSizeLoose(dna, dna_length);
        
    // update histogram
    data->spanPolyAs[spanPolyTailA]++;
    data->spanHeadTs[spanPolyHeadT]++;
    
    return 0;
}

int ParsePhredQualities(Accumulators *data, PHRED *score, boolean FirstRecord)
/* count each phred quality score per DNA sequence */
{
    size_t score_length = strlen(score);
    int    base_score;
    double curr_avg_score;
    double curr_var_score;
    
    /* loop through each base */
    for (int i = 0; i < score_length; i++)
    {
        /* current base score */
        base_score = score[i];
        
        /* correct range of base_score */
        if (base_score < 0)
            base_score = 0;
        if (base_score > PHRED_MAX)
            base_score = PHRED_MAX;
        
        /* update base histogram */
        data->countPhred[base_score]++;
        
        /* calculate per length statistics */
        if (FirstRecord == TRUE)
        {
            data->avgPhredScore[i] = (double)base_score;
            data->varPhredScore[i] = (double)(0.0);
        }
        else
        {
            curr_avg_score = data->avgPhredScore[i] + (base_score - data->avgPhredScore[i])/data->countAny[i];
            curr_var_score = data->varPhredScore[i] + (base_score - data->avgPhredScore[i])*(base_score - curr_avg_score);
            
            data->avgPhredScore[i] = curr_avg_score;
            data->varPhredScore[i] = curr_var_score;
        }
    }
    
    return 0;
}