//
//  utilities.c
//  fastqSeqStats
//
//  Created by Georgi Tushev on 22/12/14.
//  Copyright (c) 2014 Georgi Tushev. All rights reserved.
//

#include "utilities.h"

const char VERSION[] = "v1.01";
const char PROGRAM_NAME[] = "fastqSeqStats";

int parseParameters(int argc, char *argv[], Parameters *config)
{
    // initialize default values
    config->source = NULL;
    config->destination = NULL;
    config->polytail = FALSE;
    config->seqlength = MAX_LENGTH;
    
    if(argc < 5)
    {
        fprintf(stderr,"Error:: provide input & output stream!\n");
        return 1;
    }
    
    // do some parsing
    for(int i = 1; i < argc; i++)
    {
        // current paramter length
        int parameterLength = (int)strlen(argv[i]);
        
        // check each parameter and its value
        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
           (PARAMETER_CHECK("--help", 6, parameterLength)))
        {
            help();
        }
        else if((PARAMETER_CHECK("-v", 2, parameterLength)) ||
                (PARAMETER_CHECK("--version", 9, parameterLength)))
        {
            version();
        }
        else if((PARAMETER_CHECK("--polyTail", 10, parameterLength)))
        {
            config->polytail = TRUE;
        }
        else if((PARAMETER_CHECK("-l", 2, parameterLength)) ||
                (PARAMETER_CHECK("--length", 8, parameterLength)))
        {
            i += 1;
            config->seqlength = atoi(argv[i]);
            
            if(config->seqlength <= 0)
            {
                fprintf(stderr,"Error:: sequence length should be >0!\n");
                return 1;
            }
            
        }
        else if((PARAMETER_CHECK("-i", 2, parameterLength)) ||
                (PARAMETER_CHECK("--ifastq", 8, parameterLength)))
        {
            i += 1;
            
            if(strcmp("-", argv[i]) == 0)
            {
                config->source = stdin;
            }
            else
            {
                config->source = fopen(argv[i], "r");
            }
            
            if (config->source == NULL)
            {
                fprintf(stderr,"Error:: could not open file input stream!\n");
                return 1;
            }
        }
        else if((PARAMETER_CHECK("-o", 2, parameterLength)) ||
                (PARAMETER_CHECK("--otxt", 6, parameterLength)))
        {
            i += 1;
            
            if(strcmp("-", argv[i]) == 0)
            {
                config->destination = stdout;
            }
            else
            {
                config->destination = fopen(argv[i], "w");
            }
            
            if (config->destination == NULL)
            {
                fprintf(stderr, "Error:: could not open file output stream!\n");
                return 1;
            }
            
        }
        else
        {
            fprintf(stderr, "Error:: Unknown parameter %s\n", argv[i]);
            return 1;
        }
        
        
    }

    return 0;
}

int version(void)
{
    fprintf(stderr, "fastqSeqStats %s\n", VERSION);
    return 0;
}

int help(void)
{
    
    fprintf(stderr, "\nTool: fastqSeqStats\n");
    fprintf(stderr, "Version: %s\n", VERSION);
    fprintf(stderr, "Summary:\n");
    fprintf(stderr, "\tCreates histogram of base quality score.\n");
    fprintf(stderr, "\tBase composition and quality score per nucleotide are accumulated throughout sequence length.\n");
    fprintf(stderr, "\tCreates histogram of loose count for poly(A) tail or poly(T) head.\n");
    fprintf(stderr, "\nUsage:\t%s [OPTIONS] --ifastq <fastq/stdin> --otxt <txt/stdout>\n", PROGRAM_NAME);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t--polyTail\tcreate histogram of loose poly Tail/Head count, based on Kent et.al. 2002 algorithm\n");
    fprintf(stderr, "\t-h/--help\tprint help message\n");
    fprintf(stderr, "\t-v/--version\tprint current version\n");
    return 1;
}

int chompNewLine(char *line_buf)
{
    size_t  line_length;
    line_length = strlen(line_buf) - 1;
    if(line_buf[line_length] == '\n')
        line_buf[line_length] = '\0';
    
    return 0;
}
