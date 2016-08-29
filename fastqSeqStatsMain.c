//
//  fastqSeqStatsMain.c
//  fastqSeqStats
//
//  Created by Georgi Tushev on 22/12/14.
//  Copyright (c) 2014 Georgi Tushev. All rights reserved.
//

#include <stdio.h>
#include "utilities.h"
#include "fastqCounter.h"

int main(int argc, char *argv[])
{
    // configuration variables
    Parameters config;
    Accumulators data;
    
    
    // parse input parameters
    if(parseParameters(argc, argv, &config) != 0)
    {
        help();
        return 0;
    }
    
    // initialize counters
    InitializeAccumulators(&data, &config);
    
    // parsa FASTQ stream
    ParseFastqStream(&data, &config);
    
    // print stats out
    PrintResultStats(&data, &config);
    
    // destroy counters
    DestroyAccumulators(&data, &config);
    
    fprintf(stderr,"line: %lu\nrecords %lu\n", data.countLines, data.countRecords);
    
    return 0;
}