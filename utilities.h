//
//  utilities.h
//  fastqSeqStats
//
//  Created by Georgi Tushev on 22/12/14.
//  Copyright (c) 2014 Georgi Tushev. All rights reserved.
//

#ifndef __fastqSeqStats__utilities__
#define __fastqSeqStats__utilities__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// define some constants
#define PHRED_MAX 255
#define MAX_LENGTH 1024

// define MIN & MAX macros
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, MIN(actualLen, paramLen))== 0) && (actualLen == paramLen)

// define C BOOLEAN type
typedef enum{FALSE, TRUE} boolean;
typedef char DNA;
typedef char PHRED;

// define Parameters Type
typedef struct
{
    FILE *source;
    FILE *destination;
    int seqlength;
    boolean polytail;
} Parameters;

// function declarations
int chompNewLine(char *);
int parseParameters(int, char **, Parameters *);
int help(void);
int version(void);

#endif /* defined(__fastqSeqStats__utilities__) */