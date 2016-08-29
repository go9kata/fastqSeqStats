# fastqSeqStats
calculate base frequency and quality distribution from FASTQ file

### Install
C-compiler is required.

```
git clone https://github.molgen.mpg.de/MPIBR-Bioinformatics/fastqSeqStats.git
cd fastqSeqStats
make
```

### Summary
fastqSeqStats creaes histogram of base quality score. Base composition and quality score per nucleotide are accumulated throughout sequence length. Creates histogram of loose count for poly(A) tail or poly(T) head. 

### Usage
`fastqSeqStats [OPTIONS] --ifastq <fastq/stdin> --otxt <txt/stdout>`

### Options
--polyTail :: create histogram of loose poly Tail/Head count based on Kent et.al. 2002 algorithm
-h/--help :: print help message
-v/--version :: print current version
