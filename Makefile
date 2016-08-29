# define the C compiler to use
CC=gcc

# define any compile-time flags
CFLAGS=-Wall -g -std=c99
CLIBS=-lm

# define the C source files
SRCS=fastqSeqStatsMain.c utilities.c fastqCounter.c dnaPolyTail.c

# define the C object files
OBJS=$(SRCS:.c=.o)

# define the executable file
MAIN=fastqSeqStats

all:	$(MAIN)

$(MAIN):	$(OBJS)
	$(CC) $(CFLAGS) $(CLIBS) -o $(MAIN) $(OBJS)

.c.o:
	$(CC) $(CFLAGS) $(CLIBS) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o $(MAIN)
