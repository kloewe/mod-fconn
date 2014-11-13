#------------------------------------------------------------------------------
# File     : makefile
# Contents : build fcmat test program
# Author   : Kristian Loewe
#------------------------------------------------------------------------------
SHELL      = /bin/bash
CPUINFODIR = ../../cpuinfo/src
CORRDIR    = ../../corr/src

CC         = gcc -std=c99 -march=native
CFBASE     = -Wall -Wextra -Wno-unused-parameter -Wconversion -pedantic
INCS       = -I$(CPUINFODIR) -I$(CORRDIR)
CFLAGS     = $(CFBASE) -DNDEBUG -O3 -DREAL=float $(INCS)

LD         = gcc
LDFLAGS    =
LIBS       = -lm -lpthread

OBJS       = $(CPUINFODIR)/cpuinfo.o \
             $(CORRDIR)/pcc.o \
             $(CORRDIR)/tetracc.o \
             $(CORRDIR)/binarize.o
PRGS       = fcmat

#------------------------------------------------------------------------------
# Build Program
#------------------------------------------------------------------------------
all:    $(PRGS)

fcmat:  $(OBJS) fcmat.o makefile
	$(LD) $(LDFLAGS) $(OBJS) fcmat.o $(LIBS) -o $@

fcmat.o: $(OBJS)
	$(CC) $(CFLAGS) -DFCMAT_MAIN -c fcmat.c -o $@
