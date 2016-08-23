#-----------------------------------------------------------------------------
# File     : makefile
# Contents : build fcmat test program
# Author   : Kristian Loewe
#-----------------------------------------------------------------------------
SHELL      = /bin/bash
CPUINFODIR = ../../cpuinfo/src
CORRDIR    = ../../corr/src
STATSDIR    = ../../stats/src

CC         = gcc -std=c99 -march=native
CFBASE     = -Wall -Wextra -Wno-unused-parameter -Wconversion -Wshadow \
             -pedantic $(ADDFLAGS)
CFLAGS     = $(CFBASE) -DNDEBUG -O3
#CFLAGS     = $(CFBASE) -O3
#CFLAGS     = $(CFBASE) -O1 -g

REAL        ?= float
DREAL       =
ifneq "$(strip $(REAL))" ""
  DREAL       = -DREAL=$(REAL)
endif
DEFS        = $(DREAL)

CFLAGS     += $(DEFS)

INCS       = -I$(STATSDIR) -I$(CPUINFODIR) -I$(CORRDIR)

LD         = gcc
LDFLAGS    =
LIBS       = -lm -lrt -lpthread

OUTDIR     = ../bin
HDRS       = $(CPUINFODIR)/cpuinfo.h \
             $(CORRDIR)/clamp.h \
             $(CORRDIR)/binarize.h \
             $(STATSDIR)/stats.h \
             $(CORRDIR)/pcc.h \
             $(CORRDIR)/tetracc.h
OBJS       = $(CPUINFODIR)/cpuinfo.o \
             $(CORRDIR)/binarize.o \
             $(STATSDIR)/stats.o \
             $(CORRDIR)/pcc.o \
             $(CORRDIR)/tetracc.o
PRGS       = test_fcmat test_fcmat1 test_fcmat2 test_fcmat3

#-----------------------------------------------------------------------------
# Build Programs
#-----------------------------------------------------------------------------
all: $(PRGS)

$(PRGS): | $(OUTDIR)

$(OUTDIR):
	mkdir -p $(OUTDIR)

# all-in-one
test_fcmat:  ../bin/test_fcmat
	

# on-demand
test_fcmat1: ../bin/test_fcmat1
	

# cache-based
test_fcmat2: ../bin/test_fcmat2
	

# half-stored
test_fcmat3: ../bin/test_fcmat3
	

# all-in-one
../bin/test_fcmat:  $(OBJS) fcmat.o test_fcmat.o makefile
	$(LD) $(LDFLAGS) $(OBJS) fcmat.o test_fcmat.o $(LIBS) -o $@

# on-demand
../bin/test_fcmat1: $(OBJS) fcmat1.o test_fcmat.o makefile
	$(LD) $(LDFLAGS) $(OBJS) fcmat1.o test_fcmat.o $(LIBS) -o $@

# cache-based
../bin/test_fcmat2: $(OBJS) fcmat2.o test_fcmat.o makefile
	$(LD) $(LDFLAGS) $(OBJS) fcmat2.o test_fcmat.o $(LIBS) -o $@

# half-stored
../bin/test_fcmat3: $(OBJS) fcmat3.o test_fcmat.o makefile
	$(LD) $(LDFLAGS) $(OBJS) fcmat3.o test_fcmat.o $(LIBS) -o $@

#-----------------------------------------------------------------------------
# Test Program
#-----------------------------------------------------------------------------
test_fcmat.o:  fcmat.h $(HDRS)
test_fcmat.o:  test_fcmat.c makefile
	$(CC) $(CFLAGS) $(INCS) -c test_fcmat.c -o $@

#-----------------------------------------------------------------------------
# Modules
#-----------------------------------------------------------------------------
fcmat.o:  fcmat.h fcmat1.h fcmat3.h $(HDRS)
fcmat.o:  fcmat.c makefile
	$(CC) $(CFLAGS) $(INCS) -c fcmat.c -o $@

fcmat1.o:     fcmat.h fcmat1.h $(HDRS)
fcmat1.o:     fcmat1.c makefile
	$(CC) $(CFLAGS) $(INCS) -c fcmat1.c -o $@

fcmat2.o:     fcmat.h $(HDRS)
fcmat2.o:     fcmat2.c makefile
	$(CC) $(CFLAGS) $(INCS) -c fcmat2.c -o $@

fcmat3.o:     fcmat.h fcmat3.h $(HDRS)
fcmat3.o:     fcmat3.c makefile
	$(CC) $(CFLAGS) $(INCS) -c fcmat3.c -o $@

#-----------------------------------------------------------------------------
# External Modules
#-----------------------------------------------------------------------------
$(CPUINFODIR)/cpuinfo.o:
	cd $(CPUINFODIR); $(MAKE) cpuinfo.o

$(CORRDIR)/pcc.o:
	cd $(CORRDIR);    $(MAKE) pcc.o      ADDFLAGS="$(DREAL)"

$(CORRDIR)/pccbnch.o:
	cd $(CORRDIR);    $(MAKE) pccbnch.o  ADDFLAGS="$(DREAL)"

$(CORRDIR)/tetracc.o:
	cd $(CORRDIR);    $(MAKE) tetracc.o  ADDFLAGS="$(DREAL)"

$(CORRDIR)/tccbnch.o:
	cd $(CORRDIR);    $(MAKE) tccbnch.o  ADDFLAGS="$(DREAL)"

$(CORRDIR)/binarize.o:
	cd $(CORRDIR);    $(MAKE) binarize.o ADDFLAGS="$(DREAL)"

$(STATSDIR)/stats.o:
	cd $(STATSDIR);   $(MAKE) stats.o    ADDFLAGS="$(DREAL)"
