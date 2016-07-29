#-----------------------------------------------------------------------------
# File     : makefile
# Contents : build fcmat test program
# Author   : Kristian Loewe
#-----------------------------------------------------------------------------
SHELL      = /bin/bash
CPUINFODIR = ../../cpuinfo/src
CORRDIR    = ../../corr/src

CC         = gcc -std=c99 -march=native
CFBASE     = -Wall -Wextra -Wno-unused-parameter -Wconversion -Wshadow \
             -pedantic $(ADDFLAGS)
CFLAGS     = $(CFBASE) -DNDEBUG -O3
#CFLAGS     = $(CFBASE) -O1 -g

REAL        ?= float
DREAL       =
ifneq "$(strip $(REAL))" ""
  DREAL       = -DREAL=$(REAL)
endif
DEFS        = $(DREAL)

CFLAGS     += $(DEFS)

INCS       = -I$(CPUINFODIR) -I$(CORRDIR)

LD         = gcc
LDFLAGS    =
LIBS       = -lm -lrt -lpthread

OUTDIR     = ../bin
HDRS       = $(CPUINFODIR)/cpuinfo.h \
             $(CORRDIR)/clamp.h \
             $(CORRDIR)/binarize.h \
             $(CORRDIR)/pcc.h \
             $(CORRDIR)/tetracc.h
OBJS       = $(CPUINFODIR)/cpuinfo.o \
             $(CORRDIR)/binarize.o \
             $(CORRDIR)/pcc.o \
             $(CORRDIR)/tetracc.o
PRGS       = $(OUTDIR)/test_fcmat

#-----------------------------------------------------------------------------
# Build Program
#-----------------------------------------------------------------------------
all: $(PRGS)

$(PRGS): | $(OUTDIR)

$(OUTDIR):
	mkdir -p $(OUTDIR)

test_fcmat:  ../bin/test_fcmat
	

../bin/test_fcmat:  $(OBJS) fcmat.o test_fcmat.o makefile
	$(LD) $(LDFLAGS) $(OBJS) fcmat.o test_fcmat.o $(LIBS) -o $@

#-----------------------------------------------------------------------------
# Test Programs
#-----------------------------------------------------------------------------
test_fcmat.o:  fcmat.h $(HDRS)
test_fcmat.o:  test_fcmat.c makefile
	$(CC) $(CFLAGS) $(INCS) -c test_fcmat.c -o $@

#-----------------------------------------------------------------------------
# Main Modules
#-----------------------------------------------------------------------------
fcmat.o:  fcmat.h $(HDRS)
fcmat.o:  fcmat.c makefile
	$(CC) $(CFLAGS) $(INCS) -c fcmat.c -o $@

#-----------------------------------------------------------------------------
# External Modules
#-----------------------------------------------------------------------------
$(CPUINFODIR)/cpuinfo.o:
	cd $(CPUINFODIR); $(MAKE) cpuinfo.o

$(CORRDIR)/pcc.o:
	cd $(CORRDIR); $(MAKE) pcc.o      ADDFLAGS="$(DREAL)"

$(CORRDIR)/pccbnch.o:
	cd $(CORRDIR); $(MAKE) pccbnch.o  ADDFLAGS="$(DREAL)"

$(CORRDIR)/tetracc.o:
	cd $(CORRDIR); $(MAKE) tetracc.o  ADDFLAGS="$(DREAL)"

$(CORRDIR)/tccbnch.o:
	cd $(CORRDIR); $(MAKE) tccbnch.o  ADDFLAGS="$(DREAL)"

$(CORRDIR)/binarize.o:
	cd $(CORRDIR); $(MAKE) binarize.o ADDFLAGS="$(DREAL)"

