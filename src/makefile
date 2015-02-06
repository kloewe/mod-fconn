#-----------------------------------------------------------------------
# File     : makefile
# Contents : build fcmat test program
# Author   : Kristian Loewe
#-----------------------------------------------------------------------
SHELL      = /bin/bash
CPUINFODIR = ../../cpuinfo/src
CORRDIR    = ../../corr/src

CC         = gcc -std=c99 -march=native
CFBASE     = -Wall -Wextra -Wno-unused-parameter -Wconversion \
             -pedantic $(ADDFLAGS)
CFLAGS     = $(CFBASE) -DNDEBUG -O3
INCS       = -I$(CPUINFODIR) -I$(CORRDIR)

LD         = gcc
LDFLAGS    =
LIBS       = -lm -lrt -lpthread

HDRS       = $(CPUINFODIR)/cpuinfo.h \
             $(CORRDIR)/clamp.h \
             $(CORRDIR)/pcc.h \
             $(CORRDIR)/binarize.h \
             $(CORRDIR)/tetracc.h
OBJS       = $(CPUINFODIR)/cpuinfo.o \
             $(CORRDIR)/pcc.o \
             $(CORRDIR)/binarize.o \
             $(CORRDIR)/tetracc.o
PRGS       = ../bin/fcmat

#-----------------------------------------------------------------------
# Build Program
#-----------------------------------------------------------------------
all:    $(PRGS)

fcmat:      ../bin/fcmat
	echo ""

../bin/fcmat:  $(OBJS) fcm_main.o makefile
	$(LD) $(LDFLAGS) $(OBJS) fcm_main.o $(LIBS) -o $@

#-----------------------------------------------------------------------
# Main Programs
#-----------------------------------------------------------------------
fcm_main.o:   fcmat.h $(HDRS)
fcm_main.o:   fcmat.c makefile
	$(CC) $(CFLAGS) -funroll-loops $(INCS) -DFCMAT_MAIN \
              -c fcmat.c -o $@

#-----------------------------------------------------------------------
# Main Modules
#-----------------------------------------------------------------------
fcmat.o:    fcmat.h $(HDRS)
fcmat.o:    fcmat.c makefile
	$(CC) $(CFLAGS) -funroll-loops $(INCS) -c fcmat.c -o $@

#-----------------------------------------------------------------------
# External Modules
#-----------------------------------------------------------------------
$(CPUINFODIR)/cpuinfo.o:
	cd $(CPUINFODIR); $(MAKE) cpuinfo.o

$(CORRDIR)/pcc.o:
	cd $(CORRDIR); $(MAKE) pcc.o

$(CORRDIR)/tetracc.o:
	cd $(CORRDIR); $(MAKE) tetracc.o

$(CORRDIR)/binarize.o:
	cd $(CORRDIR); $(MAKE) binarize.o
