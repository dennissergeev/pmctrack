#F90    = ifort
#FFLAGS =  -convert little_endian -assume byterecl 
F90    = gfortran
FFLAGS = -O3 -frecord-marker=4 

.SUFFIXES: .o .c .f90 # .f
.PHONY: all debug clean

PROJNAME = pmctrack
TARGET = track.exe
OUTDIR = output
SRCDIR = $(PROJNAME)/src
OBJDIR = $(PROJNAME)/src/_precc
OBJ = \
$(OBJDIR)/const.o \
$(OBJDIR)/params.o \
$(OBJDIR)/vor_partition.o \
$(OBJDIR)/cf_synop_check.o \
$(OBJDIR)/synop_check.o \
$(OBJDIR)/min_z.o \
$(OBJDIR)/steering_wind.o \
$(OBJDIR)/link2.o \
$(OBJDIR)/link.o \
$(OBJDIR)/track_check.o \
$(OBJDIR)/smth.o \
$(OBJDIR)/tracking_main.o \
$(OBJDIR)/interface.o

all: $(TARGET)

debug: FFLAGS += -fcheck=all -fbacktrace -g
debug: $(TARGET)

help:
	@echo 'Makefile for the tracking code                '
	@echo '                                              '
	@echo 'Usage:                                        '
	@echo '    make all     	Compile tracking code'
	@echo '    make clean	        Clean the directory  '
	@echo '                                              '

$(TARGET) : $(OBJ)
	@mkdir -p $(OUTDIR)
	$(F90) $(FFLAGS) -o $@ $(OBJ)

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	@mkdir -p $(OBJDIR)
	$(F90) $(FFLAGS) -c $< -o $@ 

#.f.o :
#	$(F90) $(FFLAGS)  -c $<


clean :
	-rm -f $(OBJDIR)/*[.o,.mod]
	-rm -f $(TARGET)
	-rm -f $(OUTDIR)/vor*[.txt,.dat]
