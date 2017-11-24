#F90    = ifort
#FFLAGS =  -convert little_endian -assume byterecl 
F90         = gfortran
FFLAGS      = -cpp -frecord-marker=4 #-O3
NETCDF_LIB  = $(shell nc-config --flibs)
NETCDF_INC  = $(shell nc-config --fflags)

INCS = ${NETCDF_INC}
LIBS = ${NETCDF_LIB} 

.SUFFIXES: .o .c .f90 # .f
.PHONY: all debug clean

PROJNAME = pmctrack
TARGET = track.out
OUTDIR = output
SRCDIR = $(PROJNAME)/src
OBJDIR = $(PROJNAME)/src/_precc
OBJ = \
$(OBJDIR)/datetime.o \
$(OBJDIR)/types.o \
$(OBJDIR)/nc_io.o \
$(OBJDIR)/params.o \
$(OBJDIR)/constants.o \
$(OBJDIR)/utils.o \
$(OBJDIR)/vor_partition.o \
$(OBJDIR)/cf_synop_check.o \
$(OBJDIR)/synop_check.o \
$(OBJDIR)/min_z.o \
$(OBJDIR)/steering_wind.o \
$(OBJDIR)/linkin_vort2.o \
$(OBJDIR)/linkin_vort.o \
$(OBJDIR)/track_check.o \
$(OBJDIR)/smth.o \
$(OBJDIR)/tracking_main.o \
$(OBJDIR)/main.o

all: $(TARGET)

run: $(TARGET)
	./$(TARGET)

debug: FFLAGS += -g -fcheck=all -fbacktrace -Ddebug -Wall
debug: $(TARGET)

help:
	@echo 'Makefile for the tracking code                '
	@echo '                                              '
	@echo 'Usage:                                        '
	@echo '    make all     	Compile tracking code'
	@echo '    make clean	        Clean the directory  '
	@echo '                                              '

$(TARGET): $(OBJ)
	@mkdir -p $(OUTDIR)
	$(F90) $(FFLAGS) -o $@ $(OBJ) ${INCS} ${LIBS}

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@mkdir -p $(OBJDIR)
	$(F90) $(FFLAGS) ${INCS} -c $< -o $@ 

#.f.o:
#	$(F90) $(FFLAGS)  -c $<


clean:
	-rm -f $(OBJDIR)/*[.o,.mod]
	-rm -f $(TARGET)
	-rm -f $(OUTDIR)/vor*[.txt,.dat]
