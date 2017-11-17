#F90    = ifort
#FFLAGS =  -convert little_endian -assume byterecl 
F90    = gfortran
FFLAGS = -O3 -frecord-marker=4 -cpp 

.SUFFIXES: .o .c .f90 # .f
.PHONY: all debug clean

PROJNAME = pmctrack
TARGET = track.out
OUTDIR = output
# QUICKVIEW = tools/quickview.py
SRCDIR = $(PROJNAME)/src
OBJDIR = $(PROJNAME)/src/_precc
OBJ = \
$(OBJDIR)/types.o \
$(OBJDIR)/const.o \
$(OBJDIR)/params.o \
$(OBJDIR)/util.o \
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
$(OBJDIR)/interface.o

all : $(TARGET)

run : $(TARGET)
	./$(TARGET)

debug : FFLAGS += -g -fcheck=all -fbacktrace -Ddebug # -Wall
debug : $(TARGET)

help :
	@echo 'Makefile for the tracking code                                 '
	@echo '                                                               '
	@echo 'Usage:                                                         '
	@echo '    make all                              Compile tracking code'
	@echo '    make debug                     Compile with debugging flags'
	@echo '    make run                         Run the tracking algorithm'
	@echo '    make clean                Clean the output and object files'
	@echo '                                                               '

$(TARGET) : $(OBJ)
	@mkdir -p $(OUTDIR)
	$(F90) $(FFLAGS) -o $@ $(OBJ)

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	@mkdir -p $(OBJDIR)
	$(F90) $(FFLAGS) -c $< -o $@ 

clean :
	-rm -f $(OBJDIR)/*[.o,.mod]
	-rm -f $(TARGET)
	-rm -f $(OUTDIR)/vor*[.txt,.dat]
