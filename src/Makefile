#F90    = ifort
#FFLAGS =  -convert little_endian -assume byterecl 
F90    = gfortran
FFLAGS = -O3 -frecord-marker=4 


.SUFFIXES: .o .c .f .f90

TARGET = tracking
OBJ =     const.o    vor_partition.o cf_synop_check.o synop_check.o min_z.o steering_wind.o link2.o link.o  track_check.o  smth.o tracking_main.o tracking_sample.o

$(TARGET) : $(OBJ)
	$(F90) $(FFLAGS) -o $@ $(OBJ) 

.f90.o :
	$(F90) $(FFLAGS) -c $< 

.f.o :
	$(F90) $(FFLAGS)  -c $<


clean :
	-rm -f *.o *.mod $(TARGET)
