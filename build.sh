#!/usr/bin/env bash
PLATFORM=${1:-ubuntu}
MODE=${2:-quick}
export PROJNAME=pmctrack
export TARGET=track.x
export OUTDIR=output
export MODDIR=${PROJNAME}/modules
export SRCDIR=${PROJNAME}/src

export OBJDIR=${PROJNAME}/src/_precc
export OBJ="\
${OBJDIR}/datetime.o \
${OBJDIR}/types.o \
${OBJDIR}/nc_io.o \
${OBJDIR}/constants.o \
${OBJDIR}/params.o \
${OBJDIR}/utils.o \
${OBJDIR}/vor_partition.o \
${OBJDIR}/cf_synop_check.o \
${OBJDIR}/synop_check.o \
${OBJDIR}/min_z.o \
${OBJDIR}/steering_wind.o \
${OBJDIR}/link_vort_rad.o \
${OBJDIR}/smth.o \
${OBJDIR}/main.o"
#${OBJDIR}/check_track.o \
#${OBJDIR}/linkin_vort.o \
#${OBJDIR}/linkin_vort2.o \
# ${OBJDIR}/tracking_main.o \

# if { $#argv == 0 } then
#   echo "install.sh [lib|core|goodies|links|all|docu|clean] "
#   exit 0
# endif
ACTION=BUILD
if [ $PLATFORM == "ubuntu" ]; then
    export FORTRAN=gfortran
    NETCDF_INC="-I/usr/include"
    NETCDF_LIB="-L/usr/lib -lnetcdff -lnetcdf"
elif [ $PLATFORM == "jasmin" ]; then
    module load intel/14.0
    module load netcdff/intel/14.0/4.2 
    export FORTRAN=ifort
    NETCDF_INC="-I/apps/libs/netCDF/intel14/fortran/4.2/include"
    NETCDF_LIB="-L/apps/libs/netCDF/intel14/fortran/4.2/lib -lnetcdff"
elif [ $PLATFORM == "archer" ]; then
    module swap PrgEnv-cray PrgEnv-intel
    module load intel/16.0.2.181 # otherwise libifcore.so.5 is not found
    module load cray-netcdf
    export FORTRAN=ifort
    netcdf_prefix=/opt/cray/netcdf/4.4.1.1/INTEL/15.0
    NETCDF_INC="-I${netcdf_prefix}/include"
    NETCDF_LIB="-L${netcdf_prefix}/lib -lnetcdff -lnetcdf"
elif [ $PLATFORM == "clean" ]; then
    ACTION=OTHER
    make -f Makefile clean
else
    echo "build.sh [ubuntu|jasmin|archer|clean]"
    exit 0
fi

if [ $ACTION == "BUILD" ]; then
    if [ $FORTRAN == "gfortran" ]; then
        FFLAGS="-cpp -frecord-marker=4 -Dgnu"
        if [ $MODE == "quick" ]; then
            FFLAGS+=" -O3"
        elif [ $MODE == "debug" ]; then
            FFLAGS+=" -O0 -g -fcheck=all -fbacktrace -Ddebug -Wall"
        fi
        FFLAGS+=" -J ${MODDIR}"
    elif [ $FORTRAN == "ifort" ]; then
        FFLAGS="-convert little_endian -assume byterecl"
        if [ $MODE == "quick" ]; then
            FFLAGS+=" -O3"
        elif [ $MODE == "debug" ]; then
            FFLAGS+=" -O0 -g -fcheck=all -fbacktrace -Ddebug -Wall"
        fi
        FFLAGS+=" -I${MODDIR}"
    else
        echo "[gfortran|ifort]"
        exit 0
    fi
    export INCS=${NETCDF_INC}
    export LIBS=${NETCDF_LIB}
    export FFLAGS
    
    make -f Makefile all
fi
