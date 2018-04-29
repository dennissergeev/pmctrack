#!/usr/bin/env bash
PLATFORM=${1:-linux}
BUILD_TYPE=${2:-Release}

if [ $PLATFORM == "linux" ]; then
    FORTRAN=gfortran
    INCS=/usr/include
    LIBS="-L/usr/lib -lnetcdff -lnetcdf"
elif [ $PLATFORM == "monsoon" ]; then
    module swap PrgEnv-cray PrgEnv-gnu
    # module load gcc/4.8.1
    module load cray-netcdf
    netcdf_prefix="/opt/cray/netcdf/4.3.2/GNU/49"
    hdf_prefix="/opt/cray/hdf5/1.8.13/GNU/49"
    FORTRAN=gfortran
    INCS=${netcdf_prefix}/include
    LIBS="-L${netcdf_prefix}/lib -lnetcdff -L${hdf_prefix}/lib -lnetcdf"
elif [ $PLATFORM == "jasmin" ]; then
    module load intel/14.0
    module load netcdff/intel/14.0/4.2 
    FORTRAN=ifort
    INCS=/apps/libs/netCDF/intel14/fortran/4.2/include
    LIBS="-L/apps/libs/netCDF/intel14/fortran/4.2/lib -lnetcdff"
elif [ $PLATFORM == "archer" ]; then
    module swap PrgEnv-cray PrgEnv-intel
    module load intel/16.0.2.181 # otherwise libifcore.so.5 is not found
    module load cray-netcdf
    FORTRAN=ifort
    netcdf_prefix=/opt/cray/netcdf/4.4.1.1/INTEL/15.0
    INCS=${netcdf_prefix}/include
    LIBS="-L${netcdf_prefix}/lib -lnetcdff -lnetcdf"
else
    echo "No options specified for '${PLATFORM}' platform"
    exit
fi

mkdir -p build
cd build
cmake -DFC=${FORTRAN} -DINCS=${INCS} -DLIBS="${LIBS}" -DCMAKE_BUILD_TYPE=${BUILD_TYPE} ..
make
cd ..
