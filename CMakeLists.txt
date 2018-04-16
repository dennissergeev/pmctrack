cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)
project (PMCTRACK)
enable_language(Fortran)

set(EXE track.x)
set(SRCDIR src)

if(NOT PLATFORM)
    set(PLATFORM linux)
endif()

if(PLATFORM MATCHES "linux")
    set(CMAKE_Fortran_COMPILER gfortran)
    set(INCS /usr/include)
    set(LIBS "-L/usr/lib -lnetcdff -lnetcdf")
elseif(PLATFORM MATCHES "monsoon")
    set(CMAKE_Fortran_COMPILER gfortran)
    set(MODULE_CMDS
    "module swap PrgEnv-cray PrgEnv-gnu"
    "module load cray-netcdf"
    )
    set(netcdf_prefix "/opt/cray/netcdf/4.3.2/GNU/49")
    set(hdf_prefix "/opt/cray/hdf5/1.8.13/GNU/49")
    set(INCS ${netcdf_prefix}/include)
    set(LIBS "-L${netcdf_prefix}/lib -lnetcdff -L${hdf_prefix}/lib -lnetcdf")
elseif(PLATFORM MATCHES "archer")
    set(CMAKE_Fortran_COMPILER ifort)
    set(MODULE_CMDS
    "module swap PrgEnv-cray PrgEnv-intel"
    "module load intel/16.0.2.181" # otherwise libifcore.so.5 is not found
    "module load cray-netcdf"
    )
    set(netcdf_prefix "/opt/cray/netcdf/4.4.1.1/INTEL/15.0")
    set(INCS ${netcdf_prefix}/include)
    set(LIBS "-L${netcdf_prefix}/lib -lnetcdff -lnetcdf")
else()
    MESSAGE(FATAL_ERROR "Uknown option PLATFORM=${PLATFORM}")
endif()

# this is where we will place the Fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

function(hpc_env module_cmd)
    add_custom_target(${module_cmd} ALL
                      COMMAND echo ${module_cmd}
                      COMMENT
                      VERBATIM)
endfunction()

foreach(module_cmd ${MODULE_CMDS})
    hpc_env(${module_cmd})
endforeach()

# we default to Release build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    message(FATAL_ERROR "Fortran compiler does not support F90")
endif(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS         "-cpp -frecord-marker=4 -Dgnu")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fcheck=all -fbacktrace -Ddebug -Wall")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES INTEL)
    set(CMAKE_Fortran_FLAGS         "-convert little_endian -assume byterecl")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fcheck=all -fbacktrace -Ddebug -Wall")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
else()
    MESSAGE(FATAL_ERROR "Unknown compiler")
endif()

include_directories(${INCS})
add_executable(
    ${EXE}
    ${SRCDIR}/datetime.f90
    ${SRCDIR}/types.f90
    ${SRCDIR}/nc_io.f90
    ${SRCDIR}/constants.f90
    ${SRCDIR}/params.f90
    ${SRCDIR}/utils.f90
    ${SRCDIR}/vor_partition.f90
    ${SRCDIR}/cf_synop_check.f90
    ${SRCDIR}/synop_check.f90
    ${SRCDIR}/min_z.f90
    ${SRCDIR}/steering_wind.f90
    ${SRCDIR}/link_vort_rad.f90
    ${SRCDIR}/smth.f90
    ${SRCDIR}/main.f90
)
target_link_libraries(${EXE} ${LIBS})
