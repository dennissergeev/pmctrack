cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)
project (PMCTRACK)
enable_language(Fortran)

set(EXE track.x)

if(PLATFORM MATCHES "ubuntu")
    set(CMAKE_Fortran_COMPILER gfortran)
    set(INCS "/usr/include")
    set(LIBS "/usr/lib -lnetcdff -lnetcdf")
else()
    set(MODULE_CMDS
    "module load foo"
    "module load bar"
    )
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


MESSAGE("${CMAKE_Fortran_COMPILER}")

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
link_directories (${LIBS})
add_executable(
    ${EXE}
    src/datetime.f90
    src/types.f90
    src/nc_io.f90
    src/constants.f90
    src/params.f90
    src/utils.f90
    src/vor_partition.f90
    src/cf_synop_check.f90
    src/synop_check.f90
    src/min_z.f90
    src/steering_wind.f90
    src/link_vort_rad.f90
    src/smth.f90
    src/main.f90
)
