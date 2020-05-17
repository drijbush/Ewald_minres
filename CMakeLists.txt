# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(Ewald_minres Fortran)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${Ewald_minres_BINARY_DIR})
find_package(MPI REQUIRED)


# Add OMP
#find_package(OpenMP)
if (OpenMP_Fortran_FOUND)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()
find_library(fftw3 REQUIRED)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O -g -std=f2008 -Wall -Wextra -Wimplicit-interface -Wuse-without-only -fcheck=all -finit-real=snan" CACHE STRING "Flags used in debug compiles" FORCE)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops -std=f2008" CACHE STRING "Flags used in release compiles" FORCE)
endif()


add_subdirectory("src")
