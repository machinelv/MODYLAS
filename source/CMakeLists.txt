cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(MODYLAS NONE)
enable_language(Fortran)

option(GROEXT "Use Gromacs Extension" OFF)
if(GROEXT)
  enable_language(CXX)
  message(STATUS "Gromacs Extension ON : CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
  add_definitions(-DGROEXT)
endif(GROEXT)

option(FJ_RDMA      "Use Fujitsu RDMA"      OFF)
option(PROFILE_COMM "Profile Communication" OFF)
if(FJ_RDMA)
  enable_language(C)
  set(CMAKE_C_COMPILER "mpifccpx")
  message(STATUS "RDMA ON for Fujitsu system (K, FX100): CMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
  add_definitions(-DFJ_RDMA)
endif(FJ_RDMA)

if(PROFILE_COMM)
  message(STATUS "PROFILE Communication ON")
  add_definitions(-DPROFILE_COMM)
endif(PROFILE_COMM)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

option(DISABLE_OPENMP "Disable OpenMP Parallelization" OFF)
option(DISABLE_MPI "Disable MPI Parallelization" OFF)

# OpenMP
if(NOT DISABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
#    set(MODYLAS_USE_MKL_PARALLEL TRUE)
    if(GROEXT)
      if(${CMAKE_Fortran_COMPILER} MATCHES "frtpx")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      endif(${CMAKE_Fortran_COMPILER} MATCHES "frtpx")
    endif(GROEXT)
  else(OPENMP_FOUND)
    set(DISABLE_OPENMP ON)
  endif(OPENMP_FOUND)
endif(NOT DISABLE_OPENMP)
if(NOT DISABLE_OPENMP)
  set(MODYLAS_HAVE_OPENMP True)
else(NOT DISABLE_OPENMP)
  set(MODYLAS_HAVE_OPENMP False)
endif(NOT DISABLE_OPENMP)

# MPI library
if(NOT DISABLE_MPI)
  find_package(MPI)
  if(MPI_Fortran_FOUND)
    set(CMAKE_EXE_LINKER_FLAGS ${MPI_CXX_LINK_FLAGS})
  else(MPI_Fortran_FOUND)
    set(DISABLE_MPI ON)
  endif(MPI_Fortran_FOUND)
endif(NOT DISABLE_MPI)
if(NOT DISABLE_MPI)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DMPIPARA")
endif(NOT DISABLE_MPI)

# standard include directories
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
include_directories(${MPI_CXX_INCLUDE_PATH})

# timer
option(USE_MAPROF "Using MiniApp Profiler" OFF)
message(STATUS "Using MiniApp Profiler: USE_MAPROF=${USE_MAPROF}")

# inline option for Fujitsu compiler
option(INLINE_FUJITSU "Adding inline options for Fujitsu compiler" OFF)
message(STATUS "Adding inline options for Fujitsu compiler: INLINE_FUJITSU=${INLINE_FUJITSU}")

# Extracting #ifdef variables in CMAKE_Fortran_FLAGS
string(REPLACE " " ";" LIST_CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS}) # splitting string by whitespace into strings
foreach(name ${LIST_CMAKE_Fortran_FLAGS})
  string(REGEX MATCH "-D.*" TMP_BOOL_D ${name})
  if(TMP_BOOL_D)
    string(REPLACE "-D" "" PREPROCESSOR_DEF ${name})  # Trimming "-D"
    list(APPEND LIST_PREPROCESSOR_DEF ${PREPROCESSOR_DEF})
  endif(TMP_BOOL_D)
endforeach(name)

add_subdirectory(src)
##add_subdirectory(test)
##add_subdirectory(benchmark)
