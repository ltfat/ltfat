cmake_minimum_required (VERSION 2.8)
project (libltfat C CXX)

option(NOBLASLAPACK
    "Disable BLAS and LAPACK dependency" OFF)

if (MSVC)
    set(USECPP 1)
else (MSVC)
    option(USECPP
        "Use C++ compiler" OFF)
endif (MSVC)

FIND_LIBRARY(FFTW3_LIB NAMES fftw3 libfftw3)
FIND_LIBRARY(FFTW3F_LIB NAMES fftw3f libfftw3f)

if (NOT NOBLASLAPACK)
     FIND_LIBRARY(BLAS_LIB NAMES blas libblas)
     FIND_LIBRARY(LAPACK_LIB NAMES lapack liblapack)
endif (NOT NOBLASLAPACK)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)

add_subdirectory(src)

