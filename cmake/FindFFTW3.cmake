# FindFFTW3.cmake - Find FFTW3 library
#
# This module defines:
#  FFTW3_FOUND - True if FFTW3 is found
#  FFTW3_INCLUDE_DIRS - The FFTW3 include directories
#  FFTW3_LIBRARIES - The FFTW3 libraries for linking

# Look for the header file
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h
          HINTS ENV FFTW3_DIR ENV FFTW3_ROOT
          PATH_SUFFIXES include)

# Look for the library
find_library(FFTW3_LIBRARY NAMES fftw3 libfftw3
             HINTS ENV FFTW3_DIR ENV FFTW3_ROOT
             PATH_SUFFIXES lib lib64)

# Look for the OpenMP library (optional)
find_library(FFTW3_OMP_LIBRARY NAMES fftw3_omp libfftw3_omp
             HINTS ENV FFTW3_DIR ENV FFTW3_ROOT
             PATH_SUFFIXES lib lib64)

# Look for the MPI library (optional)
find_library(FFTW3_MPI_LIBRARY NAMES fftw3_mpi libfftw3_mpi
             HINTS ENV FFTW3_DIR ENV FFTW3_ROOT
             PATH_SUFFIXES lib lib64)

# Try to determine the version
if(FFTW3_INCLUDE_DIR)
  file(READ "${FFTW3_INCLUDE_DIR}/fftw3.h" _fftw3_version_header)
  string(REGEX MATCH "#define FFTW_VERSION[ \t]+\"([^\"]+)\"" _fftw3_version_match "${_fftw3_version_header}")
  if(_fftw3_version_match)
    set(FFTW3_VERSION "${CMAKE_MATCH_1}")
  endif()
endif()

# Set variables for standard handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3 
                                  REQUIRED_VARS FFTW3_LIBRARY FFTW3_INCLUDE_DIR
                                  VERSION_VAR FFTW3_VERSION)

if(FFTW3_FOUND)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
  
  # Add OpenMP library if found
  if(FFTW3_OMP_LIBRARY)
    list(APPEND FFTW3_LIBRARIES ${FFTW3_OMP_LIBRARY})
  endif()
  
  # Add MPI library if found
  if(FFTW3_MPI_LIBRARY)
    list(APPEND FFTW3_LIBRARIES ${FFTW3_MPI_LIBRARY})
  endif()
  
  set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
endif()

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY FFTW3_OMP_LIBRARY FFTW3_MPI_LIBRARY)
