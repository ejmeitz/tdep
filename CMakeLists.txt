cmake_minimum_required(VERSION 3.12)

# Project declaration
project(TDEP 
        VERSION 1.0.0
        LANGUAGES Fortran C CXX)

# RPATH settings for shared libraries
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Get Git information
find_package(Git QUIET)
if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_REVISION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()

# Options
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(USE_CGAL "Enable CGAL support" OFF)
option(ENABLE_OPENMP "Enable OpenMP support" ON)
option(GENERATE_MANPAGES "Generate manpages" ON)

# Set position independent code for all targets when building shared libraries
if(BUILD_SHARED_LIBS)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Create directories
set(TDEP_BIN_DIR ${CMAKE_BINARY_DIR}/bin)
set(TDEP_LIB_DIR ${CMAKE_BINARY_DIR}/lib)
set(TDEP_INC_DIR ${CMAKE_BINARY_DIR}/inc)
set(TDEP_MAN_DIR ${CMAKE_BINARY_DIR}/man)
set(TDEP_DOC_DIR ${CMAKE_BINARY_DIR}/doc)
set(TDEP_BUILD_DIR ${CMAKE_BINARY_DIR}/build)

file(MAKE_DIRECTORY ${TDEP_BIN_DIR})
file(MAKE_DIRECTORY ${TDEP_LIB_DIR})
file(MAKE_DIRECTORY ${TDEP_INC_DIR})
file(MAKE_DIRECTORY ${TDEP_MAN_DIR})
file(MAKE_DIRECTORY ${TDEP_MAN_DIR}/man1)
file(MAKE_DIRECTORY ${TDEP_DOC_DIR})
file(MAKE_DIRECTORY ${TDEP_BUILD_DIR})

# Configure Fortran compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -std=f2008 -cpp -fallow-argument-mismatch")
    if(ENABLE_OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
    endif()
    set(Fortran_MODULE_FLAG "-J")
    set(DEFAULT_REAL_8_FLAG "-fdefault-real-8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free -std08 -fpp")
    if(ENABLE_OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")
    endif()
    set(Fortran_MODULE_FLAG "-module")
    set(DEFAULT_REAL_8_FLAG "-r8")
endif()

# Set default optimization level
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbacktrace -fcheck=all")

# Set precompiler flags
set(PRECOMPILER_FLAGS "-DGPpng -Dgfortranprogressbar -DAGRESSIVE_SANITY")
if(USE_CGAL)
    set(PRECOMPILER_FLAGS "${PRECOMPILER_FLAGS} -Dusecgal")
endif()

# Include modules for finding dependencies
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Find external dependencies
find_package(LAPACK REQUIRED)
find_package(HDF5 REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(MPI COMPONENTS Fortran REQUIRED)

if(USE_CGAL)
    message(STATUS "CGAL support is enabled, but you will need to provide CGAL paths if needed")
endif()

# Create git information file for libolle
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/gitinformation.in
    ${CMAKE_CURRENT_BINARY_DIR}/gitinformation
    @ONLY
)

# Create bashrc_tdep
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/bashrc_tdep.in
    ${CMAKE_CURRENT_BINARY_DIR}/bashrc_tdep
    @ONLY
)

# Add subdirectories containing libraries and executables
add_subdirectory(src/libolle)
add_subdirectory(src/libflap)

# Define list of codes to build
set(TDEP_CODES
    anharmonic_free_energy
    phonon_dispersion_relations
    crystal_structure_info
    generate_structure
    canonical_configuration
    lineshape
    samples_from_md
    extract_forceconstants
    atomic_distribution
    pack_simulation
    refine_structure
    thermal_conductivity
    thermal_conductivity_2023
)

# Add CGAL-dependent code if CGAL is enabled
if(USE_CGAL)
    list(APPEND TDEP_CODES phasespace_surface)
endif()

# Function to add a TDEP program
function(add_tdep_program program_name)
    # Parse additional arguments for source files
    set(sources ${ARGN})
    if(NOT sources)
        # Default source files if none specified
        set(sources
            main.f90
            options.f90
        )
    endif()
    
    # Set build directory
    set(PROG_BUILD_DIR ${TDEP_BUILD_DIR}/${program_name})
    file(MAKE_DIRECTORY ${PROG_BUILD_DIR})
    
    # Create executable
    add_executable(${program_name} ${sources})
    
    # Set module output directory
    set_target_properties(${program_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROG_BUILD_DIR}
        Fortran_MODULE_DIRECTORY ${PROG_BUILD_DIR}
    )
    
    # Include directories
    target_include_directories(${program_name} PRIVATE
        ${TDEP_INC_DIR}/libolle
        ${TDEP_INC_DIR}/libflap
        ${HDF5_INCLUDE_DIRS}
        ${MPI_Fortran_INCLUDE_PATH}
    )
    
    # Link libraries
    target_link_libraries(${program_name} PRIVATE
        olle
        flap
        ${LAPACK_LIBRARIES}
        ${BLAS_LIBRARIES}
	${HDF5_FORTRAN_LIBRARIES}
        ${HDF5_LIBRARIES}
        MPI::MPI_Fortran
    )
    
    # Add custom command to create manpage if enabled
    if(GENERATE_MANPAGES)
        add_custom_command(
            TARGET ${program_name} POST_BUILD
            COMMAND ${PROG_BUILD_DIR}/${program_name} --manpage || true
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${program_name}.1 ${TDEP_MAN_DIR}/man1/ || true
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${program_name}.md ${TDEP_MAN_DIR}/ || true
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/${program_name}
            COMMENT "Generating manpage for ${program_name}"
        )
    endif()
    
    # Create symbolic link to bin directory
    add_custom_command(
        TARGET ${program_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink 
        ${PROG_BUILD_DIR}/${program_name}
        ${TDEP_BIN_DIR}/${program_name}
        COMMENT "Creating symbolic link for ${program_name}"
    )
    
    # Install targets
    install(TARGETS ${program_name}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endfunction()

# Add all programs 
foreach(code ${TDEP_CODES})
    file(GLOB code_sources src/${code}/*.f90)
    add_tdep_program(${code} ${code_sources})
endforeach()

# Installation
include(GNUInstallDirs)

install(DIRECTORY ${TDEP_BIN_DIR}/
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        USE_SOURCE_PERMISSIONS
        FILES_MATCHING PATTERN "*")

install(DIRECTORY ${TDEP_LIB_DIR}/
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        FILES_MATCHING PATTERN "*.a" PATTERN "*.so*")

install(DIRECTORY ${TDEP_MAN_DIR}/
        DESTINATION ${CMAKE_INSTALL_MANDIR}
        FILES_MATCHING PATTERN "*")

# Summary
message(STATUS "")
message(STATUS "TDEP Build Configuration:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  Git branch: ${GIT_BRANCH}")
message(STATUS "  Git revision: ${GIT_REVISION}")
message(STATUS "  Building shared libraries: ${BUILD_SHARED_LIBS}")
message(STATUS "  CGAL support: ${USE_CGAL}")
message(STATUS "  OpenMP support: ${ENABLE_OPENMP}")
message(STATUS "  Generate manpages: ${GENERATE_MANPAGES}")
message(STATUS "")