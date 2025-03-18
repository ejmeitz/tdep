# Set module output directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${TDEP_INC_DIR}/libflap)

# Ensure position-independent code for shared libraries
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Define source files IN THE EXACT ORDER they should be compiled
set(FLAP_SOURCES
    penf.F90                              # 1st
    flap_object_t.f90                     # 2nd
    flap_utils_m.f90                      # 3rd
    flap_command_line_argument_t.F90      # 4th
    flap_command_line_arguments_group_t.f90  # 5th
    flap_command_line_interface_t.F90     # 6th
    flap.f90                              # 7th
)

# Create the library
add_library(flap SHARED ${FLAP_SOURCES})

# Set library properties
set_target_properties(flap PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    OUTPUT_NAME "flap"
    LIBRARY_OUTPUT_DIRECTORY ${TDEP_LIB_DIR}
    ARCHIVE_OUTPUT_DIRECTORY ${TDEP_LIB_DIR}
    POSITION_INDEPENDENT_CODE ON
)

# Add explicit module dependencies to ensure correct compilation order
set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/flap_object_t.f90
    PROPERTIES OBJECT_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/penf.F90")

set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/flap_utils_m.f90
    PROPERTIES OBJECT_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/penf.F90")

set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/flap_command_line_argument_t.F90
    PROPERTIES OBJECT_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/penf.F90;${CMAKE_CURRENT_SOURCE_DIR}/flap_object_t.f90;${CMAKE_CURRENT_SOURCE_DIR}/flap_utils_m.f90")

set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/flap_command_line_arguments_group_t.f90
    PROPERTIES OBJECT_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/penf.F90;${CMAKE_CURRENT_SOURCE_DIR}/flap_command_line_argument_t.F90")

set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/flap_command_line_interface_t.F90
    PROPERTIES OBJECT_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/penf.F90;${CMAKE_CURRENT_SOURCE_DIR}/flap_command_line_argument_t.F90;${CMAKE_CURRENT_SOURCE_DIR}/flap_command_line_arguments_group_t.f90")

set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/flap.f90
    PROPERTIES OBJECT_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/flap_command_line_argument_t.F90;${CMAKE_CURRENT_SOURCE_DIR}/flap_command_line_arguments_group_t.f90;${CMAKE_CURRENT_SOURCE_DIR}/flap_command_line_interface_t.F90")

# Disable parallel builds for this target to ensure correct module dependency handling
set_property(TARGET flap PROPERTY JOB_POOL_COMPILE sequential_job_pool)
set_property(GLOBAL PROPERTY JOB_POOLS sequential_job_pool=1)

# Include directories
target_include_directories(flap PUBLIC
    ${TDEP_INC_DIR}/libolle
    ${TDEP_INC_DIR}/libflap
)

# Compiler definitions
target_compile_definitions(flap PRIVATE ${PRECOMPILER_FLAGS})

# Add explicit -fPIC flag for Fortran compiler
target_compile_options(flap PRIVATE -fPIC)

# Install targets
install(TARGETS flap
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})