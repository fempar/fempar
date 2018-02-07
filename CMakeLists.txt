#################################################################
# HEADER
#################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

PROJECT(FEMPAR Fortran C)

SET(${PROJECT_NAME}_VERSION 0.0.1)
SET(${PROJECT_NAME}_SOVERSION 1)

SET(CMAKE_VERBOSE_MAKEFILE TRUE)

#################################################################
# DEFINE PATHS (IN)
#################################################################

SET(ROOT_PATH ${CMAKE_SOURCE_DIR})             # Fempar root path
SET(SRC_PATH ${ROOT_PATH}/Sources)             # Fempar Sources path

SET(LIBRARY_PATH ${SRC_PATH}/Lib)              # Fempar Lib path

SET(TESTS_PATH ${SRC_PATH}/Tests)              # Fempar Tests path
SET(SERIAL_TESTS_PATH ${TESTS_PATH}/Serial)    # Fempar serial tests path
SET(PAR_TESTS_PATH ${TESTS_PATH}/Par)          # Fempar parallel tests path
SET(SCRIPTS_TESTS_PATH ${TESTS_PATH}/Scripts)  # Fempar tests scripts path
SET(DATA_TESTS_PATH ${TESTS_PATH}/Data)        # Fempar tests data path

SET(INCLUDE_PATH ${SRC_PATH}/Include)          # Fempar include path
SET(GENERIC_PATH ${LIBRARY_PATH}/Generic)      # Fempar generic path

SET(CMAKE_PATH "${ROOT_PATH}/CMake")           # Fempar CMake path
SET(CMAKE_MODULE_PATH "${CMAKE_PATH}/Modules") # Fempar CMake Modules path

#################################################################
# BUILD PATHS (OUT)
#################################################################

SET(BINARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})                     # Fempar binary root path
SET(THIRDPARTY_OUTPUT_PATH ${BINARY_OUTPUT_PATH}/ThirdParty)      # Fempar binary thirdparty path
SET(LIBRARY_OUTPUT_PATH ${BINARY_OUTPUT_PATH}/lib)                # Fempar binary lib path
SET(TESTS_OUTPUT_PATH ${BINARY_OUTPUT_PATH}/tests)                # Fempar binary tests path
SET(EXECUTABLE_OUTPUT_PATH ${BINARY_OUTPUT_PATH}/bin)             # Fempar binary executables path
SET(SCRIPTS_OUTPUT_PATH ${BINARY_OUTPUT_PATH}/scripts)            # Fempar binary scripts path
SET(CMAKE_Fortran_MODULE_DIRECTORY ${BINARY_OUTPUT_PATH}/modules) # Fempar binary fortran modules path

#################################################################
# ADD INCLUDE DIRS
#################################################################

SET(${PROJECT_NAME}_INCLUDES ${INCLUDE_PATH} ${CMAKE_Fortran_MODULE_DIRECTORY})
INCLUDE_DIRECTORIES(${${PROJECT_NAME}_INCLUDES} ${GENERIC_PATH})
SET(${PROJECT_NAME}_INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS} ${${PROJECT_NAME}_INCLUDES})
UNSET(${PROJECT_NAME}_DEFINITIONS CACHE)

#################################################################
# LIBRARIES TO LINK
#################################################################

UNSET(${PROJECT_NAME}_EXTERNAL_LIBRARIES CACHE)

#################################################################
# CONFIGURATION TYPES & BUILD MODE & TESTING
#################################################################

SET(CMAKE_CONFIGURATION_TYPES DEBUG RELEASE) # Specifies the available build types
IF(NOT CMAKE_BUILD_TYPE)                     # If build type was not specified, set DEBUG as default
  SET(CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are: NONE DEBUG RELEASE"
      FORCE)

  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS NONE DEBUG RELEASE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(NOT ${PROJECT_NAME}_ENABLE_TESTS)         # If FEMPAR_ENABLE_TESTS was not defined, set OFF as default
    OPTION(${PROJECT_NAME}_ENABLE_TESTS "Enable/disable tests" OFF)
ENDIF(NOT ${PROJECT_NAME}_ENABLE_TESTS)

#################################################################
# List of external libraries
#################################################################

# List of external libraries
# LIBRARY ORDER IS IMPORTANT IN THE CASE OF MKL, BLAS AND LAPACK
# IF MKL WAS FOUND BLAS AND LAPACK MACRO DEFITIONS ARE ACTIVE 
# BUT LIBRARIES ARE NOT LINKED
SET(EXTERNAL_LIBRARIES MKL BLAS LAPACK GIDPOST METIS WSMP HSL_MI20 HSL_MA87 P4EST UMFPACK ZOLTAN ZLIB)
SET(MKL_USE_STATIC_LIBS ON)
MESSAGE(MKL_USE_STATIC_LIBS: ${MKL_USE_STATIC_LIBS})

#################################################################
# Find REQUIRED MPI library and check MPI.mod compatibility
#################################################################

FIND_PACKAGE(MPI REQUIRED) # Look for MPI library. This requirement is mandatory required

IF(MPI_FOUND)
    # Add MPI directories
    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH}) # Add MPI includes to fempar includes
    SET(${PROJECT_NAME}_INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS} ${MPI_Fortran_INCLUDE_PATH})

    FIND_PATH(MPI_MOD_FOUND mpi.mod  PATHS ${MPI_Fortran_INCLUDE_PATH}) # Look for mpi.mod file
    IF(MPI_MOD_FOUND)
        # Test mpi.mod compatibility with current CMAKE_Fortran_COMPILER
        TRY_COMPILE(FORTRAN_COMPILER_MATCH_MPI_MOD 
                    ${PROJECT_BINARY_DIR} 
                    ${TESTS_PATH}/try_compile_mpi_mod.f90	
                    CMAKE_FLAGS
                        "-DINCLUDE_DIRECTORIES:STRING=${MPI_Fortran_INCLUDE_PATH}"
                        "-DLINK_LIBRARIES:STRING=${MPI_Fortran_LIBRARIES}")
        IF(FORTRAN_COMPILER_MATCH_MPI_MOD)
            ADD_DEFINITIONS(-DMPI_MOD) # If mpi.mod is compatible add MPI_MOD preprocessor definition
        ELSE()
            MESSAGE(WARNING "MPI: mpi.mod file is NOT COMPATIBLE with this compiler! mpif.h is going be included instead")
            ADD_DEFINITIONS(-DMPI_H)   # If mpi.mod is NOT compatible add MPI_H preprocessor definition
        ENDIF()
    ELSE()
        MESSAGE(WARNING "MPI: mpi.mod file NOT FOUND! mpif.h is going to be included instead")
        ADD_DEFINITIONS(-DMPI_H)       # If mpi.mod was NOT found add MPI_H preprocessor definition
    ENDIF()
    SET(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARIES} ${MPI_LIBRARIES}) # Append MPI_LIBRARIES to the list of FEMPAR_LIBRARIES
ENDIF()

#################################################################
# Find external libraries
#################################################################

# Loop over the list of external libraries 
FOREACH(EXT_LIB ${EXTERNAL_LIBRARIES})
    # If external library is enabled but not found yet
    IF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND NOT ${EXT_LIB}_FOUND)
        IF(MKL_FOUND AND (${EXT_LIB} STREQUAL "BLAS" OR ${EXT_LIB} STREQUAL "LAPACK"))
            # IF MKL WAS FOUND DO NOT LINK WITH BLAS AND LAPACK LIBRARIES
            OPTION(${PROJECT_NAME}_ENABLE_${EXT_LIB} "Enable/disable ${EXT_LIB} library" ON)
            SET(${EXT_LIB}_FOUND TRUE)
            SET(${EXT_LIB}_INCLUDES "")
            SET(${EXT_LIB}_LIBS "")
            SET(${EXT_LIB}_INCLUDE_DIR "")
            SET(${EXT_LIB}_LIBRARIES "")
            ADD_DEFINITIONS(-DENABLE_${EXT_LIB})#${${EXT_LIB}_MAJOR_VERSION})
        ELSE()
            IF(${EXT_LIB} STREQUAL "HDF5")
                FIND_PACKAGE( ${EXT_LIB} COMPONENTS C Fortran HL Fortran_HL)
            ELSE()
                FIND_PACKAGE( ${EXT_LIB} )
            ENDIF()
    
            IF (${EXT_LIB}_FOUND)
                IF (DEFINED MKL_FOUND AND MKL_FOUND AND (${EXT_LIB} STREQUAL "BLAS" OR ${EXT_LIB} STREQUAL "LAPACK"))
                # IF MKL WAS FOUND DO NOT LINK WITH BLAS AND LAPACK LIBRARIES
                ELSE()
                    ADD_DEFINITIONS(-DENABLE_${EXT_LIB})
                    IF(IS_DIRECTORY ${${EXT_LIB}_INCLUDE_DIR})
                        INCLUDE_DIRECTORIES(${${EXT_LIB}_INCLUDE_DIR})
                        SET(${PROJECT_NAME}_INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS} ${${EXT_LIB}_INCLUDE_DIR})
                    ENDIF()
                ENDIF()
            ENDIF()
        ENDIF()

    # If external library is disabled
    ELSEIF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND NOT ${PROJECT_NAME}_ENABLE_${EXT_LIB})
        SET(${EXT_LIB}_INCLUDES "")
        SET(${EXT_LIB}_LIBS "")
        SET(${EXT_LIB}_INCLUDE_DIR "")
        SET(${EXT_LIB}_LIBRARIES "")

    # If user do not define if the external library is enabled or disabled, look for it as default
    ELSE(NOT DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB})
        IF(MKL_FOUND AND (${EXT_LIB} STREQUAL "BLAS" OR ${EXT_LIB} STREQUAL "LAPACK"))
            # IF MKL WAS FOUND DO NOT LINK WITH BLAS AND LAPACK LIBRARIES
            OPTION(${PROJECT_NAME}_ENABLE_${EXT_LIB} "Enable/disable ${EXT_LIB} library" ON)
            SET(${EXT_LIB}_FOUND TRUE)
            SET(${EXT_LIB}_INCLUDES "")
            SET(${EXT_LIB}_LIBS "")
            SET(${EXT_LIB}_INCLUDE_DIR "")
            SET(${EXT_LIB}_LIBRARIES "")
            ADD_DEFINITIONS(-DENABLE_${EXT_LIB})#${${EXT_LIB}_MAJOR_VERSION})
        ELSE()
            IF(${EXT_LIB} STREQUAL "HDF5")
                FIND_PACKAGE( ${EXT_LIB} COMPONENTS C Fortran HL Fortran_HL)
            ELSE()
                FIND_PACKAGE( ${EXT_LIB} )
            ENDIF()

#            STRING( TOUPPER ${EXT_LIB}_FOUND UPPER_STRING)

#            IF (${${UPPER_STRING}})
            IF (${${EXT_LIB}_FOUND})
                SET(NO_LIBS FALSE)
                OPTION(${PROJECT_NAME}_ENABLE_${EXT_LIB} "Enable/disable ${EXT_LIB} library" ON)
                ADD_DEFINITIONS(-DENABLE_${EXT_LIB})
                IF(IS_DIRECTORY ${${EXT_LIB}_INCLUDE_DIR})
                    INCLUDE_DIRECTORIES(${${EXT_LIB}_INCLUDE_DIR})
                    SET(${PROJECT_NAME}_INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS} ${${EXT_LIB}_INCLUDE_DIR})
                ENDIF()
            ELSE()
                OPTION(${PROJECT_NAME}_ENABLE_${EXT_LIB} "Enable/disable ${EXT_LIB} library" OFF)
                SET(${EXT_LIB}_INCLUDES "")
                SET(${EXT_LIB}_LIBS "")
                SET(${EXT_LIB}_INCLUDE_DIR "")
                SET(${EXT_LIB}_LIBRARIES "")
            ENDIF()
        ENDIF()
    ENDIF()
    SET(${PROJECT_NAME}_EXTERNAL_LIBRARIES ${${PROJECT_NAME}_EXTERNAL_LIBRARIES} ${${EXT_LIB}_LIBRARIES})
ENDFOREACH ()

SET(OPENMP_FLAGS "")
IF(DEFINED ${PROJECT_NAME}_ENABLE_OPENMP AND ${PROJECT_NAME}_ENABLE_OPENMP) # If OPENMP is ENABLED
    FIND_PACKAGE( OpenMP )
    IF(${OPENMP_FOUND})
        SET(OPENMP_FLAGS ${OpenMP_Fortran_FLAGS})
    ENDIF()
ENDIF()

#################################################################
# SET BUILD NAME (Before enable testing)
#################################################################

# Append the commit ID to the BUILDNAME if Git is found
SET(BUILDNAME ${CMAKE_Fortran_COMPILER_ID}_${CMAKE_BUILD_TYPE}_MKL=${${PROJECT_NAME}_ENABLE_MKL} CACHE STRING "" )

#################################################################
# ENABLE TESTING
#################################################################

IF(${PROJECT_NAME}_ENABLE_TESTS)
    ENABLE_TESTING()
    INCLUDE(CTest)
    SET(CTEST_COMMAND ${CMAKE_COMMAND} --build . --target test)
ELSE()
    SET(BUILD_TESTING OFF)
    UNSET(CTEST_COMMAND CACHE)
ENDIF()

#################################################################
# Add external projects
#################################################################

UNSET(${PROJECT_NAME}_EXTERNAL_PROJECTS CACHE)
SET(EXTERNAL_PROJECTS FPL F3PHT FLAP VTKIO SISL XH5For) # List of Fempar external projects

# Loop over all (required) external projects
# ${EXT_PROJ}_DIR is used as a clue to find every external project
FOREACH(EXT_PROJ ${EXTERNAL_PROJECTS})
    IF(NOT ${EXT_PROJ}_DIR AND ${PROJECT_NAME}_THIRDPARTY_DIR)
        # If the directory of the external project was not defined, build it from ThirdParty directory 
        SET(${EXT_PROJ}_DIR ${CMAKE_BINARY_DIR}/${${PROJECT_NAME}_THIRDPARTY_DIR}/${EXT_PROJ})
        IF(NOT EXISTS ${${EXT_PROJ}_DIR} OR NOT IS_DIRECTORY ${${EXT_PROJ}_DIR})
            GET_FILENAME_COMPONENT(${EXT_PROJ}_DIR ${${PROJECT_NAME}_THIRDPARTY_DIR}/${EXT_PROJ} ABSOLUTE BASE_DIR ${CMAKE_BINARY_DIR})
        ENDIF()
    ENDIF() 
    SET(${EXT_PROJ}_DIR ${${EXT_PROJ}_DIR} CACHE STRING "${EXT_PROJ} root directory")
    FIND_PACKAGE(${EXT_PROJ} REQUIRED)                    # Look for the external project
    SET(${PROJECT_NAME}_EXTERNAL_PROJECTS ${${PROJECT_NAME}_EXTERNAL_PROJECTS} ${${EXT_PROJ}_LIBS})    
    INCLUDE_DIRECTORIES(${${EXT_PROJ}_INCLUDES})          # Include directory of the external project
    SET(${PROJECT_NAME}_INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS} ${${EXT_PROJ}_INCLUDES})
    IF(EXT_PROJ STREQUAL "XH5For" AND XH5For_ENABLE_HDF5) # If XH5For is enabled add ENABLE_HDF5 preprocessor definition
        ADD_DEFINITIONS(-DENABLE_HDF5)
    ENDIF()
ENDFOREACH ()


#################################################################
# FFLAGS depend on the compiler and the build type
#################################################################

GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
#   SET(MACROS "${MACROS} -DDEBUG -Dmemcheck")
    ADD_DEFINITIONS(-DDEBUG)
    ADD_DEFINITIONS(-Dmemcheck)
ENDIF()

message(STATUS "COMPILER INFO: ${CMAKE_Fortran_COMPILER_ID} - ${Fortran_COMPILER_NAME}")

IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR Fortran_COMPILER_NAME MATCHES "gfortran*")
  # flags for GNU Fortran compiler
  set (FORTRAN_EXTRA_FLAGS "-fdefault-real-8 -ffree-line-length-0 -cpp -Wimplicit-interface ${OPENMP_FLAGS} ${FORTRAN_EXTRA_FLAGS}")
  set (CMAKE_Fortran_FLAGS "${FORTRAN_EXTRA_FLAGS} ${MACROS} ${INCLUDES} " CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -fbacktrace -fbounds-check -fprofile-arcs -ftest-coverage -Wimplicit-interface " CACHE STRING "" FORCE)
ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR Fortran_COMPILER_NAME MATCHES "ifort*")
  # flags for Intel Fortran compiler
  set (FORTRAN_EXTRA_FLAGS "-r8 -fpp -W1 -fPIC -assume realloc-lhs ${OPENMP_FLAGS}  ${FORTRAN_EXTRA_FLAGS}")
  set (CMAKE_Fortran_FLAGS "${FORTRAN_EXTRA_FLAGS} ${MACROS} ${INCLUDES}" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -traceback -g -debug all -check all,noarg_temp_created -ftrapuv -warn nointerfaces" CACHE STRING "" FORCE)
# A partir de CMake 3.1
# -prof-gen:srcpos -prof-dir${PROJECT_BINARY_DIR}
ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "XL" OR Fortran_COMPILER_NAME MATCHES "xlf*")
  # Flags for IBM XLF compiler
  set(FORTRAN_EXTRA_FLAGS "-q64 -qrealsize=8 -qsuffix=f=f90:cpp=f90 ${OPENMP_FLAGS} ${FORTRAN_EXTRA_FLAGS} ")
  set (CMAKE_Fortran_FLAGS "${FORTRAN_EXTRA_FLAGS} ${MACROS} ${INCLUDES}" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -qstrict" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -qfullpath -qkeepparm" CACHE STRING "" FORCE)
ELSE ()
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
ENDIF ()

SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--allow-multiple-definition")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition")

message (STATUS "CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
message (STATUS "CMAKE_Fortran_FLAGS: " ${CMAKE_Fortran_FLAGS})
message (STATUS "CMAKE_Fortran_FLAGS_RELEASE: " ${CMAKE_Fortran_FLAGS_RELEASE})
message (STATUS "CMAKE_Fortran_FLAGS_DEBUG: " ${CMAKE_Fortran_FLAGS_DEBUG})

#################################################################
# STATIC LIBRARIES
#################################################################

# If shared libs not forced, prepend ".a" extension to library suffixes. 
# Useful to find first static libraries in future calls to FIND_PACKAGE
IF(NOT ${BUILD_SHARED_LIBS})
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a ${CMAKE_FIND_LIBRARY_SUFFIXES}")
ENDIF()

#################################################################
# ADD LIBRARY SOURCE SUBDIRS
#################################################################

ADD_SUBDIRECTORY(${LIBRARY_PATH})   # Add the Lib subdirectory

#################################################################
# ADD TESTS (IF ENABLED) AND DRIVERS SOURCE SUBDIRS
#################################################################

IF(${PROJECT_NAME}_ENABLE_TESTS)
    ADD_SUBDIRECTORY(${TESTS_PATH}) # Add the Tests subdirectory
ENDIF()

#################################################################
# CREATE TAGS
# SOUCERS LIST PICKED FROM OTHER CMAKELISTS
#################################################################

FIND_PROGRAM(ETAGS etags)
IF(ETAGS)
    ADD_CUSTOM_TARGET(TAGS ALL COMMAND ${ETAGS} ${LIBRARY_SRC} ${TESTS_SERIAL_SRC} ${TESTS_PAR_SRC})
ENDIF()

#################################################################
# CREATE CONFIG FILE
#################################################################

CONFIGURE_FILE(${CMAKE_PATH}/${PROJECT_NAME}Config.cmake.in 
               ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
               @ONLY)
