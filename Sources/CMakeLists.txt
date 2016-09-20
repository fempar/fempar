#################################################################
# HEADER
#################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

PROJECT(FEMPAR Fortran C)

SET(${PROJECT_NAME}_VERSION 0.0.1)
SET(${PROJECT_NAME}_SOVERSION 1)
SET(${PROJECT_NAME}_USER "" CACHE STRING "")

SET(CMAKE_VERBOSE_MAKEFILE TRUE)

#################################################################
# DEFINE PATHS
#################################################################

SET(ROOT_PATH ${CMAKE_SOURCE_DIR}/..)
SET(SRC_PATH ${CMAKE_SOURCE_DIR})

SET(LIBRARY_PATH ${SRC_PATH}/Lib)

SET(TESTS_PATH ${SRC_PATH}/Tests)
SET(SERIAL_TESTS_PATH ${TESTS_PATH}/Serial)
SET(PAR_TESTS_PATH ${TESTS_PATH}/Par)
SET(SCRIPTS_TESTS_PATH ${TESTS_PATH}/Scripts)
SET(DATA_TESTS_PATH ${TESTS_PATH}/Data)

SET(INCLUDE_PATH ${SRC_PATH}/Include)
SET(GENERIC_PATH ${LIBRARY_PATH}/Generic)

SET(THIRDPARTY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/ThirdParty)

SET(CMAKE_MODULE_PATH "${ROOT_PATH}/CMake/Modules/")

#################################################################
# BUILD PATHS
#################################################################

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(SCRIPTS_OUTPUT_PATH ${PROJECT_BINARY_DIR}/scripts)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

#################################################################
# ADD INCLUDE DIRS
#################################################################

INCLUDE_DIRECTORIES(${INCLUDE_PATH} ${GENERIC_PATH})

#################################################################
# CONFIGURATION TYPES & BUILD MODE & TESTING
#################################################################

SET(CMAKE_CONFIGURATION_TYPES DEBUG RELEASE)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are: NONE DEBUG RELEASE"
      FORCE)

  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS NONE DEBUG RELEASE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(NOT ${PROJECT_NAME}_ENABLE_TESTS)
    OPTION(${PROJECT_NAME}_ENABLE_TESTS "Enable/disable tests" OFF)
ENDIF(NOT ${PROJECT_NAME}_ENABLE_TESTS)

#################################################################
# Add external projects
#################################################################

SET(EXT_PROJS FPL F3PHT FLAP VTKIO)
include(ExternalProject)

FOREACH(EXT_PROJ ${EXT_PROJS})

    SET(${EXT_PROJ}_SRC_PATH ${ROOT_PATH}/ThirdParty/${EXT_PROJ})
    SET(${EXT_PROJ}_BINARY_PATH ${THIRDPARTY_OUTPUT_PATH}/${EXT_PROJ})
    EXTERNALPROJECT_Add(${EXT_PROJ}
        PREFIX ${EXT_PROJ}
        DOWNLOAD_COMMAND ""
        SOURCE_DIR ${${EXT_PROJ}_SRC_PATH}
        BINARY_DIR ${${EXT_PROJ}_BINARY_PATH}
        INSTALL_DIR ${INSTALL_PREFIX}
        CONFIGURE_COMMAND cmake -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -D${EXT_PROJ}_ENABLE_TESTS=${${PROJECT_NAME}_ENABLE_TESTS} ${${EXT_PROJ}_SRC_PATH}
        BUILD_COMMAND cmake --build .
        TEST_COMMAND ${CTEST_COMMAND}
        INSTALL_COMMAND ""
    )
    INCLUDE_DIRECTORIES(${${EXT_PROJ}_BINARY_PATH}/modules)
    SET(EXT_PROJ_LIBS ${EXT_PROJ_LIBS} ${${EXT_PROJ}_BINARY_PATH}/lib/lib${EXT_PROJ}.a)
    IF(${EXT_PROJ} STREQUAL "FLAP")
        SET(EXT_PROJ_LIBS ${EXT_PROJ_LIBS} ${${EXT_PROJ}_BINARY_PATH}/third_party/PENF/lib/libPENF.a)
    ENDIF()

ENDFOREACH ()

#################################################################
# SISL needs a different treatment because it's a C project not owned by us
SET(SISL_LIB SISL)
SET(EXT_PROJS ${EXT_PROJS} ${SISL_LIB} )

SET(${SISL_LIB}_SRC_PATH ${ROOT_PATH}/ThirdParty/${SISL_LIB})
SET(${SISL_LIB}_BINARY_PATH ${THIRDPARTY_OUTPUT_PATH}/${SISL_LIB})
EXTERNALPROJECT_Add(${SISL_LIB}
    PREFIX ${SISL_LIB}
    DOWNLOAD_COMMAND ""
    SOURCE_DIR ${${SISL_LIB}_SRC_PATH}
    BINARY_DIR ${${SISL_LIB}_BINARY_PATH}
    INSTALL_DIR ${INSTALL_PREFIX}
    CONFIGURE_COMMAND cmake -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ${${SISL_LIB}_SRC_PATH}
    BUILD_COMMAND cmake --build .
    TEST_COMMAND ""
    INSTALL_COMMAND ""
)
INCLUDE_DIRECTORIES(${${SISL_LIB}_SRC_PATH}/include)
SET(EXT_PROJ_LIBS ${EXT_PROJ_LIBS} ${${SISL_LIB}_BINARY_PATH}/libsisl.a)

#################################################################
# List of external libraries
#################################################################

# LIBRARY ORDER IS IMPORTANT IN THE CASE OF MKL, BLAS AND LAPACK
# IF MKL WAS FOUND BLAS AND LAPACK MACRO DEFITIONS ARE ACTIVE 
# BUT LIBRARIES ARE NOT LINKED
SET(EXT_LIBS MKL BLAS LAPACK GIDPOST METIS WSMP HSL_MI20 HSL_MA87 P4EST UMFPACK ZOLTAN ZLIB )
SET(MKL_FIND_QUIETLY TRUE) 

#################################################################
# Find REQUIRED MPI library and check MPI.mod compatibility
#################################################################

FIND_PACKAGE(MPI REQUIRED)

IF(MPI_FOUND)
    # Add MPI directories
    ADD_DEFINITIONS(${MPI_Fortran_COMPILE_FLAGS})
    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})

    # Find mpi.mod
    FIND_PATH(MPI_MOD_FOUND mpi.mod  PATHS ${MPI_Fortran_INCLUDE_PATH})
    IF(MPI_MOD_FOUND)
        TRY_COMPILE(FORTRAN_COMPILER_MATCH_MPI_MOD 
                    ${PROJECT_BINARY_DIR} 
                    ${TESTS_PATH}/try_compile_mpi_mod.f90	
                    CMAKE_FLAGS
                        "-DINCLUDE_DIRECTORIES:STRING=${MPI_Fortran_INCLUDE_PATH}"
                        "-DLINK_LIBRARIES:STRING=${MPI_Fortran_LIBRARIES}")
        IF(FORTRAN_COMPILER_MATCH_MPI_MOD)
            ADD_DEFINITIONS(-DMPI_MOD)
        ELSE()
            MESSAGE(WARNING "MPI: mpi.mod file is NOT COMPATIBLE with this compiler! mpif.h is going be included instead")
            ADD_DEFINITIONS(-DMPI_H)            
        ENDIF()
    ELSE()
        MESSAGE(WARNING "MPI: mpi.mod file NOT FOUND! mpif.h is going to be included instead")
        ADD_DEFINITIONS(-DMPI_H)
    ENDIF()
ENDIF()

#################################################################
# Find external libraries
#################################################################

FOREACH(EXT_LIB ${EXT_LIBS})
    IF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND NOT ${EXT_LIB}_FOUND)
        IF(MKL_FOUND AND (${EXT_LIB} STREQUAL "BLAS" OR ${EXT_LIB} STREQUAL "LAPACK"))
#           IF MKL WAS FOUND DO NOT LINK WITH BLAS AND LAPACK LIBRARIES
            OPTION(${PROJECT_NAME}_ENABLE_${EXT_LIB} "Enable/disable ${EXT_LIB} library" ON)
            SET(${EXT_LIB}_FOUND TRUE)
            ADD_DEFINITIONS(-DENABLE_${EXT_LIB})#${${EXT_LIB}_MAJOR_VERSION})
        ELSE()
            FIND_PACKAGE( ${EXT_LIB} )
    
            IF (${EXT_LIB}_FOUND)
                IF (DEFINED MKL_FOUND AND MKL_FOUND AND (${EXT_LIB} STREQUAL "BLAS" OR ${EXT_LIB} STREQUAL "LAPACK"))
#                   IF MKL WAS FOUND DO NOT LINK WITH BLAS AND LAPACK LIBRARIES
                ELSE()
                   ADD_DEFINITIONS(-DENABLE_${EXT_LIB})
                   IF(IS_DIRECTORY ${${EXT_LIB}_INCLUDE_DIR})
                       INCLUDE_DIRECTORIES(${${EXT_LIB}_INCLUDE_DIR})
                   ENDIF()
                ENDIF()
            ENDIF()
        ENDIF()


    ELSEIF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND NOT ${PROJECT_NAME}_ENABLE_${EXT_LIB})
        UNSET(${EXT_LIB}_INCLUDES CACHE)
        UNSET(${EXT_LIB}_LIBS CACHE)
        UNSET(${EXT_LIB}_INCLUDE_DIR CACHE)
        UNSET(${EXT_LIB}_LIBRARIES CACHE)

    ELSE(NOT DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB})
        IF(MKL_FOUND AND (${EXT_LIB} STREQUAL "BLAS" OR ${EXT_LIB} STREQUAL "LAPACK"))
#           IF MKL WAS FOUND DO NOT LINK WITH BLAS AND LAPACK LIBRARIES
            OPTION(${PROJECT_NAME}_ENABLE_${EXT_LIB} "Enable/disable ${EXT_LIB} library" ON)
            SET(${EXT_LIB}_FOUND TRUE)
            ADD_DEFINITIONS(-DENABLE_${EXT_LIB})#${${EXT_LIB}_MAJOR_VERSION})
        ELSE()
            FIND_PACKAGE( ${EXT_LIB} )

#            STRING( TOUPPER ${EXT_LIB}_FOUND UPPER_STRING)

#            IF (${${UPPER_STRING}})
            IF (${${EXT_LIB}_FOUND})
               SET(NO_LIBS FALSE)
               OPTION(${PROJECT_NAME}_ENABLE_${EXT_LIB} "Enable/disable ${EXT_LIB} library" ON)
               ADD_DEFINITIONS(-DENABLE_${EXT_LIB})
               IF(IS_DIRECTORY ${${EXT_LIB}_INCLUDE_DIR})
                   INCLUDE_DIRECTORIES(${${EXT_LIB}_INCLUDE_DIR})
               ENDIF()
            ELSE()
               OPTION(${PROJECT_NAME}_ENABLE_${EXT_LIB} "Enable/disable ${EXT_LIB} library" OFF)
               UNSET(${EXT_LIB}_INCLUDES CACHE)
               UNSET(${EXT_LIB}_LIBS CACHE)
               UNSET(${EXT_LIB}_INCLUDE_DIR CACHE)
               UNSET(${EXT_LIB}_LIBRARIES CACHE)
            ENDIF()
        ENDIF()
    ENDIF()
ENDFOREACH ()

FIND_PACKAGE( OpenMP )
IF(${OPENMP_FOUND})
  SET(EXTRA_FLAGS ${OpenMP_Fortran_FLAGS})
ELSE()
  SET(EXTRA_FLAGS "")
ENDIF()

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
  # gfortran 
  set(FORTRAN_FLAGS "-fdefault-real-8 -ffree-line-length-0 -cpp -Wimplicit-interface ${EXTRA_FLAGS} ")
  set (CMAKE_Fortran_FLAGS "${FORTRAN_FLAGS} ${MACROS} ${INCLUDES} " CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -fbacktrace -fbounds-check -fprofile-arcs -ftest-coverage -Wimplicit-interface ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR Fortran_COMPILER_NAME MATCHES "ifort*")
  # ifort
  set(FORTRAN_FLAGS "-r8 -fpp -W1 -fPIC -assume realloc-lhs ${EXTRA_FLAGS} ")
  set (CMAKE_Fortran_FLAGS "${FORTRAN_FLAGS} ${MACROS} ${INCLUDES}" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -traceback -g -debug all -check all,noarg_temp_created -ftrapuv -warn nointerfaces ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
# A partir de CMake 3.1
# -prof-gen:srcpos -prof-dir${PROJECT_BINARY_DIR}
ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "XL" OR Fortran_COMPILER_NAME MATCHES "xlf*")
  # xlf
  set(FORTRAN_FLAGS "-q64 -qrealsize=8 -qsuffix=f=f90:cpp=f90 ${EXTRA_FLAGS} ")
  set (CMAKE_Fortran_FLAGS "${FORTRAN_FLAGS} ${MACROS} ${INCLUDES}" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -qstrict ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -qfullpath -qkeepparm ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
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
# SET BUILD NAME
#################################################################

FIND_PACKAGE( Git )
IF(${GIT_FOUND})
    # Get the current working branch
    EXECUTE_PROCESS(
        COMMAND git rev-parse --short HEAD
        WORKING_DIRECTORY ${ROOT_PATH}
        OUTPUT_VARIABLE GIT_COMMIT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
ENDIF()

SET(BUILDNAME ${CMAKE_Fortran_COMPILER_ID}_${CMAKE_BUILD_TYPE}_MKL=${${PROJECT_NAME}_ENABLE_MKL}@${GIT_COMMIT} CACHE STRING "" )

#################################################################
# ENABLE TESTING
#################################################################

IF(${PROJECT_NAME}_ENABLE_TESTS)
    ENABLE_TESTING()
    INCLUDE(CTest)
    SET(CTEST_COMMAND ctest)
ELSE()
    SET(BUILD_TESTING OFF)
    UNSET(CTEST_COMMAND CACHE)
ENDIF()

#################################################################
# STATIC LIBRARIES
#################################################################
# Try to search first static libraries
IF(NOT ${BUILD_SHARED_LIBS})
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a ${CMAKE_FIND_LIBRARY_SUFFIXES}")
ENDIF()

#################################################################
# ADD LIBRARY SOURCE SUBDIRS
#################################################################

ADD_SUBDIRECTORY(${LIBRARY_PATH})

#################################################################
# ADD TESTS (IF ENABLED) AND DRIVERS SOURCE SUBDIRS
#################################################################

IF(${PROJECT_NAME}_ENABLE_TESTS)
    ADD_SUBDIRECTORY(${SERIAL_TESTS_PATH})
    ADD_SUBDIRECTORY(${PAR_TESTS_PATH})
ENDIF()

#################################################################
# CREATE TAGS
# SOUCERS LIST PICKED FROM OTHER CMAKELISTS
#################################################################

FIND_PROGRAM(ETAGS etags)
IF(ETAGS)
    ADD_CUSTOM_TARGET(TAGS ALL COMMAND ${ETAGS} ${LIBRARY_SRC} ${TESTS_SERIAL_SRC} ${TESTS_PAR_SRC})
ENDIF()

