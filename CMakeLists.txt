#################################################################
# HEADER
#################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

PROJECT(FEMPAR Fortran C)

SET(${PROJECT_NAME}_VERSION 0.0.1)
SET(${PROJECT_NAME}_SOVERSION 1)
SET(FEM_LIB fem) # FEM library
SET(PAR_LIB par) # PAR library

set(CTEST_PROJECT_NAME "Test_victorsndvg")
set(CTEST_NIGHTLY_START_TIME "00:00:00 EST")

set(CTEST_DROP_METHOD "http" )
set(CTEST_DROP_SITE "my.cdash.org" )
set(CTEST_DROP_LOCATION "/submit.php?project=Test_victorsndvg" )
set(CTEST_DROP_SITE_CDASH TRUE )
#set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full" )

SET(CMAKE_VERBOSE_MAKEFILE TRUE)

message ("-----------------------------------------")

#################################################################
# DEFINE PATHS
#################################################################

SET(SRC_PATH ${CMAKE_SOURCE_DIR}/Sources)
SET(FEM_PATH ${SRC_PATH}/Fem)
SET(PAR_PATH ${SRC_PATH}/Par)
SET(DRIVERS_PATH ${SRC_PATH}/Drivers)
SET(GENERIC_PATH ${SRC_PATH}/Generic)
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/Modules/")


#################################################################
# BUILD PATHS
#################################################################

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

#################################################################
# ADD INCLUDE DIRS
#################################################################

include_directories(${GENERIC_PATH})

#################################################################
# CONFIGURATION TYPES & BUILD MODE
#################################################################

SET(CMAKE_CONFIGURATION_TYPES DEBUG RELEASE)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are: NONE DEBUG RELEASE"
      FORCE)

  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS NONE DEBUG RELEASE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#################################################################
# List of external libraries
#################################################################

# LIBRARY ORDER IS IMPORTANT IN THE CASE OF MKL, BLAS AND LAPACK
# IF MKL WAS FOUND BLAS AND LAPACK MACRO DEFITIONS ARE ACTIVE 
# BUT LIBRARIES ARE NOT LINKED
SET(EXT_LIBS MKL BLAS LAPACK GIDPOST METIS WSMP HSL_MI20 HSL_MA87 P4EST UMFPACK ZOLTAN)
SET(MKL_FIND_QUIETLY TRUE) 

#################################################################
# Find packages
#################################################################

SET(NO_LIBS TRUE)

FOREACH(EXT_LIB ${EXT_LIBS})
    IF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND NOT ${EXT_LIB}_FOUND)
        IF(MKL_FOUND AND (${EXT_LIB} STREQUAL "BLAS" OR ${EXT_LIB} STREQUAL "LAPACK"))
#           IF MKL WAS FOUND DO NOT LINK WITH BLAS AND LAPACK LIBRARIES
            OPTION(${PROJECT_NAME}_ENABLE_${EXT_LIB} "Enable/disable ${EXT_LIB} library" ON)
            SET(${EXT_LIB}_FOUND TRUE)
            ADD_DEFINITIONS(-DENABLE_${EXT_LIB}${${EXT_LIB}_MAJOR_VERSION})
        ELSE()
            FIND_PACKAGE( ${EXT_LIB} )
    
            IF (${EXT_LIB}_FOUND)
                SET(NO_LIBS FALSE)
                IF (MKL_FOUND AND (${EXT_LIB} STREQUAL "BLAS" OR ${EXT_LIB} STREQUAL "LAPACK"))
#                   IF MKL WAS FOUND DO NOT LINK WITH BLAS AND LAPACK LIBRARIES
                ELSE()
#                  SET(MACROS "${MACROS} -DENABLE_${EXT_LIB}")
                   ADD_DEFINITIONS(-DENABLE_${EXT_LIB}${${EXT_LIB}_MAJOR_VERSION})
                   IF(IS_DIRECTORY ${${EXT_LIB}_INCLUDE_DIR})
                       INCLUDE_DIRECTORIES(${${EXT_LIB}_INCLUDE_DIR})
                   ENDIF()
#                  SET(INCLUDES "${INCLUDES} -I${${EXT_LIB}_INCLUDE_DIR}" )
                   LINK_DIRECTORIES(${${EXT_LIB}_LIBRARIES})
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
            ADD_DEFINITIONS(-DENABLE_${EXT_LIB}${${EXT_LIB}_MAJOR_VERSION})
        ELSE()
            FIND_PACKAGE( ${EXT_LIB} )

            IF (${EXT_LIB}_FOUND)
               SET(NO_LIBS FALSE)
               OPTION(${PROJECT_NAME}_ENABLE_${EXT_LIB} "Enable/disable ${EXT_LIB} library" ON)
               MESSAGE(STATUS "External ${EXT_LIB} library Found!")
#              SET(MACROS "${MACROS} -DENABLE_${EXT_LIB}")
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

#################################################################
# FFLAGS depend on the compiler and the build type
#################################################################

GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
   SET(MACROS "${MACROS} -DDEBUG -Dmemcheck")
ENDIF()

IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR Fortran_COMPILER_NAME MATCHES "gfortran*")
  # gfortran 
  set(FORTRAN_FLAGS "-fdefault-real-8 -ffree-line-length-0 -cpp -Wimplicit-interface ")
  set (CMAKE_Fortran_FLAGS "${FORTRAN_FLAGS} ${MACROS} ${INCLUDES}" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -fbacktrace -fbounds-check -fprofile-arcs -ftest-coverage" CACHE STRING "" FORCE)
ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR Fortran_COMPILER_NAME MATCHES "ifort*")
  # ifort (untested)
  set(FORTRAN_FLAGS "-r8 -fpp -W1 -warn nointerfaces")
  set (CMAKE_Fortran_FLAGS "${FORTRAN_FLAGS} ${MACROS} ${INCLUDES}" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -traceback -g -debug all -check all -ftrapuv" CACHE STRING "" FORCE)
ELSE ()
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
ENDIF ()

message (STATUS "CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
message (STATUS "CMAKE_Fortran_FLAGS: " ${CMAKE_Fortran_FLAGS})
message (STATUS "CMAKE_Fortran_FLAGS_RELEASE: " ${CMAKE_Fortran_FLAGS_RELEASE})
message (STATUS "CMAKE_Fortran_FLAGS_DEBUG: " ${CMAKE_Fortran_FLAGS_DEBUG})

#################################################################
# ENABLE TESTING
#################################################################

STRING(TIMESTAMP DATE "%Y-%m-%d_%H:%M")
SET(BUILDNAME ${DATE}_${CMAKE_Fortran_COMPILER_ID}_${CMAKE_BUILD_TYPE}_NOLIBS=${NO_LIBS} CACHE STRING "" )
ENABLE_TESTING()
INCLUDE(CTest)

#################################################################
# ADD SOURCE SUBDIRS
#################################################################


ADD_SUBDIRECTORY(${FEM_PATH})
ADD_SUBDIRECTORY(${PAR_PATH})
ADD_SUBDIRECTORY(${DRIVERS_PATH})


message ("-----------------------------------------")
