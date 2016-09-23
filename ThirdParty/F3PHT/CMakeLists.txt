#################################################################
# HEADER
#################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7)
PROJECT(F3PHT Fortran)

SET(${PROJECT_NAME}_VERSION 0.0.1)
SET(${PROJECT_NAME}_SOVERSION 1)
SET(LIB ${PROJECT_NAME})

#SET(CMAKE_VERBOSE_MAKEFILE TRUE)

#################################################################
# DEFINE PATHS
#################################################################

SET(CMAKE_PATH ${CMAKE_SOURCE_DIR}/CMake)
SET(SRC_PATH ${CMAKE_SOURCE_DIR}/src)
SET(LIB_PATH ${SRC_PATH}/lib)
SET(TESTS_PATH ${SRC_PATH}/tests)

#SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/Modules/")


#################################################################
# BUILD PATHS
#################################################################

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
SET(MODULE_OUTPUT_PATH ${CMAKE_Fortran_MODULE_DIRECTORY})
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

#################################################################
# CONFIGURATION TYPES & BUILD MODE & BUILD_TESTS
#################################################################

SET(CMAKE_CONFIGURATION_TYPES DEBUG RELEASE)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are: NONE DEBUG RELEASE"
      FORCE)

  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS NONE DEBUG RELEASE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(NOT ${PROJECT_NAME}_ENABLE_TESTS)
    OPTION(${PROJECT_NAME}_ENABLE_TESTS "Enable/disable tests compilation" OFF)
ENDIF(NOT ${PROJECT_NAME}_ENABLE_TESTS)

#################################################################
# List of external libraries
#################################################################

SET(EXT_LIBS  )

#################################################################
# Find packages
#################################################################

FOREACH(EXT_LIB ${EXT_LIBS})
    IF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND NOT ${EXT_LIB}_FOUND)
        FIND_PACKAGE( ${EXT_LIB} )
    
        IF (${EXT_LIB}_FOUND)
            ADD_DEFINITIONS(-DENABLE_${EXT_LIB})
            IF(IS_DIRECTORY ${${EXT_LIB}_INCLUDE_DIR})
               INCLUDE_DIRECTORIES(${${EXT_LIB}_INCLUDE_DIR})
            ENDIF()
            LINK_DIRECTORIES(${${EXT_LIB}_LIBRARIES})
        ENDIF()

    ELSEIF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND NOT ${PROJECT_NAME}_ENABLE_${EXT_LIB})
        UNSET(${EXT_LIB}_INCLUDES CACHE)
        UNSET(${EXT_LIB}_LIBS CACHE)
        UNSET(${EXT_LIB}_INCLUDE_DIR CACHE)
        UNSET(${EXT_LIB}_LIBRARIES CACHE)

    ELSE(NOT DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB})
        FIND_PACKAGE( ${EXT_LIB} )

        IF (${${EXT_LIB}_FOUND})
           SET(NO_LIBS FALSE)
           OPTION(${PROJECT_NAME}_ENABLE_${EXT_LIB} "Enable/disable ${EXT_LIB} library" ON)
           ADD_DEFINITIONS(-DENABLE_${EXT_LIB})
           IF(IS_DIRECTORY ${${EXT_LIB}_INCLUDE_DIR})
               INCLUDE_DIRECTORIES(${${EXT_LIB}_INCLUDE_DIR})
           ENDIF()
            LINK_DIRECTORIES(${${EXT_LIB}_LIBRARIES})
        ELSE()
           OPTION(${PROJECT_NAME}_ENABLE_${EXT_LIB} "Enable/disable ${EXT_LIB} library" OFF)
           UNSET(${EXT_LIB}_INCLUDES CACHE)
           UNSET(${EXT_LIB}_LIBS CACHE)
           UNSET(${EXT_LIB}_INCLUDE_DIR CACHE)
           UNSET(${EXT_LIB}_LIBRARIES CACHE)
        ENDIF()
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

ADD_DEFINITIONS(-D${CMAKE_Fortran_COMPILER_ID})

message(STATUS "COMPILER INFO: ${CMAKE_Fortran_COMPILER_ID} - ${Fortran_COMPILER_NAME}")

IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR Fortran_COMPILER_NAME MATCHES "gfortran*")
  # gfortran 
  set(FORTRAN_FLAGS "-fdefault-real-8 -ffree-line-length-0 -cpp -Wimplicit-interface ${EXTRA_FLAGS} ")
  set (CMAKE_Fortran_FLAGS "${FORTRAN_FLAGS} ${MACROS} ${INCLUDES} " CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -fbacktrace -fbounds-check -fprofile-arcs -ftest-coverage -Wimplicit-interface ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR Fortran_COMPILER_NAME MATCHES "ifort*")
  # ifort (untested)
  set(FORTRAN_FLAGS "-r8 -fpp -W1 ${EXTRA_FLAGS} ")
  set (CMAKE_Fortran_FLAGS "${FORTRAN_FLAGS} ${MACROS} ${INCLUDES}" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -traceback -g -debug all -check all -ftrapuv -warn nointerfaces ${EXTRA_FLAGS} " CACHE STRING "" FORCE)
# A partir de CMake 3.1
# -prof-gen:srcpos -prof-dir${PROJECT_BINARY_DIR}

ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL "XL" OR Fortran_COMPILER_NAME MATCHES "xlf*")
  # xlf (untested)
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
# STATIC LIBRARIES
#################################################################
# Try to search first static libraries
IF(NOT ${BUILD_SHARED_LIBS})
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a ${CMAKE_FIND_LIBRARY_SUFFIXES}")
ENDIF()

#################################################################
# ENABLE TESTING
#################################################################

SET(BUILDNAME ${CMAKE_Fortran_COMPILER_ID}_${CMAKE_BUILD_TYPE}_MKL=${${PROJECT_NAME}_ENABLE_MKL} CACHE STRING "" )
IF(${PROJECT_NAME}_ENABLE_TESTS)
    ENABLE_TESTING()
    INCLUDE(CTest)
ENDIF()

#################################################################
# ADD SOURCE SUBDIRS
#################################################################

ADD_SUBDIRECTORY(${LIB_PATH})

#################################################################
# CREATE CONFIG FILE
#################################################################

CONFIGURE_FILE(${CMAKE_PATH}/${PROJECT_NAME}Config.cmake.in 
               ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
               @ONLY)

