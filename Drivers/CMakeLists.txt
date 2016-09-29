#################################################################
# HEADER
#################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

PROJECT(FEMPAR_DRIVER Fortran C)

SET(CMAKE_VERBOSE_MAKEFILE TRUE)

#################################################################
# BUILD PATHS
#################################################################

SET(BINARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})
SET(EXECUTABLE_OUTPUT_PATH ${BINARY_OUTPUT_PATH}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${BINARY_OUTPUT_PATH}/modules)

#################################################################
# FFLAGS depend on the compiler and the build type
#################################################################

GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
#   SET(MACROS "${MACROS} -DDEBUG -Dmemcheck")
    ADD_DEFINITIONS(-DDEBUG)
    ADD_DEFINITIONS(-Dmemcheck)
ENDIF()

FIND_PACKAGE( OpenMP )
IF(${OPENMP_FOUND})
  SET(EXTRA_FLAGS ${OpenMP_Fortran_FLAGS})
ELSE()
  SET(EXTRA_FLAGS "")
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
# Find Fempar
#################################################################

FIND_PACKAGE(FEMPAR REQUIRED)
INCLUDE_DIRECTORIES(${FEMPAR_INCLUDES})

#################################################################
# Drivers
#################################################################

IF(FEMPAR_DRIVER)

    SET(DRIVERS_PATH ${CMAKE_SOURCE_DIR})
    SET(FEMPAR_DRIVER ${FEMPAR_DRIVER} CACHE STRING "Name of a directory in ${DRIVERS_PATH} containing the driver")

    FILE(GLOB DRIVERS_SRC ${DRIVERS_PATH}/${FEMPAR_DRIVER}/*.f90 ${DRIVERS_PATH}/${FEMPAR_DRIVER}/*.F90)
    SET(DRIVERS_SRC ${DRIVERS_SRC})

    IF(NOT DRIVERS_SRC)
        MESSAGE(FATAL_ERROR "Fempar Driver '${FEMPAR_DRIVER}' not found!")
    ENDIF()

    FOREACH(DRIVER_SRC ${DRIVERS_SRC})
        GET_FILENAME_COMPONENT(EXE_NAME ${DRIVER_SRC} NAME_WE)
        UNSET(${EXE_NAME}_MODULES_SRC CACHE)
        SET(${EXE_NAME}_MODULES_DIR ${DRIVERS_PATH}/${FEMPAR_DRIVER}/modules)
        IF(EXISTS "${${EXE_NAME}_MODULES_DIR}" AND IS_DIRECTORY "${${EXE_NAME}_MODULES_DIR}")
            FILE(GLOB ${EXE_NAME}_MODULES_SRC ${${EXE_NAME}_MODULES_DIR}/*.f90 ${${EXE_NAME}_MODULES_DIR}/*.F90)
        ENDIF()
        ADD_EXECUTABLE(${EXE_NAME} ${DRIVER_SRC} ${${EXE_NAME}_MODULES_SRC})
        TARGET_LINK_LIBRARIES(${EXE_NAME} ${FEMPAR_LIBS})
    ENDFOREACH()

ELSE()
    MESSAGE(FATAL_ERROR "Fempar Driver name must be passed to cmake with '-DFEMPAR_DRIVER=fempar_driver_name'")
ENDIF()



