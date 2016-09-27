#################################################################
# HEADER
#################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

PROJECT(FEMPAR Fortran C)

SET(${PROJECT_NAME}_VERSION 0.0.1)
SET(${PROJECT_NAME}_SOVERSION 1)

#################################################################
# DEFINE PATHS
#################################################################

SET(SRC_PATH ${CMAKE_SOURCE_DIR})
SET(FEMPAR_SRC_PATH ${SRC_PATH}/Sources)
SET(FEMPAR_TESTS_PATH ${FEMPAR_SRC_PATH}/Tests)
SET(FEMPAR_DRIVERS_PATH ${SRC_PATH}/Drivers)
SET(FEMPAR_THIRDPARTY_PATH ${SRC_PATH}/ThirdParty)

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
# SET BUILD NAME
#################################################################

FIND_PACKAGE( Git )
IF(${GIT_FOUND})
    # Get the current working branch
    EXECUTE_PROCESS(
        COMMAND git rev-parse --short HEAD
        WORKING_DIRECTORY ${SRC_PATH}
        OUTPUT_VARIABLE GIT_COMMIT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
ENDIF()

SET(BUILDNAME ${CMAKE_Fortran_COMPILER_ID}_${CMAKE_BUILD_TYPE}_MKL=${${PROJECT_NAME}_ENABLE_MKL}@${GIT_COMMIT} CACHE STRING "" )

#################################################################
# STATIC LIBRARIES
#################################################################
# Try to search first static libraries
IF(NOT ${BUILD_SHARED_LIBS})
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a ${CMAKE_FIND_LIBRARY_SUFFIXES}")
ENDIF()

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
                    ${FEMPAR_TESTS_PATH}/try_compile_mpi_mod.f90	
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
    SET(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARIES} ${MPI_LIBRARIES})
ENDIF()

#################################################################
# Add external projects
#################################################################

SET(EXT_PROJS FPL F3PHT FLAP VTKIO SISL XH5For FEMPAR)
include(ExternalProject)

FOREACH(EXT_PROJ ${EXT_PROJS})

    SET(${EXT_PROJ}_BINARY_PATH ${PROJECT_BINARY_DIR}/${EXT_PROJ})

    
    IF(${EXT_PROJ} STREQUAL "FEMPAR")
        EXTERNALPROJECT_Add(${EXT_PROJ}
            PREFIX ${EXT_PROJ}
            DEPENDS FPL F3PHT FLAP VTKIO SISL XH5For
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${${EXT_PROJ}_SRC_PATH}
            BINARY_DIR ${${EXT_PROJ}_BINARY_PATH}
            INSTALL_DIR ${INSTALL_PREFIX}
            CONFIGURE_COMMAND cmake -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE} -DMPI_Fortran_COMPILER=${MPI_Fortran_COMPILER} -DMPI_Fortran_COMPILER=${MPI_Fortran_COMPILER} -DMPI_C_COMPILER=${MPI_C_COMPILER} -DMPIEXEC=${MPIEXEC} -DFPL_DIR=${FPL_DIR} -DF3PHT_DIR=${F3PHT_DIR} -DFLAP_DIR=${FLAP_DIR} -DVTKIO_DIR=${VTKIO_DIR} -DSISL_DIR=${SISL_DIR} -DXH5For_DIR=${XH5For_DIR} -DFEMPAR_ENABLE_TESTS=${${PROJECT_NAME}_ENABLE_TESTS} ${${EXT_PROJ}_SRC_PATH}
            BUILD_COMMAND ${CMAKE_COMMAND} --build .
            TEST_COMMAND ""
            INSTALL_COMMAND ""
        )
    ELSEIF(${EXT_PROJ} STREQUAL "XH5For")
        SET(${EXT_PROJ}_SRC_PATH ${SRC_PATH}/ThirdParty/${EXT_PROJ})
        EXTERNALPROJECT_Add(${EXT_PROJ}
            PREFIX ${EXT_PROJ}
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${${EXT_PROJ}_SRC_PATH}
            BINARY_DIR ${${EXT_PROJ}_BINARY_PATH}
            INSTALL_DIR ${INSTALL_PREFIX}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DMPI_Fortran_COMPILER=${MPI_Fortran_COMPILER} -DMPI_C_COMPILER=${MPI_C_COMPILER} -DMPIEXEC=${MPIEXEC} -D${EXT_PROJ}_ENABLE_HDF5=${${PROJECT_NAME}_ENABLE_HDF5} -D${EXT_PROJ}_ENABLE_TESTS=${${PROJECT_NAME}_ENABLE_TESTS} ${${EXT_PROJ}_SRC_PATH}
            BUILD_COMMAND ${CMAKE_COMMAND} --build .
            TEST_COMMAND ""
            INSTALL_COMMAND ""
        )
    ELSEIF(${EXT_PROJ} STREQUAL "SISL")
        SET(${EXT_PROJ}_SRC_PATH ${SRC_PATH}/ThirdParty/${EXT_PROJ})
        EXTERNALPROJECT_Add(${EXT_PROJ}
            PREFIX ${EXT_PROJ}
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${${EXT_PROJ}_SRC_PATH}
            BINARY_DIR ${${EXT_PROJ}_BINARY_PATH}
            INSTALL_DIR ${INSTALL_PREFIX}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ${${EXT_PROJ}_SRC_PATH}
            BUILD_COMMAND ${CMAKE_COMMAND} --build .
            TEST_COMMAND ""
            INSTALL_COMMAND ""
        )
    ELSE()
        SET(${EXT_PROJ}_SRC_PATH ${SRC_PATH}/ThirdParty/${EXT_PROJ})
        EXTERNALPROJECT_Add(${EXT_PROJ}
            PREFIX ${EXT_PROJ}
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${${EXT_PROJ}_SRC_PATH}
            BINARY_DIR ${${EXT_PROJ}_BINARY_PATH}
            INSTALL_DIR ${INSTALL_PREFIX}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -D${EXT_PROJ}_ENABLE_TESTS=${${PROJECT_NAME}_ENABLE_TESTS} ${${EXT_PROJ}_SRC_PATH}
            BUILD_COMMAND ${CMAKE_COMMAND} --build .
            TEST_COMMAND ""
            INSTALL_COMMAND ""
        )
    ENDIF()


    EXTERNALPROJECT_Add_Step(${EXT_PROJ} preconfigure
        COMMAND ${CMAKE_COMMAND} -E echo "Force configure of ${EXT_PROJ}"
        DEPENDERS configure
        ALWAYS 1
    )

    IF(NOT ${EXT_PROJ} STREQUAL "SISL")
        EXTERNALPROJECT_Add_Step(${EXT_PROJ} customtest
            COMMAND ${CTEST_COMMAND}
            DEPENDERS test
            BYPRODUCTS ${${EXT_PROJ}_BINARY_PATH}/${EXT_PROJ}Config.cmake
            WORKING_DIRECTORY ${${EXT_PROJ}_BINARY_PATH}
            EXCLUDE_FROM_MAIN 1
        )
    ENDIF()

    SET(${EXT_PROJ}_DIR ${${EXT_PROJ}_BINARY_PATH})

ENDFOREACH()


