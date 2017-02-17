#################################################################
# HEADER
#################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

PROJECT(FEMPAR_THIRDPARTY Fortran C)

SET(${PROJECT_NAME}_VERSION 0.0.1)
SET(${PROJECT_NAME}_SOVERSION 1)

#################################################################
# DEFINE PATHS
#################################################################

SET(ROOT_PATH ${CMAKE_SOURCE_DIR}) 

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

#################################################################
# ENABLE TESTING
#################################################################

IF(NOT ${PROJECT_NAME}_ENABLE_TESTS)         # Disable tests as default
    OPTION(${PROJECT_NAME}_ENABLE_TESTS "Enable/disable tests" OFF)
ENDIF(NOT ${PROJECT_NAME}_ENABLE_TESTS)

IF(${PROJECT_NAME}_ENABLE_TESTS)             # If ENABLE_TESTS flag was passed, enable CTest
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
        WORKING_DIRECTORY ${ROOT_PATH}
        OUTPUT_VARIABLE GIT_COMMIT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
ENDIF()

# Append the commit ID to the BUILDNAME if Git is found
SET(BUILDNAME ${CMAKE_Fortran_COMPILER_ID}_${CMAKE_BUILD_TYPE}_MKL=${${PROJECT_NAME}_ENABLE_MKL}@${GIT_COMMIT} CACHE STRING "" )

#################################################################
# STATIC LIBRARIES
#################################################################

# If shared libs not forced, prepend ".a" extension to library suffixes. 
# Useful to find first static libraries in future calls to FIND_PACKAGE
IF(NOT ${BUILD_SHARED_LIBS}) 
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a ${CMAKE_FIND_LIBRARY_SUFFIXES}")
ENDIF()

#################################################################
# Add external projects
#################################################################

include(ExternalProject)

# List of external projects to be managed from this CMakeLists.txt
SET(EXTERNAL_PROJECTS FPL F3PHT FLAP VTKIO SISL XH5For)

# Loop over the list of external projects
FOREACH(EXT_PROJ ${EXTERNAL_PROJECTS})

    SET(${EXT_PROJ}_BINARY_PATH ${PROJECT_BINARY_DIR}/${EXT_PROJ})
    SET(${EXT_PROJ}_SRC_PATH ${ROOT_PATH}/${EXT_PROJ})
    
    IF(${EXT_PROJ} STREQUAL "XH5For")
        # Add XH5For external project
        EXTERNALPROJECT_Add(${EXT_PROJ}
            PREFIX ${EXT_PROJ}
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${${EXT_PROJ}_SRC_PATH}
            BINARY_DIR ${${EXT_PROJ}_BINARY_PATH}
            INSTALL_DIR ${INSTALL_PREFIX}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DMPI_Fortran_COMPILER=${MPI_Fortran_COMPILER} -DMPI_C_COMPILER=${MPI_C_COMPILER} -DMPIEXEC=${MPIEXEC} -D${EXT_PROJ}_ENABLE_TESTS=${${PROJECT_NAME}_ENABLE_TESTS} -D${EXT_PROJ}_ENABLE_EXAMPLES=${${PROJECT_NAME}_ENABLE_TESTS} ${${EXT_PROJ}_SRC_PATH}
            BUILD_COMMAND ${CMAKE_COMMAND} --build .
            TEST_COMMAND ${CTEST_COMMAND}
            INSTALL_COMMAND ""
        )
    ELSEIF(${EXT_PROJ} STREQUAL "SISL")
        # Add SISL external project
        EXTERNALPROJECT_Add(${EXT_PROJ}
            PREFIX ${EXT_PROJ}
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${${EXT_PROJ}_SRC_PATH}
            BINARY_DIR ${${EXT_PROJ}_BINARY_PATH}
            INSTALL_DIR ${INSTALL_PREFIX}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} ${${EXT_PROJ}_SRC_PATH}
            BUILD_COMMAND ${CMAKE_COMMAND} --build .
            TEST_COMMAND ""
            INSTALL_COMMAND ""
        )
    ELSEIF(${EXT_PROJ} STREQUAL "FPL" OR ${EXT_PROJ} STREQUAL "FLAP")
        # Add FPL external project
        EXTERNALPROJECT_Add(${EXT_PROJ}
            PREFIX ${EXT_PROJ}
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${${EXT_PROJ}_SRC_PATH}
            BINARY_DIR ${${EXT_PROJ}_BINARY_PATH}
            INSTALL_DIR ${INSTALL_PREFIX}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -D${EXT_PROJ}_ENABLE_OPENMP=${${PROJECT_NAME}_ENABLE_OPENMP} -D${EXT_PROJ}_ENABLE_TESTS=${${PROJECT_NAME}_ENABLE_TESTS} ${${EXT_PROJ}_SRC_PATH}
            BUILD_COMMAND ${CMAKE_COMMAND} --build .
            TEST_COMMAND ${CTEST_COMMAND}
            INSTALL_COMMAND ""
        )
    ELSE()
        # Add custom external project
        EXTERNALPROJECT_Add(${EXT_PROJ}
            PREFIX ${EXT_PROJ}
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${${EXT_PROJ}_SRC_PATH}
            BINARY_DIR ${${EXT_PROJ}_BINARY_PATH}
            INSTALL_DIR ${INSTALL_PREFIX}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -D${EXT_PROJ}_ENABLE_TESTS=${${PROJECT_NAME}_ENABLE_TESTS} ${${EXT_PROJ}_SRC_PATH}
            BUILD_COMMAND ${CMAKE_COMMAND} --build .
            TEST_COMMAND ${CTEST_COMMAND}
            INSTALL_COMMAND ""
        )
    ENDIF()

    SET(${EXT_PROJ}_DIR ${${EXT_PROJ}_BINARY_PATH})


    # Custom clean target for every external project
    ADD_CUSTOM_TARGET(clean-${EXT_PROJ} 
        COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target clean 
        COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH} --target clean 
    )


ENDFOREACH()


