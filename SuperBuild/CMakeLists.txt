#################################################################
# HEADER
#################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

PROJECT(SUPERBUILD Fortran C)

SET(${PROJECT_NAME}_VERSION 0.0.1)
SET(${PROJECT_NAME}_SOVERSION 1)

#################################################################
# DEFINE PATHS
#################################################################

SET(ROOT_PATH ${CMAKE_SOURCE_DIR}/..)                                        # Fempar root path 
SET(SRC_PATH ${ROOT_PATH}/Sources)                                           # Fempar Sources path 
SET(TESTS_PATH ${SRC_PATH}/Tests)                                            # Fempar Tests path 
SET(DRIVERS_PATH /Users/ocolomes/Progs/Fempar_GitLab/NSI-fempar/src/Drivers) # Fempar Drivers path 
SET(THIRDPARTY_PATH ${ROOT_PATH}/ThirdParty)                                 # Fempar ThirdParty path 

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
# STATIC LIBRARIES
#################################################################

# If shared libs not forced, prepend ".a" extension to library suffixes. 
# Useful to find first static libraries in future calls to FIND_PACKAGE
IF(NOT ${BUILD_SHARED_LIBS})
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a ${CMAKE_FIND_LIBRARY_SUFFIXES}")
ENDIF()

#################################################################
# Find REQUIRED MPI library and check MPI.mod compatibility
#################################################################

FIND_PACKAGE(MPI REQUIRED) # Look for MPI library. This requirement is mandatory required

IF(MPI_FOUND)
    # Add MPI directories
    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH}) # Add MPI includes to fempar includes

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
# Detect FEMPAR-related variables (ENABLES & External libraries)
#################################################################

# Set of fempar variable prefixes
SET(VARIABLE_PREFIXES FEMPAR_ENABLE "MPI[A-Za-z0-9_]*(EXEC|COMPILER|LIBRARIES|INCLUDES)" MKL BLAS LAPACK GIDPOST METIS WSMP HSL_MI20 HSL_MA87 P4EST UMFPACK ZOLTAN ZLIB)

# Look into CMakeCache for all variables starting with the prefixes defined in VARIABLE_PREFIXES and store them into FEMPAR_VARIABLES
# FEMPAR_VARIABLES is a list of variables that are going to be passed to FEMPAR external project configuration step
UNSET(FEMPAR_VARIABLES)
GET_CMAKE_PROPERTY(ALL_VARIABLES VARIABLES)
FOREACH(VARIABLE_PREFIX ${VARIABLE_PREFIXES})
    STRING (REGEX MATCHALL "(^|;)${VARIABLE_PREFIX}[A-Za-z0-9_]*" MATCHED_VARIABLES "${ALL_VARIABLES}")
    FOREACH (VARIABLE ${MATCHED_VARIABLES})
            SET(${VARIABLE} ${${VARIABLE}})
            SET(FEMPAR_VARIABLES ${FEMPAR_VARIABLES} "-D${VARIABLE}=${${VARIABLE}}")
    ENDFOREACH()
ENDFOREACH()

#################################################################
# Add external projects
#################################################################

SET(EXT_PROJS THIRDPARTY FEMPAR) # List of external projects to be managed by the SuperBuild
IF(FEMPAR_DRIVER)
    SET(EXT_PROJS ${EXT_PROJS} DRIVERS)
ENDIF()
include(ExternalProject)

# Iterate over the list of external projects and Add them to the SuperBuild
FOREACH(EXT_PROJ ${EXT_PROJS})

    SET(${EXT_PROJ}_BINARY_PATH ${PROJECT_BINARY_DIR}/${EXT_PROJ})                      # Define default binary path for external projects

    IF(${EXT_PROJ} STREQUAL "DRIVERS")
        SET(${EXT_PROJ}_SRC_PATH ${DRIVERS_PATH})                                       # Define Drivers source path (IN)
        SET(${EXT_PROJ}_BINARY_PATH ${PROJECT_BINARY_DIR}/${EXT_PROJ}/${FEMPAR_DRIVER}) # Define Drivers binary path (OUT)
        # Add FEMPAR_DRIVER external project
        EXTERNALPROJECT_Add(${EXT_PROJ}
            PREFIX ${EXT_PROJ}/${FEMPAR_DRIVER}
            STAMP_DIR ${EXT_PROJ}/${FEMPAR_DRIVER}/src/${FEMPAR_DRIVER}-stamp
            DEPENDS FEMPAR
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${${EXT_PROJ}_SRC_PATH}
            BINARY_DIR ${${EXT_PROJ}_BINARY_PATH}
            INSTALL_DIR ${INSTALL_PREFIX}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DFEMPAR_DIR=${FEMPAR_DIR} -DFEMPAR_DRIVER=${FEMPAR_DRIVER} ${${EXT_PROJ}_SRC_PATH}
            BUILD_COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH}
            TEST_COMMAND ""
            INSTALL_COMMAND ""
        )
        # Add Preconfigure step to force Configure-Build-Test-Install steps for Drivers (p.e., every time we call to $ make)
        EXTERNALPROJECT_Add_Step(${EXT_PROJ} PreConfigure
    	    COMMAND ${CMAKE_COMMAND} -E touch ${EXT_PROJ}/${FEMPAR_DRIVER}/src/${FEMPAR_DRIVER}-stamp/${FEMPAR_DRIVER}-PreConfigure-Fake.stamp
            DEPENDEES update
            DEPENDERS configure
            ALWAYS 1
        )
        # Add a custom target (DRIVERS-clean) to clean the current FEMPAR_DRIVER
        ADD_CUSTOM_TARGET(DRIVERS-clean
            COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target clean 
            COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH} --target clean 
        )
    ELSEIF(${EXT_PROJ} STREQUAL "FEMPAR")
        SET(${EXT_PROJ}_SRC_PATH ${ROOT_PATH}) # Define Fempar source path (IN)
        # Add FEMPAR external project
        EXTERNALPROJECT_Add(${EXT_PROJ}
            PREFIX ${EXT_PROJ}
            STAMP_DIR ${EXT_PROJ}/src/${EXT_PROJ}-stamp
            DEPENDS THIRDPARTY
            STEP_TARGETS configure build 
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${${EXT_PROJ}_SRC_PATH}
            BINARY_DIR ${${EXT_PROJ}_BINARY_PATH}
            INSTALL_DIR ${INSTALL_PREFIX}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE} ${FEMPAR_VARIABLES} -DFEMPAR_THIRDPARTY_DIR=${THIRDPARTY_DIR} ${${EXT_PROJ}_SRC_PATH}
            BUILD_COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH}
            TEST_COMMAND ""
            INSTALL_COMMAND ""
        )
        # Add Preconfigure step to force Configure-Build-Test-Install steps for Fempar (p.e., every time we call to $ make)
        EXTERNALPROJECT_Add_Step(${EXT_PROJ} PreConfigure
	        COMMAND ${CMAKE_COMMAND} -E touch ${EXT_PROJ}/src/${EXT_PROJ}-stamp/${EXT_PROJ}-PreConfigure-Fake.stamp
            DEPENDEES update
            DEPENDERS configure
            ALWAYS 1
        )
        IF(FEMPAR_ENABLE_TESTS)
            # Add a custom target (FEMPAR-test) to launch Fempar tests
            ADD_CUSTOM_TARGET(${EXT_PROJ}-test
                COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH} --target test
            )
            # Add a custom target (FEMPAR-clean-tests) to clean Fempar tests
            ADD_CUSTOM_TARGET(${EXT_PROJ}-clean-tests
                COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH} --target clean-tests
            )
        ENDIF()
        # Add a custom target (FEMPAR-clean) to clean Fempar
        ADD_CUSTOM_TARGET(${EXT_PROJ}-clean
            COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target clean 
            COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH} --target clean 
        )
    ELSEIF(${EXT_PROJ} STREQUAL "THIRDPARTY")
        SET(${EXT_PROJ}_SRC_PATH ${THIRDPARTY_PATH}) # Define ThirdParty source/root path (IN)
        # Add FEMPAR external project
        EXTERNALPROJECT_Add(${EXT_PROJ}
            PREFIX ${EXT_PROJ}
            STEP_TARGETS configure build
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${${EXT_PROJ}_SRC_PATH}
            BINARY_DIR ${${EXT_PROJ}_BINARY_PATH}
            INSTALL_DIR ${INSTALL_PREFIX}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DMPI_Fortran_COMPILER=${MPI_Fortran_COMPILER} -DMPI_C_COMPILER=${MPI_C_COMPILER} -DMPIEXEC=${MPIEXEC} -DFEMPAR_THIRDPARTY_ENABLE_OPENMP=${FEMPAR_ENABLE_OPENMP} -DFEMPAR_THIRDPARTY_ENABLE_TESTS=${FEMPAR_ENABLE_TESTS} ${${EXT_PROJ}_SRC_PATH}
            BUILD_COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH}
            TEST_COMMAND ""
            INSTALL_COMMAND ""
        )
        # Add a custom target (THIRDPARTY-clean) to clean ThirdParty libraries
        ADD_CUSTOM_TARGET(${EXT_PROJ}-clean
            COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target clean 
            COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH} --target clean 
        )
    ENDIF()

    # Set ${EXT_PROJ}_DIR (THIRDPARTY_DIR and FEMPAR_DIR) variables 
    # this variables are used to communicate the binary directories 
    # between the different external projects when needed.
    # (P.e., Fempar must know where ThirdParty libraries are)
    SET(${EXT_PROJ}_DIR ${${EXT_PROJ}_BINARY_PATH})

ENDFOREACH()


