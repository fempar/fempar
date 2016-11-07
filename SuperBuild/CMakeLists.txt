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

SET(SRC_PATH ${CMAKE_SOURCE_DIR}/..)
SET(ROOT_PATH ${SRC_PATH})
SET(SRC_PATH ${SRC_PATH}/Sources)
SET(TESTS_PATH ${SRC_PATH}/Tests)
SET(DRIVERS_PATH ${ROOT_PATH}/Drivers)
SET(THIRDPARTY_PATH ${ROOT_PATH}/ThirdParty)

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
    SET(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARIES} ${MPI_LIBRARIES})
ENDIF()


#################################################################
# Detect FEMPAR-related variables (ENABLES & External libraries)
#################################################################

SET(VARIABLE_PREFIXES FEMPAR_ENABLE "MPI[A-Za-z0-9_]*(EXEC|COMPILER|LIBRARIES|INCLUDES)" MKL BLAS LAPACK GIDPOST METIS WSMP HSL_MI20 HSL_MA87 P4EST UMFPACK ZOLTAN ZLIB)

UNSET(FEMPAR_VARIABLES)
GET_CMAKE_PROPERTY(ALL_VARIABLES VARIABLES)
FOREACH(VARIABLE_PREFIX ${VARIABLE_PREFIXES})
    STRING (REGEX MATCHALL "(^|;)${VARIABLE_PREFIX}[A-Za-z0-9_]*" MATCHED_VARIABLES "${ALL_VARIABLES}")
    FOREACH (VARIABLE ${MATCHED_VARIABLES})
            MESSAGE(${VARIABLE}=${${VARIABLE}})
            SET(${VARIABLE} ${${VARIABLE}})
            SET(FEMPAR_VARIABLES ${FEMPAR_VARIABLES} "-D${VARIABLE}=${${VARIABLE}}")
    ENDFOREACH()
ENDFOREACH()

#################################################################
# Add external projects
#################################################################

SET(EXT_PROJS THIRDPARTY FEMPAR)
IF(FEMPAR_DRIVER)
    SET(EXT_PROJS ${EXT_PROJS} DRIVERS)
ENDIF()
include(ExternalProject)

FOREACH(EXT_PROJ ${EXT_PROJS})

    SET(${EXT_PROJ}_BINARY_PATH ${PROJECT_BINARY_DIR}/${EXT_PROJ})

    IF(${EXT_PROJ} STREQUAL "DRIVERS")
        SET(${EXT_PROJ}_SRC_PATH ${DRIVERS_PATH})
        SET(${EXT_PROJ}_BINARY_PATH ${PROJECT_BINARY_DIR}/${EXT_PROJ}/${FEMPAR_DRIVER})
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
        EXTERNALPROJECT_Add_Step(${EXT_PROJ} PreConfigure
    	    COMMAND ${CMAKE_COMMAND} -E touch ${EXT_PROJ}/${FEMPAR_DRIVER}/src/${FEMPAR_DRIVER}-stamp/${FEMPAR_DRIVER}-PreConfigure-Fake.stamp
            DEPENDEES update
            DEPENDERS configure
            ALWAYS 1
        )
        ADD_CUSTOM_TARGET(DRIVERS-clean
            COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target clean 
            COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH} --target clean 
        )
    ELSEIF(${EXT_PROJ} STREQUAL "FEMPAR")
        SET(${EXT_PROJ}_SRC_PATH ${ROOT_PATH})
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
        EXTERNALPROJECT_Add_Step(${EXT_PROJ} force-build
	        COMMAND ${CMAKE_COMMAND} -E touch ${EXT_PROJ}/src/${EXT_PROJ}-stamp/${EXT_PROJ}-custom-force-build.stamp
            DEPENDEES update
            DEPENDERS configure
            ALWAYS 1
        )
        IF(FEMPAR_ENABLE_TESTS)
            ADD_CUSTOM_TARGET(${EXT_PROJ}-test
                COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH} --target test
            )
            ADD_CUSTOM_TARGET(${EXT_PROJ}-clean-tests
                COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH} --target clean-tests
            )
        ENDIF()
        ADD_CUSTOM_TARGET(${EXT_PROJ}-clean
            COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target clean 
            COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH} --target clean 
        )
    ELSEIF(${EXT_PROJ} STREQUAL "THIRDPARTY")
        SET(${EXT_PROJ}_SRC_PATH ${THIRDPARTY_PATH})
        EXTERNALPROJECT_Add(${EXT_PROJ}
            PREFIX ${EXT_PROJ}
            STEP_TARGETS configure build
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${${EXT_PROJ}_SRC_PATH}
            BINARY_DIR ${${EXT_PROJ}_BINARY_PATH}
            INSTALL_DIR ${INSTALL_PREFIX}
            CONFIGURE_COMMAND ${CMAKE_COMMAND} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DMPI_Fortran_COMPILER=${MPI_Fortran_COMPILER} -DMPI_C_COMPILER=${MPI_C_COMPILER} -DMPIEXEC=${MPIEXEC} -DFEMPAR_THIRDPARTY_ENABLE_TESTS=${FEMPAR_ENABLE_TESTS} ${${EXT_PROJ}_SRC_PATH}
            BUILD_COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH}
            TEST_COMMAND ""
            INSTALL_COMMAND ""
        )
        ADD_CUSTOM_TARGET(${EXT_PROJ}-clean
            COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target clean 
            COMMAND ${CMAKE_COMMAND} --build ${${EXT_PROJ}_BINARY_PATH} --target clean 
        )
    ENDIF()

    SET(${EXT_PROJ}_DIR ${${EXT_PROJ}_BINARY_PATH})

ENDFOREACH()


