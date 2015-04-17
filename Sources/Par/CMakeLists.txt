FIND_PACKAGE(MPI)

IF(MPI_Fortran_COMPILER)

    INCLUDE(CMakeForceCompiler)
    CMAKE_FORCE_Fortran_COMPILER(${MPI_Fortran_COMPILER} ${CMAKE_Fortran_COMPILER_ID})


    ADD_DEFINITIONS(${MPI_Fortran_COMPILE_FLAGS})
    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
    LINK_DIRECTORIES(${MPI_Fortran_LIBRARIES})

#################################################################
# Search F90 files recursively in all subdirs
#################################################################

    FILE(GLOB_RECURSE PAR_SRC *.f90)
    SET(PAR_SRC ${PAR_SRC} PARENT_SCOPE)


#################################################################
# Par library target
#################################################################

    ADD_LIBRARY(${PAR_LIB} STATIC ${PAR_SRC})
    TARGET_LINK_LIBRARIES(${PAR_LIB} ${FEM_LIB})
    TARGET_LINK_LIBRARIES(${PAR_LIB} ${MPI_Fortran_LIBRARIES})
    FOREACH(EXT_LIB ${EXT_LIBS})
        IF(${EXT_LIB}_FOUND)
            TARGET_LINK_LIBRARIES(${PAR_LIB} ${${EXT_LIB}_LIBRARIES})
        ENDIF()
    ENDFOREACH()
    #SET_TARGET_PROPERTIES(${PAR_LIB} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS} -DMPI_MOD")
    ADD_DEFINITIONS(-DMPI_MOD)
    SET_TARGET_PROPERTIES(${PAR_LIB} PROPERTIES VERSION ${FEMPAR_VERSION} SOVERSION ${FEMPAR_SOVERSION})

ENDIF()

