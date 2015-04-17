#################################################################
# Executable targets
#################################################################
#SET(EXAMPLE_DRIVERS_PATH ${SRC_PATH}/Drivers/examples)
#ADD_EXECUTABLE(part_new ${EXAMPLE_DRIVERS_PATH}/part_new.f90)
#TARGET_LINK_LIBRARIES(part_new ${FEM_LIB})



#################################################################
# Parallel drivers
#################################################################

IF(MPI_Fortran_COMPILER)


    IF(${PROJECT_NAME}_USER)

        INCLUDE(CMakeForceCompiler)
        CMAKE_FORCE_Fortran_COMPILER(${MPI_Fortran_COMPILER} ${CMAKE_Fortran_COMPILER_ID})

        FILE(GLOB_RECURSE DRIVERS_PAR_SRC ${${PROJECT_NAME}_USER}/*.f90)
        SET(DRIVERS_PAR_SRC ${DRIVERS_PAR_SRC} PARENT_SCOPE)

        FOREACH(DRIVER_SRC ${DRIVERS_PAR_SRC})
            GET_FILENAME_COMPONENT(EXE_NAME ${DRIVER_SRC} NAME_WE)
            ADD_EXECUTABLE(${EXE_NAME} ${DRIVER_SRC})
            TARGET_LINK_LIBRARIES(${EXE_NAME} ${FEM_LIB})
            TARGET_LINK_LIBRARIES(${EXE_NAME} ${PAR_LIB})
#           TARGET_LINK_LIBRARIES(${EXE_NAME} ${MPI_Fortran_LIBRARIES})
            FOREACH (EXT_LIB ${EXT_LIBS})
                IF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${EXT_LIB}_FOUND)
                    TARGET_LINK_LIBRARIES(${EXE_NAME} ${${EXT_LIB}_LIBRARIES})
#                   GET_FILENAME_COMPONENT(LIB_PATH ${${EXT_LIB}_LIBRARIES} DIRECTORY)
#                   MESSAGE("LIB_PATH: ${LIB_PATH}")
#                   LINK_DIRECTORIES(${LIB_PATH})
                ENDIF()
            ENDFOREACH()
        ENDFOREACH()
    ENDIF()
ENDIF()
