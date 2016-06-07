#################################################################
# Parallel drivers
#################################################################

IF(MPI_FOUND)

    IF(${PROJECT_NAME}_USER)

        FILE(GLOB_RECURSE DRIVERS_PAR_SRC ${${PROJECT_NAME}_USER}/*.f90  ${${PROJECT_NAME}_USER}/*.F90)
        SET(DRIVERS_PAR_SRC ${DRIVERS_PAR_SRC} PARENT_SCOPE)

        FOREACH(DRIVER_SRC ${DRIVERS_PAR_SRC})
            GET_FILENAME_COMPONENT(EXE_NAME ${DRIVER_SRC} NAME_WE)
            ADD_EXECUTABLE(${EXE_NAME} ${DRIVER_SRC})
            TARGET_LINK_LIBRARIES(${EXE_NAME} ${PROJECT_NAME})
            TARGET_LINK_LIBRARIES(${EXE_NAME} ${MPI_Fortran_LIBRARIES})
            FOREACH (EXT_LIB ${EXT_LIBS})
                IF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${EXT_LIB}_FOUND)
                    TARGET_LINK_LIBRARIES(${EXE_NAME} ${${EXT_LIB}_LIBRARIES})
                ENDIF()
            ENDFOREACH()
        ENDFOREACH()

    ENDIF()

ENDIF()
