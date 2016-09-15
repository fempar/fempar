#################################################################
# Drivers
#################################################################

IF(MPI_FOUND)

    IF(FEMPAR_DRIVER)

        FILE(GLOB DRIVERS_SRC ${DRIVERS_PATH}/${FEMPAR_DRIVER}/*.f90 ${DRIVERS_PATH}/${FEMPAR_DRIVER}/*.F90)
        SET(DRIVERS_SRC ${DRIVERS_SRC} PARENT_SCOPE)

        FOREACH(DRIVER_SRC ${DRIVERS_SRC})
            GET_FILENAME_COMPONENT(EXE_NAME ${DRIVER_SRC} NAME_WE)
            UNSET(${EXE_NAME}_MODULES_SRC CACHE)
            SET(${EXE_NAME}_MODULES_DIR ${DRIVERS_PATH}/${FEMPAR_DRIVER}/modules)
            IF(EXISTS "${${EXE_NAME}_MODULES_DIR}" AND IS_DIRECTORY "${${EXE_NAME}_MODULES_DIR}")
                FILE(GLOB ${EXE_NAME}_MODULES_SRC ${${EXE_NAME}_MODULES_DIR}/*.f90 ${${EXE_NAME}_MODULES_DIR}/*.F90)
            ENDIF()
            ADD_EXECUTABLE(${EXE_NAME} ${DRIVER_SRC} ${${EXE_NAME}_MODULES_SRC})
            TARGET_LINK_LIBRARIES(${EXE_NAME} ${EXT_PROJ_LIBS})
            TARGET_LINK_LIBRARIES(${EXE_NAME} ${MPI_Fortran_LIBRARIES})
            FOREACH (EXT_LIB ${EXT_LIBS})
                IF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${EXT_LIB}_FOUND)
                    TARGET_LINK_LIBRARIES(${EXE_NAME} ${${EXT_LIB}_LIBRARIES})
                ENDIF()
            ENDFOREACH()
        ENDFOREACH()

    ENDIF()

ENDIF()



