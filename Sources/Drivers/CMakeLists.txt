#################################################################
# Executable targets
#################################################################
#SET(EXAMPLE_DRIVERS_PATH ${SRC_PATH}/Drivers/examples)
#ADD_EXECUTABLE(part_new ${EXAMPLE_DRIVERS_PATH}/part_new.f90)
#TARGET_LINK_LIBRARIES(part_new ${FEM_LIB})



#################################################################
# Serial drivers
#################################################################

FILE(GLOB_RECURSE EXAMPLE_DRIVERS_FEM_SRC Fem/*.f90)

FOREACH(DRIVER_SRC ${EXAMPLE_DRIVERS_FEM_SRC})
    GET_FILENAME_COMPONENT(EXE_NAME ${DRIVER_SRC} NAME_WE)
    ADD_EXECUTABLE(${EXE_NAME} ${DRIVER_SRC})
    TARGET_LINK_LIBRARIES(${EXE_NAME} ${FEM_LIB})
#    TARGET_LINK_LIBRARIES(${EXE_NAME} ${PAR_LIB})
    FOREACH (EXT_LIB ${EXT_LIBS})
        IF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${EXT_LIB}_FOUND)
            TARGET_LINK_LIBRARIES(${EXE_NAME} ${${EXT_LIB}_LIBRARIES})
#            GET_FILENAME_COMPONENT(LIB_PATH ${${EXT_LIB}_LIBRARIES} DIRECTORY)
#            MESSAGE("LIB_PATH: ${LIB_PATH}")
#            LINK_DIRECTORIES(${LIB_PATH})
        ENDIF()
    ENDFOREACH()
ENDFOREACH()


#################################################################
# Parallel drivers
#################################################################

FILE(GLOB_RECURSE EXAMPLE_DRIVERS_PAR_SRC Par/*.f90)

FOREACH(DRIVER_SRC ${EXAMPLE_DRIVERS_PAR_SRC})
    GET_FILENAME_COMPONENT(EXE_NAME ${DRIVER_SRC} NAME_WE)
    ADD_EXECUTABLE(${EXE_NAME} ${DRIVER_SRC})
    TARGET_LINK_LIBRARIES(${EXE_NAME} ${FEM_LIB})
    TARGET_LINK_LIBRARIES(${EXE_NAME} ${PAR_LIB})
    FOREACH (EXT_LIB ${EXT_LIBS})
        IF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${EXT_LIB}_FOUND)
            TARGET_LINK_LIBRARIES(${EXE_NAME} ${${EXT_LIB}_LIBRARIES})
#            GET_FILENAME_COMPONENT(LIB_PATH ${${EXT_LIB}_LIBRARIES} DIRECTORY)
#            MESSAGE("LIB_PATH: ${LIB_PATH}")
#            LINK_DIRECTORIES(${LIB_PATH})
        ENDIF()
    ENDFOREACH()
ENDFOREACH()
