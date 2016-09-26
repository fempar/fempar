#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE LIBRARY_SRC *.f90 *.F90 *.c)
SET(LIBRARY_SRC ${LIBRARY_SRC} PARENT_SCOPE)

#################################################################
# ADD FEMPAR LIBRARY TARGET
#################################################################

ADD_LIBRARY(${PROJECT_NAME} ${LIBRARY_SRC})

FOREACH(EXT_LIB ${EXTERNAL_LIBRARIES})
    IF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${EXT_LIB}_FOUND)
        TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${${EXT_LIB}_LIBRARIES})
    ENDIF()
ENDFOREACH()

FOREACH(EXT_PROJ ${EXTERNAL_PROJECTS})
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${${EXT_PROJ}_LIBS})
ENDFOREACH()
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})
EXPORT(TARGETS ${PROJECT_NAME} APPEND FILE ${BINARY_OUTPUT_PATH}/${PROJECT_NAME}Targets.cmake)
