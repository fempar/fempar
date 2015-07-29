#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE SERIAL_SRC *.f90 *.F90 *.c)
SET(SERIAL_SRC ${SERIAL_SRC} PARENT_SCOPE)

#################################################################
# Fem library target
#################################################################
ADD_LIBRARY(${SERIAL_LIB} ${SERIAL_SRC})

FOREACH(EXT_LIB ${EXT_LIBS})
    IF (${EXT_LIB}_FOUND)
        TARGET_LINK_LIBRARIES(${SERIAL_LIB} ${${EXT_LIB}_LIBRARIES})
    ENDIF()
ENDFOREACH()

SET_TARGET_PROPERTIES(${SERIAL_LIB} PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})

