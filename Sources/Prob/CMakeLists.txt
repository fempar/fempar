#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE PROB_SRC *.f90 *.F90)
set(PROB_SRC ${PROB_SRC} PARENT_SCOPE)

#################################################################
# Prob library target
#################################################################
ADD_LIBRARY(${PROB_LIB} ${PROB_SRC})
TARGET_LINK_LIBRARIES(${PROB_LIB} ${FEM_LIB})

FOREACH(EXT_LIB ${EXT_LIBS})
    IF (${EXT_LIB}_FOUND)
        TARGET_LINK_LIBRARIES(${PROB_LIB} ${${EXT_LIB}_LIBRARIES})
    ENDIF()
ENDFOREACH()

SET_TARGET_PROPERTIES(${PROB_LIB} PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})

