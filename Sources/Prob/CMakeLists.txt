#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE PROB_SRC *.f90)

#################################################################
# Prob library target
#################################################################
ADD_LIBRARY(${PROB_LIB} STATIC ${PROB_SRC})
TARGET_LINK_LIBRARIES(${PROB_LIB} ${FEM_LIB})

FOREACH(EXT_LIB ${EXT_LIBS})
    IF (${EXT_LIB}_FOUND)
        TARGET_LINK_LIBRARIES(${PROB_LIB} ${${EXT_LIB}_LIBRARIES})
    ENDIF()
ENDFOREACH()

SET_TARGET_PROPERTIES(${FEM_LIB} PROPERTIES VERSION ${FEMPAR_VERSION} SOVERSION ${FEMPAR_SOVERSION})

