#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE FEM_SRC *.f90)
SET(FEM_SRC ${FEM_SRC} PARENT_SCOPE)

#################################################################
# Fem library target
#################################################################
ADD_LIBRARY(${FEM_LIB} STATIC ${FEM_SRC})

FOREACH(EXT_LIB ${EXT_LIBS})
    IF (${EXT_LIB}_FOUND)
        TARGET_LINK_LIBRARIES(${FEM_LIB} ${${EXT_LIB}_LIBRARIES})
    ENDIF()
ENDFOREACH()

SET_TARGET_PROPERTIES(${FEM_LIB} PROPERTIES VERSION ${FEMPAR_VERSION} SOVERSION ${FEMPAR_SOVERSION})

