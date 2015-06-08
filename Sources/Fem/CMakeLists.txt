#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE FEM_SRC *.f90 *.F90)
SET(FEM_SRC ${FEM_SRC} PARENT_SCOPE)

#################################################################
# Fem library target
#################################################################
ADD_LIBRARY(${FEM_LIB} ${FEM_SRC})

FOREACH(EXT_LIB ${EXT_LIBS})
    IF (${EXT_LIB}_FOUND)
        TARGET_LINK_LIBRARIES(${FEM_LIB} ${${EXT_LIB}_LIBRARIES})
    ENDIF()
ENDFOREACH()

SET_TARGET_PROPERTIES(${FEM_LIB} PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})

