#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE LIB_SRC *.f90 *.F90 *.c)
SET(LIB_SRC ${LIB_SRC} PARENT_SCOPE)

#################################################################
# Library target
#################################################################
ADD_LIBRARY(${LIB} ${LIB_SRC})

#################################################################
# External libraries
#################################################################
FOREACH(EXT_LIB ${EXT_LIBS})
    IF (${EXT_LIB}_FOUND)
        TARGET_LINK_LIBRARIES(${LIB} ${${EXT_LIB}_LIBRARIES})
    ENDIF()
ENDFOREACH()

SET_TARGET_PROPERTIES(${LIB} PROPERTIES VERSION ${${LIB}_VERSION} SOVERSION ${${LIB}_SOVERSION})
EXPORT(TARGETS ${LIB} APPEND FILE ${PROJECT_BINARY_DIR}/${LIB}Targets.cmake)

