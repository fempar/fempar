#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE LIBRARY_SRC *.f90 *.F90 *.c)
SET(LIBRARY_SRC ${LIBRARY_SRC} PARENT_SCOPE)

#################################################################
# ADD FEMPAR LIBRARY TARGET
#################################################################

ADD_LIBRARY(${PROJECT_NAME} ${LIBRARY_SRC})
ADD_DEPENDENCIES(${PROJECT_NAME} ${EXT_PROJS})
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})


