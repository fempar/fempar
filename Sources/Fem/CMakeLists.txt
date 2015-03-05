#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE FEM_SRC *.f90)

#################################################################
# Fem library target
#################################################################

ADD_LIBRARY(${FEM_LIB} STATIC ${FEM_SRC})
SET_TARGET_PROPERTIES(${FEM_LIB} PROPERTIES VERSION ${FEMPAR_VERSION} SOVERSION ${FEMPAR_SOVERSION})

