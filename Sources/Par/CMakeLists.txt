#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE PAR_SRC *.f90 *.F90 *.c)
SET(PAR_SRC ${PAR_SRC} PARENT_SCOPE)


