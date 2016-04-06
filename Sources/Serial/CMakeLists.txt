#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE SERIAL_SRC *.f90 *.F90 *.c)
SET(SERIAL_SRC ${SERIAL_SRC} PARENT_SCOPE)


