# CMake script to detect P4EST library: 
#
# Usage example:
#   find_package(P4EST)
#   if (P4EST_FOUND)
#      include_directories(${P4EST_INCLUDE_DIRS})
#      link_directories(${P4EST_LIBRARY_DIRS})
#      add_executable(foo foo.cc)
#      target_link_libraries(foo ${P4EST_LIBRARIES})
#   endif()
#
# On return this will define:
#   P4EST_FOUND                   Indicates whether P4EST was found (True/False)
#   P4EST_INCLUDE_DIRS            P4EST include folder
#   P4EST_LIBRARIES               P4EST libraries names
#
####################################################################################

SET(DEFAULT_INC_PATHS /include /usr/include /usr/local/include)
SET(DEFAULT_LIB_PATHS /lib /lib32 /lib64 /usr/lib /usr/lib32 /usr/lib64 /usr/local/lib)

SET(LIB_PREFIX P4EST)
SET(${LIB_PREFIX}_LIB_NAME p4est)
SET(${LIB_PREFIX}_INC_FILENAME p4est.h)
SET(${LIB_PREFIX}_INC_SUFFIXES "include" "p4est" "include/p4est")
SET(${LIB_PREFIX}_LIB_SUFFIXES "lib" "p4est" "lib/p4est" "lib/${CMAKE_LIBRARY_ARCHITECTURE}" "${CMAKE_LIBRARY_ARCHITECTURE}")

SET(${LIB_PREFIX}_FOUND FALSE)

FIND_PATH(${LIB_PREFIX}_INCLUDES 
    NAMES ${${LIB_PREFIX}_INC_FILENAME}
    PATHS ${DEFAULT_INC_PATHS} ENV PATH
    PATH_SUFFIXES ${${LIB_PREFIX}_INC_SUFFIXES}
    DOC "The path to ${LIB_PREFIX} header files"
)

FIND_LIBRARY(${LIB_PREFIX}_LIBS
    NAMES ${${LIB_PREFIX}_LIB_NAME}
    PATHS ${DEFAULT_LIB_PATHS}  ENV PATH LD_LIBRARY_PATH
    PATH_SUFFIXES ${${LIB_PREFIX}_LIB_SUFFIXES}
    DOC "The path to ${LIB_PREFIX} library"
)

IF(${LIB_PREFIX}_LIBS)
    SET(${LIB_PREFIX}_INCLUDE_DIR ${${LIB_PREFIX}_INCLUDES})
    SET(${LIB_PREFIX}_LIBRARIES ${${LIB_PREFIX}_LIBS})
    SET(${LIB_PREFIX}_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(${LIB_PREFIX}_INCLUDES)
MARK_AS_ADVANCED(${LIB_PREFIX}_LIBS)
