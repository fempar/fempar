# CMake script to detect GIDPOST library: 
#
# Usage example:
#   find_package(GIDPOST)
#   if (GIDPOST_FOUND)
#      include_directories(${GIDPOST_INCLUDE_DIRS})
#      link_directories(${GIDPOST_LIBRARY_DIRS})
#      add_executable(foo foo.cc)
#      target_link_libraries(foo ${GIDPOST_LIBRARIES})
#   endif()
#
# On return this will define:
#   GIDPOST_FOUND                   Indicates whether GIDPOST was found (True/False)
#   GIDPOST_INCLUDE_DIRS            GIDPOST include folder
#   GIDPOST_LIBRARIES               GIDPOST libraries names
#
####################################################################################

SET(DEFAULT_INC_PATHS /include /usr/include /usr/local/include)
SET(DEFAULT_LIB_PATHS /lib /lib32 /lib64 /usr/lib /usr/lib32 /usr/lib64 /usr/local/lib)

SET(LIB_PREFIX GIDPOST)
SET(${LIB_PREFIX}_LIB_NAME gidpost)
SET(${LIB_PREFIX}_INC_FILENAME gidpost.h)
SET(${LIB_PREFIX}_INC_SUFFIXES "include" "gidpost" "include/gidpost")
SET(${LIB_PREFIX}_LIB_SUFFIXES "lib" "gidpost" "lib/gidpost" "lib/${CMAKE_LIBRARY_ARCHITECTURE}" "${CMAKE_LIBRARY_ARCHITECTURE}")

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
