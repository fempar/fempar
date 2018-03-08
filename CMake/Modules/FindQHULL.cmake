# CMake script to detect QHULL library: 
#
# Usage example:
#   find_package(QHULL)
#   if (QHULL_FOUND)
#      include_directories(${QHULL_INCLUDE_DIRS})
#      link_directories(${QHULL_LIBRARY_DIRS})
#      add_executable(foo foo.cc)
#      target_link_libraries(foo ${QHULL_LIBRARIES})
#   endif()
#
# On return this will define:
#   QHULL_FOUND                   Indicates whether QHULL was found (True/False)
#   QHULL_INCLUDE_DIRS            QHULL include folder
#   QHULL_LIBRARIES               QHULL libraries names
#
####################################################################################

SET(DEFAULT_INC_PATHS /include /usr/include /usr/local/include)
SET(DEFAULT_LIB_PATHS /lib /lib32 /lib64 /usr/lib /usr/lib32 /usr/lib64 /usr/local/lib)

SET(LIB_PREFIX QHULL)
SET(${LIB_PREFIX}_LIB_NAME qhullstatic)
SET(${LIB_PREFIX}_INC_FILENAME libqhull.h)
SET(${LIB_PREFIX}_INC_SUFFIXES "include" "libqhull" "include/libqhull")
SET(${LIB_PREFIX}_LIB_SUFFIXES "lib" "libqhull" "lib/libqhull" "lib/${CMAKE_LIBRARY_ARCHITECTURE}" "${CMAKE_LIBRARY_ARCHITECTURE}")

SET(${LIB_PREFIX}_FOUND FALSE)

FIND_PATH(${LIB_PREFIX}_INCLUDES 
    NAMES ${${LIB_PREFIX}_INC_FILENAME}
    PATHS ${DEFAULT_INC_PATHS} ENV PATH
    PATH_SUFFIXES ${${LIB_PREFIX}_INC_SUFFIXES}
    DOC "The path to ${LIB_PREFIX} header files"
)

FIND_LIBRARY(${LIB_PREFIX}_LIBS
    NAMES ${${LIB_PREFIX}_LIB_NAME}
    PATHS ${DEFAULT_LIB_PATHS}  ENV LD_LIBRARY_PATH
    PATH_SUFFIXES ${${LIB_PREFIX}_LIB_SUFFIXES}
    DOC "The path to ${LIB_PREFIX} library"
)

IF(${LIB_PREFIX}_INCLUDES)
  IF(${LIB_PREFIX}_LIBS)
    SET(${LIB_PREFIX}_FOUND TRUE)
    SET(${LIB_PREFIX}_INCLUDE_DIR ${${LIB_PREFIX}_INCLUDES})
    SET(${LIB_PREFIX}_LIBRARIES ${${LIB_PREFIX}_LIBS})
  ENDIF()
ENDIF()

MESSAGE(STATUS "QHULL_FOUND:" ${${LIB_PREFIX}_FOUND})
IF(${${LIB_PREFIX}_FOUND})
  MESSAGE(STATUS "QHULL_INCLUDES:" ${${LIB_PREFIX}_INCLUDE_DIR})
  MESSAGE(STATUS "QHULL_LIBS:" ${${LIB_PREFIX}_LIBRARIES})
ENDIF()

MARK_AS_ADVANCED(${LIB_PREFIX}_INCLUDES)
MARK_AS_ADVANCED(${LIB_PREFIX}_LIBS)
