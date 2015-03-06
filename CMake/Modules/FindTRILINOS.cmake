# CMake script to detect TRILINOS library: 
#
# Usage example:
#   find_package(TRILINOS)
#   if (TRILINOS_FOUND)
#      include_directories(${TRILINOS_INCLUDE_DIRS})
#      link_directories(${TRILINOS_LIBRARY_DIRS})
#      add_executable(foo foo.cc)
#      target_link_libraries(foo ${TRILINOS_LIBRARIES})
#   endif()
#
# On return this will define:
#   TRILINOS_FOUND                   Indicates whether TRILINOS was found (True/False)
#   TRILINOS_INCLUDE_DIRS            TRILINOS include folder
#   TRILINOS_LIBRARIES               TRILINOS libraries names
#
####################################################################################

SET(DEFAULT_INC_PATHS /include /usr/include /usr/local/include)
SET(DEFAULT_LIB_PATHS /lib /lib32 /lib64 /usr/lib /usr/lib32 /usr/lib64 /usr/local/lib)

SET(LIB_PREFIX TRILINOS)
SET(${LIB_PREFIX}_LIB_NAME trilinos)
SET(${LIB_PREFIX}_INC_FILENAME trilinos.h)
SET(${LIB_PREFIX}_INC_SUFFIXES "include" "trilinos" "include/trilinos")
SET(${LIB_PREFIX}_LIB_SUFFIXES "lib" "trilinos" "lib/trilinos" "lib/${CMAKE_LIBRARY_ARCHITECTURE}" "${CMAKE_LIBRARY_ARCHITECTURE}")

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
