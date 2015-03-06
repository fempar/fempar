# CMake script to detect ZOLTAN library: 
#
# Usage example:
#   find_package(ZOLTAN)
#   if (ZOLTAN_FOUND)
#      include_directories(${ZOLTAN_INCLUDE_DIRS})
#      link_directories(${ZOLTAN_LIBRARY_DIRS})
#      add_executable(foo foo.cc)
#      target_link_libraries(foo ${ZOLTAN_LIBRARIES})
#   endif()
#
# On return this will define:
#   ZOLTAN_FOUND                   Indicates whether ZOLTAN was found (True/False)
#   ZOLTAN_INCLUDE_DIRS            ZOLTAN include folder
#   ZOLTAN_LIBRARIES               ZOLTAN libraries names
#
####################################################################################

SET(DEFAULT_INC_PATHS /include /usr/include /usr/local/include)
SET(DEFAULT_LIB_PATHS /lib /lib32 /lib64 /usr/lib /usr/lib32 /usr/lib64 /usr/local/lib)

SET(LIB_PREFIX ZOLTAN)
SET(${LIB_PREFIX}_LIB_NAME zoltan)
SET(${LIB_PREFIX}_INC_FILENAME zoltan.h)
SET(${LIB_PREFIX}_INC_SUFFIXES "include" "zoltan" "include/zoltan")
SET(${LIB_PREFIX}_LIB_SUFFIXES "lib" "zoltan" "lib/zoltan" "lib/${CMAKE_LIBRARY_ARCHITECTURE}" "${CMAKE_LIBRARY_ARCHITECTURE}")

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
