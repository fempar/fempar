# CMake script to detect M library: 
#
# Usage example:
#   find_package(M)
#   if (M_FOUND)
#      include_directories(${M_INCLUDE_DIRS})
#      link_directories(${M_LIBRARY_DIRS})
#      add_executable(foo foo.c)
#      target_link_libraries(foo ${M_LIBRARIES})
#   endif()
#
# On return this will define:
#   M_FOUND                   Indicates whether M was found (True/False)
#   M_INCLUDE_DIRS            M include folder
#   M_LIBRARIES               M libraries names
#
####################################################################################

SET(DEFAULT_INC_PATHS /include /usr/include /usr/local/include)
SET(DEFAULT_LIB_PATHS /lib /lib32 /lib64 /usr/lib /usr/lib32 /usr/lib64 /usr/local/lib)

SET(DEP_PREFIX M)
SET(${DEP_PREFIX}_DEP_NAME m)
SET(${DEP_PREFIX}_INC_FILENAME math.h)
SET(${DEP_PREFIX}_INC_SUFFIXES "include")
SET(${DEP_PREFIX}_DEP_SUFFIXES "lib" "lib/${CMAKE_LIBRARY_ARCHITECTURE}" "${CMAKE_LIBRARY_ARCHITECTURE}")

SET(${DEP_PREFIX}_FOUND FALSE)

FIND_PATH(${DEP_PREFIX}_INCLUDES 
    NAMES ${${DEP_PREFIX}_INC_FILENAME}
    PATHS ${DEFAULT_INC_PATHS} ENV PATH
    PATH_SUFFIXES ${${DEP_PREFIX}_INC_SUFFIXES}
    DOC "The path to ${DEP_PREFIX} header files"
)

FIND_LIBRARY(${DEP_PREFIX}_LIBS
    NAMES ${${DEP_PREFIX}_DEP_NAME}
    PATHS ${DEFAULT_LIB_PATHS}  ENV PATH LD_LIBRARY_PATH
    PATH_SUFFIXES ${${DEP_PREFIX}_DEP_SUFFIXES}
    DOC "The path to ${DEP_PREFIX} library"
)

IF(${DEP_PREFIX}_LIBS)
    SET(${DEP_PREFIX}_INCLUDE_DIR ${${DEP_PREFIX}_INCLUDES})
    SET(${DEP_PREFIX}_LIBRARIES ${${DEP_PREFIX}_LIBS})
    SET(${DEP_PREFIX}_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(${DEP_PREFIX}_INCLUDES)
MARK_AS_ADVANCED(${DEP_PREFIX}_LIBS)
