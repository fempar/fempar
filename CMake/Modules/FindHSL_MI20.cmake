# CMake script to detect HSL_MI20 library: 
#
# Usage example:
#   find_package(HSL_MI20)
#   if (HSL_MI20_FOUND)
#      include_directories(${HSL_MI20_INCLUDE_DIRS})
#      link_directories(${HSL_MI20_LIBRARY_DIRS})
#      add_executable(foo foo.cc)
#      target_link_libraries(foo ${HSL_MI20_LIBRARIES})
#   endif()
#
# On return this will define:
#   HSL_MI20_FOUND                   Indicates whether HSL_MI20 was found (True/False)
#   HSL_MI20_INCLUDE_DIRS            HSL_MI20 include folder
#   HSL_MI20_LIBRARIES               HSL_MI20 libraries names
#
####################################################################################

SET(DEFAULT_INC_PATHS /include /usr/include /usr/local/include)
SET(DEFAULT_LIB_PATHS /lib /lib32 /lib64 /usr/lib /usr/lib32 /usr/lib64 /usr/local/lib)

SET(LIB_PREFIX HSL_MI20)
SET(${LIB_PREFIX}_LIB_NAME hsl_mi20)
SET(${LIB_PREFIX}_INC_FILENAME hsl_mi20_double.mod hsl_mi20_single.mod)
SET(${LIB_PREFIX}_INC_SUFFIXES "include" "hsl_mi20" "include/hsl_mi20")
SET(${LIB_PREFIX}_LIB_SUFFIXES "lib" "hsl_mi20" "lib/hsl_mi20" "lib/${CMAKE_LIBRARY_ARCHITECTURE}" "${CMAKE_LIBRARY_ARCHITECTURE}")

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
