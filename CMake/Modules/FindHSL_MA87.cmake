# CMake script to detect HSL_MA87 library: 
#
# Usage example:
#   find_package(HSL_MA87)
#   if (HSL_MA87_FOUND)
#      include_directories(${HSL_MA87_INCLUDE_DIRS})
#      link_directories(${HSL_MA87_LIBRARY_DIRS})
#      add_executable(foo foo.cc)
#      target_link_libraries(foo ${HSL_MA87_LIBRARIES})
#   endif()
#
# On return this will define:
#   HSL_MA87_FOUND                   Indicates whether HSL_MA87 was found (True/False)
#   HSL_MA87_INCLUDE_DIRS            HSL_MA87 include folder
#   HSL_MA87_LIBRARIES               HSL_MA87 libraries names
#
####################################################################################

SET(DEFAULT_INC_PATHS /include /usr/include /usr/local/include)
SET(DEFAULT_LIB_PATHS /lib /lib32 /lib64 /usr/lib /usr/lib32 /usr/lib64 /usr/local/lib)

SET(LIB_PREFIX HSL_MA87)
SET(${LIB_PREFIX}_LIB_NAME hsl_ma87)
SET(${LIB_PREFIX}_INC_FILENAME hsl_ma87_double.mod hsl_ma87_single.mod hsl_ma87_complex.mod hsl_ma87_double_complex.mod)
SET(${LIB_PREFIX}_INC_SUFFIXES "include" "hsl_ma87" "include/hsl_ma87")
SET(${LIB_PREFIX}_LIB_SUFFIXES "lib" "hsl_ma87" "lib/hsl_ma87" "lib/${CMAKE_LIBRARY_ARCHITECTURE}" "${CMAKE_LIBRARY_ARCHITECTURE}")

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
