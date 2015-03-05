# CMake script to detect WSMP library: 
#
# Usage example:
#   find_package(WSMP)
#   if (WSMP_FOUND)
#      include_directories(${WSMP_INCLUDE_DIRS})
#      link_directories(${WSMP_LIBRARY_DIRS})
#      add_executable(foo foo.cc)
#      target_link_libraries(foo ${WSMP_LIBRARIES})
#   endif()
#
# On return this will define:
#   WSMP_FOUND                   Indicates whether WSMP was found (True/False)
#   WSMP_INCLUDE_DIRS            WSMP include folder
#   WSMP_LIBRARIES               WSMP libraries names
#
# WSMP package can contains the following libraries:
#   libwsmp.a:      32-bit serial/multithreaded library
#   libwsmp64.a:    64-bit serial/multithreaded library
#   libwsmp8_8.a:   64-bit serial/multithreaded library w/ 8-byte integers
#   libwsmp64gpu.a:   GPU-accelerated multithreaded library; needs Qccel_WSMP library from Nvidia
#   libpwsmp64.a:   64-bit MPI-based library
#
####################################################################################

SET(DEFAULT_INC_PATHS /include /usr/include /usr/local/include)
SET(DEFAULT_LIB_PATHS /lib /lib32 /lib64 /usr/lib /usr/lib32 /usr/lib64 /usr/local/lib)

SET(LIB_PREFIX WSMP)
SET(${LIB_PREFIX}_LIB_NAME wsmp wsmp64 wsmp8_8 wsmp64gpu pwsmp64)
SET(${LIB_PREFIX}_INC_FILENAME wsmp.h)
SET(${LIB_PREFIX}_INC_SUFFIXES "include" "wsmp" "include/wsmp")
SET(${LIB_PREFIX}_LIB_SUFFIXES "lib" "wsmp" "lib/wsmp" "lib/${CMAKE_LIBRARY_ARCHITECTURE}" "${CMAKE_LIBRARY_ARCHITECTURE}")

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
