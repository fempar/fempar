# CMake script to detect Metis library: 
#
# Usage example:
#   find_package(METIS)
#   if (METIS_FOUND)
#      include_directories(${METIS_INCLUDE_DIRS})
#      link_directories(${METIS_LIBRARY_DIRS})
#      add_executable(foo foo.cc)
#      target_link_libraries(foo ${METIS_LIBRARIES})
#   endif()
#
# On return this will define:
#   METIS_FOUND                   Indicates whether METIS was found (True/False)
#   METIS_INCLUDE_DIRS            METIS include folder
#   METIS_LIBRARIES               METIS libraries names
#
####################################################################################

SET(DEFAULT_INC_PATHS /include /usr/include /usr/local/include)
SET(DEFAULT_LIB_PATHS /lib /lib32 /lib64 /usr/lib /usr/lib32 /usr/lib64 /usr/local/lib)

SET(LIB_PREFIX UMFPACK)
SET(${LIB_PREFIX}_LIB_NAME umfpack)
SET(${LIB_PREFIX}_INC_FILENAME umfpack.h)
SET(${LIB_PREFIX}_INC_SUFFIXES "include" "umfpack" "include/umfpack")
SET(${LIB_PREFIX}_LIB_SUFFIXES "lib" "umfpack" "lib/umfpack" "lib/${CMAKE_LIBRARY_ARCHITECTURE}" "${CMAKE_LIBRARY_ARCHITECTURE}")

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

IF(${LIB_PREFIX}_LIBS)
    SET(${LIB_PREFIX}_INCLUDE_DIR ${${LIB_PREFIX}_INCLUDES})
    SET(${LIB_PREFIX}_LIBRARIES ${${LIB_PREFIX}_LIBS})
    SET(${LIB_PREFIX}_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(${LIB_PREFIX}_INCLUDES)
MARK_AS_ADVANCED(${LIB_PREFIX}_LIBS)
