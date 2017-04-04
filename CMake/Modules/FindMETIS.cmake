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

INCLUDE(CheckLibraryExists)

SET(DEFAULT_INC_PATHS /include /usr/include /usr/local/include)
SET(DEFAULT_LIB_PATHS /lib /lib32 /lib64 /usr/lib /usr/lib32 /usr/lib64 /usr/local/lib)

SET(LIB_PREFIX METIS)
SET(${LIB_PREFIX}_LIB_NAME metis)
SET(${LIB_PREFIX}_INC_FILENAME metis.h)
SET(${LIB_PREFIX}_INC_SUFFIXES "include" "libmetis" "include/metis")
SET(${LIB_PREFIX}_LIB_SUFFIXES "lib" "libmetis" "lib/metis" "lib/${CMAKE_LIBRARY_ARCHITECTURE}" "${CMAKE_LIBRARY_ARCHITECTURE}")
SET(${LIB_PREFIX}_EXT_DEPS M )

SET(${LIB_PREFIX}_FOUND FALSE)

FIND_PATH(${LIB_PREFIX}_INCLUDES 
    NAMES ${${LIB_PREFIX}_INC_FILENAME}
    PATHS ${DEFAULT_INC_PATHS} ENV PATH
    PATH_SUFFIXES ${${LIB_PREFIX}_INC_SUFFIXES}
    DOC "The path to METIS header files"
)

FIND_LIBRARY(${LIB_PREFIX}_LIBS
    NAMES ${${LIB_PREFIX}_LIB_NAME}
    PATHS ${DEFAULT_LIB_PATHS} ENV LD_LIBRARY_PATH
    PATH_SUFFIXES ${${LIB_PREFIX}_LIB_SUFFIXES}
    DOC "The path to ${LIB_PREFIX} library"
)


IF(${LIB_PREFIX}_LIBS)
    UNSET(HAVE_${LIB_PREFIX}_SETDEFAULTOPTIONS CACHE)
    GET_FILENAME_COMPONENT(${LIB_PREFIX}_PATH ${${LIB_PREFIX}_LIBS} PATH) #PATH legacy alias for DIRECTORY
    GET_FILENAME_COMPONENT(${LIB_PREFIX}_EXT ${${LIB_PREFIX}_LIBS} EXT)

    IF (${${LIB_PREFIX}_EXT} MATCHES ${CMAKE_STATIC_LIBRARY_SUFFIX})
        FOREACH (EXT_DEP ${${LIB_PREFIX}_EXT_DEPS})
            FIND_PACKAGE( ${EXT_DEP} )
            IF(${EXT_DEP}_FOUND)
                SET(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${${EXT_DEP}_LIBS} )
                SET(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${${EXT_DEP}_INCLUDES})
#               IF(IS_DIRECTORY ${${EXT_DEP}_INCLUDE_DIR})
#                   SET(${LIB_PREFIX}_INCLUDES ${${LIB_PREFIX}_INCLUDES} ${${EXT_DEP}_INCLUDES})
#               ENDIF()
#               SET(${LIB_PREFIX}_LIBRARIES ${${LIB_PREFIX}_LIBRARIES} ${${EXT_DEP}_LIBRARIES})
            ENDIF()
            UNSET(${EXT_DEP}_INCLUDES CACHE)
            UNSET(${EXT_DEP}_LIBS CACHE)
            UNSET(${EXT_DEP}_INCLUDE_DIR CACHE)
            UNSET(${EXT_DEP}_LIBRARIES CACHE)
        ENDFOREACH()
    ELSE()

    ENDIF()

    CHECK_LIBRARY_EXISTS(${${LIB_PREFIX}_LIB_NAME} METIS_SetDefaultOptions ${${LIB_PREFIX}_PATH} HAVE_${LIB_PREFIX}_SETDEFAULTOPTIONS)
    IF(HAVE_${LIB_PREFIX}_SETDEFAULTOPTIONS)
        SET(${LIB_PREFIX}_MAJOR_VERSION 5)
        SET(${LIB_PREFIX}_INCLUDE_DIR ${${LIB_PREFIX}_INCLUDES})
        SET(${LIB_PREFIX}_LIBRARIES ${${LIB_PREFIX}_LIBS})
        SET(${LIB_PREFIX}_FOUND TRUE)
    ELSE()
#        SET(${LIB_PREFIX}_MAJOR_VERSION 4)
    ENDIF()


ENDIF()

MARK_AS_ADVANCED(${LIB_PREFIX}_INCLUDES)
MARK_AS_ADVANCED(${LIB_PREFIX}_LIBS)

