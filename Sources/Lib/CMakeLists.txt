#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE LIBRARY_SRC *.f90 *.F90 *.c)
SET(LIBRARY_SRC ${LIBRARY_SRC} PARENT_SCOPE)

#################################################################
# ADD FEMPAR LIBRARY TARGET
#################################################################

ADD_LIBRARY(${PROJECT_NAME} ${LIBRARY_SRC})

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${${PROJECT_NAME}_EXTERNAL_LIBRARIES} ${${PROJECT_NAME}_EXTERNAL_PROJECTS})

SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})
EXPORT(TARGETS ${PROJECT_NAME}  FILE ${BINARY_OUTPUT_PATH}/${PROJECT_NAME}Targets.cmake)
