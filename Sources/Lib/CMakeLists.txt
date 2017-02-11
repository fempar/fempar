#################################################################
# Search F90 files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE LIBRARY_SRC *.f90 *.F90 *.c) # Recursively search files matching the Regex and store them in LIBRARY_SRC
SET(LIBRARY_SRC ${LIBRARY_SRC} PARENT_SCOPE)   # Publish the list of files in the parent scope (main CMakeLists.txt)

#################################################################
# ADD FEMPAR LIBRARY TARGET
#################################################################

# Create a library target (FEMPAR) with the sources in LIBRARY_SRC
ADD_LIBRARY(${PROJECT_NAME} ${LIBRARY_SRC})

# Link FEMPAR target with extenal libraries and projects
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${${PROJECT_NAME}_EXTERNAL_LIBRARIES} ${${PROJECT_NAME}_EXTERNAL_PROJECTS})

# Set version properties to library
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} SOVERSION ${${PROJECT_NAME}_SOVERSION})

# Export FEMPARTargets.cmake file
EXPORT(TARGETS ${PROJECT_NAME}  FILE ${BINARY_OUTPUT_PATH}/${PROJECT_NAME}Targets.cmake)
