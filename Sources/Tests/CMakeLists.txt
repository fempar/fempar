#################################################################
# Serial and Parallel tests
#################################################################

# Change the module output directory for tests
SET(TEMP_Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Copy test scripts to build dir
FILE(COPY ${SCRIPTS_TESTS_PATH} DESTINATION ${TESTS_OUTPUT_PATH} USE_SOURCE_PERMISSIONS)
# Copy test data to build dir
FILE(COPY ${DATA_TESTS_PATH} DESTINATION ${TESTS_OUTPUT_PATH} USE_SOURCE_PERMISSIONS)

# Add serial tests subdirectory
ADD_SUBDIRECTORY(${SERIAL_TESTS_PATH})
# Add parallel tests subdirectory
ADD_SUBDIRECTORY(${PAR_TESTS_PATH})
# Restore the module output directory
SET(CMAKE_Fortran_MODULE_DIRECTORY ${TEMP_Fortran_MODULE_DIRECTORY})

#################################################################
# Custom clean-Tests target
#################################################################

# Create clean-tests custom target
ADD_CUSTOM_TARGET(clean-tests
    COMMAND ${CMAKE_MAKE_PROGRAM} clean 
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
