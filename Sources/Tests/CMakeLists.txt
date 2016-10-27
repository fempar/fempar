#################################################################
# Serial and Parallel tests
#################################################################

SET(TEMP_Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})
ADD_SUBDIRECTORY(${SERIAL_TESTS_PATH})
ADD_SUBDIRECTORY(${PAR_TESTS_PATH})
SET(CMAKE_Fortran_MODULE_DIRECTORY ${TEMP_Fortran_MODULE_DIRECTORY})

#################################################################
# Custom clean-Tests target
#################################################################

ADD_CUSTOM_TARGET(clean-tests
    COMMAND ${CMAKE_MAKE_PROGRAM} clean 
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
