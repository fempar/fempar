#################################################################
# Serial and Parallel tests
#################################################################

ADD_SUBDIRECTORY(${SERIAL_TESTS_PATH})
ADD_SUBDIRECTORY(${PAR_TESTS_PATH})

#################################################################
# Custom clean-Tests target
#################################################################

ADD_CUSTOM_TARGET(clean-tests
    COMMAND ${CMAKE_MAKE_PROGRAM} clean 
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
