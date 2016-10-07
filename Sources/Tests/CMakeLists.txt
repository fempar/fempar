#################################################################
# Serial and Parallel tests
#################################################################

ADD_SUBDIRECTORY(${SERIAL_TESTS_PATH})
ADD_SUBDIRECTORY(${PAR_TESTS_PATH})

#################################################################
# Custom clean-Tests target
#################################################################

ADD_CUSTOM_TARGET(clean-Tests
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR} --target clean 
)
