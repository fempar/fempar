#################################################################
# Serial tests
#################################################################

IF(MPI_FOUND)

    FILE(GLOB DIRECTORIES *)

    FOREACH( DIR ${DIRECTORIES})
        FILE(GLOB TESTS_SERIAL_SRC ${DIR}/*.f90 ${DIR}/*.F90)

        FOREACH(TEST_SRC ${TESTS_SERIAL_SRC})
            GET_FILENAME_COMPONENT(EXE_NAME ${TEST_SRC} NAME_WE)
            SET(${EXE_NAME}_MODULES_DIR "${DIR}/modules")
            UNSET(${EXE_NAME}_MODULES_SRC)
            IF(IS_DIRECTORY ${${EXE_NAME}_MODULES_DIR})
                FILE(GLOB ${EXE_NAME}_MODULES_SRC ${${EXE_NAME}_MODULES_DIR}/*.f90 ${${EXE_NAME}_MODULES_DIR}/*.F90)
            ENDIF()
            ADD_EXECUTABLE(${EXE_NAME} ${TEST_SRC} ${${EXE_NAME}_MODULES_SRC})
            TARGET_LINK_LIBRARIES(${EXE_NAME} ${PROJECT_NAME})
            TARGET_LINK_LIBRARIES(${EXE_NAME} ${MPI_Fortran_LIBRARIES})
            TARGET_LINK_LIBRARIES(${EXE_NAME} ${${PROJECT_NAME}_EXTERNAL_LIBRARIES} ${${PROJECT_NAME}_EXTERNAL_PROJECTS})
            IF(EXISTS ${SCRIPTS_TESTS_PATH}/${EXE_NAME})
                ADD_TEST(${EXE_NAME}_TEST ${SCRIPTS_TESTS_PATH}/${EXE_NAME} ${EXECUTABLE_OUTPUT_PATH}/${EXE_NAME} ${DATA_TESTS_PATH}/${EXE_NAME})
            ENDIF()
        ENDFOREACH()
    ENDFOREACH()

ENDIF()



