#################################################################
# Executable targets
#################################################################
#SET(EXAMPLE_TESTS_PATH ${SRC_PATH}/Tests/examples)
#ADD_EXECUTABLE(part_new ${EXAMPLE_TESTS_PATH}/part_new.f90)
#TARGET_LINK_LIBRARIES(part_new ${FEM_LIB})



#################################################################
# Parallel tests
#################################################################

IF(MPI_Fortran_COMPILER)

    INCLUDE(CMakeForceCompiler)
    CMAKE_FORCE_Fortran_COMPILER(${MPI_Fortran_COMPILER} ${CMAKE_Fortran_COMPILER_ID})

    FILE(GLOB_RECURSE EXAMPLE_TESTS_PAR_SRC *.f90)

    FOREACH(TEST_SRC ${EXAMPLE_TESTS_PAR_SRC})
        GET_FILENAME_COMPONENT(EXE_NAME ${TEST_SRC} NAME_WE)
        ADD_EXECUTABLE(${EXE_NAME} ${TEST_SRC})
        TARGET_LINK_LIBRARIES(${EXE_NAME} ${FEM_LIB})
        TARGET_LINK_LIBRARIES(${EXE_NAME} ${PAR_LIB})
#        TARGET_LINK_LIBRARIES(${EXE_NAME} ${MPI_Fortran_LIBRARIES})
        FOREACH (EXT_LIB ${EXT_LIBS})
            IF(DEFINED ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${PROJECT_NAME}_ENABLE_${EXT_LIB} AND ${EXT_LIB}_FOUND)
                TARGET_LINK_LIBRARIES(${EXE_NAME} ${${EXT_LIB}_LIBRARIES})
    #            GET_FILENAME_COMPONENT(LIB_PATH ${${EXT_LIB}_LIBRARIES} DIRECTORY)
    #            MESSAGE("LIB_PATH: ${LIB_PATH}")
    #            LINK_DIRECTORIES(${LIB_PATH})
            ENDIF()
        ENDFOREACH()
        ADD_TEST(${EXE_NAME}_TEST 
             ${EXECUTABLE_OUTPUT_PATH}/${EXE_NAME} ${CMAKE_SOURCE_DIR}/Sources/Tests/Data/${EXE_NAME}.txt)
    ENDFOREACH()
ENDIF()