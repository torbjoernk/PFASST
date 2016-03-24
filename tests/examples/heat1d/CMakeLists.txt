# Building and Running Tests
include_directories(
    ${3rdparty_INCLUDES}
    ${TESTS_3rdparty_INCLUDES}
    ${pfasst_INCLUDES}
    ${PROJECT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FFTW3_INCLUDE_PATH}
)

set(TESTS
    test_heat1d_sweeper
    test_heat1d_convergence
)

set(TESTS_MPI
)

if(NOT ${pfasst_WITH_MPI})
    foreach(test ${TESTS})
        message(STATUS "  ${test}")
        add_executable(${test} ${test}.cpp)
        add_dependencies(${test} ${pfasst_DEPENDEND_TARGETS})
        target_link_libraries(${test}
            ${3rdparty_DEPENDEND_LIBS}
            ${TESTS_3rdparty_DEPENDEND_LIBS}
            ${pfasst_DEPENDEND_LIBS}
            ${FFTW3_LIBRARIES}
        )
        if(pfasst_WITH_GCC_PROF AND ${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
            set_target_properties(${test}
                PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -ftest-coverage -fprofile-arcs"
                           LINK_FLAGS "-fprofile-arcs"
            )
        endif()
        add_test(NAME ${test}
            COMMAND ${CMAKE_BINARY_DIR}/tests/examples/heat1d/${test} --gtest_output=xml:${test}_out.xml
        )
    endforeach(test)
else()
    include_directories(${MPI_CXX_INCLUDE_PATH})
    foreach(test ${TESTS_MPI})
        message(STATUS "  ${test}")
        add_executable(${test} ${test}.cpp)
        add_dependencies(${test} ${pfasst_DEPENDEND_TARGETS})
        target_link_libraries(${test}
            ${3rdparty_DEPENDEND_LIBS}
            ${TESTS_3rdparty_DEPENDEND_LIBS}
            ${pfasst_DEPENDEND_LIBS}
            ${FFTW3_LIBRARIES}
        )
        if(pfasst_WITH_GCC_PROF AND ${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
            set_target_properties(${test}
                PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -ftest-coverage -fprofile-arcs"
                           LINK_FLAGS "-fprofile-arcs"
            )
        endif()
        add_test(NAME ${test}
            COMMAND ${CMAKE_BINARY_DIR}/tests/examples/heat1d/${test} --gtest_output=xml:${test}_out.xml
        )
    endforeach(test)
endif()
