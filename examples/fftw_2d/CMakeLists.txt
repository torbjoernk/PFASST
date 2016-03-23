set(examples_to_install ${examples_to_install})

message(STATUS "  FFT Test")
include_directories(
    ${3rdparty_INCLUDES}
    ${pfasst_INCLUDES}
    ${FFTW3_INCLUDE_DIR}
)

set(fftw_examples
    fftw_test
)

set(all_fftw_examples ${fftw_examples})


foreach(example ${fftw_examples})
    add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
    add_dependencies(${example} ${pfasst_DEPENDEND_TARGETS})
    target_link_libraries(${example}
        fftw3_lib
        ${3rdparty_DEPENDEND_LIBS}
        ${pfasst_DEPENDEND_LIBS}
    )
    if(pfasst_INSTALL_EXAMPLES)
        install(TARGETS ${example} RUNTIME DESTINATION bin)
    endif()
endforeach(example)

set(examples_to_install ${examples_to_install} ${all_fftw_examples} PARENT_SCOPE)
