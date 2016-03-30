set(examples_to_install ${examples_to_install})

message(STATUS "  Heat 1D")
include_directories(
    ${3rdparty_INCLUDES}
    ${pfasst_INCLUDES}
    ${FFTW3_INCLUDE_PATH}
)

set(heat_examples
    heat1d_sdc
    heat1d_mlsdc
)

set(heat_mpi_examples
    heat1d_pfasst
)

if(${pfasst_WITH_MPI})
    set(all_heat_examples ${heat_mpi_examples})
else()
    set(all_heat_examples ${heat_examples})
endif()


if(NOT ${pfasst_WITH_MPI})
    foreach(example ${heat_examples})
        add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
        add_dependencies(${example} ${pfasst_DEPENDEND_TARGETS})
        target_link_libraries(${example}
            ${3rdparty_DEPENDEND_LIBS}
            ${pfasst_DEPENDEND_LIBS}
            ${FFTW3_LIBRARIES}
        )
        target_compile_features(${example} PRIVATE ${pfasst_CXX_FEATURES})
        if(pfasst_INSTALL_EXAMPLES)
            install(TARGETS ${example} RUNTIME DESTINATION bin)
        endif()
    endforeach(example)
else()
    include_directories(${MPI_CXX_INCLUDE_PATH})
    foreach(example ${heat_mpi_examples})
        add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
        add_dependencies(${example} ${pfasst_DEPENDEND_TARGETS})
        target_link_libraries(${example}
            ${3rdparty_DEPENDEND_LIBS}
            ${pfasst_DEPENDEND_LIBS}
            ${FFTW3_LIBRARIES}
        )
        target_compile_features(${example} PRIVATE ${pfasst_CXX_FEATURES})
        if(pfasst_INSTALL_EXAMPLES)
            install(TARGETS ${example} RUNTIME DESTINATION bin)
        endif()
    endforeach(example)
endif()

set(examples_to_install ${examples_to_install} ${all_heat_examples} PARENT_SCOPE)
