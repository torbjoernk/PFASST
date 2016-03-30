set(examples_to_install ${examples_to_install})

message(STATUS "  Heat 2D")
include_directories(
    ${3rdparty_INCLUDES}
    ${pfasst_INCLUDES}
    ${FFTW3_INCLUDE_PATH}
)

set(heat2d_examples
    heat2d_sdc
    heat2d_mlsdc
)

set(heat2d_mpi_examples
    heat2d_pfasst
)

if(${pfasst_WITH_MPI})
    set(all_heat2d_examples ${heat2d_mpi_examples})
else()
    set(all_heat2d_examples ${heat2d_examples})
endif()


if(NOT ${pfasst_WITH_MPI})
    foreach(example ${heat2d_examples})
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
    foreach(example ${heat2d_mpi_examples})
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

set(examples_to_install ${examples_to_install} ${all_heat2d_examples} PARENT_SCOPE)
