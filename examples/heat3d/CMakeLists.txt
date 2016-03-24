set(examples_to_install ${examples_to_install})

message(STATUS "  Heat 3D")
include_directories(
    ${3rdparty_INCLUDES}
    ${pfasst_INCLUDES}
    ${FFTW3_INCLUDE_PATH}
)

set(heat3d_examples
    heat3d_sdc
    heat3d_mlsdc
)

set(heat3d_mpi_examples
#     heat3d_pfasst
)

if(${pfasst_WITH_MPI})
    set(all_heat3d_examples ${heat3d_mpi_examples})
else()
    set(all_heat3d_examples ${heat3d_examples})
endif()


if(NOT ${pfasst_WITH_MPI})
    foreach(example ${heat3d_examples})
        add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
        add_dependencies(${example} ${pfasst_DEPENDEND_TARGETS})
        target_link_libraries(${example}
            ${3rdparty_DEPENDEND_LIBS}
            ${pfasst_DEPENDEND_LIBS}
            ${FFTW3_LIBRARIES}
        )
        if(pfasst_INSTALL_EXAMPLES)
            install(TARGETS ${example} RUNTIME DESTINATION bin)
        endif()
    endforeach(example)
else()
    include_directories(${MPI_CXX_INCLUDE_PATH})
    foreach(example ${heat3d_mpi_examples})
        add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
        add_dependencies(${example} ${pfasst_DEPENDEND_TARGETS})
        target_link_libraries(${example}
            ${3rdparty_DEPENDEND_LIBS}
            ${pfasst_DEPENDEND_LIBS}
            ${FFTW3_LIBRARIES}
        )
        if(pfasst_INSTALL_EXAMPLES)
            install(TARGETS ${example} RUNTIME DESTINATION bin)
        endif()
    endforeach(example)
endif()

set(examples_to_install ${examples_to_install} ${all_heat3d_examples} PARENT_SCOPE)
