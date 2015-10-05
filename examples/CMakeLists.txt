set(examples_to_install ${examples_to_install})

add_subdirectory(advec_diff)
add_subdirectory(heat1d)
add_subdirectory(heat2d)
add_subdirectory(heat3d)
add_subdirectory(fftw_2d)

set(examples_to_install ${examples_to_install} PARENT_SCOPE)
