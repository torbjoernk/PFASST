message(STATUS "trying to use system's Dune")

find_package(dune-common QUIET CONFIG)
find_package(dune-istl QUIET CONFIG)
find_package(dune-localfunctions QUIET CONFIG)
find_package(dune-grid QUIET CONFIG)

if(dune-common_DIR AND dune-istl_DIR AND dune-localfunctions_DIR AND dune-grid_DIR)
    message(STATUS "  using components:")
    message(STATUS "    dune-common")
    message(STATUS "    dune-istl")
    message(STATUS "    dune-local")
    message(STATUS "    dune-grid")

    set(dune_INCLUDE_PATHS
        ${dune-common_INCLUDE_DIRS}
        ${dune-istl_INCLUDE_DIRS}
        ${dune-localfunctions_INCLUDE_DIRS}
        ${dune-grid_INCLUDE_DIRS}
    )
    list(REMOVE_DUPLICATES dune_INCLUDE_PATHS)
    set(dune_LIBRARIES
        ${dune-common_LIBRARIES}
        ${dune-istl_LIBRARIES}
        ${dune-localfunctions_LIBRARIES}
        ${dune-grid_LIBRARIES}
    )
    set(dune_FOUND TRUE)

    message(STATUS "  include path: ${dune_INCLUDE_PATHS}")
    message(STATUS "  libraries: ${dune_LIBRARIES}")
else()
    message(STATUS "!! not found")
    set(dune_FOUND FALSE)
endif()

set(dune_FOUND ${dune_FOUND} PARENT_SCOPE)
set(dune_INCLUDE_PATHS ${dune_INCLUDE_PATHS} PARENT_SCOPE)
set(dune_LIBRARIES ${dune_LIBRARIES} PARENT_SCOPE)
