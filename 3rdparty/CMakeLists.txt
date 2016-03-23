# Managing 3rd-party libraries and tools
set(3rdparty_INCLUDES ${3rdparty_INCLUDES})
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS})
set(pfasst_DEPENDEND_TARGETS ${pfasst_DEPENDEND_TARGETS})
set(pfasst_TESTS_DEPENDEND_TARGETS ${pfasst_TESTS_DEPENDEND_TARGETS})

if(3rdparty_SOURCES)
    # TODO: only define if dirs present
    if(NOT Boost_SOURCES)
        set(Boost_SOURCES ${3rdparty_SOURCES}/boost)
    endif()
    if(NOT Eigen3_SOURCES)
        set(Eigen3_SOURCES ${3rdparty_SOURCES}/eigen3)
    endif()
    if(NOT fftw3_SOURCES)
        set(fftw3_SOURCES ${3rdparty_SOURCES}/fftw3)
    endif()
    if(NOT gtest_SOURCES)
        set(gtest_SOURCES ${3rdparty_SOURCES}/gtest)
    endif()
    if(NOT gmock_SOURCES)
        set(gmock_SOURCES ${3rdparty_SOURCES}/gmock)
    endif()
endif()
set(3rdparty_DOWNLOAD_TIMEOUT "60")


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Leathers")
list(APPEND 3rdparty_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/leathers)


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Better Enums")
list(APPEND 3rdparty_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/better-enums)


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Boost")

set(Boost_SOURCE_NAME "boost_1_60_0.tar.bz2")
set(Boost_SOURCE_MD5 "65a840e1a0b13a558ff19eeb2c4f0cbe")
message(STATUS "  using version: 1.60.0  (file: ${Boost_SOURCE_NAME})")

if(NOT Boost_SOURCES)
    set(Boost_SOURCE_LOCATION "https://downloads.sourceforge.net/project/boost/boost/1.60.0")
    message(STATUS "  going to download from: ${Boost_SOURCE_LOCATION}")
    message(STATUS "    HINT: specify -DBoost_SOURCES=/local/directory/with/${Boost_SOURCE_NAME} to avoid redownloading")
else()
    set(Boost_SOURCE_LOCATION ${Boost_SOURCES})
    message(STATUS "  using pre-downloaded sources at: ${Boost_SOURCE_LOCATION}")
endif()

set(Boost_COMPONENTS "program_options")

set(Boost_PREFIX "${3rdparty_PREFIX}/Boost")
set(Boost_TMP_PREFIX "${3rdparty_TMP_PREFIX}/Boost")

set(Boost_TMP_DIR "${Boost_TMP_PREFIX}/tmp")
set(Boost_STAMP_DIR "${Boost_TMP_PREFIX}/stamp")
set(Boost_DOWNLOAD_DIR "${Boost_PREFIX}/download")
set(Boost_SOURCE_DIR "${Boost_PREFIX}/src")
set(Boost_BINARY_DIR "${Boost_TMP_PREFIX}/build")
set(Boost_INSTALL_DIR "${Boost_PREFIX}/install")

set(Boost_build_CMD "${Boost_SOURCE_DIR}/bootstrap.sh")
list(APPEND Boost_build_CMD
    --prefix=${Boost_INSTALL_DIR}
    --with-libraries=${Boost_COMPONENTS}
)
if(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
    list(APPEND Boost_build_CMD toolset=clang)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    list(APPEND Boost_build_CMD toolset=gcc)
endif()

ExternalProject_Add(boost_build
    LIST_SEPARATOR " "
    TMP_DIR "${Boost_TMP_DIR}"
    STAMP_DIR "${Boost_STAMP_DIR}"
    DOWNLOAD_DIR "${Boost_DOWNLOAD_DIR}"
    SOURCE_DIR "${Boost_SOURCE_DIR}"
    # BINARY_DIR "${Boost_BINARY_DIR}"  # not used as in-source build is specified
    INSTALL_DIR "${Boost_INSTALL_DIR}"

    EXCLUDE_FROM_ALL ON

    URL ${Boost_SOURCE_LOCATION}/${Boost_SOURCE_NAME}
    URL_MD5 ${Boost_SOURCE_MD5}
    TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}

    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${Boost_build_CMD}
    TEST_COMMAND ""
    INSTALL_COMMAND ""

    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

set(BOOST_LIBRARIES)
set(BOOST_LIBRARY_DIR ${Boost_SOURCE_DIR}/stage/lib)
set(BOOST_INCLUDE_DIR ${Boost_SOURCE_DIR})

message(STATUS "  using components:")
foreach(component ${Boost_COMPONENTS})
    message(STATUS "    ${component}")

    set(Boost_${component}_BUILD_CMD ${Boost_SOURCE_DIR}/b2)
    list(APPEND Boost_${component}_BUILD_CMD
        --build-dir=${Boost_BINARY_DIR}
        --with-${component}
        address-model=64
        cxxflags=-std=c++11
        cxxflags=-fPIC
        link=static
        threading=multi
        runtime-link=shared
        variant=release
    )
    if(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
        list(APPEND Boost_${component}_BUILD_CMD toolset=clang)
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        list(APPEND Boost_${component}_BUILD_CMD toolset=gcc)
    endif()

    ExternalProject_Add(boost_${component}
        LIST_SEPARATOR " "
        TMP_DIR "${Boost_TMP_DIR}"
        STAMP_DIR "${Boost_STAMP_DIR}"
        DOWNLOAD_DIR "${Boost_DOWNLOAD_DIR}"
        SOURCE_DIR "${Boost_SOURCE_DIR}"
        # BINARY_DIR "${Boost_BINARY_DIR}"  # not used as in-source build is specified
        INSTALL_DIR "${Boost_INSTALL_DIR}"

        EXCLUDE_FROM_ALL ON
        DEPENDS boost_build

        DOWNLOAD_COMMAND ""
        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        BUILD_IN_SOURCE ON
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${Boost_${component}_BUILD_CMD} -j4 stage
        TEST_COMMAND ""
        INSTALL_COMMAND ""

        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON
        LOG_INSTALL ON
    )

    list(APPEND pfasst_DEPENDEND_TARGETS boost_${component})
    list(APPEND BOOST_LIBRARIES ${BOOST_LIBRARY_DIR}/libboost_${component}.a)
endforeach()

message(STATUS "  Boost include path: ${BOOST_INCLUDE_DIR}")
message(STATUS "  Boost library path: ${BOOST_LIBRARY_DIR}")
list(APPEND 3rdparty_INCLUDES ${BOOST_INCLUDE_DIR})
list(APPEND 3rdparty_DEPENDEND_LIBS ${BOOST_LIBRARIES})


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Eigen3")

set(Eigen3_SOURCE_NAME "3.2.5.tar.bz2")
set(Eigen3_SOURCE_MD5 "21a928f6e0f1c7f24b6f63ff823593f5")
message(STATUS "  using version: 3.2.5  (file: ${Eigen3_SOURCE_NAME})")

if(NOT Eigen3_SOURCES)
    set(Eigen3_SOURCE_LOCATION "https://bitbucket.org/eigen/eigen/get/")
    message(STATUS "  going to download from: ${Eigen3_SOURCE_LOCATION}")
    message(STATUS "    HINT: specify -DEigen3_SOURCES=/local/directory/with/${Eigen3_SOURCE_NAME} to avoid redownloading")
else()
    set(Eigen3_SOURCE_LOCATION ${Eigen3_SOURCES})
    message(STATUS "  using pre-downloaded sources at: ${Eigen3_SOURCE_LOCATION}")
endif()

set(Eigen3_PREFIX "${3rdparty_PREFIX}/Eigen3")
set(Eigen3_TMP_PREFIX "${3rdparty_TMP_PREFIX}/Eigen3")

set(Eigen3_TMP_DIR "${Eigen3_TMP_PREFIX}/tmp")
set(Eigen3_STAMP_DIR "${Eigen3_TMP_PREFIX}/stamp")
set(Eigen3_DOWNLOAD_DIR "${Eigen3_PREFIX}/download")
set(Eigen3_SOURCE_DIR "${Eigen3_PREFIX}/src")
set(Eigen3_BINARY_DIR "${Eigen3_TMP_PREFIX}/build")
set(Eigen3_INSTALL_DIR "${Eigen3_PREFIX}/install")

ExternalProject_Add(Eigen3
    LIST_SEPARATOR " "
    TMP_DIR "${Eigen3_TMP_DIR}"
    STAMP_DIR "${Eigen3_STAMP_DIR}"
    DOWNLOAD_DIR "${Eigen3_DOWNLOAD_DIR}"
    SOURCE_DIR "${Eigen3_SOURCE_DIR}"
    # BINARY_DIR "${Eigen3_BINARY_DIR}"  # not used as in-source build is specified
    INSTALL_DIR "${Eigen3_INSTALL_DIR}"

    EXCLUDE_FROM_ALL ON

    URL ${Eigen3_SOURCE_LOCATION}/${Eigen3_SOURCE_NAME}
    URL_MD5 ${Eigen3_SOURCE_MD5}
    TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}

    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    TEST_COMMAND ""
    INSTALL_COMMAND ""

    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
    LOG_INSTALL ON
)

# Specify include dir
set(Eigen3_INCLUDE_DIR "${Eigen3_SOURCE_DIR}")
message(STATUS "  include path: ${Eigen3_INCLUDE_DIR}")

list(APPEND pfasst_DEPENDEND_TARGETS Eigen3)
list(APPEND 3rdparty_INCLUDES "${Eigen3_INCLUDE_DIR}")


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "FFTW3")

set(fftw3_SOURCE_NAME "fftw-3.3.4.tar.gz")
set(fftw3_SOURCE_MD5 "2edab8c06b24feeb3b82bbb3ebf3e7b3")
message(STATUS "  using version: 3.3.4  (file: ${fftw3_SOURCE_NAME})")

if(NOT fftw3_SOURCES)
    set(fftw3_SOURCE_LOCATION "http://www.fftw.org/")
    message(STATUS "  going to download from: ${fftw3_SOURCE_LOCATION}")
    message(STATUS "    HINT: specify -Dfftw3_SOURCES=/local/directory/with/${fftw3_SOURCE_NAME} to avoid redownloading")
else()
    set(fftw3_SOURCE_LOCATION ${fftw3_SOURCES})
    message(STATUS "  using pre-downloaded sources at: ${fftw3_SOURCE_LOCATION}")
endif()

set(fftw3_PREFIX "${3rdparty_PREFIX}/fftw3")
set(fftw3_TMP_PREFIX "${3rdparty_TMP_PREFIX}/fftw3")

set(fftw3_TMP_DIR "${fftw3_TMP_PREFIX}/tmp")
set(fftw3_STAMP_DIR "${fftw3_TMP_PREFIX}/stamp")
set(fftw3_DOWNLOAD_DIR "${fftw3_PREFIX}/download")
set(fftw3_SOURCE_DIR "${fftw3_PREFIX}/src")
set(fftw3_BINARY_DIR "${fftw3_TMP_PREFIX}/build")
set(fftw3_INSTALL_DIR "${fftw3_PREFIX}/install")

ExternalProject_Add(fftw3
    LIST_SEPARATOR " "
    TMP_DIR "${fftw3_TMP_DIR}"
    STAMP_DIR "${fftw3_STAMP_DIR}"
    DOWNLOAD_DIR "${fftw3_DOWNLOAD_DIR}"
    SOURCE_DIR "${fftw3_SOURCE_DIR}"
    # BINARY_DIR "${fftw3_BINARY_DIR}"  # not used as in-source build is specified
    INSTALL_DIR "${fftw3_INSTALL_DIR}"

    EXCLUDE_FROM_ALL ON

    URL ${fftw3_SOURCE_LOCATION}/${fftw3_SOURCE_NAME}
    URL_MD5 ${fftw3_SOURCE_MD5}
    TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}

    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ${fftw3_SOURCE_DIR}/configure
        CC=${CMAKE_C_COMPILER}
        CXX=${CMAKE_CXX_COMPILER}
        CXXFLAGS=${CMAKE_CXX_FLAGS}
        --prefix=${fftw3_INSTALL_DIR}
        --libdir=${fftw3_INSTALL_DIR}/lib
    BUILD_COMMAND make -j4
    TEST_COMMAND ""
    INSTALL_COMMAND make install

    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
    LOG_INSTALL ON
)

set(FFTW3_INCLUDE_DIR ${fftw3_INSTALL_DIR}/include)
set(FFTW3_LIBRARY_DIR ${fftw3_INSTALL_DIR}/lib)
add_library(fftw3_lib STATIC IMPORTED GLOBAL)
set_property(TARGET fftw3_lib
    PROPERTY IMPORTED_LOCATION ${FFTW3_LIBRARY_DIR}/libfftw3.a
)
add_dependencies(fftw3_lib fftw3)
message(STATUS "  include path: ${FFTW3_INCLUDE_DIR}")
message(STATUS "  library path: ${FFTW3_LIBRARY_DIR}")

set(FFTW3_INCLUDE_DIR ${FFTW3_INCLUDE_DIR} PARENT_SCOPE)


if(pfasst_BUILD_TESTS)
    set(TESTS_3rdparty_INCLUDES ${TESTS_3rdparty_INCLUDES})
    set(TESTS_3rdparty_DEPENDEND_LIBS ${TESTS_3rdparty_DEPENDEND_LIBS})

    message(STATUS "--------------------------------------------------------------------------------")
    message(STATUS "Google Testing Framework")

    message(STATUS "  googletest")
    set(gtest_SOURCE_NAME "release-1.7.0.tar.gz")
    set(gtest_SOURCE_MD5 "4ff6353b2560df0afecfbda3b2763847")
    message(STATUS "    using version: 1.7.0  (file: ${gtest_SOURCE_NAME})")

    if(NOT gtest_SOURCES)
        set(gtest_SOURCE_LOCATION "https://github.com/google/googletest/archive")
        message(STATUS "    going to download from: ${gtest_SOURCE_LOCATION}")
        message(STATUS "      HINT: specify -Dgtest_SOURCES=/local/directory/with/${gtest_SOURCE_NAME} to avoid redownloading")
    else()
        set(gtest_SOURCE_LOCATION ${gtest_SOURCES})
        message(STATUS "    using pre-downloaded sources at: ${gtest_SOURCE_LOCATION}")
    endif()


    message(STATUS "  googlemock")
    set(gmock_SOURCE_NAME "release-1.7.0.tar.gz")
    set(gmock_SOURCE_MD5 "13c3b4a57ad575763deb73fc0ad96e07")
    message(STATUS "    using version: 1.7.0  (file: ${gmock_SOURCE_NAME})")

    if(NOT gmock_SOURCES)
        set(gmock_SOURCE_LOCATION "https://github.com/google/googlemock/archive")
        message(STATUS "    going to download from: ${gmock_SOURCE_LOCATION}")
        message(STATUS "      HINT: specify -Dgmock_SOURCES=/local/directory/with/${gmock_SOURCE_NAME} to avoid redownloading")
    else()
        set(gmock_SOURCE_LOCATION ${gmock_SOURCES})
        message(STATUS "    using pre-downloaded sources at: ${gmock_SOURCE_LOCATION}")
    endif()

    set(gmock_PREFIX "${3rdparty_PREFIX}/gmock")
    set(gmock_TMP_PREFIX "${3rdparty_TMP_PREFIX}/gmock")

    set(gmock_TMP_DIR "${gmock_TMP_PREFIX}/tmp")
    set(gmock_STAMP_DIR "${gmock_TMP_PREFIX}/stamp")
    set(gmock_DOWNLOAD_DIR "${gmock_PREFIX}/download")
    set(gmock_SOURCE_DIR "${gmock_PREFIX}/src")
    set(gmock_BINARY_DIR "${gmock_TMP_PREFIX}/build")

    set(gtest_PREFIX "${3rdparty_PREFIX}/gtest")
    set(gtest_TMP_PREFIX "${3rdparty_TMP_PREFIX}/gtest")

    set(gtest_TMP_DIR "${gtest_TMP_PREFIX}/tmp")
    set(gtest_STAMP_DIR "${gtest_TMP_PREFIX}/stamp")
    set(gtest_DOWNLOAD_DIR "${gtest_PREFIX}/download")
    # gtest sources should be in a subdir of gmock sources
    set(gtest_SOURCE_DIR "${gmock_PREFIX}/gtest")
    # and gtest is build by gmock
    set(gtest_BINARY_DIR "${gmock_BINARY_DIR}/gtest")

    # Add gtest -- only download
    ExternalProject_Add(gtest
        TMP_DIR "${gtest_TMP_DIR}"
        STAMP_DIR "${gtest_STAMP_DIR}"
        DOWNLOAD_DIR "${gtest_DOWNLOAD_DIR}"
        SOURCE_DIR "${gtest_SOURCE_DIR}"
        BINARY_DIR "${gtest_BINARY_DIR}"
        INSTALL_DIR ""

        EXCLUDE_FROM_ALL ON

        URL ${gtest_SOURCE_LOCATION}/${gtest_SOURCE_NAME}
        URL_MD5 ${gtest_SOURCE_MD5}
        TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}

        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""

        LOG_DOWNLOAD ON
    )

    # Add gmock
    ExternalProject_Add(gmock
        TMP_DIR "${gmock_TMP_DIR}"
        STAMP_DIR "${gmock_STAMP_DIR}"
        DOWNLOAD_DIR "${gmock_DOWNLOAD_DIR}"
        SOURCE_DIR "${gmock_SOURCE_DIR}"
        BINARY_DIR "${gmock_BINARY_DIR}"
        INSTALL_DIR ""

        DEPENDS gtest
        EXCLUDE_FROM_ALL ON

        URL ${gmock_SOURCE_LOCATION}/${gmock_SOURCE_NAME}
        URL_MD5 ${gmock_SOURCE_MD5}
        TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}

        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_C_COMPILE=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
            -Dgtest_USE_OWN_TR1_TUPLE=ON
            -Dgtest_force_shared_crt=ON
            -Dgtest_build_tests=OFF
        BUILD_COMMAND make -j4
        INSTALL_COMMAND ""

        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON
    )

    set(gtest_INCLUDE_DIR ${gtest_SOURCE_DIR}/include)
    set(gtest_LIBRARY_DIR ${gtest_BINARY_DIR})
    set(gmock_INCLUDE_DIR ${gmock_SOURCE_DIR}/include)
    set(gmock_LIBRARY_DIR ${gmock_BINARY_DIR})

    add_library(gtest_lib STATIC IMPORTED GLOBAL)
    set_property(TARGET gtest_lib
        PROPERTY IMPORTED_LOCATION ${gtest_LIBRARY_DIR}/libgtest.a
    )
    add_dependencies(gtest_lib gmock)

    add_library(gmock_lib STATIC IMPORTED GLOBAL)
    set_property(TARGET gmock_lib
        PROPERTY IMPORTED_LOCATION ${gmock_LIBRARY_DIR}/libgmock.a
    )
    add_dependencies(gmock_lib gmock)

    message(STATUS "  googletest")
    message(STATUS "    include path: ${gtest_INCLUDE_DIR}")
    message(STATUS "    library path: ${gtest_LIBRARY_DIR}")
    message(STATUS "  googlemock")
    message(STATUS "    include path: ${gmock_INCLUDE_DIR}")
    message(STATUS "    library path: ${gmock_LIBRARY_DIR}")

    set(GOOGLETEST_INCLUDE_DIRS ${gtest_INCLUDE_DIR} ${gmock_INCLUDE_DIR})
    set(GOOGLETEST_LIBRARY_DIRS ${gtest_LIBRARY_DIR} ${gmock_LIBRARY_DIR})

    set(CMAKE_THREAD_PREFER_PTHREAD)
    find_package(Threads QUIET)

    list(APPEND TESTS_3rdparty_INCLUDES ${GOOGLETEST_INCLUDE_DIRS})
    list(APPEND TESTS_3rdparty_DEPENDEND_LIBS gtest_lib gmock_lib ${CMAKE_THREAD_LIBS_INIT})
endif()

# propagate include lists to parent directory
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS} PARENT_SCOPE)
set(3rdparty_INCLUDES ${3rdparty_INCLUDES} PARENT_SCOPE)
set(pfasst_DEPENDEND_TARGETS ${pfasst_DEPENDEND_TARGETS} PARENT_SCOPE)
set(pfasst_TESTS_DEPENDEND_TARGETS ${pfasst_TESTS_DEPENDEND_TARGETS} PARENT_SCOPE)
set(TESTS_3rdparty_DEPENDEND_LIBS ${TESTS_3rdparty_DEPENDEND_LIBS} PARENT_SCOPE)
set(TESTS_3rdparty_INCLUDES ${TESTS_3rdparty_INCLUDES} PARENT_SCOPE)
