# Managing 3rd-party libraries and tools
set(3rdparty_INCLUDES ${3rdparty_INCLUDES})
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS})
set(pfasst_DEPENDEND_TARGETS ${pfasst_DEPENDEND_TARGETS})
set(pfasst_TESTS_DEPENDEND_TARGETS ${pfasst_TESTS_DEPENDEND_TARGETS})

set(3rdparty_DOWNLOAD_TIMEOUT "60")

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

        DEPENDS boost_build

        DOWNLOAD_COMMAND ""
        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        BUILD_IN_SOURCE ON
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${Boost_${component}_BUILD_CMD} stage
        TEST_COMMAND ""
        INSTALL_COMMAND ""

        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON
    )

    list(APPEND pfasst_DEPENDEND_TARGETS boost_${component})
    list(APPEND BOOST_LIBRARIES ${BOOST_LIBRARY_DIR}/libboost_${component}.a)
endforeach()

message(STATUS "Boost include path: ${BOOST_INCLUDE_DIR}")
message(STATUS "Boost library path: ${BOOST_LIBRARY_DIR}")
list(APPEND 3rdparty_INCLUDES ${BOOST_INCLUDE_DIR})
list(APPEND 3rdparty_DEPENDEND_LIBS ${BOOST_LIBRARIES})


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Leathers")
list(APPEND 3rdparty_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/leathers)


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Better Enums")
list(APPEND 3rdparty_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/better-enums)


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "Eigen3")

if(pfasst_SYSTEM_EIGEN)
    find_package(Eigen3)
endif()
set_package_properties(Eigen3
    PROPERTIES
      DESCRIPTION "linear algebra library"
      URL "http://eigen.tuxfamily.org/"
      TYPE OPTIONAL
      PURPOSE "matrix and vector types and arithmetics"
)

if(NOT EIGEN3_FOUND)
    if(pfasst_SYSTEM_EIGEN)
        message(FATAL "Eigen3 was not found but requested to only use a system's installation.")
    else()
        set(Eigen3_SOURCE_URL "http://bitbucket.org/eigen/eigen/get/3.2.2.tar.bz2")
        set(Eigen3_SOURCE_MD5 "fc2e814ae449d16b331f7e1f4e272bd3")

        message(STATUS "going to download from: ${Eigen3_SOURCE_URL}")

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

            URL ${Eigen3_SOURCE_URL}
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
            LOG_CONFIGURE OFF
            LOG_BUILD OFF
        )

        # Specify include dir
        set(Eigen3_INCLUDE_PATH "${Eigen3_SOURCE_DIR}")

        list(APPEND pfasst_DEPENDEND_TARGETS Eigen3)

        msg_not_installed("Eigen3")
    endif()
endif()

list(APPEND 3rdparty_INCLUDES "${Eigen3_INCLUDE_PATH}")


message(STATUS "--------------------------------------------------------------------------------")
message(STATUS "FFTW3")

if(pfasst_SYSTEM_FFTW)
    find_package(FFTW)
endif()
set_package_properties(FFTW
    PROPERTIES
      DESCRIPTION "FFTW3"
      URL "http://fftw.org/"
      TYPE OPTIONAL
      PURPOSE "fast Fourier transformation"
)

if(NOT FFTW_FOUND)
    if(pfasst_SYSTEM_FFTW)
        message(FATAL "FFTW was not found but requested to only use a system's installation.")
    else()
        set(fftw3_SOURCE_URL "http://fftw.org/fftw-3.3.4.tar.gz")
        set(fftw3_SOURCE_MD5 "2edab8c06b24feeb3b82bbb3ebf3e7b3")

        message(STATUS "going to download and compile from sources retrieved from: ${fftw3_SOURCE_URL}")

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

            URL ${fftw3_SOURCE_URL}
            URL_MD5 ${fftw3_SOURCE_MD5}
            TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}

            UPDATE_COMMAND ""
            PATCH_COMMAND ""
            BUILD_IN_SOURCE ON
            CONFIGURE_COMMAND ${fftw3_SOURCE_DIR}/configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CXXFLAGS=${CMAKE_CXX_FLAGS} --prefix=${fftw3_INSTALL_DIR} --libdir=${fftw3_INSTALL_DIR}/lib
            BUILD_COMMAND make
            TEST_COMMAND ""
            INSTALL_COMMAND make install

            LOG_DOWNLOAD ON
            LOG_CONFIGURE ON
            LOG_BUILD ON
        )
        # Specify include dir
        set(FFTW_INCLUDE_PATH ${fftw3_INSTALL_DIR}/include)
        set(FFTW_LIBRARIES ${fftw3_INSTALL_DIR}/lib/libfftw3.a)

        list(APPEND pfasst_DEPENDEND_LIBS fftw3)

        msg_not_installed("FFTW3")
    endif()
else()
    set(FFTW_INCLUDE_PATH ${FFTW_INCLUDE_PATH} PARENT_SCOPE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} PARENT_SCOPE)
endif()
set(FFTW_FOUND ${FFTW_FOUND} PARENT_SCOPE)

list(APPEND 3rdparty_INCLUDES "${FFTW_INCLUDE_PATH}")
list(APPEND 3rdparty_DEPENDEND_LIBS "${FFTW_LIBRARIES}")


if(pfasst_BUILD_TESTS)
    set(TESTS_3rdparty_INCLUDES ${TESTS_3rdparty_INCLUDES})
    set(TESTS_3rdparty_DEPENDEND_LIBS ${TESTS_3rdparty_DEPENDEND_LIBS})

    message(STATUS "--------------------------------------------------------------------------------")
    message(STATUS "Google Testing Framework (gtest & gmock)")

    set_package_properties(GMock
        PROPERTIES
          DESCRIPTION "GMock"
          URL "https://code.google.com/p/googlemock/"
          TYPE OPTIONAL
          PURPOSE "Google testing and mocking framework"
    )

    set(gmock_SOURCE_URL "http://googlemock.googlecode.com/files/gmock-1.7.0.zip")
    set(gmock_SOURCE_MD5 "073b984d8798ea1594f5e44d85b20d66")

    message(STATUS "going to download and compile from sources retrieved from: ${gmock_SOURCE_URL}")

    set(gmock_PREFIX "${3rdparty_PREFIX}/gmock")
    set(gmock_TMP_PREFIX "${3rdparty_TMP_PREFIX}/gmock")

    set(gmock_TMP_DIR "${gmock_TMP_PREFIX}/tmp")
    set(gmock_STAMP_DIR "${gmock_TMP_PREFIX}/stamp")
    set(gmock_DOWNLOAD_DIR "${gmock_PREFIX}/download")
    set(gmock_SOURCE_DIR "${gmock_PREFIX}/src")
    set(gmock_BINARY_DIR "${gmock_TMP_PREFIX}/build")
    set(gmock_INSTALL_DIR "${gmock_PREFIX}/install")

    # Add gmock
    ExternalProject_Add(gmock
        TMP_DIR "${gmock_TMP_DIR}"
        STAMP_DIR "${gmock_STAMP_DIR}"
        DOWNLOAD_DIR "${gmock_DOWNLOAD_DIR}"
        SOURCE_DIR "${gmock_SOURCE_DIR}"
        BINARY_DIR "${gmock_BINARY_DIR}"
        INSTALL_DIR "${gmock_INSTALL_DIR}"
        
        URL ${gmock_SOURCE_URL}
        URL_MD5 ${gmock_SOURCE_MD5}
        TIMEOUT ${3rdparty_DOWNLOAD_TIMEOUT}

        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        CMAKE_ARGS "-DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILE=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=ON -Dgtest_force_shared_crt=ON -Dgmock_build_tests=OFF"

        INSTALL_COMMAND ""

        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON
    )

    # Specify include dir
    list(APPEND TESTS_3rdparty_INCLUDES ${gmock_SOURCE_DIR}/include ${gmock_SOURCE_DIR}/gtest/include)

    set(Suffix ".a")
    set(Pthread "-pthread")
    list(APPEND TESTS_3rdparty_DEPENDEND_LIBS ${gmock_BINARY_DIR}/${CMAKE_FIND_LIBRARY_PREFIXES}gmock${Suffix})
    list(APPEND TESTS_3rdparty_DEPENDEND_LIBS ${Pthread})
    list(APPEND pfasst_TESTS_DEPENDEND_TARGETS gmock)
endif()

# propagate include lists to parent directory
set(3rdparty_DEPENDEND_LIBS ${3rdparty_DEPENDEND_LIBS} PARENT_SCOPE)
set(3rdparty_INCLUDES ${3rdparty_INCLUDES} PARENT_SCOPE)
set(pfasst_DEPENDEND_TARGETS ${pfasst_DEPENDEND_TARGETS} PARENT_SCOPE)
set(pfasst_TESTS_DEPENDEND_TARGETS ${pfasst_TESTS_DEPENDEND_TARGETS} PARENT_SCOPE)
set(TESTS_3rdparty_DEPENDEND_LIBS ${TESTS_3rdparty_DEPENDEND_LIBS} PARENT_SCOPE)
set(TESTS_3rdparty_INCLUDES ${TESTS_3rdparty_INCLUDES} PARENT_SCOPE)
