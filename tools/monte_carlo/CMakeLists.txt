cmake_minimum_required(VERSION 3.16)

project(
    QMOMKitMonteCarlo
    VERSION 0.1.0
    LANGUAGES CXX
)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(
        CMAKE_INSTALL_PREFIX
        "${CMAKE_SOURCE_DIR}/install"
        CACHE PATH
        "Install path prefix"
        FORCE
    )
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
option(BUILD_TESTING "Build QMOMKitMonteCarlo tests" OFF)
include(CTest)

find_package(QMOMKit CONFIG REQUIRED)

add_library(
    qmomkit_monte_carlo
    src/birth.cpp
    src/evolution.cpp
    src/moments.cpp
    src/output.cpp
    src/sampling.cpp
    src/sampled_initialization.cpp
)
add_library(qmomkit_monte_carlo::qmomkit_monte_carlo ALIAS qmomkit_monte_carlo)

target_compile_features(qmomkit_monte_carlo PUBLIC cxx_std_17)
target_include_directories(
    qmomkit_monte_carlo
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/qmomkit_monte_carlo>
)
target_link_libraries(
    qmomkit_monte_carlo
    PUBLIC
        qmomkit::qmomkit
)

if(MSVC)
    target_compile_options(qmomkit_monte_carlo PRIVATE /W4 /permissive-)
else()
    target_compile_options(qmomkit_monte_carlo PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(BUILD_TESTING)
    # Each logical case is registered with CTest on its own, so a failure names
    # the behavior that broke instead of the file it lives in. Selectors are
    # discovered from the `void test_*()` definitions, so adding a case to a
    # source file is enough to register it.
    function(monte_carlo_add_test_executable target_name source_file)
        add_executable(${target_name} ${source_file})
        target_link_libraries(${target_name} PRIVATE qmomkit_monte_carlo qmomkit::qmomkit)
        target_include_directories(${target_name} PRIVATE "${CMAKE_CURRENT_LIST_DIR}/tests")

        string(REGEX REPLACE "^monte_carlo_(.*)_tests$" "\\1" test_subsystem "${target_name}")
        file(STRINGS "${source_file}" case_definitions REGEX "^[ \t]*(static[ \t]+)?void[ \t]+test_[A-Za-z0-9_]+[ \t]*\\([ \t]*\\)")
        if(case_definitions STREQUAL "")
            message(FATAL_ERROR "No test cases discovered in ${source_file}")
        endif()
        foreach(definition IN LISTS case_definitions)
            string(REGEX REPLACE "^[ \t]*(static[ \t]+)?void[ \t]+test_([A-Za-z0-9_]+)[ \t]*\\(.*" "\\2" selector "${definition}")
            add_test(NAME monte_carlo_${test_subsystem}_${selector} COMMAND ${target_name} ${selector})
            set_tests_properties(monte_carlo_${test_subsystem}_${selector} PROPERTIES LABELS monte_carlo)
        endforeach()
    endfunction()

    monte_carlo_add_test_executable(monte_carlo_population_birth_tests tests/population_birth.cpp)
    monte_carlo_add_test_executable(monte_carlo_sampling_moment_verification_tests tests/sampling_moment_verification.cpp)
    monte_carlo_add_test_executable(monte_carlo_ensemble_evolution_tests tests/ensemble_evolution.cpp)
    monte_carlo_add_test_executable(monte_carlo_weighted_point_evolution_tests tests/weighted_point_evolution.cpp)
    monte_carlo_add_test_executable(monte_carlo_sampled_initialization_tests tests/sampled_initialization.cpp)
    monte_carlo_add_test_executable(monte_carlo_selective_moments_tests tests/selective_moments.cpp)
    monte_carlo_add_test_executable(monte_carlo_output_tests tests/output.cpp)
endif()

install(
    TARGETS qmomkit_monte_carlo
    EXPORT QMOMKitMonteCarloTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
    DIRECTORY include/monte_carlo
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qmomkit_monte_carlo
    FILES_MATCHING PATTERN "*.hpp"
)

install(
    EXPORT QMOMKitMonteCarloTargets
    FILE QMOMKitMonteCarloTargets.cmake
    NAMESPACE qmomkit_monte_carlo::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QMOMKitMonteCarlo
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/QMOMKitMonteCarloConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/QMOMKitMonteCarloConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QMOMKitMonteCarlo
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/QMOMKitMonteCarloConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/QMOMKitMonteCarloConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/QMOMKitMonteCarloConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QMOMKitMonteCarlo
)
