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
    src/comparison.cpp
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

add_executable(monte_carlo_sample_based_initialization_example examples/sample_based_initialization.cpp)
target_link_libraries(monte_carlo_sample_based_initialization_example PRIVATE qmomkit_monte_carlo qmomkit::qmomkit)

if(BUILD_TESTING)
    add_executable(monte_carlo_mc1_tests tests/mc1_validation.cpp)
    target_link_libraries(monte_carlo_mc1_tests PRIVATE qmomkit_monte_carlo qmomkit::qmomkit)
    add_test(NAME monte_carlo_mc1_validation COMMAND monte_carlo_mc1_tests)
    set_tests_properties(monte_carlo_mc1_validation PROPERTIES LABELS monte_carlo)

    add_executable(monte_carlo_mc2a_tests tests/mc2a_evolution.cpp)
    target_link_libraries(monte_carlo_mc2a_tests PRIVATE qmomkit_monte_carlo qmomkit::qmomkit)
    add_test(NAME monte_carlo_mc2a_evolution COMMAND monte_carlo_mc2a_tests)
    set_tests_properties(monte_carlo_mc2a_evolution PROPERTIES LABELS monte_carlo)

    add_executable(monte_carlo_mc2b_tests tests/mc2b_comparison.cpp)
    target_link_libraries(monte_carlo_mc2b_tests PRIVATE qmomkit_monte_carlo qmomkit::qmomkit)
    add_test(NAME monte_carlo_mc2b_comparison COMMAND monte_carlo_mc2b_tests)
    set_tests_properties(monte_carlo_mc2b_comparison PROPERTIES LABELS monte_carlo)

    add_executable(monte_carlo_sampled_initialization_tests tests/sampled_initialization.cpp)
    target_link_libraries(monte_carlo_sampled_initialization_tests PRIVATE qmomkit_monte_carlo qmomkit::qmomkit)
    add_test(NAME monte_carlo_sampled_initialization COMMAND monte_carlo_sampled_initialization_tests)
    set_tests_properties(monte_carlo_sampled_initialization PROPERTIES LABELS monte_carlo)

    add_executable(monte_carlo_node_vs_ensemble_example examples/comparison/node_vs_ensemble.cpp)
    target_link_libraries(monte_carlo_node_vs_ensemble_example PRIVATE qmomkit_monte_carlo qmomkit::qmomkit)
    add_test(NAME monte_carlo_node_vs_ensemble_example COMMAND monte_carlo_node_vs_ensemble_example)
    set_tests_properties(monte_carlo_node_vs_ensemble_example PROPERTIES LABELS monte_carlo)

    add_test(
        NAME monte_carlo_sample_based_initialization_qmom
        COMMAND
            monte_carlo_sample_based_initialization_example
            ${CMAKE_CURRENT_SOURCE_DIR}/configuration/qmom_distribution_lognormal.qmom
            ${CMAKE_CURRENT_BINARY_DIR}/sample_based_qmom
    )
    set_tests_properties(monte_carlo_sample_based_initialization_qmom PROPERTIES LABELS monte_carlo)
    add_test(
        NAME monte_carlo_sample_based_initialization_cqmom
        COMMAND
            monte_carlo_sample_based_initialization_example
            ${CMAKE_CURRENT_SOURCE_DIR}/configuration/cqmom_independent_product_2d.qmom
            ${CMAKE_CURRENT_BINARY_DIR}/sample_based_cqmom
    )
    set_tests_properties(monte_carlo_sample_based_initialization_cqmom PROPERTIES LABELS monte_carlo)
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
