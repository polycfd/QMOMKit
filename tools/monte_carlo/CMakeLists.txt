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
    add_executable(monte_carlo_population_birth_tests tests/population_birth.cpp)
    target_link_libraries(monte_carlo_population_birth_tests PRIVATE qmomkit_monte_carlo qmomkit::qmomkit)
    add_test(NAME monte_carlo_population_birth COMMAND monte_carlo_population_birth_tests)
    set_tests_properties(monte_carlo_population_birth PROPERTIES LABELS monte_carlo)

    add_executable(monte_carlo_sampling_moment_validation_tests tests/sampling_moment_validation.cpp)
    target_link_libraries(monte_carlo_sampling_moment_validation_tests PRIVATE qmomkit_monte_carlo qmomkit::qmomkit)
    add_test(NAME monte_carlo_sampling_moment_validation COMMAND monte_carlo_sampling_moment_validation_tests)
    set_tests_properties(monte_carlo_sampling_moment_validation PROPERTIES LABELS monte_carlo)

    add_executable(monte_carlo_ensemble_evolution_tests tests/ensemble_evolution.cpp)
    target_link_libraries(monte_carlo_ensemble_evolution_tests PRIVATE qmomkit_monte_carlo qmomkit::qmomkit)
    add_test(NAME monte_carlo_ensemble_evolution COMMAND monte_carlo_ensemble_evolution_tests)
    set_tests_properties(monte_carlo_ensemble_evolution PROPERTIES LABELS monte_carlo)

    add_executable(monte_carlo_quadrature_ensemble_comparison_tests tests/quadrature_ensemble_comparison.cpp)
    target_link_libraries(monte_carlo_quadrature_ensemble_comparison_tests PRIVATE qmomkit_monte_carlo qmomkit::qmomkit)
    add_test(NAME monte_carlo_quadrature_ensemble_comparison COMMAND monte_carlo_quadrature_ensemble_comparison_tests)
    set_tests_properties(monte_carlo_quadrature_ensemble_comparison PROPERTIES LABELS monte_carlo)

    add_executable(monte_carlo_sampled_initialization_tests tests/sampled_initialization.cpp)
    target_link_libraries(monte_carlo_sampled_initialization_tests PRIVATE qmomkit_monte_carlo qmomkit::qmomkit)
    add_test(NAME monte_carlo_sampled_initialization COMMAND monte_carlo_sampled_initialization_tests)
    set_tests_properties(monte_carlo_sampled_initialization PROPERTIES LABELS monte_carlo)

    add_executable(monte_carlo_selective_moments_tests tests/selective_moments.cpp)
    target_link_libraries(monte_carlo_selective_moments_tests PRIVATE qmomkit_monte_carlo qmomkit::qmomkit)
    add_test(NAME monte_carlo_selective_moments COMMAND monte_carlo_selective_moments_tests)
    set_tests_properties(monte_carlo_selective_moments PROPERTIES LABELS monte_carlo)

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
