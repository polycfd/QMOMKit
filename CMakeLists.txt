cmake_minimum_required(VERSION 3.16)

project(QMOMKit VERSION 0.1.0 LANGUAGES CXX)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "Install path prefix" FORCE)
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION "Enable internal QMOMKit performance instrumentation counters" OFF)

add_library(
    qmomkit
    src/conditional_moments.cpp
    src/conditioning_path.cpp
    src/configuration.cpp
    src/continuous_moment_rates.cpp
    src/distributions.cpp
    src/internal/linear_algebra.cpp
    src/internal/vandermonde.cpp
    src/moment_source.cpp
    src/node_advection.cpp
    src/output.cpp
    src/public_api.cpp
    src/moment_reconstruction.cpp
    src/realizability.cpp
    src/recursive_quadrature.cpp
    src/runge_kutta.cpp
    src/transport_integrator.cpp
    src/wheeler.cpp
)
add_library(qmomkit::qmomkit ALIAS qmomkit)

target_compile_features(qmomkit PUBLIC cxx_std_17)
target_include_directories(
    qmomkit
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_compile_definitions(
    qmomkit
    PUBLIC
        $<$<BOOL:${QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION}>:QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION>
)

if(MSVC)
    target_compile_options(qmomkit PRIVATE /W4 /permissive-)
else()
    target_compile_options(qmomkit PRIVATE -Wall -Wextra -Wpedantic)
endif()

install(
    TARGETS qmomkit
    EXPORT QMOMKitTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    FILES
        include/qmomkit/conditional_moments.hpp
        include/qmomkit/conditioning_path.hpp
        include/qmomkit/configuration.hpp
        include/qmomkit/cqmom.hpp
        include/qmomkit/distributions.hpp
        include/qmomkit/moment_set.hpp
        include/qmomkit/moment_source.hpp
        include/qmomkit/multi_index.hpp
        include/qmomkit/node_advection.hpp
        include/qmomkit/output.hpp
        include/qmomkit/phase_space_hierarchy.hpp
        include/qmomkit/qmom.hpp
        include/qmomkit/qmomkit.hpp
        include/qmomkit/quadrature.hpp
        include/qmomkit/moment_reconstruction.hpp
        include/qmomkit/realizability.hpp
        include/qmomkit/realizability_policy.hpp
        include/qmomkit/runge_kutta.hpp
        include/qmomkit/transport_integrator.hpp
        include/qmomkit/vandermonde_solver.hpp
        include/qmomkit/wheeler.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qmomkit
)

install(
    EXPORT QMOMKitTargets
    FILE QMOMKitTargets.cmake
    NAMESPACE qmomkit::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QMOMKit
)

configure_package_config_file(
    cmake/QMOMKitConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/QMOMKitConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QMOMKit
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/QMOMKitConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/QMOMKitConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/QMOMKitConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/QMOMKit
)

message(STATUS "QMOMKit build type: ${CMAKE_BUILD_TYPE}")
