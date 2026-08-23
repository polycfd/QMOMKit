cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmomSphericalBreakupQuadrature4dExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

set(CQMOM_SPHERICAL_BREAKUP_QUADRATURE_4D_SOURCES
    physics.cpp
    sampling.cpp
    utilities.cpp
)

add_executable(cqmom_spherical_breakup_quadrature_4d_example main.cpp ${CQMOM_SPHERICAL_BREAKUP_QUADRATURE_4D_SOURCES})
target_link_libraries(cqmom_spherical_breakup_quadrature_4d_example PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)

add_executable(cqmom_spherical_breakup_quadrature_4d_mc_example main_mc.cpp ${CQMOM_SPHERICAL_BREAKUP_QUADRATURE_4D_SOURCES})
target_link_libraries(cqmom_spherical_breakup_quadrature_4d_mc_example PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)

configure_file(configuration.qmom configuration.qmom COPYONLY)
