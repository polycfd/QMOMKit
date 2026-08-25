cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom4dSphericalBreakupQuadratureExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

set(CQMOM_4D_SPHERICAL_BREAKUP_QUADRATURE_SOURCES
    physics.cpp
    sampling.cpp
    utilities.cpp
)

add_executable(cqmom_4d_spherical_breakup_quadrature_example main.cpp ${CQMOM_4D_SPHERICAL_BREAKUP_QUADRATURE_SOURCES})
target_link_libraries(cqmom_4d_spherical_breakup_quadrature_example PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)

add_executable(cqmom_4d_spherical_breakup_quadrature_mc_example main_mc.cpp ${CQMOM_4D_SPHERICAL_BREAKUP_QUADRATURE_SOURCES})
target_link_libraries(cqmom_4d_spherical_breakup_quadrature_mc_example PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)

configure_file(configuration.qmom configuration.qmom COPYONLY)
