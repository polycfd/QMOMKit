cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmomSphericalBreakupMoments5dExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

set(CQMOM_SPHERICAL_BREAKUP_MOMENTS_5D_SOURCES
    physics.cpp
    sampling.cpp
    utilities.cpp
)

add_executable(cqmom_spherical_breakup_moments_5d_example main.cpp ${CQMOM_SPHERICAL_BREAKUP_MOMENTS_5D_SOURCES})
target_link_libraries(cqmom_spherical_breakup_moments_5d_example PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)

add_executable(cqmom_spherical_breakup_moments_5d_mc_example main_mc.cpp ${CQMOM_SPHERICAL_BREAKUP_MOMENTS_5D_SOURCES})
target_link_libraries(cqmom_spherical_breakup_moments_5d_mc_example PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)

configure_file(configuration.qmom configuration.qmom COPYONLY)
