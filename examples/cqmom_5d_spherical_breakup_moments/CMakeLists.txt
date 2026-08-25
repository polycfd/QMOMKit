cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom5dSphericalBreakupMomentsExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

set(CQMOM_5D_SPHERICAL_BREAKUP_MOMENTS_SOURCES
    physics.cpp
    sampling.cpp
    utilities.cpp
)

add_executable(cqmom_5d_spherical_breakup_moments_example main.cpp ${CQMOM_5D_SPHERICAL_BREAKUP_MOMENTS_SOURCES})
target_link_libraries(cqmom_5d_spherical_breakup_moments_example PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)

add_executable(cqmom_5d_spherical_breakup_moments_mc_example main_mc.cpp ${CQMOM_5D_SPHERICAL_BREAKUP_MOMENTS_SOURCES})
target_link_libraries(cqmom_5d_spherical_breakup_moments_mc_example PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)

configure_file(configuration.qmom configuration.qmom COPYONLY)
