cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmomThermalSphericalBreakupExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

set(CQMOM_THERMAL_SPHERICAL_BREAKUP_SOURCES
    physics.cpp
    sampling.cpp
    utilities.cpp
)

add_executable(cqmom_thermal_spherical_breakup_example main.cpp ${CQMOM_THERMAL_SPHERICAL_BREAKUP_SOURCES})
target_link_libraries(cqmom_thermal_spherical_breakup_example PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)

add_executable(cqmom_thermal_spherical_breakup_mc_example main_mc.cpp ${CQMOM_THERMAL_SPHERICAL_BREAKUP_SOURCES})
target_link_libraries(cqmom_thermal_spherical_breakup_mc_example PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)

configure_file(configuration.qmom configuration.qmom COPYONLY)
