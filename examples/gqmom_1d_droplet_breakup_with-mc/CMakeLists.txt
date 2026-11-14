cmake_minimum_required(VERSION 3.16)

project(QMOMKitGqmom1DDropletBreakupWithMcExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

add_executable(gqmom_1d_droplet_breakup_with_mc_example main.cpp physics.cpp sampling.cpp utilities.cpp)
target_link_libraries(
    gqmom_1d_droplet_breakup_with_mc_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(gqmom_1d_droplet_breakup_with_mc_mc_example main_mc.cpp particles_mc.cpp physics.cpp sampling.cpp utilities.cpp)
target_link_libraries(
    gqmom_1d_droplet_breakup_with_mc_mc_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

configure_file(configuration.qmom configuration.qmom COPYONLY)
