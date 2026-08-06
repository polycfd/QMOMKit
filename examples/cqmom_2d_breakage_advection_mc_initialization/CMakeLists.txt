cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom2DAdvectionBreakageMcInitializationExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

add_executable(cqmom_2d_breakage_advection_mc_initialization_example main.cpp physics.cpp sampling.cpp utilities.cpp)
target_link_libraries(
    cqmom_2d_breakage_advection_mc_initialization_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(cqmom_2d_breakage_advection_mc_initialization_mc_example main_mc.cpp particles_mc.cpp physics.cpp sampling.cpp utilities.cpp)
target_link_libraries(
    cqmom_2d_breakage_advection_mc_initialization_mc_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

configure_file(configuration.qmom configuration.qmom COPYONLY)
