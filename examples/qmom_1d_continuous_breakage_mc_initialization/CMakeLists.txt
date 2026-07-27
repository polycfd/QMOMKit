cmake_minimum_required(VERSION 3.16)

project(QMOMKitQmom1DContinuousBreakageMcInitializationExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

add_executable(qmom_1d_continuous_breakage_mc_initialization_example main.cpp physics.cpp sampling.cpp utilities.cpp)
target_link_libraries(
    qmom_1d_continuous_breakage_mc_initialization_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(qmom_1d_continuous_breakage_mc_initialization_mc_example main_mc.cpp particles_mc.cpp physics.cpp sampling.cpp utilities.cpp)
target_link_libraries(
    qmom_1d_continuous_breakage_mc_initialization_mc_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

configure_file(configuration.qmom configuration.qmom COPYONLY)
