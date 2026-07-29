cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom2DVanishingVarianceExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

add_executable(cqmom_2d_vanishing_variance_example main.cpp physics.cpp sampling.cpp utilities.cpp)
target_link_libraries(
    cqmom_2d_vanishing_variance_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(cqmom_2d_vanishing_variance_mc_example main_mc.cpp particles_mc.cpp physics.cpp sampling.cpp utilities.cpp)
target_link_libraries(
    cqmom_2d_vanishing_variance_mc_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

configure_file(configuration.qmom configuration.qmom COPYONLY)
