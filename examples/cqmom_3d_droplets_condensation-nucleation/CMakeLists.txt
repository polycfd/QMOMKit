cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom3DDropletsCondensationNucleationExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

add_executable(cqmom_3d_droplets_condensation_nucleation_example main.cpp sampling.cpp utilities.cpp droplet_model.cpp)
target_link_libraries(
    cqmom_3d_droplets_condensation_nucleation_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(cqmom_3d_droplets_condensation_nucleation_mc_example main_mc.cpp sampling.cpp utilities.cpp droplet_model.cpp)
target_link_libraries(
    cqmom_3d_droplets_condensation_nucleation_mc_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

configure_file(configuration.qmom configuration.qmom COPYONLY)
