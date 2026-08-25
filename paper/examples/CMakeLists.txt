cmake_minimum_required(VERSION 3.16)

project(QMOMKitPaperExamples LANGUAGES CXX)

find_package(QMOMKit CONFIG REQUIRED)
find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

add_executable(
    qmomkit_paper_1d_aerosols_evaporation
    1d_aerosols_evaporation/main.cpp
    1d_aerosols_evaporation/physics.cpp
)
target_link_libraries(
    qmomkit_paper_1d_aerosols_evaporation
    PRIVATE
        qmomkit::qmomkit
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(
    qmomkit_paper_2d_bubbles_cavitation
    2d_bubbles_cavitation/main.cpp
    2d_bubbles_cavitation/physics.cpp
    2d_bubbles_cavitation/particles_mc.cpp
    2d_bubbles_cavitation/sampling.cpp
    2d_bubbles_cavitation/utilities.cpp
)
target_link_libraries(
    qmomkit_paper_2d_bubbles_cavitation
    PRIVATE
        qmomkit::qmomkit
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(
    qmomkit_paper_3d_droplets_condensation_nucleation
    3d_droplets_condensation-nucleation/main.cpp
    3d_droplets_condensation-nucleation/droplet_model.cpp
    3d_droplets_condensation-nucleation/sampling.cpp
    3d_droplets_condensation-nucleation/utilities.cpp
)
target_link_libraries(
    qmomkit_paper_3d_droplets_condensation_nucleation
    PRIVATE
        qmomkit::qmomkit
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(
    qmomkit_paper_4d_particles_breakup_recoil
    4d_particles_breakup_recoil/main.cpp
    4d_particles_breakup_recoil/physics.cpp
    4d_particles_breakup_recoil/sampling.cpp
    4d_particles_breakup_recoil/utilities.cpp
)
target_link_libraries(
    qmomkit_paper_4d_particles_breakup_recoil
    PRIVATE
        qmomkit::qmomkit
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(
    qmomkit_paper_5d_particles_breakup_heat_transfer
    5d_particles_breakup_heat-transfer/main.cpp
    5d_particles_breakup_heat-transfer/physics.cpp
    5d_particles_breakup_heat-transfer/sampling.cpp
    5d_particles_breakup_heat-transfer/utilities.cpp
)
target_link_libraries(
    qmomkit_paper_5d_particles_breakup_heat_transfer
    PRIVATE
        qmomkit::qmomkit
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

configure_file(1d_aerosols_evaporation/configuration.qmom 1d_aerosols_evaporation/configuration.qmom COPYONLY)
configure_file(2d_bubbles_cavitation/configuration.qmom 2d_bubbles_cavitation/configuration.qmom COPYONLY)
configure_file(3d_droplets_condensation-nucleation/configuration.qmom 3d_droplets_condensation-nucleation/configuration.qmom COPYONLY)
configure_file(4d_particles_breakup_recoil/configuration.qmom 4d_particles_breakup_recoil/configuration.qmom COPYONLY)
configure_file(5d_particles_breakup_heat-transfer/configuration.qmom 5d_particles_breakup_heat-transfer/configuration.qmom COPYONLY)
