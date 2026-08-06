cmake_minimum_required(VERSION 3.16)

project(QMOMKitPaperExamples LANGUAGES CXX)

find_package(QMOMKit CONFIG REQUIRED)
find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

add_executable(
    qmomkit_paper_non_spherical_particles
    non_spherical_particles/main.cpp
    non_spherical_particles/particle_model.cpp
    non_spherical_particles/breakage_model.cpp
    non_spherical_particles/particles_mc.cpp
    non_spherical_particles/sampling.cpp
    non_spherical_particles/utilities.cpp
)
target_link_libraries(
    qmomkit_paper_non_spherical_particles
    PRIVATE
        qmomkit::qmomkit
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(
    qmomkit_paper_condensing_droplets
    condensing_droplets/main.cpp
    condensing_droplets/droplet_model.cpp
    condensing_droplets/sampling.cpp
    condensing_droplets/utilities.cpp
)
target_link_libraries(
    qmomkit_paper_condensing_droplets
    PRIVATE
        qmomkit::qmomkit
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(
    qmomkit_paper_polydisperse_cavitation
    polydisperse_cavitation/main.cpp
    polydisperse_cavitation/physics.cpp
    polydisperse_cavitation/particles_mc.cpp
    polydisperse_cavitation/sampling.cpp
    polydisperse_cavitation/utilities.cpp
)
target_link_libraries(
    qmomkit_paper_polydisperse_cavitation
    PRIVATE
        qmomkit::qmomkit
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

configure_file(non_spherical_particles/configuration.qmom non_spherical_particles/configuration.qmom COPYONLY)
configure_file(condensing_droplets/configuration.qmom condensing_droplets/configuration.qmom COPYONLY)
configure_file(polydisperse_cavitation/configuration.qmom polydisperse_cavitation/configuration.qmom COPYONLY)
