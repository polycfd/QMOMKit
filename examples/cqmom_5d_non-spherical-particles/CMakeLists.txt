cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom5dNonSphericalParticlesExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

add_executable(
    cqmom_5d_non_spherical_particles_example
    main.cpp
    particle_model.cpp
    breakage_model.cpp
    sampling.cpp
    utilities.cpp
)
target_link_libraries(
    cqmom_5d_non_spherical_particles_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(
    cqmom_5d_non_spherical_particles_mc_example
    main_mc.cpp
    particle_model.cpp
    breakage_model.cpp
    particles_mc.cpp
    sampling.cpp
    utilities.cpp
)
target_link_libraries(
    cqmom_5d_non_spherical_particles_mc_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)
