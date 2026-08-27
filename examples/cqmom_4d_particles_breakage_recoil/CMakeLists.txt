cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom4DParticlesBreakageRecoilExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

set(CQMOM_4D_PARTICLES_BREAKAGE_RECOIL_SOURCES
    physics.cpp
    sampling.cpp
    utilities.cpp
)

add_executable(cqmom_4d_particles_breakage_recoil_example main.cpp ${CQMOM_4D_PARTICLES_BREAKAGE_RECOIL_SOURCES})
target_link_libraries(cqmom_4d_particles_breakage_recoil_example PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)

add_executable(cqmom_4d_particles_breakage_recoil_mc_example main_mc.cpp ${CQMOM_4D_PARTICLES_BREAKAGE_RECOIL_SOURCES})
target_link_libraries(cqmom_4d_particles_breakage_recoil_mc_example PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo)

configure_file(configuration.qmom configuration.qmom COPYONLY)
