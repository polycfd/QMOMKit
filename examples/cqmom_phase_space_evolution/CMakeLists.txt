cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmomPhaseSpaceEvolutionExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(cqmom_phase_space_evolution_example main.cpp physics.cpp utilities.cpp)
target_link_libraries(cqmom_phase_space_evolution_example PRIVATE qmomkit::qmomkit)

configure_file(configuration.qmom configuration.qmom COPYONLY)
