cmake_minimum_required(VERSION 3.16)

project(QMOMKitQmom1DAerosolEvaporationExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(qmom_1d_aerosol_evaporation_example main.cpp physics.cpp)
target_link_libraries(qmom_1d_aerosol_evaporation_example PRIVATE qmomkit::qmomkit)

configure_file(configuration.qmom configuration.qmom COPYONLY)
