cmake_minimum_required(VERSION 3.16)

project(QMOMKitHyqmom1DDragBgkCollisionsExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(hyqmom_1d_drag_bgk_collisions_example main.cpp physics.cpp)
target_link_libraries(hyqmom_1d_drag_bgk_collisions_example PRIVATE qmomkit::qmomkit)

configure_file(configuration.qmom configuration.qmom COPYONLY)
