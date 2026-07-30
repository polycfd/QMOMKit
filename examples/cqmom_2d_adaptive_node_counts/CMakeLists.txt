cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom2DAdaptiveNodeCountsExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(cqmom_2d_adaptive_node_counts_example main.cpp physics.cpp utilities.cpp)
target_link_libraries(
    cqmom_2d_adaptive_node_counts_example
    PRIVATE
        qmomkit::qmomkit
)

configure_file(configuration.qmom configuration.qmom COPYONLY)
