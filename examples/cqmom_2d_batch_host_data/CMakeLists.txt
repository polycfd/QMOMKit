cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom2DBatchHostDataExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(cqmom_2d_batch_host_data_example main.cpp host_data.cpp physics.cpp utilities.cpp)
target_link_libraries(
    cqmom_2d_batch_host_data_example
    PRIVATE
        qmomkit::qmomkit
)

configure_file(configuration.qmom configuration.qmom COPYONLY)
