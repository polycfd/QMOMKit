cmake_minimum_required(VERSION 3.16)

project(QMOMKitQmom1DMinimalExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(qmom_1d_minimal_example main.cpp)
target_link_libraries(qmom_1d_minimal_example PRIVATE qmomkit::qmomkit)

configure_file(configuration.qmom configuration.qmom COPYONLY)
