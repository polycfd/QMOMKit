cmake_minimum_required(VERSION 3.16)

project(QMOMKitQmom1DExplicitStrangSplittingExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(qmom_1d_explicit_strang_splitting_example main.cpp utilities.cpp)
target_link_libraries(qmom_1d_explicit_strang_splitting_example PRIVATE qmomkit::qmomkit)

configure_file(configuration.qmom configuration.qmom COPYONLY)
