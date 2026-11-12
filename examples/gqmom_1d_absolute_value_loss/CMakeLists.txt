cmake_minimum_required(VERSION 3.16)

project(QMOMKitGqmom1DAbsoluteValueLossExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(gqmom_1d_absolute_value_loss_example main.cpp physics.cpp)
target_link_libraries(gqmom_1d_absolute_value_loss_example PRIVATE qmomkit::qmomkit)

configure_file(configuration.qmom configuration.qmom COPYONLY)
