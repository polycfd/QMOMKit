cmake_minimum_required(VERSION 3.16)

project(QMOMKitHyqmomCharacteristicSpeedsExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(hyqmom_characteristic_speeds_example main.cpp)
target_link_libraries(hyqmom_characteristic_speeds_example PRIVATE qmomkit::qmomkit)
