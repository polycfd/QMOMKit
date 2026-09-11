cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom2DMinimalExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(cqmom_2d_minimal_example main.cpp)
target_link_libraries(cqmom_2d_minimal_example PRIVATE qmomkit::qmomkit)
