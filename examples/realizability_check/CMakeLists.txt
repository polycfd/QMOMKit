cmake_minimum_required(VERSION 3.16)

project(QMOMKitRealizabilityCheckExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(
    realizability_check
    main.cpp
)
target_link_libraries(
    realizability_check
    PRIVATE qmomkit::qmomkit
)
