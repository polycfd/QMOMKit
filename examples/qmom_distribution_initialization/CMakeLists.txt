cmake_minimum_required(VERSION 3.16)

project(QMOMKitQmomDistributionInitializationExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(
    qmom_distribution_initialization_example
    main.cpp
)
target_link_libraries(
    qmom_distribution_initialization_example
    PRIVATE qmomkit::qmomkit
)
