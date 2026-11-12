cmake_minimum_required(VERSION 3.16)

project(QMOMKitGqmomTrimodalDistributionExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(
    gqmom_trimodal_distribution_example
    main.cpp
)
target_link_libraries(
    gqmom_trimodal_distribution_example
    PRIVATE qmomkit::qmomkit
)

configure_file(configuration.qmom configuration.qmom COPYONLY)
