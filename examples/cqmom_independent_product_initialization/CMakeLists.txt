cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmomIndependentProductInitializationExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(cqmom_independent_product_initialization_example main.cpp)
target_link_libraries(
    cqmom_independent_product_initialization_example
    PRIVATE qmomkit::qmomkit
)

configure_file(configuration.qmom configuration.qmom COPYONLY)
