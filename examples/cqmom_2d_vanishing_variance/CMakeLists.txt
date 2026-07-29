cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom2DVanishingVarianceExample LANGUAGES CXX)

find_package(QMOMKit CONFIG REQUIRED)

add_executable(cqmom_2d_vanishing_variance_example main.cpp physics.cpp utilities.cpp)
target_link_libraries(
    cqmom_2d_vanishing_variance_example
    PRIVATE
        qmomkit::qmomkit
)

configure_file(configuration.qmom configuration.qmom COPYONLY)
