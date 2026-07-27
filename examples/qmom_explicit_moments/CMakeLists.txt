cmake_minimum_required(VERSION 3.16)

project(QMOMKitQmomExplicitMomentsExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(
    qmom_explicit_moments_example
    main.cpp
)
target_link_libraries(
    qmom_explicit_moments_example
    PRIVATE qmomkit::qmomkit
)
