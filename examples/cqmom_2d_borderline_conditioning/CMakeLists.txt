cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom2DBorderlineConditioningExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(cqmom_2d_borderline_conditioning_example main.cpp conditioning_workflow.cpp conditioning_report.cpp)
target_link_libraries(
    cqmom_2d_borderline_conditioning_example
    PRIVATE
        qmomkit::qmomkit
)
