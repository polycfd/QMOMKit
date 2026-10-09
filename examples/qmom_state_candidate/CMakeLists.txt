cmake_minimum_required(VERSION 3.16)

project(QMOMKitQmomStateCandidateExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(
    qmom_state_candidate_example
    main.cpp
)
target_link_libraries(
    qmom_state_candidate_example
    PRIVATE qmomkit::qmomkit
)
