cmake_minimum_required(VERSION 3.16)

project(qmomkit_consumer LANGUAGES CXX)

find_package(QMOMKit CONFIG REQUIRED)

add_executable(consumer main.cpp)

target_link_libraries(
    consumer
    PRIVATE qmomkit::qmomkit
)

# Building this library is the check: it links the installed static library
# into a shared one.
add_library(shared_library_consumer SHARED shared_library.cpp)

target_link_libraries(
    shared_library_consumer
    PRIVATE qmomkit::qmomkit
)
