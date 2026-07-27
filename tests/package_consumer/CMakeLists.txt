cmake_minimum_required(VERSION 3.16)

project(qmomkit_consumer LANGUAGES CXX)

find_package(QMOMKit CONFIG REQUIRED)

add_executable(consumer main.cpp)

target_link_libraries(
    consumer
    PRIVATE qmomkit::qmomkit
)
