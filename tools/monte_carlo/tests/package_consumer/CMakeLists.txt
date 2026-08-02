cmake_minimum_required(VERSION 3.16)

project(qmomkit_monte_carlo_consumer LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

add_executable(consumer main.cpp)

target_link_libraries(
    consumer
    PRIVATE qmomkit_monte_carlo::qmomkit_monte_carlo
)
