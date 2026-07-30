cmake_minimum_required(VERSION 3.16)

project(QMOMKitMomentRateElementOperator LANGUAGES CXX)

find_package(QMOMKit CONFIG REQUIRED)

add_executable(moment_rate_element_operator_check correctness.cpp)
add_executable(moment_rate_element_operator_benchmark benchmark.cpp)

target_compile_features(moment_rate_element_operator_check PRIVATE cxx_std_17)
target_compile_features(moment_rate_element_operator_benchmark PRIVATE cxx_std_17)

target_link_libraries(moment_rate_element_operator_check PRIVATE qmomkit::qmomkit)
target_link_libraries(moment_rate_element_operator_benchmark PRIVATE qmomkit::qmomkit)

if(CMAKE_BUILD_TYPE)
    set(MOMENT_RATE_ELEMENT_OPERATOR_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
else()
    set(MOMENT_RATE_ELEMENT_OPERATOR_BUILD_TYPE "unknown")
endif()

target_compile_definitions(moment_rate_element_operator_benchmark
    PRIVATE
        MOMENT_RATE_ELEMENT_OPERATOR_COMPILER_ID="${CMAKE_CXX_COMPILER_ID}"
        MOMENT_RATE_ELEMENT_OPERATOR_COMPILER_VERSION="${CMAKE_CXX_COMPILER_VERSION}"
        MOMENT_RATE_ELEMENT_OPERATOR_BUILD_TYPE="${MOMENT_RATE_ELEMENT_OPERATOR_BUILD_TYPE}"
        MOMENT_RATE_ELEMENT_OPERATOR_SYSTEM_NAME="${CMAKE_SYSTEM_NAME}")

if(MSVC)
    target_compile_options(moment_rate_element_operator_check PRIVATE /W4 /permissive-)
    target_compile_options(moment_rate_element_operator_benchmark PRIVATE /W4 /permissive-)
else()
    target_compile_options(moment_rate_element_operator_check PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(moment_rate_element_operator_benchmark PRIVATE -Wall -Wextra -Wpedantic)
endif()
