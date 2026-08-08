cmake_minimum_required(VERSION 3.16)

project(QMOMKitMomentRateElementOperator LANGUAGES CXX)

find_package(QMOMKit CONFIG REQUIRED)

add_executable(moment_rate_element_operator main.cpp)

target_compile_features(moment_rate_element_operator PRIVATE cxx_std_17)
target_link_libraries(moment_rate_element_operator PRIVATE qmomkit::qmomkit)

if(MSVC)
    target_compile_options(moment_rate_element_operator PRIVATE /W4 /permissive-)
else()
    target_compile_options(moment_rate_element_operator PRIVATE -Wall -Wextra -Wpedantic)
endif()
