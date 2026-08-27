cmake_minimum_required(VERSION 3.16)

project(QMOMKitGalerkinElementCoupling LANGUAGES CXX)

find_package(QMOMKit CONFIG REQUIRED)

add_executable(galerkin_element_coupling main.cpp)

target_compile_features(galerkin_element_coupling PRIVATE cxx_std_17)
target_link_libraries(galerkin_element_coupling PRIVATE qmomkit::qmomkit)

if(MSVC)
    target_compile_options(galerkin_element_coupling PRIVATE /W4 /permissive-)
else()
    target_compile_options(galerkin_element_coupling PRIVATE -Wall -Wextra -Wpedantic)
endif()
