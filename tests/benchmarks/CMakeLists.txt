cmake_minimum_required(VERSION 3.16)

project(QMOMKitBenchmarks LANGUAGES CXX)

find_package(QMOMKit CONFIG REQUIRED)

set(QMOMKIT_REPOSITORY_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")

include(FetchContent)

find_package(benchmark CONFIG QUIET)

if(NOT TARGET benchmark::benchmark_main)
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_INSTALL_DOCS OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
        google_benchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG v1.8.5
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(google_benchmark)
endif()

function(qmomkit_add_benchmark target_name source_file)
    add_executable(${target_name} ${source_file})
    target_include_directories(${target_name} PRIVATE "${QMOMKIT_REPOSITORY_ROOT}/src")
    target_link_libraries(${target_name} PRIVATE qmomkit::qmomkit benchmark::benchmark_main)
    set_target_properties(${target_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()

qmomkit_add_benchmark(wheeler_benchmark wheeler_benchmark.cpp)
qmomkit_add_benchmark(vandermonde_benchmark vandermonde_benchmark.cpp)
qmomkit_add_benchmark(cqmom_benchmark cqmom_benchmark.cpp)
qmomkit_add_benchmark(performance_instrumentation_benchmark performance_instrumentation_benchmark.cpp)
