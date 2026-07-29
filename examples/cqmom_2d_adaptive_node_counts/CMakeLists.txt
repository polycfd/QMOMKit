cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom2DAdaptiveNodeCountsExample LANGUAGES CXX)

include(CTest)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(cqmom_2d_adaptive_node_counts_example main.cpp physics.cpp utilities.cpp)
target_link_libraries(
    cqmom_2d_adaptive_node_counts_example
    PRIVATE
        qmomkit::qmomkit
)

configure_file(configuration.qmom configuration.qmom COPYONLY)

if(QMOMKIT_BUILD_TESTS OR BUILD_TESTING)
    set(adaptive_node_counts_test_configs default configuration.qmom)
    while(adaptive_node_counts_test_configs)
        list(POP_FRONT adaptive_node_counts_test_configs test_name test_config)
        set(test_directory "${CMAKE_CURRENT_BINARY_DIR}/test_${test_name}")
        file(MAKE_DIRECTORY "${test_directory}")
        configure_file("${test_config}" "test_${test_name}/configuration.qmom" COPYONLY)
        add_test(
            NAME cqmom_2d_adaptive_node_counts_${test_name}
            COMMAND cqmom_2d_adaptive_node_counts_example
            WORKING_DIRECTORY "${test_directory}"
        )
    endwhile()
endif()
