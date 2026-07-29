cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom2DVarianceReinjectionExample LANGUAGES CXX)

if(NOT TARGET qmomkit::qmomkit)
    find_package(QMOMKit CONFIG REQUIRED)
endif()

add_executable(cqmom_2d_variance_reinjection_example main.cpp physics.cpp utilities.cpp)
target_link_libraries(
    cqmom_2d_variance_reinjection_example
    PRIVATE
        qmomkit::qmomkit
)

configure_file(configuration.qmom configuration.qmom COPYONLY)
configure_file(configuration_swapped.qmom configuration_swapped.qmom COPYONLY)
configure_file(configuration_velocity3.qmom configuration_velocity3.qmom COPYONLY)
configure_file(configuration_velocity3_swapped.qmom configuration_velocity3_swapped.qmom COPYONLY)

if(QMOMKIT_BUILD_TESTS OR BUILD_TESTING)
    set(variance_reinjection_test_configs
        default configuration.qmom
        swapped configuration_swapped.qmom
        velocity3 configuration_velocity3.qmom
        velocity3_swapped configuration_velocity3_swapped.qmom
    )
    while(variance_reinjection_test_configs)
        list(POP_FRONT variance_reinjection_test_configs test_name test_config)
        set(test_directory "${CMAKE_CURRENT_BINARY_DIR}/test_${test_name}")
        file(MAKE_DIRECTORY "${test_directory}")
        configure_file("${test_config}" "test_${test_name}/configuration.qmom" COPYONLY)
        add_test(
            NAME cqmom_2d_variance_reinjection_${test_name}
            COMMAND cqmom_2d_variance_reinjection_example
            WORKING_DIRECTORY "${test_directory}"
        )
    endwhile()
endif()
