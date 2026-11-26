cmake_minimum_required(VERSION 3.16)

project(QMOMKitCqmom2DHyperbolicDragDispersionWithMcExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

# The QMOMKit driver is built twice. Each executable reads the configuration
# file named by its compile definition, and the two files differ only in the
# conditioning order.
add_executable(cqmom_2d_hyperbolic_drag_dispersion_with_mc_velocity_first_example main.cpp physics.cpp sampling.cpp utilities.cpp)
target_compile_definitions(cqmom_2d_hyperbolic_drag_dispersion_with_mc_velocity_first_example PRIVATE CONFIGURATION_FILE="velocity_first.qmom")
target_link_libraries(
    cqmom_2d_hyperbolic_drag_dispersion_with_mc_velocity_first_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(cqmom_2d_hyperbolic_drag_dispersion_with_mc_size_first_example main.cpp physics.cpp sampling.cpp utilities.cpp)
target_compile_definitions(cqmom_2d_hyperbolic_drag_dispersion_with_mc_size_first_example PRIVATE CONFIGURATION_FILE="size_first.qmom")
target_link_libraries(
    cqmom_2d_hyperbolic_drag_dispersion_with_mc_size_first_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(cqmom_2d_hyperbolic_drag_dispersion_with_mc_mc_example main_mc.cpp particles_mc.cpp physics.cpp sampling.cpp utilities.cpp)
target_link_libraries(
    cqmom_2d_hyperbolic_drag_dispersion_with_mc_mc_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

configure_file(velocity_first.qmom velocity_first.qmom COPYONLY)
configure_file(size_first.qmom size_first.qmom COPYONLY)
