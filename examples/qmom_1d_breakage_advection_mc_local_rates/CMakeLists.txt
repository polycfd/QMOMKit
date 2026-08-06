cmake_minimum_required(VERSION 3.16)

project(QMOMKitQmom1DAdvectionBreakageMcLocalRatesExample LANGUAGES CXX)

find_package(QMOMKitMonteCarlo CONFIG REQUIRED)

add_executable(qmom_1d_breakage_advection_mc_local_rates_example main.cpp local_rate_workflow.cpp physics.cpp sampling.cpp utilities.cpp)
target_compile_features(qmom_1d_breakage_advection_mc_local_rates_example PRIVATE cxx_std_17)
target_link_libraries(
    qmom_1d_breakage_advection_mc_local_rates_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(qmom_1d_breakage_advection_mc_local_rates_mc_example main_mc.cpp particles_mc.cpp physics.cpp sampling.cpp utilities.cpp)
target_compile_features(qmom_1d_breakage_advection_mc_local_rates_mc_example PRIVATE cxx_std_17)
target_link_libraries(
    qmom_1d_breakage_advection_mc_local_rates_mc_example
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

add_executable(qmom_1d_breakage_advection_mc_local_rates_nucleation_balance_check nucleation_balance_check.cpp particles_mc.cpp physics.cpp)
target_compile_features(qmom_1d_breakage_advection_mc_local_rates_nucleation_balance_check PRIVATE cxx_std_17)
target_link_libraries(
    qmom_1d_breakage_advection_mc_local_rates_nucleation_balance_check
    PRIVATE
        qmomkit_monte_carlo::qmomkit_monte_carlo
)

configure_file(configuration.qmom configuration.qmom COPYONLY)
