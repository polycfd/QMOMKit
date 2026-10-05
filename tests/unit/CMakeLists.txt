set(unit_test_sources
    abscissa_bounds.cpp
    wheeler.cpp
    runge_kutta.cpp
    recursive_quadrature.cpp
    bjorck_pereyra.cpp
    householder_qr.cpp
    jacobi_svd.cpp
    linear_algebra.cpp
    vandermonde.cpp
    public_api.cpp
    configuration.cpp
    output.cpp
    distributions.cpp
    realizability.cpp
    realizability_policy.cpp
    node_advection.cpp
    moment_source.cpp
    moment_reconstruction.cpp
    phase_space_evolution.cpp
    phase_space_evolution_batch.cpp
    performance_instrumentation.cpp
)

foreach(source_file IN LISTS unit_test_sources)
    get_filename_component(test_name "${source_file}" NAME_WE)
    set(target_name "test_${test_name}")
    qmomkit_add_test_executable(${target_name} "${source_file}")
    qmomkit_register_source_test_cases(${target_name} unit "${CMAKE_CURRENT_LIST_DIR}/${source_file}")
endforeach()
