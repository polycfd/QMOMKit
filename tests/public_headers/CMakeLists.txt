add_executable(
    test_public_headers
    main.cpp
    include_conditional_moments.cpp
    include_conditioning_path.cpp
    include_configuration.cpp
    include_cqmom.cpp
    include_distributions.cpp
    include_moment_set.cpp
    include_moment_rates.cpp
    include_moment_source.cpp
    include_qmomkit.cpp
    include_multi_index.cpp
    include_node_advection.cpp
    include_output.cpp
    include_phase_space_hierarchy.cpp
    include_qmom.cpp
    include_quadrature.cpp
    include_moment_reconstruction.cpp
    include_realizability.cpp
    include_realizability_policy.cpp
    include_runge_kutta.cpp
    include_transport_integrator.cpp
    include_vandermonde_solver.cpp
    include_wheeler.cpp
)
target_link_libraries(test_public_headers PRIVATE qmomkit::qmomkit)
add_test(NAME public_headers COMMAND test_public_headers)
set_tests_properties(public_headers PROPERTIES LABELS public_headers)
