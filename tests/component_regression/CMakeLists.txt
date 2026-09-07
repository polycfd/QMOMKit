set(component_case_sources
    moment_source.cpp
    node_advection.cpp
    phase_space_evolution.cpp
    cqmom_case_moment_replacement.cpp
)

foreach(source_file IN LISTS component_case_sources)
    get_filename_component(test_name "${source_file}" NAME_WE)
    set(target_name "component_regression_${test_name}")
    qmomkit_add_test_executable(${target_name} "${source_file}")
    qmomkit_register_source_test_cases(${target_name} component_regression "${CMAKE_CURRENT_LIST_DIR}/${source_file}")
endforeach()

qmomkit_add_test_executable(component_regression_cqmom_compute_failures cqmom_compute_failures.cpp)
foreach(case_name IN ITEMS nonrealizable_public_workflow bounded_abscissa_violation dense_extent_rank dense_value_count dense_nonfinite_moment dense_missing_required_moment lean_missing_required_value lean_nonfinite_value lean_incompatible_layout adaptive_policy)
    qmomkit_register_test_case(component_regression_cqmom_compute_failures "component_regression_cqmom_compute_${case_name}" component_regression "${case_name}")
endforeach()

qmomkit_add_test_executable(component_regression_cqmom_dense_lean_parity_3d cqmom_dense_lean_parity_3d.cpp)
foreach(case_name IN ITEMS fixed_identity fixed_reverse_order adaptive_branch_specific near_degenerate_conditional_branch physical_analytical_contraction failure_status phase_space_evolution)
    qmomkit_register_test_case(component_regression_cqmom_dense_lean_parity_3d "component_regression_cqmom_dense_lean_${case_name}" component_regression "${case_name}")
endforeach()

qmomkit_add_test_executable(component_regression_cqmom_physical_transport_probe cqmom_physical_transport_probe.cpp)
qmomkit_register_test_case(component_regression_cqmom_physical_transport_probe component_regression_cqmom_physical_transport component_regression physical_transport)
