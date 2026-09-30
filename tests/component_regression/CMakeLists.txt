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

# The reference-value and dense/lean parity matrix tests share one case matrix,
# so their CTest case names are listed once here.
set(cqmom_matrix_cases
    static_2d_fixed_identity
    static_3d_fixed_identity
    static_5d_fixed_identity
    static_2d_fixed_reverse_order
    static_3d_fixed_reverse_order
    static_2d_fixed_unequal_counts
    static_3d_adaptive_collapsed_branch
    static_5d_adaptive_collapsed_branch
    static_3d_bounded_support
    static_3d_solver_automatic
    static_3d_solver_dense
    static_3d_solver_pivoted_qr
    static_3d_solver_svd
    static_3d_solver_moore_penrose
    static_3d_solver_bjorck_pereyra
    evolution_3d_strang
    evolution_3d_advsrc
    evolution_3d_srcadv
    evolution_3d_adaptive_strang
    evolution_5d_strang
)

qmomkit_add_test_executable(component_regression_cqmom_reference_values cqmom_reference_values.cpp cqmom_case_matrix.cpp)
foreach(case_name IN LISTS cqmom_matrix_cases)
    qmomkit_register_test_case(component_regression_cqmom_reference_values "component_regression_cqmom_reference_${case_name}" component_regression "${case_name}")
endforeach()

qmomkit_add_test_executable(component_regression_cqmom_dense_lean_parity_matrix cqmom_dense_lean_parity_matrix.cpp cqmom_case_matrix.cpp)
foreach(case_name IN LISTS cqmom_matrix_cases)
    qmomkit_register_test_case(component_regression_cqmom_dense_lean_parity_matrix "component_regression_cqmom_dense_lean_matrix_${case_name}" component_regression "${case_name}")
endforeach()
