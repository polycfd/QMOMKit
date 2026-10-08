qmomkit_add_test_executable(configuration_failures_verification configuration_failures.cpp)

set(configuration_failure_cases
    malformed_syntax.qmom
    qmom_invalid_adaptive.qmom
    qmom_missing_moments.qmom
    qmom_insufficient_moments.qmom
    cqmom_invalid_vandermonde_solver.qmom
    cqmom_removed_adaptive_branch_counts.qmom
    cqmom_missing_phase_space.qmom
    cqmom_invalid_conditioning_order.qmom
    cqmom_incomplete_tensor.qmom
)

foreach(case_name IN LISTS configuration_failure_cases)
    get_filename_component(case_stem "${case_name}" NAME_WE)
    qmomkit_register_test_case(configuration_failures_verification "verification_configuration_${case_stem}" verification
        "${CMAKE_CURRENT_LIST_DIR}/${case_name}")
endforeach()
