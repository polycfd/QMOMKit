add_executable(configuration_failures_verification configuration_failures.cpp)
target_link_libraries(configuration_failures_verification PRIVATE qmomkit)

set(configuration_failure_cases
    malformed_syntax.qmom
    qmom_invalid_adaptive.qmom
    qmom_missing_moments.qmom
    qmom_insufficient_moments.qmom
    cqmom_invalid_vandermonde_solver.qmom
    cqmom_invalid_adaptive_branch_counts.qmom
    cqmom_missing_phase_space.qmom
    cqmom_invalid_conditioning_order.qmom
    cqmom_incomplete_tensor.qmom
)

list(TRANSFORM configuration_failure_cases PREPEND "${CMAKE_CURRENT_LIST_DIR}/")

add_test(
    NAME verification_configuration_failure_paths
    COMMAND
        configuration_failures_verification
        ${configuration_failure_cases}
)
set_tests_properties(verification_configuration_failure_paths PROPERTIES LABELS verification)
