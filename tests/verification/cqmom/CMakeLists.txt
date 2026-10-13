qmomkit_add_test_executable(test_recursive_quadrature_manufactured recursive_quadrature.cpp)
qmomkit_register_source_logical_cases(test_recursive_quadrature_manufactured verification "${CMAKE_CURRENT_LIST_DIR}/recursive_quadrature.cpp")
qmomkit_add_test_executable(test_cqmom_quadrature_bounds quadrature_bounds.cpp)
qmomkit_register_source_logical_cases(test_cqmom_quadrature_bounds verification "${CMAKE_CURRENT_LIST_DIR}/quadrature_bounds.cpp")
qmomkit_add_test_executable(test_cqmom_adaptive_branch_growth adaptive_branch_growth.cpp)
qmomkit_register_source_logical_cases(test_cqmom_adaptive_branch_growth verification "${CMAKE_CURRENT_LIST_DIR}/adaptive_branch_growth.cpp")
qmomkit_add_test_executable(test_cqmom_power_of_two_scaling power_of_two_scaling.cpp)
qmomkit_register_source_logical_cases(test_cqmom_power_of_two_scaling verification "${CMAKE_CURRENT_LIST_DIR}/power_of_two_scaling.cpp")
qmomkit_add_test_executable(test_recursive_validation recursive_cqmom.cpp)
foreach(case_name IN ITEMS 1d 2d_identity 2d_reversed 3d_non_identity non_uniform_fixed_counts 4d_adversarial_non_identity 4d_adversarial_adaptive_non_identity dimension_sweep_D2 dimension_sweep_D3 dimension_sweep_D4 dimension_sweep_D5 dimension_sweep_D6)
    qmomkit_register_test_case(test_recursive_validation "verification_cqmom_recursive_${case_name}" verification "${case_name}")
endforeach()

qmomkit_add_test_executable(cqmom_unit_invariance_verification unit_invariance.cpp)
qmomkit_register_test_case(cqmom_unit_invariance_verification verification_cqmom_unit_invariance verification unit-invariance)
qmomkit_register_test_case(cqmom_unit_invariance_verification verification_cqmom_disparate_units verification disparate-units)

qmomkit_add_test_executable(moment_reconstruction_verification moment_reconstruction.cpp)
foreach(case_name IN ITEMS qmom_all_orders cqmom_selected cqmom_required_layout failure_paths)
    qmomkit_register_test_case(moment_reconstruction_verification "verification_cqmom_moment_reconstruction_${case_name}" verification
        "${CMAKE_CURRENT_LIST_DIR}/independent_product.qmom;${case_name}")
endforeach()

add_executable(cqmom_public_workflows_verification cqmom_public_workflows.cpp)
target_link_libraries(cqmom_public_workflows_verification PRIVATE qmomkit)
add_test(
    NAME verification_cqmom_public_workflows
    COMMAND
        cqmom_public_workflows_verification
        ${CMAKE_CURRENT_LIST_DIR}/independent_product.qmom
)
set_tests_properties(verification_cqmom_public_workflows PROPERTIES LABELS verification TIMEOUT ${QMOMKIT_TEST_TIMEOUT_SECONDS})

add_executable(cqmom_2d_analytical_verification cqmom_2d_analytical.cpp)
target_link_libraries(cqmom_2d_analytical_verification PRIVATE qmomkit)
foreach(case_name IN ITEMS independent_product deterministic_affine branch_dependent)
    add_test(
        NAME verification_cqmom_${case_name}
        COMMAND
            cqmom_2d_analytical_verification
            ${CMAKE_CURRENT_LIST_DIR}/${case_name}.qmom
            ${case_name}
    )
    set_tests_properties(verification_cqmom_${case_name} PROPERTIES LABELS verification TIMEOUT ${QMOMKIT_TEST_TIMEOUT_SECONDS})
endforeach()

add_executable(cqmom_high_dimensional_analytical cqmom_high_dimensional_analytical.cpp)
target_link_libraries(cqmom_high_dimensional_analytical PRIVATE qmomkit)
foreach(case_name IN ITEMS affine_chain_5d independent_product_3d independent_product_4d independent_product_5d
                           independent_product_6d mixed_affine_independent_4d)
    add_test(
        NAME verification_cqmom_${case_name}
        COMMAND
            cqmom_high_dimensional_analytical
            ${CMAKE_CURRENT_LIST_DIR}/${case_name}.qmom
            ${case_name}
    )
    set_tests_properties(verification_cqmom_${case_name} PROPERTIES LABELS verification TIMEOUT ${QMOMKIT_TEST_TIMEOUT_SECONDS})
endforeach()
