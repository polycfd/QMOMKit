qmomkit_add_ctest(test_recursive_quadrature_manufactured verification_cqmom_recursive_quadrature_manufactured verification recursive_quadrature.cpp)
qmomkit_add_ctest(test_cqmom_quadrature_bounds verification_cqmom_quadrature_bounds verification quadrature_bounds.cpp)
qmomkit_add_ctest(test_recursive_validation verification_cqmom_recursive_validation verification recursive_cqmom.cpp)
qmomkit_add_ctest(test_conditional_moments_analytical verification_cqmom_conditional_moments_analytical verification conditional_moments.cpp)
qmomkit_add_ctest(test_cqmom_dense_lean_parity_3d verification_cqmom_dense_lean_parity_3d verification dense_lean_parity_3d.cpp)
qmomkit_add_ctest(test_cqmom_physical_transport_probe verification_cqmom_physical_transport_probe verification physical_transport_probe.cpp)

add_executable(moment_reconstruction_verification moment_reconstruction.cpp)
target_link_libraries(moment_reconstruction_verification PRIVATE qmomkit)
add_test(
    NAME verification_cqmom_moment_reconstruction_selected
    COMMAND
        moment_reconstruction_verification
        ${CMAKE_CURRENT_LIST_DIR}/independent_product.qmom
)
set_tests_properties(verification_cqmom_moment_reconstruction_selected PROPERTIES LABELS verification)

add_executable(cqmom_public_workflows_verification cqmom_public_workflows.cpp)
target_link_libraries(cqmom_public_workflows_verification PRIVATE qmomkit)
add_test(
    NAME verification_cqmom_public_workflows
    COMMAND
        cqmom_public_workflows_verification
        ${CMAKE_CURRENT_LIST_DIR}/independent_product.qmom
)
set_tests_properties(verification_cqmom_public_workflows PROPERTIES LABELS verification)

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
    set_tests_properties(verification_cqmom_${case_name} PROPERTIES LABELS verification)
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
    set_tests_properties(verification_cqmom_${case_name} PROPERTIES LABELS verification)
endforeach()
