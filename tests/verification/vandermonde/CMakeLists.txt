qmomkit_add_test_executable(vandermonde_solver_selection_verification solver_selection.cpp)
foreach(case_name IN ITEMS automatic dense pivoted_qr svd moore_penrose bjorck_pereyra)
    qmomkit_register_test_case(vandermonde_solver_selection_verification "verification_vandermonde_solver_${case_name}" verification
        "${CMAKE_CURRENT_LIST_DIR};${case_name}")
endforeach()
qmomkit_register_test_case(vandermonde_solver_selection_verification verification_vandermonde_solver_dense_matches_bjorck_pereyra verification
    "${CMAKE_CURRENT_LIST_DIR};dense_matches_bjorck_pereyra")

qmomkit_add_test_executable(clustered_solver_policy_verification clustered_solver_policy.cpp)
foreach(case_name IN ITEMS automatic dense pivoted_qr svd moore_penrose bjorck_pereyra)
    qmomkit_register_test_case(clustered_solver_policy_verification "verification_vandermonde_clustered_${case_name}" verification
        "${CMAKE_CURRENT_LIST_DIR}/clustered_policy_automatic.qmom;${case_name}")
endforeach()

qmomkit_add_test_executable(test_vandermonde_diagnostics diagnostics_invariance.cpp)
qmomkit_register_source_logical_cases(test_vandermonde_diagnostics verification "${CMAKE_CURRENT_LIST_DIR}/diagnostics_invariance.cpp")
