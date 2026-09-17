qmomkit_add_test_executable(test_wheeler_analytical wheeler_analytical.cpp)
qmomkit_register_source_logical_cases(test_wheeler_analytical verification "${CMAKE_CURRENT_LIST_DIR}/wheeler_analytical.cpp")
qmomkit_add_test_executable(test_adaptive_qmom adaptive_qmom.cpp)
qmomkit_register_source_logical_cases(test_adaptive_qmom verification "${CMAKE_CURRENT_LIST_DIR}/adaptive_qmom.cpp")
qmomkit_add_test_executable(test_qmom_moment_space_boundary moment_space_boundary.cpp)
qmomkit_register_source_logical_cases(test_qmom_moment_space_boundary verification "${CMAKE_CURRENT_LIST_DIR}/moment_space_boundary.cpp")
qmomkit_add_test_executable(test_qmom_quadrature_bounds quadrature_bounds.cpp)
qmomkit_register_source_logical_cases(test_qmom_quadrature_bounds verification "${CMAKE_CURRENT_LIST_DIR}/quadrature_bounds.cpp")
qmomkit_add_test_executable(test_qmom_compute_failures qmom_compute_failures.cpp)
foreach(case_name IN ITEMS negative_zeroth_moment negative_variance insufficient_moments zero_node_count)
    string(REPLACE "_" " " case_label "${case_name}")
    qmomkit_register_test_case(test_qmom_compute_failures "verification_qmom_compute_failure_${case_name}" verification "${case_label}")
endforeach()

add_executable(qmom_analytical_verification qmom_analytical.cpp)
target_link_libraries(qmom_analytical_verification PRIVATE qmomkit)
foreach(case_name IN ITEMS single_dirac two_node_discrete three_node_discrete uniform_n2 uniform_n3 exponential_n2 exponential_n3
                           gaussian_n2 gaussian_n3)
    add_test(
        NAME verification_qmom_analytical_${case_name}
        COMMAND
            qmom_analytical_verification
            ${CMAKE_CURRENT_LIST_DIR}/${case_name}.qmom
            ${case_name}
    )
    set_tests_properties(verification_qmom_analytical_${case_name} PROPERTIES LABELS verification TIMEOUT ${QMOMKIT_TEST_TIMEOUT_SECONDS})
endforeach()
