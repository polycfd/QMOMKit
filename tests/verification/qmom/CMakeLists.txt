qmomkit_add_ctest(test_wheeler_analytical verification_qmom_wheeler_analytical verification wheeler_analytical.cpp)
qmomkit_add_ctest(test_qmom_compute_failures verification_qmom_compute_failures verification qmom_compute_failures.cpp)

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
    set_tests_properties(verification_qmom_analytical_${case_name} PROPERTIES LABELS verification)
endforeach()
