add_executable(moment_source_analytical_verification moment_source_analytical.cpp)
target_link_libraries(moment_source_analytical_verification PRIVATE qmomkit)
set(moment_source_cases
    order_dependent_relaxation_qmom order_dependent_relaxation
    common_rate_relaxation_cqmom common_rate_relaxation
    linear_coupled_source_qmom linear_coupled
)
list(LENGTH moment_source_cases moment_source_case_count)
math(EXPR moment_source_last_index "${moment_source_case_count} - 1")
foreach(index RANGE 0 ${moment_source_last_index} 2)
    math(EXPR case_index "${index} + 1")
    list(GET moment_source_cases ${index} config_stem)
    list(GET moment_source_cases ${case_index} case_name)
    add_test(
        NAME verification_moment_source_${case_name}
        COMMAND
            moment_source_analytical_verification
            ${CMAKE_CURRENT_LIST_DIR}/${config_stem}.qmom
            ${case_name}
    )
    set_tests_properties(verification_moment_source_${case_name} PROPERTIES LABELS verification TIMEOUT ${QMOMKIT_TEST_TIMEOUT_SECONDS})
endforeach()

qmomkit_add_test_executable(test_moment_source_cqmom_update cqmom_source_update.cpp)
qmomkit_register_source_logical_cases(test_moment_source_cqmom_update verification "${CMAKE_CURRENT_LIST_DIR}/cqmom_source_update.cpp")
