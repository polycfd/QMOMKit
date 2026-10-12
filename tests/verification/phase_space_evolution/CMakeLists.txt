add_executable(phase_space_evolution_analytical_verification phase_space_evolution_analytical.cpp)
target_link_libraries(phase_space_evolution_analytical_verification PRIVATE qmomkit)
add_executable(qmom_linear_advection_manufactured_verification qmom_linear_advection_manufactured.cpp)
target_link_libraries(qmom_linear_advection_manufactured_verification PRIVATE qmomkit)
foreach(case_name IN ITEMS direct_advection phase_space_source node_source)
    if(case_name STREQUAL direct_advection)
        set(test_name verification_qmom_advection_matches_exponential_phase_evolution)
    elseif(case_name STREQUAL node_source)
        set(test_name verification_qmom_node_source_matches_exact_translation)
    else()
        set(test_name verification_qmom_phase_space_source_evolution_matches_manufactured_solution)
    endif()
    add_test(
        NAME ${test_name}
        COMMAND qmom_linear_advection_manufactured_verification ${case_name}
    )
    set_tests_properties(${test_name} PROPERTIES LABELS verification TIMEOUT ${QMOMKIT_TEST_TIMEOUT_SECONDS})
endforeach()
foreach(case_name IN ITEMS commuting_qmom qmom_boundary_contribution qmom_phase_space_evolution_failures qmom_phase_space_evolution_inversion_failures noncommuting_qmom
                           noncommuting_cqmom commuting_cqmom_2d cqmom_boundary_contribution commuting_cqmom_state_adaptive_2d cqmom_native_source
                           cqmom_combined_native_source cqmom_native_source_failures
                           cqmom_phase_space_evolution_failures cqmom_phase_space_evolution_inversion_failures cqmom_case_phase_space_evolution_failures)
    if(case_name STREQUAL noncommuting_qmom)
        set(test_name verification_phase_space_evolution_noncommuting_convergence)
    elseif(case_name STREQUAL noncommuting_cqmom)
        set(test_name verification_phase_space_evolution_noncommuting_cqmom_convergence)
    else()
        set(test_name verification_phase_space_evolution_${case_name})
    endif()
    set(configuration_case_name ${case_name})
    if(case_name STREQUAL qmom_boundary_contribution)
        set(configuration_case_name commuting_qmom)
    elseif(case_name STREQUAL noncommuting_cqmom)
        set(configuration_case_name noncommuting_cqmom_2d)
    elseif(case_name STREQUAL cqmom_boundary_contribution)
        set(configuration_case_name commuting_cqmom_2d)
    elseif(case_name STREQUAL cqmom_case_phase_space_evolution_failures)
        set(configuration_case_name cqmom_phase_space_evolution_failures)
    endif()
    add_test(
        NAME ${test_name}
        COMMAND
            phase_space_evolution_analytical_verification
            ${CMAKE_CURRENT_LIST_DIR}/${configuration_case_name}.qmom
            ${case_name}
    )
    set_tests_properties(${test_name} PROPERTIES LABELS verification TIMEOUT ${QMOMKIT_TEST_TIMEOUT_SECONDS})
endforeach()

add_test(
    NAME verification_phase_space_evolution_cqmom_state_is_unit_invariant
    COMMAND phase_space_evolution_analytical_verification cqmom_state_is_unit_invariant
)
set_tests_properties(verification_phase_space_evolution_cqmom_state_is_unit_invariant PROPERTIES LABELS verification TIMEOUT ${QMOMKIT_TEST_TIMEOUT_SECONDS})

qmomkit_add_test_executable(test_phase_space_evolution_adaptive_qmom adaptive_qmom_evolution.cpp)
qmomkit_register_source_logical_cases(test_phase_space_evolution_adaptive_qmom verification "${CMAKE_CURRENT_LIST_DIR}/adaptive_qmom_evolution.cpp")
qmomkit_add_test_executable(test_phase_space_evolution_source_retries source_step_retries.cpp)
qmomkit_register_source_logical_cases(test_phase_space_evolution_source_retries verification "${CMAKE_CURRENT_LIST_DIR}/source_step_retries.cpp")
qmomkit_add_test_executable(test_phase_space_evolution_time_dependent_splitting time_dependent_splitting.cpp)
qmomkit_register_source_logical_cases(test_phase_space_evolution_time_dependent_splitting verification
                                      "${CMAKE_CURRENT_LIST_DIR}/time_dependent_splitting.cpp")
