add_executable(phase_space_evolution_analytical_verification phase_space_evolution_analytical.cpp)
target_link_libraries(phase_space_evolution_analytical_verification PRIVATE qmomkit)
add_executable(qmom_linear_advection_manufactured_verification qmom_linear_advection_manufactured.cpp)
target_link_libraries(qmom_linear_advection_manufactured_verification PRIVATE qmomkit)
foreach(case_name IN ITEMS direct_advection phase_space_source)
    if(case_name STREQUAL direct_advection)
        set(test_name verification_qmom_advection_matches_exponential_phase_evolution)
    else()
        set(test_name verification_qmom_phase_space_source_evolution_matches_manufactured_solution)
    endif()
    add_test(
        NAME ${test_name}
        COMMAND qmom_linear_advection_manufactured_verification ${case_name}
    )
    set_tests_properties(${test_name} PROPERTIES LABELS verification)
endforeach()
foreach(case_name IN ITEMS commuting_qmom qmom_boundary_contribution qmom_phase_space_evolution_failures qmom_phase_space_evolution_inversion_failures noncommuting_qmom
                           noncommuting_cqmom commuting_cqmom_2d cqmom_boundary_contribution commuting_cqmom_state_adaptive_2d cqmom_native_source
                           cqmom_combined_native_source cqmom_native_source_failures
                           cqmom_phase_space_evolution_failures cqmom_phase_space_evolution_inversion_failures cqmom_dense_phase_space_evolution_failures
                           cqmom_dense_phase_space_evolution_inversion_failures cqmom_lean_phase_space_evolution_failures cqmom_lean_phase_space_evolution_inversion_failures)
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
    elseif(case_name STREQUAL cqmom_dense_phase_space_evolution_failures OR case_name STREQUAL cqmom_lean_phase_space_evolution_failures)
        set(configuration_case_name cqmom_phase_space_evolution_failures)
    elseif(case_name STREQUAL cqmom_dense_phase_space_evolution_inversion_failures OR case_name STREQUAL cqmom_lean_phase_space_evolution_inversion_failures)
        set(configuration_case_name cqmom_phase_space_evolution_inversion_failures)
    endif()
    add_test(
        NAME ${test_name}
        COMMAND
            phase_space_evolution_analytical_verification
            ${CMAKE_CURRENT_LIST_DIR}/${configuration_case_name}.qmom
            ${case_name}
    )
    set_tests_properties(${test_name} PROPERTIES LABELS verification)
endforeach()

add_test(
    NAME verification_phase_space_evolution_cqmom_automatic_scaling_is_unit_invariant
    COMMAND phase_space_evolution_analytical_verification cqmom_automatic_scaling_is_unit_invariant
)
set_tests_properties(verification_phase_space_evolution_cqmom_automatic_scaling_is_unit_invariant PROPERTIES LABELS verification)

add_custom_target(
    generate_phase_space_evolution_verification_outputs
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/phase_space_evolution_verification_outputs
    COMMAND
        phase_space_evolution_analytical_verification ${CMAKE_CURRENT_LIST_DIR}/commuting_qmom.qmom commuting_qmom --output
        ${CMAKE_CURRENT_BINARY_DIR}/phase_space_evolution_verification_outputs/commuting_qmom_refinement.txt
    COMMAND
        phase_space_evolution_analytical_verification ${CMAKE_CURRENT_LIST_DIR}/noncommuting_qmom.qmom noncommuting_qmom --output
        ${CMAKE_CURRENT_BINARY_DIR}/phase_space_evolution_verification_outputs/noncommuting_qmom_refinement.txt
    DEPENDS phase_space_evolution_analytical_verification
)
