add_executable(node_advection_analytical_verification node_advection_analytical.cpp)
target_link_libraries(node_advection_analytical_verification PRIVATE qmomkit)
foreach(case_name IN ITEMS constant_velocity linear_growth affine_evolution linear_system_2d)
    if(case_name STREQUAL affine_evolution)
        set(test_name verification_node_advection_affine)
    else()
        set(test_name verification_node_advection_${case_name})
    endif()
    add_test(
        NAME ${test_name}
        COMMAND
            node_advection_analytical_verification
            ${CMAKE_CURRENT_LIST_DIR}/${case_name}.qmom
            ${case_name}
    )
    set_tests_properties(${test_name} PROPERTIES LABELS verification TIMEOUT ${QMOMKIT_TEST_TIMEOUT_SECONDS})
endforeach()

qmomkit_add_test_executable(test_node_advection_ssp_rk3_stage_times ssp_rk3_stage_times.cpp)
qmomkit_register_source_logical_cases(test_node_advection_ssp_rk3_stage_times verification "${CMAKE_CURRENT_LIST_DIR}/ssp_rk3_stage_times.cpp")
