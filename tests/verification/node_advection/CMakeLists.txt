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
