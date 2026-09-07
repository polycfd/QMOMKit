qmomkit_add_test_executable(output_workflows_verification output_workflows.cpp)
foreach(case_name IN ITEMS qmom cqmom)
    qmomkit_register_test_case(output_workflows_verification "verification_output_writer_${case_name}" verification
        "${CMAKE_CURRENT_LIST_DIR}/../qmom/two_node_discrete.qmom;${CMAKE_CURRENT_LIST_DIR}/../cqmom/independent_product.qmom;${CMAKE_CURRENT_BINARY_DIR}/output_writer_${case_name};${case_name}")
endforeach()
