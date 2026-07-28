add_executable(output_workflows_verification output_workflows.cpp)
target_link_libraries(output_workflows_verification PRIVATE qmomkit)

add_test(
    NAME verification_output_writer_workflows
    COMMAND
        output_workflows_verification
        ${CMAKE_CURRENT_LIST_DIR}/../qmom/two_node_discrete.qmom
        ${CMAKE_CURRENT_LIST_DIR}/../cqmom/independent_product.qmom
        ${CMAKE_CURRENT_BINARY_DIR}/output_writer_workflows
)
set_tests_properties(verification_output_writer_workflows PROPERTIES LABELS verification)
