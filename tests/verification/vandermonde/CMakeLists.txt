add_executable(vandermonde_solver_selection_verification solver_selection.cpp)
target_include_directories(vandermonde_solver_selection_verification PRIVATE "${QMOMKIT_REPOSITORY_ROOT}/src")
target_link_libraries(vandermonde_solver_selection_verification PRIVATE qmomkit)

add_test(
    NAME verification_vandermonde_solver_selection
    COMMAND
        vandermonde_solver_selection_verification
        ${CMAKE_CURRENT_LIST_DIR}
)
set_tests_properties(verification_vandermonde_solver_selection PROPERTIES LABELS verification)
