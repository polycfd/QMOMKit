cmake_minimum_required(VERSION 3.16)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message(FATAL_ERROR "Configure QMOMKit tests from the repository root with -DQMOMKIT_BUILD_TESTS=ON.")
endif()

set(QMOMKIT_REPOSITORY_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")

function(qmomkit_add_ctest target_name test_name label)
    add_executable(${target_name} ${ARGN})
    target_include_directories(${target_name} PRIVATE "${QMOMKIT_REPOSITORY_ROOT}/src")
    target_link_libraries(${target_name} PRIVATE qmomkit)
    add_test(NAME ${test_name} COMMAND ${target_name})
    set_tests_properties(${test_name} PROPERTIES LABELS ${label})
endfunction()

add_subdirectory(unit)
add_subdirectory(component_regression)
add_subdirectory(verification)
add_subdirectory(public_headers)
