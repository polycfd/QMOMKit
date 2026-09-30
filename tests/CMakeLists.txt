cmake_minimum_required(VERSION 3.16)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message(FATAL_ERROR "Configure QMOMKit tests from the repository root with -DQMOMKIT_BUILD_TESTS=ON.")
endif()

set(QMOMKIT_REPOSITORY_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")

# Registered tests are fast; the slowest takes a few seconds in a Debug build.
# This bound is generous against that and exists so that a test whose failure
# mode is non-termination reports a failure instead of stalling for CTest's much
# larger default timeout.
set(QMOMKIT_TEST_TIMEOUT_SECONDS 300)

# Hold the tests to the same warning level as the library. This is set for the
# whole directory rather than per target because several test executables are
# declared with add_executable directly instead of through the helper below.
#
# Debug only: the tests assert through C/C++ assert(), which NDEBUG removes, so
# a Release build leaves the parameters and locals that feed those assertions
# unused and reports hundreds of warnings that say nothing about the code. The
# test suite is run in Debug, which is also where these warnings are actionable.
#
# Aggregate initialization in the tests deliberately leaves trailing fields at
# their default member initializers, which is all that
# -Wmissing-field-initializers reports here, so only that warning is disabled.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options($<$<CONFIG:Debug>:-Wall> $<$<CONFIG:Debug>:-Wextra> $<$<CONFIG:Debug>:-Wno-missing-field-initializers>)
endif()

function(qmomkit_add_test_executable target_name)
    add_executable(${target_name} ${ARGN})
    target_include_directories(${target_name} PRIVATE "${QMOMKIT_REPOSITORY_ROOT}/src" "${QMOMKIT_REPOSITORY_ROOT}/tests")
    target_link_libraries(${target_name} PRIVATE qmomkit)
endfunction()

function(qmomkit_register_test_case target_name test_name label selector)
    add_test(NAME ${test_name} COMMAND ${target_name} ${selector})
    set_tests_properties(${test_name} PROPERTIES LABELS ${label} TIMEOUT ${QMOMKIT_TEST_TIMEOUT_SECONDS})
endfunction()

function(qmomkit_register_source_test_cases target_name label source_file)
    file(STRINGS "${source_file}" test_case_definitions REGEX "^[ \t]*(static[ \t]+)?void[ \t]+test_[A-Za-z0-9_]+[ \t]*\\(")
    string(REGEX REPLACE "^(test_|component_regression_)" "" test_subsystem "${target_name}")
    foreach(definition IN LISTS test_case_definitions)
        string(REGEX REPLACE "^[ \t]*(static[ \t]+)?void[ \t]+test_([A-Za-z0-9_]+)[ \t]*\\(.*" "\\2" selector "${definition}")
        qmomkit_register_test_case(${target_name} "${label}_${test_subsystem}_${selector}" ${label} ${selector})
    endforeach()
endfunction()

function(qmomkit_register_source_logical_cases target_name label source_file)
    file(STRINGS "${source_file}" test_case_definitions REGEX "^[ \t]*(static[ \t]+)?void[ \t]+(test|verify)_[A-Za-z0-9_]+[ \t]*\\([ \t]*\\)")
    string(REGEX REPLACE "^(test_|component_regression_)" "" test_subsystem "${target_name}")
    foreach(definition IN LISTS test_case_definitions)
        string(REGEX REPLACE "^[ \t]*(static[ \t]+)?void[ \t]+(test|verify)_([A-Za-z0-9_]+)[ \t]*\\([ \t]*\\).*" "\\3" selector "${definition}")
        qmomkit_register_test_case(${target_name} "${label}_${test_subsystem}_${selector}" ${label} ${selector})
    endforeach()
endfunction()

add_subdirectory(unit)
add_subdirectory(component_regression)
add_subdirectory(verification)
add_subdirectory(public_headers)
