cmake_minimum_required(VERSION 3.16)

project(QMOMKitPaperVerification LANGUAGES CXX)

find_package(QMOMKit CONFIG REQUIRED)

add_executable(qmomkit_paper_quadrature_1d quadrature_1d/main.cpp)
add_executable(qmomkit_paper_cqmom_affine_chain cqmom_affine_chain/main.cpp)
add_executable(qmomkit_paper_adaptive_branch_quadrature adaptive_branch_quadrature/main.cpp)
add_executable(qmomkit_paper_time_convergence time_convergence/main.cpp)

foreach(target_name IN ITEMS qmomkit_paper_quadrature_1d qmomkit_paper_cqmom_affine_chain qmomkit_paper_adaptive_branch_quadrature
                             qmomkit_paper_time_convergence)
    target_compile_features(${target_name} PRIVATE cxx_std_17)
    target_link_libraries(${target_name} PRIVATE qmomkit::qmomkit)

    if(MSVC)
        target_compile_options(${target_name} PRIVATE /W4 /permissive-)
    else()
        target_compile_options(${target_name} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endforeach()

configure_file(adaptive_branch_quadrature/configuration.qmom configuration.qmom COPYONLY)
