#pragma once

#include "qmomkit/configuration.hpp"
#include "qmomkit/cqmom.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace cqmom_2d_minimal
{

void print_configuration_errors(const std::string& error);

void ensure_output_directory(const std::filesystem::path& path);

void print_result_summary(const qmomkit::SelectedRawMomentResult& initial_moments, const qmomkit::SelectedRawMomentResult& final_moments,
                          const qmomkit::CqmomQuadrature& final_quadrature, const std::filesystem::path& output_directory, bool wrote_moments,
                          bool wrote_quadrature);

} // namespace cqmom_2d_minimal
