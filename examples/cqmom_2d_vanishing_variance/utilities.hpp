#pragma once

#include "qmomkit/cqmom.hpp"
#include "qmomkit/moment_reconstruction.hpp"

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

namespace cqmom_2d_vanishing_variance
{

void print_configuration_errors(const std::string& error);

void ensure_output_directory(const std::filesystem::path& path);

void print_result_summary(const qmomkit::SelectedRawMomentResult& initial_moments, const qmomkit::SelectedRawMomentResult& final_moments,
                          const std::vector<std::pair<double, qmomkit::CqmomQuadrature>>& quadrature_history,
                          const qmomkit::CqmomQuadrature& final_quadrature, const std::filesystem::path& output_directory,
                          const std::vector<std::size_t>& conditioning_order, double time_step, double final_time, bool wrote_moments,
                          bool wrote_quadrature);

} // namespace cqmom_2d_vanishing_variance
