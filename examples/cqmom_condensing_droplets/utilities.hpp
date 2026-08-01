#pragma once

#include "monte_carlo/ensemble.hpp"
#include "qmomkit/configuration.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/node_advection.hpp"
#include "qmomkit/transport_integrator.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace cqmom_condensing_droplets
{

void print_configuration_errors(const std::string& error);

void ensure_output_directory(const std::filesystem::path& path);

void print_cqmom_run_configuration(const qmomkit::CqmomCase& cqmom_case, const qmomkit::CqmomTransportOptions& transport_options);

void print_monte_carlo_run_configuration(const qmomkit::CqmomCase& cqmom_case, const qmomkit::NodeAdvectionOptions& node_options);

[[nodiscard]] qmomkit::SelectedRawMomentResult compute_ensemble_moments(const qmomkit::monte_carlo::Ensemble& ensemble, double particle_weight,
                                                                        const qmomkit::RawMomentRequest& request);

void print_result_summary(const qmomkit::SelectedRawMomentResult& initial_moments, const qmomkit::SelectedRawMomentResult& final_moments,
                          const qmomkit::CqmomQuadrature& final_quadrature, const std::filesystem::path& output_directory, std::size_t sample_count,
                          unsigned long long initialization_seed, double time_step, double final_time, bool wrote_moments, bool wrote_quadrature);

void print_monte_carlo_result_summary(const qmomkit::SelectedRawMomentResult& initial_moments, const qmomkit::SelectedRawMomentResult& final_moments,
                                      const std::filesystem::path& output_directory, std::size_t final_particle_count, std::size_t sample_count,
                                      unsigned long long initialization_seed, double time_step, double final_time, bool wrote_moments);

} // namespace cqmom_condensing_droplets
