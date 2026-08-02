#pragma once

#include "example_workflow.hpp"
#include "unbounded_diagnostic_types.hpp"

#include <filesystem>
#include <string>

namespace cqmom_non_spherical_particles
{

void write_failure_context(const std::filesystem::path& directory, const qmomkit::CqmomQuadrature& last_quadrature,
                           const PreparedCqmomExample& example, double time, const std::string& message);

void write_report_file(const std::filesystem::path& directory, const DiagnosticReport& report, const PreparedCqmomExample& unbounded,
                       const PreparedCqmomExample& bounded);

} // namespace cqmom_non_spherical_particles
