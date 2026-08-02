#pragma once

#include "conditioning_types.hpp"

#include <filesystem>
#include <vector>

namespace cqmom_2d_borderline_conditioning
{

void write_diagnostics_file(const std::vector<RunRecord>& records, const std::filesystem::path& filename);

void write_markdown_report(const std::vector<RunRecord>& records, const std::filesystem::path& filename);

void validate_report_file(const std::filesystem::path& filename);

void validate_successful_rows(const std::vector<RunRecord>& records);

void validate_report_inputs(const std::vector<RunRecord>& records);

void print_terminal_summary(const std::vector<RunRecord>& records);

} // namespace cqmom_2d_borderline_conditioning
