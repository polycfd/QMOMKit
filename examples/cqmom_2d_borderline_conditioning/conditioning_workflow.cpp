#include "conditioning_workflow.hpp"

#include "conditioning_report.hpp"
#include "conditioning_types.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

constexpr std::size_t dimensions = 2;
constexpr std::size_t configured_node_count = 2;
constexpr std::size_t maximum_order = 3;
constexpr double rank_tolerance = 1.0e-14;
constexpr double population_weight_tolerance = 32.0 * std::numeric_limits<double>::epsilon();

struct Atom
{
    double x = 0.0;
    double y = 0.0;
    double weight = 0.0;
};

using cqmom_2d_borderline_conditioning::delta_values;
using cqmom_2d_borderline_conditioning::MatrixDiagnostics;
using cqmom_2d_borderline_conditioning::RunRecord;
using cqmom_2d_borderline_conditioning::solver_choices;
using cqmom_2d_borderline_conditioning::SolverChoice;

const char* status_name(const qmomkit::CqmomStatus status)
{
    switch (status)
    {
    case qmomkit::CqmomStatus::success:
        return "success";
    case qmomkit::CqmomStatus::invalid_input:
        return "invalid_input";
    case qmomkit::CqmomStatus::qmom_failed:
        return "qmom_failed";
    case qmomkit::CqmomStatus::conditional_moment_failed:
        return "conditional_moment_failed";
    }
    return "unknown";
}

const char* solver_name(const qmomkit::VandermondeSolver solver)
{
    for (const SolverChoice& choice : solver_choices)
    {
        if (choice.solver == solver)
        {
            return choice.label;
        }
    }
    if (solver == qmomkit::VandermondeSolver::automatic)
    {
        return "automatic";
    }
    return "unknown";
}

std::string format_path(const qmomkit::ConditioningPath& path)
{
    if (path.empty())
    {
        return "root";
    }

    std::ostringstream stream;
    stream << '[';
    for (std::size_t entry = 0; entry < path.size(); ++entry)
    {
        if (entry > 0)
        {
            stream << ',';
        }
        stream << path.node_index_at_level(entry);
    }
    stream << ']';
    return stream.str();
}

std::string sanitize_message(std::string message)
{
    for (char& character : message)
    {
        if (std::isspace(static_cast<unsigned char>(character)))
        {
            character = '_';
        }
    }
    return message.empty() ? "-" : message;
}

std::vector<Atom> make_population(const double delta)
{
    return {{-delta / 2.0, -1.0, 0.10}, {-delta / 2.0, 1.0, 0.30}, {delta / 2.0, -1.0, 0.15}, {delta / 2.0, 1.0, 0.45}};
}

qmomkit::RawMomentTensor exact_moments(const std::vector<Atom>& population)
{
    qmomkit::RawMomentTensor moments{{maximum_order + 1, maximum_order + 1}, std::vector<double>((maximum_order + 1) * (maximum_order + 1), 0.0)};
    for (std::size_t p = 0; p <= maximum_order; ++p)
    {
        for (std::size_t q = 0; q <= maximum_order; ++q)
        {
            double value = 0.0;
            for (const Atom& atom : population)
            {
                value += atom.weight * std::pow(atom.x, static_cast<int>(p)) * std::pow(atom.y, static_cast<int>(q));
            }
            moments.values[p * (maximum_order + 1) + q] = value;
        }
    }
    return moments;
}

qmomkit::RawMomentTensor exact_moments(const double delta)
{
    return exact_moments(make_population(delta));
}

struct Matrix2x2
{
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    double d = 0.0;
};

Matrix2x2 diagnostic_matrix(const double delta, const std::vector<std::size_t>& conditioning_order)
{
    if (conditioning_order == std::vector<std::size_t>{0, 1})
    {
        return {0.4, 0.6, -0.2 * delta, 0.3 * delta};
    }
    return {0.25, 0.75, -0.25, 0.75};
}

MatrixDiagnostics analyse_matrix(const Matrix2x2& matrix)
{
    const double p = std::hypot(matrix.a + matrix.d, matrix.b - matrix.c);
    const double q = std::hypot(matrix.a - matrix.d, matrix.b + matrix.c);
    const double sigma_max = (p + q) / 2.0;
    const double sigma_min = std::abs(p - q) / 2.0;
    const double ratio = sigma_max > 0.0 ? sigma_min / sigma_max : 0.0;

    std::size_t rank = 0;
    if (sigma_max > 0.0)
    {
        rank += sigma_max > rank_tolerance * sigma_max ? 1 : 0;
        rank += sigma_min > rank_tolerance * sigma_max ? 1 : 0;
    }

    std::optional<double> condition_number;
    if (rank == 2 && sigma_min > 0.0)
    {
        condition_number = sigma_max / sigma_min;
    }

    return {sigma_max, sigma_min, ratio, rank, condition_number};
}

std::string order_value(const std::vector<std::size_t>& conditioning_order)
{
    return conditioning_order == std::vector<std::size_t>{0, 1} ? "(0,1)" : "(1,0)";
}

std::string active_counts(const qmomkit::CqmomResult& result)
{
    if (result.diagnostics.empty())
    {
        return "-";
    }

    std::ostringstream stream;
    bool first = true;
    for (const qmomkit::CqmomDiagnostics& diagnostic : result.diagnostics)
    {
        if (diagnostic.actual_node_count == 0)
        {
            continue;
        }
        if (!first)
        {
            stream << ';';
        }
        first = false;
        stream << 'L' << diagnostic.conditioning_level << '(' << format_path(diagnostic.parent_path) << ")=" << diagnostic.actual_node_count << '/'
               << diagnostic.requested_node_count;
    }
    return first ? "-" : stream.str();
}

std::optional<qmomkit::CqmomDiagnostics> conditional_diagnostic(const qmomkit::CqmomResult& result)
{
    for (const qmomkit::CqmomDiagnostics& diagnostic : result.diagnostics)
    {
        if (diagnostic.qmom_status == qmomkit::QmomStatus::success && diagnostic.conditioning_level + 1 < dimensions)
        {
            return diagnostic;
        }
    }
    return std::nullopt;
}

std::optional<double> moment_error(const qmomkit::CqmomResult& result, const qmomkit::RawMomentTensor& exact)
{
    if (!result.ok())
    {
        return std::nullopt;
    }

    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(result.quadrature, exact.extents);
    double maximum_error = 0.0;
    for (std::size_t entry = 0; entry < exact.values.size(); ++entry)
    {
        if (!std::isfinite(reconstructed.values[entry]))
        {
            throw std::runtime_error("Successful reconstruction produced a non-finite raw moment.");
        }
        maximum_error =
            std::max(maximum_error, std::abs(reconstructed.values[entry] - exact.values[entry]) / std::max(1.0, std::abs(exact.values[entry])));
    }
    return maximum_error;
}

void validate_successful_quadrature(const qmomkit::CqmomResult& result)
{
    if (!result.ok())
    {
        return;
    }
    for (const double weight : result.quadrature.weights())
    {
        if (!std::isfinite(weight))
        {
            throw std::runtime_error("Successful reconstruction produced a non-finite quadrature weight.");
        }
    }
    for (const std::vector<double>& dimension : result.quadrature.abscissas())
    {
        for (const double abscissa : dimension)
        {
            if (!std::isfinite(abscissa))
            {
                throw std::runtime_error("Successful reconstruction produced a non-finite quadrature abscissa.");
            }
        }
    }
}

RunRecord run_reconstruction(const double delta, const std::vector<std::size_t>& conditioning_order, const SolverChoice& choice,
                             const bool adaptive_branch_counts, const std::string& mode)
{
    const qmomkit::RawMomentTensor moments = exact_moments(delta);
    const qmomkit::PhaseSpaceHierarchy hierarchy({configured_node_count, configured_node_count}, conditioning_order);
    const qmomkit::CqmomDenseCase case_{hierarchy, moments, {}, {}, adaptive_branch_counts, qmomkit::CqmomDiagnosticLevel::detailed, choice.solver};
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(case_);
    validate_successful_quadrature(result);

    const std::optional<qmomkit::CqmomDiagnostics> diagnostic = conditional_diagnostic(result);
    RunRecord record;
    record.mode = mode;
    record.delta = delta;
    record.conditioning_order = order_value(conditioning_order);
    record.requested_solver = choice.label;
    record.selected_solver = diagnostic.has_value() ? solver_name(diagnostic->conditional_vandermonde_selected_solver) : "-";
    record.success = result.ok();
    record.status = status_name(result.status);
    record.fallback = diagnostic.has_value() && diagnostic->dense_vandermonde_fallback_occurred;
    record.active_counts = active_counts(result);
    record.diagnostic_level = diagnostic.has_value() ? std::to_string(diagnostic->conditioning_level) : "-";
    record.parent_path = diagnostic.has_value() ? format_path(diagnostic->parent_path) : "-";
    record.condition_estimate_inf = diagnostic.has_value() ? diagnostic->conditional_moment_diagnostics.condition_estimate : std::nullopt;
    record.matrix = analyse_matrix(diagnostic_matrix(delta, conditioning_order));
    record.residual_norm = diagnostic.has_value() ? diagnostic->conditional_moment_diagnostics.residual_norm : std::nullopt;
    record.moment_error = moment_error(result, moments);
    record.failed_dimension = result.ok() ? "-" : std::to_string(result.failed_phase_space_dimension);
    record.message = sanitize_message(result.message);
    if (result.ok() && ((!record.condition_estimate_inf.has_value() || std::isfinite(*record.condition_estimate_inf)) &&
                        (!record.residual_norm.has_value() || std::isfinite(*record.residual_norm)) &&
                        (!record.moment_error.has_value() || std::isfinite(*record.moment_error))))
    {
        return record;
    }
    if (result.ok())
    {
        throw std::runtime_error("Successful reconstruction produced a non-finite diagnostic.");
    }
    return record;
}

void validate_population_and_moments(const double delta)
{
    const std::vector<Atom> population = make_population(delta);
    double weight_sum = 0.0;
    for (const Atom& atom : population)
    {
        weight_sum += atom.weight;
    }
    if (std::abs(weight_sum - 1.0) > population_weight_tolerance)
    {
        throw std::runtime_error("Prescribed population weights do not sum to one.");
    }

    const qmomkit::RawMomentTensor moments = exact_moments(population);
    for (const double value : moments.values)
    {
        if (!std::isfinite(value))
        {
            throw std::runtime_error("Prescribed raw moments are not finite.");
        }
    }
}

} // namespace

int cqmom_2d_borderline_conditioning::run_conditioning_workflow()
{
    try
    {
        for (const double delta : delta_values)
        {
            validate_population_and_moments(delta);
        }

        std::vector<RunRecord> records;
        records.reserve(delta_values.size() * 2 * solver_choices.size() + 2);
        const std::array<std::vector<std::size_t>, 2> conditioning_orders{{std::vector<std::size_t>{0, 1}, std::vector<std::size_t>{1, 0}}};

        std::cout << "2D CQMOM borderline-conditioning sweep\n";
        std::cout << "Running fixed-count solver sweep and adaptive delta=0 references...\n";
        for (const double delta : delta_values)
        {
            for (const std::vector<std::size_t>& conditioning_order : conditioning_orders)
            {
                for (const SolverChoice& choice : solver_choices)
                {
                    records.push_back(run_reconstruction(delta, conditioning_order, choice, false, "fixed"));
                }
            }
        }
        for (const std::vector<std::size_t>& conditioning_order : conditioning_orders)
        {
            records.push_back(
                run_reconstruction(0.0, conditioning_order, {qmomkit::VandermondeSolver::automatic, "automatic"}, true, "adaptive_reference"));
        }

        cqmom_2d_borderline_conditioning::validate_successful_rows(records);
        cqmom_2d_borderline_conditioning::validate_report_inputs(records);

        const std::filesystem::path output_directory = "output";
        std::error_code directory_error;
        std::filesystem::create_directories(output_directory, directory_error);
        if (directory_error)
        {
            throw std::runtime_error("Could not create output directory: " + directory_error.message());
        }

        cqmom_2d_borderline_conditioning::write_diagnostics_file(records, output_directory / "diagnostics.txt");
        cqmom_2d_borderline_conditioning::write_markdown_report(records, output_directory / "report.md");
        cqmom_2d_borderline_conditioning::validate_report_file(output_directory / "report.md");
        cqmom_2d_borderline_conditioning::print_terminal_summary(records);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
