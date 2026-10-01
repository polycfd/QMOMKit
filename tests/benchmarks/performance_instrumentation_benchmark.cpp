#include "internal/cqmom_case.hpp"
#include "internal/performance.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <benchmark/benchmark.h>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

namespace perf = qmomkit::internal::performance;

struct CqmomWorkload
{
    qmomkit::PhaseSpaceHierarchy hierarchy;
    qmomkit::RawMomentTensor dense_moments;
    qmomkit::CqmomMomentData lean_moments;
    qmomkit::CqmomLeanCase lean_case;
    qmomkit::CqmomCase case_;
};

struct BenchmarkSummary
{
    std::string name;
    std::string workload_kind;
    std::vector<std::size_t> node_counts;
    std::string solver_policy;
    std::string phase_space_evolution_scheme;
    std::string node_advection_scheme;
    std::string instrumentation_mode;
    std::size_t state_count = 0;
    std::size_t steps_per_state = 0;
    std::size_t operations_per_iteration = 0;
    std::string operation_label;
    std::size_t iterations = 0;
    double total_runtime_seconds = 0.0;
    std::size_t requested_terminal_nodes = 0;
    std::size_t dense_moments = 0;
    std::size_t lean_moments = 0;
    std::size_t reconstructed_moments = 0;
    std::size_t diagnostic_node_advection_steps = 0;
    std::size_t diagnostic_node_advection_rejected_steps = 0;
    std::size_t diagnostic_node_advection_rhs_evaluations = 0;
    bool has_trajectory_error = false;
    double max_trajectory_error = 0.0;
    perf::PerformanceCounters counters;
};

std::vector<BenchmarkSummary>& benchmark_summaries()
{
    static auto* summaries = new std::vector<BenchmarkSummary>();
    return *summaries;
}

std::string join_node_counts(const std::vector<std::size_t>& node_counts)
{
    std::ostringstream stream;
    for (std::size_t dimension = 0; dimension < node_counts.size(); ++dimension)
    {
        if (dimension > 0)
        {
            stream << " x ";
        }
        stream << node_counts[dimension];
    }
    return stream.str();
}

std::vector<std::size_t> nodes_2x2()
{
    return {2, 2};
}

std::vector<std::size_t> nodes_3x3()
{
    return {3, 3};
}

std::vector<std::size_t> nodes_4x4()
{
    return {4, 4};
}

std::vector<std::size_t> nodes_2x2x2()
{
    return {2, 2, 2};
}

std::vector<std::size_t> nodes_2x2x2x2()
{
    return {2, 2, 2, 2};
}

std::vector<std::size_t> nodes_2x2x2x2x2()
{
    return {2, 2, 2, 2, 2};
}

double safe_ratio(const std::size_t numerator, const std::size_t denominator)
{
    return denominator == 0 ? 0.0 : static_cast<double>(numerator) / static_cast<double>(denominator);
}

double percentage(const std::size_t numerator, const std::size_t denominator)
{
    return 100.0 * safe_ratio(numerator, denominator);
}

std::string seconds_to_text(const double seconds)
{
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(3);
    if (seconds >= 1.0)
    {
        stream << seconds << " s";
    }
    else if (seconds >= 1.0e-3)
    {
        stream << seconds * 1.0e3 << " ms";
    }
    else if (seconds >= 1.0e-6)
    {
        stream << seconds * 1.0e6 << " us";
    }
    else
    {
        stream << seconds * 1.0e9 << " ns";
    }
    return stream.str();
}

void print_counter_line(const char* label, const std::size_t value)
{
    std::cout << "    " << std::left << std::setw(40) << label << value << '\n';
}

void print_metric_line(const char* label, const double value, const char* suffix = "")
{
    std::cout << "    " << std::left << std::setw(40) << label << std::fixed << std::setprecision(3) << value << suffix << '\n';
}

void print_scientific_metric_line(const char* label, const double value)
{
    std::cout << "    " << std::left << std::setw(40) << label << std::scientific << std::setprecision(6) << value << std::fixed << '\n';
}

bool has_reconstruction_counters(const perf::PerformanceCounters& counters)
{
    return counters.reconstruction.qmom_reconstructions != 0 || counters.reconstruction.cqmom_reconstructions != 0 ||
           counters.reconstruction.cqmom_branch_solves != 0 || counters.reconstruction.wheeler_calls != 0 ||
           counters.reconstruction.wheeler_adaptive_reductions != 0 || counters.reconstruction.terminal_nodes_assembled != 0 ||
           counters.reconstruction.moment_layout_builds != 0 || counters.reconstruction.cqmom_prepared_metadata_builds != 0 ||
           counters.reconstruction.cqmom_prepared_metadata_uses != 0 || counters.reconstruction.cqmom_dynamic_metadata_builds != 0;
}

bool has_linear_algebra_counters(const perf::PerformanceCounters& counters)
{
    return counters.linear_algebra.vandermonde_solves != 0 || counters.linear_algebra.bjorck_pereyra_solves != 0 ||
           counters.linear_algebra.pivoted_qr_solves != 0 || counters.linear_algebra.svd_solves != 0 ||
           counters.linear_algebra.dense_linear_solves != 0 || counters.linear_algebra.residual_checks != 0 ||
           counters.linear_algebra.condition_estimates != 0;
}

bool has_phase_space_evolution_counters(const perf::PerformanceCounters& counters)
{
    return counters.phase_space_evolution.phase_space_evolution_macro_steps != 0 || counters.phase_space_evolution.advection_operator_calls != 0 ||
           counters.phase_space_evolution.moment_source_operator_calls != 0 || counters.phase_space_evolution.moment_reconstruction_calls != 0 ||
           counters.phase_space_evolution.post_source_inversions != 0 || counters.time_integration.source_rk_stages != 0 ||
           counters.time_integration.source_rk_retries != 0 || counters.time_integration.node_rhs_evaluations != 0 ||
           counters.time_integration.embedded_rk_accepted_steps != 0 || counters.time_integration.embedded_rk_rejected_steps != 0;
}

bool has_cqmom_topology_counters(const perf::PerformanceCounters& counters)
{
    return counters.cqmom_topology.active_count_derivations != 0 || counters.cqmom_topology.active_count_validations != 0 ||
           counters.cqmom_topology.active_count_transition_builds != 0;
}

void record_summary(BenchmarkSummary summary)
{
    std::vector<BenchmarkSummary>& summaries = benchmark_summaries();
    for (BenchmarkSummary& existing : summaries)
    {
        if (existing.name == summary.name)
        {
            existing = std::move(summary);
            return;
        }
    }
    summaries.push_back(std::move(summary));
}

class SummaryPrinter
{
public:
    ~SummaryPrinter()
    {
        const std::vector<BenchmarkSummary>& summaries = benchmark_summaries();
        if (summaries.empty())
        {
            return;
        }

        std::cout << "\nQMOMKit performance instrumentation report\n";
        std::cout << "==========================================\n";
        std::cout << "Timing tells how expensive the benchmark body was; counters explain the workload shape.\n";
        std::cout << "Use these summaries to form profiling hypotheses, not as optimization conclusions.\n";

        for (const BenchmarkSummary& summary : summaries)
        {
            const std::size_t total_operations = summary.operations_per_iteration * summary.iterations;
            const std::size_t reconstructions = summary.counters.reconstruction.cqmom_reconstructions;
            const std::size_t vandermonde_solves = summary.counters.linear_algebra.vandermonde_solves;
            const std::size_t phase_space_evolution_ops =
                summary.counters.phase_space_evolution.advection_operator_calls + summary.counters.phase_space_evolution.moment_source_operator_calls;
            const bool has_reconstruction = has_reconstruction_counters(summary.counters);
            const bool has_linear_algebra = has_linear_algebra_counters(summary.counters);
            const bool has_phase_space_evolution = has_phase_space_evolution_counters(summary.counters);
            const bool has_cqmom_topology = has_cqmom_topology_counters(summary.counters);

            std::cout << "\n" << summary.name << "\n";
            std::cout << std::string(summary.name.size(), '-') << "\n";
            std::cout << "  Workload\n";
            std::cout << "    type: " << summary.workload_kind << '\n';
            std::cout << "    phase-space dimensions: " << summary.node_counts.size() << '\n';
            std::cout << "    node counts: " << join_node_counts(summary.node_counts) << '\n';
            std::cout << "    requested terminal nodes: " << summary.requested_terminal_nodes << '\n';
            std::cout << "    dense/lean moment counts: " << summary.dense_moments << " / " << summary.lean_moments << '\n';
            if (summary.reconstructed_moments > 0)
            {
                std::cout << "    reconstructed moments: " << summary.reconstructed_moments << '\n';
            }
            if (!summary.solver_policy.empty())
            {
                std::cout << "    solver policy: " << summary.solver_policy << '\n';
            }
            if (!summary.phase_space_evolution_scheme.empty())
            {
                std::cout << "    phase-space evolution scheme: " << summary.phase_space_evolution_scheme << '\n';
            }
            if (!summary.node_advection_scheme.empty())
            {
                std::cout << "    node-advection scheme: " << summary.node_advection_scheme << '\n';
            }
            if (!summary.instrumentation_mode.empty())
            {
                std::cout << "    instrumentation mode: " << summary.instrumentation_mode << '\n';
            }
            if (summary.state_count > 0)
            {
                std::cout << "    states: " << summary.state_count << '\n';
            }
            if (summary.steps_per_state > 0)
            {
                std::cout << "    steps per state: " << summary.steps_per_state << '\n';
            }

            std::cout << "  Timing\n";
            std::cout << "    total measured benchmark-body runtime: " << seconds_to_text(summary.total_runtime_seconds) << '\n';
            if (total_operations > 0)
            {
                std::cout << "    average runtime per " << summary.operation_label << ": "
                          << seconds_to_text(summary.total_runtime_seconds / static_cast<double>(total_operations)) << '\n';
            }
            if (summary.state_count > 0 && summary.steps_per_state > 0)
            {
                std::cout << "    average runtime per state update: "
                          << seconds_to_text(summary.total_runtime_seconds / static_cast<double>(total_operations)) << '\n';
            }

            if (has_reconstruction)
            {
                std::cout << "  Reconstruction counters\n";
                print_counter_line("CQMOM reconstructions", summary.counters.reconstruction.cqmom_reconstructions);
                print_counter_line("QMOM reconstructions", summary.counters.reconstruction.qmom_reconstructions);
                print_counter_line("Wheeler calls", summary.counters.reconstruction.wheeler_calls);
                print_counter_line("CQMOM branch solves", summary.counters.reconstruction.cqmom_branch_solves);
                print_counter_line("terminal nodes assembled", summary.counters.reconstruction.terminal_nodes_assembled);
                print_counter_line("moment layout builds", summary.counters.reconstruction.moment_layout_builds);
                print_counter_line("prepared metadata builds", summary.counters.reconstruction.cqmom_prepared_metadata_builds);
                print_counter_line("prepared metadata uses", summary.counters.reconstruction.cqmom_prepared_metadata_uses);
                print_counter_line("dynamic metadata builds", summary.counters.reconstruction.cqmom_dynamic_metadata_builds);
            }

            if (has_linear_algebra)
            {
                std::cout << "  Linear algebra counters\n";
                print_counter_line("Vandermonde solves", summary.counters.linear_algebra.vandermonde_solves);
                print_counter_line("Bjorck-Pereyra solves", summary.counters.linear_algebra.bjorck_pereyra_solves);
                print_counter_line("pivoted QR solves", summary.counters.linear_algebra.pivoted_qr_solves);
                print_counter_line("SVD solves", summary.counters.linear_algebra.svd_solves);
                print_counter_line("dense linear solves", summary.counters.linear_algebra.dense_linear_solves);
                print_counter_line("residual checks", summary.counters.linear_algebra.residual_checks);
                print_counter_line("condition estimates", summary.counters.linear_algebra.condition_estimates);
            }

            if (has_phase_space_evolution)
            {
                std::cout << "  Phase-space evolution/time-integration counters\n";
                print_counter_line("advection operators", summary.counters.phase_space_evolution.advection_operator_calls);
                print_counter_line("moment-source operators", summary.counters.phase_space_evolution.moment_source_operator_calls);
                print_counter_line("moment reconstructions", summary.counters.phase_space_evolution.moment_reconstruction_calls);
                print_counter_line("post-source inversions", summary.counters.phase_space_evolution.post_source_inversions);
                print_counter_line("moment-source SSP-RK stages", summary.counters.time_integration.source_rk_stages);
                print_counter_line("moment-source retries", summary.counters.time_integration.source_rk_retries);
                print_counter_line("node RHS evaluations", summary.counters.time_integration.node_rhs_evaluations);
                print_counter_line("adaptive embedded-RK accepted/rejected", summary.counters.time_integration.embedded_rk_accepted_steps +
                                                                                 summary.counters.time_integration.embedded_rk_rejected_steps);
            }

            if (has_cqmom_topology)
            {
                std::cout << "  CQMOM topology counters\n";
                print_counter_line("active-count derivations", summary.counters.cqmom_topology.active_count_derivations);
                print_counter_line("active-count validations", summary.counters.cqmom_topology.active_count_validations);
                print_counter_line("active-count transition builds", summary.counters.cqmom_topology.active_count_transition_builds);
            }

            if (summary.diagnostic_node_advection_steps != 0 || summary.diagnostic_node_advection_rejected_steps != 0 ||
                summary.diagnostic_node_advection_rhs_evaluations != 0 || summary.has_trajectory_error)
            {
                std::cout << "  Node-advection diagnostics\n";
                print_counter_line("node-advection steps", summary.diagnostic_node_advection_steps);
                print_counter_line("node-advection rejected steps", summary.diagnostic_node_advection_rejected_steps);
                print_counter_line("node-advection RHS evaluations", summary.diagnostic_node_advection_rhs_evaluations);
                if (summary.has_trajectory_error)
                {
                    print_scientific_metric_line("max trajectory error", summary.max_trajectory_error);
                }
            }

            if (has_reconstruction || has_linear_algebra || has_phase_space_evolution || has_cqmom_topology)
            {
                std::cout << "  Derived metrics\n";
                if (has_reconstruction)
                {
                    print_metric_line("Wheeler calls / reconstruction", safe_ratio(summary.counters.reconstruction.wheeler_calls, reconstructions));
                    print_metric_line("layout builds / operation",
                                      safe_ratio(summary.counters.reconstruction.moment_layout_builds, total_operations));
                }
                if (has_linear_algebra)
                {
                    print_metric_line("Vandermonde solves / reconstruction", safe_ratio(vandermonde_solves, reconstructions));
                    print_metric_line("residual checks / Vandermonde solve",
                                      safe_ratio(summary.counters.linear_algebra.residual_checks, vandermonde_solves));
                }
                if (has_phase_space_evolution && total_operations > 0)
                {
                    print_metric_line("phase-space evolution operators / update", safe_ratio(phase_space_evolution_ops, total_operations));
                    print_metric_line("node RHS evaluations / update",
                                      safe_ratio(summary.counters.time_integration.node_rhs_evaluations, total_operations));
                    print_metric_line("moment-source SSP-RK stages / update",
                                      safe_ratio(summary.counters.time_integration.source_rk_stages, total_operations));
                }
                if (has_cqmom_topology && total_operations > 0)
                {
                    print_metric_line("active-count derivations / update",
                                      safe_ratio(summary.counters.cqmom_topology.active_count_derivations, total_operations));
                    print_metric_line("active-count validations / update",
                                      safe_ratio(summary.counters.cqmom_topology.active_count_validations, total_operations));
                    print_metric_line("active-count transition builds / update",
                                      safe_ratio(summary.counters.cqmom_topology.active_count_transition_builds, total_operations));
                }
            }
        }
        std::cout << '\n';
    }
};

SummaryPrinter summary_printer;

std::vector<double> separated_nodes(const std::size_t node_count, const std::size_t dimension)
{
    std::vector<double> nodes(node_count, 0.0);
    if (node_count == 1)
    {
        nodes.front() = 0.05 * static_cast<double>(dimension + 1);
        return nodes;
    }

    for (std::size_t node = 0; node < node_count; ++node)
    {
        const double fraction = static_cast<double>(node) / static_cast<double>(node_count - 1);
        nodes[node] = -0.9 + 1.8 * fraction + 0.04 * static_cast<double>(dimension);
    }
    return nodes;
}

std::vector<double> normalized_weights(const std::size_t node_count, const std::size_t dimension)
{
    std::vector<double> weights(node_count, 0.0);
    double sum = 0.0;
    for (std::size_t node = 0; node < node_count; ++node)
    {
        weights[node] = 1.0 + static_cast<double>((node + 1) * (dimension + 2));
        sum += weights[node];
    }
    for (double& weight : weights)
    {
        weight /= sum;
    }
    return weights;
}

std::vector<double> one_dimensional_raw_moments(const std::vector<double>& nodes, const std::vector<double>& weights, const std::size_t extent)
{
    std::vector<double> moments(extent, 0.0);
    for (std::size_t order = 0; order < extent; ++order)
    {
        double value = 0.0;
        for (std::size_t node = 0; node < nodes.size(); ++node)
        {
            value += weights[node] * std::pow(nodes[node], static_cast<int>(order));
        }
        moments[order] = value;
    }
    return moments;
}

qmomkit::RawMomentTensor independent_product_moments(const std::vector<std::size_t>& node_counts)
{
    std::vector<std::size_t> extents;
    extents.reserve(node_counts.size());
    for (const std::size_t node_count : node_counts)
    {
        extents.push_back(2 * node_count);
    }

    std::vector<std::vector<double>> moments_by_dimension;
    moments_by_dimension.reserve(node_counts.size());
    for (std::size_t dimension = 0; dimension < node_counts.size(); ++dimension)
    {
        moments_by_dimension.push_back(one_dimensional_raw_moments(separated_nodes(node_counts[dimension], dimension),
                                                                   normalized_weights(node_counts[dimension], dimension), extents[dimension]));
    }

    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat = 0; flat < moments.size(); ++flat)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat, extents);
        double value = 1.0;
        for (std::size_t dimension = 0; dimension < node_counts.size(); ++dimension)
        {
            value *= moments_by_dimension[dimension][moment_index[dimension]];
        }
        moments[flat] = value;
    }

    return {std::move(extents), std::move(moments)};
}

CqmomWorkload make_workload(const std::vector<std::size_t>& node_counts,
                            const qmomkit::VandermondeSolver solver = qmomkit::VandermondeSolver::automatic,
                            const qmomkit::CqmomDiagnosticLevel diagnostic_level = qmomkit::CqmomDiagnosticLevel::basic)
{
    qmomkit::PhaseSpaceHierarchy hierarchy(node_counts);
    qmomkit::RawMomentTensor dense_moments = independent_product_moments(node_counts);
    qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::CqmomMomentData lean_moments = qmomkit::extract_cqmom_moments(dense_moments, layout);
    qmomkit::CqmomLeanCase lean_case{lean_moments, {}, {}, false, diagnostic_level, solver};
    qmomkit::CqmomCase case_ = qmomkit::internal::CqmomConstructionAccess::create_case(lean_case);
    return {std::move(hierarchy), std::move(dense_moments), std::move(lean_moments), std::move(lean_case), std::move(case_)};
}

qmomkit::RawMomentTensor ragged_phase_space_evolution_moments()
{
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.5};
    const std::vector<std::vector<double>> abscissas{{-1.0, 1.0, 1.0}, {0.0, -0.5, 0.5}};
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][node], static_cast<int>(moment_index[dimension]));
            }
            moments[flat_moment] += contribution;
        }
    }
    return {extents, std::move(moments)};
}

qmomkit::RawMomentTensor ragged_phase_space_evolution_moments_3d()
{
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.25, 0.25, 0.5};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0}, {-1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}};
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][node], static_cast<int>(moment_index[dimension]));
            }
            moments[flat_moment] += contribution;
        }
    }
    return {extents, std::move(moments)};
}

CqmomWorkload make_phase_space_evolution_workload(const std::vector<std::size_t>& node_counts, const bool adaptive)
{
    qmomkit::PhaseSpaceHierarchy hierarchy(node_counts);
    qmomkit::RawMomentTensor dense_moments = independent_product_moments(node_counts);
    if (adaptive)
    {
        if (node_counts == nodes_2x2())
        {
            dense_moments = ragged_phase_space_evolution_moments();
        }
        else if (node_counts == nodes_2x2x2())
        {
            dense_moments = ragged_phase_space_evolution_moments_3d();
        }
        else
        {
            throw std::invalid_argument("No adaptive CQMOM phase-space evolution benchmark fixture is available for these node counts.");
        }
    }
    qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::CqmomMomentData lean_moments = qmomkit::extract_cqmom_moments(dense_moments, layout);
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = adaptive;
    qmom_options.rmin = 1.0e-10;
    qmom_options.eabs = 1.0e-12;
    qmom_options.cutoff = 1.0e-12;
    qmomkit::CqmomLeanCase lean_case{lean_moments, {}, qmom_options, adaptive, qmomkit::CqmomDiagnosticLevel::basic};
    qmomkit::CqmomCase case_ = qmomkit::internal::CqmomConstructionAccess::create_case(lean_case);
    return {std::move(hierarchy), std::move(dense_moments), std::move(lean_moments), std::move(lean_case), std::move(case_)};
}

std::vector<qmomkit::CqmomState> make_states(const qmomkit::CqmomCase& case_, const std::size_t state_count)
{
    std::vector<qmomkit::CqmomState> states;
    states.reserve(state_count);
    for (std::size_t state = 0; state < state_count; ++state)
    {
        states.push_back(case_.initial_state());
    }
    return states;
}

qmomkit::CqmomState state_without_accepted_quadrature(const qmomkit::CqmomState& state)
{
    const auto& prepared_metadata = qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(state);
    const qmomkit::CqmomLeanCase* lean_case = qmomkit::internal::CqmomConstructionAccess::state_lean_case(state);
    if (lean_case == nullptr)
    {
        throw std::logic_error("CQMOM phase-space evolution benchmark state has unknown storage.");
    }
    return qmomkit::internal::CqmomConstructionAccess::create_state(*lean_case, state.active_branch_counts(), prepared_metadata);
}

std::vector<qmomkit::CqmomState> make_phase_space_evolution_states(const qmomkit::CqmomCase& case_, const std::size_t state_count,
                                                                   const bool cache_miss)
{
    std::vector<qmomkit::CqmomState> states = make_states(case_, state_count);
    if (cache_miss)
    {
        for (qmomkit::CqmomState& state : states)
        {
            state = state_without_accepted_quadrature(state);
        }
    }
    return states;
}

enum class NodeTrajectory
{
    constant,
    affine
};

// Written in the packed form the examples use, so the benchmark measures the
// same velocity path as a real consumer. The per-node adapter allocates one
// vector per node per stage, which is a cost no shipped example pays;
// allocation_benchmark.cpp compares the two forms deliberately.
qmomkit::PhaseSpaceVelocity velocity_model()
{
    return [](const double, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
    {
        for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
        {
            const std::size_t first = node * dimensions;
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
            {
                velocities[first + dimension] = 0.015 * static_cast<double>(dimension + 1) + 0.01 * coordinates[first + dimension];
            }
        }
    };
}

qmomkit::Quadrature1D make_qmom_node_advection_fixture(const std::size_t node_count)
{
    return {separated_nodes(node_count, 0), normalized_weights(node_count, 0), qmomkit::QuadratureStatus::success, ""};
}

[[nodiscard]] double constant_velocity_component(const std::size_t dimension)
{
    return 0.12 + 0.04 * static_cast<double>(dimension + 1);
}

[[nodiscard]] double affine_offset_component(const std::size_t dimension)
{
    return 0.05 + 0.03 * static_cast<double>(dimension + 1);
}

[[nodiscard]] double affine_rate_component(const std::size_t dimension)
{
    return 0.18 + 0.04 * static_cast<double>(dimension + 1);
}

qmomkit::PhaseSpaceVelocity node_advection_velocity_model(const NodeTrajectory trajectory)
{
    return [trajectory](const double, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
    {
        for (std::size_t node = 0; node < coordinates.size() / dimensions; ++node)
        {
            const std::size_t first = node * dimensions;
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
            {
                if (trajectory == NodeTrajectory::constant)
                {
                    velocities[first + dimension] = constant_velocity_component(dimension);
                }
                else
                {
                    velocities[first + dimension] =
                        affine_offset_component(dimension) + affine_rate_component(dimension) * coordinates[first + dimension];
                }
            }
        }
    };
}

[[nodiscard]] double exact_node_coordinate(const NodeTrajectory trajectory, const std::size_t dimension, const double initial_coordinate,
                                           const double time_step)
{
    if (trajectory == NodeTrajectory::constant)
    {
        return initial_coordinate + constant_velocity_component(dimension) * time_step;
    }

    const double offset = affine_offset_component(dimension);
    const double rate = affine_rate_component(dimension);
    return (initial_coordinate + offset / rate) * std::exp(rate * time_step) - offset / rate;
}

[[nodiscard]] double max_qmom_trajectory_error(const qmomkit::Quadrature1D& initial, const qmomkit::Quadrature1D& advanced,
                                               const NodeTrajectory trajectory, const double time_step)
{
    double error = 0.0;
    for (std::size_t node = 0; node < initial.abscissas.size(); ++node)
    {
        const double expected = exact_node_coordinate(trajectory, 0, initial.abscissas[node], time_step);
        error = std::max(error, std::abs(advanced.abscissas[node] - expected));
    }
    return error;
}

[[nodiscard]] double max_cqmom_trajectory_error(const qmomkit::CqmomQuadrature& initial, const qmomkit::CqmomQuadrature& advanced,
                                                const NodeTrajectory trajectory, const double time_step)
{
    double error = 0.0;
    for (std::size_t dimension = 0; dimension < initial.number_of_dimensions(); ++dimension)
    {
        for (std::size_t node = 0; node < initial.weights().size(); ++node)
        {
            const double expected = exact_node_coordinate(trajectory, dimension, initial.abscissas()[dimension][node], time_step);
            error = std::max(error, std::abs(advanced.abscissas()[dimension][node] - expected));
        }
    }
    return error;
}

qmomkit::NodeAdvectionOptions node_advection_benchmark_options(const qmomkit::NodeAdvectionScheme scheme, const double time_step)
{
    qmomkit::NodeAdvectionOptions options;
    options.scheme = scheme;
    options.adaptive_options.relative_tolerance = 1.0e-8;
    options.adaptive_options.absolute_tolerance = 1.0e-10;
    options.adaptive_options.initial_step = time_step;
    options.adaptive_options.max_step = time_step;
    return options;
}

qmomkit::CqmomMomentSource source_model()
{
    return [](const double, const qmomkit::CqmomSourceContext& context)
    {
        qmomkit::CqmomMomentRate rate;
        rate.values.assign(context.size(), 0.0);
        for (std::size_t moment = 0; moment < context.size(); ++moment)
        {
            std::size_t order_sum = 0;
            for (const std::size_t order : context.index(moment))
            {
                order_sum += order;
            }
            const double damping = 0.001 * static_cast<double>(order_sum + 1);
            rate.values[moment] = -damping * context.value(moment);
        }
        return rate;
    };
}

qmomkit::CqmomMomentBoundaryContribution zero_boundary_model()
{
    return [](const double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };
}

[[nodiscard]] bool instrumentation_is_available(benchmark::State& state)
{
    if (!perf::enabled())
    {
        state.SkipWithError("Rebuild QMOMKit with QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION=ON.");
        return false;
    }
    return true;
}

enum class CounterMode
{
    enabled_session,
    enabled_no_session,
    timing_only
};

enum class CqmomPhaseSpaceEvolutionWorkflow
{
    advection,
    moment_source,
    advsrc,
    srcadv,
    strang
};

std::string phase_space_evolution_workflow_name(const CqmomPhaseSpaceEvolutionWorkflow workflow)
{
    switch (workflow)
    {
    case CqmomPhaseSpaceEvolutionWorkflow::advection:
        return "advection";
    case CqmomPhaseSpaceEvolutionWorkflow::moment_source:
        return "moment_source";
    case CqmomPhaseSpaceEvolutionWorkflow::advsrc:
        return "advsrc";
    case CqmomPhaseSpaceEvolutionWorkflow::srcadv:
        return "srcadv";
    case CqmomPhaseSpaceEvolutionWorkflow::strang:
        return "strang";
    }
    return "unknown";
}

std::string instrumentation_mode_name(const CounterMode mode)
{
    if (!perf::enabled())
    {
        return "compile-time disabled";
    }
    switch (mode)
    {
    case CounterMode::enabled_session:
        return "enabled session";
    case CounterMode::enabled_no_session:
        return "enabled build, no active session";
    case CounterMode::timing_only:
        return "timing only";
    }
    return "timing only";
}

std::string solver_name(const qmomkit::VandermondeSolver solver)
{
    switch (solver)
    {
    case qmomkit::VandermondeSolver::automatic:
        return "automatic";
    case qmomkit::VandermondeSolver::dense:
        return "dense";
    case qmomkit::VandermondeSolver::pivoted_qr:
        return "pivoted_qr";
    case qmomkit::VandermondeSolver::svd:
        return "svd";
    case qmomkit::VandermondeSolver::moore_penrose:
        return "moore_penrose";
    case qmomkit::VandermondeSolver::bjorck_pereyra:
        return "bjorck_pereyra";
    }
    return "unknown";
}

std::string diagnostic_name(const qmomkit::CqmomDiagnosticLevel level)
{
    switch (level)
    {
    case qmomkit::CqmomDiagnosticLevel::none:
        return "none";
    case qmomkit::CqmomDiagnosticLevel::basic:
        return "basic";
    case qmomkit::CqmomDiagnosticLevel::detailed:
        return "detailed";
    }
    return "unknown";
}

std::string phase_space_evolution_scheme_name(const qmomkit::PhaseSpaceEvolutionScheme scheme)
{
    switch (scheme)
    {
    case qmomkit::PhaseSpaceEvolutionScheme::strang:
        return "strang";
    case qmomkit::PhaseSpaceEvolutionScheme::advsrc:
        return "advsrc";
    case qmomkit::PhaseSpaceEvolutionScheme::srcadv:
        return "srcadv";
    }
    return "unknown";
}

std::string node_advection_scheme_name(const qmomkit::NodeAdvectionScheme scheme)
{
    switch (scheme)
    {
    case qmomkit::NodeAdvectionScheme::ssp_rk3:
        return "ssp_rk3";
    case qmomkit::NodeAdvectionScheme::rk23:
        return "rk23";
    case qmomkit::NodeAdvectionScheme::rk45:
        return "rk45";
    }
    return "unknown";
}

std::string node_trajectory_name(const NodeTrajectory trajectory)
{
    switch (trajectory)
    {
    case NodeTrajectory::constant:
        return "constant";
    case NodeTrajectory::affine:
        return "affine";
    }
    return "unknown";
}

void publish_counter(benchmark::State& state, const char* name, const std::size_t value)
{
    state.counters[name] = benchmark::Counter(static_cast<double>(value), benchmark::Counter::kAvgIterations);
}

void publish_workload_counters(benchmark::State& state, const perf::PerformanceCounters& counters)
{
    if (has_reconstruction_counters(counters))
    {
        publish_counter(state, "cqmom_reconstructions", counters.reconstruction.cqmom_reconstructions);
        publish_counter(state, "branch_solves", counters.reconstruction.cqmom_branch_solves);
        publish_counter(state, "wheeler_calls", counters.reconstruction.wheeler_calls);
        publish_counter(state, "terminal_nodes", counters.reconstruction.terminal_nodes_assembled);
        publish_counter(state, "moment_layout_builds", counters.reconstruction.moment_layout_builds);
        publish_counter(state, "prepared_metadata_builds", counters.reconstruction.cqmom_prepared_metadata_builds);
        publish_counter(state, "prepared_metadata_uses", counters.reconstruction.cqmom_prepared_metadata_uses);
        publish_counter(state, "dynamic_metadata_builds", counters.reconstruction.cqmom_dynamic_metadata_builds);
    }
    if (has_linear_algebra_counters(counters))
    {
        publish_counter(state, "vandermonde_solves", counters.linear_algebra.vandermonde_solves);
        publish_counter(state, "bp_solves", counters.linear_algebra.bjorck_pereyra_solves);
        publish_counter(state, "qr_solves", counters.linear_algebra.pivoted_qr_solves);
        publish_counter(state, "svd_solves", counters.linear_algebra.svd_solves);
        publish_counter(state, "dense_solves", counters.linear_algebra.dense_linear_solves);
        publish_counter(state, "residual_checks", counters.linear_algebra.residual_checks);
        publish_counter(state, "condition_estimates", counters.linear_algebra.condition_estimates);
    }
    if (has_phase_space_evolution_counters(counters))
    {
        publish_counter(state, "advection_ops", counters.phase_space_evolution.advection_operator_calls);
        publish_counter(state, "moment_source_ops", counters.phase_space_evolution.moment_source_operator_calls);
        publish_counter(state, "moment_reconstructions", counters.phase_space_evolution.moment_reconstruction_calls);
        publish_counter(state, "post_source_inversions", counters.phase_space_evolution.post_source_inversions);
        publish_counter(state, "source_rk_stages", counters.time_integration.source_rk_stages);
        publish_counter(state, "moment_source_retries", counters.time_integration.source_rk_retries);
        publish_counter(state, "node_rhs_evals", counters.time_integration.node_rhs_evaluations);
        publish_counter(state, "embedded_rk_accepted", counters.time_integration.embedded_rk_accepted_steps);
        publish_counter(state, "embedded_rk_rejected", counters.time_integration.embedded_rk_rejected_steps);
    }
    if (has_cqmom_topology_counters(counters))
    {
        publish_counter(state, "active_count_derivations", counters.cqmom_topology.active_count_derivations);
        publish_counter(state, "active_count_validations", counters.cqmom_topology.active_count_validations);
        publish_counter(state, "active_count_transition_builds", counters.cqmom_topology.active_count_transition_builds);
    }
}

void publish_node_advection_diagnostics(benchmark::State& state, const std::size_t step_count, const std::size_t rejected_step_count,
                                        const std::size_t rhs_evaluations, const double max_trajectory_error)
{
    publish_counter(state, "diagnostic_node_steps", step_count);
    publish_counter(state, "diagnostic_node_rejected", rejected_step_count);
    publish_counter(state, "diagnostic_node_rhs", rhs_evaluations);
    state.counters["max_trajectory_error"] = max_trajectory_error;
}

void publish_normalized_counters(benchmark::State& state, const perf::PerformanceCounters& counters, const std::size_t operations_per_iteration)
{
    if (operations_per_iteration == 0 || (!has_reconstruction_counters(counters) && !has_linear_algebra_counters(counters) &&
                                          !has_phase_space_evolution_counters(counters) && !has_cqmom_topology_counters(counters)))
    {
        return;
    }
    const double operations = static_cast<double>(operations_per_iteration);
    if (has_reconstruction_counters(counters))
    {
        state.counters["layout_builds_per_op"] =
            benchmark::Counter(static_cast<double>(counters.reconstruction.moment_layout_builds) / operations, benchmark::Counter::kAvgIterations);
        state.counters["prepared_metadata_uses_per_op"] = benchmark::Counter(
            static_cast<double>(counters.reconstruction.cqmom_prepared_metadata_uses) / operations, benchmark::Counter::kAvgIterations);
        state.counters["dynamic_metadata_builds_per_op"] = benchmark::Counter(
            static_cast<double>(counters.reconstruction.cqmom_dynamic_metadata_builds) / operations, benchmark::Counter::kAvgIterations);
        state.counters["wheeler_per_op"] =
            benchmark::Counter(static_cast<double>(counters.reconstruction.wheeler_calls) / operations, benchmark::Counter::kAvgIterations);
    }
    if (has_linear_algebra_counters(counters))
    {
        state.counters["vandermonde_per_op"] =
            benchmark::Counter(static_cast<double>(counters.linear_algebra.vandermonde_solves) / operations, benchmark::Counter::kAvgIterations);
        state.counters["residual_per_vandermonde"] =
            counters.linear_algebra.vandermonde_solves == 0
                ? 0.0
                : static_cast<double>(counters.linear_algebra.residual_checks) / static_cast<double>(counters.linear_algebra.vandermonde_solves);
    }
    if (has_phase_space_evolution_counters(counters))
    {
        state.counters["node_rhs_per_op"] =
            benchmark::Counter(static_cast<double>(counters.time_integration.node_rhs_evaluations) / operations, benchmark::Counter::kAvgIterations);
        state.counters["source_rk_stages_per_op"] =
            benchmark::Counter(static_cast<double>(counters.time_integration.source_rk_stages) / operations, benchmark::Counter::kAvgIterations);
    }
    if (has_cqmom_topology_counters(counters))
    {
        state.counters["active_count_derivations_per_op"] = benchmark::Counter(
            static_cast<double>(counters.cqmom_topology.active_count_derivations) / operations, benchmark::Counter::kAvgIterations);
        state.counters["active_count_validations_per_op"] = benchmark::Counter(
            static_cast<double>(counters.cqmom_topology.active_count_validations) / operations, benchmark::Counter::kAvgIterations);
        state.counters["active_count_transition_builds_per_op"] = benchmark::Counter(
            static_cast<double>(counters.cqmom_topology.active_count_transition_builds) / operations, benchmark::Counter::kAvgIterations);
    }
}

void publish_geometry_counters(benchmark::State& state, const std::vector<std::size_t>& node_counts, const CqmomWorkload& workload)
{
    state.counters["dimensions"] = static_cast<double>(node_counts.size());
    state.counters["requested_terminal_nodes"] = static_cast<double>(workload.hierarchy.total_nodes());
    state.counters["dense_moments"] = static_cast<double>(workload.dense_moments.values.size());
    state.counters["lean_moments"] = static_cast<double>(workload.lean_moments.values.size());
}

BenchmarkSummary make_summary(const std::string& name, const std::string& workload_kind, const std::vector<std::size_t>& node_counts,
                              const CqmomWorkload& workload, const perf::PerformanceCounters& counters, const double measured_seconds,
                              const std::size_t iterations, const std::size_t operations_per_iteration, const std::string& operation_label,
                              const CounterMode counter_mode)
{
    BenchmarkSummary summary;
    summary.name = name;
    summary.workload_kind = workload_kind;
    summary.node_counts = node_counts;
    summary.instrumentation_mode = instrumentation_mode_name(counter_mode);
    summary.operations_per_iteration = operations_per_iteration;
    summary.operation_label = operation_label;
    summary.iterations = iterations;
    summary.total_runtime_seconds = measured_seconds;
    summary.requested_terminal_nodes = workload.hierarchy.total_nodes();
    summary.dense_moments = workload.dense_moments.values.size();
    summary.lean_moments = workload.lean_moments.values.size();
    summary.counters = counters;
    return summary;
}

void bm_cqmom_phase_space_evolution_workload(benchmark::State& state, const qmomkit::PhaseSpaceEvolutionScheme scheme, const CounterMode counter_mode,
                                             const std::string name)
{
    if (counter_mode == CounterMode::enabled_session && !instrumentation_is_available(state))
    {
        return;
    }

    const std::size_t state_count = static_cast<std::size_t>(state.range(0));
    const std::size_t steps_per_state = static_cast<std::size_t>(state.range(1));
    const CqmomWorkload workload = make_workload({2, 2});
    const qmomkit::PhaseSpaceVelocity velocity = velocity_model();
    const qmomkit::CqmomMomentSource source = source_model();
    const qmomkit::CqmomMomentBoundaryContribution boundary = zero_boundary_model();
    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = scheme;
    options.node_options.adaptive_options.initial_step = 0.0025;
    options.node_options.adaptive_options.max_step = 0.0025;

    perf::CounterSession session;
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        state.PauseTiming();
        std::vector<qmomkit::CqmomState> states = make_states(workload.case_, state_count);
        state.ResumeTiming();

        const auto run_steps = [&]() -> bool
        {
            for (std::size_t step = 0; step < steps_per_state; ++step)
            {
                const double time = 0.0025 * static_cast<double>(step);
                for (qmomkit::CqmomState& cqmom_state : states)
                {
                    qmomkit::CqmomPhaseSpaceEvolutionResult result =
                        qmomkit::advance_cqmom_phase_space(cqmom_state, time, 0.0025, velocity, source, boundary, options);
                    if (!result.ok())
                    {
                        state.SkipWithError(result.message.empty() ? "CQMOM phase-space evolution workload failed." : result.message.c_str());
                        return false;
                    }
                    cqmom_state = std::move(result.state);
                    benchmark::DoNotOptimize(cqmom_state.hierarchy().number_of_dimensions());
                }
            }
            return true;
        };

        const auto start = std::chrono::steady_clock::now();
        bool ok = true;
        if (counter_mode == CounterMode::enabled_session)
        {
            perf::CounterScope scope(session);
            ok = run_steps();
        }
        else
        {
            ok = run_steps();
        }
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
        if (!ok)
        {
            return;
        }
    }

    const std::size_t state_steps = state_count * steps_per_state;
    const perf::PerformanceCounters counters = session.snapshot();
    publish_counter(state, "state_steps", state_steps * state.iterations());
    publish_geometry_counters(state, {2, 2}, workload);
    publish_workload_counters(state, counters);
    publish_normalized_counters(state, counters, state_steps);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations() * state_steps));
    BenchmarkSummary summary = make_summary(name, "CQMOM phase-space evolution", {2, 2}, workload, counters, measured_seconds,
                                            static_cast<std::size_t>(state.iterations()), state_steps, "state-step", counter_mode);
    summary.phase_space_evolution_scheme = phase_space_evolution_scheme_name(scheme);
    summary.state_count = state_count;
    summary.steps_per_state = steps_per_state;
    record_summary(std::move(summary));
}

void bm_cqmom_phase_space_evolution_ownership_workload(benchmark::State& state, const CqmomPhaseSpaceEvolutionWorkflow workflow,
                                                       const std::vector<std::size_t>& node_counts, const bool adaptive, const bool cache_miss,
                                                       const std::string& name)
{
    if (!instrumentation_is_available(state))
    {
        return;
    }

    const std::size_t state_count = static_cast<std::size_t>(state.range(0));
    const std::size_t steps_per_state = static_cast<std::size_t>(state.range(1));
    const CqmomWorkload workload = make_phase_space_evolution_workload(node_counts, adaptive);
    const qmomkit::PhaseSpaceVelocity velocity = velocity_model();
    const qmomkit::CqmomMomentSource source = source_model();
    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.node_options.adaptive_options.initial_step = 0.0025;
    options.node_options.adaptive_options.max_step = 0.0025;
    switch (workflow)
    {
    case CqmomPhaseSpaceEvolutionWorkflow::advsrc:
        options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;
        break;
    case CqmomPhaseSpaceEvolutionWorkflow::srcadv:
        options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
        break;
    case CqmomPhaseSpaceEvolutionWorkflow::strang:
    case CqmomPhaseSpaceEvolutionWorkflow::advection:
    case CqmomPhaseSpaceEvolutionWorkflow::moment_source:
        options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
        break;
    }

    perf::CounterSession session;
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        state.PauseTiming();
        std::vector<qmomkit::CqmomState> states = make_phase_space_evolution_states(workload.case_, state_count, cache_miss);
        state.ResumeTiming();

        const auto run_steps = [&]() -> bool
        {
            for (std::size_t step = 0; step < steps_per_state; ++step)
            {
                const double time = 0.0025 * static_cast<double>(step);
                for (qmomkit::CqmomState& cqmom_state : states)
                {
                    qmomkit::CqmomPhaseSpaceEvolutionResult result = [&]()
                    {
                        switch (workflow)
                        {
                        case CqmomPhaseSpaceEvolutionWorkflow::advection:
                            return qmomkit::advance_cqmom_advection(cqmom_state, time, 0.0025, velocity, options);
                        case CqmomPhaseSpaceEvolutionWorkflow::moment_source:
                            return qmomkit::advance_cqmom_moment_source(cqmom_state, time, 0.0025, source, options);
                        case CqmomPhaseSpaceEvolutionWorkflow::advsrc:
                        case CqmomPhaseSpaceEvolutionWorkflow::srcadv:
                        case CqmomPhaseSpaceEvolutionWorkflow::strang:
                            return qmomkit::advance_cqmom_phase_space(cqmom_state, time, 0.0025, velocity, source, options);
                        }
                        throw std::logic_error("Unsupported CQMOM phase-space evolution benchmark workflow.");
                    }();
                    if (!result.ok())
                    {
                        state.SkipWithError(result.message.empty() ? "CQMOM phase-space evolution ownership workload failed."
                                                                   : result.message.c_str());
                        return false;
                    }
                    cqmom_state = std::move(result.state);
                    benchmark::DoNotOptimize(cqmom_state.hierarchy().number_of_dimensions());
                }
            }
            return true;
        };

        const auto start = std::chrono::steady_clock::now();
        bool ok = false;
        {
            perf::CounterScope scope(session);
            ok = run_steps();
        }
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
        if (!ok)
        {
            return;
        }
    }

    const std::size_t state_steps = state_count * steps_per_state;
    const perf::PerformanceCounters counters = session.snapshot();
    publish_counter(state, "state_steps", state_steps * state.iterations());
    publish_geometry_counters(state, node_counts, workload);
    publish_workload_counters(state, counters);
    publish_normalized_counters(state, counters, state_steps);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations() * state_steps));
    BenchmarkSummary summary =
        make_summary(name, "CQMOM accepted-quadrature ownership phase-space evolution", node_counts, workload, counters, measured_seconds,
                     static_cast<std::size_t>(state.iterations()), state_steps, "state-step", CounterMode::enabled_session);
    summary.phase_space_evolution_scheme = phase_space_evolution_workflow_name(workflow);
    summary.state_count = state_count;
    summary.steps_per_state = steps_per_state;
    record_summary(std::move(summary));
}

void bm_cqmom_reconstruction_scaling(benchmark::State& state, const std::vector<std::size_t> node_counts, const qmomkit::VandermondeSolver solver,
                                     const qmomkit::CqmomDiagnosticLevel diagnostic_level, const CounterMode counter_mode, const std::string name)
{
    if (counter_mode == CounterMode::enabled_session && !instrumentation_is_available(state))
    {
        return;
    }

    const CqmomWorkload workload = make_workload(node_counts, solver, diagnostic_level);
    perf::CounterSession session;
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        const auto start = std::chrono::steady_clock::now();
        const auto compute = [&]() { return qmomkit::compute_cqmom(workload.lean_case); };
        std::optional<qmomkit::CqmomResult> result;
        if (counter_mode == CounterMode::enabled_session)
        {
            perf::CounterScope scope(session);
            result.emplace(compute());
        }
        else
        {
            result.emplace(compute());
        }
        if (!result->ok())
        {
            state.SkipWithError(result->message.empty() ? "CQMOM reconstruction workload failed." : result->message.c_str());
            return;
        }
        benchmark::DoNotOptimize(result->quadrature.weights().data());
        benchmark::DoNotOptimize(result->quadrature.abscissas().data());
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
    }

    const perf::PerformanceCounters counters = session.snapshot();
    publish_geometry_counters(state, node_counts, workload);
    publish_workload_counters(state, counters);
    publish_normalized_counters(state, counters, 1);
    BenchmarkSummary summary = make_summary(name, "CQMOM reconstruction", node_counts, workload, counters, measured_seconds,
                                            static_cast<std::size_t>(state.iterations()), 1, "reconstruction", counter_mode);
    summary.solver_policy = solver_name(solver) + ", diagnostics=" + diagnostic_name(diagnostic_level);
    record_summary(std::move(summary));
}

void bm_cqmom_normal_prepared_reconstruction_scaling(benchmark::State& state, const std::vector<std::size_t>& node_counts, const std::string name)
{
    if (!instrumentation_is_available(state))
    {
        return;
    }

    const CqmomWorkload workload = make_workload(node_counts);
    perf::CounterSession session;
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        const auto start = std::chrono::steady_clock::now();
        std::optional<qmomkit::CqmomResult> result;
        {
            perf::CounterScope scope(session);
            result.emplace(qmomkit::compute_cqmom(workload.case_));
        }
        if (!result->ok())
        {
            state.SkipWithError(result->message.empty() ? "Normal prepared CQMOM reconstruction workload failed." : result->message.c_str());
            return;
        }
        benchmark::DoNotOptimize(result->quadrature.weights().data());
        benchmark::DoNotOptimize(result->quadrature.abscissas().data());
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
    }

    const perf::PerformanceCounters counters = session.snapshot();
    publish_geometry_counters(state, node_counts, workload);
    publish_workload_counters(state, counters);
    publish_normalized_counters(state, counters, 1);
    BenchmarkSummary summary = make_summary(name, "CQMOM normal prepared reconstruction", node_counts, workload, counters, measured_seconds,
                                            static_cast<std::size_t>(state.iterations()), 1, "reconstruction", CounterMode::enabled_session);
    summary.solver_policy = "automatic, diagnostics=basic";
    record_summary(std::move(summary));
}

void bm_cqmom_fixed_case_preparation(benchmark::State& state, const std::vector<std::size_t> node_counts, const std::string name)
{
    if (!instrumentation_is_available(state))
    {
        return;
    }

    const qmomkit::PhaseSpaceHierarchy hierarchy(node_counts);
    const qmomkit::RawMomentTensor dense_moments = independent_product_moments(node_counts);
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData lean_moments = qmomkit::extract_cqmom_moments(dense_moments, layout);
    const CqmomWorkload workload = make_workload(node_counts);
    perf::CounterSession session;
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        const auto start = std::chrono::steady_clock::now();
        std::optional<qmomkit::CqmomCase> prepared_case;
        {
            perf::CounterScope scope(session);
            prepared_case.emplace(qmomkit::internal::CqmomConstructionAccess::create_case(
                qmomkit::CqmomLeanCase{lean_moments, {}, {}, false, qmomkit::CqmomDiagnosticLevel::basic}));
        }
        benchmark::DoNotOptimize(prepared_case->hierarchy().number_of_dimensions());
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
    }

    const perf::PerformanceCounters counters = session.snapshot();
    publish_geometry_counters(state, node_counts, workload);
    publish_workload_counters(state, counters);
    publish_normalized_counters(state, counters, 1);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
    BenchmarkSummary summary = make_summary(name, "CQMOM programmatic fixed case construction", node_counts, workload, counters, measured_seconds,
                                            static_cast<std::size_t>(state.iterations()), 1, "case construction", CounterMode::enabled_session);
    record_summary(std::move(summary));
}

void bm_cqmom_initial_state_prepared(benchmark::State& state, const std::vector<std::size_t> node_counts, const std::string name)
{
    if (!instrumentation_is_available(state))
    {
        return;
    }

    const CqmomWorkload workload = make_workload(node_counts);
    perf::CounterSession session;
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        const auto start = std::chrono::steady_clock::now();
        std::optional<qmomkit::CqmomState> initial_state;
        {
            perf::CounterScope scope(session);
            initial_state.emplace(workload.case_.initial_state());
        }
        benchmark::DoNotOptimize(initial_state->moment_values().data());
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
    }

    const perf::PerformanceCounters counters = session.snapshot();
    publish_geometry_counters(state, node_counts, workload);
    publish_workload_counters(state, counters);
    publish_normalized_counters(state, counters, 1);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
    BenchmarkSummary summary = make_summary(name, "CQMOM prepared initial-state construction", node_counts, workload, counters, measured_seconds,
                                            static_cast<std::size_t>(state.iterations()), 1, "initial state", CounterMode::enabled_session);
    record_summary(std::move(summary));
}

void bm_cqmom_explicit_reconstruction_control(benchmark::State& state, const std::string name)
{
    if (!instrumentation_is_available(state))
    {
        return;
    }

    const CqmomWorkload workload = make_workload(nodes_2x2());
    perf::CounterSession session;
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        const auto start = std::chrono::steady_clock::now();
        std::optional<qmomkit::CqmomResult> result;
        {
            perf::CounterScope scope(session);
            result.emplace(qmomkit::compute_cqmom(workload.lean_case));
        }
        if (!result->ok())
        {
            state.SkipWithError(result->message.empty() ? "Explicit CQMOM reconstruction control failed." : result->message.c_str());
            return;
        }
        benchmark::DoNotOptimize(result->quadrature.weights().data());
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
    }

    const perf::PerformanceCounters counters = session.snapshot();
    publish_geometry_counters(state, nodes_2x2(), workload);
    publish_workload_counters(state, counters);
    publish_normalized_counters(state, counters, 1);
    BenchmarkSummary summary =
        make_summary(name, "CQMOM standalone explicit reconstruction control", nodes_2x2(), workload, counters, measured_seconds,
                     static_cast<std::size_t>(state.iterations()), 1, "reconstruction", CounterMode::enabled_session);
    summary.solver_policy = "explicit";
    record_summary(std::move(summary));
}

void bm_cqmom_adaptive_reconstruction_control(benchmark::State& state, const std::string name)
{
    if (!instrumentation_is_available(state))
    {
        return;
    }

    const CqmomWorkload workload = make_phase_space_evolution_workload(nodes_2x2(), true);
    perf::CounterSession session;
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        const auto start = std::chrono::steady_clock::now();
        std::optional<qmomkit::CqmomResult> result;
        {
            perf::CounterScope scope(session);
            result.emplace(qmomkit::compute_cqmom(workload.case_));
        }
        if (!result->ok())
        {
            state.SkipWithError(result->message.empty() ? "Adaptive CQMOM reconstruction control failed." : result->message.c_str());
            return;
        }
        benchmark::DoNotOptimize(result->quadrature.weights().data());
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
    }

    const perf::PerformanceCounters counters = session.snapshot();
    publish_geometry_counters(state, nodes_2x2(), workload);
    publish_workload_counters(state, counters);
    publish_normalized_counters(state, counters, 1);
    BenchmarkSummary summary = make_summary(name, "CQMOM adaptive normal reconstruction control", nodes_2x2(), workload, counters, measured_seconds,
                                            static_cast<std::size_t>(state.iterations()), 1, "reconstruction", CounterMode::enabled_session);
    record_summary(std::move(summary));
}

void bm_cqmom_source_integration_workload(benchmark::State& state, const CounterMode counter_mode, const std::vector<std::size_t> node_counts,
                                          const std::string name, const std::string workload_kind)
{
    if (counter_mode == CounterMode::enabled_session && !instrumentation_is_available(state))
    {
        return;
    }

    const std::size_t state_count = static_cast<std::size_t>(state.range(0));
    const std::size_t steps_per_state = static_cast<std::size_t>(state.range(1));
    const CqmomWorkload workload = make_workload(node_counts);
    const qmomkit::CqmomMomentSource source = source_model();
    qmomkit::MomentSourceUpdateOptions options;
    options.ssp_rk3_options.estimate_error = true;

    perf::CounterSession session;
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        state.PauseTiming();
        std::vector<qmomkit::CqmomMomentData> moments;
        moments.reserve(state_count);
        for (std::size_t state_index = 0; state_index < state_count; ++state_index)
        {
            moments.push_back(workload.lean_moments);
        }
        state.ResumeTiming();

        const auto run_updates = [&]() -> bool
        {
            for (std::size_t step = 0; step < steps_per_state; ++step)
            {
                const double time = 0.001 * static_cast<double>(step);
                for (qmomkit::CqmomMomentData& moment_data : moments)
                {
                    qmomkit::CqmomMomentSourceUpdateResult result =
                        qmomkit::advance_cqmom_moments_with_source(moment_data, time, 0.001, source, options);
                    if (!result.ok())
                    {
                        state.SkipWithError(result.message.empty() ? "CQMOM moment-source integration workload failed." : result.message.c_str());
                        return false;
                    }
                    moment_data = std::move(result.moments);
                    benchmark::DoNotOptimize(moment_data.values.data());
                }
            }
            return true;
        };

        const auto start = std::chrono::steady_clock::now();
        bool ok = true;
        if (counter_mode == CounterMode::enabled_session)
        {
            perf::CounterScope scope(session);
            ok = run_updates();
        }
        else
        {
            ok = run_updates();
        }
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
        if (!ok)
        {
            return;
        }
    }

    const std::size_t source_steps = state_count * steps_per_state;
    const perf::PerformanceCounters counters = session.snapshot();
    publish_counter(state, "source_steps", source_steps * state.iterations());
    publish_geometry_counters(state, node_counts, workload);
    publish_workload_counters(state, counters);
    publish_normalized_counters(state, counters, source_steps);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations() * source_steps));
    BenchmarkSummary summary = make_summary(name, workload_kind, node_counts, workload, counters, measured_seconds,
                                            static_cast<std::size_t>(state.iterations()), source_steps, "moment-source update", counter_mode);
    summary.state_count = state_count;
    summary.steps_per_state = steps_per_state;
    record_summary(std::move(summary));
}

void bm_cqmom_trusted_moment_source_state_workload(benchmark::State& state, const std::vector<std::size_t> node_counts, const std::string name)
{
    if (!instrumentation_is_available(state))
    {
        return;
    }

    const std::size_t state_count = static_cast<std::size_t>(state.range(0));
    const std::size_t steps_per_state = static_cast<std::size_t>(state.range(1));
    const CqmomWorkload workload = make_workload(node_counts);
    const qmomkit::CqmomMomentSource source = source_model();
    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.source_options.ssp_rk3_options.estimate_error = true;

    perf::CounterSession session;
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        state.PauseTiming();
        std::vector<qmomkit::CqmomState> states = make_states(workload.case_, state_count);
        state.ResumeTiming();

        const auto run_updates = [&]() -> bool
        {
            for (std::size_t step = 0; step < steps_per_state; ++step)
            {
                const double time = 0.001 * static_cast<double>(step);
                for (qmomkit::CqmomState& cqmom_state : states)
                {
                    qmomkit::CqmomPhaseSpaceEvolutionResult result = qmomkit::advance_cqmom_moment_source(cqmom_state, time, 0.001, source, options);
                    if (!result.ok())
                    {
                        state.SkipWithError(result.message.empty() ? "Trusted CQMOM moment-source state workload failed." : result.message.c_str());
                        return false;
                    }
                    cqmom_state = std::move(result.state);
                    benchmark::DoNotOptimize(cqmom_state.hierarchy().number_of_dimensions());
                }
            }
            return true;
        };

        const auto start = std::chrono::steady_clock::now();
        bool ok = false;
        {
            perf::CounterScope scope(session);
            ok = run_updates();
        }
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
        if (!ok)
        {
            return;
        }
    }

    const std::size_t source_steps = state_count * steps_per_state;
    const perf::PerformanceCounters counters = session.snapshot();
    publish_counter(state, "source_steps", source_steps * state.iterations());
    publish_geometry_counters(state, node_counts, workload);
    publish_workload_counters(state, counters);
    publish_normalized_counters(state, counters, source_steps);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations() * source_steps));
    BenchmarkSummary summary =
        make_summary(name, "CQMOM trusted moment-source state phase-space evolution", node_counts, workload, counters, measured_seconds,
                     static_cast<std::size_t>(state.iterations()), source_steps, "source state-step", CounterMode::enabled_session);
    summary.state_count = state_count;
    summary.steps_per_state = steps_per_state;
    record_summary(std::move(summary));
}

void bm_cqmom_source_callback_only(benchmark::State& state)
{
    const std::size_t state_count = static_cast<std::size_t>(state.range(0));
    const std::size_t steps_per_state = static_cast<std::size_t>(state.range(1));
    const CqmomWorkload workload = make_workload({3, 3});
    const qmomkit::CqmomMomentSource source = source_model();
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t step = 0; step < steps_per_state; ++step)
        {
            const double time = 0.001 * static_cast<double>(step);
            const qmomkit::CqmomSourceContext context(time, workload.hierarchy, workload.lean_moments.layout.required_indices,
                                                      workload.lean_moments.values);
            for (std::size_t state_index = 0; state_index < state_count; ++state_index)
            {
                qmomkit::CqmomMomentRate rate = source(time, context);
                benchmark::DoNotOptimize(rate.values.data());
            }
        }
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
    }
    const std::size_t callback_evaluations = state_count * steps_per_state;
    publish_counter(state, "callback_evals", callback_evaluations * state.iterations());
    publish_geometry_counters(state, {3, 3}, workload);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations() * callback_evaluations));
    BenchmarkSummary summary =
        make_summary("performance_source_callback_only_2d_3x3", "CQMOM source callback only", {3, 3}, workload, {}, measured_seconds,
                     static_cast<std::size_t>(state.iterations()), callback_evaluations, "callback", CounterMode::timing_only);
    summary.state_count = state_count;
    summary.steps_per_state = steps_per_state;
    record_summary(std::move(summary));
}

void bm_cqmom_moment_data_copy(benchmark::State& state)
{
    const std::size_t copies_per_iteration = static_cast<std::size_t>(state.range(0));
    const CqmomWorkload workload = make_workload({3, 3});
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t copy = 0; copy < copies_per_iteration; ++copy)
        {
            qmomkit::CqmomMomentData moment_data = workload.lean_moments;
            benchmark::DoNotOptimize(moment_data.values.data());
        }
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
    }
    publish_counter(state, "copies", copies_per_iteration * state.iterations());
    publish_geometry_counters(state, {3, 3}, workload);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations() * copies_per_iteration));
    record_summary(make_summary("performance_cqmom_moment_data_copy_2d_3x3", "CqmomMomentData copy/layout preparation", {3, 3}, workload, {},
                                measured_seconds, static_cast<std::size_t>(state.iterations()), copies_per_iteration, "copy",
                                CounterMode::timing_only));
}

void bm_cqmom_ssp_source_existing_layout(benchmark::State& state)
{
    if (!instrumentation_is_available(state))
    {
        return;
    }
    const std::size_t state_count = static_cast<std::size_t>(state.range(0));
    const std::size_t steps_per_state = static_cast<std::size_t>(state.range(1));
    const CqmomWorkload workload = make_workload({3, 3});
    const qmomkit::CqmomMomentSource source = source_model();
    qmomkit::SspRk3Options options;
    options.estimate_error = true;
    perf::CounterSession session;
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        state.PauseTiming();
        std::vector<std::vector<double>> values(state_count, workload.lean_moments.values);
        state.ResumeTiming();
        const auto start = std::chrono::steady_clock::now();
        perf::CounterScope scope(session);
        for (std::size_t step = 0; step < steps_per_state; ++step)
        {
            const double time = 0.001 * static_cast<double>(step);
            for (std::vector<double>& value_state : values)
            {
                qmomkit::SspRk3Result result = qmomkit::ssp_rk3_step(
                    value_state, time, 0.001,
                    [&](const double rhs_time, const std::vector<double>& state_values, std::vector<double>& derivative)
                    {
                        const qmomkit::CqmomSourceContext context(rhs_time, workload.hierarchy, workload.lean_moments.layout.required_indices,
                                                                  state_values);
                        const qmomkit::CqmomMomentRate rate = source(rhs_time, context);
                        derivative.assign(rate.values.begin(), rate.values.end());
                    },
                    options);
                if (!result.ok())
                {
                    state.SkipWithError("SSP-RK3 source existing-layout workload failed.");
                    return;
                }
                value_state = std::move(result.state);
                benchmark::DoNotOptimize(value_state.data());
            }
        }
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
    }
    const std::size_t source_steps = state_count * steps_per_state;
    const perf::PerformanceCounters counters = session.snapshot();
    publish_counter(state, "source_steps", source_steps * state.iterations());
    publish_geometry_counters(state, {3, 3}, workload);
    publish_workload_counters(state, counters);
    publish_normalized_counters(state, counters, source_steps);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations() * source_steps));
    BenchmarkSummary summary = make_summary("performance_ssp_source_existing_layout_2d_3x3", "SSP-RK3 source with existing layout", {3, 3}, workload,
                                            counters, measured_seconds, static_cast<std::size_t>(state.iterations()), source_steps,
                                            "moment-source update", CounterMode::enabled_session);
    summary.state_count = state_count;
    summary.steps_per_state = steps_per_state;
    record_summary(std::move(summary));
}

void bm_qmom_node_advection_scheme(benchmark::State& state, const qmomkit::NodeAdvectionScheme scheme, const NodeTrajectory trajectory,
                                   const std::string name)
{
    if (!instrumentation_is_available(state))
    {
        return;
    }

    constexpr std::size_t node_count = 4;
    const double time_step = trajectory == NodeTrajectory::constant ? 0.25 : 0.75;
    const qmomkit::Quadrature1D initial = make_qmom_node_advection_fixture(node_count);
    const qmomkit::PhaseSpaceVelocity velocity = node_advection_velocity_model(trajectory);
    const qmomkit::NodeAdvectionOptions options = node_advection_benchmark_options(scheme, time_step);

    perf::CounterSession session;
    double measured_seconds = 0.0;
    std::size_t diagnostic_steps = 0;
    std::size_t diagnostic_rejected_steps = 0;
    std::size_t diagnostic_rhs_evaluations = 0;
    double max_trajectory_error = 0.0;
    for (auto _ : state)
    {
        const auto start = std::chrono::steady_clock::now();
        qmomkit::QmomNodeAdvectionResult result;
        {
            perf::CounterScope scope(session);
            result = qmomkit::advance_nodes(initial, 0.0, time_step, velocity, options);
        }
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
        if (!result.ok())
        {
            state.SkipWithError(result.message.empty() ? "QMOM node advection failed." : result.message.c_str());
            return;
        }
        diagnostic_steps += result.diagnostics.step_count;
        diagnostic_rejected_steps += result.diagnostics.rejected_step_count;
        diagnostic_rhs_evaluations += result.diagnostics.rhs_evaluations;
        max_trajectory_error = std::max(max_trajectory_error, max_qmom_trajectory_error(initial, result.quadrature, trajectory, time_step));
        benchmark::DoNotOptimize(result.quadrature.abscissas.data());
    }

    const perf::PerformanceCounters counters = session.snapshot();
    state.counters["dimensions"] = 1.0;
    state.counters["requested_terminal_nodes"] = static_cast<double>(node_count);
    publish_counter(state, "node_updates", static_cast<std::size_t>(state.iterations()));
    publish_workload_counters(state, counters);
    publish_node_advection_diagnostics(state, diagnostic_steps, diagnostic_rejected_steps, diagnostic_rhs_evaluations, max_trajectory_error);
    publish_normalized_counters(state, counters, 1);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));

    BenchmarkSummary summary;
    summary.name = name;
    summary.workload_kind = "QMOM node advection, " + node_trajectory_name(trajectory) + " trajectory";
    summary.node_counts = {node_count};
    summary.node_advection_scheme = node_advection_scheme_name(scheme);
    summary.instrumentation_mode = instrumentation_mode_name(CounterMode::enabled_session);
    summary.operations_per_iteration = 1;
    summary.operation_label = "node update";
    summary.iterations = static_cast<std::size_t>(state.iterations());
    summary.total_runtime_seconds = measured_seconds;
    summary.requested_terminal_nodes = node_count;
    summary.diagnostic_node_advection_steps = diagnostic_steps;
    summary.diagnostic_node_advection_rejected_steps = diagnostic_rejected_steps;
    summary.diagnostic_node_advection_rhs_evaluations = diagnostic_rhs_evaluations;
    summary.has_trajectory_error = true;
    summary.max_trajectory_error = max_trajectory_error;
    summary.counters = counters;
    record_summary(std::move(summary));
}

void bm_cqmom_node_advection_scheme(benchmark::State& state, const qmomkit::NodeAdvectionScheme scheme, const NodeTrajectory trajectory,
                                    const std::string name)
{
    if (!instrumentation_is_available(state))
    {
        return;
    }
    const CqmomWorkload workload = make_workload({2, 2});
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(workload.lean_case);
    if (!initial.ok())
    {
        state.SkipWithError("CQMOM setup inversion failed.");
        return;
    }
    const double time_step = trajectory == NodeTrajectory::constant ? 0.25 : 0.75;
    const qmomkit::PhaseSpaceVelocity velocity = node_advection_velocity_model(trajectory);
    const qmomkit::NodeAdvectionOptions options = node_advection_benchmark_options(scheme, time_step);
    perf::CounterSession session;
    double measured_seconds = 0.0;
    std::size_t diagnostic_steps = 0;
    std::size_t diagnostic_rejected_steps = 0;
    std::size_t diagnostic_rhs_evaluations = 0;
    double max_trajectory_error = 0.0;
    for (auto _ : state)
    {
        const auto start = std::chrono::steady_clock::now();
        std::optional<qmomkit::CqmomNodeAdvectionResult> result;
        {
            perf::CounterScope scope(session);
            result.emplace(qmomkit::advance_nodes(initial.quadrature, 0.0, time_step, velocity, options));
        }
        if (!result->ok())
        {
            state.SkipWithError(result->message.empty() ? "CQMOM node advection failed." : result->message.c_str());
            return;
        }
        diagnostic_steps += result->diagnostics.step_count;
        diagnostic_rejected_steps += result->diagnostics.rejected_step_count;
        diagnostic_rhs_evaluations += result->diagnostics.rhs_evaluations;
        max_trajectory_error =
            std::max(max_trajectory_error, max_cqmom_trajectory_error(initial.quadrature, result->quadrature, trajectory, time_step));
        benchmark::DoNotOptimize(result->quadrature.abscissas().data());
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
    }
    const perf::PerformanceCounters counters = session.snapshot();
    publish_geometry_counters(state, {2, 2}, workload);
    publish_workload_counters(state, counters);
    publish_node_advection_diagnostics(state, diagnostic_steps, diagnostic_rejected_steps, diagnostic_rhs_evaluations, max_trajectory_error);
    publish_normalized_counters(state, counters, 1);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()));
    BenchmarkSummary summary =
        make_summary(name, "CQMOM node advection, " + node_trajectory_name(trajectory) + " trajectory", {2, 2}, workload, counters, measured_seconds,
                     static_cast<std::size_t>(state.iterations()), 1, "advection update", CounterMode::enabled_session);
    summary.node_advection_scheme = node_advection_scheme_name(scheme);
    summary.reconstructed_moments = workload.lean_moments.values.size();
    summary.diagnostic_node_advection_steps = diagnostic_steps;
    summary.diagnostic_node_advection_rejected_steps = diagnostic_rejected_steps;
    summary.diagnostic_node_advection_rhs_evaluations = diagnostic_rhs_evaluations;
    summary.has_trajectory_error = true;
    summary.max_trajectory_error = max_trajectory_error;
    record_summary(std::move(summary));
}

void bm_cqmom_moment_reconstruction_only(benchmark::State& state)
{
    const CqmomWorkload workload = make_workload({2, 2});
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(workload.lean_case);
    if (!initial.ok())
    {
        state.SkipWithError("CQMOM setup inversion failed.");
        return;
    }
    double measured_seconds = 0.0;
    for (auto _ : state)
    {
        const auto start = std::chrono::steady_clock::now();
        qmomkit::CqmomMomentData reconstructed = qmomkit::reconstruct_cqmom_moments(initial.quadrature, workload.lean_moments.layout);
        benchmark::DoNotOptimize(reconstructed.values.data());
        const auto stop = std::chrono::steady_clock::now();
        measured_seconds += std::chrono::duration<double>(stop - start).count();
    }
    publish_geometry_counters(state, {2, 2}, workload);
    state.counters["reconstructed_moments"] = static_cast<double>(workload.lean_moments.values.size());
    BenchmarkSummary summary =
        make_summary("performance_phase_space_evolution_component_moment_reconstruction_2d_2x2", "CQMOM raw-moment reconstruction", {2, 2}, workload,
                     {}, measured_seconds, static_cast<std::size_t>(state.iterations()), 1, "moment reconstruction", CounterMode::timing_only);
    summary.reconstructed_moments = workload.lean_moments.values.size();
    record_summary(std::move(summary));
}

void bm_cqmom_inversion_only(benchmark::State& state)
{
    bm_cqmom_reconstruction_scaling(state, {2, 2}, qmomkit::VandermondeSolver::automatic, qmomkit::CqmomDiagnosticLevel::basic,
                                    CounterMode::enabled_session, "performance_phase_space_evolution_component_inversion_2d_2x2");
}

BENCHMARK_CAPTURE(bm_cqmom_phase_space_evolution_workload, strang_2d_2x2, qmomkit::PhaseSpaceEvolutionScheme::strang, CounterMode::enabled_session,
                  "performance_cqmom_phase_space_evolution_strang_2d_2x2")
    ->Name("performance_cqmom_phase_space_evolution_strang_2d_2x2")
    ->Args({32, 4});
BENCHMARK_CAPTURE(bm_cqmom_phase_space_evolution_workload, advsrc_2d_2x2, qmomkit::PhaseSpaceEvolutionScheme::advsrc, CounterMode::enabled_session,
                  "performance_cqmom_phase_space_evolution_advsrc_2d_2x2")
    ->Name("performance_cqmom_phase_space_evolution_advsrc_2d_2x2")
    ->Args({32, 4});
BENCHMARK_CAPTURE(bm_cqmom_phase_space_evolution_workload, srcadv_2d_2x2, qmomkit::PhaseSpaceEvolutionScheme::srcadv, CounterMode::enabled_session,
                  "performance_cqmom_phase_space_evolution_srcadv_2d_2x2")
    ->Name("performance_cqmom_phase_space_evolution_srcadv_2d_2x2")
    ->Args({32, 4});

BENCHMARK_CAPTURE(bm_cqmom_phase_space_evolution_workload, overhead_no_session_2d_2x2, qmomkit::PhaseSpaceEvolutionScheme::strang,
                  CounterMode::enabled_no_session, "performance_instrumentation_overhead_phase_space_evolution_no_session_2d_2x2")
    ->Name("performance_instrumentation_overhead_phase_space_evolution_no_session_2d_2x2")
    ->Args({32, 4});

const bool ownership_phase_space_evolution_benchmarks_registered = []
{
    struct Fixture
    {
        const char* name;
        std::vector<std::size_t> node_counts;
        bool adaptive;
        bool cache_miss;
    };
    const std::vector<Fixture> fixtures{{"fixed_2d_2x2", nodes_2x2(), false, false},
                                        {"fixed_2d_3x3", nodes_3x3(), false, false},
                                        {"fixed_3d_2x2x2", nodes_2x2x2(), false, false},
                                        {"adaptive_ragged_2d", nodes_2x2(), true, false},
                                        {"adaptive_ragged_3d", nodes_2x2x2(), true, false},
                                        {"adaptive_ragged_2d_cache_miss", nodes_2x2(), true, true},
                                        {"adaptive_ragged_3d_cache_miss", nodes_2x2x2(), true, true}};
    // Isolated operators plus one composite. The first-order splitting orders
    // (advsrc/srcadv) are covered by the dedicated
    // performance_cqmom_phase_space_evolution_{advsrc,srcadv} benchmarks; repeating
    // both orders across every ownership fixture would re-ask an answered question.
    const std::vector<CqmomPhaseSpaceEvolutionWorkflow> workflows{
        CqmomPhaseSpaceEvolutionWorkflow::advection, CqmomPhaseSpaceEvolutionWorkflow::moment_source, CqmomPhaseSpaceEvolutionWorkflow::strang};
    for (const Fixture& fixture : fixtures)
    {
        for (const CqmomPhaseSpaceEvolutionWorkflow workflow : workflows)
        {
            const std::string name = "performance_cqmom_phase_space_evolution_ownership_" + std::string(fixture.name) + "_" +
                                     phase_space_evolution_workflow_name(workflow);
            benchmark::RegisterBenchmark(name.c_str(),
                                         [fixture, workflow, name](benchmark::State& state)
                                         {
                                             bm_cqmom_phase_space_evolution_ownership_workload(state, workflow, fixture.node_counts, fixture.adaptive,
                                                                                               fixture.cache_miss, name);
                                         })
                ->Args({32, fixture.cache_miss ? 1 : 4});
        }
    }
    return true;
}();

BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, auto_basic_2d_2x2, nodes_2x2(), qmomkit::VandermondeSolver::automatic,
                  qmomkit::CqmomDiagnosticLevel::basic, CounterMode::enabled_session, "performance_vandermonde_auto_basic_2d_2x2")
    ->Name("performance_vandermonde_auto_basic_2d_2x2");
BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, bp_basic_2d_2x2, nodes_2x2(), qmomkit::VandermondeSolver::bjorck_pereyra,
                  qmomkit::CqmomDiagnosticLevel::basic, CounterMode::enabled_session, "performance_vandermonde_bp_basic_2d_2x2")
    ->Name("performance_vandermonde_bp_basic_2d_2x2");
BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, dense_basic_2d_2x2, nodes_2x2(), qmomkit::VandermondeSolver::dense,
                  qmomkit::CqmomDiagnosticLevel::basic, CounterMode::enabled_session, "performance_vandermonde_dense_basic_2d_2x2")
    ->Name("performance_vandermonde_dense_basic_2d_2x2");
BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, auto_detailed_2d_2x2, nodes_2x2(), qmomkit::VandermondeSolver::automatic,
                  qmomkit::CqmomDiagnosticLevel::detailed, CounterMode::enabled_session, "performance_vandermonde_auto_detailed_2d_2x2")
    ->Name("performance_vandermonde_auto_detailed_2d_2x2");

BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, auto_basic_2d_4x4, nodes_4x4(), qmomkit::VandermondeSolver::automatic,
                  qmomkit::CqmomDiagnosticLevel::basic, CounterMode::enabled_session, "performance_vandermonde_auto_basic_2d_4x4")
    ->Name("performance_vandermonde_auto_basic_2d_4x4");
BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, bp_basic_2d_4x4, nodes_4x4(), qmomkit::VandermondeSolver::bjorck_pereyra,
                  qmomkit::CqmomDiagnosticLevel::basic, CounterMode::enabled_session, "performance_vandermonde_bp_basic_2d_4x4")
    ->Name("performance_vandermonde_bp_basic_2d_4x4");
BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, dense_basic_2d_4x4, nodes_4x4(), qmomkit::VandermondeSolver::dense,
                  qmomkit::CqmomDiagnosticLevel::basic, CounterMode::enabled_session, "performance_vandermonde_dense_basic_2d_4x4")
    ->Name("performance_vandermonde_dense_basic_2d_4x4");
BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, auto_detailed_2d_4x4, nodes_4x4(), qmomkit::VandermondeSolver::automatic,
                  qmomkit::CqmomDiagnosticLevel::detailed, CounterMode::enabled_session, "performance_vandermonde_auto_detailed_2d_4x4")
    ->Name("performance_vandermonde_auto_detailed_2d_4x4");

BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, auto_basic_3d_2x2x2, nodes_2x2x2(), qmomkit::VandermondeSolver::automatic,
                  qmomkit::CqmomDiagnosticLevel::basic, CounterMode::enabled_session, "performance_vandermonde_auto_basic_3d_2x2x2")
    ->Name("performance_vandermonde_auto_basic_3d_2x2x2");
BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, bp_basic_3d_2x2x2, nodes_2x2x2(), qmomkit::VandermondeSolver::bjorck_pereyra,
                  qmomkit::CqmomDiagnosticLevel::basic, CounterMode::enabled_session, "performance_vandermonde_bp_basic_3d_2x2x2")
    ->Name("performance_vandermonde_bp_basic_3d_2x2x2");
BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, dense_basic_3d_2x2x2, nodes_2x2x2(), qmomkit::VandermondeSolver::dense,
                  qmomkit::CqmomDiagnosticLevel::basic, CounterMode::enabled_session, "performance_vandermonde_dense_basic_3d_2x2x2")
    ->Name("performance_vandermonde_dense_basic_3d_2x2x2");
BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, auto_detailed_3d_2x2x2, nodes_2x2x2(), qmomkit::VandermondeSolver::automatic,
                  qmomkit::CqmomDiagnosticLevel::detailed, CounterMode::enabled_session, "performance_vandermonde_auto_detailed_3d_2x2x2")
    ->Name("performance_vandermonde_auto_detailed_3d_2x2x2");

BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, auto_basic_2d_3x3, nodes_3x3(), qmomkit::VandermondeSolver::automatic,
                  qmomkit::CqmomDiagnosticLevel::basic, CounterMode::enabled_session, "performance_vandermonde_auto_basic_2d_3x3")
    ->Name("performance_vandermonde_auto_basic_2d_3x3");
BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, auto_basic_4d_2x2x2x2, nodes_2x2x2x2(), qmomkit::VandermondeSolver::automatic,
                  qmomkit::CqmomDiagnosticLevel::basic, CounterMode::enabled_session, "performance_vandermonde_auto_basic_4d_2x2x2x2")
    ->Name("performance_vandermonde_auto_basic_4d_2x2x2x2");
BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, auto_basic_5d_2x2x2x2x2, nodes_2x2x2x2x2(), qmomkit::VandermondeSolver::automatic,
                  qmomkit::CqmomDiagnosticLevel::basic, CounterMode::enabled_session, "performance_vandermonde_auto_basic_5d_2x2x2x2x2")
    ->Name("performance_vandermonde_auto_basic_5d_2x2x2x2x2");

BENCHMARK_CAPTURE(bm_cqmom_normal_prepared_reconstruction_scaling, normal_prepared_2d_2x2, nodes_2x2(), "performance_cqmom_normal_prepared_2d_2x2")
    ->Name("performance_cqmom_normal_prepared_2d_2x2");
BENCHMARK_CAPTURE(bm_cqmom_normal_prepared_reconstruction_scaling, normal_prepared_2d_3x3, nodes_3x3(), "performance_cqmom_normal_prepared_2d_3x3")
    ->Name("performance_cqmom_normal_prepared_2d_3x3");
BENCHMARK_CAPTURE(bm_cqmom_normal_prepared_reconstruction_scaling, normal_prepared_3d_2x2x2, nodes_2x2x2(),
                  "performance_cqmom_normal_prepared_3d_2x2x2")
    ->Name("performance_cqmom_normal_prepared_3d_2x2x2");
BENCHMARK_CAPTURE(bm_cqmom_normal_prepared_reconstruction_scaling, normal_prepared_4d_2x2x2x2, nodes_2x2x2x2(),
                  "performance_cqmom_normal_prepared_4d_2x2x2x2")
    ->Name("performance_cqmom_normal_prepared_4d_2x2x2x2");
BENCHMARK_CAPTURE(bm_cqmom_normal_prepared_reconstruction_scaling, normal_prepared_5d_2x2x2x2x2, nodes_2x2x2x2x2(),
                  "performance_cqmom_normal_prepared_5d_2x2x2x2x2")
    ->Name("performance_cqmom_normal_prepared_5d_2x2x2x2x2");

BENCHMARK_CAPTURE(bm_cqmom_reconstruction_scaling, overhead_no_session_2d_2x2, nodes_2x2(), qmomkit::VandermondeSolver::automatic,
                  qmomkit::CqmomDiagnosticLevel::basic, CounterMode::enabled_no_session,
                  "performance_instrumentation_overhead_reconstruction_no_session_2d_2x2")
    ->Name("performance_instrumentation_overhead_reconstruction_no_session_2d_2x2");

BENCHMARK_CAPTURE(bm_cqmom_fixed_case_preparation, cold_2d_2x2, nodes_2x2(), "performance_cqmom_programmatic_case_construction_2d_2x2")
    ->Name("performance_cqmom_programmatic_case_construction_2d_2x2");
BENCHMARK_CAPTURE(bm_cqmom_fixed_case_preparation, cold_2d_3x3, nodes_3x3(), "performance_cqmom_programmatic_case_construction_2d_3x3")
    ->Name("performance_cqmom_programmatic_case_construction_2d_3x3");
BENCHMARK_CAPTURE(bm_cqmom_fixed_case_preparation, cold_3d_2x2x2, nodes_2x2x2(), "performance_cqmom_programmatic_case_construction_3d_2x2x2")
    ->Name("performance_cqmom_programmatic_case_construction_3d_2x2x2");
BENCHMARK_CAPTURE(bm_cqmom_fixed_case_preparation, cold_4d_2x2x2x2, nodes_2x2x2x2(), "performance_cqmom_programmatic_case_construction_4d_2x2x2x2")
    ->Name("performance_cqmom_programmatic_case_construction_4d_2x2x2x2");
BENCHMARK_CAPTURE(bm_cqmom_fixed_case_preparation, cold_5d_2x2x2x2x2, nodes_2x2x2x2x2(),
                  "performance_cqmom_programmatic_case_construction_5d_2x2x2x2x2")
    ->Name("performance_cqmom_programmatic_case_construction_5d_2x2x2x2x2");

BENCHMARK_CAPTURE(bm_cqmom_initial_state_prepared, fixed_2d_2x2, nodes_2x2(), "performance_cqmom_initial_state_prepared_2d_2x2")
    ->Name("performance_cqmom_initial_state_prepared_2d_2x2");

BENCHMARK_CAPTURE(bm_cqmom_explicit_reconstruction_control, lean_2d_2x2, "performance_cqmom_reconstruction_explicit_lean_2d_2x2")
    ->Name("performance_cqmom_reconstruction_explicit_lean_2d_2x2");
BENCHMARK_CAPTURE(bm_cqmom_adaptive_reconstruction_control, normal_lean_2d_2x2, "performance_cqmom_reconstruction_adaptive_normal_lean_2d_2x2")
    ->Name("performance_cqmom_reconstruction_adaptive_normal_lean_2d_2x2");

BENCHMARK_CAPTURE(bm_cqmom_source_integration_workload, source_normal_2d_3x3, CounterMode::enabled_session, nodes_3x3(),
                  "performance_source_normal_2d_3x3", "CQMOM moment-source integration")
    ->Name("performance_source_normal_2d_3x3")
    ->Args({64, 8});
BENCHMARK_CAPTURE(bm_cqmom_source_integration_workload, source_no_session_2d_3x3, CounterMode::enabled_no_session, nodes_3x3(),
                  "performance_instrumentation_overhead_source_no_session_2d_3x3", "CQMOM moment-source integration")
    ->Name("performance_instrumentation_overhead_source_no_session_2d_3x3")
    ->Args({64, 8});
BENCHMARK_CAPTURE(bm_cqmom_trusted_moment_source_state_workload, trusted_state_2d_2x2, nodes_2x2(), "performance_source_trusted_state_2d_2x2")
    ->Name("performance_source_trusted_state_2d_2x2")
    ->Args({32, 4});
BENCHMARK_CAPTURE(bm_cqmom_trusted_moment_source_state_workload, trusted_state_2d_3x3, nodes_3x3(), "performance_source_trusted_state_2d_3x3")
    ->Name("performance_source_trusted_state_2d_3x3")
    ->Args({64, 8});
BENCHMARK_CAPTURE(bm_cqmom_trusted_moment_source_state_workload, trusted_state_3d_2x2x2, nodes_2x2x2(), "performance_source_trusted_state_3d_2x2x2")
    ->Name("performance_source_trusted_state_3d_2x2x2")
    ->Args({16, 2});
BENCHMARK(bm_cqmom_source_callback_only)->Name("performance_source_callback_only_2d_3x3")->Args({64, 8});
BENCHMARK(bm_cqmom_moment_data_copy)->Name("performance_cqmom_moment_data_copy_2d_3x3")->Arg(512);
BENCHMARK(bm_cqmom_ssp_source_existing_layout)->Name("performance_ssp_source_existing_layout_2d_3x3")->Args({64, 8});

BENCHMARK_CAPTURE(bm_qmom_node_advection_scheme, constant_ssp_rk3_1d_4, qmomkit::NodeAdvectionScheme::ssp_rk3, NodeTrajectory::constant,
                  "performance_node_advection_scheme_qmom_constant_1d_4_ssp_rk3")
    ->Name("performance_node_advection_scheme_qmom_constant_1d_4_ssp_rk3");
BENCHMARK_CAPTURE(bm_qmom_node_advection_scheme, constant_rk23_1d_4, qmomkit::NodeAdvectionScheme::rk23, NodeTrajectory::constant,
                  "performance_node_advection_scheme_qmom_constant_1d_4_rk23")
    ->Name("performance_node_advection_scheme_qmom_constant_1d_4_rk23");
BENCHMARK_CAPTURE(bm_qmom_node_advection_scheme, affine_ssp_rk3_1d_4, qmomkit::NodeAdvectionScheme::ssp_rk3, NodeTrajectory::affine,
                  "performance_node_advection_scheme_qmom_affine_1d_4_ssp_rk3")
    ->Name("performance_node_advection_scheme_qmom_affine_1d_4_ssp_rk3");
BENCHMARK_CAPTURE(bm_qmom_node_advection_scheme, affine_rk23_1d_4, qmomkit::NodeAdvectionScheme::rk23, NodeTrajectory::affine,
                  "performance_node_advection_scheme_qmom_affine_1d_4_rk23")
    ->Name("performance_node_advection_scheme_qmom_affine_1d_4_rk23");

BENCHMARK_CAPTURE(bm_cqmom_node_advection_scheme, constant_ssp_rk3_2d_2x2, qmomkit::NodeAdvectionScheme::ssp_rk3, NodeTrajectory::constant,
                  "performance_node_advection_scheme_cqmom_constant_2d_2x2_ssp_rk3")
    ->Name("performance_node_advection_scheme_cqmom_constant_2d_2x2_ssp_rk3");
BENCHMARK_CAPTURE(bm_cqmom_node_advection_scheme, constant_rk23_2d_2x2, qmomkit::NodeAdvectionScheme::rk23, NodeTrajectory::constant,
                  "performance_node_advection_scheme_cqmom_constant_2d_2x2_rk23")
    ->Name("performance_node_advection_scheme_cqmom_constant_2d_2x2_rk23");
BENCHMARK_CAPTURE(bm_cqmom_node_advection_scheme, affine_ssp_rk3_2d_2x2, qmomkit::NodeAdvectionScheme::ssp_rk3, NodeTrajectory::affine,
                  "performance_node_advection_scheme_cqmom_affine_2d_2x2_ssp_rk3")
    ->Name("performance_node_advection_scheme_cqmom_affine_2d_2x2_ssp_rk3");
BENCHMARK_CAPTURE(bm_cqmom_node_advection_scheme, affine_rk23_2d_2x2, qmomkit::NodeAdvectionScheme::rk23, NodeTrajectory::affine,
                  "performance_node_advection_scheme_cqmom_affine_2d_2x2_rk23")
    ->Name("performance_node_advection_scheme_cqmom_affine_2d_2x2_rk23");
BENCHMARK_CAPTURE(bm_cqmom_source_integration_workload, phase_space_evolution_component_source_2d_2x2, CounterMode::enabled_session, nodes_2x2(),
                  "performance_phase_space_evolution_component_source_update_2d_2x2", "CQMOM moment-source update")
    ->Name("performance_phase_space_evolution_component_source_update_2d_2x2")
    ->Args({32, 4});
BENCHMARK(bm_cqmom_moment_reconstruction_only)->Name("performance_phase_space_evolution_component_moment_reconstruction_2d_2x2");
BENCHMARK(bm_cqmom_inversion_only)->Name("performance_phase_space_evolution_component_inversion_2d_2x2");

} // namespace
