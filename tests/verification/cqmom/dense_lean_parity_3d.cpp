#include "../verification_helpers.hpp"
#include "internal/cqmom_case.hpp"
#include "internal/linear_algebra.hpp"
#include "internal/recursive_quadrature.hpp"
#include "internal/vandermonde.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/moment_source.hpp"
#include "qmomkit/multi_index.hpp"
#include "qmomkit/phase_space_evolution.hpp"
#include "qmomkit/wheeler.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace
{

constexpr double value_abs_tol = 2.0e-8;
constexpr double value_rel_tol = 2.0e-10;
constexpr double quadrature_abs_tol = 2.0e-8;
constexpr double quadrature_rel_tol = 2.0e-10;

struct TerminalNode
{
    double weight = 0.0;
    std::array<double, 3> coordinate{};
    qmomkit::ConditioningPath path;
};

struct BranchSigmaDiagnostic
{
    std::vector<double> moments;
    double sigma22 = 0.0;
    double tolerance = 0.0;
};

struct Recurrence2Diagnostic
{
    double sigma11 = 0.0;
    double sigma12 = 0.0;
    double sigma13 = 0.0;
    double sigma22 = 0.0;
    double sigma23 = 0.0;
    double a0 = 0.0;
    double a1 = 0.0;
    double b0 = 0.0;
    double b1 = 0.0;
    double variance = 0.0;
};

struct BranchDirectTrace
{
    std::vector<std::size_t> flat_terminal_indices;
    std::vector<TerminalNode> terminal_nodes;
    double branch_mass = 0.0;
    std::vector<double> conditional_weights;
    std::vector<double> normalized_moments;
    std::vector<double> unnormalized_moments;
    std::vector<long double> normalized_moments_long;
    long double variance_long = 0.0L;
};

struct VandermondeStageTrace
{
    std::vector<double> parent_abscissas;
    std::vector<double> parent_weights;
    qmomkit::internal::DenseMatrix rhs{0, 0};
    qmomkit::internal::VandermondeSolveResult solve{qmomkit::internal::DenseMatrix(0, 0)};
    double residual = 0.0;
    double condition = 0.0;
};

struct RecursiveBranchTrace
{
    std::vector<double> d_moments;
    std::vector<double> d_nodes;
    std::vector<double> d_weights;
    VandermondeStageTrace d_to_chi_u;
    std::vector<qmomkit::MultiIndex> d0_child_indices;
    std::vector<double> d0_child_values;
    std::vector<double> d0_chi_moments;
    std::vector<double> d0_chi_nodes;
    std::vector<double> d0_chi_weights;
    VandermondeStageTrace chi_to_u;
    std::vector<qmomkit::MultiIndex> d0_chi1_child_indices;
    std::vector<double> normalized_moments;
};

struct FailureProbe
{
    qmomkit::RawMomentTensor dense_moments;
    std::vector<TerminalNode> terminal_nodes;
    std::vector<TerminalNode> previous_terminal_nodes;
    qmomkit::CqmomResult dense_result;
    qmomkit::CqmomResult lean_result;
    std::size_t macro_step = 0;
    double time_step = 0.0;
    double fluid_velocity_base = 0.0;
    double fluid_chi_coefficient = 0.0;
    double fluid_d_coefficient = 0.0;
};

struct DragTransformParameters
{
    double time_step = 0.0;
    double fluid_velocity_base = 0.0;
    double fluid_chi_coefficient = 0.0;
    double fluid_d_coefficient = 0.0;
    std::string label;
};

struct DragTransformTrace
{
    qmomkit::ConditioningPath path;
    double u_old = 0.0;
    double fluid_velocity = 0.0;
    double d = 0.0;
    double chi = 0.0;
    double lambda = 0.0;
    double time_step = 0.0;
    double lambda_dt = 0.0;
    double decay = 0.0;
    double u_new = 0.0;
};

BranchDirectTrace direct_trace_for_branch(const std::vector<TerminalNode>& nodes, const qmomkit::ConditioningPath& branch_path, std::size_t dimension,
                                          std::size_t moment_count);

std::string index_string(const qmomkit::MultiIndex& index)
{
    std::ostringstream stream;
    stream << '(';
    for (std::size_t i = 0; i < index.size(); ++i)
    {
        if (i > 0)
        {
            stream << ',';
        }
        stream << index[i];
    }
    stream << ')';
    return stream.str();
}

std::string path_string(const qmomkit::ConditioningPath& path)
{
    return index_string(path.node_indices_by_level());
}

double integer_power(const double base, const std::size_t exponent)
{
    double value = 1.0;
    for (std::size_t i = 0; i < exponent; ++i)
    {
        value *= base;
    }
    return value;
}

long double integer_power_long(const long double base, const std::size_t exponent)
{
    long double value = 1.0L;
    for (std::size_t i = 0; i < exponent; ++i)
    {
        value *= base;
    }
    return value;
}

qmomkit::RawMomentTensor raw_moments_from_terminals(const std::vector<TerminalNode>& nodes, const std::vector<std::size_t>& extents)
{
    std::vector<double> values(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat = 0; flat < values.size(); ++flat)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat, extents);
        double moment = 0.0;
        for (const TerminalNode& node : nodes)
        {
            double contribution = node.weight;
            for (std::size_t dimension = 0; dimension < moment_index.size(); ++dimension)
            {
                contribution *= integer_power(node.coordinate[dimension], moment_index[dimension]);
            }
            moment += contribution;
        }
        values[flat] = moment;
    }
    return qmomkit::RawMomentTensor{extents, std::move(values)};
}

qmomkit::CqmomMomentData lean_moments_from_dense(const qmomkit::RawMomentTensor& moments, const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    return qmomkit::extract_cqmom_moments(moments, layout);
}

qmomkit::CqmomResult compute_dense(const qmomkit::RawMomentTensor& moments, const qmomkit::PhaseSpaceHierarchy& hierarchy,
                                   const qmomkit::QmomOptions& qmom_options, const bool adaptive_branch_counts)
{
    qmomkit::CqmomDenseCase dense_case{hierarchy, moments};
    dense_case.qmom_options = qmom_options;
    dense_case.adaptive_branch_counts = adaptive_branch_counts;
    dense_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::detailed;
    dense_case.vandermonde_solver = qmomkit::VandermondeSolver::automatic;
    return qmomkit::compute_cqmom(dense_case);
}

qmomkit::CqmomResult compute_lean(const qmomkit::RawMomentTensor& moments, const qmomkit::PhaseSpaceHierarchy& hierarchy,
                                  const qmomkit::QmomOptions& qmom_options, const bool adaptive_branch_counts)
{
    qmomkit::CqmomLeanCase lean_case{lean_moments_from_dense(moments, hierarchy)};
    lean_case.qmom_options = qmom_options;
    lean_case.adaptive_branch_counts = adaptive_branch_counts;
    lean_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::detailed;
    lean_case.vandermonde_solver = qmomkit::VandermondeSolver::automatic;
    return qmomkit::compute_cqmom(lean_case);
}

bool path_less(const TerminalNode& lhs, const TerminalNode& rhs)
{
    return lhs.path < rhs.path;
}

std::vector<TerminalNode> terminal_nodes_from_quadrature(const qmomkit::CqmomQuadrature& quadrature)
{
    std::vector<TerminalNode> nodes;
    nodes.reserve(quadrature.weights().size());
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        TerminalNode terminal;
        terminal.weight = quadrature.weights()[node];
        for (std::size_t dimension = 0; dimension < terminal.coordinate.size(); ++dimension)
        {
            terminal.coordinate[dimension] = quadrature.abscissas()[dimension][node];
        }
        terminal.path = quadrature.terminal_paths()[node];
        nodes.push_back(terminal);
    }
    std::sort(nodes.begin(), nodes.end(), path_less);
    return nodes;
}

void compare_quadratures(const qmomkit::CqmomQuadrature& dense, const qmomkit::CqmomQuadrature& lean)
{
    QMOMKIT_VERIFY_CHECK(dense.hierarchy().node_counts() == lean.hierarchy().node_counts());
    QMOMKIT_VERIFY_CHECK(dense.hierarchy().conditioning_order() == lean.hierarchy().conditioning_order());
    QMOMKIT_VERIFY_CHECK(dense.total_nodes() == lean.total_nodes());
    QMOMKIT_VERIFY_CHECK(dense.has_rectangular_storage() == lean.has_rectangular_storage());
    QMOMKIT_VERIFY_CHECK(dense.adaptive_branch_counts_enabled() == lean.adaptive_branch_counts_enabled());

    const std::vector<TerminalNode> dense_nodes = terminal_nodes_from_quadrature(dense);
    const std::vector<TerminalNode> lean_nodes = terminal_nodes_from_quadrature(lean);
    QMOMKIT_VERIFY_CHECK(dense_nodes.size() == lean_nodes.size());
    for (std::size_t node = 0; node < dense_nodes.size(); ++node)
    {
        QMOMKIT_VERIFY_CHECK(dense_nodes[node].path == lean_nodes[node].path);
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(lean_nodes[node].weight, dense_nodes[node].weight, quadrature_abs_tol, quadrature_rel_tol);
        for (std::size_t dimension = 0; dimension < dense_nodes[node].coordinate.size(); ++dimension)
        {
            QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(lean_nodes[node].coordinate[dimension], dense_nodes[node].coordinate[dimension], quadrature_abs_tol,
                                                quadrature_rel_tol);
        }
    }
}

void compare_diagnostics(const std::vector<qmomkit::CqmomDiagnostics>& dense, const std::vector<qmomkit::CqmomDiagnostics>& lean)
{
    QMOMKIT_VERIFY_CHECK(dense.size() == lean.size());
    for (std::size_t i = 0; i < dense.size(); ++i)
    {
        QMOMKIT_VERIFY_CHECK(dense[i].conditioning_level == lean[i].conditioning_level);
        QMOMKIT_VERIFY_CHECK(dense[i].parent_path == lean[i].parent_path);
        QMOMKIT_VERIFY_CHECK(dense[i].phase_space_dimension == lean[i].phase_space_dimension);
        QMOMKIT_VERIFY_CHECK(dense[i].requested_node_count == lean[i].requested_node_count);
        QMOMKIT_VERIFY_CHECK(dense[i].actual_node_count == lean[i].actual_node_count);
        QMOMKIT_VERIFY_CHECK(dense[i].adaptive_branch_counts_enabled == lean[i].adaptive_branch_counts_enabled);
        QMOMKIT_VERIFY_CHECK(dense[i].adaptive_node_count_reduced == lean[i].adaptive_node_count_reduced);
        QMOMKIT_VERIFY_CHECK(dense[i].qmom_status == lean[i].qmom_status);
        QMOMKIT_VERIFY_CHECK(dense[i].conditional_moment_status == lean[i].conditional_moment_status);
        QMOMKIT_VERIFY_CHECK(dense[i].conditional_vandermonde_requested_solver == lean[i].conditional_vandermonde_requested_solver);
        QMOMKIT_VERIFY_CHECK(dense[i].conditional_vandermonde_selected_solver == lean[i].conditional_vandermonde_selected_solver);
        QMOMKIT_VERIFY_CHECK(dense[i].dense_vandermonde_fallback_occurred == lean[i].dense_vandermonde_fallback_occurred);
        QMOMKIT_VERIFY_CHECK(dense[i].message == lean[i].message);
    }
}

std::vector<qmomkit::CqmomNodeCountTransition> transitions_for_stage(const std::vector<qmomkit::CqmomNodeCountTransition>& transitions,
                                                                     const qmomkit::PhaseSpaceEvolutionStage stage)
{
    std::vector<qmomkit::CqmomNodeCountTransition> selected;
    for (const qmomkit::CqmomNodeCountTransition& transition : transitions)
    {
        if (transition.stage == stage)
        {
            selected.push_back(transition);
        }
    }
    return selected;
}

void compare_node_count_transitions(const std::vector<qmomkit::CqmomNodeCountTransition>& dense,
                                    const std::vector<qmomkit::CqmomNodeCountTransition>& lean)
{
    QMOMKIT_VERIFY_CHECK(dense.size() == lean.size());
    for (std::size_t i = 0; i < dense.size(); ++i)
    {
        QMOMKIT_VERIFY_CHECK(dense[i].stage == lean[i].stage);
        QMOMKIT_VERIFY_CHECK(dense[i].conditioning_level == lean[i].conditioning_level);
        QMOMKIT_VERIFY_CHECK(dense[i].parent_path == lean[i].parent_path);
        QMOMKIT_VERIFY_CHECK(dense[i].phase_space_dimension == lean[i].phase_space_dimension);
        QMOMKIT_VERIFY_CHECK(dense[i].configured_maximum_node_count == lean[i].configured_maximum_node_count);
        QMOMKIT_VERIFY_CHECK(dense[i].active_node_count_before == lean[i].active_node_count_before);
        QMOMKIT_VERIFY_CHECK(dense[i].active_node_count_after == lean[i].active_node_count_after);
        QMOMKIT_VERIFY_CHECK(dense[i].kind == lean[i].kind);
        QMOMKIT_VERIFY_CHECK(dense[i].rank_reduced == lean[i].rank_reduced);
        QMOMKIT_VERIFY_CHECK(dense[i].rank_grew == lean[i].rank_grew);
    }
}

bool diagnostics_request_solver(const std::vector<qmomkit::CqmomDiagnostics>& diagnostics, const std::size_t dimensions,
                                const qmomkit::VandermondeSolver solver)
{
    bool saw_conditional_branch = false;
    for (const qmomkit::CqmomDiagnostics& diagnostic : diagnostics)
    {
        if (diagnostic.conditioning_level + 1 >= dimensions)
        {
            continue;
        }
        saw_conditional_branch = true;
        if (diagnostic.conditional_vandermonde_requested_solver != solver || diagnostic.conditional_vandermonde_selected_solver != solver)
        {
            return false;
        }
    }
    return saw_conditional_branch;
}

bool all_diagnostics_enable_adaptive_branch_counts(const std::vector<qmomkit::CqmomDiagnostics>& diagnostics)
{
    if (diagnostics.empty())
    {
        return false;
    }
    for (const qmomkit::CqmomDiagnostics& diagnostic : diagnostics)
    {
        if (!diagnostic.adaptive_branch_counts_enabled)
        {
            return false;
        }
    }
    return true;
}

void compare_common_reconstructed_moments(const qmomkit::CqmomResult& dense, const qmomkit::CqmomResult& lean,
                                          const qmomkit::CqmomMomentLayout& layout)
{
    QMOMKIT_VERIFY_CHECK(dense.ok());
    QMOMKIT_VERIFY_CHECK(lean.ok());
    const qmomkit::CqmomMomentData dense_common = qmomkit::reconstruct_cqmom_moments(dense.quadrature, layout);
    const qmomkit::CqmomMomentData lean_common = qmomkit::reconstruct_cqmom_moments(lean.quadrature, layout);
    QMOMKIT_VERIFY_CHECK(dense_common.values.size() == lean_common.values.size());
    for (std::size_t i = 0; i < dense_common.values.size(); ++i)
    {
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(lean_common.values[i], dense_common.values[i], value_abs_tol, value_rel_tol);
    }
}

void compare_results(const std::string& name, const qmomkit::CqmomResult& dense, const qmomkit::CqmomResult& lean,
                     const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    std::cout << "[dense/lean parity] " << name << ": dense_status=" << static_cast<int>(dense.status)
              << " lean_status=" << static_cast<int>(lean.status) << '\n';
    QMOMKIT_VERIFY_CHECK(dense.status == lean.status);
    QMOMKIT_VERIFY_CHECK(dense.ok() == lean.ok());
    QMOMKIT_VERIFY_CHECK(dense.failed_phase_space_dimension == lean.failed_phase_space_dimension);
    QMOMKIT_VERIFY_CHECK(dense.message == lean.message);
    compare_diagnostics(dense.diagnostics, lean.diagnostics);
    if (dense.ok())
    {
        compare_quadratures(dense.quadrature, lean.quadrature);
        compare_common_reconstructed_moments(dense, lean, qmomkit::create_cqmom_moment_layout(hierarchy));
    }
}

std::vector<TerminalNode> deterministic_rectangular_nodes()
{
    const std::array<double, 2> d_nodes{0.82, 1.66};
    const std::array<double, 2> d_weights{0.43, 0.57};
    const std::array<std::array<double, 2>, 2> chi_nodes{{{{0.30, 0.94}}, {{0.48, 1.32}}}};
    const std::array<std::array<double, 2>, 2> chi_weights{{{{0.37, 0.63}}, {{0.59, 0.41}}}};
    const std::array<std::array<std::array<double, 2>, 2>, 2> u_nodes{{{{{{-0.35, 1.15}}, {{0.15, 1.80}}}}, {{{{-0.85, 0.90}}, {{0.42, 2.35}}}}}};
    const std::array<std::array<std::array<double, 2>, 2>, 2> u_weights{{{{{{0.62, 0.38}}, {{0.44, 0.56}}}}, {{{{0.71, 0.29}}, {{0.52, 0.48}}}}}};

    std::vector<TerminalNode> nodes;
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            for (std::size_t k = 0; k < 2; ++k)
            {
                nodes.push_back(TerminalNode{d_weights[i] * chi_weights[i][j] * u_weights[i][j][k],
                                             {d_nodes[i], chi_nodes[i][j], u_nodes[i][j][k]},
                                             qmomkit::ConditioningPath{{i, j, k}}});
            }
        }
    }
    return nodes;
}

std::vector<TerminalNode> independent_rectangular_nodes()
{
    const std::array<double, 2> d_nodes{0.76, 1.58};
    const std::array<double, 2> d_weights{0.46, 0.54};
    const std::array<double, 2> chi_nodes{0.28, 1.18};
    const std::array<double, 2> chi_weights{0.39, 0.61};
    const std::array<double, 2> u_nodes{-0.52, 1.44};
    const std::array<double, 2> u_weights{0.57, 0.43};

    std::vector<TerminalNode> nodes;
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            for (std::size_t k = 0; k < 2; ++k)
            {
                nodes.push_back(TerminalNode{
                    d_weights[i] * chi_weights[j] * u_weights[k], {d_nodes[i], chi_nodes[j], u_nodes[k]}, qmomkit::ConditioningPath{{i, j, k}}});
            }
        }
    }
    return nodes;
}

std::vector<TerminalNode> adaptive_rectangular_nodes()
{
    std::vector<TerminalNode> nodes = deterministic_rectangular_nodes();
    nodes.erase(std::remove_if(nodes.begin(), nodes.end(),
                               [](const TerminalNode& node) { return node.path.node_indices_by_level() == qmomkit::MultiIndex{0, 0, 1}; }),
                nodes.end());
    const double removed_weight = 0.43 * 0.37 * 0.38;
    for (TerminalNode& node : nodes)
    {
        if (node.path.node_indices_by_level() == qmomkit::MultiIndex{0, 0, 0})
        {
            node.weight += removed_weight;
            node.coordinate[2] = -0.35;
        }
    }
    return nodes;
}

std::vector<TerminalNode> near_degenerate_u_nodes()
{
    std::vector<TerminalNode> nodes = deterministic_rectangular_nodes();
    for (TerminalNode& node : nodes)
    {
        if (node.path.node_indices_by_level() == qmomkit::MultiIndex{1, 1, 0})
        {
            node.coordinate[2] = 0.7500000;
        }
        if (node.path.node_indices_by_level() == qmomkit::MultiIndex{1, 1, 1})
        {
            node.coordinate[2] = 0.7500015;
        }
    }
    return nodes;
}

void run_static_parity_case(const std::string& name, const qmomkit::PhaseSpaceHierarchy& hierarchy, const std::vector<TerminalNode>& nodes,
                            const bool adaptive_branch_counts)
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = adaptive_branch_counts;
    qmom_options.rmin = 1.0e-12;
    qmom_options.eabs = 1.0e-12;
    const qmomkit::RawMomentTensor moments = raw_moments_from_terminals(nodes, hierarchy.moment_extents());
    const qmomkit::CqmomResult dense = compute_dense(moments, hierarchy, qmom_options, adaptive_branch_counts);
    const qmomkit::CqmomResult lean = compute_lean(moments, hierarchy, qmom_options, adaptive_branch_counts);
    compare_results(name, dense, lean, hierarchy);
}

std::set<qmomkit::MultiIndex> index_set(const std::vector<qmomkit::MultiIndex>& indices)
{
    return std::set<qmomkit::MultiIndex>(indices.begin(), indices.end());
}

std::vector<qmomkit::MultiIndex> omitted_dense_indices(const qmomkit::CqmomMomentLayout& layout)
{
    const std::set<qmomkit::MultiIndex> required = index_set(layout.required_indices);
    std::vector<qmomkit::MultiIndex> omitted;
    for (const qmomkit::MultiIndex& index : qmomkit::enumerate_multi_indices(layout.hierarchy.moment_extents()))
    {
        if (required.find(index) == required.end())
        {
            omitted.push_back(index);
        }
    }
    return omitted;
}

double drag_lambda(const double d, const double chi)
{
    const double k_chi = 1.0 + 0.35 * chi + 0.08 * chi * chi;
    return k_chi / (d * d);
}

double drag_fluid_velocity(const double d, const double chi, const DragTransformParameters& parameters)
{
    return parameters.fluid_velocity_base + parameters.fluid_chi_coefficient * chi + parameters.fluid_d_coefficient * d;
}

void assert_valid_contraction(const DragTransformTrace& trace)
{
    QMOMKIT_VERIFY_CHECK(trace.lambda >= 0.0);
    QMOMKIT_VERIFY_CHECK(trace.time_step >= 0.0);
    QMOMKIT_VERIFY_CHECK(trace.decay >= 0.0);
    QMOMKIT_VERIFY_CHECK(trace.decay <= 1.0);
    const double lower = std::min(trace.u_old, trace.fluid_velocity);
    const double upper = std::max(trace.u_old, trace.fluid_velocity);
    const double tolerance = 32.0 * std::numeric_limits<double>::epsilon() * std::max({1.0, std::abs(trace.u_old), std::abs(trace.fluid_velocity)});
    QMOMKIT_VERIFY_CHECK(trace.u_new >= lower - tolerance);
    QMOMKIT_VERIFY_CHECK(trace.u_new <= upper + tolerance);
}

DragTransformTrace drag_transform_trace(const TerminalNode& node, const DragTransformParameters& parameters)
{
    DragTransformTrace trace;
    trace.path = node.path;
    trace.u_old = node.coordinate[2];
    trace.d = node.coordinate[0];
    trace.chi = node.coordinate[1];
    trace.lambda = drag_lambda(trace.d, trace.chi);
    trace.time_step = parameters.time_step;
    trace.lambda_dt = trace.lambda * trace.time_step;
    trace.decay = std::exp(-trace.lambda_dt);
    trace.fluid_velocity = drag_fluid_velocity(trace.d, trace.chi, parameters);
    trace.u_new = trace.fluid_velocity + (trace.u_old - trace.fluid_velocity) * trace.decay;
    assert_valid_contraction(trace);
    return trace;
}

std::vector<TerminalNode> transformed_nodes(const std::vector<TerminalNode>& nodes, const DragTransformParameters& parameters,
                                            std::vector<DragTransformTrace>* traces = nullptr)
{
    std::vector<TerminalNode> transformed = nodes;
    for (TerminalNode& node : transformed)
    {
        const DragTransformTrace trace = drag_transform_trace(node, parameters);
        node.coordinate[2] = trace.u_new;
        if (traces != nullptr)
        {
            traces->push_back(trace);
        }
    }
    return transformed;
}

DragTransformParameters physical_drag_parameters()
{
    return DragTransformParameters{0.7, 0.35, 0.25, -0.10, "physical-order-one analytical contraction"};
}

std::vector<DragTransformParameters> large_offset_stress_drag_parameters()
{
    return {
        {8.0, 1.0e8, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {14.0, 1.0e8, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {20.0, 1.0e8, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {28.0, 1.0e8, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {36.0, 1.0e8, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {8.0, 1.0e10, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {14.0, 1.0e10, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {20.0, 1.0e10, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {28.0, 1.0e10, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {36.0, 1.0e10, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {8.0, 1.0e12, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {14.0, 1.0e12, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {20.0, 1.0e12, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {28.0, 1.0e12, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
        {36.0, 1.0e12, 25000.0, -10000.0, "large-offset tiny-spread raw-moment stress"},
    };
}

void run_physical_analytical_contraction_case(const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = false;
    qmom_options.rmin = 1.0e-12;
    qmom_options.eabs = 1.0e-12;

    const DragTransformParameters parameters = physical_drag_parameters();
    std::vector<DragTransformTrace> traces;
    std::vector<TerminalNode> initial = deterministic_rectangular_nodes();
    std::sort(initial.begin(), initial.end(), path_less);
    const std::vector<TerminalNode> transformed = transformed_nodes(initial, parameters, &traces);
    const qmomkit::RawMomentTensor moments = raw_moments_from_terminals(transformed, hierarchy.moment_extents());
    const qmomkit::CqmomResult dense = compute_dense(moments, hierarchy, qmom_options, false);
    const qmomkit::CqmomResult lean = compute_lean(moments, hierarchy, qmom_options, false);
    compare_results("physical analytical contraction 3D 2,2,2", dense, lean, hierarchy);
    QMOMKIT_VERIFY_CHECK(lean.ok());

    const qmomkit::ConditioningPath branch{{0, 1}};
    const BranchDirectTrace direct = direct_trace_for_branch(transformed, branch, 2, 4);
    const qmomkit::WheelerOptions wheeler_options{qmom_options.adaptive, qmom_options.rmin, qmom_options.eabs, qmom_options.cutoff};
    const qmomkit::Quadrature1D direct_wheeler = qmomkit::compute_quadrature(qmomkit::MomentSet(direct.normalized_moments), 2, wheeler_options);
    QMOMKIT_VERIFY_CHECK(direct_wheeler.ok());
    std::cout << "[physical contraction] branch=(0,1) direct Wheeler status=" << static_cast<int>(direct_wheeler.status)
              << " initial_u_order=O(1) fluid_coefficients=(" << std::scientific << std::setprecision(17) << parameters.fluid_velocity_base << ", "
              << parameters.fluid_chi_coefficient << ", " << parameters.fluid_d_coefficient << ") time_step=" << parameters.time_step << '\n';
}

std::optional<FailureProbe> find_large_offset_failure_probe(const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = false;
    qmom_options.rmin = 1.0e-12;
    qmom_options.eabs = 1.0e-12;

    for (const DragTransformParameters& parameters : large_offset_stress_drag_parameters())
    {
        const qmomkit::RawMomentTensor initial_moments = raw_moments_from_terminals(deterministic_rectangular_nodes(), hierarchy.moment_extents());
        const qmomkit::CqmomResult initial = compute_lean(initial_moments, hierarchy, qmom_options, false);
        QMOMKIT_VERIFY_CHECK(initial.ok());
        const std::vector<TerminalNode> public_initial_nodes = terminal_nodes_from_quadrature(initial.quadrature);
        std::vector<TerminalNode> current = deterministic_rectangular_nodes();
        std::sort(current.begin(), current.end(), path_less);
        QMOMKIT_VERIFY_CHECK(current.size() == public_initial_nodes.size());
        for (std::size_t node = 0; node < current.size(); ++node)
        {
            QMOMKIT_VERIFY_CHECK(current[node].path == public_initial_nodes[node].path);
        }
        for (std::size_t macro_step = 1; macro_step <= 30; ++macro_step)
        {
            const std::vector<TerminalNode> previous = current;
            const std::vector<TerminalNode> next = transformed_nodes(current, parameters);
            const qmomkit::RawMomentTensor dense_moments = raw_moments_from_terminals(next, hierarchy.moment_extents());
            const qmomkit::CqmomResult dense = compute_dense(dense_moments, hierarchy, qmom_options, false);
            const qmomkit::CqmomResult lean = compute_lean(dense_moments, hierarchy, qmom_options, false);
            compare_results("large-offset transformed stress probe step", dense, lean, hierarchy);
            if (!lean.ok())
            {
                return FailureProbe{dense_moments,
                                    next,
                                    previous,
                                    dense,
                                    lean,
                                    macro_step,
                                    parameters.time_step,
                                    parameters.fluid_velocity_base,
                                    parameters.fluid_chi_coefficient,
                                    parameters.fluid_d_coefficient};
            }
            current = terminal_nodes_from_quadrature(lean.quadrature);
        }
    }
    return std::nullopt;
}

bool path_has_prefix(const qmomkit::ConditioningPath& path, const qmomkit::ConditioningPath& prefix)
{
    const qmomkit::MultiIndex& path_indices = path.node_indices_by_level();
    const qmomkit::MultiIndex& prefix_indices = prefix.node_indices_by_level();
    if (prefix_indices.size() > path_indices.size())
    {
        return false;
    }
    return std::equal(prefix_indices.begin(), prefix_indices.end(), path_indices.begin());
}

BranchSigmaDiagnostic branch_sigma_direct_from_terminals(const std::vector<TerminalNode>& nodes, const qmomkit::CqmomDiagnostics& diagnostic,
                                                         const qmomkit::WheelerOptions& options)
{
    const std::size_t moment_count = 2 * diagnostic.requested_node_count;
    BranchSigmaDiagnostic branch;
    branch.moments.assign(moment_count, 0.0);
    double parent_mass = 0.0;
    for (const TerminalNode& node : nodes)
    {
        if (path_has_prefix(node.path, diagnostic.parent_path))
        {
            parent_mass += node.weight;
        }
    }
    for (const TerminalNode& node : nodes)
    {
        if (path_has_prefix(node.path, diagnostic.parent_path))
        {
            const double local_weight = node.weight / parent_mass;
            for (std::size_t order = 0; order < moment_count; ++order)
            {
                branch.moments[order] += local_weight * integer_power(node.coordinate[diagnostic.phase_space_dimension], order);
            }
        }
    }
    branch.sigma22 = branch.moments[2] - branch.moments[1] * branch.moments[1] / branch.moments[0];
    double scale = 1.0;
    for (const double moment : branch.moments)
    {
        scale = std::max(scale, std::abs(moment));
    }
    branch.tolerance = std::max(options.cutoff, 100.0 * std::numeric_limits<double>::epsilon() * scale);
    return branch;
}

std::array<double, 2> solve_conditional_vandermonde_system_2(const std::array<double, 2>& nodes, const std::array<double, 2>& weights,
                                                             const std::array<double, 2>& rhs)
{
    const double separation = nodes[1] - nodes[0];
    return {(rhs[0] * nodes[1] - rhs[1]) / (weights[0] * separation), (rhs[1] - rhs[0] * nodes[0]) / (weights[1] * separation)};
}

double dense_moment(const qmomkit::RawMomentTensor& moments, const qmomkit::MultiIndex& index)
{
    return moments.values[qmomkit::flatten_multi_index(index, moments.extents)];
}

Recurrence2Diagnostic recurrence2(const std::vector<double>& moments)
{
    QMOMKIT_VERIFY_CHECK(moments.size() >= 4);
    Recurrence2Diagnostic diagnostic;
    diagnostic.sigma11 = moments[0];
    diagnostic.sigma12 = moments[1];
    diagnostic.sigma13 = moments[2];
    diagnostic.a0 = moments[1] / moments[0];
    diagnostic.b0 = 0.0;
    diagnostic.sigma22 = moments[2] - diagnostic.a0 * moments[1];
    diagnostic.sigma23 = moments[3] - diagnostic.a0 * moments[2];
    diagnostic.b1 = diagnostic.sigma22 / diagnostic.sigma11;
    diagnostic.a1 = diagnostic.sigma23 / diagnostic.sigma22 - diagnostic.sigma12 / diagnostic.sigma11;
    diagnostic.variance = moments[2] / moments[0] - (moments[1] / moments[0]) * (moments[1] / moments[0]);
    return diagnostic;
}

BranchDirectTrace direct_trace_for_branch(const std::vector<TerminalNode>& nodes, const qmomkit::ConditioningPath& branch_path,
                                          const std::size_t dimension, const std::size_t moment_count)
{
    BranchDirectTrace trace;
    trace.normalized_moments.assign(moment_count, 0.0);
    trace.unnormalized_moments.assign(moment_count, 0.0);
    trace.normalized_moments_long.assign(moment_count, 0.0L);

    for (std::size_t flat = 0; flat < nodes.size(); ++flat)
    {
        if (path_has_prefix(nodes[flat].path, branch_path))
        {
            trace.flat_terminal_indices.push_back(flat);
            trace.terminal_nodes.push_back(nodes[flat]);
            trace.branch_mass += nodes[flat].weight;
        }
    }
    QMOMKIT_VERIFY_CHECK(trace.branch_mass > 0.0);

    long double branch_mass_long = 0.0L;
    for (const TerminalNode& node : trace.terminal_nodes)
    {
        branch_mass_long += static_cast<long double>(node.weight);
    }

    long double mean_long = 0.0L;
    for (const TerminalNode& node : trace.terminal_nodes)
    {
        const double conditional_weight = node.weight / trace.branch_mass;
        trace.conditional_weights.push_back(conditional_weight);
        mean_long += static_cast<long double>(conditional_weight) * static_cast<long double>(node.coordinate[dimension]);
        for (std::size_t order = 0; order < moment_count; ++order)
        {
            const double power = integer_power(node.coordinate[dimension], order);
            trace.normalized_moments[order] += conditional_weight * power;
            trace.unnormalized_moments[order] += node.weight * power;
            trace.normalized_moments_long[order] += (static_cast<long double>(node.weight) / branch_mass_long) *
                                                    integer_power_long(static_cast<long double>(node.coordinate[dimension]), order);
        }
    }
    for (std::size_t local = 0; local < trace.terminal_nodes.size(); ++local)
    {
        const long double delta = static_cast<long double>(trace.terminal_nodes[local].coordinate[dimension]) - mean_long;
        trace.variance_long += static_cast<long double>(trace.conditional_weights[local]) * delta * delta;
    }
    return trace;
}

double direct_conditioned_child_moment(const std::vector<TerminalNode>& nodes, const qmomkit::ConditioningPath& branch_path,
                                       const std::size_t chi_order, const std::size_t u_order)
{
    double branch_mass = 0.0;
    for (const TerminalNode& node : nodes)
    {
        if (path_has_prefix(node.path, branch_path))
        {
            branch_mass += node.weight;
        }
    }
    double moment = 0.0;
    for (const TerminalNode& node : nodes)
    {
        if (path_has_prefix(node.path, branch_path))
        {
            moment += (node.weight / branch_mass) * integer_power(node.coordinate[1], chi_order) * integer_power(node.coordinate[2], u_order);
        }
    }
    return moment;
}

VandermondeStageTrace solve_stage(const qmomkit::internal::BranchMomentSlice& parent_slice, const std::vector<double>& parent_abscissas,
                                  const std::vector<double>& parent_weights, const qmomkit::ConditionalMomentOptions& conditional_options)
{
    VandermondeStageTrace trace;
    trace.parent_abscissas = parent_abscissas;
    trace.parent_weights = parent_weights;
    trace.rhs = qmomkit::internal::make_branch_rhs_matrix(parent_slice, parent_abscissas.size());
    const qmomkit::internal::VandermondeSolveOptions solve_options{qmomkit::VandermondeSolver::automatic, conditional_options.singular_tolerance,
                                                                   true, 1.0e-10};
    trace.solve = qmomkit::internal::solve_conditional_vandermonde_system(parent_abscissas, parent_weights, trace.rhs, solve_options,
                                                                          "Test-local weighted Vandermonde system is singular.");
    QMOMKIT_VERIFY_CHECK(trace.solve.success);
    const qmomkit::internal::DenseMatrix system = qmomkit::internal::make_conditional_vandermonde_matrix(parent_abscissas, parent_weights);
    trace.residual = qmomkit::internal::residual_norm(system, trace.solve.solution, trace.rhs);
    trace.condition = qmomkit::internal::condition_estimate(system, conditional_options.singular_tolerance);
    return trace;
}

RecursiveBranchTrace recursive_trace_for_identity_branch(const qmomkit::RawMomentTensor& moments, const qmomkit::WheelerOptions& wheeler_options,
                                                         const qmomkit::ConditionalMomentOptions& conditional_options)
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {0, 1, 2});
    const qmomkit::CqmomMomentLayout canonical = qmomkit::create_cqmom_moment_layout(hierarchy);
    const std::shared_ptr<const qmomkit::internal::BranchMomentLayout> root_layout =
        qmomkit::internal::create_branch_moment_layout(hierarchy.moment_extents(), canonical.required_indices, false);
    qmomkit::internal::BranchMomentSlice root_slice;
    root_slice.layout = root_layout.get();
    root_slice.values.reserve(canonical.required_indices.size());
    for (const qmomkit::MultiIndex& index : canonical.required_indices)
    {
        root_slice.values.push_back(dense_moment(moments, index));
    }

    RecursiveBranchTrace trace;
    trace.d_moments = qmomkit::internal::extract_wheeler_moments(root_slice);
    const qmomkit::Quadrature1D d_quadrature = qmomkit::compute_quadrature(qmomkit::MomentSet(trace.d_moments), 2, wheeler_options);
    QMOMKIT_VERIFY_CHECK(d_quadrature.ok());
    trace.d_nodes = d_quadrature.abscissas;
    trace.d_weights = d_quadrature.weights;

    trace.d_to_chi_u = solve_stage(root_slice, trace.d_nodes, trace.d_weights, conditional_options);
    const qmomkit::internal::BranchMomentChildLayout& root_child_case = qmomkit::internal::child_layout_for_parent_node_count(*root_slice.layout, 2);
    qmomkit::internal::BranchMomentSlice d0_slice{root_child_case.child_layout.get(),
                                                  std::vector<double>(root_child_case.child_layout->indices.size(), 0.0)};
    trace.d0_child_indices = root_child_case.child_layout->indices;
    for (std::size_t col = 0; col < d0_slice.values.size(); ++col)
    {
        d0_slice.values[col] = trace.d_to_chi_u.solve.solution(0, col);
    }
    trace.d0_child_values = d0_slice.values;

    trace.d0_chi_moments = qmomkit::internal::extract_wheeler_moments(d0_slice);
    const qmomkit::Quadrature1D chi_quadrature = qmomkit::compute_quadrature(qmomkit::MomentSet(trace.d0_chi_moments), 2, wheeler_options);
    QMOMKIT_VERIFY_CHECK(chi_quadrature.ok());
    trace.d0_chi_nodes = chi_quadrature.abscissas;
    trace.d0_chi_weights = chi_quadrature.weights;

    trace.chi_to_u = solve_stage(d0_slice, trace.d0_chi_nodes, trace.d0_chi_weights, conditional_options);
    const qmomkit::internal::BranchMomentChildLayout& chi_child_case = qmomkit::internal::child_layout_for_parent_node_count(*d0_slice.layout, 2);
    qmomkit::internal::BranchMomentSlice d0_chi1_slice{chi_child_case.child_layout.get(),
                                                       std::vector<double>(chi_child_case.child_layout->indices.size(), 0.0)};
    trace.d0_chi1_child_indices = chi_child_case.child_layout->indices;
    for (std::size_t col = 0; col < d0_chi1_slice.values.size(); ++col)
    {
        d0_chi1_slice.values[col] = trace.chi_to_u.solve.solution(1, col);
    }
    trace.normalized_moments = qmomkit::internal::extract_wheeler_moments(d0_chi1_slice);
    return trace;
}

BranchSigmaDiagnostic branch_sigma_reconstructed_for_identity_order(const qmomkit::RawMomentTensor& moments,
                                                                    const qmomkit::CqmomDiagnostics& diagnostic,
                                                                    const qmomkit::WheelerOptions& options)
{
    BranchSigmaDiagnostic branch;
    branch.moments.assign(2 * diagnostic.requested_node_count, 0.0);
    if (diagnostic.conditioning_level != 2 || diagnostic.parent_path.node_indices_by_level().size() != 2)
    {
        return branch;
    }

    std::vector<double> d_moments;
    for (std::size_t order = 0; order < 4; ++order)
    {
        d_moments.push_back(dense_moment(moments, {order, 0, 0}));
    }
    const qmomkit::Quadrature1D d_quadrature = qmomkit::compute_quadrature(qmomkit::MomentSet(d_moments), 2, options);
    if (!d_quadrature.ok())
    {
        return branch;
    }
    const std::array<double, 2> d_nodes{d_quadrature.abscissas[0], d_quadrature.abscissas[1]};
    const std::array<double, 2> d_weights{d_quadrature.weights[0], d_quadrature.weights[1]};

    double child_moments[2][4][4]{};
    for (std::size_t chi_order = 0; chi_order < 4; ++chi_order)
    {
        for (std::size_t u_order = 0; u_order < 4; ++u_order)
        {
            const std::array<double, 2> rhs{dense_moment(moments, {0, chi_order, u_order}), dense_moment(moments, {1, chi_order, u_order})};
            const std::array<double, 2> values = solve_conditional_vandermonde_system_2(d_nodes, d_weights, rhs);
            child_moments[0][chi_order][u_order] = values[0];
            child_moments[1][chi_order][u_order] = values[1];
        }
    }

    const std::size_t d_branch = diagnostic.parent_path.node_indices_by_level()[0];
    const std::size_t chi_branch = diagnostic.parent_path.node_indices_by_level()[1];
    std::vector<double> chi_moments;
    for (std::size_t order = 0; order < 4; ++order)
    {
        chi_moments.push_back(child_moments[d_branch][order][0]);
    }
    const qmomkit::Quadrature1D chi_quadrature = qmomkit::compute_quadrature(qmomkit::MomentSet(chi_moments), 2, options);
    if (!chi_quadrature.ok())
    {
        return branch;
    }
    const std::array<double, 2> chi_nodes{chi_quadrature.abscissas[0], chi_quadrature.abscissas[1]};
    const std::array<double, 2> chi_weights{chi_quadrature.weights[0], chi_quadrature.weights[1]};

    for (std::size_t u_order = 0; u_order < branch.moments.size(); ++u_order)
    {
        const std::array<double, 2> rhs{child_moments[d_branch][0][u_order], child_moments[d_branch][1][u_order]};
        const std::array<double, 2> values = solve_conditional_vandermonde_system_2(chi_nodes, chi_weights, rhs);
        branch.moments[u_order] = values[chi_branch];
    }
    branch.sigma22 = branch.moments[2] - branch.moments[1] * branch.moments[1] / branch.moments[0];
    double scale = 1.0;
    for (const double moment : branch.moments)
    {
        scale = std::max(scale, std::abs(moment));
    }
    branch.tolerance = std::max(options.cutoff, 100.0 * std::numeric_limits<double>::epsilon() * scale);
    return branch;
}

std::optional<qmomkit::CqmomDiagnostics> first_failed_branch(const qmomkit::CqmomResult& result)
{
    for (const qmomkit::CqmomDiagnostics& diagnostic : result.diagnostics)
    {
        if (diagnostic.qmom_status != qmomkit::QmomStatus::success ||
            diagnostic.conditional_moment_status != qmomkit::ConditionalMomentStatus::success)
        {
            return diagnostic;
        }
    }
    return std::nullopt;
}

void print_vector(const std::string& label, const std::vector<double>& values)
{
    std::cout << label;
    for (const double value : values)
    {
        std::cout << ' ' << std::scientific << std::setprecision(17) << value;
    }
    std::cout << '\n';
}

void print_recurrence(const std::string& label, const Recurrence2Diagnostic& recurrence)
{
    std::cout << label << " sigma11=" << std::scientific << std::setprecision(17) << recurrence.sigma11 << " sigma12=" << recurrence.sigma12
              << " sigma13=" << recurrence.sigma13 << " sigma22=" << recurrence.sigma22 << " sigma23=" << recurrence.sigma23
              << " variance=" << recurrence.variance << " a=[" << recurrence.a0 << ", " << recurrence.a1 << "] b=[" << recurrence.b0 << ", "
              << recurrence.b1 << "]\n";
}

void print_stage(const std::string& label, const VandermondeStageTrace& trace)
{
    std::cout << label << " parent_nodes=";
    for (const double node : trace.parent_abscissas)
    {
        std::cout << ' ' << std::scientific << std::setprecision(17) << node;
    }
    std::cout << " parent_weights=";
    for (const double weight : trace.parent_weights)
    {
        std::cout << ' ' << std::scientific << std::setprecision(17) << weight;
    }
    std::cout << " solver_requested=" << static_cast<int>(trace.solve.requested_solver)
              << " solver_selected=" << static_cast<int>(trace.solve.selected_solver) << " algorithm=" << static_cast<int>(trace.solve.algorithm_used)
              << " residual=" << std::scientific << std::setprecision(17) << trace.residual << " condition=" << trace.condition << '\n';
}

void print_drag_inputs_for_branch(const FailureProbe& probe, const qmomkit::ConditioningPath& branch_path)
{
    const DragTransformParameters parameters{probe.time_step, probe.fluid_velocity_base, probe.fluid_chi_coefficient, probe.fluid_d_coefficient,
                                             "large-offset tiny-spread raw-moment stress"};
    std::cout << "[large-offset contraction inputs] branch=" << path_string(branch_path) << " formula=u_f + (u_old - u_f) * exp(-lambda * dt)"
              << " fluid_base=" << std::scientific << std::setprecision(17) << parameters.fluid_velocity_base
              << " fluid_chi_coefficient=" << parameters.fluid_chi_coefficient << " fluid_d_coefficient=" << parameters.fluid_d_coefficient << '\n';
    for (const TerminalNode& previous : probe.previous_terminal_nodes)
    {
        if (path_has_prefix(previous.path, branch_path))
        {
            const DragTransformTrace trace = drag_transform_trace(previous, parameters);
            std::cout << "[large-offset contraction input] path=" << path_string(trace.path) << " u_old=" << std::scientific << std::setprecision(17)
                      << trace.u_old << " u_f=" << trace.fluid_velocity << " d=" << trace.d << " chi=" << trace.chi << " lambda=" << trace.lambda
                      << " dt=" << trace.time_step << " lambda_dt=" << trace.lambda_dt << " exp(-lambda_dt)=" << trace.decay
                      << " u_new=" << trace.u_new << '\n';
        }
    }
}

void audit_failing_branch_01(const FailureProbe& probe, const qmomkit::CqmomDiagnostics& failed, const qmomkit::WheelerOptions& wheeler_options,
                             const qmomkit::ConditionalMomentOptions& conditional_options)
{
    const qmomkit::ConditioningPath requested_branch{{0, 1}};
    QMOMKIT_VERIFY_CHECK(failed.parent_path == requested_branch);
    print_drag_inputs_for_branch(probe, requested_branch);
    const BranchDirectTrace direct =
        direct_trace_for_branch(probe.terminal_nodes, failed.parent_path, failed.phase_space_dimension, 2 * failed.requested_node_count);
    QMOMKIT_VERIFY_CHECK(direct.terminal_nodes.size() == 2);
    QMOMKIT_VERIFY_CHECK(direct.conditional_weights.size() == direct.terminal_nodes.size());

    double conditional_weight_sum = 0.0;
    for (const double weight : direct.conditional_weights)
    {
        conditional_weight_sum += weight;
        QMOMKIT_VERIFY_CHECK(weight > 0.0);
    }
    QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(conditional_weight_sum, 1.0, 1.0e-15, 1.0e-15);

    std::cout << "[branch (0,1) terminals] count=" << direct.terminal_nodes.size() << " branch_mass=" << std::scientific << std::setprecision(17)
              << direct.branch_mass << " conditional_weight_sum=" << conditional_weight_sum << '\n';
    for (std::size_t i = 0; i < direct.terminal_nodes.size(); ++i)
    {
        const TerminalNode& node = direct.terminal_nodes[i];
        std::cout << "[branch (0,1) terminal] flat=" << direct.flat_terminal_indices[i] << " path=" << path_string(node.path)
                  << " global_weight=" << std::scientific << std::setprecision(17) << node.weight
                  << " conditional_weight=" << direct.conditional_weights[i] << " d=" << node.coordinate[0] << " chi=" << node.coordinate[1]
                  << " transformed_u=" << node.coordinate[2] << '\n';
    }

    print_vector("[branch (0,1) direct normalized moments]", direct.normalized_moments);
    print_vector("[branch (0,1) direct unnormalized moments]", direct.unnormalized_moments);
    std::cout << "[branch (0,1) direct long-double normalized moments]";
    for (const long double value : direct.normalized_moments_long)
    {
        std::cout << ' ' << std::scientific << std::setprecision(21) << value;
    }
    std::cout << " central_variance=" << std::scientific << std::setprecision(21) << direct.variance_long << '\n';

    const qmomkit::Quadrature1D direct_wheeler =
        qmomkit::compute_quadrature(qmomkit::MomentSet(direct.normalized_moments), failed.requested_node_count, wheeler_options);
    std::cout << "[branch (0,1) direct Wheeler] status=" << static_cast<int>(direct_wheeler.status) << " message=\"" << direct_wheeler.message
              << "\" node_count=" << direct_wheeler.weights.size() << '\n';
    print_vector("[branch (0,1) direct Wheeler weights]", direct_wheeler.weights);
    print_vector("[branch (0,1) direct Wheeler abscissas]", direct_wheeler.abscissas);

    qmomkit::WheelerOptions adaptive_wheeler_options = wheeler_options;
    adaptive_wheeler_options.adaptive = true;
    const qmomkit::Quadrature1D adaptive_direct_wheeler =
        qmomkit::compute_quadrature(qmomkit::MomentSet(direct.normalized_moments), failed.requested_node_count, adaptive_wheeler_options);
    const double node_separation = std::abs(direct.terminal_nodes[1].coordinate[failed.phase_space_dimension] -
                                            direct.terminal_nodes[0].coordinate[failed.phase_space_dimension]);
    const double node_scale = std::max(std::abs(direct.terminal_nodes[0].coordinate[failed.phase_space_dimension]),
                                       std::abs(direct.terminal_nodes[1].coordinate[failed.phase_space_dimension]));
    const double relative_node_separation = node_separation / std::max(1.0, node_scale);
    const double ulp_spacing = std::nextafter(node_scale, std::numeric_limits<double>::infinity()) - node_scale;
    const Recurrence2Diagnostic direct_recurrence = recurrence2(direct.normalized_moments);
    std::cout << "[branch (0,1) adaptive collapse] node_separation=" << std::scientific << std::setprecision(17) << node_separation
              << " relative_node_separation=" << relative_node_separation << " ulp_spacing=" << ulp_spacing
              << " long_double_true_variance=" << std::setprecision(21) << direct.variance_long << " double_raw_variance=" << std::setprecision(17)
              << direct_recurrence.variance << " eabs_spacing_threshold=" << adaptive_wheeler_options.eabs
              << " spacing_test_would_treat_distinct=" << (node_separation > adaptive_wheeler_options.eabs ? "true" : "false")
              << " adaptive_direct_status=" << static_cast<int>(adaptive_direct_wheeler.status)
              << " adaptive_direct_node_count=" << adaptive_direct_wheeler.weights.size() << '\n';

    const RecursiveBranchTrace recursive = recursive_trace_for_identity_branch(probe.dense_moments, wheeler_options, conditional_options);
    print_vector("[branch (0,1) recursive d moments]", recursive.d_moments);
    print_vector("[branch (0,1) recursive d nodes]", recursive.d_nodes);
    print_vector("[branch (0,1) recursive d weights]", recursive.d_weights);
    print_stage("[branch (0,1) stage d -> (chi,u)]", recursive.d_to_chi_u);
    print_vector("[branch (0,1) recursive chi|d0 moments]", recursive.d0_chi_moments);
    print_vector("[branch (0,1) recursive chi|d0 nodes]", recursive.d0_chi_nodes);
    print_vector("[branch (0,1) recursive chi|d0 weights]", recursive.d0_chi_weights);
    print_stage("[branch (0,1) stage chi -> u]", recursive.chi_to_u);
    print_vector("[branch (0,1) recursive normalized u moments]", recursive.normalized_moments);

    std::cout << "[branch (0,1) moment comparison] order direct_normalized recursive_normalized absolute_error relative_error unnormalized_error\n";
    for (std::size_t order = 0; order < direct.normalized_moments.size(); ++order)
    {
        const double direct_value = direct.normalized_moments[order];
        const double recursive_value = recursive.normalized_moments[order];
        const double absolute_error = std::abs(recursive_value - direct_value);
        const double relative_error = absolute_error / std::max(1.0, std::abs(direct_value));
        const double unnormalized_error = std::abs(recursive_value * direct.branch_mass - direct.unnormalized_moments[order]);
        std::cout << "[branch (0,1) moment comparison] " << order << ' ' << std::scientific << std::setprecision(17) << direct_value << ' '
                  << recursive_value << ' ' << absolute_error << ' ' << relative_error << ' ' << unnormalized_error << '\n';
    }

    const Recurrence2Diagnostic recursive_recurrence = recurrence2(recursive.normalized_moments);
    print_recurrence("[branch (0,1) direct recurrence]", direct_recurrence);
    print_recurrence("[branch (0,1) recursive recurrence]", recursive_recurrence);

    double max_stage1_error = 0.0;
    double max_stage1_relative_error = 0.0;
    qmomkit::MultiIndex max_stage1_index;
    for (const qmomkit::MultiIndex& local_index : recursive.d0_child_indices)
    {
        const double direct_value =
            direct_conditioned_child_moment(probe.terminal_nodes, qmomkit::ConditioningPath{{0}}, local_index[0], local_index[1]);
        const auto found = std::find(recursive.d0_child_indices.begin(), recursive.d0_child_indices.end(), local_index);
        const std::size_t offset = static_cast<std::size_t>(std::distance(recursive.d0_child_indices.begin(), found));
        const double error = std::abs(recursive.d0_child_values[offset] - direct_value);
        if (error > max_stage1_error)
        {
            max_stage1_error = error;
            max_stage1_relative_error = error / std::max(1.0, std::abs(direct_value));
            max_stage1_index = local_index;
        }
    }
    double max_stage2_error = 0.0;
    double max_stage2_relative_error = 0.0;
    qmomkit::MultiIndex max_stage2_index;
    for (std::size_t order = 0; order < recursive.normalized_moments.size(); ++order)
    {
        const double error = std::abs(recursive.normalized_moments[order] - direct.normalized_moments[order]);
        if (error > max_stage2_error)
        {
            max_stage2_error = error;
            max_stage2_relative_error = error / std::max(1.0, std::abs(direct.normalized_moments[order]));
            max_stage2_index = qmomkit::MultiIndex{order};
        }
    }
    std::cout << "[branch (0,1) divergence] stage1 max_abs_error=" << std::scientific << std::setprecision(17) << max_stage1_error
              << " max_relative_error=" << max_stage1_relative_error << " at local(chi,u)=" << index_string(max_stage1_index) << '\n';
    std::cout << "[branch (0,1) divergence] stage2 max_abs_error=" << std::scientific << std::setprecision(17) << max_stage2_error
              << " max_relative_error=" << max_stage2_relative_error << " at local(u)=" << index_string(max_stage2_index) << '\n';
    if (direct_wheeler.status != qmomkit::QuadratureStatus::success)
    {
        std::cout << "[branch (0,1) conclusion] Direct double-precision raw moments are already rejected by Wheeler; long-double central variance "
                     "from the same positive nodes is positive, so the failure is a raw-moment scaling/cancellation problem before any "
                     "dense-versus-lean distinction.\n";
    }
    else if (max_stage2_error > 1.0e-6 * std::max(1.0, std::abs(direct.normalized_moments[3])))
    {
        std::cout << "[branch (0,1) conclusion] Direct Wheeler succeeds, but recursive nested Vandermonde reconstruction materially changes the "
                     "u moment sequence.\n";
    }
    else
    {
        std::cout << "[branch (0,1) conclusion] Direct and recursive moments agree at the tested scale; any failure is in common Wheeler "
                     "realizability handling for the supplied raw sequence.\n";
    }
}

std::vector<qmomkit::MultiIndex> expanded_indices(const qmomkit::CqmomMomentLayout& canonical, const std::vector<qmomkit::MultiIndex>& additions)
{
    std::vector<qmomkit::MultiIndex> indices = canonical.required_indices;
    for (const qmomkit::MultiIndex& index : additions)
    {
        if (std::find(indices.begin(), indices.end(), index) == indices.end())
        {
            indices.push_back(index);
        }
    }
    return indices;
}

qmomkit::internal::RecursiveQuadratureResult compute_with_expanded_layout(const qmomkit::RawMomentTensor& moments,
                                                                          const qmomkit::CqmomMomentLayout& layout,
                                                                          const qmomkit::QmomOptions& qmom_options)
{
    std::vector<double> values;
    values.reserve(layout.required_indices.size());
    for (const qmomkit::MultiIndex& index : layout.required_indices)
    {
        values.push_back(dense_moment(moments, index));
    }

    const qmomkit::internal::IndexedMomentSet indexed(moments.extents, layout.required_indices, std::move(values));
    const qmomkit::internal::CqmomMomentView view(indexed);
    qmomkit::internal::RecursiveQuadratureMomentViewRequest request{layout.hierarchy, view};
    request.wheeler_options = qmomkit::WheelerOptions{qmom_options.adaptive, qmom_options.rmin, qmom_options.eabs, qmom_options.cutoff};
    request.diagnostic_level = qmomkit::CqmomDiagnosticLevel::detailed;
    request.conditional_vandermonde_solver = qmomkit::VandermondeSolver::automatic;
    request.validated_layout = &layout;
    return qmomkit::internal::compute_recursive_quadrature_from_moment_view(request);
}

void report_omitted_moments(const qmomkit::RawMomentTensor& before, const qmomkit::RawMomentTensor& after, const qmomkit::CqmomMomentLayout& layout)
{
    const std::vector<qmomkit::MultiIndex> omitted = omitted_dense_indices(layout);
    std::cout << "[dense/lean layout] canonical required moments=" << layout.required_indices.size() << " omitted_dense_moments=" << omitted.size()
              << '\n';

    double largest_change = 0.0;
    qmomkit::MultiIndex largest_change_index;
    for (const qmomkit::MultiIndex& index : omitted)
    {
        const double change = std::abs(dense_moment(after, index) - dense_moment(before, index));
        if (change > largest_change)
        {
            largest_change = change;
            largest_change_index = index;
        }
    }
    std::cout << "[dense/lean layout] largest omitted transformed change " << index_string(largest_change_index) << " = " << std::setprecision(17)
              << largest_change << '\n';

    const std::vector<qmomkit::MultiIndex> requested{{2, 1, 1}, {1, 2, 1}, {2, 2, 1}, {2, 1, 2}, {1, 2, 2}, {2, 2, 2}};
    for (const qmomkit::MultiIndex& index : requested)
    {
        QMOMKIT_VERIFY_CHECK(std::find(omitted.begin(), omitted.end(), index) != omitted.end());
        std::cout << "[dense/lean layout] omitted " << index_string(index) << " before=" << std::setprecision(17) << dense_moment(before, index)
                  << " after=" << dense_moment(after, index) << " delta=" << dense_moment(after, index) - dense_moment(before, index) << '\n';
    }
}

void run_expanded_layout_experiments(const qmomkit::RawMomentTensor& moments, const qmomkit::CqmomMomentLayout& canonical,
                                     const qmomkit::QmomOptions& qmom_options)
{
    const std::vector<std::pair<std::string, std::vector<qmomkit::MultiIndex>>> experiments{
        {"add (2,1,1)", {{2, 1, 1}}},
        {"add first-order u cross terms", {{2, 1, 1}, {1, 2, 1}, {2, 2, 1}}},
        {"add second-order u cross terms", {{2, 1, 2}, {1, 2, 2}, {2, 2, 2}}},
        {"add requested omitted set", {{2, 1, 1}, {1, 2, 1}, {2, 2, 1}, {2, 1, 2}, {1, 2, 2}, {2, 2, 2}}},
    };

    bool found_successful_expansion = false;
    for (const auto& experiment : experiments)
    {
        qmomkit::CqmomMomentLayout expanded{canonical.hierarchy, expanded_indices(canonical, experiment.second)};
        try
        {
            const qmomkit::internal::RecursiveQuadratureResult result = compute_with_expanded_layout(moments, expanded, qmom_options);
            std::cout << "[expanded layout] " << experiment.first << ": status=" << static_cast<int>(result.status)
                      << " diagnostics=" << result.diagnostics.size() << '\n';
            found_successful_expansion = found_successful_expansion || result.ok();
        }
        catch (const std::exception& error)
        {
            std::cout << "[expanded layout] " << experiment.first << ": rejected by internal layout construction: " << error.what() << '\n';
        }
    }
    if (!found_successful_expansion)
    {
        std::cout << "[expanded layout] no tested diagnostic expansion avoided the observed failure.\n";
    }
}

void run_failure_probe()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {0, 1, 2});
    const qmomkit::RawMomentTensor initial = raw_moments_from_terminals(deterministic_rectangular_nodes(), hierarchy.moment_extents());
    const std::optional<FailureProbe> probe = find_large_offset_failure_probe(hierarchy);
    if (!probe.has_value())
    {
        std::cout << "[large-offset transformed stress probe] no failure found in deterministic no-RK sweep; dense and lean remained equivalent.\n";
        return;
    }

    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData lean_from_dense = qmomkit::extract_cqmom_moments(probe->dense_moments, layout);
    const qmomkit::CqmomMomentData lean_from_nodes =
        qmomkit::extract_cqmom_moments(raw_moments_from_terminals(probe->terminal_nodes, hierarchy.moment_extents()), layout);
    for (std::size_t i = 0; i < lean_from_dense.values.size(); ++i)
    {
        QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(lean_from_dense.values[i], lean_from_nodes.values[i], value_abs_tol, value_rel_tol);
    }

    std::cout << "[large-offset transformed stress probe] reproduced non-success without RK at step=" << probe->macro_step
              << " dt=" << std::setprecision(17) << probe->time_step << " fluid_velocity_base=" << probe->fluid_velocity_base
              << " fluid_chi_coefficient=" << probe->fluid_chi_coefficient << " fluid_d_coefficient=" << probe->fluid_d_coefficient
              << " status=" << static_cast<int>(probe->lean_result.status) << '\n';

    report_omitted_moments(initial, probe->dense_moments, layout);

    const std::optional<qmomkit::CqmomDiagnostics> failed = first_failed_branch(probe->lean_result);
    QMOMKIT_VERIFY_CHECK(failed.has_value());
    std::cout << "[failure branch] dim=" << failed->phase_space_dimension << " level=" << failed->conditioning_level
              << " path=" << path_string(failed->parent_path) << " requested=" << failed->requested_node_count
              << " qmom_status=" << static_cast<int>(failed->qmom_status) << " message=\"" << failed->message << "\"\n";

    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = false;
    qmom_options.rmin = 1.0e-12;
    qmom_options.eabs = 1.0e-12;
    const qmomkit::WheelerOptions wheeler_options{qmom_options.adaptive, qmom_options.rmin, qmom_options.eabs, qmom_options.cutoff};
    const BranchSigmaDiagnostic direct = branch_sigma_direct_from_terminals(probe->terminal_nodes, *failed, wheeler_options);
    const BranchSigmaDiagnostic recursive = branch_sigma_reconstructed_for_identity_order(probe->dense_moments, *failed, wheeler_options);
    std::cout << "[failure branch] direct moments:";
    for (const double moment : direct.moments)
    {
        std::cout << ' ' << std::setprecision(17) << moment;
    }
    std::cout << " sigma22=" << direct.sigma22 << " tolerance=" << direct.tolerance << '\n';
    std::cout << "[failure branch] recursive moments:";
    for (const double moment : recursive.moments)
    {
        std::cout << ' ' << std::setprecision(17) << moment;
    }
    std::cout << " sigma22=" << recursive.sigma22 << " tolerance=" << recursive.tolerance << '\n';

    audit_failing_branch_01(*probe, *failed, wheeler_options, qmomkit::ConditionalMomentOptions{});
    run_expanded_layout_experiments(probe->dense_moments, layout, qmom_options);
}

void run_explicit_phase_space_evolution_parity()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2}, {0, 1, 2});
    const qmomkit::RawMomentTensor moments = raw_moments_from_terminals(deterministic_rectangular_nodes(), hierarchy.moment_extents());
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData lean_moments = qmomkit::extract_cqmom_moments(moments, layout);

    const auto velocity = [](double, const std::vector<double>& coordinate) { return std::vector<double>{0.03, -0.02, 0.04 + 0.01 * coordinate[0]}; };
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = true;
    qmom_options.rmin = 1.0e-12;
    qmom_options.eabs = 1.0e-12;
    const qmomkit::CqmomDenseCase dense_case{
        hierarchy, moments, {}, qmom_options, true, qmomkit::CqmomDiagnosticLevel::basic, qmomkit::VandermondeSolver::dense};
    const qmomkit::CqmomLeanCase lean_case{
        lean_moments, {}, qmom_options, true, qmomkit::CqmomDiagnosticLevel::basic, qmomkit::VandermondeSolver::dense};
    const qmomkit::CqmomState dense_state = qmomkit::internal::CqmomConstructionAccess::create_state(dense_case);
    const qmomkit::CqmomState lean_state = qmomkit::internal::CqmomConstructionAccess::create_state(lean_case);

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;

    const qmomkit::CqmomPhaseSpaceEvolutionResult dense = qmomkit::advance_cqmom_phase_space(dense_state, 0.0, 0.05, velocity, zero_source, options);
    const qmomkit::CqmomPhaseSpaceEvolutionResult lean = qmomkit::advance_cqmom_phase_space(lean_state, 0.0, 0.05, velocity, zero_source, options);

    std::cout << "[explicit phase-space evolution parity] dense_status=" << static_cast<int>(dense.status)
              << " lean_status=" << static_cast<int>(lean.status)
              << " remaining_limits=explicit active topology is one-call only; runtime storage selection is not exposed\n";
    QMOMKIT_VERIFY_CHECK(dense.status == lean.status);
    QMOMKIT_VERIFY_CHECK(dense.ok() == lean.ok());
    QMOMKIT_VERIFY_CHECK(dense.diagnostics.failed_stage == lean.diagnostics.failed_stage);
    QMOMKIT_VERIFY_CHECK(dense.diagnostics.cqmom_moment_storage == qmomkit::CqmomMomentStorage::dense);
    QMOMKIT_VERIFY_CHECK(lean.diagnostics.cqmom_moment_storage == qmomkit::CqmomMomentStorage::lean);
    QMOMKIT_VERIFY_CHECK(dense.diagnostics.cqmom_diagnostic_level == qmomkit::CqmomDiagnosticLevel::basic);
    QMOMKIT_VERIFY_CHECK(lean.diagnostics.cqmom_diagnostic_level == qmomkit::CqmomDiagnosticLevel::basic);
    QMOMKIT_VERIFY_CHECK(dense.diagnostics.cqmom_vandermonde_solver == qmomkit::VandermondeSolver::dense);
    QMOMKIT_VERIFY_CHECK(lean.diagnostics.cqmom_vandermonde_solver == qmomkit::VandermondeSolver::dense);
    QMOMKIT_VERIFY_CHECK(all_diagnostics_enable_adaptive_branch_counts(dense.diagnostics.initial_cqmom_diagnostics));
    QMOMKIT_VERIFY_CHECK(all_diagnostics_enable_adaptive_branch_counts(lean.diagnostics.initial_cqmom_diagnostics));
    QMOMKIT_VERIFY_CHECK(
        diagnostics_request_solver(dense.diagnostics.initial_cqmom_diagnostics, hierarchy.number_of_dimensions(), qmomkit::VandermondeSolver::dense));
    QMOMKIT_VERIFY_CHECK(
        diagnostics_request_solver(lean.diagnostics.initial_cqmom_diagnostics, hierarchy.number_of_dimensions(), qmomkit::VandermondeSolver::dense));
    compare_node_count_transitions(
        transitions_for_stage(dense.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::initial_inversion),
        transitions_for_stage(lean.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::initial_inversion));
    if (dense.ok() && lean.ok())
    {
        QMOMKIT_VERIFY_CHECK(dense.quadrature.has_value());
        QMOMKIT_VERIFY_CHECK(lean.quadrature.has_value());
        if (dense.quadrature.has_value() && lean.quadrature.has_value())
        {
            QMOMKIT_VERIFY_CHECK(dense.quadrature->adaptive_branch_counts_enabled());
            QMOMKIT_VERIFY_CHECK(lean.quadrature->adaptive_branch_counts_enabled());
        }
        const qmomkit::RawMomentTensor dense_moments{moments.extents, dense.state.moment_values()};
        const qmomkit::CqmomMomentData dense_common = qmomkit::extract_cqmom_moments(dense_moments, layout);
        QMOMKIT_VERIFY_CHECK(dense_common.values.size() == lean.state.moment_values().size());
        for (std::size_t i = 0; i < dense_common.values.size(); ++i)
        {
            QMOMKIT_VERIFY_ABSOLUTE_OR_RELATIVE(lean.state.moment_values()[i], dense_common.values[i], 5.0e-8, 5.0e-10);
        }
    }
}

} // namespace

int main()
{
    const qmomkit::PhaseSpaceHierarchy identity({2, 2, 2}, {0, 1, 2});
    const qmomkit::PhaseSpaceHierarchy order_210({2, 2, 2}, {2, 1, 0});

    run_static_parity_case("fixed 3D 2,2,2 identity order", identity, deterministic_rectangular_nodes(), false);
    run_static_parity_case("fixed 3D 2,2,2 reverse order", order_210, independent_rectangular_nodes(), false);
    run_static_parity_case("adaptive branch-specific 3D 2,2,2 identity order", identity, adaptive_rectangular_nodes(), true);
    run_static_parity_case("near-degenerate conditional u branch", identity, near_degenerate_u_nodes(), false);
    run_physical_analytical_contraction_case(identity);
    run_failure_probe();
    run_explicit_phase_space_evolution_parity();

    std::cout << "[dense/lean parity] Stage 1 diagnostic verification completed.\n";
    return 0;
}
