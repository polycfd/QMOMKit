#include "internal/cqmom_case.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

const char* fixed_cqmom_case_text = R"(
[phase_space]
names = x, y
node_counts = 2, 2

[moments]
extents = 4, 4
M[0,0] = 1.0
M[0,1] = 0.90000000000000013
M[0,2] = 1.6499999999999999
M[0,3] = 2.3250000000000002
M[1,0] = 0.5
M[1,1] = 0.75000000000000022
M[1,2] = 1.125
M[1,3] = 1.6875000000000002
M[2,0] = 2.5
M[2,1] = 2.5500000000000003
M[2,2] = 4.4249999999999998
M[2,3] = 6.3375000000000004
M[3,0] = 3.5
M[3,1] = 4.0500000000000007
M[3,2] = 6.6750000000000007
M[3,3] = 9.7125000000000004
)";

template <typename T, typename = void> struct HasPublicMomentsMember : std::false_type
{
};

template <typename T> struct HasPublicMomentsMember<T, std::void_t<decltype(std::declval<T&>().moments)>> : std::true_type
{
};

void assert_near(const double actual, const double expected, const double tolerance)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void assert_vector_near(const std::vector<double>& actual, const std::vector<double>& expected, const double tolerance)
{
    assert(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        assert_near(actual[i], expected[i], tolerance);
    }
}

void assert_cqmom_quadrature_near(const qmomkit::CqmomQuadrature& actual, const qmomkit::CqmomQuadrature& expected, const double tolerance)
{
    assert(actual.hierarchy().node_counts() == expected.hierarchy().node_counts());
    assert(actual.hierarchy().conditioning_order() == expected.hierarchy().conditioning_order());
    assert(actual.total_nodes() == expected.total_nodes());
    assert(actual.has_rectangular_storage() == expected.has_rectangular_storage());
    assert(actual.terminal_paths() == expected.terminal_paths());
    assert_vector_near(actual.weights(), expected.weights(), tolerance);
    assert(actual.abscissas().size() == expected.abscissas().size());
    for (std::size_t dimension = 0; dimension < actual.abscissas().size(); ++dimension)
    {
        assert_vector_near(actual.abscissas()[dimension], expected.abscissas()[dimension], tolerance);
    }
}

std::vector<double> raw_moments_1d(const std::vector<double>& abscissas, const std::vector<double>& weights, const std::size_t moment_count)
{
    std::vector<double> moments(moment_count, 0.0);
    for (std::size_t order = 0; order < moment_count; ++order)
    {
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(abscissas[node], static_cast<int>(order));
        }
    }
    return moments;
}

std::vector<double> raw_moments_from_terminals(const std::vector<std::size_t>& extents, const std::vector<double>& weights,
                                               const std::vector<std::vector<double>>& abscissas)
{
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
    return moments;
}

qmomkit::RawMomentTensor fixed_cqmom_moments()
{
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 2.0, 2.0}, {-0.5, 1.5, -0.5, 1.5}};
    return {extents, raw_moments_from_terminals(extents, weights, abscissas)};
}

qmomkit::RawMomentTensor adaptive_cqmom_moments()
{
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.1, 0.3, 0.6};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0}, {-2.0, 2.0, 3.0}};
    return {extents, raw_moments_from_terminals(extents, weights, abscissas)};
}

qmomkit::RawMomentTensor branch_specific_ragged_cqmom_moments()
{
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.5};
    const std::vector<std::vector<double>> abscissas{{-1.0, 1.0, 1.0}, {0.0, -0.5, 0.5}};
    return {extents, raw_moments_from_terminals(extents, weights, abscissas)};
}

qmomkit::RawMomentTensor collapsed_branch_cqmom_moments()
{
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.8};
    const std::vector<std::vector<double>> abscissas{{-1.0, 1.0}, {0.0, 0.125}};
    return {extents, raw_moments_from_terminals(extents, weights, abscissas)};
}

qmomkit::RawMomentTensor root_collapsed_child_spread_cqmom_moments()
{
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.5, 0.5};
    const std::vector<std::vector<double>> abscissas{{0.0, 0.0}, {-1.0, 1.0}};
    return {extents, raw_moments_from_terminals(extents, weights, abscissas)};
}

qmomkit::RawMomentTensor second_level_split_3d_cqmom_moments()
{
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.25, 0.25, 0.5};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0}, {-1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}};
    return {extents, raw_moments_from_terminals(extents, weights, abscissas)};
}

qmomkit::RawMomentTensor second_level_collapsed_grandchild_spread_3d_cqmom_moments()
{
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.25, 0.25, 0.5};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, 1.0}, {0.0, 0.0, 0.0}, {-1.0, 1.0, 0.0}};
    return {extents, raw_moments_from_terminals(extents, weights, abscissas)};
}

qmomkit::RawMomentTensor nonidentity_ragged_cqmom_moments()
{
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.5};
    const std::vector<std::vector<double>> abscissas{{-1.0, 1.0, 1.0}, {-1.0, -1.0, 1.0}};
    return {extents, raw_moments_from_terminals(extents, weights, abscissas)};
}

qmomkit::NodeAdvectionOptions tight_node_options()
{
    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    options.rk23_options.relative_tolerance = 1.0e-9;
    options.rk23_options.absolute_tolerance = 1.0e-12;
    return options;
}

qmomkit::CqmomLeanTransportState lean_state_from_dense(const qmomkit::RawMomentTensor& moments, const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    return qmomkit::CqmomLeanTransportState{qmomkit::extract_cqmom_moments(moments, layout)};
}

qmomkit::CqmomCase fixed_cqmom_case()
{
    const auto parsed = qmomkit::config::parse_cqmom(fixed_cqmom_case_text, "fixed-cqmom-case");
    assert(parsed.ok());
    return *parsed.value;
}

qmomkit::CqmomCase fixed_cqmom_case_with_vandermonde_solver(const char* const solver)
{
    std::string text = fixed_cqmom_case_text;
    text += "\n[vandermonde]\nsolver = ";
    text += solver;
    text += "\n";
    const auto parsed = qmomkit::config::parse_cqmom(text, "fixed-cqmom-case-with-vandermonde-solver");
    assert(parsed.ok());
    return *parsed.value;
}

std::string cqmom_case_text_from_moments_with_hierarchy(const qmomkit::RawMomentTensor& moments, const std::vector<std::size_t>& node_counts,
                                                        const std::vector<std::size_t>& conditioning_order, const bool adaptive_branch_counts,
                                                        const bool adaptive_wheeler)
{
    std::ostringstream text;
    text.precision(17);
    text << "[phase_space]\n";
    text << "names = ";
    for (std::size_t dimension = 0; dimension < node_counts.size(); ++dimension)
    {
        if (dimension > 0)
        {
            text << ", ";
        }
        text << "x" << dimension;
    }
    text << "\nnode_counts = ";
    for (std::size_t dimension = 0; dimension < node_counts.size(); ++dimension)
    {
        if (dimension > 0)
        {
            text << ", ";
        }
        text << node_counts[dimension];
    }
    text << "\nconditioning_order = ";
    for (std::size_t level = 0; level < conditioning_order.size(); ++level)
    {
        if (level > 0)
        {
            text << ", ";
        }
        text << conditioning_order[level];
    }
    text << "\n";
    text << "adaptive_branch_counts = " << (adaptive_branch_counts ? "true" : "false") << "\n\n";
    text << "[qmom]\n";
    text << "adaptive = " << (adaptive_wheeler ? "true" : "false") << "\n\n";
    text << "[moments]\n";
    text << "extents = ";
    for (std::size_t dimension = 0; dimension < moments.extents.size(); ++dimension)
    {
        if (dimension > 0)
        {
            text << ", ";
        }
        text << moments.extents[dimension];
    }
    text << "\n";
    for (std::size_t flat = 0; flat < moments.values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, moments.extents);
        text << "M[";
        for (std::size_t dimension = 0; dimension < index.size(); ++dimension)
        {
            if (dimension > 0)
            {
                text << ",";
            }
            text << index[dimension];
        }
        text << "] = " << moments.values[flat] << "\n";
    }
    return text.str();
}

std::string cqmom_case_text_from_moments(const qmomkit::RawMomentTensor& moments, const bool adaptive_branch_counts, const bool adaptive_wheeler)
{
    return cqmom_case_text_from_moments_with_hierarchy(moments, {2, 2}, {0, 1}, adaptive_branch_counts, adaptive_wheeler);
}

qmomkit::CqmomCase adaptive_cqmom_case()
{
    const std::string text = cqmom_case_text_from_moments(adaptive_cqmom_moments(), true, true);
    const auto parsed = qmomkit::config::parse_cqmom(text, "adaptive-cqmom-case");
    assert(parsed.ok());
    return *parsed.value;
}

qmomkit::CqmomCase adaptive_cqmom_case_from_moments(const qmomkit::RawMomentTensor& moments)
{
    const std::string text = cqmom_case_text_from_moments(moments, true, true);
    const auto parsed = qmomkit::config::parse_cqmom(text, "adaptive-cqmom-case-from-moments");
    assert(parsed.ok());
    return *parsed.value;
}

qmomkit::CqmomCase adaptive_cqmom_case_from_moments_with_hierarchy(const qmomkit::RawMomentTensor& moments,
                                                                   const std::vector<std::size_t>& node_counts,
                                                                   const std::vector<std::size_t>& conditioning_order)
{
    const std::string text = cqmom_case_text_from_moments_with_hierarchy(moments, node_counts, conditioning_order, true, true);
    const auto parsed = qmomkit::config::parse_cqmom(text, "adaptive-cqmom-case-from-moments-with-hierarchy");
    assert(parsed.ok());
    return *parsed.value;
}

bool diagnostics_contain_reduced_branch(const std::vector<qmomkit::CqmomDiagnostics>& diagnostics)
{
    for (const qmomkit::CqmomDiagnostics& item : diagnostics)
    {
        if (item.adaptive_branch_counts_enabled && item.adaptive_node_count_reduced && item.actual_node_count < item.requested_node_count)
        {
            return true;
        }
    }
    return false;
}

bool diagnostics_all_use_requested_node_counts(const std::vector<qmomkit::CqmomDiagnostics>& diagnostics)
{
    bool saw_branch = false;
    for (const qmomkit::CqmomDiagnostics& item : diagnostics)
    {
        if (!item.adaptive_branch_counts_enabled)
        {
            return false;
        }
        saw_branch = true;
        if (item.actual_node_count != item.requested_node_count || item.adaptive_node_count_reduced)
        {
            return false;
        }
    }
    return saw_branch;
}

bool active_counts_contain(const std::vector<qmomkit::CqmomActiveBranchCount>& counts, const std::size_t level, const qmomkit::ConditioningPath& path,
                           const std::size_t active_count, const std::optional<std::size_t> phase_space_dimension = std::nullopt,
                           const std::optional<std::size_t> configured_node_count = std::nullopt)
{
    for (const qmomkit::CqmomActiveBranchCount& count : counts)
    {
        if (count.conditioning_level == level && count.parent_path == path && count.active_node_count == active_count &&
            (!phase_space_dimension.has_value() || count.phase_space_dimension == *phase_space_dimension) &&
            (!configured_node_count.has_value() || count.configured_node_count == *configured_node_count))
        {
            return true;
        }
    }
    return false;
}

bool transitions_contain(const std::vector<qmomkit::CqmomNodeCountTransition>& transitions, const qmomkit::TransportStage stage,
                         const std::size_t level, const qmomkit::ConditioningPath& path, const std::optional<std::size_t> before,
                         const std::optional<std::size_t> after, const qmomkit::CqmomNodeCountEvolution policy,
                         const qmomkit::CqmomNodeCountTransitionKind kind)
{
    for (const qmomkit::CqmomNodeCountTransition& transition : transitions)
    {
        if (transition.stage == stage && transition.conditioning_level == level && transition.parent_path == path &&
            transition.active_node_count_before == before && transition.active_node_count_after == after &&
            transition.node_count_evolution == policy && transition.kind == kind)
        {
            return true;
        }
    }
    return false;
}

bool active_counts_equal(const std::vector<qmomkit::CqmomActiveBranchCount>& lhs, const std::vector<qmomkit::CqmomActiveBranchCount>& rhs)
{
    if (lhs.size() != rhs.size())
    {
        return false;
    }
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        if (lhs[i].conditioning_level != rhs[i].conditioning_level || lhs[i].parent_path != rhs[i].parent_path ||
            lhs[i].phase_space_dimension != rhs[i].phase_space_dimension || lhs[i].configured_node_count != rhs[i].configured_node_count ||
            lhs[i].active_node_count != rhs[i].active_node_count)
        {
            return false;
        }
    }
    return true;
}

qmomkit::CqmomMomentSource source_from_to(const qmomkit::CqmomMomentData& initial, const qmomkit::CqmomMomentData& target)
{
    return [initial, target](double, const qmomkit::CqmomSourceContext& moments)
    {
        assert(moments.size() == target.values.size());
        std::vector<double> rate(moments.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            rate[i] = target.values[i] - initial.values[i];
        }
        return qmomkit::CqmomMomentRate{rate};
    };
}

void compare_lean_transport_to_dense_layout_closed_source(const qmomkit::TransportScheme scheme)
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::CqmomLeanTransportState state = lean_state_from_dense(moments, hierarchy);
    const auto velocity = [](double, const std::vector<double>& xi)
    {
        assert(xi.size() == 2);
        return std::vector<double>{0.2 * xi[0] + 0.1, -0.15 * xi[1] + 0.05};
    };
    constexpr double lambda = -0.35;
    const auto dense_source = [lambda](double, const std::vector<double>& state_values)
    {
        std::vector<double> rate(state_values.size(), 0.0);
        for (std::size_t i = 0; i < state_values.size(); ++i)
        {
            rate[i] = lambda * state_values[i];
        }
        return rate;
    };
    const auto lean_source = [lambda](double, const qmomkit::CqmomSourceContext& state_values)
    {
        std::vector<double> rate(state_values.size(), 0.0);
        for (std::size_t i = 0; i < state_values.size(); ++i)
        {
            rate[i] = lambda * state_values.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };

    qmomkit::CqmomDenseTransportOptions dense_options;
    dense_options.node_options = tight_node_options();
    dense_options.scheme = scheme;
    const qmomkit::CqmomDenseTransportResult dense =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.2, velocity, dense_source, dense_options);

    qmomkit::CqmomLeanTransportOptions lean_options;
    lean_options.node_options = tight_node_options();
    lean_options.scheme = scheme;
    const qmomkit::CqmomLeanTransportResult lean = qmomkit::advance_transport_cqmom(state, 0.0, 0.2, velocity, lean_source, lean_options);

    assert(dense.ok());
    assert(lean.ok());
    assert(lean.quadrature.has_value());
    assert(lean.quadrature->has_rectangular_storage());
    assert(lean.diagnostics.scheme == scheme);

    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments(dense.moments, state.moments.layout);
    assert_vector_near(lean.state.moments.values, expected.values, 5.0e-9);
}

void test_continuous_only_matches_direct_operator_qmom()
{
    const std::vector<double> moments = raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4);
    const qmomkit::QmomResult initial = qmomkit::compute_qmom(qmomkit::MomentSet(moments), 2);
    assert(initial.ok());

    const auto velocity = [](double, const std::vector<double>& xi) { return std::vector<double>{0.5 * xi[0]}; };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };

    qmomkit::QmomTransportOptions options;
    options.node_options = tight_node_options();
    options.scheme = qmomkit::TransportScheme::strang;
    const qmomkit::QmomTransportResult transport = qmomkit::advance_transport_qmom(moments, 2, 0.0, 0.3, velocity, zero_source, options);
    const qmomkit::QmomNodeAdvectionResult direct = qmomkit::advance_nodes(initial, 0.0, 0.3, velocity, options.node_options);

    assert(transport.ok());
    assert(direct.ok());
    const std::vector<double> expected = qmomkit::reconstruct_raw_moments(direct.quadrature, moments.size());
    assert_vector_near(transport.moments, expected, 5.0e-9);
}

void test_source_only_matches_direct_operator_qmom()
{
    const std::vector<double> moments{1.0, 2.0};
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto source = [](double, const std::vector<double>& state)
    {
        assert(state.size() == 2);
        return std::vector<double>{0.0, 3.0};
    };

    const qmomkit::QmomTransportResult transport = qmomkit::advance_transport_qmom(moments, 1, 0.0, 0.25, zero_velocity, source);
    const qmomkit::MomentSourceUpdateResult direct = qmomkit::advance_moments_with_source(moments, 0.125, 0.25, source);

    assert(transport.ok());
    assert(direct.ok());
    assert_vector_near(transport.moments, direct.moments, 1.0e-13);
}

void test_source_plus_boundary_matches_direct_operator_qmom()
{
    const std::vector<double> moments{1.0, 2.0};
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::MomentSource source = [](double, const std::vector<double>& state)
    {
        assert(state.size() == 2);
        return std::vector<double>{0.0, 1.0};
    };
    const qmomkit::MomentBoundaryContribution boundary = [](double, const std::vector<double>& state)
    {
        assert(state.size() == 2);
        return qmomkit::MomentBoundaryRate{std::vector<double>{0.25, -0.5}};
    };

    const qmomkit::QmomTransportResult transport = qmomkit::advance_transport_qmom(moments, 1, 0.0, 0.25, zero_velocity, source, boundary, {});
    const qmomkit::MomentSourceUpdateResult direct = qmomkit::advance_moments_with_source(moments, 0.125, 0.25, source, boundary, {});

    assert(transport.ok());
    assert(direct.ok());
    assert_vector_near(transport.moments, direct.moments, 1.0e-13);
}

void test_adaptive_branch_specific_cqmom_case()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = adaptive_cqmom_moments();
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };

    qmomkit::CqmomDenseTransportOptions options;
    options.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    options.adaptive_branch_counts = true;
    const qmomkit::CqmomDenseTransportResult result =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(result.ok());
    assert(result.quadrature.has_value());
    assert(result.quadrature->adaptive_branch_counts_enabled());
    assert(!result.quadrature->has_rectangular_storage());
    assert(result.quadrature->total_nodes() == 3);
    assert_vector_near(result.moments.values, moments.values, 1.0e-11);
}

void test_adaptive_state_initialization_preserves_policy()
{
    const qmomkit::CqmomCase cqmom_case = adaptive_cqmom_case();
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_case);
    const qmomkit::CqmomState state = cqmom_case.initial_state();

    assert(cqmom_case.adaptive_branch_counts());
    assert(cqmom_case.qmom_options().adaptive);
    assert(state.adaptive_branch_counts());
    assert(state.qmom_options().adaptive);
    assert(initial.ok());
    assert(initial.quadrature.adaptive_branch_counts_enabled());
    assert(!initial.quadrature.has_rectangular_storage());
    assert(initial.quadrature.total_nodes() == 3);
    assert(diagnostics_contain_reduced_branch(initial.diagnostics));
    assert(!state.active_branch_counts().empty());
    assert(active_counts_contain(state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, 2));
    assert(active_counts_contain(state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 2, 1, 2));
    assert(active_counts_contain(state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 1, 1, 2));
}

void test_adaptive_state_continuous_transport_matches_dense()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = adaptive_cqmom_moments();
    const qmomkit::CqmomState state = adaptive_cqmom_case().initial_state();
    const auto velocity = [](double, const std::vector<double>& xi)
    {
        assert(xi.size() == 2);
        return std::vector<double>{0.15 * xi[0] + 0.05, -0.2 * xi[1] + 0.1};
    };
    const auto dense_zero_source = [](double, const std::vector<double>& state_values) { return std::vector<double>(state_values.size(), 0.0); };
    const auto lean_zero_source = [](double, const qmomkit::CqmomSourceContext& state_values)
    { return qmomkit::CqmomMomentRate{std::vector<double>(state_values.size(), 0.0)}; };

    qmomkit::CqmomDenseTransportOptions dense_options;
    dense_options.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    dense_options.adaptive_branch_counts = true;
    dense_options.node_options = tight_node_options();
    dense_options.scheme = qmomkit::TransportScheme::cd;
    const qmomkit::CqmomDenseTransportResult dense =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.2, velocity, dense_zero_source, dense_options);

    qmomkit::CqmomStateTransportOptions state_options;
    state_options.node_options = dense_options.node_options;
    state_options.scheme = dense_options.scheme;
    const qmomkit::CqmomStateTransportResult transported =
        qmomkit::advance_transport_cqmom(state, 0.0, 0.2, velocity, lean_zero_source, state_options);

    assert(dense.ok());
    assert(transported.ok());
    assert(dense.quadrature.has_value());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->adaptive_branch_counts_enabled());
    assert(!transported.quadrature->has_rectangular_storage());
    assert_cqmom_quadrature_near(*transported.quadrature, *dense.quadrature, 5.0e-9);
    assert(transported.diagnostics.first_reconstruction_diagnostics.success);
}

void test_adaptive_state_source_transport_changes_active_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = adaptive_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = fixed_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const qmomkit::CqmomState state = adaptive_cqmom_case().initial_state();
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto source = [initial_required, target_required](double, const qmomkit::CqmomSourceContext& moments)
    {
        assert(moments.size() == target_required.values.size());
        std::vector<double> rate(moments.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            rate[i] = target_required.values[i] - initial_required.values[i];
        }
        return qmomkit::CqmomMomentRate{rate};
    };

    qmomkit::CqmomStateTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    const qmomkit::CqmomStateTransportResult transported = qmomkit::advance_transport_cqmom(state, 0.0, 1.0, zero_velocity, source, options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->adaptive_branch_counts_enabled());
    assert(!transported.quadrature->has_rectangular_storage());
    assert(transported.quadrature->total_nodes() == 4);
    assert(diagnostics_all_use_requested_node_counts(transported.diagnostics.post_source_cqmom_diagnostics));
    assert(transported.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(transported.diagnostics.post_source_cqmom_status == qmomkit::CqmomStatus::success);
}

void test_adaptive_state_preserves_branch_counts_across_steps()
{
    const qmomkit::CqmomCase case_ = adaptive_cqmom_case_from_moments(branch_specific_ragged_cqmom_moments());
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomStateTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.source_node_count_evolution = qmomkit::CqmomNodeCountEvolution::preserve_active_counts;
    qmomkit::CqmomStateTransportResult first = qmomkit::advance_transport_cqmom(case_.initial_state(), 0.0, 0.1, zero_velocity, zero_source, options);
    assert(first.ok());
    assert(first.quadrature.has_value());
    assert(first.quadrature->total_nodes() == 3);
    assert(active_counts_contain(first.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1));
    assert(active_counts_contain(first.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2));
    assert(first.diagnostics.cqmom_node_count_transitions.empty());

    const qmomkit::CqmomStateTransportResult second =
        qmomkit::advance_transport_cqmom(first.state, first.final_time, 0.1, zero_velocity, zero_source, options);
    assert(second.ok());
    assert(second.quadrature.has_value());
    assert(second.quadrature->terminal_paths() == first.quadrature->terminal_paths());
    assert(active_counts_contain(second.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1));
    assert(active_counts_contain(second.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2));
}

void test_adaptive_state_preserved_source_can_reduce_active_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = collapsed_branch_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    qmomkit::CqmomStateTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.source_node_count_evolution = qmomkit::CqmomNodeCountEvolution::preserve_active_counts;
    const qmomkit::CqmomStateTransportResult transported =
        qmomkit::advance_transport_cqmom(adaptive_cqmom_case_from_moments(initial_moments).initial_state(), 0.0, 1.0, zero_velocity,
                                         source_from_to(initial_required, target_required), options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->total_nodes() == 2);
    assert(active_counts_contain(transported.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 1));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{1}}, 2, 1, qmomkit::CqmomNodeCountEvolution::preserve_active_counts,
                               qmomkit::CqmomNodeCountTransitionKind::reduced));
}

void test_adaptive_state_source_policy_allows_growth_to_configured_maximum()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = collapsed_branch_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    qmomkit::CqmomStateTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    assert(options.source_node_count_evolution == qmomkit::CqmomNodeCountEvolution::allow_growth_to_configured_maximum);
    const qmomkit::CqmomStateTransportResult transported =
        qmomkit::advance_transport_cqmom(adaptive_cqmom_case_from_moments(initial_moments).initial_state(), 0.0, 1.0, zero_velocity,
                                         source_from_to(initial_required, target_required), options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->total_nodes() == 3);
    assert(active_counts_contain(transported.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{1}}, 1, 2, qmomkit::CqmomNodeCountEvolution::allow_growth_to_configured_maximum,
                               qmomkit::CqmomNodeCountTransitionKind::grew));
}

void test_adaptive_state_preserved_source_invalidates_descendants_after_root_collapse()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = root_collapsed_child_spread_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    qmomkit::CqmomStateTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.source_node_count_evolution = qmomkit::CqmomNodeCountEvolution::preserve_active_counts;
    const qmomkit::CqmomStateTransportResult transported =
        qmomkit::advance_transport_cqmom(adaptive_cqmom_case_from_moments(initial_moments).initial_state(), 0.0, 1.0, zero_velocity,
                                         source_from_to(initial_required, target_required), options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->total_nodes() == 2);
    assert(active_counts_contain(transported.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 1));
    assert(active_counts_contain(transported.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 2));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 0,
                               qmomkit::ConditioningPath{}, 2, 1, qmomkit::CqmomNodeCountEvolution::preserve_active_counts,
                               qmomkit::CqmomNodeCountTransitionKind::reduced));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{0}}, 1, std::nullopt, qmomkit::CqmomNodeCountEvolution::preserve_active_counts,
                               qmomkit::CqmomNodeCountTransitionKind::removed));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{0}}, std::nullopt, 2, qmomkit::CqmomNodeCountEvolution::preserve_active_counts,
                               qmomkit::CqmomNodeCountTransitionKind::created));
}

void test_adaptive_state_preserved_source_invalidates_grandchildren_after_second_level_collapse()
{
    const std::vector<std::size_t> node_counts{2, 2, 2};
    const std::vector<std::size_t> conditioning_order{0, 1, 2};
    const qmomkit::PhaseSpaceHierarchy hierarchy(node_counts, conditioning_order);
    const qmomkit::RawMomentTensor initial_moments = second_level_split_3d_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = second_level_collapsed_grandchild_spread_3d_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    qmomkit::CqmomStateTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.source_node_count_evolution = qmomkit::CqmomNodeCountEvolution::preserve_active_counts;
    const qmomkit::CqmomStateTransportResult transported = qmomkit::advance_transport_cqmom(
        adaptive_cqmom_case_from_moments_with_hierarchy(initial_moments, node_counts, conditioning_order).initial_state(), 0.0, 1.0, zero_velocity,
        source_from_to(initial_required, target_required), options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->total_nodes() == 3);
    assert(active_counts_contain(transported.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2));
    assert(active_counts_contain(transported.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1));
    assert(active_counts_contain(transported.state.active_branch_counts(), 2, qmomkit::ConditioningPath{{0, 0}}, 2));
    assert(active_counts_contain(transported.state.active_branch_counts(), 2, qmomkit::ConditioningPath{{1, 0}}, 1));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{0}}, 2, 1, qmomkit::CqmomNodeCountEvolution::preserve_active_counts,
                               qmomkit::CqmomNodeCountTransitionKind::reduced));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 2,
                               qmomkit::ConditioningPath{{0, 0}}, 1, std::nullopt, qmomkit::CqmomNodeCountEvolution::preserve_active_counts,
                               qmomkit::CqmomNodeCountTransitionKind::removed));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 2,
                               qmomkit::ConditioningPath{{0, 0}}, std::nullopt, 2, qmomkit::CqmomNodeCountEvolution::preserve_active_counts,
                               qmomkit::CqmomNodeCountTransitionKind::created));
}

void test_adaptive_state_initial_counts_preserve_nonidentity_conditioning_order()
{
    const std::vector<std::size_t> node_counts{2, 2};
    const std::vector<std::size_t> conditioning_order{1, 0};
    const qmomkit::RawMomentTensor moments = nonidentity_ragged_cqmom_moments();
    const qmomkit::CqmomState state = adaptive_cqmom_case_from_moments_with_hierarchy(moments, node_counts, conditioning_order).initial_state();

    assert(state.adaptive_branch_counts());
    assert(active_counts_contain(state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 1, 2));
    assert(active_counts_contain(state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 2, 0, 2));
    assert(active_counts_contain(state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 1, 0, 2));
}

void test_fixed_state_transport_does_not_store_adaptive_branch_counts()
{
    const qmomkit::CqmomState state = fixed_cqmom_case().initial_state();
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomStateTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.source_node_count_evolution = qmomkit::CqmomNodeCountEvolution::preserve_active_counts;
    const qmomkit::CqmomStateTransportResult transported = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->has_rectangular_storage());
    assert(transported.quadrature->total_nodes() == 4);
    assert(transported.state.active_branch_counts().empty());
    assert(transported.diagnostics.cqmom_node_count_transitions.empty());
}

void test_adaptive_state_rejects_duplicate_active_count_metadata()
{
    const qmomkit::CqmomState good_state = adaptive_cqmom_case().initial_state();
    std::vector<qmomkit::CqmomActiveBranchCount> bad_counts = good_state.active_branch_counts();
    assert(!bad_counts.empty());
    bad_counts.push_back(bad_counts.front());
    const qmomkit::CqmomState bad_state =
        qmomkit::internal::create_cqmom_state(qmomkit::internal::cqmom_state_case(good_state), std::move(bad_counts));
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomStateTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.source_node_count_evolution = qmomkit::CqmomNodeCountEvolution::preserve_active_counts;
    const qmomkit::CqmomStateTransportResult result = qmomkit::advance_transport_cqmom(bad_state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::initial_inversion);
    assert(result.diagnostics.initial_cqmom_status == qmomkit::CqmomStatus::invalid_input);
}

void test_adaptive_state_rejects_missing_active_count_metadata()
{
    const qmomkit::CqmomState good_state = adaptive_cqmom_case().initial_state();
    std::vector<qmomkit::CqmomActiveBranchCount> bad_counts = good_state.active_branch_counts();
    assert(!bad_counts.empty());
    bad_counts.pop_back();
    const qmomkit::CqmomState bad_state =
        qmomkit::internal::create_cqmom_state(qmomkit::internal::cqmom_state_case(good_state), std::move(bad_counts));
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomStateTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.source_node_count_evolution = qmomkit::CqmomNodeCountEvolution::preserve_active_counts;
    const qmomkit::CqmomStateTransportResult result = qmomkit::advance_transport_cqmom(bad_state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::initial_inversion);
    assert(result.diagnostics.initial_cqmom_status == qmomkit::CqmomStatus::invalid_input);
}

void test_adaptive_state_failure_preserves_active_counts_for_retry()
{
    const qmomkit::CqmomState state = adaptive_cqmom_case().initial_state();
    const std::vector<qmomkit::CqmomActiveBranchCount> original_counts = state.active_branch_counts();
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto bad_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomStateTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.source_node_count_evolution = qmomkit::CqmomNodeCountEvolution::preserve_active_counts;
    const qmomkit::CqmomStateTransportResult failed = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, bad_source, options);
    assert(!failed.ok());
    assert(active_counts_equal(state.active_branch_counts(), original_counts));
    assert(active_counts_equal(failed.state.active_branch_counts(), original_counts));

    const qmomkit::CqmomStateTransportResult retry = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, zero_source, options);
    assert(retry.ok());
    assert(active_counts_equal(retry.state.active_branch_counts(), original_counts));
}

void test_adaptive_state_post_source_failure_preserves_active_counts_for_retry()
{
    const qmomkit::CqmomState state = adaptive_cqmom_case().initial_state();
    const std::vector<qmomkit::CqmomActiveBranchCount> original_counts = state.active_branch_counts();
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto bad_realizability_source = [](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        rate.front() = -20.0;
        return qmomkit::CqmomMomentRate{rate};
    };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomStateTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.source_node_count_evolution = qmomkit::CqmomNodeCountEvolution::preserve_active_counts;
    const qmomkit::CqmomStateTransportResult failed =
        qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, bad_realizability_source, options);
    assert(!failed.ok());
    assert(failed.status == qmomkit::TransportStatus::inversion_failed);
    assert(active_counts_equal(state.active_branch_counts(), original_counts));
    assert(active_counts_equal(failed.state.active_branch_counts(), original_counts));

    const qmomkit::CqmomStateTransportResult retry = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, zero_source, options);
    assert(retry.ok());
    assert(active_counts_equal(retry.state.active_branch_counts(), original_counts));
}

void test_adaptive_state_strang_transport_matches_dense()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = adaptive_cqmom_moments();
    const qmomkit::CqmomState state = adaptive_cqmom_case().initial_state();
    const auto velocity = [](double, const std::vector<double>& xi)
    {
        assert(xi.size() == 2);
        return std::vector<double>{0.1 + 0.05 * xi[0], -0.08 * xi[1]};
    };
    const auto dense_zero_source = [](double, const std::vector<double>& state_values) { return std::vector<double>(state_values.size(), 0.0); };
    const auto lean_zero_source = [](double, const qmomkit::CqmomSourceContext& state_values)
    { return qmomkit::CqmomMomentRate{std::vector<double>(state_values.size(), 0.0)}; };

    qmomkit::CqmomDenseTransportOptions dense_options;
    dense_options.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    dense_options.adaptive_branch_counts = true;
    dense_options.node_options = tight_node_options();
    dense_options.scheme = qmomkit::TransportScheme::strang;
    const qmomkit::CqmomDenseTransportResult dense =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.2, velocity, dense_zero_source, dense_options);

    qmomkit::CqmomStateTransportOptions state_options;
    state_options.node_options = dense_options.node_options;
    state_options.scheme = dense_options.scheme;
    const qmomkit::CqmomStateTransportResult transported =
        qmomkit::advance_transport_cqmom(state, 0.0, 0.2, velocity, lean_zero_source, state_options);

    assert(dense.ok());
    assert(transported.ok());
    assert(dense.quadrature.has_value());
    assert(transported.quadrature.has_value());
    assert_cqmom_quadrature_near(*transported.quadrature, *dense.quadrature, 5.0e-9);
    assert(transported.diagnostics.first_reconstruction_diagnostics.success);
    assert(transported.diagnostics.final_reconstruction_diagnostics.success);
}

void test_adaptive_state_transport_reports_initial_failure()
{
    qmomkit::RawMomentTensor moments = adaptive_cqmom_moments();
    std::fill(moments.values.begin(), moments.values.end(), 0.0);
    const std::string text = cqmom_case_text_from_moments(moments, true, true);
    const auto parsed = qmomkit::config::parse_cqmom(text, "adaptive-failing-cqmom-case");
    assert(parsed.ok());
    bool threw = false;
    try
    {
        static_cast<void>(parsed.value->initial_state());
    }
    catch (const std::runtime_error& error)
    {
        threw = true;
        assert(std::string(error.what()).find("active-count seeding failed") != std::string::npos);
    }
    assert(threw);
}

void test_lean_cqmom_strang_matches_dense_layout_closed_source()
{
    compare_lean_transport_to_dense_layout_closed_source(qmomkit::TransportScheme::strang);
}

void test_lean_cqmom_cd_matches_dense_layout_closed_source()
{
    compare_lean_transport_to_dense_layout_closed_source(qmomkit::TransportScheme::cd);
}

void test_lean_cqmom_dc_matches_dense_layout_closed_source()
{
    compare_lean_transport_to_dense_layout_closed_source(qmomkit::TransportScheme::dc);
}

void test_state_cqmom_transport_matches_lean_fixed_layout()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::CqmomLeanTransportState lean_state = lean_state_from_dense(moments, hierarchy);
    const qmomkit::CqmomState state = fixed_cqmom_case().initial_state();
    const auto velocity = [](double, const std::vector<double>& xi)
    {
        assert(xi.size() == 2);
        return std::vector<double>{0.2 * xi[0] + 0.1, -0.15 * xi[1] + 0.05};
    };
    constexpr double lambda = -0.35;
    const auto source = [lambda](double, const qmomkit::CqmomSourceContext& state_values)
    {
        std::vector<double> rate(state_values.size(), 0.0);
        for (std::size_t i = 0; i < state_values.size(); ++i)
        {
            rate[i] = lambda * state_values.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };

    qmomkit::CqmomStateTransportOptions state_options;
    state_options.node_options = tight_node_options();
    state_options.scheme = qmomkit::TransportScheme::strang;
    const qmomkit::CqmomStateTransportResult state_result = qmomkit::advance_transport_cqmom(state, 0.0, 0.2, velocity, source, state_options);

    qmomkit::CqmomLeanTransportOptions lean_options;
    lean_options.node_options = state_options.node_options;
    lean_options.scheme = state_options.scheme;
    const qmomkit::CqmomLeanTransportResult lean_result = qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.2, velocity, source, lean_options);

    assert(state_result.ok());
    assert(lean_result.ok());
    assert(state_result.quadrature.has_value());
    assert(lean_result.quadrature.has_value());
    assert_cqmom_quadrature_near(*state_result.quadrature, *lean_result.quadrature, 5.0e-9);
    assert(state_result.diagnostics.scheme == lean_result.diagnostics.scheme);
    assert(state_result.final_time == lean_result.final_time);
    assert(state_result.state.hierarchy().node_counts() == state.hierarchy().node_counts());
    assert(state_result.state.hierarchy().conditioning_order() == state.hierarchy().conditioning_order());
}

void test_state_cqmom_transport_preserves_vandermonde_solver_policy()
{
    const qmomkit::CqmomCase case_ = fixed_cqmom_case_with_vandermonde_solver("dense");
    const qmomkit::CqmomState state = case_.initial_state();
    assert(state.vandermonde_solver() == qmomkit::VandermondeSolver::dense);

    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(case_);
    assert(initial.ok());

    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomStateTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    const qmomkit::CqmomStateTransportResult transported = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.state.vandermonde_solver() == qmomkit::VandermondeSolver::dense);
    assert(transported.state.vandermonde_solver() == state.vandermonde_solver());
    assert_cqmom_quadrature_near(*transported.quadrature, initial.quadrature, 5.0e-10);
}

void test_state_cqmom_source_evolution_matches_lean()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomLeanTransportState lean_state = lean_state_from_dense(fixed_cqmom_moments(), hierarchy);
    const qmomkit::CqmomState state = fixed_cqmom_case().initial_state();
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto source = [](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            rate[i] = 0.1 * moments.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };

    qmomkit::CqmomStateTransportOptions state_options;
    state_options.scheme = qmomkit::TransportScheme::dc;
    const qmomkit::CqmomStateTransportResult state_result = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, source, state_options);

    qmomkit::CqmomLeanTransportOptions lean_options;
    lean_options.scheme = state_options.scheme;
    const qmomkit::CqmomLeanTransportResult lean_result = qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.1, zero_velocity, source, lean_options);

    assert(state_result.ok());
    assert(lean_result.ok());
    assert(state_result.quadrature.has_value());
    assert(lean_result.quadrature.has_value());
    assert_cqmom_quadrature_near(*state_result.quadrature, *lean_result.quadrature, 5.0e-10);
    assert(state_result.diagnostics.source_status == lean_result.diagnostics.source_status);
}

void test_state_cqmom_failure_propagation_from_source_operator()
{
    const qmomkit::CqmomState state = fixed_cqmom_case().initial_state();
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto bad_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; };

    qmomkit::CqmomStateTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    const qmomkit::CqmomStateTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, bad_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.message.find("Lean CQMOM") == std::string::npos);
    assert(result.diagnostics.message.find("Lean CQMOM") == std::string::npos);
}

void test_state_cqmom_transport_api_does_not_expose_moments()
{
    static_assert(!HasPublicMomentsMember<qmomkit::CqmomStateTransportOptions>::value);
    static_assert(!HasPublicMomentsMember<qmomkit::CqmomStateTransportResult>::value);
    static_assert(std::is_same_v<decltype(std::declval<qmomkit::CqmomStateTransportResult&>().state), qmomkit::CqmomState>);
}

void test_first_order_schemes_differ_for_noncommuting_qmom()
{
    const std::vector<double> moments = raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4);
    const auto velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.35}; };
    const auto source = [](double, const std::vector<double>& state)
    {
        std::vector<double> rhs(state.size(), 0.0);
        for (std::size_t order = 1; order < state.size(); ++order)
        {
            rhs[order] = -0.2 * static_cast<double>(order) * state[order];
        }
        return rhs;
    };

    qmomkit::QmomTransportOptions cd_options;
    cd_options.node_options = tight_node_options();
    cd_options.scheme = qmomkit::TransportScheme::cd;
    const qmomkit::QmomTransportResult cd = qmomkit::advance_transport_qmom(moments, 2, 0.0, 0.2, velocity, source, cd_options);

    qmomkit::QmomTransportOptions dc_options = cd_options;
    dc_options.scheme = qmomkit::TransportScheme::dc;
    const qmomkit::QmomTransportResult dc = qmomkit::advance_transport_qmom(moments, 2, 0.0, 0.2, velocity, source, dc_options);

    assert(cd.ok());
    assert(dc.ok());
    assert(cd.diagnostics.scheme == qmomkit::TransportScheme::cd);
    assert(dc.diagnostics.scheme == qmomkit::TransportScheme::dc);
    assert(std::abs(cd.moments[1] - dc.moments[1]) > 1.0e-4);
}

void test_diagnostics_are_aggregated()
{
    const std::vector<double> moments = raw_moments_1d({1.0}, {2.0}, 2);
    const auto velocity = [](double, const std::vector<double>& xi) { return std::vector<double>{xi[0]}; };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };

    qmomkit::QmomTransportOptions options;
    options.node_options = tight_node_options();
    options.source_options.ssp_rk3_options.estimate_error = true;
    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(moments, 1, 0.0, 0.5, velocity, zero_source, options);

    assert(result.ok());
    assert(result.diagnostics.initial_qmom_status == qmomkit::QmomStatus::success);
    assert(result.diagnostics.post_source_qmom_status == qmomkit::QmomStatus::success);
    assert(result.diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::success);
    assert(result.diagnostics.second_continuous_status == qmomkit::NodeAdvectionStatus::success);
    assert(result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(result.diagnostics.source_diagnostics.error_estimate.has_value());
    assert(result.diagnostics.total_node_advection_steps ==
           result.diagnostics.first_continuous_diagnostics.step_count + result.diagnostics.second_continuous_diagnostics.step_count);
    assert(result.diagnostics.total_node_advection_rhs_evaluations > 0);
}

void test_qmom_transport_node_advection_scheme_selection()
{
    const std::vector<double> moments = raw_moments_1d({-1.0, 2.0}, {0.25, 0.75}, 4);
    const auto velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.5}; };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };

    qmomkit::QmomTransportOptions ssp_options;
    ssp_options.scheme = qmomkit::TransportScheme::cd;
    ssp_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::QmomTransportResult ssp = qmomkit::advance_transport_qmom(moments, 2, 0.0, 0.2, velocity, zero_source, ssp_options);
    assert(ssp.ok());
    assert(ssp.diagnostics.first_continuous_diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    assert(ssp.diagnostics.first_continuous_diagnostics.step_count == 1);
    assert(ssp.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(ssp.diagnostics.total_node_advection_rejected_steps == 0);

    qmomkit::QmomTransportOptions rk_options = ssp_options;
    rk_options.node_options.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    rk_options.node_options.rk23_options.relative_tolerance = 1.0e-10;
    rk_options.node_options.rk23_options.absolute_tolerance = 1.0e-12;
    const qmomkit::QmomTransportResult rk = qmomkit::advance_transport_qmom(moments, 2, 0.0, 0.2, velocity, zero_source, rk_options);
    assert(rk.ok());
    assert(rk.diagnostics.first_continuous_diagnostics.scheme == qmomkit::NodeAdvectionScheme::adaptive_rk23);
    assert(rk.diagnostics.first_continuous_diagnostics.adaptive_rk23_status == qmomkit::AdaptiveRk23Status::success);
    assert(rk.diagnostics.total_node_advection_steps > 0);

    assert_vector_near(ssp.moments, rk.moments, 1.0e-10);
}

void test_cqmom_state_transport_node_advection_scheme_selection()
{
    const qmomkit::CqmomCase cqmom_case = fixed_cqmom_case();
    qmomkit::CqmomState state = cqmom_case.initial_state();
    const auto velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.25); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    qmomkit::CqmomStateTransportOptions ssp_options;
    ssp_options.scheme = qmomkit::TransportScheme::cd;
    ssp_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomStateTransportResult ssp = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, velocity, zero_source, ssp_options);
    assert(ssp.ok());
    assert(ssp.quadrature.has_value());
    assert(ssp.diagnostics.first_continuous_diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    assert(ssp.diagnostics.first_continuous_diagnostics.step_count == 1);
    assert(ssp.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::success);

    qmomkit::CqmomStateTransportOptions rk_options = ssp_options;
    rk_options.node_options.scheme = qmomkit::NodeAdvectionScheme::adaptive_rk23;
    rk_options.node_options.rk23_options.relative_tolerance = 1.0e-10;
    rk_options.node_options.rk23_options.absolute_tolerance = 1.0e-12;
    const qmomkit::CqmomStateTransportResult rk = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, velocity, zero_source, rk_options);
    assert(rk.ok());
    assert(rk.quadrature.has_value());
    assert(rk.diagnostics.first_continuous_diagnostics.scheme == qmomkit::NodeAdvectionScheme::adaptive_rk23);
    assert(rk.diagnostics.first_continuous_diagnostics.adaptive_rk23_status == qmomkit::AdaptiveRk23Status::success);
}

void test_failure_propagation_from_continuous_operator()
{
    const std::vector<double> moments{1.0, 2.0};
    const auto bad_velocity = [](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0}; };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };

    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(moments, 1, 0.0, 0.1, bad_velocity, zero_source);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::continuous_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::first_continuous_half_step);
    assert(result.diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::velocity_failed);
}

void test_failure_propagation_from_source_operator()
{
    const std::vector<double> moments{1.0, 2.0};
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto bad_source = [](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0, 3.0}; };

    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(moments, 1, 0.0, 0.1, zero_velocity, bad_source);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
}

void test_failure_propagation_from_boundary_contribution()
{
    const std::vector<double> moments{1.0, 2.0};
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };
    const qmomkit::MomentBoundaryContribution bad_boundary = [](double, const std::vector<double>& state)
    { return qmomkit::MomentBoundaryRate{std::vector<double>(state.size() + 1, 0.0)}; };

    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(moments, 1, 0.0, 0.1, zero_velocity, zero_source, bad_boundary, {});

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.message.find("boundary contribution") != std::string::npos);
}

void test_boundary_contribution_post_source_inversion_failure()
{
    const std::vector<double> moments{1.0, 0.0};
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };
    const qmomkit::MomentBoundaryContribution removal = [](double, const std::vector<double>& state)
    {
        std::vector<double> rate(state.size(), 0.0);
        rate.front() = -20.0;
        return qmomkit::MomentBoundaryRate{rate};
    };

    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(moments, 1, 0.0, 0.1, zero_velocity, zero_source, removal, {});

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::post_source_inversion);
    assert(result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(result.diagnostics.post_source_qmom_status != qmomkit::QmomStatus::success);
}

void test_failure_propagation_from_inversion()
{
    const std::vector<double> moments{1.0, 0.0, -1.0, 0.0};
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };

    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(moments, 2, 0.0, 0.1, zero_velocity, zero_source);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::initial_inversion);
    assert(result.diagnostics.initial_qmom_status != qmomkit::QmomStatus::success);
}

void test_lean_cqmom_failure_propagation_from_inversion()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomLeanTransportState state{qmomkit::CqmomMomentData{layout, std::vector<double>(layout.required_indices.size(), 0.0)}};
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    const qmomkit::CqmomLeanTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::initial_inversion);
    assert(result.diagnostics.initial_cqmom_status != qmomkit::CqmomStatus::success);
}

void test_lean_cqmom_transport_preserves_invalid_state_messages()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    qmomkit::CqmomLeanTransportState missing_values = lean_state_from_dense(moments, hierarchy);
    missing_values.moments.values.pop_back();
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    const qmomkit::CqmomLeanTransportResult missing = qmomkit::advance_transport_cqmom(missing_values, 0.0, 0.1, zero_velocity);
    assert(!missing.ok());
    assert(missing.status == qmomkit::TransportStatus::invalid_input);
    assert(missing.message == "Lean CQMOM transport moments must match the fixed layout size.");
    assert(missing.diagnostics.message == missing.message);

    qmomkit::CqmomLeanTransportState wrong_layout = lean_state_from_dense(moments, hierarchy);
    wrong_layout.moments.layout.required_indices.pop_back();
    const qmomkit::CqmomLeanTransportResult layout = qmomkit::advance_transport_cqmom(wrong_layout, 0.0, 0.1, zero_velocity);
    assert(!layout.ok());
    assert(layout.status == qmomkit::TransportStatus::invalid_input);
    assert(layout.message == "Lean CQMOM transport requires the fixed CQMOM moment layout.");
    assert(layout.diagnostics.message == layout.message);
}

void test_lean_cqmom_failure_propagation_from_reconstruction()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::CqmomLeanTransportState state = lean_state_from_dense(moments, hierarchy);
    const auto huge_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 1.0e200); };

    qmomkit::CqmomLeanTransportOptions options;
    options.node_options.rk23_options.max_step = 1.0;
    options.scheme = qmomkit::TransportScheme::cd;
    const qmomkit::CqmomLeanTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, 1.0, huge_velocity, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::moment_reconstruction_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::first_moment_reconstruction);
    assert(!result.diagnostics.first_reconstruction_diagnostics.success);
}

void test_lean_cqmom_failure_propagation_from_source_operator()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::CqmomLeanTransportState state = lean_state_from_dense(moments, hierarchy);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto bad_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; };

    qmomkit::CqmomLeanTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    const qmomkit::CqmomLeanTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, bad_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
}

void test_lean_cqmom_failure_propagation_from_post_source_inversion()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::CqmomLeanTransportState state = lean_state_from_dense(moments, hierarchy);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto bad_realizability_source = [](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        rate.front() = -20.0;
        return qmomkit::CqmomMomentRate{rate};
    };

    qmomkit::CqmomLeanTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    const qmomkit::CqmomLeanTransportResult result =
        qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, bad_realizability_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::post_source_inversion);
    assert(result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(result.diagnostics.post_source_cqmom_status != qmomkit::CqmomStatus::success);
}

} // namespace

int main()
{
    test_continuous_only_matches_direct_operator_qmom();
    test_source_only_matches_direct_operator_qmom();
    test_source_plus_boundary_matches_direct_operator_qmom();
    test_adaptive_branch_specific_cqmom_case();
    test_adaptive_state_initialization_preserves_policy();
    test_adaptive_state_continuous_transport_matches_dense();
    test_adaptive_state_source_transport_changes_active_counts();
    test_adaptive_state_preserves_branch_counts_across_steps();
    test_adaptive_state_preserved_source_can_reduce_active_counts();
    test_adaptive_state_source_policy_allows_growth_to_configured_maximum();
    test_adaptive_state_preserved_source_invalidates_descendants_after_root_collapse();
    test_adaptive_state_preserved_source_invalidates_grandchildren_after_second_level_collapse();
    test_adaptive_state_initial_counts_preserve_nonidentity_conditioning_order();
    test_fixed_state_transport_does_not_store_adaptive_branch_counts();
    test_adaptive_state_rejects_duplicate_active_count_metadata();
    test_adaptive_state_rejects_missing_active_count_metadata();
    test_adaptive_state_failure_preserves_active_counts_for_retry();
    test_adaptive_state_post_source_failure_preserves_active_counts_for_retry();
    test_adaptive_state_strang_transport_matches_dense();
    test_adaptive_state_transport_reports_initial_failure();
    test_lean_cqmom_strang_matches_dense_layout_closed_source();
    test_lean_cqmom_cd_matches_dense_layout_closed_source();
    test_lean_cqmom_dc_matches_dense_layout_closed_source();
    test_state_cqmom_transport_matches_lean_fixed_layout();
    test_state_cqmom_transport_preserves_vandermonde_solver_policy();
    test_state_cqmom_source_evolution_matches_lean();
    test_state_cqmom_failure_propagation_from_source_operator();
    test_state_cqmom_transport_api_does_not_expose_moments();
    test_first_order_schemes_differ_for_noncommuting_qmom();
    test_diagnostics_are_aggregated();
    test_qmom_transport_node_advection_scheme_selection();
    test_cqmom_state_transport_node_advection_scheme_selection();
    test_failure_propagation_from_continuous_operator();
    test_failure_propagation_from_source_operator();
    test_failure_propagation_from_boundary_contribution();
    test_boundary_contribution_post_source_inversion_failure();
    test_failure_propagation_from_inversion();
    test_lean_cqmom_failure_propagation_from_inversion();
    test_lean_cqmom_transport_preserves_invalid_state_messages();
    test_lean_cqmom_failure_propagation_from_reconstruction();
    test_lean_cqmom_failure_propagation_from_source_operator();
    test_lean_cqmom_failure_propagation_from_post_source_inversion();
    return 0;
}
