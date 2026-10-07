#include "internal/cqmom_active_counts.hpp"
#include "internal/cqmom_case.hpp"
#include "internal/performance.hpp"
#include "internal/recursive_quadrature.hpp"
#include "qmomkit/qmomkit.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
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
adaptive_branch_counts = false

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

template <typename T, typename = void> struct HasPublicQuadratureMember : std::false_type
{
};

template <typename T> struct HasPublicQuadratureMember<T, std::void_t<decltype(std::declval<T&>().quadrature)>> : std::true_type
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

void assert_qmom_quadrature_near(const qmomkit::QmomResult& actual, const qmomkit::QmomResult& expected, const double tolerance)
{
    assert(actual.weights.size() == expected.weights.size());
    assert(actual.abscissas.size() == expected.abscissas.size());
    assert_vector_near(actual.weights, expected.weights, tolerance);
    assert_vector_near(actual.abscissas, expected.abscissas, tolerance);
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

qmomkit::QmomState qmom_state_from_moments(const std::vector<double>& moments, const std::size_t node_count, const qmomkit::QmomOptions& options = {},
                                           std::optional<qmomkit::AbscissaBounds> abscissa_bounds = {})
{
    return qmomkit::QmomCase{qmomkit::MomentSet(moments), node_count, options, std::move(abscissa_bounds)}.initial_state();
}

qmomkit::QmomMomentSource qmom_source_from_to(const std::vector<double>& initial, const std::vector<double>& target)
{
    return [initial, target](double, const qmomkit::QmomSourceContext& state)
    {
        assert(state.size() == initial.size());
        std::vector<double> rate(state.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            rate[i] = target[i] - initial[i];
        }
        return qmomkit::QmomMomentRate{std::move(rate)};
    };
}

void assert_qmom_state_canonical(const qmomkit::QmomState& state, const double tolerance)
{
    const qmomkit::Quadrature1D quadrature{state.quadrature().abscissas, state.quadrature().weights, qmomkit::QuadratureStatus::success, ""};
    const std::vector<double> reconstructed = qmomkit::reconstruct_raw_moments(quadrature, state.moment_values().size());
    assert_vector_near(state.moment_values(), reconstructed, tolerance);
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

qmomkit::RawMomentTensor near_collapsed_adaptive_cqmom_moments()
{
    const std::vector<std::size_t> extents{4, 4};
    constexpr double secondary_weight = 1.0e-12;
    const std::vector<double> weights{1.0 - secondary_weight, secondary_weight};
    const std::vector<std::vector<double>> abscissas{{0.0, 1.0}, {0.0, 0.0}};
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

qmomkit::RawMomentTensor one_dimensional_cqmom_moments(const std::size_t configured_node_count, const std::vector<double>& abscissas,
                                                       const std::vector<double>& weights)
{
    const std::vector<std::size_t> extents{2 * configured_node_count};
    return {extents, raw_moments_1d(abscissas, weights, extents.front())};
}

qmomkit::CqmomState one_dimensional_adaptive_cqmom_state(const qmomkit::CqmomMomentData& moments)
{
    qmomkit::CqmomCase cqmom_case{moments, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, true, qmomkit::CqmomDiagnosticLevel::basic};
    return qmomkit::internal::CqmomConstructionAccess::create_state(cqmom_case);
}

qmomkit::NodeAdvectionOptions tight_node_options()
{
    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    options.adaptive_options.relative_tolerance = 1.0e-9;
    options.adaptive_options.absolute_tolerance = 1.0e-12;
    options.adaptive_options.initial_step = 0.01;
    return options;
}

} // namespace

namespace qmomkit
{

struct TestCqmomPhaseSpaceEvolutionState
{
    CqmomMomentData moments;
};

} // namespace qmomkit

qmomkit::TestCqmomPhaseSpaceEvolutionState cqmom_state_from_tensor(const qmomkit::RawMomentTensor& moments,
                                                                   const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    return qmomkit::TestCqmomPhaseSpaceEvolutionState{qmomkit::extract_cqmom_moments(moments, layout)};
}

namespace qmomkit
{

struct TestCqmomPhaseSpaceEvolutionOptions : CqmomPhaseSpaceEvolutionOptions
{
    ConditionalMomentOptions conditional_options;
    QmomOptions qmom_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::none;
    VandermondeSolver vandermonde_solver = VandermondeSolver::automatic;
};

struct TestCqmomPhaseSpaceEvolutionResult
{
    TestCqmomPhaseSpaceEvolutionState state;
    std::optional<CqmomQuadrature> quadrature;
    PhaseSpaceEvolutionStatus status = PhaseSpaceEvolutionStatus::success;
    PhaseSpaceEvolutionDiagnostics diagnostics;
    double final_time = 0.0;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == PhaseSpaceEvolutionStatus::success;
    }
};

TestCqmomPhaseSpaceEvolutionResult advance_cqmom_phase_space(const TestCqmomPhaseSpaceEvolutionState& test_state, const double t,
                                                             const double time_step, const PhaseSpaceVelocity& velocity,
                                                             const CqmomMomentSource& source, const TestCqmomPhaseSpaceEvolutionOptions& options = {})
{
    if (test_state.moments.values.size() != test_state.moments.layout.required_indices.size())
    {
        const std::string message = "CQMOM phase-space evolution moments must match the fixed layout size.";
        PhaseSpaceEvolutionDiagnostics diagnostics;
        diagnostics.failed_stage = PhaseSpaceEvolutionStage::none;
        diagnostics.message = message;
        return TestCqmomPhaseSpaceEvolutionResult{test_state, std::nullopt, PhaseSpaceEvolutionStatus::invalid_input, std::move(diagnostics),
                                                  t,          message};
    }
    CqmomCase cqmom_case{test_state.moments,       options.conditional_options, options.qmom_options, options.adaptive_branch_counts,
                         options.diagnostic_level, options.vandermonde_solver};
    CqmomState state = internal::CqmomConstructionAccess::create_state(cqmom_case);
    const CqmomPhaseSpaceEvolutionResult result =
        advance_cqmom_phase_space(state, t, time_step, velocity, source, static_cast<const CqmomPhaseSpaceEvolutionOptions&>(options));
    const CqmomState& result_state = result.ok() ? result.state : state;
    return TestCqmomPhaseSpaceEvolutionResult{
        TestCqmomPhaseSpaceEvolutionState{CqmomMomentData{test_state.moments.layout, result_state.moment_values()}},
        result.ok() ? std::optional<CqmomQuadrature>{result.state.quadrature()} : result.diagnostic_quadrature,
        result.status,
        result.diagnostics,
        result.final_time,
        result.message};
}

TestCqmomPhaseSpaceEvolutionResult advance_cqmom_phase_space(const TestCqmomPhaseSpaceEvolutionState& test_state, const double t,
                                                             const double time_step, const PhaseSpaceVelocity& velocity,
                                                             const CqmomMomentSource& source, const CqmomMomentBoundaryContribution& boundary,
                                                             const TestCqmomPhaseSpaceEvolutionOptions& options)
{
    if (test_state.moments.values.size() != test_state.moments.layout.required_indices.size())
    {
        const std::string message = "CQMOM phase-space evolution moments must match the fixed layout size.";
        PhaseSpaceEvolutionDiagnostics diagnostics;
        diagnostics.failed_stage = PhaseSpaceEvolutionStage::none;
        diagnostics.message = message;
        return TestCqmomPhaseSpaceEvolutionResult{test_state, std::nullopt, PhaseSpaceEvolutionStatus::invalid_input, std::move(diagnostics),
                                                  t,          message};
    }
    CqmomCase cqmom_case{test_state.moments,       options.conditional_options, options.qmom_options, options.adaptive_branch_counts,
                         options.diagnostic_level, options.vandermonde_solver};
    CqmomState state = internal::CqmomConstructionAccess::create_state(cqmom_case);
    const CqmomPhaseSpaceEvolutionResult result =
        advance_cqmom_phase_space(state, t, time_step, velocity, source, boundary, static_cast<const CqmomPhaseSpaceEvolutionOptions&>(options));
    const CqmomState& result_state = result.ok() ? result.state : state;
    return TestCqmomPhaseSpaceEvolutionResult{
        TestCqmomPhaseSpaceEvolutionState{CqmomMomentData{test_state.moments.layout, result_state.moment_values()}},
        result.ok() ? std::optional<CqmomQuadrature>{result.state.quadrature()} : result.diagnostic_quadrature,
        result.status,
        result.diagnostics,
        result.final_time,
        result.message};
}

TestCqmomPhaseSpaceEvolutionResult advance_cqmom_phase_space(const TestCqmomPhaseSpaceEvolutionState& test_state, const double t,
                                                             const double time_step, const PhaseSpaceVelocity& velocity,
                                                             const TestCqmomPhaseSpaceEvolutionOptions& options = {})
{
    const CqmomMomentSource zero_source = [](double, const CqmomSourceContext& context)
    { return CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };
    return advance_cqmom_phase_space(test_state, t, time_step, velocity, zero_source, options);
}

} // namespace qmomkit

namespace
{

qmomkit::CqmomPreparedCase fixed_cqmom_prepared_case()
{
    const auto parsed = qmomkit::config::parse_cqmom_prepared(fixed_cqmom_case_text, "fixed-cqmom-case");
    assert(parsed.ok());
    return parsed.value();
}

qmomkit::CqmomPreparedCase fixed_cqmom_prepared_case_with_vandermonde_solver(const char* const solver)
{
    std::string text = fixed_cqmom_case_text;
    text += "\n[vandermonde]\nsolver = ";
    text += solver;
    text += "\n";
    const auto parsed = qmomkit::config::parse_cqmom_prepared(text, "fixed-cqmom-case-with-vandermonde-solver");
    assert(parsed.ok());
    return parsed.value();
}

void test_fixed_cqmom_phase_space_evolution_successors_share_prepared_metadata()
{
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = fixed_cqmom_prepared_case();
    const auto& prepared_metadata = qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(cqmom_prepared_case);
    assert(prepared_metadata != nullptr);
    const qmomkit::CqmomState initial_state = cqmom_prepared_case.initial_state();
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(initial_state).get() == prepared_metadata.get());

    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& coordinate) { return std::vector<double>(coordinate.size(), 0.0); });
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    const qmomkit::CqmomPhaseSpaceEvolutionResult advection = qmomkit::advance_cqmom_advection(initial_state, 0.0, 0.1, zero_velocity);
    assert(advection.ok());
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(advection.state).get() == prepared_metadata.get());

    const qmomkit::CqmomPhaseSpaceEvolutionResult moment_source = qmomkit::advance_cqmom_moment_source(initial_state, 0.0, 0.1, zero_source);
    assert(moment_source.ok());
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(moment_source.state).get() == prepared_metadata.get());

    for (const qmomkit::PhaseSpaceEvolutionScheme scheme :
         {qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv, qmomkit::PhaseSpaceEvolutionScheme::strang})
    {
        qmomkit::CqmomPhaseSpaceEvolutionOptions options;
        options.scheme = scheme;
        const qmomkit::CqmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_cqmom_phase_space(initial_state, 0.0, 0.1, zero_velocity, zero_source, options);
        assert(result.ok());
        assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(result.state).get() == prepared_metadata.get());
    }

    const qmomkit::CqmomMomentSource failing_source = [](double, const qmomkit::CqmomSourceContext&) -> qmomkit::CqmomMomentRate
    { throw std::runtime_error("prepared metadata rollback failure"); };
    const qmomkit::CqmomPhaseSpaceEvolutionResult failed = qmomkit::advance_cqmom_moment_source(initial_state, 0.0, 0.1, failing_source);
    assert(!failed.ok());
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(failed.state).get() == prepared_metadata.get());
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

qmomkit::CqmomPreparedCase adaptive_cqmom_prepared_case()
{
    const std::string text = cqmom_case_text_from_moments(adaptive_cqmom_moments(), true, true) + "\n[cqmom]\ndiagnostic_level = basic\n";
    const auto parsed = qmomkit::config::parse_cqmom_prepared(text, "adaptive-cqmom-case");
    assert(parsed.ok());
    return parsed.value();
}

qmomkit::CqmomPreparedCase adaptive_cqmom_prepared_case_from_moments(const qmomkit::RawMomentTensor& moments)
{
    const std::string text = cqmom_case_text_from_moments(moments, true, true);
    const auto parsed = qmomkit::config::parse_cqmom_prepared(text, "adaptive-cqmom-case-from-moments");
    assert(parsed.ok());
    return parsed.value();
}

qmomkit::CqmomPreparedCase adaptive_cqmom_prepared_case_from_moments_with_fixed_qmom(const qmomkit::RawMomentTensor& moments)
{
    const std::string text = cqmom_case_text_from_moments(moments, true, false) + "\n[cqmom]\ndiagnostic_level = basic\n";
    const auto parsed = qmomkit::config::parse_cqmom_prepared(text, "adaptive-cqmom-case-from-moments-with-fixed-qmom");
    assert(parsed.ok());
    return parsed.value();
}

qmomkit::CqmomPreparedCase adaptive_cqmom_prepared_case_from_moments_with_hierarchy(const qmomkit::RawMomentTensor& moments,
                                                                                    const std::vector<std::size_t>& node_counts,
                                                                                    const std::vector<std::size_t>& conditioning_order)
{
    const std::string text = cqmom_case_text_from_moments_with_hierarchy(moments, node_counts, conditioning_order, true, true);
    const auto parsed = qmomkit::config::parse_cqmom_prepared(text, "adaptive-cqmom-case-from-moments-with-hierarchy");
    assert(parsed.ok());
    return parsed.value();
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

bool diagnostics_contain_branch_count(const std::vector<qmomkit::CqmomDiagnostics>& diagnostics, const std::size_t level,
                                      const qmomkit::ConditioningPath& path, const std::size_t requested_count, const std::size_t actual_count)
{
    for (const qmomkit::CqmomDiagnostics& item : diagnostics)
    {
        if (item.conditioning_level == level && item.parent_path == path && item.requested_node_count == requested_count &&
            item.actual_node_count == actual_count)
        {
            return true;
        }
    }
    return false;
}

bool transitions_contain(const std::vector<qmomkit::CqmomNodeCountTransition>& transitions, const qmomkit::PhaseSpaceEvolutionStage stage,
                         const std::size_t level, const qmomkit::ConditioningPath& path, const std::optional<std::size_t> before,
                         const std::optional<std::size_t> after, const qmomkit::CqmomNodeCountTransitionKind kind)
{
    for (const qmomkit::CqmomNodeCountTransition& transition : transitions)
    {
        if (transition.stage == stage && transition.conditioning_level == level && transition.parent_path == path &&
            transition.active_node_count_before == before && transition.active_node_count_after == after && transition.kind == kind)
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

void assert_cqmom_moments_near(const qmomkit::CqmomState& lhs, const qmomkit::CqmomState& rhs, const double tolerance)
{
    assert(lhs.hierarchy().node_counts() == rhs.hierarchy().node_counts());
    assert(lhs.hierarchy().conditioning_order() == rhs.hierarchy().conditioning_order());
    assert_vector_near(lhs.moment_values(), rhs.moment_values(), tolerance);
}

void assert_cqmom_state_moments_match_quadrature(const qmomkit::CqmomState& state, const qmomkit::CqmomQuadrature& quadrature, const double tolerance)
{
    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(quadrature, state.hierarchy().moment_extents());
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(state.hierarchy());
    const qmomkit::CqmomMomentData reconstructed_required = qmomkit::extract_cqmom_moments(reconstructed, layout);
    assert_vector_near(state.moment_values(), reconstructed_required.values, tolerance);
}

void assert_cqmom_state_active_counts_match_accepted_quadrature(const qmomkit::CqmomState& state)
{
    const std::optional<qmomkit::CqmomQuadrature>& accepted = qmomkit::internal::CqmomConstructionAccess::state_accepted_quadrature(state);
    assert(accepted.has_value());
    assert(active_counts_equal(state.active_branch_counts(), qmomkit::internal::active_branch_counts_from_quadrature(*accepted)));
}

// A case promoted to state must carry exactly the quadrature that reconstructing
// that case accepts, together with the active branch counts and moments implied
// by it. This is the invariant every evolution entry point relies on.
void assert_cqmom_initial_state_stores_accepted_quadrature(const qmomkit::CqmomPreparedCase& cqmom_prepared_case, const double tolerance)
{
    const qmomkit::CqmomResult expected = qmomkit::compute_cqmom(cqmom_prepared_case);
    assert(expected.ok());

    const qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
    const std::optional<qmomkit::CqmomQuadrature>& accepted = qmomkit::internal::CqmomConstructionAccess::state_accepted_quadrature(state);

    assert(accepted.has_value());
    assert_cqmom_quadrature_near(*accepted, expected.quadrature, tolerance);
    assert_cqmom_state_active_counts_match_accepted_quadrature(state);
    assert_cqmom_state_moments_match_quadrature(state, *accepted, tolerance);
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

// An adaptive case prepares metadata as a fixed one does, so its successors
// must share the same immutable tree rather than rebuild one per step. Without
// this, an adaptive lineage could silently fall back to rebuilding and only a
// benchmark would notice.
void test_adaptive_cqmom_phase_space_evolution_successors_share_prepared_metadata()
{
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = adaptive_cqmom_prepared_case();
    const auto& prepared_metadata = qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(cqmom_prepared_case);
    assert(prepared_metadata != nullptr);
    assert(prepared_metadata->adaptive_branch_counts);
    const qmomkit::CqmomState initial_state = cqmom_prepared_case.initial_state();
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(initial_state).get() == prepared_metadata.get());

    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& coordinate) { return std::vector<double>(coordinate.size(), 0.0); });
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    const qmomkit::CqmomPhaseSpaceEvolutionResult advection = qmomkit::advance_cqmom_advection(initial_state, 0.0, 0.1, zero_velocity);
    assert(advection.ok());
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(advection.state).get() == prepared_metadata.get());

    const qmomkit::CqmomPhaseSpaceEvolutionResult moment_source = qmomkit::advance_cqmom_moment_source(initial_state, 0.0, 0.1, zero_source);
    assert(moment_source.ok());
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(moment_source.state).get() == prepared_metadata.get());

    for (const qmomkit::PhaseSpaceEvolutionScheme scheme :
         {qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv, qmomkit::PhaseSpaceEvolutionScheme::strang})
    {
        qmomkit::CqmomPhaseSpaceEvolutionOptions options;
        options.scheme = scheme;
        const qmomkit::CqmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_cqmom_phase_space(initial_state, 0.0, 0.1, zero_velocity, zero_source, options);
        assert(result.ok());
        assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(result.state).get() == prepared_metadata.get());
    }
}

void test_explicit_cqmom_adaptive_branch_options_are_honored()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = adaptive_cqmom_moments();
    const qmomkit::TestCqmomPhaseSpaceEvolutionState state = cqmom_state_from_tensor(moments, hierarchy);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& values)
    { return qmomkit::CqmomMomentRate{std::vector<double>(values.size(), 0.0)}; };

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options;
    options.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    options.adaptive_branch_counts = true;
    options.diagnostic_level = qmomkit::CqmomDiagnosticLevel::detailed;
    options.vandermonde_solver = qmomkit::VandermondeSolver::dense;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::TestCqmomPhaseSpaceEvolutionResult evolved =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(evolved.ok());
    assert(evolved.quadrature.has_value());
    assert(evolved.quadrature->adaptive_branch_counts_enabled());
    assert(!evolved.quadrature->has_rectangular_storage());
    assert(evolved.diagnostics.cqmom_diagnostic_level == qmomkit::CqmomDiagnosticLevel::detailed);
    assert(evolved.diagnostics.cqmom_vandermonde_solver == qmomkit::VandermondeSolver::dense);
    assert(diagnostics_request_solver(evolved.diagnostics.initial_cqmom_diagnostics, hierarchy.number_of_dimensions(),
                                      qmomkit::VandermondeSolver::dense));
}

void test_explicit_cqmom_diagnostic_level_none_is_propagated()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestCqmomPhaseSpaceEvolutionState state = cqmom_state_from_tensor(moments, hierarchy);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& values)
    { return qmomkit::CqmomMomentRate{std::vector<double>(values.size(), 0.0)}; };

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options;
    options.diagnostic_level = qmomkit::CqmomDiagnosticLevel::none;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::TestCqmomPhaseSpaceEvolutionResult evolved =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(evolved.ok());
    assert(evolved.diagnostics.cqmom_diagnostic_level == qmomkit::CqmomDiagnosticLevel::none);
    assert(evolved.diagnostics.initial_cqmom_diagnostics.empty());
    assert(evolved.diagnostics.post_source_cqmom_diagnostics.empty());
}

void test_explicit_cqmom_adaptive_options_are_honored()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestCqmomPhaseSpaceEvolutionState state = cqmom_state_from_tensor(moments, hierarchy);
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& xi)
        {
            assert(xi.size() == 2);
            return std::vector<double>{0.25, -0.1};
        });
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& values)
    { return qmomkit::CqmomMomentRate{std::vector<double>(values.size(), 0.0)}; };

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    options.node_options.adaptive_options.initial_step = 0.01;
    options.node_options.adaptive_options.max_step = 0.01;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;

    const qmomkit::TestCqmomPhaseSpaceEvolutionResult evolved = qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, velocity, zero_source, options);

    assert(evolved.ok());
    assert(evolved.quadrature.has_value());
    assert(evolved.diagnostics.first_advection_diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
}

void test_explicit_cqmom_source_options_are_honored()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestCqmomPhaseSpaceEvolutionState state = cqmom_state_from_tensor(moments, hierarchy);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& values)
    { return qmomkit::CqmomMomentRate{std::vector<double>(values.size(), 0.0)}; };

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    options.source_options.ssp_rk3_options.max_retries = 2;
    options.source_options.ssp_rk3_options.accept = [](const std::vector<double>&) { return false; };

    const qmomkit::TestCqmomPhaseSpaceEvolutionResult evolved =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(!evolved.ok());
    assert(evolved.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(evolved.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
    assert(evolved.diagnostics.moment_source_retries == 2);
}

// The assertions that matter here are inside the source callback: they check
// that the context QMOMKit hands a CQMOM-native source is well formed - the
// stage time, the dimension and moment counts, agreement between the index and
// value tables, indices inside the tensor extents, and the indexed lookup
// agreeing with the positional one. The count at the end proves the callback
// was invoked at all, without which every check above would pass vacuously.
void test_explicit_cqmom_native_source_context_contract()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestCqmomPhaseSpaceEvolutionState state = cqmom_state_from_tensor(moments, hierarchy);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    std::size_t inspected_context_count = 0;

    const qmomkit::CqmomMomentSource source = [&](const double time, const qmomkit::CqmomSourceContext& context)
    {
        assert_near(time, context.time(), 1.0e-15);
        assert(context.dimension_count() == hierarchy.number_of_dimensions());
        if (context.indices() == state.moments.layout.required_indices)
        {
            ++inspected_context_count;
            assert(context.moment_count() == state.moments.values.size());
            assert(context.indices().size() == context.values().size());
            for (std::size_t entry = 0; entry < context.size(); ++entry)
            {
                const qmomkit::MultiIndexView index = context.index(entry);
                assert(index.size() == moments.extents.size());
                for (std::size_t dimension = 0; dimension < moments.extents.size(); ++dimension)
                {
                    assert(index[dimension] < moments.extents[dimension]);
                }
                const std::optional<double> indexed_value = context.value(index);
                assert(indexed_value.has_value());
                assert_near(*indexed_value, context.value(entry), 1.0e-15);
            }
        }

        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            rate[entry] = -0.05 * context.value(entry);
        }
        return qmomkit::CqmomMomentRate{rate};
    };

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::TestCqmomPhaseSpaceEvolutionResult evolved = qmomkit::advance_cqmom_phase_space(state, 0.25, 0.02, zero_velocity, source, options);

    assert(evolved.ok());
    assert(inspected_context_count > 0);
}

void test_explicit_cqmom_native_source_failures_are_reported()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestCqmomPhaseSpaceEvolutionState state = cqmom_state_from_tensor(moments, hierarchy);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const qmomkit::CqmomMomentSource bad_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size() + 1, 0.0)}; };

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::TestCqmomPhaseSpaceEvolutionResult evolved =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, bad_source, options);

    assert(!evolved.ok());
    assert(evolved.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(evolved.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    assert(evolved.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
}

void test_explicit_cqmom_native_source_invalid_values_are_reported()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestCqmomPhaseSpaceEvolutionState state = cqmom_state_from_tensor(moments, hierarchy);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::CqmomMomentSource nonfinite_source = [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        rate.front() = std::numeric_limits<double>::infinity();
        return qmomkit::CqmomMomentRate{rate};
    };
    const qmomkit::TestCqmomPhaseSpaceEvolutionResult nonfinite =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, nonfinite_source, options);
    assert(!nonfinite.ok());
    assert(nonfinite.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(nonfinite.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    assert(nonfinite.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::source_failed);

    const qmomkit::CqmomMomentSource throwing_source = [](double, const qmomkit::CqmomSourceContext&) -> qmomkit::CqmomMomentRate
    { throw std::runtime_error("native CQMOM source failure"); };
    const qmomkit::TestCqmomPhaseSpaceEvolutionResult throwing =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, throwing_source, options);
    assert(!throwing.ok());
    assert(throwing.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(throwing.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    assert(throwing.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(throwing.message.find("native CQMOM source failure") != std::string::npos);
}

void test_explicit_cqmom_source_growth_is_automatic()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = collapsed_branch_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const qmomkit::TestCqmomPhaseSpaceEvolutionState state{initial_required};
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions growth_options;
    growth_options.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    growth_options.adaptive_branch_counts = true;
    growth_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    growth_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::TestCqmomPhaseSpaceEvolutionResult growth =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 1.0, zero_velocity, source_from_to(initial_required, target_required), growth_options);

    // An accepted moment-source stage may grow an active branch toward the
    // configured maximum, which is what raises the terminal node count here.
    assert(growth.ok());
    assert(growth.quadrature.has_value());
    assert(growth.quadrature->total_nodes() == 3);
}

void test_public_qmom_advection_stage_advances_canonical_state()
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = true;
    qmom_options.rmin = 1.0e-10;
    qmom_options.eabs = 1.0e-12;
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.75, 1.25}, {0.4, 0.6}, 6), 3, qmom_options);
    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>{0.2 * xi[0] + 0.1}; });
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options = tight_node_options();

    const qmomkit::QmomPhaseSpaceEvolutionResult stage = qmomkit::advance_qmom_advection(state, 0.0, 0.2, velocity, options);
    const qmomkit::QmomNodeAdvectionResult direct = qmomkit::advance_nodes(state.quadrature(), 0.0, 0.2, velocity, options.node_options);

    assert(stage.ok());
    assert(direct.ok());
    assert(stage.final_time == 0.2);
    assert(stage.state.active_node_count() == state.active_node_count());
    assert(stage.state.active_node_count() == 2);
    assert_qmom_state_canonical(stage.state, 5.0e-10);
    assert_qmom_quadrature_near(stage.state.quadrature(), qmomkit::QmomResult{direct.quadrature.abscissas, direct.quadrature.weights}, 5.0e-10);
    assert_vector_near(stage.state.moment_values(), qmomkit::reconstruct_raw_moments(direct.quadrature, state.moment_values().size()), 5.0e-10);
    assert(stage.diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success);
    assert(stage.diagnostics.first_reconstruction_diagnostics.success);
}

void test_public_qmom_advection_stage_failure_rolls_back()
{
    const qmomkit::QmomState state = qmom_state_from_moments({1.0, 2.0}, 1);
    const qmomkit::PhaseSpaceVelocity bad_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0}; });

    const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_advection(state, 0.0, 0.1, bad_velocity);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::advection_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::first_advection_half_step);
    assert_vector_near(result.state.moment_values(), state.moment_values(), 0.0);
    assert_qmom_quadrature_near(result.state.quadrature(), state.quadrature(), 0.0);
}

void test_qmom_bounded_state_advection_and_source_reconstruction_semantics()
{
    const qmomkit::AbscissaBounds bounds{0.0, 1.0};
    const qmomkit::QmomState initial = qmom_state_from_moments({1.0, 0.5}, 1, {}, bounds);
    const qmomkit::PhaseSpaceVelocity outward_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{2.0}; });

    const qmomkit::QmomPhaseSpaceEvolutionResult advected = qmomkit::advance_qmom_advection(initial, 0.0, 1.0, outward_velocity);
    assert(advected.ok());
    assert(advected.state.abscissa_bounds().has_value());
    assert(advected.state.abscissa_bounds()->lower == bounds.lower);
    assert(advected.state.abscissa_bounds()->upper == bounds.upper);
    assert(advected.state.quadrature().abscissas[0] > bounds.upper);

    const qmomkit::QmomMomentSource untracked_moment_source = [](double, const qmomkit::QmomSourceContext& context)
    {
        static_cast<void>(context.moment(2));
        return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)};
    };
    const qmomkit::QmomPhaseSpaceEvolutionResult lazy_reconstruction =
        qmomkit::advance_qmom_moment_source(advected.state, 1.0, 0.1, untracked_moment_source);
    assert(!lazy_reconstruction.ok());
    assert(lazy_reconstruction.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(lazy_reconstruction.message.find("source-stage quadrature reconstruction failed") != std::string::npos);

    const qmomkit::QmomMomentSource move_outside_bounds = [](double, const qmomkit::QmomSourceContext&)
    { return qmomkit::QmomMomentRate{std::vector<double>{0.0, 1.0}}; };
    const qmomkit::QmomPhaseSpaceEvolutionResult post_source = qmomkit::advance_qmom_moment_source(initial, 0.0, 1.0, move_outside_bounds);
    assert(!post_source.ok());
    assert(post_source.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    assert(post_source.diagnostics.post_source_qmom_status == qmomkit::QmomStatus::not_realizable);
    assert_vector_near(post_source.state.moment_values(), initial.moment_values(), 0.0);
}

void test_qmom_bounded_post_source_automatic_scaling_is_unit_covariant()
{
    // The rule must be sufficiently separated that its double-precision raw
    // moments retain a two-node oracle at the asserted 1e-10 tolerance.  A
    // 1e-4 relative separation makes the moment inversion ill-conditioned
    // before automatic scaling is applied; scaling changes units, not that
    // intrinsic conditioning.
    const std::vector<double> reference_nodes{1.0, 1.05};
    const std::vector<double> weights{0.5, 0.5};
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = true;
    qmom_options.automatic_scaling = true;

    for (const double factor : std::vector<double>{1.0e-9, 1.0e6})
    {
        std::vector<double> physical_nodes = reference_nodes;
        for (double& node : physical_nodes)
        {
            node *= factor;
        }
        const qmomkit::AbscissaBounds bounds{0.5 * factor, 1.5 * factor};
        const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d(physical_nodes, weights, 4), 2, qmom_options, bounds);
        const qmomkit::QmomMomentSource zero_source = [](double, const qmomkit::QmomSourceContext& context)
        { return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

        qmomkit::QmomPhaseSpaceEvolutionOptions scaled_options;
        const qmomkit::QmomPhaseSpaceEvolutionResult scaled = qmomkit::advance_qmom_moment_source(state, 0.0, 1.0, zero_source, scaled_options);
        assert(scaled.ok());
        assert(scaled.state.active_node_count() == 2);
        for (std::size_t node = 0; node < reference_nodes.size(); ++node)
        {
            assert_near(scaled.state.quadrature().weights[node], weights[node], 1.0e-10);
            assert_near(scaled.state.quadrature().abscissas[node] / factor, reference_nodes[node], 1.0e-8);
        }

        qmomkit::QmomOptions unscaled_qmom_options = qmom_options;
        unscaled_qmom_options.automatic_scaling = false;
        const qmomkit::QmomState unscaled_state =
            qmom_state_from_moments(raw_moments_1d(physical_nodes, weights, 4), 2, unscaled_qmom_options, bounds);
        qmomkit::QmomPhaseSpaceEvolutionOptions unscaled_options = scaled_options;
        const qmomkit::QmomPhaseSpaceEvolutionResult unscaled =
            qmomkit::advance_qmom_moment_source(unscaled_state, 0.0, 1.0, zero_source, unscaled_options);
        assert(unscaled.ok());
        // Adaptive reduction does not depend on the unit, so the representation
        // keeps both nodes without scaling as well.
        assert(unscaled.state.active_node_count() == 2);
    }
}

void test_public_qmom_moment_source_stage_fixed_mode_is_strict()
{
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.0}, {0.4, 0.6}, 4), 2, qmomkit::QmomOptions{false});
    const std::vector<double> one_node_moments = raw_moments_1d({0.0}, {1.0}, 4);
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    const qmomkit::QmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_qmom_moment_source(state, 0.0, 1.0, qmom_source_from_to(state.moment_values(), one_node_moments), options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion);
    assert(result.diagnostics.post_source_qmom_status != qmomkit::QmomStatus::success);
    assert_vector_near(result.state.moment_values(), state.moment_values(), 0.0);
    assert_qmom_quadrature_near(result.state.quadrature(), state.quadrature(), 0.0);
}

void test_public_qmom_moment_source_stage_adaptive_policy()
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = true;
    qmom_options.rmin = 1.0e-10;
    qmom_options.eabs = 1.0e-12;
    qmomkit::QmomPhaseSpaceEvolutionOptions options;

    const qmomkit::QmomState two_node_state = qmom_state_from_moments(raw_moments_1d({-0.75, 1.25}, {0.4, 0.6}, 6), 3, qmom_options);
    const qmomkit::QmomPhaseSpaceEvolutionResult retained = qmomkit::advance_qmom_moment_source(
        two_node_state, 0.0, 0.2,
        [](double, const qmomkit::QmomSourceContext& moments) { return qmomkit::QmomMomentRate{std::vector<double>(moments.size(), 0.0)}; }, options);
    assert(retained.ok());
    assert(retained.state.active_node_count() == two_node_state.active_node_count());
    assert_qmom_state_canonical(retained.state, 1.0e-12);

    const std::vector<double> one_node_moments = raw_moments_1d({0.5}, {1.0}, 6);
    const qmomkit::QmomPhaseSpaceEvolutionResult reduction =
        qmomkit::advance_qmom_moment_source(two_node_state, 0.0, 1.0, qmom_source_from_to(two_node_state.moment_values(), one_node_moments), options);
    assert(reduction.ok());
    assert(reduction.state.active_node_count() == 1);
    assert_qmom_state_canonical(reduction.state, 1.0e-12);

    const qmomkit::QmomState one_node_state = qmom_state_from_moments(one_node_moments, 3, qmom_options);
    const std::vector<double> two_node_moments = raw_moments_1d({-0.25, 1.5}, {0.35, 0.65}, 6);
    const qmomkit::QmomPhaseSpaceEvolutionResult growth =
        qmomkit::advance_qmom_moment_source(one_node_state, 0.0, 1.0, qmom_source_from_to(one_node_state.moment_values(), two_node_moments), options);
    assert(growth.ok());
    assert(growth.state.active_node_count() == 2);
    assert_qmom_state_canonical(growth.state, 1.0e-12);

    const std::vector<double> three_node_moments = raw_moments_1d({-1.0, 0.25, 1.5}, {0.2, 0.3, 0.5}, 6);
    const qmomkit::QmomPhaseSpaceEvolutionResult capped_growth = qmomkit::advance_qmom_moment_source(
        one_node_state, 0.0, 1.0, qmom_source_from_to(one_node_state.moment_values(), three_node_moments), options);
    assert(capped_growth.ok());
    assert(capped_growth.state.active_node_count() == 2);
    assert_qmom_state_canonical(capped_growth.state, 1.0e-12);
}

void test_public_qmom_moment_source_stage_failures_roll_back()
{
    const qmomkit::QmomState state = qmom_state_from_moments({1.0, 2.0}, 1);
    const auto bad_source = [](double, const qmomkit::QmomSourceContext&) { return qmomkit::QmomMomentRate{std::vector<double>{1.0, 2.0, 3.0}}; };
    const qmomkit::QmomPhaseSpaceEvolutionResult source_failure = qmomkit::advance_qmom_moment_source(state, 0.0, 0.1, bad_source);
    assert(!source_failure.ok());
    assert(source_failure.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(source_failure.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    assert_vector_near(source_failure.state.moment_values(), state.moment_values(), 0.0);

    const auto nonrealizable_source = [](double, const qmomkit::QmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        rate.front() = -20.0;
        return qmomkit::QmomMomentRate{std::move(rate)};
    };
    const qmomkit::QmomPhaseSpaceEvolutionResult inversion_failure = qmomkit::advance_qmom_moment_source(state, 0.0, 0.1, nonrealizable_source);
    assert(!inversion_failure.ok());
    assert(inversion_failure.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    assert(inversion_failure.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion);
    assert_vector_near(inversion_failure.state.moment_values(), state.moment_values(), 0.0);
    assert_qmom_quadrature_near(inversion_failure.state.quadrature(), state.quadrature(), 0.0);
}

void test_qmom_source_context_uses_current_rk_stage_and_caches_quadrature()
{
    const std::vector<double> initial_moments = raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4);
    qmomkit::QmomOptions qmom_options;
    qmom_options.automatic_scaling = true;
    const qmomkit::QmomState state = qmom_state_from_moments(initial_moments, 2, qmom_options);
    const qmomkit::QmomResult initial_quadrature = qmomkit::compute_qmom(qmomkit::MomentSet(initial_moments), 2, qmom_options);
    assert(initial_quadrature.ok());
    const qmomkit::Quadrature1D quadrature{initial_quadrature.abscissas, initial_quadrature.weights, qmomkit::QuadratureStatus::success, ""};
    const std::vector<double> initial_extended = qmomkit::reconstruct_raw_moments(quadrature, 6);
    constexpr double lambda = 0.1;
    constexpr double time_step = 0.2;
    std::size_t calls = 0;

    const qmomkit::QmomMomentSource source = [&](double, const qmomkit::QmomSourceContext& context)
    {
        assert_near(context.moment(0), context.value(0), 0.0);
        const double density = context.moment(0);
        assert_near(context.moment(4), initial_extended[4] * density / initial_extended[0], 2.0e-12);
        assert_near(context.moment(5), initial_extended[5] * density / initial_extended[0], 2.0e-12);
        ++calls;
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t order = 0; order < rate.size(); ++order)
        {
            rate[order] = lambda * context.value(order);
        }
        return qmomkit::QmomMomentRate{std::move(rate)};
    };

    qmomkit::internal::performance::CounterSession session;
    qmomkit::QmomPhaseSpaceEvolutionOptions evolution_options;
    std::optional<qmomkit::QmomPhaseSpaceEvolutionResult> result;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        result.emplace(qmomkit::advance_qmom_moment_source(state, 0.0, time_step, source, evolution_options));
    }

    assert(result->ok());
    assert(calls == 3);
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(session.snapshot().reconstruction.qmom_reconstructions == 4);
#else
    assert(session.snapshot().reconstruction.qmom_reconstructions == 0);
#endif
}

// The boundary callback receives the same context as the source callback, so it
// can read the current stage and derive moments from it rather than being handed
// only the raw moment vector. This pins that the two callbacks observe the same
// stage within one invocation.
void test_qmom_boundary_context_matches_source_context()
{
    const std::vector<double> abscissas{1.0, 4.0};
    const std::vector<double> weights{0.25, 0.75};
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d(abscissas, weights, 6), 2);

    std::vector<double> source_view;
    std::vector<double> boundary_view;
    double source_generalized = 0.0;
    double boundary_generalized = 0.0;
    bool boundary_called = false;

    const qmomkit::QmomMomentSource source = [&](double, const qmomkit::QmomSourceContext& context)
    {
        source_view.assign(context.values().begin(), context.values().end());
        source_generalized = context.generalized_moment(0.5);
        return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)};
    };
    const qmomkit::QmomMomentBoundaryContribution boundary = [&](double, const qmomkit::QmomSourceContext& context)
    {
        boundary_view.assign(context.values().begin(), context.values().end());
        // Available only because the boundary callback now receives the context.
        boundary_generalized = context.generalized_moment(0.5);
        assert_near(context.generalized_moment(0.0), 1.0, 2.0e-12);
        boundary_called = true;
        return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)};
    };

    const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_moment_source(state, 0.0, 0.1, source, boundary, {});
    assert(result.ok());
    assert(boundary_called);
    assert(boundary_view.size() == source_view.size());
    for (std::size_t i = 0; i < source_view.size(); ++i)
    {
        assert_near(boundary_view[i], source_view[i], 0.0);
    }
    assert_near(boundary_generalized, source_generalized, 0.0);
    assert_near(boundary_generalized, 0.25 * std::sqrt(1.0) + 0.75 * std::sqrt(4.0), 2.0e-12);
}

void test_qmom_source_context_generalized_moment_semantics()
{
    const std::vector<double> positive_abscissas{1.0, 4.0};
    const std::vector<double> positive_weights{0.25, 0.75};
    const qmomkit::QmomState positive_state = qmom_state_from_moments(raw_moments_1d(positive_abscissas, positive_weights, 6), 2);
    bool positive_checked = false;
    const qmomkit::QmomMomentSource positive_source = [&positive_checked](double, const qmomkit::QmomSourceContext& context)
    {
        assert_near(context.generalized_moment(0.0), 1.0, 2.0e-12);
        assert_near(context.generalized_moment(2.0), context.value(2), 2.0e-10);
        assert_near(context.generalized_moment(0.5), 0.25 * std::sqrt(1.0) + 0.75 * std::sqrt(4.0), 2.0e-12);
        assert_near(context.generalized_moment(1.5), 0.25 * std::pow(1.0, 1.5) + 0.75 * std::pow(4.0, 1.5), 2.0e-12);

        const auto expect_invalid_order = [&context](const double order)
        {
            bool rejected = false;
            try
            {
                static_cast<void>(context.generalized_moment(order));
            }
            catch (const std::invalid_argument&)
            {
                rejected = true;
            }
            assert(rejected);
        };
        expect_invalid_order(-0.5);
        expect_invalid_order(std::numeric_limits<double>::quiet_NaN());
        expect_invalid_order(std::numeric_limits<double>::infinity());
        expect_invalid_order(-std::numeric_limits<double>::infinity());
        positive_checked = true;
        return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)};
    };
    const qmomkit::QmomPhaseSpaceEvolutionResult positive_result = qmomkit::advance_qmom_moment_source(positive_state, 0.0, 0.1, positive_source);
    assert(positive_result.ok());
    assert(positive_checked);

    const qmomkit::QmomState negative_state = qmom_state_from_moments(raw_moments_1d({-1.0, 2.0}, {0.4, 0.6}, 6), 2);
    bool negative_checked = false;
    const qmomkit::QmomMomentSource negative_source = [&negative_checked](double, const qmomkit::QmomSourceContext& context)
    {
        bool rejected_fractional_power = false;
        try
        {
            static_cast<void>(context.generalized_moment(0.5));
        }
        catch (const std::domain_error&)
        {
            rejected_fractional_power = true;
        }
        assert(rejected_fractional_power);
        assert_near(context.generalized_moment(3.0), context.value(3), 2.0e-10);
        negative_checked = true;
        return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)};
    };
    const qmomkit::QmomPhaseSpaceEvolutionResult negative_result = qmomkit::advance_qmom_moment_source(negative_state, 0.0, 0.1, negative_source);
    assert(negative_result.ok());
    assert(negative_checked);
}

void test_qmom_source_context_generalized_moment_uses_current_rk_stage()
{
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({1.0, 4.0}, {0.25, 0.75}, 6), 2);
    constexpr double coefficient = 0.5;
    constexpr double time_step = 0.2;
    constexpr double initial_generalized_moment = 1.75;
    std::vector<double> observed_generalized_moments;
    std::vector<double> observed_densities;
    const qmomkit::QmomMomentSource source = [&](double, const qmomkit::QmomSourceContext& context)
    {
        const double generalized_moment = context.generalized_moment(0.5);
        observed_generalized_moments.push_back(generalized_moment);
        observed_densities.push_back(context.value(0));
        // The source scales every raw moment equally, so the stage quadrature
        // keeps its abscissas and scales its weights by the stage density.
        assert_near(generalized_moment, initial_generalized_moment * context.value(0), 2.0e-10);
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t order = 0; order < rate.size(); ++order)
        {
            rate[order] = coefficient * context.value(order);
        }
        return qmomkit::QmomMomentRate{std::move(rate)};
    };

    const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_moment_source(state, 0.0, time_step, source);
    assert(result.ok());
    assert(observed_generalized_moments.size() == 3);
    assert(observed_densities.size() == 3);
    assert_near(observed_generalized_moments[0], initial_generalized_moment, 2.0e-10);
    assert_near(observed_generalized_moments[1], initial_generalized_moment * 1.1, 2.0e-10);
    assert_near(observed_generalized_moments[2], initial_generalized_moment * 1.0525, 2.0e-10);
    assert(observed_generalized_moments[1] != observed_generalized_moments[0]);
    assert(observed_generalized_moments[2] != observed_generalized_moments[0]);
}

// The node visitor sees the quadrature of each SSP-RK3 stage. It shares that
// quadrature with generalized_moment(...), so a stage is reconstructed once
// however many derived requests the source makes.
void test_qmom_source_context_visits_nodes_from_each_rk_stage()
{
    const std::vector<double> abscissas{1.0, 4.0};
    const std::vector<double> weights{0.25, 0.75};
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d(abscissas, weights, 4), 2);
    constexpr double coefficient = 0.5;
    constexpr double time_step = 0.2;
    std::vector<double> observed_densities;
    const qmomkit::QmomMomentSource source = [&](double, const qmomkit::QmomSourceContext& context)
    {
        std::vector<double> node_weights;
        std::vector<double> node_abscissas;
        context.for_each_node(
            [&node_weights, &node_abscissas](const double weight, const double abscissa)
            {
                node_weights.push_back(weight);
                node_abscissas.push_back(abscissa);
            });

        // The source scales every raw moment equally, so each stage quadrature
        // keeps the abscissas and scales the weights by the stage density.
        const double density = context.value(0);
        assert(node_weights.size() == abscissas.size());
        double half_order_moment = 0.0;
        for (std::size_t node = 0; node < node_weights.size(); ++node)
        {
            assert_near(node_abscissas[node], abscissas[node], 2.0e-12);
            assert_near(node_weights[node], weights[node] * density, 2.0e-12);
            half_order_moment += node_weights[node] * std::sqrt(node_abscissas[node]);
        }
        assert_near(context.generalized_moment(0.5), half_order_moment, 1.0e-14);

        std::size_t second_visit_nodes = 0;
        context.for_each_node([&second_visit_nodes](double, double) { ++second_visit_nodes; });
        assert(second_visit_nodes == node_weights.size());

        observed_densities.push_back(density);
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t order = 0; order < rate.size(); ++order)
        {
            rate[order] = coefficient * context.value(order);
        }
        return qmomkit::QmomMomentRate{std::move(rate)};
    };

    qmomkit::internal::performance::CounterSession session;
    std::optional<qmomkit::QmomPhaseSpaceEvolutionResult> result;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        result.emplace(qmomkit::advance_qmom_moment_source(state, 0.0, time_step, source));
    }

    assert(result->ok());
    assert(observed_densities.size() == 3);
    assert_near(observed_densities[0], 1.0, 1.0e-15);
    assert_near(observed_densities[1], 1.1, 1.0e-15);
    assert_near(observed_densities[2], 1.0525, 1.0e-15);
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    // One reconstruction per stage and one for the final state.
    assert(session.snapshot().reconstruction.qmom_reconstructions == 4);
#else
    assert(session.snapshot().reconstruction.qmom_reconstructions == 0);
#endif
}

void test_configured_qmom_source_realizability_policy_matches_programmatic_adapter()
{
    const qmomkit::QmomState state = qmom_state_from_moments({1.0, 0.0, 1.0, 0.0}, 2, qmomkit::QmomOptions{false});
    const auto source_with_negative_second_moment = [](std::size_t& calls)
    {
        return [&calls](double, const qmomkit::QmomSourceContext& moments)
        {
            ++calls;
            std::vector<double> rate(moments.size(), 0.0);
            rate[2] = -10.0;
            return qmomkit::QmomMomentRate{std::move(rate)};
        };
    };
    const auto parse_options = [](const std::string& check_mode)
    {
        const auto parsed = qmomkit::config::parse_qmom_phase_space_evolution_options("[time_integration.realizability]\ncheck_mode = " + check_mode +
                                                                                      "\nsupport = real_line\n");
        assert(parsed.ok());
        return parsed.value();
    };

    std::size_t disabled_calls = 0;
    const qmomkit::QmomPhaseSpaceEvolutionResult disabled =
        qmomkit::advance_qmom_moment_source(state, 0.0, 0.2, source_with_negative_second_moment(disabled_calls), parse_options("none"));
    assert(!disabled.ok());
    assert(disabled.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    assert(disabled.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(disabled_calls == 3);

    std::size_t final_only_calls = 0;
    const qmomkit::QmomPhaseSpaceEvolutionOptions final_only_options = parse_options("final_only");
    const qmomkit::QmomPhaseSpaceEvolutionResult final_only =
        qmomkit::advance_qmom_moment_source(state, 0.0, 0.2, source_with_negative_second_moment(final_only_calls), final_only_options);
    assert(!final_only.ok());
    assert(final_only.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(final_only.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
    assert(final_only_calls == 3);
    assert_vector_near(final_only.state.moment_values(), state.moment_values(), 0.0);

    std::size_t boundary_calls = 0;
    const qmomkit::QmomMomentSource zero_source = [](double, const qmomkit::QmomSourceContext& moments)
    { return qmomkit::QmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    const qmomkit::QmomMomentBoundaryContribution boundary_source = [&boundary_calls](double, const qmomkit::QmomSourceContext& moments)
    {
        ++boundary_calls;
        std::vector<double> rate(moments.size(), 0.0);
        rate[2] = -10.0;
        return qmomkit::QmomMomentRate{std::move(rate)};
    };
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& location) { return std::vector<double>(location.size(), 0.0); });
    const qmomkit::QmomPhaseSpaceEvolutionResult boundary_phase_space =
        qmomkit::advance_qmom_phase_space(state, 0.0, 0.2, zero_velocity, zero_source, boundary_source, final_only_options);
    assert(!boundary_phase_space.ok());
    assert(boundary_phase_space.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(boundary_phase_space.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
    assert(boundary_calls == 3);
    assert_vector_near(boundary_phase_space.state.moment_values(), state.moment_values(), 0.0);

    qmomkit::QmomPhaseSpaceEvolutionOptions programmatic_options;
    programmatic_options.source_options.ssp_rk3_options.accept = qmomkit::make_realizability_acceptance_predicate(
        final_only_options.realizability_policy, [](const std::vector<double>& values) { return qmomkit::MomentSet(values); });
    std::size_t programmatic_calls = 0;
    const qmomkit::QmomPhaseSpaceEvolutionResult programmatic =
        qmomkit::advance_qmom_moment_source(state, 0.0, 0.2, source_with_negative_second_moment(programmatic_calls), programmatic_options);
    assert(!programmatic.ok());
    assert(programmatic.status == final_only.status);
    assert(programmatic.diagnostics.moment_source_status == final_only.diagnostics.moment_source_status);
    assert(programmatic_calls == final_only_calls);
    assert_vector_near(programmatic.state.moment_values(), final_only.state.moment_values(), 0.0);

    std::size_t every_stage_calls = 0;
    const qmomkit::QmomPhaseSpaceEvolutionResult every_stage =
        qmomkit::advance_qmom_moment_source(state, 0.0, 0.2, source_with_negative_second_moment(every_stage_calls), parse_options("every_stage"));
    assert(!every_stage.ok());
    assert(every_stage.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(every_stage.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
    assert(every_stage_calls == 1);
    assert_vector_near(every_stage.state.moment_values(), state.moment_values(), 0.0);
}

// A source that translates the distribution by velocity v,
// dm_k/dt = k v m_{k-1}, keeps every moment sequence realizable on the real
// line, so a configured policy can reject it only through its support. For
// four moments the generator of this linear system is nilpotent, A^4 = 0, so
// its exact propagator is the cubic Taylor polynomial that SSP-RK3 applies,
// and an accepted state holds the exactly translated moments.
void test_configured_qmom_source_realizability_policy_enforces_its_support()
{
    const std::vector<double> points{0.5, 1.5};
    const std::vector<double> weights{0.5, 0.5};
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d(points, weights, 4), 2);
    const auto translation = [](const double velocity)
    {
        return [velocity](double, const qmomkit::QmomSourceContext& moments)
        {
            std::vector<double> rate(moments.size(), 0.0);
            for (std::size_t order = 1; order < moments.size(); ++order)
            {
                rate[order] = static_cast<double>(order) * velocity * moments.moment(order - 1);
            }
            return qmomkit::QmomMomentRate{std::move(rate)};
        };
    };
    const auto parse_options = [](const std::string& support)
    {
        const auto parsed =
            qmomkit::config::parse_qmom_phase_space_evolution_options("[time_integration.realizability]\ncheck_mode = final_only\n" + support);
        assert(parsed.ok());
        return parsed.value();
    };
    const auto translated_moments = [&points, &weights](const double shift)
    { return raw_moments_1d({points[0] + shift, points[1] + shift}, weights, 4); };
    const auto assert_rejected = [&state](const qmomkit::QmomPhaseSpaceEvolutionResult& result)
    {
        assert(!result.ok());
        assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
        assert(result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
        assert_vector_near(result.state.moment_values(), state.moment_values(), 0.0);
    };

    // Moving by -2 places both points below zero.
    const qmomkit::QmomPhaseSpaceEvolutionOptions real_line = parse_options("support = real_line\n");
    const qmomkit::QmomPhaseSpaceEvolutionOptions nonnegative = parse_options("support = nonnegative\n");
    const qmomkit::QmomPhaseSpaceEvolutionResult left_on_real_line =
        qmomkit::advance_qmom_moment_source(state, 0.0, 0.2, translation(-10.0), real_line);
    assert(left_on_real_line.ok());
    assert_vector_near(left_on_real_line.state.moment_values(), translated_moments(-2.0), 1.0e-12);
    assert_rejected(qmomkit::advance_qmom_moment_source(state, 0.0, 0.2, translation(-10.0), nonnegative));

    // On [0, 2], moving by 0.2 stays inside and moving by 2 leaves it.
    const qmomkit::QmomPhaseSpaceEvolutionOptions interval = parse_options("support = bounded_interval\nlower = 0.0\nupper = 2.0\n");
    const qmomkit::QmomPhaseSpaceEvolutionResult inside = qmomkit::advance_qmom_moment_source(state, 0.0, 0.2, translation(1.0), interval);
    assert(inside.ok());
    assert_vector_near(inside.state.moment_values(), translated_moments(0.2), 1.0e-12);
    assert_rejected(qmomkit::advance_qmom_moment_source(state, 0.0, 0.2, translation(10.0), interval));
    assert(qmomkit::advance_qmom_moment_source(state, 0.0, 0.2, translation(10.0), real_line).ok());
}

void test_public_qmom_moment_source_stage_boundary_matches_source_update()
{
    const qmomkit::QmomState state = qmom_state_from_moments({1.0, 2.0}, 1);
    const qmomkit::QmomMomentSource source = [](double, const qmomkit::QmomSourceContext& moments)
    {
        assert(moments.size() == 2);
        return qmomkit::QmomMomentRate{std::vector<double>{0.0, 1.0}};
    };
    const qmomkit::QmomMomentBoundaryContribution boundary = [](double, const qmomkit::QmomSourceContext& moments)
    {
        assert(moments.size() == 2);
        return qmomkit::QmomMomentRate{std::vector<double>{0.25, -0.5}};
    };

    const qmomkit::QmomPhaseSpaceEvolutionResult stage = qmomkit::advance_qmom_moment_source(state, 0.125, 0.25, source, boundary, {});
    const qmomkit::QmomMomentSourceUpdateResult direct =
        qmomkit::advance_qmom_moments_with_source(state.moment_values(), 0.125, 0.25, source, boundary, {});

    assert(stage.ok());
    assert(direct.ok());
    assert_vector_near(stage.state.moment_values(), direct.moments, 1.0e-13);
    assert_qmom_state_canonical(stage.state, 1.0e-13);
}

void assert_qmom_phase_space_evolution_results_match(const qmomkit::QmomPhaseSpaceEvolutionResult& actual,
                                                     const qmomkit::QmomPhaseSpaceEvolutionResult& expected, const double tolerance)
{
    assert(actual.ok());
    assert(expected.ok());
    assert_near(actual.final_time, expected.final_time, tolerance);
    assert(actual.state.active_node_count() == expected.state.active_node_count());
    assert_vector_near(actual.state.moment_values(), expected.state.moment_values(), tolerance);
    assert_qmom_quadrature_near(actual.state.quadrature(), expected.state.quadrature(), tolerance);
}

void test_public_qmom_stage_composition_matches_strang_phase_space_evolution()
{
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4), 2);
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity(
        [](const double time, const std::vector<double>& xi) { return std::vector<double>{0.35 + 0.1 * time - 0.05 * xi[0]}; });
    const auto source = [](const double time, const qmomkit::QmomSourceContext& moments)
    {
        std::vector<double> rhs(moments.size(), 0.0);
        for (std::size_t order = 1; order < moments.size(); ++order)
        {
            rhs[order] = -(0.2 + 0.05 * time) * static_cast<double>(order) * moments.value(order);
        }
        return qmomkit::QmomMomentRate{std::move(rhs)};
    };
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options = tight_node_options();
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    constexpr double dt = 0.2;

    const qmomkit::QmomPhaseSpaceEvolutionResult high_level = qmomkit::advance_qmom_phase_space(state, 0.0, dt, velocity, source, options);
    const qmomkit::QmomPhaseSpaceEvolutionResult first = qmomkit::advance_qmom_advection(state, 0.0, 0.5 * dt, velocity, options);
    assert(first.ok());
    const qmomkit::QmomPhaseSpaceEvolutionResult moment_source = qmomkit::advance_qmom_moment_source(first.state, 0.0, dt, source, options);
    assert(moment_source.ok());
    const qmomkit::QmomPhaseSpaceEvolutionResult second = qmomkit::advance_qmom_advection(moment_source.state, 0.5 * dt, 0.5 * dt, velocity, options);

    assert_qmom_phase_space_evolution_results_match(second, high_level, 5.0e-10);
    assert(high_level.diagnostics.first_advection_diagnostics.step_count == first.diagnostics.first_advection_diagnostics.step_count);
    assert(high_level.diagnostics.moment_source_status == moment_source.diagnostics.moment_source_status);
    assert(high_level.diagnostics.second_advection_diagnostics.step_count == second.diagnostics.first_advection_diagnostics.step_count);
}

void test_public_qmom_stage_composition_matches_lie_phase_space_evolution()
{
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4), 2);
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity([](const double time, const std::vector<double>& xi)
                                                                                    { return std::vector<double>{0.2 + time - 0.1 * xi[0]}; });
    const auto source = [](const double time, const qmomkit::QmomSourceContext& moments)
    {
        std::vector<double> rhs(moments.size(), 0.0);
        for (std::size_t order = 1; order < moments.size(); ++order)
        {
            rhs[order] = -(0.1 + 0.02 * time) * static_cast<double>(order) * moments.value(order);
        }
        return qmomkit::QmomMomentRate{std::move(rhs)};
    };
    constexpr double dt = 0.2;

    qmomkit::QmomPhaseSpaceEvolutionOptions advsrc_options;
    advsrc_options.node_options = tight_node_options();
    advsrc_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;
    const qmomkit::QmomPhaseSpaceEvolutionResult advsrc_high = qmomkit::advance_qmom_phase_space(state, 0.0, dt, velocity, source, advsrc_options);
    const qmomkit::QmomPhaseSpaceEvolutionResult advsrc_advection = qmomkit::advance_qmom_advection(state, 0.0, dt, velocity, advsrc_options);
    assert(advsrc_advection.ok());
    const qmomkit::QmomPhaseSpaceEvolutionResult advsrc_explicit =
        qmomkit::advance_qmom_moment_source(advsrc_advection.state, 0.0, dt, source, advsrc_options);
    assert_qmom_phase_space_evolution_results_match(advsrc_explicit, advsrc_high, 5.0e-10);

    qmomkit::QmomPhaseSpaceEvolutionOptions srcadv_options = advsrc_options;
    srcadv_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    const qmomkit::QmomPhaseSpaceEvolutionResult srcadv_high = qmomkit::advance_qmom_phase_space(state, 0.0, dt, velocity, source, srcadv_options);
    const qmomkit::QmomPhaseSpaceEvolutionResult srcadv_moment_source = qmomkit::advance_qmom_moment_source(state, 0.0, dt, source, srcadv_options);
    assert(srcadv_moment_source.ok());
    const qmomkit::QmomPhaseSpaceEvolutionResult srcadv_explicit =
        qmomkit::advance_qmom_advection(srcadv_moment_source.state, 0.0, dt, velocity, srcadv_options);
    assert_qmom_phase_space_evolution_results_match(srcadv_explicit, srcadv_high, 5.0e-10);
}

void test_qmom_phase_space_evolution_failure_after_first_stage_rolls_back_to_macro_input()
{
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4), 2);
    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{0.5}; });
    const auto bad_source = [](double, const qmomkit::QmomSourceContext& moments)
    { return qmomkit::QmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; };
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options = tight_node_options();
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;

    const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(state, 0.0, 0.2, velocity, bad_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    assert(result.diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success);
    assert_vector_near(result.state.moment_values(), state.moment_values(), 0.0);
    assert_qmom_quadrature_near(result.state.quadrature(), state.quadrature(), 0.0);
}

void test_qmom_macro_rollbacks_preserve_the_original_state_after_later_successes()
{
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4), 2);
    const auto zero_source = [](double, const qmomkit::QmomSourceContext& moments)
    { return qmomkit::QmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::QmomPhaseSpaceEvolutionOptions srcadv_options;
    srcadv_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    srcadv_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::PhaseSpaceVelocity failing_velocity = qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>&) -> std::vector<double> { throw std::runtime_error("advection failure after moment-source success"); });
    const qmomkit::QmomPhaseSpaceEvolutionResult srcadv =
        qmomkit::advance_qmom_phase_space(state, 0.0, 0.2, failing_velocity, zero_source, srcadv_options);

    assert(!srcadv.ok());
    assert(srcadv.status == qmomkit::PhaseSpaceEvolutionStatus::advection_failed);
    assert(srcadv.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(srcadv.diagnostics.second_advection_status == qmomkit::NodeAdvectionStatus::velocity_failed);
    // The srcadv scheme's single advection substep is the "second" stage, not the "first".
    assert(srcadv.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::second_advection_half_step);
    assert_vector_near(srcadv.state.moment_values(), state.moment_values(), 0.0);
    assert_qmom_quadrature_near(srcadv.state.quadrature(), state.quadrature(), 0.0);

    qmomkit::QmomPhaseSpaceEvolutionOptions strang_options;
    strang_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    strang_options.node_options = tight_node_options();
    std::size_t first_stage_velocity_calls = 0;
    const qmomkit::PhaseSpaceVelocity counting_velocity = qmomkit::make_phase_space_velocity(
        [&first_stage_velocity_calls](double, const std::vector<double>& xi)
        {
            ++first_stage_velocity_calls;
            return std::vector<double>(xi.size(), 0.2);
        });
    const qmomkit::QmomPhaseSpaceEvolutionResult first = qmomkit::advance_qmom_advection(state, 0.0, 0.1, counting_velocity, strang_options);
    assert(first.ok());
    assert(first_stage_velocity_calls > 0);

    std::size_t macro_velocity_calls = 0;
    const qmomkit::PhaseSpaceVelocity fail_final_advection = qmomkit::make_phase_space_velocity(
        [&macro_velocity_calls, first_stage_velocity_calls](double, const std::vector<double>& xi)
        {
            ++macro_velocity_calls;
            if (macro_velocity_calls > first_stage_velocity_calls)
            {
                throw std::runtime_error("final advection failure");
            }
            return std::vector<double>(xi.size(), 0.2);
        });
    const qmomkit::QmomPhaseSpaceEvolutionResult strang =
        qmomkit::advance_qmom_phase_space(state, 0.0, 0.2, fail_final_advection, zero_source, strang_options);

    assert(!strang.ok());
    assert(strang.status == qmomkit::PhaseSpaceEvolutionStatus::advection_failed);
    assert(strang.diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success);
    assert(strang.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(strang.diagnostics.second_advection_status == qmomkit::NodeAdvectionStatus::velocity_failed);
    assert(strang.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::second_advection_half_step);
    assert_vector_near(strang.state.moment_values(), state.moment_values(), 0.0);
    assert_qmom_quadrature_near(strang.state.quadrature(), state.quadrature(), 0.0);
}

// A failed step returns its input state, so final_time must be the input time
// whichever stage failed. The injected failures include stages whose own time
// differs from t: the second Strang half-step starts at t + dt/2, and moment
// reconstruction follows advection to the end of the advected interval.
void test_qmom_failed_steps_report_the_time_of_the_returned_state()
{
    constexpr double t0 = 1.0;
    constexpr double dt = 0.2;
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4), 2);

    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    // Fails only after t0 + dt/2, so in a Strang step only the second half-step fails.
    const qmomkit::PhaseSpaceVelocity late_failing_velocity = qmomkit::make_phase_space_velocity(
        [](const double time, const std::vector<double>& xi)
        {
            if (time > t0 + 0.6 * dt)
            {
                throw std::runtime_error("late advection failure");
            }
            return std::vector<double>(xi.size(), 0.2);
        });
    // Moves the nodes so far that their raw moments overflow, so the moment
    // reconstruction after the first advection stage fails.
    const qmomkit::PhaseSpaceVelocity overflowing_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 1.0e300); });
    const qmomkit::QmomMomentSource zero_source = [](double, const qmomkit::QmomSourceContext& moments)
    { return qmomkit::QmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    const qmomkit::QmomMomentSource non_finite_source = [](double, const qmomkit::QmomSourceContext& moments)
    {
        std::vector<double> rates(moments.size(), 0.0);
        rates.front() = std::numeric_limits<double>::quiet_NaN();
        return qmomkit::QmomMomentRate{std::move(rates)};
    };
    // Drives m0 negative over the step, so the post-source inversion fails.
    const qmomkit::QmomMomentSource negative_mass_source = [](double, const qmomkit::QmomSourceContext& moments)
    {
        std::vector<double> rates(moments.size(), 0.0);
        rates.front() = -100.0;
        return qmomkit::QmomMomentRate{std::move(rates)};
    };

    struct FailureCase
    {
        qmomkit::PhaseSpaceVelocity velocity;
        qmomkit::QmomMomentSource source;
        qmomkit::PhaseSpaceEvolutionStatus status;
        qmomkit::PhaseSpaceEvolutionStage strang_stage;
    };
    const std::vector<FailureCase> failure_cases{{late_failing_velocity, zero_source, qmomkit::PhaseSpaceEvolutionStatus::advection_failed,
                                                  qmomkit::PhaseSpaceEvolutionStage::second_advection_half_step},
                                                 {overflowing_velocity, zero_source, qmomkit::PhaseSpaceEvolutionStatus::moment_reconstruction_failed,
                                                  qmomkit::PhaseSpaceEvolutionStage::first_moment_reconstruction},
                                                 {zero_velocity, non_finite_source, qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed,
                                                  qmomkit::PhaseSpaceEvolutionStage::moment_source_step},
                                                 {zero_velocity, negative_mass_source, qmomkit::PhaseSpaceEvolutionStatus::inversion_failed,
                                                  qmomkit::PhaseSpaceEvolutionStage::post_source_inversion}};

    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options = tight_node_options();
    for (const qmomkit::PhaseSpaceEvolutionScheme scheme :
         {qmomkit::PhaseSpaceEvolutionScheme::strang, qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv})
    {
        options.scheme = scheme;
        for (const FailureCase& failure_case : failure_cases)
        {
            const qmomkit::QmomPhaseSpaceEvolutionResult result =
                qmomkit::advance_qmom_phase_space(state, t0, dt, failure_case.velocity, failure_case.source, options);
            assert(!result.ok());
            assert(result.status == failure_case.status);
            if (scheme == qmomkit::PhaseSpaceEvolutionScheme::strang)
            {
                assert(result.diagnostics.failed_stage == failure_case.strang_stage);
            }
            assert(result.final_time == t0);
            assert_vector_near(result.state.moment_values(), state.moment_values(), 0.0);
        }
    }

    // A standalone stage also returns its input state on failure.
    const qmomkit::QmomPhaseSpaceEvolutionResult advection = qmomkit::advance_qmom_advection(state, t0, dt, overflowing_velocity, options);
    assert(advection.status == qmomkit::PhaseSpaceEvolutionStatus::moment_reconstruction_failed);
    assert(advection.final_time == t0);
    assert_vector_near(advection.state.moment_values(), state.moment_values(), 0.0);

    const qmomkit::QmomPhaseSpaceEvolutionResult source = qmomkit::advance_qmom_moment_source(state, t0, dt, negative_mass_source, options);
    assert(source.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    assert(source.final_time == t0);
    assert_vector_near(source.state.moment_values(), state.moment_values(), 0.0);

    // A successful step reports the end of the step.
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    const qmomkit::QmomPhaseSpaceEvolutionResult success = qmomkit::advance_qmom_phase_space(state, t0, dt, zero_velocity, zero_source, options);
    assert(success.ok());
    assert(success.final_time == t0 + dt);
}

void test_cqmom_initial_state_stores_accepted_quadrature()
{
    assert_cqmom_initial_state_stores_accepted_quadrature(fixed_cqmom_prepared_case(), 5.0e-10);
    assert_cqmom_initial_state_stores_accepted_quadrature(adaptive_cqmom_prepared_case(), 5.0e-10);
}

void test_cqmom_initial_state_canonicalizes_adaptively_reduced_moments()
{
    const qmomkit::RawMomentTensor moments = near_collapsed_adaptive_cqmom_moments();
    constexpr double tolerance = 1.0e-30;

    {
        const qmomkit::CqmomState initial = adaptive_cqmom_prepared_case_from_moments(moments).initial_state();
        const std::optional<qmomkit::CqmomQuadrature>& accepted = qmomkit::internal::CqmomConstructionAccess::state_accepted_quadrature(initial);

        assert(accepted.has_value());
        assert(accepted->total_nodes() == 1);
        assert_cqmom_state_moments_match_quadrature(initial, *accepted, tolerance);

        std::size_t source_calls = 0;
        const qmomkit::CqmomMomentSource zero_source = [&source_calls](double, const qmomkit::CqmomSourceContext& context)
        {
            ++source_calls;
            assert_near(context.moment({3, 0}), context.generalized_moment({3.0, 0.0}), 1.0e-30);
            return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)};
        };

        const qmomkit::CqmomPhaseSpaceEvolutionResult advanced = qmomkit::advance_cqmom_moment_source(initial, 0.0, 0.1, zero_source);
        assert(advanced.ok());
        assert(!advanced.diagnostic_quadrature.has_value());
        assert(source_calls == 3);
        assert_vector_near(advanced.state.moment_values(), initial.moment_values(), tolerance);
        assert_cqmom_quadrature_near(advanced.state.quadrature(), *accepted, 1.0e-20);
        assert_cqmom_state_moments_match_quadrature(advanced.state, advanced.state.quadrature(), tolerance);
    }
}

void test_cqmom_first_advection_stage_uses_initial_accepted_quadrature()
{
    const qmomkit::CqmomState state = adaptive_cqmom_prepared_case_from_moments(branch_specific_ragged_cqmom_moments()).initial_state();
    const std::optional<qmomkit::CqmomQuadrature>& initial_quadrature = qmomkit::internal::CqmomConstructionAccess::state_accepted_quadrature(state);
    assert(initial_quadrature.has_value());
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.node_options = tight_node_options();
    const qmomkit::CqmomPhaseSpaceEvolutionResult stage = qmomkit::advance_cqmom_advection(state, 0.0, 0.1, zero_velocity, options);

    assert(stage.ok());
    assert(!stage.diagnostic_quadrature.has_value());
    assert(stage.diagnostics.initial_cqmom_diagnostics.empty());
    assert(stage.state.quadrature().terminal_paths() == initial_quadrature->terminal_paths());
    assert_cqmom_quadrature_near(stage.state.quadrature(), *initial_quadrature, 5.0e-10);
    assert(active_counts_equal(stage.state.active_branch_counts(), state.active_branch_counts()));
}

void assert_cqmom_phase_space_evolution_results_match(const qmomkit::CqmomPhaseSpaceEvolutionResult& actual,
                                                      const qmomkit::CqmomPhaseSpaceEvolutionResult& expected, const double tolerance)
{
    assert(actual.ok());
    assert(expected.ok());
    assert_near(actual.final_time, expected.final_time, tolerance);
    assert(!actual.diagnostic_quadrature.has_value());
    assert(!expected.diagnostic_quadrature.has_value());
    assert_cqmom_quadrature_near(actual.state.quadrature(), expected.state.quadrature(), tolerance);
    assert(active_counts_equal(actual.state.active_branch_counts(), expected.state.active_branch_counts()));
    assert_cqmom_moments_near(actual.state, expected.state, tolerance);
}

// A successful result carries no diagnostic quadrature: its accepted
// quadrature is the returned state's own, and that is the only copy.
void assert_cqmom_success_has_no_diagnostic_quadrature(const qmomkit::CqmomPhaseSpaceEvolutionResult& result)
{
    assert(result.ok());
    assert(!result.diagnostic_quadrature.has_value());
    const std::optional<qmomkit::CqmomQuadrature>& accepted = qmomkit::internal::CqmomConstructionAccess::state_accepted_quadrature(result.state);
    assert(accepted.has_value());
    assert(&result.state.quadrature() == &*accepted);
    assert_cqmom_state_active_counts_match_accepted_quadrature(result.state);
}

void assert_cqmom_states_equal(const qmomkit::CqmomState& actual, const qmomkit::CqmomState& expected)
{
    assert(actual.hierarchy().node_counts() == expected.hierarchy().node_counts());
    assert(actual.hierarchy().conditioning_order() == expected.hierarchy().conditioning_order());
    assert(actual.moment_indices() == expected.moment_indices());
    assert_vector_near(actual.moment_values(), expected.moment_values(), 0.0);
    assert(active_counts_equal(actual.active_branch_counts(), expected.active_branch_counts()));

    const std::optional<qmomkit::CqmomQuadrature>& actual_quadrature = qmomkit::internal::CqmomConstructionAccess::state_accepted_quadrature(actual);
    const std::optional<qmomkit::CqmomQuadrature>& expected_quadrature =
        qmomkit::internal::CqmomConstructionAccess::state_accepted_quadrature(expected);
    assert(actual_quadrature.has_value() == expected_quadrature.has_value());
    if (actual_quadrature.has_value())
    {
        assert_cqmom_quadrature_near(*actual_quadrature, *expected_quadrature, 0.0);
    }
}

// A state created from the open case must carry every setting of that case,
// because the case is the only place the caller set them.
void assert_cqmom_state_carries_case_settings(const qmomkit::CqmomState& state, const qmomkit::CqmomCase& cqmom_case)
{
    assert(state.qmom_options().adaptive == cqmom_case.qmom_options.adaptive);
    assert(state.qmom_options().rmin == cqmom_case.qmom_options.rmin);
    assert(state.qmom_options().eabs == cqmom_case.qmom_options.eabs);
    assert(state.qmom_options().cutoff == cqmom_case.qmom_options.cutoff);
    assert(state.qmom_options().automatic_scaling == cqmom_case.qmom_options.automatic_scaling);
    assert(state.conditional_options().singular_tolerance == cqmom_case.conditional_options.singular_tolerance);
    assert(state.conditional_options().conditioning_warning_threshold == cqmom_case.conditional_options.conditioning_warning_threshold);
    assert(state.conditional_options().moore_penrose_rcond == cqmom_case.conditional_options.moore_penrose_rcond);
    assert(state.adaptive_branch_counts() == cqmom_case.adaptive_branch_counts);
    assert(state.diagnostic_level() == cqmom_case.diagnostic_level);
    assert(state.vandermonde_solver() == cqmom_case.vandermonde_solver);
    assert(state.abscissa_bounds().size() == cqmom_case.abscissa_bounds.size());
    for (std::size_t dimension = 0; dimension < cqmom_case.abscissa_bounds.size(); ++dimension)
    {
        const std::optional<qmomkit::AbscissaBounds>& actual = state.abscissa_bounds()[dimension];
        const std::optional<qmomkit::AbscissaBounds>& expected = cqmom_case.abscissa_bounds[dimension];
        assert(actual.has_value() == expected.has_value());
        if (expected.has_value())
        {
            assert(actual->lower == expected->lower);
            assert(actual->upper == expected->upper);
        }
    }
}

// The open case and CqmomPreparedCase build their initial state the same way;
// they differ only in when the branch structure is built. The two states must
// therefore be identical, down to the last bit.
void assert_cqmom_case_initial_state_matches_prepared_case(const qmomkit::CqmomCase& cqmom_case)
{
    const qmomkit::CqmomState state = cqmom_case.initial_state();
    assert_cqmom_state_carries_case_settings(state, cqmom_case);
    assert_cqmom_states_equal(state, qmomkit::CqmomPreparedCase(cqmom_case).initial_state());
}

void test_cqmom_case_initial_state_matches_prepared_case()
{
    // Every setting except the two adaptive flags is moved away from its
    // default, so a setting the state failed to carry would show. The adaptive
    // cases below set both flags.
    qmomkit::CqmomCase fixed_cqmom_case = qmomkit::internal::CqmomConstructionAccess::case_definition(fixed_cqmom_prepared_case());
    fixed_cqmom_case.conditional_options.singular_tolerance = 1.0e-13;
    fixed_cqmom_case.conditional_options.conditioning_warning_threshold = 1.0e11;
    fixed_cqmom_case.conditional_options.moore_penrose_rcond = 1.0e-13;
    fixed_cqmom_case.qmom_options.rmin = 1.0e-11;
    fixed_cqmom_case.qmom_options.eabs = 1.0e-11;
    fixed_cqmom_case.qmom_options.cutoff = 1.0e-31;
    fixed_cqmom_case.qmom_options.automatic_scaling = true;
    fixed_cqmom_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::basic;
    fixed_cqmom_case.vandermonde_solver = qmomkit::VandermondeSolver::pivoted_qr;
    fixed_cqmom_case.abscissa_bounds = qmomkit::CqmomAbscissaBounds(fixed_cqmom_case.moments.layout.hierarchy.number_of_dimensions());
    fixed_cqmom_case.abscissa_bounds[0] = qmomkit::AbscissaBounds{-1.0e3, 1.0e3};
    assert_cqmom_case_initial_state_matches_prepared_case(fixed_cqmom_case);

    assert_cqmom_case_initial_state_matches_prepared_case(
        qmomkit::internal::CqmomConstructionAccess::case_definition(adaptive_cqmom_prepared_case()));

    // Adaptive reconstruction keeps one node here, so the stored moments are
    // replaced by those of the reduced quadrature.
    assert_cqmom_case_initial_state_matches_prepared_case(qmomkit::internal::CqmomConstructionAccess::case_definition(
        adaptive_cqmom_prepared_case_from_moments(near_collapsed_adaptive_cqmom_moments())));
}

template <class Exception, class Function> std::string thrown_message(const Function& function)
{
    try
    {
        function();
    }
    catch (const Exception& error)
    {
        return error.what();
    }
    assert(false && "expected an exception");
    return {};
}

// Both routes must throw the same exception type with the same message. The
// message is returned so that the caller can check what it says.
template <class Exception> std::string initial_state_failure_message(const qmomkit::CqmomCase& cqmom_case)
{
    const std::string message = thrown_message<Exception>([&] { static_cast<void>(cqmom_case.initial_state()); });
    assert(message == thrown_message<Exception>([&] { static_cast<void>(qmomkit::CqmomPreparedCase(cqmom_case).initial_state()); }));
    return message;
}

// The open case reports failures exactly as CqmomPreparedCase does.
void test_cqmom_case_initial_state_reports_failures()
{
    const std::string reconstruction_prefix = "CQMOM initial state quadrature reconstruction failed: ";

    qmomkit::RawMomentTensor zero_moments = adaptive_cqmom_moments();
    std::fill(zero_moments.values.begin(), zero_moments.values.end(), 0.0);
    const qmomkit::CqmomCase zero_cqmom_case =
        qmomkit::internal::CqmomConstructionAccess::case_definition(adaptive_cqmom_prepared_case_from_moments(zero_moments));
    assert(initial_state_failure_message<std::runtime_error>(zero_cqmom_case).rfind(reconstruction_prefix, 0) == 0);

    // The moment values of an open case can be changed after it is built, so
    // the values are checked again when the state is created.
    const qmomkit::CqmomCase valid_cqmom_case = qmomkit::internal::CqmomConstructionAccess::case_definition(fixed_cqmom_prepared_case());
    qmomkit::CqmomCase nonfinite_cqmom_case = valid_cqmom_case;
    nonfinite_cqmom_case.moments.values[3] = std::numeric_limits<double>::quiet_NaN();
    assert(initial_state_failure_message<std::runtime_error>(nonfinite_cqmom_case) ==
           reconstruction_prefix + "CQMOM moment data values must be finite.");

    qmomkit::CqmomCase short_cqmom_case = valid_cqmom_case;
    short_cqmom_case.moments.values.pop_back();
    assert(initial_state_failure_message<std::runtime_error>(short_cqmom_case) ==
           reconstruction_prefix + "CQMOM moment data values must match the layout size.");

    // A layout listing the required indices in another order is rejected before
    // any reconstruction.
    const qmomkit::CqmomMomentLayout fixed_layout = qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy({2, 2}));
    std::vector<qmomkit::MultiIndex> reordered_indices = fixed_layout.required_indices.to_multi_indices();
    std::reverse(reordered_indices.begin(), reordered_indices.end());
    const qmomkit::CqmomCase reordered_cqmom_case{
        qmomkit::CqmomMomentData{qmomkit::CqmomMomentLayout{fixed_layout.hierarchy, qmomkit::MultiIndexTable(reordered_indices)},
                                 std::vector<double>(fixed_layout.required_indices.size(), 1.0)}};
    assert(initial_state_failure_message<std::invalid_argument>(reordered_cqmom_case).find("canonical CQMOM moment layout") != std::string::npos);
}

void test_public_cqmom_advection_stage_preserves_topology()
{
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = adaptive_cqmom_prepared_case_from_moments(branch_specific_ragged_cqmom_moments());
    const qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_prepared_case);
    assert(initial.ok());
    const std::optional<qmomkit::CqmomQuadrature>& input_quadrature = qmomkit::internal::CqmomConstructionAccess::state_accepted_quadrature(state);
    assert(input_quadrature.has_value());
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& xi)
        {
            assert(xi.size() == 2);
            return std::vector<double>{0.1 + 0.05 * xi[0], -0.08 * xi[1]};
        });

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.node_options = tight_node_options();
    const qmomkit::CqmomPhaseSpaceEvolutionResult stage = qmomkit::advance_cqmom_advection(state, 0.0, 0.2, velocity, options);

    assert(stage.ok());
    assert(!stage.diagnostic_quadrature.has_value());
    assert(stage.final_time == 0.2);
    assert(stage.state.quadrature().hierarchy().node_counts() == initial.quadrature.hierarchy().node_counts());
    assert(stage.state.quadrature().hierarchy().conditioning_order() == initial.quadrature.hierarchy().conditioning_order());
    assert(stage.state.quadrature().terminal_paths() == initial.quadrature.terminal_paths());
    assert_vector_near(stage.state.quadrature().weights(), initial.quadrature.weights(), 0.0);
    assert(active_counts_equal(stage.state.active_branch_counts(), state.active_branch_counts()));
    assert_cqmom_state_active_counts_match_accepted_quadrature(stage.state);
    bool abscissa_changed = false;
    for (std::size_t dimension = 0; dimension < input_quadrature->abscissas().size(); ++dimension)
    {
        for (std::size_t node = 0; node < input_quadrature->abscissas()[dimension].size(); ++node)
        {
            abscissa_changed = abscissa_changed || std::abs(stage.state.quadrature().abscissas()[dimension][node] -
                                                            input_quadrature->abscissas()[dimension][node]) > 1.0e-12;
        }
    }
    assert(abscissa_changed);
    assert(stage.diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success);
    assert(stage.diagnostics.first_reconstruction_diagnostics.success);
}

void test_public_cqmom_moment_source_stage_adaptive_policy()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor retained_moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData retained_required = qmomkit::extract_cqmom_moments(retained_moments, layout);
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    const qmomkit::CqmomState retained_state = adaptive_cqmom_prepared_case_from_moments(retained_moments).initial_state();
    const qmomkit::CqmomPhaseSpaceEvolutionResult retained = qmomkit::advance_cqmom_moment_source(retained_state, 0.0, 0.1, zero_source, options);
    assert(retained.ok());
    assert(!retained.diagnostic_quadrature.has_value());
    assert(active_counts_equal(retained.state.active_branch_counts(), retained_state.active_branch_counts()));
    assert_cqmom_state_moments_match_quadrature(retained.state, retained.state.quadrature(), 5.0e-10);

    const qmomkit::RawMomentTensor reduced_moments = collapsed_branch_cqmom_moments();
    const qmomkit::CqmomMomentData reduced_required = qmomkit::extract_cqmom_moments(reduced_moments, layout);
    const qmomkit::CqmomPhaseSpaceEvolutionResult reduced =
        qmomkit::advance_cqmom_moment_source(retained_state, 0.0, 1.0, source_from_to(retained_required, reduced_required), options);
    assert(reduced.ok());
    assert(!reduced.diagnostic_quadrature.has_value());
    assert(reduced.state.quadrature().total_nodes() == 2);
    assert_cqmom_state_moments_match_quadrature(reduced.state, reduced.state.quadrature(), 5.0e-10);
    assert_cqmom_state_active_counts_match_accepted_quadrature(reduced.state);
    assert(transitions_contain(reduced.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{1}}, 2, 1, qmomkit::CqmomNodeCountTransitionKind::reduced));

    const qmomkit::CqmomState collapsed_state = adaptive_cqmom_prepared_case_from_moments(reduced_moments).initial_state();
    const qmomkit::CqmomPhaseSpaceEvolutionResult growth =
        qmomkit::advance_cqmom_moment_source(collapsed_state, 0.0, 1.0, source_from_to(reduced_required, retained_required), options);
    assert(growth.ok());
    assert(!growth.diagnostic_quadrature.has_value());
    assert(growth.state.quadrature().total_nodes() == 3);
    assert_cqmom_state_moments_match_quadrature(growth.state, growth.state.quadrature(), 5.0e-10);
    assert_cqmom_state_active_counts_match_accepted_quadrature(growth.state);
    assert(transitions_contain(growth.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{1}}, 1, 2, qmomkit::CqmomNodeCountTransitionKind::grew));
}

void test_public_cqmom_moment_source_stage_growth_is_limited_to_one_node()
{
    constexpr std::size_t configured_node_count = 3;
    const qmomkit::RawMomentTensor initial_moments = one_dimensional_cqmom_moments(configured_node_count, {0.0}, {1.0});
    const qmomkit::RawMomentTensor target_moments = one_dimensional_cqmom_moments(configured_node_count, {-1.0, 0.25, 2.0}, {0.2, 0.3, 0.5});
    const qmomkit::PhaseSpaceHierarchy hierarchy({configured_node_count});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);

    const qmomkit::CqmomPhaseSpaceEvolutionResult stage = qmomkit::advance_cqmom_moment_source(
        one_dimensional_adaptive_cqmom_state(initial_required), 0.0, 1.0, source_from_to(initial_required, target_required));

    assert(stage.ok());
    assert(!stage.diagnostic_quadrature.has_value());
    assert(stage.state.quadrature().total_nodes() == 2);
    assert(active_counts_contain(stage.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, configured_node_count));
    assert_cqmom_state_moments_match_quadrature(stage.state, stage.state.quadrature(), 5.0e-10);
}

void test_public_cqmom_moment_source_stage_canonicalizes_adaptive_reduction()
{
    constexpr std::size_t configured_node_count = 3;
    const qmomkit::PhaseSpaceHierarchy hierarchy({configured_node_count});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::RawMomentTensor initial_moments = one_dimensional_cqmom_moments(configured_node_count, {-1.0, 1.0}, {0.35, 0.65});
    qmomkit::RawMomentTensor source_candidate_moments = one_dimensional_cqmom_moments(configured_node_count, {0.25}, {1.0});
    source_candidate_moments.values[4] = 0.5;
    source_candidate_moments.values[5] = -0.25;
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData source_candidate_required = qmomkit::extract_cqmom_moments(source_candidate_moments, layout);

    const qmomkit::CqmomPhaseSpaceEvolutionResult stage = qmomkit::advance_cqmom_moment_source(
        one_dimensional_adaptive_cqmom_state(initial_required), 0.0, 1.0, source_from_to(initial_required, source_candidate_required));

    assert(stage.ok());
    assert(!stage.diagnostic_quadrature.has_value());
    assert(stage.state.quadrature().total_nodes() == 1);
    assert(active_counts_contain(stage.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 1, 0, configured_node_count));
    assert(transitions_contain(stage.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion, 0,
                               qmomkit::ConditioningPath{}, 2, 1, qmomkit::CqmomNodeCountTransitionKind::reduced));
    assert_cqmom_state_moments_match_quadrature(stage.state, stage.state.quadrature(), 5.0e-12);

    const qmomkit::RawMomentTensor reconstructed =
        qmomkit::reconstruct_raw_moments(stage.state.quadrature(), stage.state.hierarchy().moment_extents());
    const qmomkit::CqmomMomentData canonical_required = qmomkit::extract_cqmom_moments(reconstructed, layout);
    assert_vector_near(stage.state.moment_values(), canonical_required.values, 5.0e-12);

    double largest_source_candidate_residual = 0.0;
    for (std::size_t i = 0; i < source_candidate_required.values.size(); ++i)
    {
        largest_source_candidate_residual =
            std::max(largest_source_candidate_residual, std::abs(source_candidate_required.values[i] - canonical_required.values[i]));
    }
    assert(largest_source_candidate_residual > 0.49);
}

void assert_public_cqmom_stage_composition_matches_phase_space_evolution(const qmomkit::PhaseSpaceEvolutionScheme scheme)
{
    const qmomkit::CqmomState state = fixed_cqmom_prepared_case().initial_state();
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity(
        [](const double time, const std::vector<double>& xi)
        {
            assert(xi.size() == 2);
            return std::vector<double>{0.2 + 0.03 * time - 0.04 * xi[0], -0.1 + 0.02 * time + 0.05 * xi[1]};
        });
    const qmomkit::CqmomMomentSource source = [](const double time, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            rate[i] = -(0.18 + 0.04 * time) * static_cast<double>(context.index(i).size()) * context.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };
    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.node_options = tight_node_options();
    options.scheme = scheme;
    constexpr double dt = 0.2;

    const qmomkit::CqmomPhaseSpaceEvolutionResult high_level = qmomkit::advance_cqmom_phase_space(state, 0.0, dt, velocity, source, options);
    qmomkit::CqmomPhaseSpaceEvolutionResult explicit_result = high_level;
    if (scheme == qmomkit::PhaseSpaceEvolutionScheme::strang)
    {
        const qmomkit::CqmomPhaseSpaceEvolutionResult first = qmomkit::advance_cqmom_advection(state, 0.0, 0.5 * dt, velocity, options);
        assert(first.ok());
        const qmomkit::CqmomPhaseSpaceEvolutionResult moment_source = qmomkit::advance_cqmom_moment_source(first.state, 0.0, dt, source, options);
        assert(moment_source.ok());
        explicit_result = qmomkit::advance_cqmom_advection(moment_source.state, 0.5 * dt, 0.5 * dt, velocity, options);
    }
    else if (scheme == qmomkit::PhaseSpaceEvolutionScheme::advsrc)
    {
        const qmomkit::CqmomPhaseSpaceEvolutionResult advection = qmomkit::advance_cqmom_advection(state, 0.0, dt, velocity, options);
        assert(advection.ok());
        explicit_result = qmomkit::advance_cqmom_moment_source(advection.state, 0.0, dt, source, options);
    }
    else
    {
        const qmomkit::CqmomPhaseSpaceEvolutionResult moment_source = qmomkit::advance_cqmom_moment_source(state, 0.0, dt, source, options);
        assert(moment_source.ok());
        explicit_result = qmomkit::advance_cqmom_advection(moment_source.state, 0.0, dt, velocity, options);
    }

    assert_cqmom_phase_space_evolution_results_match(explicit_result, high_level, 5.0e-9);
}

void test_public_cqmom_stage_composition_matches_phase_space_evolution()
{
    assert_public_cqmom_stage_composition_matches_phase_space_evolution(qmomkit::PhaseSpaceEvolutionScheme::strang);
    assert_public_cqmom_stage_composition_matches_phase_space_evolution(qmomkit::PhaseSpaceEvolutionScheme::advsrc);
    assert_public_cqmom_stage_composition_matches_phase_space_evolution(qmomkit::PhaseSpaceEvolutionScheme::srcadv);
}

void test_cqmom_phase_space_evolution_success_returns_no_diagnostic_quadrature()
{
    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.1); });
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    {
        const qmomkit::CqmomState state = fixed_cqmom_prepared_case().initial_state();
        qmomkit::CqmomPhaseSpaceEvolutionOptions options;
        options.node_options = tight_node_options();

        assert_cqmom_success_has_no_diagnostic_quadrature(qmomkit::advance_cqmom_advection(state, 0.0, 0.1, velocity, options));
        assert_cqmom_success_has_no_diagnostic_quadrature(qmomkit::advance_cqmom_moment_source(state, 0.0, 0.1, zero_source, options));

        for (const qmomkit::PhaseSpaceEvolutionScheme scheme :
             {qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv, qmomkit::PhaseSpaceEvolutionScheme::strang})
        {
            options.scheme = scheme;
            assert_cqmom_success_has_no_diagnostic_quadrature(qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, velocity, zero_source, options));
        }
    }

    const qmomkit::CqmomPreparedCase cqmom_prepared_case = fixed_cqmom_prepared_case();
    const qmomkit::CqmomState cache_miss_state = qmomkit::internal::CqmomConstructionAccess::create_state(
        qmomkit::internal::CqmomConstructionAccess::case_definition(cqmom_prepared_case),
        qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(cqmom_prepared_case));
    assert(!qmomkit::internal::CqmomConstructionAccess::state_accepted_quadrature(cache_miss_state).has_value());
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(cache_miss_state).get() ==
           qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(cqmom_prepared_case).get());
    assert_cqmom_success_has_no_diagnostic_quadrature(qmomkit::advance_cqmom_moment_source(cache_miss_state, 0.0, 0.1, zero_source));
}

qmomkit::CqmomState adaptive_cqmom_state_without_accepted_quadrature(const qmomkit::CqmomState& state)
{
    return qmomkit::internal::CqmomConstructionAccess::create_state(qmomkit::internal::CqmomConstructionAccess::state_definition(state),
                                                                    state.active_branch_counts());
}

void test_adaptive_cqmom_phase_space_evolution_state_counts_match_accepted_quadrature()
{
    const qmomkit::RawMomentTensor moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& coordinate)
        {
            assert(coordinate.size() == 2);
            return std::vector<double>{0.05 + 0.01 * coordinate[0], -0.02 + 0.03 * coordinate[1]};
        });
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    {
        const qmomkit::CqmomState state = adaptive_cqmom_prepared_case_from_moments(moments).initial_state();
        assert_cqmom_state_active_counts_match_accepted_quadrature(state);

        qmomkit::CqmomPhaseSpaceEvolutionOptions options;
        options.node_options = tight_node_options();
        const qmomkit::CqmomPhaseSpaceEvolutionResult advection = qmomkit::advance_cqmom_advection(state, 0.0, 0.1, velocity, options);
        const qmomkit::CqmomPhaseSpaceEvolutionResult moment_source = qmomkit::advance_cqmom_moment_source(state, 0.0, 0.1, zero_source, options);
        assert(advection.ok());
        assert(moment_source.ok());
        assert_cqmom_state_active_counts_match_accepted_quadrature(advection.state);
        assert_cqmom_state_active_counts_match_accepted_quadrature(moment_source.state);

        for (const qmomkit::PhaseSpaceEvolutionScheme scheme :
             {qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv, qmomkit::PhaseSpaceEvolutionScheme::strang})
        {
            options.scheme = scheme;
            const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolution =
                qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, velocity, zero_source, options);
            assert(phase_space_evolution.ok());
            assert_cqmom_state_active_counts_match_accepted_quadrature(phase_space_evolution.state);
        }

        const qmomkit::CqmomState cache_miss_state = adaptive_cqmom_state_without_accepted_quadrature(state);
        assert(!qmomkit::internal::CqmomConstructionAccess::state_accepted_quadrature(cache_miss_state).has_value());

        const qmomkit::CqmomPhaseSpaceEvolutionResult cache_miss_advection =
            qmomkit::advance_cqmom_advection(cache_miss_state, 0.0, 0.1, velocity, options);
        const qmomkit::CqmomPhaseSpaceEvolutionResult cache_miss_moment_source =
            qmomkit::advance_cqmom_moment_source(cache_miss_state, 0.0, 0.1, zero_source, options);
        assert(cache_miss_advection.ok());
        assert(cache_miss_moment_source.ok());
        assert_cqmom_state_active_counts_match_accepted_quadrature(cache_miss_advection.state);
        assert_cqmom_state_active_counts_match_accepted_quadrature(cache_miss_moment_source.state);

        for (const qmomkit::PhaseSpaceEvolutionScheme scheme :
             {qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv, qmomkit::PhaseSpaceEvolutionScheme::strang})
        {
            options.scheme = scheme;
            const qmomkit::CqmomState cache_miss_scheme_state = adaptive_cqmom_state_without_accepted_quadrature(state);
            const qmomkit::CqmomPhaseSpaceEvolutionResult cache_miss_scheme =
                qmomkit::advance_cqmom_phase_space(cache_miss_scheme_state, 0.0, 0.1, velocity, zero_source, options);
            assert(cache_miss_scheme.ok());
            assert(!cache_miss_scheme.diagnostic_quadrature.has_value());
            assert(cache_miss_scheme.state.quadrature().adaptive_branch_counts_enabled());
            assert(!cache_miss_scheme.state.quadrature().terminal_paths().empty());
            assert_cqmom_state_active_counts_match_accepted_quadrature(cache_miss_scheme.state);
        }
    }
}

void test_public_cqmom_stage_failures_roll_back()
{
    const qmomkit::CqmomState state = fixed_cqmom_prepared_case().initial_state();
    const qmomkit::PhaseSpaceVelocity bad_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{1.0}; });
    const qmomkit::CqmomPhaseSpaceEvolutionResult advection = qmomkit::advance_cqmom_advection(state, 0.0, 0.1, bad_velocity);
    assert(!advection.ok());
    assert(advection.status == qmomkit::PhaseSpaceEvolutionStatus::advection_failed);
    assert_vector_near(advection.state.moment_values(), state.moment_values(), 0.0);
    assert(active_counts_equal(advection.state.active_branch_counts(), state.active_branch_counts()));

    const auto bad_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size() + 1, 0.0)}; };
    const qmomkit::CqmomPhaseSpaceEvolutionResult moment_source = qmomkit::advance_cqmom_moment_source(state, 0.0, 0.1, bad_source);
    assert(!moment_source.ok());
    assert(moment_source.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert_vector_near(moment_source.state.moment_values(), state.moment_values(), 0.0);
    assert(active_counts_equal(moment_source.state.active_branch_counts(), state.active_branch_counts()));
}

void test_cqmom_phase_space_evolution_failure_after_first_stage_rolls_back_to_macro_input()
{
    const qmomkit::CqmomState state = fixed_cqmom_prepared_case().initial_state();
    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.2); });
    const auto bad_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size() + 1, 0.0)}; };
    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.node_options = tight_node_options();
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;

    const qmomkit::CqmomPhaseSpaceEvolutionResult result = qmomkit::advance_cqmom_phase_space(state, 0.0, 0.2, velocity, bad_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    assert(result.diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success);
    assert_cqmom_states_equal(result.state, state);
}

void test_cqmom_macro_rollbacks_preserve_the_original_state_after_later_successes()
{
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    {
        const qmomkit::CqmomState state = fixed_cqmom_prepared_case().initial_state();

        qmomkit::CqmomPhaseSpaceEvolutionOptions srcadv_options;
        srcadv_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
        srcadv_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
        const qmomkit::PhaseSpaceVelocity failing_velocity =
            qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) -> std::vector<double>
                                               { throw std::runtime_error("advection failure after moment-source success"); });
        const qmomkit::CqmomPhaseSpaceEvolutionResult srcadv =
            qmomkit::advance_cqmom_phase_space(state, 0.0, 0.2, failing_velocity, zero_source, srcadv_options);

        assert(!srcadv.ok());
        assert(srcadv.status == qmomkit::PhaseSpaceEvolutionStatus::advection_failed);
        assert(srcadv.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);
        assert(srcadv.diagnostics.second_advection_status == qmomkit::NodeAdvectionStatus::velocity_failed);
        // The srcadv scheme's single advection substep is the "second" stage, not the "first".
        assert(srcadv.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::second_advection_half_step);
        assert_cqmom_states_equal(srcadv.state, state);

        qmomkit::CqmomPhaseSpaceEvolutionOptions strang_options;
        strang_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
        strang_options.node_options = tight_node_options();
        std::size_t first_stage_velocity_calls = 0;
        const qmomkit::PhaseSpaceVelocity counting_velocity = qmomkit::make_phase_space_velocity(
            [&first_stage_velocity_calls](double, const std::vector<double>& xi)
            {
                ++first_stage_velocity_calls;
                return std::vector<double>(xi.size(), 0.2);
            });
        const qmomkit::CqmomPhaseSpaceEvolutionResult first = qmomkit::advance_cqmom_advection(state, 0.0, 0.1, counting_velocity, strang_options);
        assert(first.ok());
        assert(first_stage_velocity_calls > 0);

        std::size_t macro_velocity_calls = 0;
        const qmomkit::PhaseSpaceVelocity fail_final_advection = qmomkit::make_phase_space_velocity(
            [&macro_velocity_calls, first_stage_velocity_calls](double, const std::vector<double>& xi)
            {
                ++macro_velocity_calls;
                if (macro_velocity_calls > first_stage_velocity_calls)
                {
                    throw std::runtime_error("final advection failure");
                }
                return std::vector<double>(xi.size(), 0.2);
            });
        const qmomkit::CqmomPhaseSpaceEvolutionResult strang =
            qmomkit::advance_cqmom_phase_space(state, 0.0, 0.2, fail_final_advection, zero_source, strang_options);

        assert(!strang.ok());
        assert(strang.status == qmomkit::PhaseSpaceEvolutionStatus::advection_failed);
        assert(strang.diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success);
        assert(strang.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);
        assert(strang.diagnostics.second_advection_status == qmomkit::NodeAdvectionStatus::velocity_failed);
        assert(strang.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::second_advection_half_step);
        assert_cqmom_states_equal(strang.state, state);
    }
}

// The CQMOM counterpart of test_qmom_failed_steps_report_the_time_of_the_returned_state.
void test_cqmom_failed_steps_report_the_time_of_the_returned_state()
{
    constexpr double t0 = 1.0;
    constexpr double dt = 0.2;
    const qmomkit::CqmomState state = fixed_cqmom_prepared_case().initial_state();

    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    // Fails only after t0 + dt/2, so in a Strang step only the second half-step fails.
    const qmomkit::PhaseSpaceVelocity late_failing_velocity = qmomkit::make_phase_space_velocity(
        [](const double time, const std::vector<double>& xi)
        {
            if (time > t0 + 0.6 * dt)
            {
                throw std::runtime_error("late advection failure");
            }
            return std::vector<double>(xi.size(), 0.2);
        });
    // Moves the nodes so far that their raw moments overflow, so the moment
    // reconstruction after the first advection stage fails.
    const qmomkit::PhaseSpaceVelocity overflowing_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 1.0e300); });
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };
    const qmomkit::CqmomMomentSource non_finite_source = [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rates(context.size(), 0.0);
        rates.front() = std::numeric_limits<double>::quiet_NaN();
        return qmomkit::CqmomMomentRate{std::move(rates)};
    };
    // Drives M[0,0] negative over the step, so the post-source inversion fails.
    const qmomkit::CqmomMomentSource negative_mass_source = [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rates(context.size(), 0.0);
        for (std::size_t index = 0; index < context.size(); ++index)
        {
            const qmomkit::MultiIndexView order = context.index(index);
            if (order[0] == 0 && order[1] == 0)
            {
                rates[index] = -100.0;
            }
        }
        return qmomkit::CqmomMomentRate{std::move(rates)};
    };

    struct FailureCase
    {
        qmomkit::PhaseSpaceVelocity velocity;
        qmomkit::CqmomMomentSource source;
        qmomkit::PhaseSpaceEvolutionStatus status;
        qmomkit::PhaseSpaceEvolutionStage strang_stage;
    };
    const std::vector<FailureCase> failure_cases{{late_failing_velocity, zero_source, qmomkit::PhaseSpaceEvolutionStatus::advection_failed,
                                                  qmomkit::PhaseSpaceEvolutionStage::second_advection_half_step},
                                                 {overflowing_velocity, zero_source, qmomkit::PhaseSpaceEvolutionStatus::moment_reconstruction_failed,
                                                  qmomkit::PhaseSpaceEvolutionStage::first_moment_reconstruction},
                                                 {zero_velocity, non_finite_source, qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed,
                                                  qmomkit::PhaseSpaceEvolutionStage::moment_source_step},
                                                 {zero_velocity, negative_mass_source, qmomkit::PhaseSpaceEvolutionStatus::inversion_failed,
                                                  qmomkit::PhaseSpaceEvolutionStage::post_source_inversion}};

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.node_options = tight_node_options();
    for (const qmomkit::PhaseSpaceEvolutionScheme scheme :
         {qmomkit::PhaseSpaceEvolutionScheme::strang, qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv})
    {
        options.scheme = scheme;
        for (const FailureCase& failure_case : failure_cases)
        {
            const qmomkit::CqmomPhaseSpaceEvolutionResult result =
                qmomkit::advance_cqmom_phase_space(state, t0, dt, failure_case.velocity, failure_case.source, options);
            assert(!result.ok());
            assert(result.status == failure_case.status);
            if (scheme == qmomkit::PhaseSpaceEvolutionScheme::strang)
            {
                assert(result.diagnostics.failed_stage == failure_case.strang_stage);
            }
            assert(result.final_time == t0);
            assert_vector_near(result.state.moment_values(), state.moment_values(), 0.0);
            // Every stage failure carries the quadrature the failing stage produced or started from.
            assert(result.diagnostic_quadrature.has_value());
        }
    }

    // A standalone stage also returns its input state on failure.
    const qmomkit::CqmomPhaseSpaceEvolutionResult advection = qmomkit::advance_cqmom_advection(state, t0, dt, overflowing_velocity, options);
    assert(advection.status == qmomkit::PhaseSpaceEvolutionStatus::moment_reconstruction_failed);
    assert(advection.final_time == t0);
    assert_vector_near(advection.state.moment_values(), state.moment_values(), 0.0);
    // The diagnostic quadrature is the advected one whose moments overflowed,
    // not the returned input state's.
    assert(advection.diagnostic_quadrature.has_value());
    assert(advection.diagnostic_quadrature->abscissas()[0][0] > 1.0e290);
    assert(advection.state.quadrature().abscissas()[0][0] < 1.0e3);

    const qmomkit::CqmomPhaseSpaceEvolutionResult source = qmomkit::advance_cqmom_moment_source(state, t0, dt, negative_mass_source, options);
    assert(source.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    assert(source.final_time == t0);
    assert_vector_near(source.state.moment_values(), state.moment_values(), 0.0);

    // Invalid input fails before any stage runs, so there is no diagnostic quadrature.
    const qmomkit::CqmomPhaseSpaceEvolutionResult invalid = qmomkit::advance_cqmom_phase_space(state, t0, -dt, zero_velocity, zero_source, options);
    assert(invalid.status == qmomkit::PhaseSpaceEvolutionStatus::invalid_input);
    assert(invalid.final_time == t0);
    assert(!invalid.diagnostic_quadrature.has_value());

    // A successful step reports the end of the step.
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    const qmomkit::CqmomPhaseSpaceEvolutionResult success = qmomkit::advance_cqmom_phase_space(state, t0, dt, zero_velocity, zero_source, options);
    assert(success.ok());
    assert(success.final_time == t0 + dt);
    assert(!success.diagnostic_quadrature.has_value());
}

void test_cqmom_macro_success_preserves_input_and_final_value_semantics()
{
    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.1); });
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    {
        const qmomkit::CqmomState state = fixed_cqmom_prepared_case().initial_state();
        const qmomkit::CqmomState original = state;
        qmomkit::CqmomPhaseSpaceEvolutionOptions options;
        options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
        options.node_options = tight_node_options();
        qmomkit::CqmomPhaseSpaceEvolutionResult result = qmomkit::advance_cqmom_phase_space(state, 0.0, 0.2, velocity, zero_source, options);

        assert(result.ok());
        assert_cqmom_states_equal(state, original);
        qmomkit::CqmomState copied_result = result.state;
        qmomkit::CqmomState moved_result = std::move(copied_result);
        assert_cqmom_states_equal(moved_result, result.state);

        const qmomkit::CqmomPhaseSpaceEvolutionResult next = qmomkit::advance_cqmom_advection(moved_result, 0.2, 0.1, velocity, options);
        assert(next.ok());
        assert_cqmom_states_equal(result.state, moved_result);
    }
}

void test_public_cqmom_boundary_moment_source_stage_matches_phase_space_evolution()
{
    const qmomkit::CqmomState state = fixed_cqmom_prepared_case().initial_state();
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const qmomkit::CqmomMomentSource source = [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            rate[i] = 0.10 * context.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };
    const qmomkit::CqmomMomentBoundaryContribution boundary = [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            rate[i] = 0.15 * context.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };
    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::CqmomPhaseSpaceEvolutionResult high_level =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, source, boundary, options);
    const qmomkit::CqmomPhaseSpaceEvolutionResult stage = qmomkit::advance_cqmom_moment_source(state, 0.0, 0.1, source, boundary, options);

    assert_cqmom_phase_space_evolution_results_match(stage, high_level, 5.0e-10);
    assert(stage.diagnostics.moment_source_status == high_level.diagnostics.moment_source_status);
}

void test_public_cqmom_stage_irrelevant_options_are_ignored()
{
    const qmomkit::CqmomState state = fixed_cqmom_prepared_case().initial_state();
    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    qmomkit::CqmomPhaseSpaceEvolutionOptions advection_options;
    advection_options.source_options.ssp_rk3_options.min_time_step = 1.0;
    advection_options.node_options.adaptive_options.initial_step = 0.1;
    const qmomkit::CqmomPhaseSpaceEvolutionResult advection = qmomkit::advance_cqmom_advection(state, 0.0, 0.1, velocity, advection_options);
    assert(advection.ok());

    qmomkit::CqmomPhaseSpaceEvolutionOptions moment_source_options;
    moment_source_options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    moment_source_options.node_options.adaptive_options.max_steps = 0;
    const qmomkit::CqmomPhaseSpaceEvolutionResult moment_source =
        qmomkit::advance_cqmom_moment_source(state, 0.0, 0.1, zero_source, moment_source_options);
    assert(moment_source.ok());
}

void test_advection_only_matches_direct_operator_qmom()
{
    const std::vector<double> moments = raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4);
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 2);

    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>{0.5 * xi[0]}; });
    const auto zero_source = [](double, const qmomkit::QmomSourceContext& state)
    { return qmomkit::QmomMomentRate{std::vector<double>(state.size(), 0.0)}; };

    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options = tight_node_options();
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    const qmomkit::QmomPhaseSpaceEvolutionResult phase_space_evolution =
        qmomkit::advance_qmom_phase_space(state, 0.0, 0.3, velocity, zero_source, options);
    const qmomkit::QmomNodeAdvectionResult direct = qmomkit::advance_nodes(state.quadrature(), 0.0, 0.3, velocity, options.node_options);

    assert(phase_space_evolution.ok());
    assert(direct.ok());
    const std::vector<double> expected = qmomkit::reconstruct_raw_moments(direct.quadrature, moments.size());
    assert_vector_near(phase_space_evolution.state.moment_values(), expected, 5.0e-9);
    assert_qmom_state_canonical(phase_space_evolution.state, 5.0e-9);
}

void test_source_only_matches_direct_operator_qmom()
{
    const std::vector<double> moments{1.0, 2.0};
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto source = [](double, const qmomkit::QmomSourceContext& state)
    {
        assert(state.size() == 2);
        return qmomkit::QmomMomentRate{std::vector<double>{0.0, 3.0}};
    };

    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::QmomPhaseSpaceEvolutionResult phase_space_evolution =
        qmomkit::advance_qmom_phase_space(state, 0.0, 0.25, zero_velocity, source, options);
    const qmomkit::QmomMomentSourceUpdateResult direct = qmomkit::advance_qmom_moments_with_source(moments, 0.125, 0.25, source);

    assert(phase_space_evolution.ok());
    assert(direct.ok());
    assert_vector_near(phase_space_evolution.state.moment_values(), direct.moments, 1.0e-13);
    assert_qmom_state_canonical(phase_space_evolution.state, 1.0e-13);
}

void test_source_plus_boundary_matches_direct_operator_qmom()
{
    const std::vector<double> moments{1.0, 2.0};
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const qmomkit::QmomMomentSource source = [](double, const qmomkit::QmomSourceContext& state)
    {
        assert(state.size() == 2);
        return qmomkit::QmomMomentRate{std::vector<double>{0.0, 1.0}};
    };
    const qmomkit::QmomMomentBoundaryContribution boundary = [](double, const qmomkit::QmomSourceContext& state)
    {
        assert(state.size() == 2);
        return qmomkit::QmomMomentRate{std::vector<double>{0.25, -0.5}};
    };

    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::QmomPhaseSpaceEvolutionResult phase_space_evolution =
        qmomkit::advance_qmom_phase_space(state, 0.0, 0.25, zero_velocity, source, boundary, options);
    const qmomkit::QmomMomentSourceUpdateResult direct = qmomkit::advance_qmom_moments_with_source(moments, 0.125, 0.25, source, boundary, {});

    assert(phase_space_evolution.ok());
    assert(direct.ok());
    assert_vector_near(phase_space_evolution.state.moment_values(), direct.moments, 1.0e-13);
    assert_qmom_state_canonical(phase_space_evolution.state, 1.0e-13);
}

void test_adaptive_branch_specific_cqmom_case()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = adaptive_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData required = qmomkit::extract_cqmom_moments(moments, layout);
    const qmomkit::TestCqmomPhaseSpaceEvolutionState state{required};
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& values)
    { return qmomkit::CqmomMomentRate{std::vector<double>(values.size(), 0.0)}; };

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options;
    options.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    options.adaptive_branch_counts = true;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::TestCqmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(result.ok());
    assert(result.quadrature.has_value());
    assert(result.quadrature->adaptive_branch_counts_enabled());
    assert(!result.quadrature->has_rectangular_storage());
    assert(result.quadrature->total_nodes() == 3);
    assert_vector_near(result.state.moments.values, required.values, 1.0e-11);
}

void test_adaptive_state_initialization_preserves_policy()
{
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = adaptive_cqmom_prepared_case();
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_prepared_case);
    const qmomkit::CqmomState state = cqmom_prepared_case.initial_state();

    assert(cqmom_prepared_case.adaptive_branch_counts());
    assert(cqmom_prepared_case.qmom_options().adaptive);
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

void test_adaptive_state_source_phase_space_evolution_changes_active_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = adaptive_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = fixed_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const qmomkit::CqmomState state = adaptive_cqmom_prepared_case().initial_state();
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
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

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 1.0, zero_velocity, source, options);

    assert(phase_space_evolved.ok());
    assert(!phase_space_evolved.diagnostic_quadrature.has_value());
    assert(phase_space_evolved.state.quadrature().adaptive_branch_counts_enabled());
    assert(!phase_space_evolved.state.quadrature().has_rectangular_storage());
    assert(phase_space_evolved.state.quadrature().total_nodes() == 4);
    assert(diagnostics_all_use_requested_node_counts(phase_space_evolved.diagnostics.post_source_cqmom_diagnostics));
    assert(phase_space_evolved.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(phase_space_evolved.diagnostics.post_source_cqmom_status == qmomkit::CqmomStatus::success);
    assert_cqmom_state_active_counts_match_accepted_quadrature(phase_space_evolved.state);
}

void test_adaptive_state_preserves_branch_counts_across_steps()
{
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = adaptive_cqmom_prepared_case_from_moments(branch_specific_ragged_cqmom_moments());
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    qmomkit::CqmomPhaseSpaceEvolutionResult first =
        qmomkit::advance_cqmom_phase_space(cqmom_prepared_case.initial_state(), 0.0, 0.1, zero_velocity, zero_source, options);
    assert(first.ok());
    assert(!first.diagnostic_quadrature.has_value());
    assert(first.state.quadrature().total_nodes() == 3);
    assert(active_counts_contain(first.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1));
    assert(active_counts_contain(first.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2));
    assert(first.diagnostics.cqmom_node_count_transitions.empty());
    assert_cqmom_state_moments_match_quadrature(first.state, first.state.quadrature(), 5.0e-10);
    assert_cqmom_state_active_counts_match_accepted_quadrature(first.state);

    const qmomkit::CqmomPhaseSpaceEvolutionResult second =
        qmomkit::advance_cqmom_phase_space(first.state, first.final_time, 0.1, zero_velocity, zero_source, options);
    assert(second.ok());
    assert(!second.diagnostic_quadrature.has_value());
    assert(second.state.quadrature().terminal_paths() == first.state.quadrature().terminal_paths());
    assert(active_counts_contain(second.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1));
    assert(active_counts_contain(second.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2));
    assert_cqmom_state_moments_match_quadrature(second.state, second.state.quadrature(), 5.0e-10);
    assert_cqmom_state_active_counts_match_accepted_quadrature(second.state);
}

void test_adaptive_state_preserved_source_can_reduce_active_counts()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = collapsed_branch_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
        qmomkit::advance_cqmom_phase_space(adaptive_cqmom_prepared_case_from_moments(initial_moments).initial_state(), 0.0, 1.0, zero_velocity,
                                           source_from_to(initial_required, target_required), options);

    assert(phase_space_evolved.ok());
    assert(!phase_space_evolved.diagnostic_quadrature.has_value());
    assert(phase_space_evolved.state.quadrature().total_nodes() == 2);
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 1));
    assert(transitions_contain(phase_space_evolved.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion,
                               1, qmomkit::ConditioningPath{{1}}, 2, 1, qmomkit::CqmomNodeCountTransitionKind::reduced));
    assert_cqmom_state_active_counts_match_accepted_quadrature(phase_space_evolved.state);
}

void test_adaptive_state_source_policy_allows_growth_to_configured_maximum()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = collapsed_branch_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
        qmomkit::advance_cqmom_phase_space(adaptive_cqmom_prepared_case_from_moments(initial_moments).initial_state(), 0.0, 1.0, zero_velocity,
                                           source_from_to(initial_required, target_required), options);

    assert(phase_space_evolved.ok());
    assert(!phase_space_evolved.diagnostic_quadrature.has_value());
    assert(phase_space_evolved.state.quadrature().total_nodes() == 3);
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1));
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2));
    assert(transitions_contain(phase_space_evolved.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion,
                               1, qmomkit::ConditioningPath{{1}}, 1, 2, qmomkit::CqmomNodeCountTransitionKind::grew));
    assert_cqmom_state_active_counts_match_accepted_quadrature(phase_space_evolved.state);
}

void test_adaptive_state_source_growth_forces_adaptive_wheeler_without_mutating_qmom_options()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = collapsed_branch_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = adaptive_cqmom_prepared_case_from_moments_with_fixed_qmom(initial_moments);
    const qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    assert(cqmom_prepared_case.adaptive_branch_counts());
    assert(!cqmom_prepared_case.qmom_options().adaptive);
    assert(state.adaptive_branch_counts());
    assert(!state.qmom_options().adaptive);
    assert(active_counts_contain(state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, 2));
    assert(active_counts_contain(state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1, 1, 2));
    assert(active_counts_contain(state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 1, 1, 2));

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 1.0, zero_velocity, source_from_to(initial_required, target_required), options);

    assert(phase_space_evolved.ok());
    assert(!phase_space_evolved.diagnostic_quadrature.has_value());
    assert(phase_space_evolved.state.quadrature().adaptive_branch_counts_enabled());
    assert(!phase_space_evolved.state.quadrature().has_rectangular_storage());
    assert(phase_space_evolved.state.quadrature().total_nodes() == 3);
    assert(!phase_space_evolved.state.qmom_options().adaptive);
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, 2));
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1, 1, 2));
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2, 1, 2));
    assert(diagnostics_contain_branch_count(phase_space_evolved.diagnostics.post_source_cqmom_diagnostics, 1, qmomkit::ConditioningPath{{1}}, 2, 2));
    assert(transitions_contain(phase_space_evolved.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion,
                               1, qmomkit::ConditioningPath{{1}}, 1, 2, qmomkit::CqmomNodeCountTransitionKind::grew));
    assert_cqmom_state_active_counts_match_accepted_quadrature(phase_space_evolved.state);
    assert(phase_space_evolved.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(phase_space_evolved.diagnostics.post_source_cqmom_status == qmomkit::CqmomStatus::success);

    const qmomkit::QmomResult fixed_qmom = qmomkit::compute_qmom(qmomkit::MomentSet({1.0, 0.0, 0.0, 0.0}), 2, qmomkit::QmomOptions{false});
    const qmomkit::QmomResult adaptive_qmom = qmomkit::compute_qmom(qmomkit::MomentSet({1.0, 0.0, 0.0, 0.0}), 2, qmomkit::QmomOptions{true});
    assert(!fixed_qmom.ok());
    assert(adaptive_qmom.ok());
    assert(adaptive_qmom.weights.size() == 1);
}

void test_adaptive_state_source_growth_canonicalizes_between_steps()
{
    constexpr std::size_t configured_node_count = 3;
    const qmomkit::RawMomentTensor initial_moments = one_dimensional_cqmom_moments(configured_node_count, {0.0}, {1.0});
    const qmomkit::RawMomentTensor target_moments = one_dimensional_cqmom_moments(configured_node_count, {-1.0, 0.25, 2.0}, {0.2, 0.3, 0.5});
    const qmomkit::PhaseSpaceHierarchy hierarchy({configured_node_count});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    qmomkit::CqmomPhaseSpaceEvolutionOptions source_last_options;
    source_last_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;
    source_last_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult first =
        qmomkit::advance_cqmom_phase_space(one_dimensional_adaptive_cqmom_state(initial_required), 0.0, 1.0, zero_velocity,
                                           source_from_to(initial_required, target_required), source_last_options);

    assert(first.ok());
    assert(!first.diagnostic_quadrature.has_value());
    assert(first.state.quadrature().total_nodes() == 2);
    assert(active_counts_contain(first.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, configured_node_count));
    assert(diagnostics_contain_branch_count(first.diagnostics.post_source_cqmom_diagnostics, 0, qmomkit::ConditioningPath{}, 2, 2));
    assert_cqmom_state_moments_match_quadrature(first.state, first.state.quadrature(), 5.0e-10);
    assert_cqmom_state_active_counts_match_accepted_quadrature(first.state);

    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    qmomkit::CqmomPhaseSpaceEvolutionOptions growth_options;
    growth_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    growth_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult second =
        qmomkit::advance_cqmom_phase_space(first.state, first.final_time, 0.1, zero_velocity, zero_source, growth_options);

    assert(second.ok());
    assert(!second.diagnostic_quadrature.has_value());
    assert(second.state.quadrature().total_nodes() == 2);
    assert(active_counts_contain(second.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, configured_node_count));
    assert(second.diagnostics.cqmom_node_count_transitions.empty());
    assert_cqmom_state_moments_match_quadrature(second.state, second.state.quadrature(), 5.0e-10);
    assert_cqmom_state_active_counts_match_accepted_quadrature(second.state);
}

void test_adaptive_state_failed_growth_trial_keeps_current_count()
{
    constexpr std::size_t configured_node_count = 3;
    const qmomkit::RawMomentTensor initial_moments = one_dimensional_cqmom_moments(configured_node_count, {-1.0, 1.0}, {0.35, 0.65});
    const qmomkit::RawMomentTensor target_moments = one_dimensional_cqmom_moments(configured_node_count, {-0.75, 1.25}, {0.4, 0.6});
    const qmomkit::PhaseSpaceHierarchy hierarchy({configured_node_count});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved = qmomkit::advance_cqmom_phase_space(
        one_dimensional_adaptive_cqmom_state(initial_required), 0.0, 1.0, zero_velocity, source_from_to(initial_required, target_required), options);

    assert(phase_space_evolved.ok());
    assert(!phase_space_evolved.diagnostic_quadrature.has_value());
    assert(phase_space_evolved.state.quadrature().total_nodes() == 2);
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, configured_node_count));
    assert(diagnostics_contain_branch_count(phase_space_evolved.diagnostics.post_source_cqmom_diagnostics, 0, qmomkit::ConditioningPath{}, 2, 2));
    assert(transitions_contain(phase_space_evolved.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::initial_inversion, 0,
                               qmomkit::ConditioningPath{}, std::nullopt, 2, qmomkit::CqmomNodeCountTransitionKind::created));
    assert_cqmom_state_active_counts_match_accepted_quadrature(phase_space_evolved.state);
}

void test_adaptive_state_bounded_support_rejects_post_source_growth()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({3, 1});
    const std::vector<std::size_t> extents{6, 2};
    const std::vector<double> initial_sequence{1.0, 0.5, 0.3125, 0.21875, 0.16015625, 0.1220703125};
    std::vector<double> initial_values(qmomkit::product_extents(extents), 0.0);
    std::vector<double> target_values(qmomkit::product_extents(extents), 0.0);
    for (std::size_t order = 0; order < initial_sequence.size(); ++order)
    {
        initial_values[qmomkit::flatten_multi_index({order, 0}, extents)] = initial_sequence[order];
        target_values[qmomkit::flatten_multi_index({order, 0}, extents)] = initial_sequence[order];
    }
    target_values[qmomkit::flatten_multi_index({4, 0}, extents)] = 10.0;
    target_values[qmomkit::flatten_multi_index({5, 0}, extents)] = 10.0;

    const qmomkit::RawMomentTensor initial_moments{extents, initial_values};
    const qmomkit::RawMomentTensor target_moments{extents, target_values};
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    qmomkit::CqmomAbscissaBounds support(2);
    support[0] = qmomkit::AbscissaBounds{0.0, 1.0};

    qmomkit::CqmomCase cqmom_case{initial_required, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, true, qmomkit::CqmomDiagnosticLevel::basic};
    cqmom_case.abscissa_bounds = support;
    const qmomkit::CqmomState state = qmomkit::internal::CqmomConstructionAccess::create_state(cqmom_case);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 1.0, zero_velocity, source_from_to(initial_required, target_required), options);

    assert(phase_space_evolved.ok());
    assert(!phase_space_evolved.diagnostic_quadrature.has_value());
    assert(phase_space_evolved.state.quadrature().total_nodes() == 2);
    assert(phase_space_evolved.state.abscissa_bounds().size() == 2);
    assert(phase_space_evolved.state.abscissa_bounds()[0].has_value());
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, 3));
    assert(!phase_space_evolved.diagnostics.post_source_cqmom_diagnostics.empty());
    const qmomkit::CqmomDiagnostics& diagnostic = phase_space_evolved.diagnostics.post_source_cqmom_diagnostics.front();
    assert(diagnostic.bounded_support_checked);
    assert(diagnostic.bounded_candidate_attempt_type == "growth");
    assert(diagnostic.bounded_candidate_node_count == 3);
    assert(diagnostic.bounded_accepted_fallback_node_count == 2);
    assert_cqmom_state_active_counts_match_accepted_quadrature(phase_space_evolved.state);
}

void test_adaptive_state_reduction_prevents_source_growth_trial()
{
    constexpr std::size_t configured_node_count = 3;
    const qmomkit::RawMomentTensor initial_moments = one_dimensional_cqmom_moments(configured_node_count, {-1.0, 1.0}, {0.35, 0.65});
    const qmomkit::RawMomentTensor target_moments = one_dimensional_cqmom_moments(configured_node_count, {0.25}, {1.0});
    const qmomkit::PhaseSpaceHierarchy hierarchy({configured_node_count});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved = qmomkit::advance_cqmom_phase_space(
        one_dimensional_adaptive_cqmom_state(initial_required), 0.0, 1.0, zero_velocity, source_from_to(initial_required, target_required), options);

    assert(phase_space_evolved.ok());
    assert(!phase_space_evolved.diagnostic_quadrature.has_value());
    assert(phase_space_evolved.state.quadrature().total_nodes() == 1);
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 1, 0, configured_node_count));
    assert(diagnostics_contain_branch_count(phase_space_evolved.diagnostics.post_source_cqmom_diagnostics, 0, qmomkit::ConditioningPath{}, 2, 1));
    assert_cqmom_state_active_counts_match_accepted_quadrature(phase_space_evolved.state);
}

void test_adaptive_state_parent_growth_creates_unconstrained_descendants()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = root_collapsed_child_spread_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = fixed_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
        qmomkit::advance_cqmom_phase_space(adaptive_cqmom_prepared_case_from_moments(initial_moments).initial_state(), 0.0, 1.0, zero_velocity,
                                           source_from_to(initial_required, target_required), options);

    assert(phase_space_evolved.ok());
    assert(!phase_space_evolved.diagnostic_quadrature.has_value());
    assert(phase_space_evolved.state.quadrature().total_nodes() == 4);
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2));
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 2));
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2));
    assert(transitions_contain(phase_space_evolved.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion,
                               0, qmomkit::ConditioningPath{}, 1, 2, qmomkit::CqmomNodeCountTransitionKind::grew));
    assert(transitions_contain(phase_space_evolved.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion,
                               1, qmomkit::ConditioningPath{{1}}, std::nullopt, 2, qmomkit::CqmomNodeCountTransitionKind::created));
    assert_cqmom_state_active_counts_match_accepted_quadrature(phase_space_evolved.state);
}

void test_adaptive_state_preserved_source_invalidates_descendants_after_root_collapse()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = root_collapsed_child_spread_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
        qmomkit::advance_cqmom_phase_space(adaptive_cqmom_prepared_case_from_moments(initial_moments).initial_state(), 0.0, 1.0, zero_velocity,
                                           source_from_to(initial_required, target_required), options);

    assert(phase_space_evolved.ok());
    assert(!phase_space_evolved.diagnostic_quadrature.has_value());
    assert(phase_space_evolved.state.quadrature().total_nodes() == 2);
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 1));
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 2));
    assert(transitions_contain(phase_space_evolved.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion,
                               0, qmomkit::ConditioningPath{}, 2, 1, qmomkit::CqmomNodeCountTransitionKind::reduced));
    assert(transitions_contain(phase_space_evolved.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion,
                               1, qmomkit::ConditioningPath{{0}}, 1, std::nullopt, qmomkit::CqmomNodeCountTransitionKind::removed));
    assert(transitions_contain(phase_space_evolved.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion,
                               1, qmomkit::ConditioningPath{{0}}, std::nullopt, 2, qmomkit::CqmomNodeCountTransitionKind::created));
    assert_cqmom_state_active_counts_match_accepted_quadrature(phase_space_evolved.state);
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
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved = qmomkit::advance_cqmom_phase_space(
        adaptive_cqmom_prepared_case_from_moments_with_hierarchy(initial_moments, node_counts, conditioning_order).initial_state(), 0.0, 1.0,
        zero_velocity, source_from_to(initial_required, target_required), options);

    assert(phase_space_evolved.ok());
    assert(!phase_space_evolved.diagnostic_quadrature.has_value());
    assert(phase_space_evolved.state.quadrature().total_nodes() == 3);
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2));
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1));
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 2, qmomkit::ConditioningPath{{0, 0}}, 2));
    assert(active_counts_contain(phase_space_evolved.state.active_branch_counts(), 2, qmomkit::ConditioningPath{{1, 0}}, 1));
    assert(transitions_contain(phase_space_evolved.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion,
                               1, qmomkit::ConditioningPath{{0}}, 2, 1, qmomkit::CqmomNodeCountTransitionKind::reduced));
    assert(transitions_contain(phase_space_evolved.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion,
                               2, qmomkit::ConditioningPath{{0, 0}}, 1, std::nullopt, qmomkit::CqmomNodeCountTransitionKind::removed));
    assert(transitions_contain(phase_space_evolved.diagnostics.cqmom_node_count_transitions, qmomkit::PhaseSpaceEvolutionStage::post_source_inversion,
                               2, qmomkit::ConditioningPath{{0, 0}}, std::nullopt, 2, qmomkit::CqmomNodeCountTransitionKind::created));
    assert_cqmom_state_active_counts_match_accepted_quadrature(phase_space_evolved.state);
}

void test_adaptive_state_initial_counts_preserve_nonidentity_conditioning_order()
{
    const std::vector<std::size_t> node_counts{2, 2};
    const std::vector<std::size_t> conditioning_order{1, 0};
    const qmomkit::RawMomentTensor moments = nonidentity_ragged_cqmom_moments();
    const qmomkit::CqmomState state =
        adaptive_cqmom_prepared_case_from_moments_with_hierarchy(moments, node_counts, conditioning_order).initial_state();

    assert(state.adaptive_branch_counts());
    assert(active_counts_contain(state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 1, 2));
    assert(active_counts_contain(state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 2, 0, 2));
    assert(active_counts_contain(state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 1, 0, 2));
}

void test_fixed_state_phase_space_evolution_does_not_store_adaptive_branch_counts()
{
    const qmomkit::CqmomState state = fixed_cqmom_prepared_case().initial_state();
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(phase_space_evolved.ok());
    assert(!phase_space_evolved.diagnostic_quadrature.has_value());
    assert(phase_space_evolved.state.quadrature().has_rectangular_storage());
    assert(phase_space_evolved.state.quadrature().total_nodes() == 4);
    assert(phase_space_evolved.state.active_branch_counts().empty());
    assert(phase_space_evolved.diagnostics.cqmom_node_count_transitions.empty());
}

void test_adaptive_state_rejects_duplicate_active_count_metadata()
{
    const qmomkit::CqmomState good_state = adaptive_cqmom_prepared_case().initial_state();
    const qmomkit::CqmomCase& cqmom_case = qmomkit::internal::CqmomConstructionAccess::state_definition(good_state);
    std::vector<qmomkit::CqmomActiveBranchCount> bad_counts = good_state.active_branch_counts();
    assert(!bad_counts.empty());
    bad_counts.push_back(bad_counts.front());
    const qmomkit::CqmomState bad_state = qmomkit::internal::CqmomConstructionAccess::create_state(cqmom_case, std::move(bad_counts));
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_cqmom_phase_space(bad_state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::initial_inversion);
    assert(result.diagnostics.initial_cqmom_status == qmomkit::CqmomStatus::invalid_input);
}

void test_adaptive_state_rejects_missing_active_count_metadata()
{
    const qmomkit::CqmomState good_state = adaptive_cqmom_prepared_case().initial_state();
    const qmomkit::CqmomCase& cqmom_case = qmomkit::internal::CqmomConstructionAccess::state_definition(good_state);
    std::vector<qmomkit::CqmomActiveBranchCount> bad_counts = good_state.active_branch_counts();
    assert(!bad_counts.empty());
    bad_counts.pop_back();
    const qmomkit::CqmomState bad_state = qmomkit::internal::CqmomConstructionAccess::create_state(cqmom_case, std::move(bad_counts));
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_cqmom_phase_space(bad_state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::initial_inversion);
    assert(result.diagnostics.initial_cqmom_status == qmomkit::CqmomStatus::invalid_input);
}

void test_adaptive_state_failure_preserves_active_counts_for_retry()
{
    const qmomkit::CqmomState state = adaptive_cqmom_prepared_case().initial_state();
    const std::vector<qmomkit::CqmomActiveBranchCount> original_counts = state.active_branch_counts();
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto bad_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult failed = qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, bad_source, options);
    assert(!failed.ok());
    assert(active_counts_equal(state.active_branch_counts(), original_counts));
    assert(active_counts_equal(failed.state.active_branch_counts(), original_counts));

    const qmomkit::CqmomPhaseSpaceEvolutionResult retry = qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, zero_source, options);
    assert(retry.ok());
    assert(active_counts_equal(retry.state.active_branch_counts(), original_counts));
}

void test_adaptive_state_post_source_failure_preserves_active_counts_for_retry()
{
    const qmomkit::CqmomState state = adaptive_cqmom_prepared_case().initial_state();
    const std::vector<qmomkit::CqmomActiveBranchCount> original_counts = state.active_branch_counts();
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto bad_realizability_source = [](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        rate.front() = -20.0;
        return qmomkit::CqmomMomentRate{rate};
    };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult failed =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, bad_realizability_source, options);
    assert(!failed.ok());
    assert(failed.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    assert(active_counts_equal(state.active_branch_counts(), original_counts));
    assert(active_counts_equal(failed.state.active_branch_counts(), original_counts));

    const qmomkit::CqmomPhaseSpaceEvolutionResult retry = qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, zero_source, options);
    assert(retry.ok());
    assert(active_counts_equal(retry.state.active_branch_counts(), original_counts));
}

void test_adaptive_state_phase_space_evolution_reports_initial_failure()
{
    qmomkit::RawMomentTensor moments = adaptive_cqmom_moments();
    std::fill(moments.values.begin(), moments.values.end(), 0.0);
    const std::string text = cqmom_case_text_from_moments(moments, true, true);
    const auto parsed = qmomkit::config::parse_cqmom_prepared(text, "adaptive-failing-cqmom-case");
    assert(parsed.ok());
    bool threw = false;
    try
    {
        static_cast<void>(parsed.value().initial_state());
    }
    catch (const std::runtime_error& error)
    {
        threw = true;
        assert(std::string(error.what()).find("quadrature reconstruction failed") != std::string::npos);
    }
    assert(threw);
}

void test_state_cqmom_phase_space_evolution_matches_cqmom_case_fixed_layout()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestCqmomPhaseSpaceEvolutionState test_state = cqmom_state_from_tensor(moments, hierarchy);
    const qmomkit::CqmomState state = fixed_cqmom_prepared_case().initial_state();
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& xi)
        {
            assert(xi.size() == 2);
            return std::vector<double>{0.2 * xi[0] + 0.1, -0.15 * xi[1] + 0.05};
        });
    constexpr double lambda = -0.35;
    const auto source = [](double, const qmomkit::CqmomSourceContext& state_values)
    {
        std::vector<double> rate(state_values.size(), 0.0);
        for (std::size_t i = 0; i < state_values.size(); ++i)
        {
            rate[i] = lambda * state_values.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };

    qmomkit::CqmomPhaseSpaceEvolutionOptions state_options;
    state_options.node_options = tight_node_options();
    state_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    const qmomkit::CqmomPhaseSpaceEvolutionResult state_result = qmomkit::advance_cqmom_phase_space(state, 0.0, 0.2, velocity, source, state_options);

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions test_options;
    test_options.node_options = state_options.node_options;
    test_options.scheme = state_options.scheme;
    const qmomkit::TestCqmomPhaseSpaceEvolutionResult test_result =
        qmomkit::advance_cqmom_phase_space(test_state, 0.0, 0.2, velocity, source, test_options);

    assert(state_result.ok());
    assert(test_result.ok());
    assert(!state_result.diagnostic_quadrature.has_value());
    assert(test_result.quadrature.has_value());
    assert_cqmom_quadrature_near(state_result.state.quadrature(), *test_result.quadrature, 5.0e-9);
    assert(state_result.diagnostics.scheme == test_result.diagnostics.scheme);
    assert(state_result.final_time == test_result.final_time);
    assert(state_result.state.hierarchy().node_counts() == state.hierarchy().node_counts());
    assert(state_result.state.hierarchy().conditioning_order() == state.hierarchy().conditioning_order());
}

void test_cqmom_source_context_reconstructs_untracked_moments_from_each_rk_stage()
{
    const qmomkit::MultiIndex tracked_index{0, 0};
    const qmomkit::MultiIndex first_untracked_index{4, 0};
    const qmomkit::MultiIndex second_untracked_index{0, 4};
    constexpr double lambda = 0.1;
    constexpr double time_step = 0.2;

    {
        const qmomkit::CqmomPreparedCase cqmom_prepared_case = fixed_cqmom_prepared_case();
        const qmomkit::CqmomResult initial_quadrature = qmomkit::compute_cqmom(cqmom_prepared_case);
        assert(initial_quadrature.ok());
        const qmomkit::SelectedRawMomentResult initial_moments = qmomkit::reconstruct_raw_moments(
            initial_quadrature.quadrature, qmomkit::RawMomentRequest{{tracked_index, first_untracked_index, second_untracked_index}});
        assert(initial_moments.ok());
        const double first_untracked_per_number_density = initial_moments.values[1] / initial_moments.values[0];
        const double second_untracked_per_number_density = initial_moments.values[2] / initial_moments.values[0];

        std::size_t calls = 0;
        const qmomkit::CqmomMomentSource source = [&](double, const qmomkit::CqmomSourceContext& context)
        {
            assert(context.position(first_untracked_index) == std::nullopt);
            assert(context.position(second_untracked_index) == std::nullopt);
            assert_near(context.moment(tracked_index), context.value(tracked_index).value(), 0.0);

            const double number_density = context.moment(tracked_index);
            assert_near(context.moment(first_untracked_index), first_untracked_per_number_density * number_density, 2.0e-10);
            assert_near(context.moment(second_untracked_index), second_untracked_per_number_density * number_density, 2.0e-10);
            ++calls;

            std::vector<double> rate(context.size(), 0.0);
            for (std::size_t entry = 0; entry < context.size(); ++entry)
            {
                rate[entry] = lambda * context.value(entry);
            }
            return qmomkit::CqmomMomentRate{rate};
        };

        const qmomkit::CqmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_cqmom_moment_source(cqmom_prepared_case.initial_state(), 0.0, time_step, source);
        assert(result.ok());
        assert(calls == 3);
    }
}

void test_cqmom_source_context_visits_terminal_nodes_from_each_rk_stage()
{
    {
        const qmomkit::CqmomPreparedCase cqmom_prepared_case = fixed_cqmom_prepared_case();
        const qmomkit::CqmomResult expected = qmomkit::compute_cqmom(cqmom_prepared_case);
        assert(expected.ok());

        std::size_t calls = 0;
        const qmomkit::CqmomMomentSource source = [&](double, const qmomkit::CqmomSourceContext& context)
        {
            std::vector<double> weights;
            std::vector<std::vector<double>> coordinates;
            context.for_each_terminal_node(
                [&weights, &coordinates](const double weight, const std::vector<double>& node_coordinates)
                {
                    weights.push_back(weight);
                    coordinates.push_back(node_coordinates);
                });

            assert_vector_near(weights, expected.quadrature.weights(), 2.0e-10);
            assert(coordinates.size() == expected.quadrature.total_nodes());
            for (std::size_t node = 0; node < coordinates.size(); ++node)
            {
                assert(coordinates[node].size() == expected.quadrature.number_of_dimensions());
                for (std::size_t dimension = 0; dimension < coordinates[node].size(); ++dimension)
                {
                    assert_near(coordinates[node][dimension], expected.quadrature.abscissas()[dimension][node], 2.0e-10);
                }
            }

            ++calls;
            return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)};
        };

        const qmomkit::CqmomPhaseSpaceEvolutionResult result =
            qmomkit::advance_cqmom_moment_source(cqmom_prepared_case.initial_state(), 0.0, 0.2, source);
        assert(result.ok());
        assert(calls == 3);
    }
}

qmomkit::CqmomPreparedCase fixed_cqmom_prepared_case_with_hierarchy(const qmomkit::RawMomentTensor& moments,
                                                                    const std::vector<std::size_t>& conditioning_order)
{
    const std::string text = cqmom_case_text_from_moments_with_hierarchy(moments, {2, 2}, conditioning_order, false, false);
    const auto parsed = qmomkit::config::parse_cqmom_prepared(text, "fixed-cqmom-case-with-hierarchy");
    assert(parsed.ok());
    return parsed.value();
}

qmomkit::CqmomPreparedCase generalized_moment_cqmom_prepared_case()
{
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.1, 0.4};
    // The first coordinate is strictly positive for fractional powers; the
    // second deliberately has negative nodes to exercise integral powers.
    const std::vector<std::vector<double>> abscissas{{1.0, 1.0, 2.0, 2.0}, {-1.0, 1.0, -1.0, 1.0}};
    return fixed_cqmom_prepared_case_with_hierarchy({extents, raw_moments_from_terminals(extents, weights, abscissas)}, {0, 1});
}

void test_cqmom_source_context_generalized_moments_use_stage_quadrature()
{
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = generalized_moment_cqmom_prepared_case();
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_prepared_case);
    assert(initial.ok());
    const std::vector<double> fractional_orders{1.0 / 3.0, 0.0};
    const std::vector<double> mixed_orders{4.0 / 3.0, 2.0};
    const std::vector<double> integral_orders{2.0, 1.0};
    const double expected_fractional = qmomkit::evaluate_generalized_moment(initial.quadrature, fractional_orders);
    const double expected_mixed = qmomkit::evaluate_generalized_moment(initial.quadrature, mixed_orders);

    std::size_t source_calls = 0;
    const qmomkit::CqmomMomentSource source = [&](double, const qmomkit::CqmomSourceContext& context)
    {
        assert_near(context.generalized_moment(fractional_orders), expected_fractional, 2.0e-10);
        assert_near(context.generalized_moment(mixed_orders), expected_mixed, 2.0e-10);
        assert_near(context.generalized_moment(integral_orders), context.moment({2, 1}), 2.0e-10);
        // Multiple requests in the same source stage must reuse the context's
        // lazy source-stage quadrature rather than needing a second path.
        assert_near(context.generalized_moment(fractional_orders), expected_fractional, 2.0e-10);
        ++source_calls;
        return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)};
    };

    const qmomkit::CqmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_cqmom_moment_source(cqmom_prepared_case.initial_state(), 0.0, 0.1, source);
    assert(result.ok());
    assert(source_calls == 3);
}

void test_cqmom_source_context_generalized_moment_validation()
{
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = generalized_moment_cqmom_prepared_case();
    bool checked = false;
    const qmomkit::CqmomMomentSource source = [&](double, const qmomkit::CqmomSourceContext& context)
    {
        const auto expect_invalid_argument = [&context](const std::vector<double>& orders)
        {
            bool rejected = false;
            try
            {
                static_cast<void>(context.generalized_moment(orders));
            }
            catch (const std::invalid_argument&)
            {
                rejected = true;
            }
            assert(rejected);
        };
        expect_invalid_argument({1.0});
        expect_invalid_argument({-0.5, 0.0});
        expect_invalid_argument({std::numeric_limits<double>::quiet_NaN(), 0.0});
        expect_invalid_argument({std::numeric_limits<double>::infinity(), 0.0});
        expect_invalid_argument({-std::numeric_limits<double>::infinity(), 0.0});

        bool rejected_fractional_negative_abscissa = false;
        try
        {
            static_cast<void>(context.generalized_moment({1.0, 0.5}));
        }
        catch (const std::domain_error&)
        {
            rejected_fractional_negative_abscissa = true;
        }
        assert(rejected_fractional_negative_abscissa);
        // The same negative coordinate is valid for a real-valued integral order.
        assert_near(context.generalized_moment({1.0, 3.0}), context.moment({1, 3}), 2.0e-10);
        checked = true;
        return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)};
    };

    const qmomkit::CqmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_cqmom_moment_source(cqmom_prepared_case.initial_state(), 0.0, 0.1, source);
    assert(result.ok());
    assert(checked);
}

void test_state_cqmom_phase_space_evolution_preserves_vandermonde_solver_policy()
{
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = fixed_cqmom_prepared_case_with_vandermonde_solver("dense");
    const qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
    assert(state.vandermonde_solver() == qmomkit::VandermondeSolver::dense);

    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(cqmom_prepared_case);
    assert(initial.ok());

    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(phase_space_evolved.ok());
    assert(!phase_space_evolved.diagnostic_quadrature.has_value());
    assert(phase_space_evolved.state.vandermonde_solver() == qmomkit::VandermondeSolver::dense);
    assert(phase_space_evolved.state.vandermonde_solver() == state.vandermonde_solver());
    assert_cqmom_quadrature_near(phase_space_evolved.state.quadrature(), initial.quadrature, 5.0e-10);
}

qmomkit::RawMomentTensor phase_space_evolution_vandermonde_policy_moments()
{
    return fixed_cqmom_moments();
}

void test_cqmom_phase_space_evolution_vandermonde_policy_independent_of_diagnostic_level()
{
    const qmomkit::RawMomentTensor moments = phase_space_evolution_vandermonde_policy_moments();
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    std::optional<qmomkit::CqmomPhaseSpaceEvolutionResult> reference;
    for (const qmomkit::CqmomDiagnosticLevel level :
         {qmomkit::CqmomDiagnosticLevel::none, qmomkit::CqmomDiagnosticLevel::basic, qmomkit::CqmomDiagnosticLevel::detailed})
    {
        const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
        const qmomkit::CqmomCase cqmom_case{qmomkit::extract_cqmom_moments(moments, layout), {}, {}, false, level,
                                            qmomkit::VandermondeSolver::bjorck_pereyra};

        const qmomkit::CqmomPhaseSpaceEvolutionResult evolved = qmomkit::advance_cqmom_phase_space(
            qmomkit::internal::CqmomConstructionAccess::create_state(cqmom_case), 0.0, 0.1, zero_velocity, zero_source, options);

        if (level == qmomkit::CqmomDiagnosticLevel::none)
        {
            assert(evolved.diagnostics.initial_cqmom_diagnostics.empty());
            assert(evolved.diagnostics.post_source_cqmom_diagnostics.empty());
        }
        if (!reference.has_value())
        {
            reference = evolved;
            continue;
        }

        assert(evolved.status == reference->status);
        assert(evolved.ok() == reference->ok());
        assert(evolved.diagnostics.failed_stage == reference->diagnostics.failed_stage);
        assert(evolved.diagnostics.initial_cqmom_status == reference->diagnostics.initial_cqmom_status);
        assert(evolved.diagnostics.post_source_cqmom_status == reference->diagnostics.post_source_cqmom_status);
        assert(evolved.diagnostic_quadrature.has_value() == reference->diagnostic_quadrature.has_value());
        if (evolved.diagnostic_quadrature.has_value())
        {
            assert_cqmom_quadrature_near(*evolved.diagnostic_quadrature, *reference->diagnostic_quadrature, 1.0e-5);
        }
        if (evolved.ok())
        {
            assert_cqmom_quadrature_near(evolved.state.quadrature(), reference->state.quadrature(), 1.0e-5);
            assert_cqmom_moments_near(evolved.state, reference->state, 1.0e-5);
        }
        assert(active_counts_equal(evolved.state.active_branch_counts(), reference->state.active_branch_counts()));
    }
}

void test_state_cqmom_source_evolution_matches_cqmom_case()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::TestCqmomPhaseSpaceEvolutionState test_state = cqmom_state_from_tensor(fixed_cqmom_moments(), hierarchy);
    const qmomkit::CqmomState state = fixed_cqmom_prepared_case().initial_state();
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto source = [](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            rate[i] = 0.1 * moments.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };

    qmomkit::CqmomPhaseSpaceEvolutionOptions state_options;
    state_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    state_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult state_result =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, source, state_options);

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions test_options;
    test_options.scheme = state_options.scheme;
    test_options.node_options = state_options.node_options;
    const qmomkit::TestCqmomPhaseSpaceEvolutionResult test_result =
        qmomkit::advance_cqmom_phase_space(test_state, 0.0, 0.1, zero_velocity, source, test_options);

    assert(state_result.ok());
    assert(test_result.ok());
    assert(!state_result.diagnostic_quadrature.has_value());
    assert(test_result.quadrature.has_value());
    assert_cqmom_quadrature_near(state_result.state.quadrature(), *test_result.quadrature, 5.0e-10);
    assert(state_result.diagnostics.moment_source_status == test_result.diagnostics.moment_source_status);
}

void test_state_cqmom_failure_propagation_from_source_operator()
{
    const qmomkit::CqmomState state = fixed_cqmom_prepared_case().initial_state();
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto bad_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; };

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult result = qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, bad_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    assert(result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.message == "CQMOM moment-source rate values must match the source context size.");
    assert(result.diagnostics.message == result.message);
}

void test_state_cqmom_phase_space_evolution_api_does_not_expose_moments()
{
    static_assert(!HasPublicMomentsMember<qmomkit::CqmomPhaseSpaceEvolutionOptions>::value);
    static_assert(!HasPublicMomentsMember<qmomkit::CqmomPhaseSpaceEvolutionResult>::value);
    static_assert(std::is_same_v<decltype(std::declval<qmomkit::CqmomPhaseSpaceEvolutionResult&>().state), qmomkit::CqmomState>);
    static_assert(std::is_copy_constructible_v<qmomkit::CqmomState>);
    static_assert(std::is_move_constructible_v<qmomkit::CqmomState>);
    static_assert(std::is_copy_assignable_v<qmomkit::CqmomState>);
    static_assert(std::is_move_assignable_v<qmomkit::CqmomState>);
}

void test_state_qmom_phase_space_evolution_api_does_not_expose_moments()
{
    static_assert(!HasPublicMomentsMember<qmomkit::QmomState>::value);
    static_assert(!HasPublicQuadratureMember<qmomkit::QmomState>::value);
    static_assert(!HasPublicMomentsMember<qmomkit::QmomPhaseSpaceEvolutionResult>::value);
    static_assert(std::is_same_v<decltype(std::declval<qmomkit::QmomPhaseSpaceEvolutionResult&>().state), qmomkit::QmomState>);
    static_assert(std::is_copy_constructible_v<qmomkit::QmomState>);
    static_assert(std::is_move_constructible_v<qmomkit::QmomState>);
    static_assert(std::is_copy_assignable_v<qmomkit::QmomState>);
    static_assert(std::is_move_assignable_v<qmomkit::QmomState>);
    static_assert(!std::is_default_constructible_v<qmomkit::QmomState>);
    static_assert(!std::is_aggregate_v<qmomkit::QmomState>);
    static_assert(!std::is_default_constructible_v<qmomkit::QmomPhaseSpaceEvolutionResult>);
    static_assert(std::is_copy_constructible_v<qmomkit::QmomPhaseSpaceEvolutionResult>);
    static_assert(std::is_move_constructible_v<qmomkit::QmomPhaseSpaceEvolutionResult>);
}

void test_qmom_state_copy_move_workflow()
{
    qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.0}, {0.4, 0.6}, 4), 2);
    qmomkit::QmomState copied = state;
    qmomkit::QmomState moved = std::move(copied);
    assert_qmom_state_canonical(moved, 1.0e-12);

    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto zero_source = [](double, const qmomkit::QmomSourceContext& moments)
    { return qmomkit::QmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(moved, 0.0, 0.1, zero_velocity, zero_source, options);
    assert(result.ok());
    state = result.state;
    assert_qmom_state_canonical(state, 1.0e-12);
    state = std::move(result.state);
    assert_qmom_state_canonical(state, 1.0e-12);
}

void test_qmom_initial_state_canonicalizes_adaptive_moments()
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = true;
    qmom_options.rmin = 1.0e-10;
    qmom_options.eabs = 1.0e-12;
    const std::vector<double> moments{1.0, 0.8, 2.8, 4.4, 0.0, 0.0};

    const qmomkit::QmomState state = qmom_state_from_moments(moments, 3, qmom_options);

    assert(state.maximum_node_count() == 3);
    assert(state.qmom_options().adaptive);
    assert_near(state.qmom_options().rmin, qmom_options.rmin, 0.0);
    assert_near(state.qmom_options().eabs, qmom_options.eabs, 0.0);
    assert(state.active_node_count() == 2);
    assert_qmom_state_canonical(state, 1.0e-12);
}

void test_qmom_adaptive_phase_space_evolution_preserves_active_count_with_zero_source()
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = true;
    qmom_options.rmin = 1.0e-10;
    qmom_options.eabs = 1.0e-12;
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.75, 1.25}, {0.4, 0.6}, 6), 3, qmom_options);
    assert(state.active_node_count() == 2);

    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>{0.2 * xi[0]}; });
    const auto zero_source = [](double, const qmomkit::QmomSourceContext& moments)
    { return qmomkit::QmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;
    options.node_options = tight_node_options();

    const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(state, 0.0, 0.2, velocity, zero_source, options);

    assert(result.ok());
    assert(result.state.maximum_node_count() == 3);
    assert(result.state.qmom_options().adaptive);
    assert_near(result.state.qmom_options().rmin, qmom_options.rmin, 0.0);
    assert_near(result.state.qmom_options().eabs, qmom_options.eabs, 0.0);
    assert(result.state.active_node_count() == 2);
    assert_qmom_state_canonical(result.state, 5.0e-10);
}

void test_qmom_adaptive_post_source_reduction_and_growth_policy()
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = true;
    qmom_options.rmin = 1.0e-10;
    qmom_options.eabs = 1.0e-12;
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    const auto source_from_to = [](const std::vector<double>& initial, const std::vector<double>& target)
    {
        return [initial, target](double, const qmomkit::QmomSourceContext& state)
        {
            assert(state.size() == initial.size());
            std::vector<double> rate(state.size(), 0.0);
            for (std::size_t i = 0; i < rate.size(); ++i)
            {
                rate[i] = target[i] - initial[i];
            }
            return qmomkit::QmomMomentRate{std::move(rate)};
        };
    };

    const qmomkit::QmomState two_node_state = qmom_state_from_moments(raw_moments_1d({-0.75, 1.25}, {0.4, 0.6}, 6), 3, qmom_options);
    const std::vector<double> one_node_moments = raw_moments_1d({0.5}, {1.0}, 6);
    const qmomkit::QmomPhaseSpaceEvolutionResult reduction = qmomkit::advance_qmom_phase_space(
        two_node_state, 0.0, 1.0, zero_velocity, source_from_to(two_node_state.moment_values(), one_node_moments), options);
    assert(reduction.ok());
    assert(reduction.state.active_node_count() == 1);
    assert_qmom_state_canonical(reduction.state, 1.0e-12);

    const qmomkit::QmomState one_node_state = qmom_state_from_moments(one_node_moments, 3, qmom_options);
    const std::vector<double> two_node_moments = raw_moments_1d({-0.25, 1.5}, {0.35, 0.65}, 6);
    const qmomkit::QmomPhaseSpaceEvolutionResult growth = qmomkit::advance_qmom_phase_space(
        one_node_state, 0.0, 1.0, zero_velocity, source_from_to(one_node_state.moment_values(), two_node_moments), options);
    assert(growth.ok());
    assert(growth.state.active_node_count() == 2);
    assert_qmom_state_canonical(growth.state, 1.0e-12);

    const std::vector<double> three_node_moments = raw_moments_1d({-1.0, 0.25, 1.5}, {0.2, 0.3, 0.5}, 6);
    const qmomkit::QmomPhaseSpaceEvolutionResult capped_growth = qmomkit::advance_qmom_phase_space(
        one_node_state, 0.0, 1.0, zero_velocity, source_from_to(one_node_state.moment_values(), three_node_moments), options);
    assert(capped_growth.ok());
    assert(capped_growth.state.active_node_count() == 2);
    assert_qmom_state_canonical(capped_growth.state, 1.0e-12);

    const std::vector<double> failed_growth_moments{1.0, 0.8, 2.8, 4.4, 0.0, 0.0};
    const qmomkit::QmomPhaseSpaceEvolutionResult failed_growth = qmomkit::advance_qmom_phase_space(
        two_node_state, 0.0, 1.0, zero_velocity, source_from_to(two_node_state.moment_values(), failed_growth_moments), options);
    assert(failed_growth.ok());
    assert(failed_growth.state.active_node_count() == 2);
    assert_qmom_state_canonical(failed_growth.state, 1.0e-12);
}

void test_qmom_fixed_post_source_requires_configured_count()
{
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.0}, {0.4, 0.6}, 4), 2, qmomkit::QmomOptions{false});
    const std::vector<double> one_node_moments = raw_moments_1d({0.0}, {1.0}, 4);
    const auto source = [initial = state.moment_values(), one_node_moments](double, const qmomkit::QmomSourceContext& moments)
    {
        assert(moments.size() == initial.size());
        std::vector<double> rate(moments.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            rate[i] = one_node_moments[i] - initial[i];
        }
        return qmomkit::QmomMomentRate{std::move(rate)};
    };
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(state, 0.0, 1.0, zero_velocity, source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion);
    assert(result.diagnostics.post_source_qmom_status != qmomkit::QmomStatus::success);
    assert_vector_near(result.state.moment_values(), state.moment_values(), 0.0);
    assert_vector_near(result.state.quadrature().weights, state.quadrature().weights, 0.0);
    assert_vector_near(result.state.quadrature().abscissas, state.quadrature().abscissas, 0.0);
}

void test_first_order_schemes_differ_for_noncommuting_qmom()
{
    const std::vector<double> moments = raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4);
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 2);
    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{0.35}; });
    const auto source = [](double, const qmomkit::QmomSourceContext& state)
    {
        std::vector<double> rhs(state.size(), 0.0);
        for (std::size_t order = 1; order < state.size(); ++order)
        {
            rhs[order] = -0.2 * static_cast<double>(order) * state.value(order);
        }
        return qmomkit::QmomMomentRate{std::move(rhs)};
    };

    qmomkit::QmomPhaseSpaceEvolutionOptions advsrc_options;
    advsrc_options.node_options = tight_node_options();
    advsrc_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;
    const qmomkit::QmomPhaseSpaceEvolutionResult advsrc = qmomkit::advance_qmom_phase_space(state, 0.0, 0.2, velocity, source, advsrc_options);

    qmomkit::QmomPhaseSpaceEvolutionOptions srcadv_options = advsrc_options;
    srcadv_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    const qmomkit::QmomPhaseSpaceEvolutionResult srcadv = qmomkit::advance_qmom_phase_space(state, 0.0, 0.2, velocity, source, srcadv_options);

    assert(advsrc.ok());
    assert(srcadv.ok());
    assert(advsrc.diagnostics.scheme == qmomkit::PhaseSpaceEvolutionScheme::advsrc);
    assert(srcadv.diagnostics.scheme == qmomkit::PhaseSpaceEvolutionScheme::srcadv);
    assert(std::abs(advsrc.state.moment_values()[1] - srcadv.state.moment_values()[1]) > 1.0e-4);
}

void test_strang_source_uses_physical_time_interval_for_nonautonomous_source()
{
    // With zero advection, dm0/dt = time has exact physical update
    // m0(t + h) = m0(t) + t h + h^2 / 2.  SSP-RK3 is exact for this
    // linear-in-time rate, so this checks the clock, not RK truncation error.
    const qmomkit::QmomState state = qmom_state_from_moments({1.0, 0.0}, 1);
    std::vector<double> callback_times;
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto source = [&callback_times](const double time, const qmomkit::QmomSourceContext&)
    {
        callback_times.push_back(time);
        return qmomkit::QmomMomentRate{std::vector<double>{time, 0.0}};
    };

    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::QmomPhaseSpaceEvolutionResult source_only = qmomkit::advance_qmom_moment_source(state, 2.0, 1.0, source, options);
    assert(source_only.ok());
    assert_near(source_only.state.moment_values()[0], 3.5, 1.0e-12);
    qmomkit::QmomPhaseSpaceEvolutionOptions advsrc_options = options;
    advsrc_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;
    const qmomkit::QmomPhaseSpaceEvolutionResult advsrc = qmomkit::advance_qmom_phase_space(state, 2.0, 1.0, zero_velocity, source, advsrc_options);
    assert(advsrc.ok());
    assert_near(advsrc.state.moment_values()[0], 3.5, 1.0e-12);
    qmomkit::QmomPhaseSpaceEvolutionOptions srcadv_options = options;
    srcadv_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    const qmomkit::QmomPhaseSpaceEvolutionResult srcadv = qmomkit::advance_qmom_phase_space(state, 2.0, 1.0, zero_velocity, source, srcadv_options);
    assert(srcadv.ok());
    assert_near(srcadv.state.moment_values()[0], 3.5, 1.0e-12);
    callback_times.clear();
    const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(state, 2.0, 1.0, zero_velocity, source, options);
    assert(result.ok());
    assert_vector_near(callback_times, {2.0, 3.0, 2.5}, 0.0);
    assert_near(result.state.moment_values()[0], 3.5, 1.0e-12);
}

void test_diagnostics_are_aggregated()
{
    const std::vector<double> moments = raw_moments_1d({1.0}, {2.0}, 2);
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>{xi[0]}; });
    const auto zero_source = [](double, const qmomkit::QmomSourceContext& state)
    { return qmomkit::QmomMomentRate{std::vector<double>(state.size(), 0.0)}; };

    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options = tight_node_options();
    options.source_options.ssp_rk3_options.estimate_error = true;
    const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(state, 0.0, 0.5, velocity, zero_source, options);

    assert(result.ok());
    assert(result.diagnostics.post_source_qmom_status == qmomkit::QmomStatus::success);
    assert(result.diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::success);
    assert(result.diagnostics.second_advection_status == qmomkit::NodeAdvectionStatus::success);
    assert(result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(result.diagnostics.moment_source_diagnostics.error_estimate.has_value());
    assert(result.diagnostics.total_node_advection_steps ==
           result.diagnostics.first_advection_diagnostics.step_count + result.diagnostics.second_advection_diagnostics.step_count);
    assert(result.diagnostics.total_node_advection_rhs_evaluations > 0);
}

void test_qmom_phase_space_evolution_node_advection_scheme_selection()
{
    const std::vector<double> moments = raw_moments_1d({-1.0, 2.0}, {0.25, 0.75}, 4);
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 2);
    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{0.5}; });
    const auto zero_source = [](double, const qmomkit::QmomSourceContext& state)
    { return qmomkit::QmomMomentRate{std::vector<double>(state.size(), 0.0)}; };

    qmomkit::QmomPhaseSpaceEvolutionOptions ssp_options;
    ssp_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;
    ssp_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::QmomPhaseSpaceEvolutionResult ssp = qmomkit::advance_qmom_phase_space(state, 0.0, 0.2, velocity, zero_source, ssp_options);
    assert(ssp.ok());
    assert(ssp.diagnostics.first_advection_diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    assert(ssp.diagnostics.first_advection_diagnostics.step_count == 1);
    assert(ssp.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(ssp.diagnostics.total_node_advection_rejected_steps == 0);

    qmomkit::QmomPhaseSpaceEvolutionOptions rk_options = ssp_options;
    rk_options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    rk_options.node_options.adaptive_options.relative_tolerance = 1.0e-10;
    rk_options.node_options.adaptive_options.absolute_tolerance = 1.0e-12;
    rk_options.node_options.adaptive_options.initial_step = 0.01;
    const qmomkit::QmomPhaseSpaceEvolutionResult rk = qmomkit::advance_qmom_phase_space(state, 0.0, 0.2, velocity, zero_source, rk_options);
    assert(rk.ok());
    assert(rk.diagnostics.first_advection_diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(rk.diagnostics.first_advection_diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::success);
    assert(rk.diagnostics.total_node_advection_steps > 0);

    assert_vector_near(ssp.state.moment_values(), rk.state.moment_values(), 1.0e-10);

    qmomkit::QmomPhaseSpaceEvolutionOptions rk45_options = rk_options;
    rk45_options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk45;
    const qmomkit::QmomPhaseSpaceEvolutionResult rk45 = qmomkit::advance_qmom_phase_space(state, 0.0, 0.2, velocity, zero_source, rk45_options);
    assert(rk45.ok());
    assert(rk45.diagnostics.first_advection_diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk45);
    assert(rk45.diagnostics.first_advection_diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::success);
    assert_vector_near(rk.state.moment_values(), rk45.state.moment_values(), 1.0e-10);
}

void test_cqmom_state_phase_space_evolution_node_advection_scheme_selection()
{
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = fixed_cqmom_prepared_case();
    qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.25); });
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    qmomkit::CqmomPhaseSpaceEvolutionOptions ssp_options;
    ssp_options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;
    ssp_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomPhaseSpaceEvolutionResult ssp = qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, velocity, zero_source, ssp_options);
    assert(ssp.ok());
    assert(!ssp.diagnostic_quadrature.has_value());
    assert(ssp.diagnostics.first_advection_diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    assert(ssp.diagnostics.first_advection_diagnostics.step_count == 1);
    assert(ssp.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);

    qmomkit::CqmomPhaseSpaceEvolutionOptions rk_options = ssp_options;
    rk_options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    rk_options.node_options.adaptive_options.relative_tolerance = 1.0e-10;
    rk_options.node_options.adaptive_options.absolute_tolerance = 1.0e-12;
    rk_options.node_options.adaptive_options.initial_step = 0.01;
    const qmomkit::CqmomPhaseSpaceEvolutionResult rk = qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, velocity, zero_source, rk_options);
    assert(rk.ok());
    assert(!rk.diagnostic_quadrature.has_value());
    assert(rk.diagnostics.first_advection_diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(rk.diagnostics.first_advection_diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::success);

    qmomkit::CqmomPhaseSpaceEvolutionOptions rk45_options = rk_options;
    rk45_options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk45;
    const qmomkit::CqmomPhaseSpaceEvolutionResult rk45 = qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, velocity, zero_source, rk45_options);
    assert(rk45.ok());
    assert(!rk45.diagnostic_quadrature.has_value());
    assert(rk45.diagnostics.first_advection_diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk45);
    assert(rk45.diagnostics.first_advection_diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::success);
    assert_cqmom_quadrature_near(rk.state.quadrature(), rk45.state.quadrature(), 1.0e-10);
}

void test_failure_propagation_from_advection_operator()
{
    const std::vector<double> moments{1.0, 2.0};
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    const qmomkit::PhaseSpaceVelocity bad_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0}; });
    const auto zero_source = [](double, const qmomkit::QmomSourceContext& state)
    { return qmomkit::QmomMomentRate{std::vector<double>(state.size(), 0.0)}; };
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(state, 0.0, 0.1, bad_velocity, zero_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::advection_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::first_advection_half_step);
    assert(result.diagnostics.first_advection_status == qmomkit::NodeAdvectionStatus::velocity_failed);
    assert_vector_near(result.state.moment_values(), state.moment_values(), 0.0);
}

void test_failure_propagation_from_source_operator()
{
    const std::vector<double> moments{1.0, 2.0};
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto bad_source = [](double, const qmomkit::QmomSourceContext&) { return qmomkit::QmomMomentRate{std::vector<double>{1.0, 2.0, 3.0}}; };
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::QmomPhaseSpaceEvolutionResult result = qmomkit::advance_qmom_phase_space(state, 0.0, 0.1, zero_velocity, bad_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    assert(result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert_vector_near(result.state.moment_values(), state.moment_values(), 0.0);
}

void test_failure_propagation_from_boundary_contribution()
{
    const std::vector<double> moments{1.0, 2.0};
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto zero_source = [](double, const qmomkit::QmomSourceContext& state)
    { return qmomkit::QmomMomentRate{std::vector<double>(state.size(), 0.0)}; };
    const qmomkit::QmomMomentBoundaryContribution bad_boundary = [](double, const qmomkit::QmomSourceContext& state)
    { return qmomkit::QmomMomentRate{std::vector<double>(state.size() + 1, 0.0)}; };
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::QmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_qmom_phase_space(state, 0.0, 0.1, zero_velocity, zero_source, bad_boundary, options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    assert(result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.message.find("boundary contribution") != std::string::npos);
    assert_vector_near(result.state.moment_values(), state.moment_values(), 0.0);
}

void test_boundary_contribution_post_source_inversion_failure()
{
    const std::vector<double> moments{1.0, 0.0};
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto zero_source = [](double, const qmomkit::QmomSourceContext& state)
    { return qmomkit::QmomMomentRate{std::vector<double>(state.size(), 0.0)}; };
    const qmomkit::QmomMomentBoundaryContribution removal = [](double, const qmomkit::QmomSourceContext& state)
    {
        std::vector<double> rate(state.size(), 0.0);
        rate.front() = -20.0;
        return qmomkit::QmomMomentRate{rate};
    };
    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::QmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_qmom_phase_space(state, 0.0, 0.1, zero_velocity, zero_source, removal, options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion);
    assert(result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(result.diagnostics.post_source_qmom_status != qmomkit::QmomStatus::success);
    assert_vector_near(result.state.moment_values(), state.moment_values(), 0.0);
}

// The packed velocity form is an alternative spelling of the same contract, so
// a full phase-space step must agree with the per-node form exactly, including
// diagnostics, for every splitting scheme.
void test_packed_velocity_phase_space_step_matches_per_node_form()
{
    const qmomkit::PhaseSpaceVelocity adapted = qmomkit::make_phase_space_velocity(
        [](const double time, const std::vector<double>& xi)
        {
            std::vector<double> velocity(xi.size(), 0.0);
            for (std::size_t dimension = 0; dimension < xi.size(); ++dimension)
            {
                velocity[dimension] = 0.2 * xi[dimension] - 0.05 * time;
            }
            return velocity;
        });
    const qmomkit::PhaseSpaceVelocity packed =
        [](const double time, const std::vector<double>& coordinates, std::vector<double>& velocities, const std::size_t dimensions)
    {
        for (std::size_t entry = 0; entry < coordinates.size(); ++entry)
        {
            velocities[entry] = 0.2 * coordinates[entry] - 0.05 * time;
        }
        (void)dimensions;
    };

    for (const qmomkit::PhaseSpaceEvolutionScheme scheme :
         {qmomkit::PhaseSpaceEvolutionScheme::strang, qmomkit::PhaseSpaceEvolutionScheme::advsrc, qmomkit::PhaseSpaceEvolutionScheme::srcadv})
    {
        const qmomkit::CqmomState cqmom_state = fixed_cqmom_prepared_case().initial_state();
        const qmomkit::CqmomMomentSource cqmom_source = [](double, const qmomkit::CqmomSourceContext& context)
        { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

        qmomkit::CqmomPhaseSpaceEvolutionOptions cqmom_options;
        cqmom_options.scheme = scheme;

        const qmomkit::CqmomPhaseSpaceEvolutionResult cqmom_adapted =
            qmomkit::advance_cqmom_phase_space(cqmom_state, 0.25, 0.05, adapted, cqmom_source, cqmom_options);
        const qmomkit::CqmomPhaseSpaceEvolutionResult cqmom_packed =
            qmomkit::advance_cqmom_phase_space(cqmom_state, 0.25, 0.05, packed, cqmom_source, cqmom_options);

        assert(cqmom_adapted.ok());
        assert(cqmom_packed.ok());
        assert(cqmom_packed.state.moment_values() == cqmom_adapted.state.moment_values());
        assert(cqmom_packed.final_time == cqmom_adapted.final_time);
        assert(cqmom_packed.diagnostics.total_node_advection_rhs_evaluations == cqmom_adapted.diagnostics.total_node_advection_rhs_evaluations);
        assert(cqmom_packed.diagnostics.total_node_advection_steps == cqmom_adapted.diagnostics.total_node_advection_steps);

        const qmomkit::QmomState qmom_state = qmom_state_from_moments({1.0, 0.5, 0.75, 0.875, 1.5, 3.0}, 3);
        const qmomkit::QmomMomentSource qmom_source = [](double, const qmomkit::QmomSourceContext& context)
        { return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

        qmomkit::QmomPhaseSpaceEvolutionOptions qmom_options;
        qmom_options.scheme = scheme;

        const qmomkit::QmomPhaseSpaceEvolutionResult qmom_adapted =
            qmomkit::advance_qmom_phase_space(qmom_state, 0.25, 0.05, adapted, qmom_source, qmom_options);
        const qmomkit::QmomPhaseSpaceEvolutionResult qmom_packed =
            qmomkit::advance_qmom_phase_space(qmom_state, 0.25, 0.05, packed, qmom_source, qmom_options);

        assert(qmom_adapted.ok());
        assert(qmom_packed.ok());
        assert(qmom_packed.state.moment_values() == qmom_adapted.state.moment_values());
        assert(qmom_packed.state.quadrature().abscissas == qmom_adapted.state.quadrature().abscissas);
        assert(qmom_packed.diagnostics.total_node_advection_rhs_evaluations == qmom_adapted.diagnostics.total_node_advection_rhs_evaluations);
    }
}

// Advection-only entry points share the same adaptation, including the
// missing-callback contract.
void test_packed_velocity_advection_entry_points_match_and_reject_empty_callback()
{
    const qmomkit::CqmomState cqmom_state = fixed_cqmom_prepared_case().initial_state();
    const qmomkit::PhaseSpaceVelocity adapted =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.1); });
    const qmomkit::PhaseSpaceVelocity packed = [](double, const std::vector<double>& coordinates, std::vector<double>& velocities, std::size_t)
    { velocities.assign(coordinates.size(), 0.1); };

    const qmomkit::CqmomPhaseSpaceEvolutionResult from_adapted = qmomkit::advance_cqmom_advection(cqmom_state, 0.0, 0.1, adapted);
    const qmomkit::CqmomPhaseSpaceEvolutionResult from_packed = qmomkit::advance_cqmom_advection(cqmom_state, 0.0, 0.1, packed);
    assert(from_adapted.ok());
    assert(from_packed.ok());
    assert(from_packed.state.moment_values() == from_adapted.state.moment_values());

    const qmomkit::CqmomPhaseSpaceEvolutionResult missing = qmomkit::advance_cqmom_advection(cqmom_state, 0.0, 0.1, qmomkit::PhaseSpaceVelocity{});
    assert(!missing.ok());
    assert(missing.status == qmomkit::PhaseSpaceEvolutionStatus::invalid_input);
}

void test_cqmom_case_failure_propagation_from_inversion()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::TestCqmomPhaseSpaceEvolutionState state{
        qmomkit::CqmomMomentData{layout, std::vector<double>(layout.required_indices.size(), 0.0)}};
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    const qmomkit::TestCqmomPhaseSpaceEvolutionResult result = qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::initial_inversion);
    assert(result.diagnostics.initial_cqmom_status != qmomkit::CqmomStatus::success);
}

void test_cqmom_case_phase_space_evolution_preserves_invalid_state_messages()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    qmomkit::TestCqmomPhaseSpaceEvolutionState missing_values = cqmom_state_from_tensor(moments, hierarchy);
    missing_values.moments.values.pop_back();
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });

    const qmomkit::TestCqmomPhaseSpaceEvolutionResult missing = qmomkit::advance_cqmom_phase_space(missing_values, 0.0, 0.1, zero_velocity);
    assert(!missing.ok());
    assert(missing.status == qmomkit::PhaseSpaceEvolutionStatus::invalid_input);
    assert(missing.message == "CQMOM phase-space evolution moments must match the fixed layout size.");
    assert(missing.diagnostics.message == missing.message);

    qmomkit::TestCqmomPhaseSpaceEvolutionState wrong_layout = cqmom_state_from_tensor(moments, hierarchy);
    wrong_layout.moments.layout.required_indices.pop_back();
    const qmomkit::TestCqmomPhaseSpaceEvolutionResult layout = qmomkit::advance_cqmom_phase_space(wrong_layout, 0.0, 0.1, zero_velocity);
    assert(!layout.ok());
    assert(layout.status == qmomkit::PhaseSpaceEvolutionStatus::invalid_input);
    assert(layout.message == "CQMOM phase-space evolution moments must match the fixed layout size.");
    assert(layout.diagnostics.message == layout.message);

    qmomkit::TestCqmomPhaseSpaceEvolutionState nonfinite_values = cqmom_state_from_tensor(moments, hierarchy);
    nonfinite_values.moments.values[2] = std::numeric_limits<double>::quiet_NaN();
    const qmomkit::TestCqmomPhaseSpaceEvolutionResult nonfinite = qmomkit::advance_cqmom_phase_space(nonfinite_values, 0.0, 0.1, zero_velocity);
    assert(!nonfinite.ok());
    assert(nonfinite.status == qmomkit::PhaseSpaceEvolutionStatus::invalid_input);
    assert(nonfinite.message == "CQMOM phase-space evolution moments must be finite.");
    assert(nonfinite.diagnostics.message == nonfinite.message);
}

void test_cqmom_case_failure_propagation_from_reconstruction()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestCqmomPhaseSpaceEvolutionState state = cqmom_state_from_tensor(moments, hierarchy);
    const qmomkit::PhaseSpaceVelocity huge_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 1.0e200); });

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options;
    options.node_options.adaptive_options.initial_step = 1.0;
    options.node_options.adaptive_options.max_step = 1.0;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::advsrc;
    const qmomkit::TestCqmomPhaseSpaceEvolutionResult result = qmomkit::advance_cqmom_phase_space(state, 0.0, 1.0, huge_velocity, options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_reconstruction_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::first_moment_reconstruction);
    assert(!result.diagnostics.first_reconstruction_diagnostics.success);
}

void test_cqmom_case_failure_propagation_from_source_operator()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestCqmomPhaseSpaceEvolutionState state = cqmom_state_from_tensor(moments, hierarchy);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto bad_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; };

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::TestCqmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, bad_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::moment_source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::moment_source_step);
    assert(result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
}

void test_cqmom_case_failure_propagation_from_post_source_inversion()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestCqmomPhaseSpaceEvolutionState state = cqmom_state_from_tensor(moments, hierarchy);
    const qmomkit::PhaseSpaceVelocity zero_velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); });
    const auto bad_realizability_source = [](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        rate.front() = -20.0;
        return qmomkit::CqmomMomentRate{rate};
    };

    qmomkit::TestCqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::srcadv;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::TestCqmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_cqmom_phase_space(state, 0.0, 0.1, zero_velocity, bad_realizability_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::PhaseSpaceEvolutionStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion);
    assert(result.diagnostics.moment_source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(result.diagnostics.post_source_cqmom_status != qmomkit::CqmomStatus::success);
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"explicit_cqmom_adaptive_branch_options_are_honored", test_explicit_cqmom_adaptive_branch_options_are_honored},
            {"explicit_cqmom_diagnostic_level_none_is_propagated", test_explicit_cqmom_diagnostic_level_none_is_propagated},
            {"explicit_cqmom_adaptive_options_are_honored", test_explicit_cqmom_adaptive_options_are_honored},
            {"explicit_cqmom_source_options_are_honored", test_explicit_cqmom_source_options_are_honored},
            {"explicit_cqmom_native_source_context_contract", test_explicit_cqmom_native_source_context_contract},
            {"explicit_cqmom_native_source_failures_are_reported", test_explicit_cqmom_native_source_failures_are_reported},
            {"explicit_cqmom_native_source_invalid_values_are_reported", test_explicit_cqmom_native_source_invalid_values_are_reported},
            {"explicit_cqmom_source_growth_is_automatic", test_explicit_cqmom_source_growth_is_automatic},
            {"fixed_cqmom_phase_space_evolution_successors_share_prepared_metadata",
             test_fixed_cqmom_phase_space_evolution_successors_share_prepared_metadata},
            {"adaptive_cqmom_phase_space_evolution_successors_share_prepared_metadata",
             test_adaptive_cqmom_phase_space_evolution_successors_share_prepared_metadata},
            {"public_qmom_advection_stage_advances_canonical_state", test_public_qmom_advection_stage_advances_canonical_state},
            {"public_qmom_advection_stage_failure_rolls_back", test_public_qmom_advection_stage_failure_rolls_back},
            {"qmom_bounded_state_advection_and_source_reconstruction_semantics",
             test_qmom_bounded_state_advection_and_source_reconstruction_semantics},
            {"qmom_bounded_post_source_automatic_scaling_is_unit_covariant", test_qmom_bounded_post_source_automatic_scaling_is_unit_covariant},
            {"public_qmom_moment_source_stage_fixed_mode_is_strict", test_public_qmom_moment_source_stage_fixed_mode_is_strict},
            {"public_qmom_moment_source_stage_adaptive_policy", test_public_qmom_moment_source_stage_adaptive_policy},
            {"public_qmom_moment_source_stage_failures_roll_back", test_public_qmom_moment_source_stage_failures_roll_back},
            {"qmom_source_context_uses_current_rk_stage_and_caches_quadrature", test_qmom_source_context_uses_current_rk_stage_and_caches_quadrature},
            {"qmom_boundary_context_matches_source_context", test_qmom_boundary_context_matches_source_context},
            {"qmom_source_context_generalized_moment_semantics", test_qmom_source_context_generalized_moment_semantics},
            {"qmom_source_context_generalized_moment_uses_current_rk_stage", test_qmom_source_context_generalized_moment_uses_current_rk_stage},
            {"qmom_source_context_visits_nodes_from_each_rk_stage", test_qmom_source_context_visits_nodes_from_each_rk_stage},
            {"configured_qmom_source_realizability_policy_matches_programmatic_adapter",
             test_configured_qmom_source_realizability_policy_matches_programmatic_adapter},
            {"configured_qmom_source_realizability_policy_enforces_its_support",
             test_configured_qmom_source_realizability_policy_enforces_its_support},
            {"public_qmom_moment_source_stage_boundary_matches_source_update", test_public_qmom_moment_source_stage_boundary_matches_source_update},
            {"public_qmom_stage_composition_matches_strang_phase_space_evolution",
             test_public_qmom_stage_composition_matches_strang_phase_space_evolution},
            {"public_qmom_stage_composition_matches_lie_phase_space_evolution", test_public_qmom_stage_composition_matches_lie_phase_space_evolution},
            {"qmom_phase_space_evolution_failure_after_first_stage_rolls_back_to_macro_input",
             test_qmom_phase_space_evolution_failure_after_first_stage_rolls_back_to_macro_input},
            {"qmom_macro_rollbacks_preserve_the_original_state_after_later_successes",
             test_qmom_macro_rollbacks_preserve_the_original_state_after_later_successes},
            {"qmom_failed_steps_report_the_time_of_the_returned_state", test_qmom_failed_steps_report_the_time_of_the_returned_state},
            {"cqmom_initial_state_stores_accepted_quadrature", test_cqmom_initial_state_stores_accepted_quadrature},
            {"cqmom_initial_state_canonicalizes_adaptively_reduced_moments", test_cqmom_initial_state_canonicalizes_adaptively_reduced_moments},
            {"cqmom_case_initial_state_matches_prepared_case", test_cqmom_case_initial_state_matches_prepared_case},
            {"cqmom_case_initial_state_reports_failures", test_cqmom_case_initial_state_reports_failures},
            {"cqmom_first_advection_stage_uses_initial_accepted_quadrature", test_cqmom_first_advection_stage_uses_initial_accepted_quadrature},
            {"public_cqmom_advection_stage_preserves_topology", test_public_cqmom_advection_stage_preserves_topology},
            {"public_cqmom_moment_source_stage_adaptive_policy", test_public_cqmom_moment_source_stage_adaptive_policy},
            {"public_cqmom_moment_source_stage_growth_is_limited_to_one_node", test_public_cqmom_moment_source_stage_growth_is_limited_to_one_node},
            {"public_cqmom_moment_source_stage_canonicalizes_adaptive_reduction",
             test_public_cqmom_moment_source_stage_canonicalizes_adaptive_reduction},
            {"public_cqmom_stage_composition_matches_phase_space_evolution", test_public_cqmom_stage_composition_matches_phase_space_evolution},
            {"cqmom_phase_space_evolution_success_returns_no_diagnostic_quadrature",
             test_cqmom_phase_space_evolution_success_returns_no_diagnostic_quadrature},
            {"adaptive_cqmom_phase_space_evolution_state_counts_match_accepted_quadrature",
             test_adaptive_cqmom_phase_space_evolution_state_counts_match_accepted_quadrature},
            {"public_cqmom_stage_failures_roll_back", test_public_cqmom_stage_failures_roll_back},
            {"cqmom_phase_space_evolution_failure_after_first_stage_rolls_back_to_macro_input",
             test_cqmom_phase_space_evolution_failure_after_first_stage_rolls_back_to_macro_input},
            {"cqmom_macro_rollbacks_preserve_the_original_state_after_later_successes",
             test_cqmom_macro_rollbacks_preserve_the_original_state_after_later_successes},
            {"cqmom_failed_steps_report_the_time_of_the_returned_state", test_cqmom_failed_steps_report_the_time_of_the_returned_state},
            {"cqmom_macro_success_preserves_input_and_final_value_semantics", test_cqmom_macro_success_preserves_input_and_final_value_semantics},
            {"public_cqmom_boundary_moment_source_stage_matches_phase_space_evolution",
             test_public_cqmom_boundary_moment_source_stage_matches_phase_space_evolution},
            {"public_cqmom_stage_irrelevant_options_are_ignored", test_public_cqmom_stage_irrelevant_options_are_ignored},
            {"advection_only_matches_direct_operator_qmom", test_advection_only_matches_direct_operator_qmom},
            {"source_only_matches_direct_operator_qmom", test_source_only_matches_direct_operator_qmom},
            {"source_plus_boundary_matches_direct_operator_qmom", test_source_plus_boundary_matches_direct_operator_qmom},
            {"adaptive_branch_specific_cqmom_case", test_adaptive_branch_specific_cqmom_case},
            {"adaptive_state_initialization_preserves_policy", test_adaptive_state_initialization_preserves_policy},
            {"adaptive_state_source_phase_space_evolution_changes_active_counts",
             test_adaptive_state_source_phase_space_evolution_changes_active_counts},
            {"adaptive_state_preserves_branch_counts_across_steps", test_adaptive_state_preserves_branch_counts_across_steps},
            {"adaptive_state_preserved_source_can_reduce_active_counts", test_adaptive_state_preserved_source_can_reduce_active_counts},
            {"adaptive_state_source_policy_allows_growth_to_configured_maximum",
             test_adaptive_state_source_policy_allows_growth_to_configured_maximum},
            {"adaptive_state_source_growth_forces_adaptive_wheeler_without_mutating_qmom_options",
             test_adaptive_state_source_growth_forces_adaptive_wheeler_without_mutating_qmom_options},
            {"adaptive_state_source_growth_canonicalizes_between_steps", test_adaptive_state_source_growth_canonicalizes_between_steps},
            {"adaptive_state_failed_growth_trial_keeps_current_count", test_adaptive_state_failed_growth_trial_keeps_current_count},
            {"adaptive_state_bounded_support_rejects_post_source_growth", test_adaptive_state_bounded_support_rejects_post_source_growth},
            {"adaptive_state_reduction_prevents_source_growth_trial", test_adaptive_state_reduction_prevents_source_growth_trial},
            {"adaptive_state_parent_growth_creates_unconstrained_descendants", test_adaptive_state_parent_growth_creates_unconstrained_descendants},
            {"adaptive_state_preserved_source_invalidates_descendants_after_root_collapse",
             test_adaptive_state_preserved_source_invalidates_descendants_after_root_collapse},
            {"adaptive_state_preserved_source_invalidates_grandchildren_after_second_level_collapse",
             test_adaptive_state_preserved_source_invalidates_grandchildren_after_second_level_collapse},
            {"adaptive_state_initial_counts_preserve_nonidentity_conditioning_order",
             test_adaptive_state_initial_counts_preserve_nonidentity_conditioning_order},
            {"fixed_state_phase_space_evolution_does_not_store_adaptive_branch_counts",
             test_fixed_state_phase_space_evolution_does_not_store_adaptive_branch_counts},
            {"adaptive_state_rejects_duplicate_active_count_metadata", test_adaptive_state_rejects_duplicate_active_count_metadata},
            {"adaptive_state_rejects_missing_active_count_metadata", test_adaptive_state_rejects_missing_active_count_metadata},
            {"adaptive_state_failure_preserves_active_counts_for_retry", test_adaptive_state_failure_preserves_active_counts_for_retry},
            {"adaptive_state_post_source_failure_preserves_active_counts_for_retry",
             test_adaptive_state_post_source_failure_preserves_active_counts_for_retry},
            {"adaptive_state_phase_space_evolution_reports_initial_failure", test_adaptive_state_phase_space_evolution_reports_initial_failure},
            {"state_cqmom_phase_space_evolution_matches_cqmom_case_fixed_layout",
             test_state_cqmom_phase_space_evolution_matches_cqmom_case_fixed_layout},
            {"cqmom_source_context_reconstructs_untracked_moments_from_each_rk_stage",
             test_cqmom_source_context_reconstructs_untracked_moments_from_each_rk_stage},
            {"cqmom_source_context_visits_terminal_nodes_from_each_rk_stage", test_cqmom_source_context_visits_terminal_nodes_from_each_rk_stage},
            {"cqmom_source_context_generalized_moments_use_stage_quadrature", test_cqmom_source_context_generalized_moments_use_stage_quadrature},
            {"cqmom_source_context_generalized_moment_validation", test_cqmom_source_context_generalized_moment_validation},
            {"state_cqmom_phase_space_evolution_preserves_vandermonde_solver_policy",
             test_state_cqmom_phase_space_evolution_preserves_vandermonde_solver_policy},
            {"cqmom_phase_space_evolution_vandermonde_policy_independent_of_diagnostic_level",
             test_cqmom_phase_space_evolution_vandermonde_policy_independent_of_diagnostic_level},
            {"state_cqmom_source_evolution_matches_cqmom_case", test_state_cqmom_source_evolution_matches_cqmom_case},
            {"state_cqmom_failure_propagation_from_source_operator", test_state_cqmom_failure_propagation_from_source_operator},
            {"state_cqmom_phase_space_evolution_api_does_not_expose_moments", test_state_cqmom_phase_space_evolution_api_does_not_expose_moments},
            {"state_qmom_phase_space_evolution_api_does_not_expose_moments", test_state_qmom_phase_space_evolution_api_does_not_expose_moments},
            {"qmom_state_copy_move_workflow", test_qmom_state_copy_move_workflow},
            {"qmom_initial_state_canonicalizes_adaptive_moments", test_qmom_initial_state_canonicalizes_adaptive_moments},
            {"qmom_adaptive_phase_space_evolution_preserves_active_count_with_zero_source",
             test_qmom_adaptive_phase_space_evolution_preserves_active_count_with_zero_source},
            {"qmom_adaptive_post_source_reduction_and_growth_policy", test_qmom_adaptive_post_source_reduction_and_growth_policy},
            {"qmom_fixed_post_source_requires_configured_count", test_qmom_fixed_post_source_requires_configured_count},
            {"first_order_schemes_differ_for_noncommuting_qmom", test_first_order_schemes_differ_for_noncommuting_qmom},
            {"strang_source_uses_physical_time_interval_for_nonautonomous_source",
             test_strang_source_uses_physical_time_interval_for_nonautonomous_source},
            {"diagnostics_are_aggregated", test_diagnostics_are_aggregated},
            {"qmom_phase_space_evolution_node_advection_scheme_selection", test_qmom_phase_space_evolution_node_advection_scheme_selection},
            {"cqmom_state_phase_space_evolution_node_advection_scheme_selection",
             test_cqmom_state_phase_space_evolution_node_advection_scheme_selection},
            {"failure_propagation_from_advection_operator", test_failure_propagation_from_advection_operator},
            {"failure_propagation_from_source_operator", test_failure_propagation_from_source_operator},
            {"failure_propagation_from_boundary_contribution", test_failure_propagation_from_boundary_contribution},
            {"boundary_contribution_post_source_inversion_failure", test_boundary_contribution_post_source_inversion_failure},
            {"packed_velocity_phase_space_step_matches_per_node_form", test_packed_velocity_phase_space_step_matches_per_node_form},
            {"packed_velocity_advection_entry_points_match_and_reject_empty_callback",
             test_packed_velocity_advection_entry_points_match_and_reject_empty_callback},
            {"cqmom_case_failure_propagation_from_inversion", test_cqmom_case_failure_propagation_from_inversion},
            {"cqmom_case_phase_space_evolution_preserves_invalid_state_messages",
             test_cqmom_case_phase_space_evolution_preserves_invalid_state_messages},
            {"cqmom_case_failure_propagation_from_reconstruction", test_cqmom_case_failure_propagation_from_reconstruction},
            {"cqmom_case_failure_propagation_from_source_operator", test_cqmom_case_failure_propagation_from_source_operator},
            {"cqmom_case_failure_propagation_from_post_source_inversion", test_cqmom_case_failure_propagation_from_post_source_inversion},
        });
}
