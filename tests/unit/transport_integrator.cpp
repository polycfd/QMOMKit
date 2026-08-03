#include "internal/cqmom_active_counts.hpp"
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

qmomkit::QmomState qmom_state_from_moments(const std::vector<double>& moments, const std::size_t node_count, const qmomkit::QmomOptions& options = {})
{
    return qmomkit::QmomCase{qmomkit::MomentSet(moments), node_count, options}.initial_state();
}

qmomkit::MomentSource qmom_source_from_to(const std::vector<double>& initial, const std::vector<double>& target)
{
    return [initial, target](double, const std::vector<double>& state)
    {
        assert(state.size() == initial.size());
        std::vector<double> rate(state.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            rate[i] = target[i] - initial[i];
        }
        return rate;
    };
}

void assert_qmom_state_canonical(const qmomkit::QmomState& state, const double tolerance)
{
    const qmomkit::Quadrature1D quadrature{state.quadrature().abscissas, state.quadrature().weights, qmomkit::QuadratureStatus::success, ""};
    const std::vector<double> reconstructed = qmomkit::reconstruct_raw_moments(quadrature, state.moments().size());
    assert_vector_near(state.moments(), reconstructed, tolerance);
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

qmomkit::RawMomentTensor one_dimensional_cqmom_moments(const std::size_t configured_node_count, const std::vector<double>& abscissas,
                                                       const std::vector<double>& weights)
{
    const std::vector<std::size_t> extents{2 * configured_node_count};
    return {extents, raw_moments_1d(abscissas, weights, extents.front())};
}

qmomkit::CqmomState one_dimensional_adaptive_cqmom_state(const qmomkit::CqmomMomentData& moments)
{
    qmomkit::CqmomLeanCase lean_case{moments, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, true};
    return qmomkit::internal::create_cqmom_state(lean_case);
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

struct TestLeanCqmomTransportState
{
    CqmomMomentData moments;
};

} // namespace qmomkit

qmomkit::TestLeanCqmomTransportState lean_state_from_dense(const qmomkit::RawMomentTensor& moments, const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    return qmomkit::TestLeanCqmomTransportState{qmomkit::extract_cqmom_moments(moments, layout)};
}

namespace qmomkit
{

struct TestDenseCqmomTransportOptions : CqmomTransportOptions
{
    ConditionalMomentOptions conditional_options;
    QmomOptions qmom_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::detailed;
    VandermondeSolver vandermonde_solver = VandermondeSolver::automatic;
};

struct TestLeanCqmomTransportOptions : CqmomTransportOptions
{
    ConditionalMomentOptions conditional_options;
    QmomOptions qmom_options;
    bool adaptive_branch_counts = false;
    CqmomDiagnosticLevel diagnostic_level = CqmomDiagnosticLevel::detailed;
    VandermondeSolver vandermonde_solver = VandermondeSolver::automatic;
};

struct TestDenseCqmomTransportResult
{
    RawMomentTensor moments;
    std::optional<CqmomQuadrature> quadrature;
    TransportStatus status = TransportStatus::success;
    TransportDiagnostics diagnostics;
    double final_time = 0.0;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == TransportStatus::success;
    }
};

struct TestLeanCqmomTransportResult
{
    TestLeanCqmomTransportState state;
    std::optional<CqmomQuadrature> quadrature;
    TransportStatus status = TransportStatus::success;
    TransportDiagnostics diagnostics;
    double final_time = 0.0;
    std::string message;

    [[nodiscard]] bool ok() const noexcept
    {
        return status == TransportStatus::success;
    }
};

RawMomentTensor dense_moments_from_state(const CqmomState& state, const std::vector<std::size_t>& extents)
{
    return RawMomentTensor{extents, state.moment_values()};
}

TestDenseCqmomTransportResult advance_transport_cqmom(const RawMomentTensor& moments, const PhaseSpaceHierarchy& hierarchy, const double t,
                                                      const double time_step, const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                      const TestDenseCqmomTransportOptions& options)
{
    CqmomDenseCase dense_case{hierarchy,
                              moments,
                              options.conditional_options,
                              options.qmom_options,
                              options.adaptive_branch_counts,
                              options.diagnostic_level,
                              options.vandermonde_solver};
    CqmomState state = internal::create_cqmom_state(dense_case);
    const CqmomTransportResult result =
        advance_transport_cqmom(state, t, time_step, velocity, source, static_cast<const CqmomTransportOptions&>(options));
    const CqmomState& result_state = result.ok() ? result.state : state;
    return TestDenseCqmomTransportResult{dense_moments_from_state(result_state, moments.extents),
                                         result.quadrature,
                                         result.status,
                                         result.diagnostics,
                                         result.final_time,
                                         result.message};
}

TestDenseCqmomTransportResult advance_transport_cqmom(const RawMomentTensor& moments, const PhaseSpaceHierarchy& hierarchy, const double t,
                                                      const double time_step, const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                      const CqmomMomentBoundaryContribution& boundary, const TestDenseCqmomTransportOptions& options)
{
    CqmomDenseCase dense_case{hierarchy,
                              moments,
                              options.conditional_options,
                              options.qmom_options,
                              options.adaptive_branch_counts,
                              options.diagnostic_level,
                              options.vandermonde_solver};
    CqmomState state = internal::create_cqmom_state(dense_case);
    const CqmomTransportResult result =
        advance_transport_cqmom(state, t, time_step, velocity, source, boundary, static_cast<const CqmomTransportOptions&>(options));
    const CqmomState& result_state = result.ok() ? result.state : state;
    return TestDenseCqmomTransportResult{dense_moments_from_state(result_state, moments.extents),
                                         result.quadrature,
                                         result.status,
                                         result.diagnostics,
                                         result.final_time,
                                         result.message};
}

TestLeanCqmomTransportResult advance_transport_cqmom(const TestLeanCqmomTransportState& lean_state, const double t, const double time_step,
                                                     const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                     const TestLeanCqmomTransportOptions& options = {})
{
    CqmomLeanCase lean_case{lean_state.moments,       options.conditional_options, options.qmom_options, options.adaptive_branch_counts,
                            options.diagnostic_level, options.vandermonde_solver};
    CqmomState state = internal::create_cqmom_state(lean_case);
    const CqmomTransportResult result =
        advance_transport_cqmom(state, t, time_step, velocity, source, static_cast<const CqmomTransportOptions&>(options));
    const CqmomState& result_state = result.ok() ? result.state : state;
    return TestLeanCqmomTransportResult{TestLeanCqmomTransportState{CqmomMomentData{lean_state.moments.layout, result_state.moment_values()}},
                                        result.quadrature,
                                        result.status,
                                        result.diagnostics,
                                        result.final_time,
                                        result.message};
}

TestLeanCqmomTransportResult advance_transport_cqmom(const TestLeanCqmomTransportState& lean_state, const double t, const double time_step,
                                                     const PhaseSpaceVelocity& velocity, const CqmomMomentSource& source,
                                                     const CqmomMomentBoundaryContribution& boundary, const TestLeanCqmomTransportOptions& options)
{
    CqmomLeanCase lean_case{lean_state.moments,       options.conditional_options, options.qmom_options, options.adaptive_branch_counts,
                            options.diagnostic_level, options.vandermonde_solver};
    CqmomState state = internal::create_cqmom_state(lean_case);
    const CqmomTransportResult result =
        advance_transport_cqmom(state, t, time_step, velocity, source, boundary, static_cast<const CqmomTransportOptions&>(options));
    const CqmomState& result_state = result.ok() ? result.state : state;
    return TestLeanCqmomTransportResult{TestLeanCqmomTransportState{CqmomMomentData{lean_state.moments.layout, result_state.moment_values()}},
                                        result.quadrature,
                                        result.status,
                                        result.diagnostics,
                                        result.final_time,
                                        result.message};
}

TestLeanCqmomTransportResult advance_transport_cqmom(const TestLeanCqmomTransportState& lean_state, const double t, const double time_step,
                                                     const PhaseSpaceVelocity& velocity, const TestLeanCqmomTransportOptions& options = {})
{
    const CqmomMomentSource zero_source = [](double, const CqmomSourceContext& context)
    { return CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };
    return advance_transport_cqmom(lean_state, t, time_step, velocity, zero_source, options);
}

} // namespace qmomkit

namespace
{

qmomkit::CqmomCase fixed_cqmom_case()
{
    const auto parsed = qmomkit::config::parse_cqmom(fixed_cqmom_case_text, "fixed-cqmom-case");
    assert(parsed.ok());
    return parsed.value();
}

qmomkit::CqmomCase fixed_cqmom_case_with_vandermonde_solver(const char* const solver)
{
    std::string text = fixed_cqmom_case_text;
    text += "\n[vandermonde]\nsolver = ";
    text += solver;
    text += "\n";
    const auto parsed = qmomkit::config::parse_cqmom(text, "fixed-cqmom-case-with-vandermonde-solver");
    assert(parsed.ok());
    return parsed.value();
}

qmomkit::CqmomCase fixed_cqmom_case_with_storage(const char* const storage)
{
    std::string text = fixed_cqmom_case_text;
    text += "\n[cqmom]\nmoment_storage = ";
    text += storage;
    text += "\n";
    const auto parsed = qmomkit::config::parse_cqmom(text, "fixed-cqmom-case-with-storage");
    assert(parsed.ok());
    return parsed.value();
}

void test_internal_cqmom_state_backend_accessors_are_explicit_and_guarded()
{
    const qmomkit::CqmomState dense_state = fixed_cqmom_case_with_storage("dense").initial_state();
    const qmomkit::CqmomState lean_state = fixed_cqmom_case_with_storage("lean").initial_state();

    assert(dense_state.moment_storage() == qmomkit::CqmomMomentStorage::dense);
    assert(lean_state.moment_storage() == qmomkit::CqmomMomentStorage::lean);

    assert(qmomkit::internal::cqmom_state_dense_case(dense_state) != nullptr);
    assert(qmomkit::internal::cqmom_state_lean_case(dense_state) == nullptr);
    assert(qmomkit::internal::cqmom_state_dense_case(lean_state) == nullptr);
    assert(qmomkit::internal::cqmom_state_lean_case(lean_state) != nullptr);
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
    return parsed.value();
}

qmomkit::CqmomCase adaptive_cqmom_case_from_moments(const qmomkit::RawMomentTensor& moments)
{
    const std::string text = cqmom_case_text_from_moments(moments, true, true);
    const auto parsed = qmomkit::config::parse_cqmom(text, "adaptive-cqmom-case-from-moments");
    assert(parsed.ok());
    return parsed.value();
}

qmomkit::CqmomCase adaptive_cqmom_case_from_moments_with_fixed_qmom(const qmomkit::RawMomentTensor& moments)
{
    const std::string text = cqmom_case_text_from_moments(moments, true, false);
    const auto parsed = qmomkit::config::parse_cqmom(text, "adaptive-cqmom-case-from-moments-with-fixed-qmom");
    assert(parsed.ok());
    return parsed.value();
}

qmomkit::CqmomCase adaptive_cqmom_case_from_moments_with_storage(const qmomkit::RawMomentTensor& moments, const char* const storage)
{
    std::string text = cqmom_case_text_from_moments(moments, true, true);
    text += "\n[cqmom]\nmoment_storage = ";
    text += storage;
    text += "\n";
    const auto parsed = qmomkit::config::parse_cqmom(text, "adaptive-cqmom-case-from-moments-with-storage");
    assert(parsed.ok());
    return parsed.value();
}

qmomkit::CqmomCase adaptive_cqmom_case_from_moments_with_hierarchy(const qmomkit::RawMomentTensor& moments,
                                                                   const std::vector<std::size_t>& node_counts,
                                                                   const std::vector<std::size_t>& conditioning_order)
{
    const std::string text = cqmom_case_text_from_moments_with_hierarchy(moments, node_counts, conditioning_order, true, true);
    const auto parsed = qmomkit::config::parse_cqmom(text, "adaptive-cqmom-case-from-moments-with-hierarchy");
    assert(parsed.ok());
    return parsed.value();
}

qmomkit::CqmomCase adaptive_cqmom_case_from_moments_with_hierarchy_and_storage(const qmomkit::RawMomentTensor& moments,
                                                                               const std::vector<std::size_t>& node_counts,
                                                                               const std::vector<std::size_t>& conditioning_order,
                                                                               const char* const storage)
{
    std::string text = cqmom_case_text_from_moments_with_hierarchy(moments, node_counts, conditioning_order, true, true);
    text += "\n[cqmom]\nmoment_storage = ";
    text += storage;
    text += "\n";
    const auto parsed = qmomkit::config::parse_cqmom(text, "adaptive-cqmom-case-from-moments-with-hierarchy-and-storage");
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

bool transitions_contain(const std::vector<qmomkit::CqmomNodeCountTransition>& transitions, const qmomkit::TransportStage stage,
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

void assert_transition_diagnostics_equal(const std::vector<qmomkit::CqmomNodeCountTransition>& lhs,
                                         const std::vector<qmomkit::CqmomNodeCountTransition>& rhs)
{
    assert(lhs.size() == rhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        assert(lhs[i].stage == rhs[i].stage);
        assert(lhs[i].conditioning_level == rhs[i].conditioning_level);
        assert(lhs[i].parent_path == rhs[i].parent_path);
        assert(lhs[i].phase_space_dimension == rhs[i].phase_space_dimension);
        assert(lhs[i].configured_maximum_node_count == rhs[i].configured_maximum_node_count);
        assert(lhs[i].active_node_count_before == rhs[i].active_node_count_before);
        assert(lhs[i].active_node_count_after == rhs[i].active_node_count_after);
        assert(lhs[i].kind == rhs[i].kind);
        assert(lhs[i].rank_reduced == rhs[i].rank_reduced);
        assert(lhs[i].rank_grew == rhs[i].rank_grew);
    }
}

void assert_common_state_moments_near(const qmomkit::CqmomState& dense_state, const qmomkit::CqmomState& lean_state, const double tolerance)
{
    assert(dense_state.moment_storage() == qmomkit::CqmomMomentStorage::dense);
    assert(lean_state.moment_storage() == qmomkit::CqmomMomentStorage::lean);
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(dense_state.hierarchy());
    const qmomkit::RawMomentTensor dense_moments{dense_state.hierarchy().moment_extents(), dense_state.moment_values()};
    const qmomkit::CqmomMomentData dense_required = qmomkit::extract_cqmom_moments(dense_moments, layout);
    assert_vector_near(dense_required.values, lean_state.moment_values(), tolerance);
}

std::vector<double> common_cqmom_moment_values(const qmomkit::CqmomState& state)
{
    if (state.moment_storage() == qmomkit::CqmomMomentStorage::lean)
    {
        return state.moment_values();
    }
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(state.hierarchy());
    const qmomkit::RawMomentTensor dense_moments{state.hierarchy().moment_extents(), state.moment_values()};
    return qmomkit::extract_cqmom_moments(dense_moments, layout).values;
}

void assert_common_cqmom_moments_near(const qmomkit::CqmomState& lhs, const qmomkit::CqmomState& rhs, const double tolerance)
{
    assert(lhs.hierarchy().node_counts() == rhs.hierarchy().node_counts());
    assert(lhs.hierarchy().conditioning_order() == rhs.hierarchy().conditioning_order());
    assert_vector_near(common_cqmom_moment_values(lhs), common_cqmom_moment_values(rhs), tolerance);
}

void assert_cqmom_state_moments_match_quadrature(const qmomkit::CqmomState& state, const qmomkit::CqmomQuadrature& quadrature, const double tolerance)
{
    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(quadrature, state.hierarchy().moment_extents());
    if (state.moment_storage() == qmomkit::CqmomMomentStorage::dense)
    {
        assert_vector_near(state.moment_values(), reconstructed.values, tolerance);
        return;
    }

    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(state.hierarchy());
    const qmomkit::CqmomMomentData reconstructed_required = qmomkit::extract_cqmom_moments(reconstructed, layout);
    assert_vector_near(state.moment_values(), reconstructed_required.values, tolerance);
}

void assert_cqmom_initial_state_stores_accepted_quadrature(const qmomkit::CqmomCase& case_, const double tolerance)
{
    const qmomkit::CqmomResult expected = qmomkit::compute_cqmom(case_);
    assert(expected.ok());

    const qmomkit::CqmomState state = case_.initial_state();
    const std::optional<qmomkit::CqmomQuadrature>& accepted = qmomkit::internal::cqmom_state_accepted_quadrature(state);

    assert(accepted.has_value());
    assert_cqmom_quadrature_near(*accepted, expected.quadrature, tolerance);
    assert(active_counts_equal(state.active_branch_counts(), qmomkit::internal::active_branch_counts_from_quadrature(*accepted)));
    assert_cqmom_state_moments_match_quadrature(state, *accepted, tolerance);
}

bool layout_contains_index(const qmomkit::CqmomMomentLayout& layout, const qmomkit::MultiIndex& index)
{
    return std::find(layout.required_indices.begin(), layout.required_indices.end(), index) != layout.required_indices.end();
}

void assert_dense_state_matches_reconstructed_quadrature(const qmomkit::CqmomTransportResult& result, const std::vector<std::size_t>& extents,
                                                         const double tolerance)
{
    assert(result.ok());
    assert(result.state.moment_storage() == qmomkit::CqmomMomentStorage::dense);
    assert(result.quadrature.has_value());
    const qmomkit::RawMomentTensor stored{extents, result.state.moment_values()};
    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(*result.quadrature, extents);
    assert_vector_near(stored.values, reconstructed.values, tolerance);
}

bool dense_only_moment_changed(const qmomkit::RawMomentTensor& before, const qmomkit::CqmomState& after,
                               const qmomkit::PhaseSpaceHierarchy& hierarchy, const double tolerance)
{
    assert(after.moment_storage() == qmomkit::CqmomMomentStorage::dense);
    assert(before.values.size() == after.moment_values().size());
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    for (std::size_t flat = 0; flat < before.values.size(); ++flat)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat, before.extents);
        if (!layout_contains_index(layout, index) && std::abs(before.values[flat] - after.moment_values()[flat]) > tolerance)
        {
            return true;
        }
    }
    return false;
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

qmomkit::CqmomMomentSource dense_source_from_to(const qmomkit::RawMomentTensor& initial, const qmomkit::RawMomentTensor& target)
{
    return [initial, target](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            const std::size_t flat = qmomkit::flatten_multi_index(moments.index(i), target.extents);
            rate[i] = target.values[flat] - initial.values[flat];
        }
        return qmomkit::CqmomMomentRate{rate};
    };
}

struct SourceContextSnapshot
{
    bool captured = false;
    double time = 0.0;
    std::size_t dimension_count = 0;
    std::vector<std::size_t> node_counts;
    std::vector<std::size_t> conditioning_order;
    std::vector<qmomkit::MultiIndex> indices;
    std::vector<double> values;
    std::vector<double> lookup_values;
    std::optional<double> m0;
};

void capture_source_context_once(SourceContextSnapshot& snapshot, const qmomkit::CqmomSourceContext& context)
{
    if (snapshot.captured)
    {
        return;
    }
    snapshot.captured = true;
    snapshot.time = context.time();
    snapshot.dimension_count = context.dimension_count();
    snapshot.node_counts = context.hierarchy().node_counts();
    snapshot.conditioning_order = context.hierarchy().conditioning_order();
    snapshot.indices = context.indices();
    snapshot.values = context.values();
    snapshot.lookup_values.reserve(context.size());
    for (const qmomkit::MultiIndex& index : context.indices())
    {
        const std::optional<double> value = context.value(index);
        assert(value.has_value());
        snapshot.lookup_values.push_back(*value);
    }
    snapshot.m0 = context.value(qmomkit::MultiIndex(context.dimension_count(), 0));
}

void assert_source_context_snapshots_near(const SourceContextSnapshot& dense, const SourceContextSnapshot& lean, const double tolerance)
{
    assert(dense.captured);
    assert(lean.captured);
    assert(dense.time == lean.time);
    assert(dense.dimension_count == lean.dimension_count);
    assert(dense.node_counts == lean.node_counts);
    assert(dense.conditioning_order == lean.conditioning_order);
    assert(dense.indices == lean.indices);
    assert(dense.m0.has_value());
    assert(lean.m0.has_value());
    assert_vector_near(dense.values, lean.values, tolerance);
    assert_vector_near(dense.lookup_values, lean.lookup_values, tolerance);
    assert_near(*dense.m0, *lean.m0, tolerance);
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

void compare_lean_transport_to_dense_layout_closed_source(const qmomkit::TransportScheme scheme)
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestLeanCqmomTransportState state = lean_state_from_dense(moments, hierarchy);
    const auto velocity = [](double, const std::vector<double>& xi)
    {
        assert(xi.size() == 2);
        return std::vector<double>{0.2 * xi[0] + 0.1, -0.15 * xi[1] + 0.05};
    };
    constexpr double lambda = -0.35;
    const qmomkit::CqmomMomentSource source = [lambda](double, const qmomkit::CqmomSourceContext& state_values)
    {
        std::vector<double> rate(state_values.size(), 0.0);
        for (std::size_t i = 0; i < state_values.size(); ++i)
        {
            rate[i] = lambda * state_values.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };

    qmomkit::TestDenseCqmomTransportOptions dense_options;
    dense_options.node_options = tight_node_options();
    dense_options.scheme = scheme;
    const qmomkit::TestDenseCqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.2, velocity, source, dense_options);

    qmomkit::TestLeanCqmomTransportOptions lean_options;
    lean_options.node_options = tight_node_options();
    lean_options.scheme = scheme;
    const qmomkit::TestLeanCqmomTransportResult lean = qmomkit::advance_transport_cqmom(state, 0.0, 0.2, velocity, source, lean_options);

    assert(dense.ok());
    assert(lean.ok());
    assert(lean.quadrature.has_value());
    assert(lean.quadrature->has_rectangular_storage());
    assert(lean.diagnostics.scheme == scheme);

    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments(dense.moments, state.moments.layout);
    assert_vector_near(lean.state.moments.values, expected.values, 5.0e-9);
}

void test_explicit_dense_and_lean_adaptive_branch_options_are_equivalent()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = adaptive_cqmom_moments();
    const qmomkit::TestLeanCqmomTransportState lean_state = lean_state_from_dense(moments, hierarchy);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& values)
    { return qmomkit::CqmomMomentRate{std::vector<double>(values.size(), 0.0)}; };

    qmomkit::TestDenseCqmomTransportOptions dense_options;
    dense_options.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    dense_options.adaptive_branch_counts = true;
    dense_options.diagnostic_level = qmomkit::CqmomDiagnosticLevel::detailed;
    dense_options.vandermonde_solver = qmomkit::VandermondeSolver::dense;
    dense_options.scheme = qmomkit::TransportScheme::dc;
    dense_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    qmomkit::TestLeanCqmomTransportOptions lean_options;
    lean_options.qmom_options = dense_options.qmom_options;
    lean_options.adaptive_branch_counts = dense_options.adaptive_branch_counts;
    lean_options.diagnostic_level = dense_options.diagnostic_level;
    lean_options.vandermonde_solver = dense_options.vandermonde_solver;
    lean_options.scheme = dense_options.scheme;
    lean_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::TestDenseCqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.1, zero_velocity, zero_source, dense_options);
    const qmomkit::TestLeanCqmomTransportResult lean =
        qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.1, zero_velocity, zero_source, lean_options);

    assert(dense.ok());
    assert(lean.ok());
    assert(dense.quadrature.has_value());
    assert(lean.quadrature.has_value());
    assert(dense.quadrature->adaptive_branch_counts_enabled());
    assert(lean.quadrature->adaptive_branch_counts_enabled());
    assert(!dense.quadrature->has_rectangular_storage());
    assert(!lean.quadrature->has_rectangular_storage());
    assert_cqmom_quadrature_near(*lean.quadrature, *dense.quadrature, 1.0e-12);
    assert(dense.diagnostics.cqmom_diagnostic_level == qmomkit::CqmomDiagnosticLevel::detailed);
    assert(lean.diagnostics.cqmom_diagnostic_level == qmomkit::CqmomDiagnosticLevel::detailed);
    assert(dense.diagnostics.cqmom_vandermonde_solver == qmomkit::VandermondeSolver::dense);
    assert(lean.diagnostics.cqmom_vandermonde_solver == qmomkit::VandermondeSolver::dense);
    assert(
        diagnostics_request_solver(dense.diagnostics.initial_cqmom_diagnostics, hierarchy.number_of_dimensions(), qmomkit::VandermondeSolver::dense));
    assert(
        diagnostics_request_solver(lean.diagnostics.initial_cqmom_diagnostics, hierarchy.number_of_dimensions(), qmomkit::VandermondeSolver::dense));
    assert_vector_near(lean.state.moments.values, qmomkit::extract_cqmom_moments(dense.moments, lean_state.moments.layout).values, 1.0e-12);
}

void test_explicit_dense_and_lean_diagnostic_level_none_is_propagated()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestLeanCqmomTransportState lean_state = lean_state_from_dense(moments, hierarchy);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& values)
    { return qmomkit::CqmomMomentRate{std::vector<double>(values.size(), 0.0)}; };

    qmomkit::TestDenseCqmomTransportOptions dense_options;
    dense_options.diagnostic_level = qmomkit::CqmomDiagnosticLevel::none;
    dense_options.scheme = qmomkit::TransportScheme::dc;
    dense_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    qmomkit::TestLeanCqmomTransportOptions lean_options;
    lean_options.diagnostic_level = dense_options.diagnostic_level;
    lean_options.scheme = dense_options.scheme;
    lean_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::TestDenseCqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.1, zero_velocity, zero_source, dense_options);
    const qmomkit::TestLeanCqmomTransportResult lean =
        qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.1, zero_velocity, zero_source, lean_options);

    assert(dense.ok());
    assert(lean.ok());
    assert(dense.diagnostics.cqmom_diagnostic_level == qmomkit::CqmomDiagnosticLevel::none);
    assert(lean.diagnostics.cqmom_diagnostic_level == qmomkit::CqmomDiagnosticLevel::none);
    assert(dense.diagnostics.initial_cqmom_diagnostics.empty());
    assert(lean.diagnostics.initial_cqmom_diagnostics.empty());
    assert(dense.diagnostics.post_source_cqmom_diagnostics.empty());
    assert(lean.diagnostics.post_source_cqmom_diagnostics.empty());
}

void test_explicit_dense_and_lean_adaptive_options_are_equivalent()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestLeanCqmomTransportState lean_state = lean_state_from_dense(moments, hierarchy);
    const auto velocity = [](double, const std::vector<double>& xi)
    {
        assert(xi.size() == 2);
        return std::vector<double>{0.25, -0.1};
    };
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& values)
    { return qmomkit::CqmomMomentRate{std::vector<double>(values.size(), 0.0)}; };

    qmomkit::TestDenseCqmomTransportOptions dense_options;
    dense_options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    dense_options.node_options.adaptive_options.initial_step = 0.01;
    dense_options.node_options.adaptive_options.max_step = 0.01;
    dense_options.scheme = qmomkit::TransportScheme::cd;
    qmomkit::TestLeanCqmomTransportOptions lean_options;
    lean_options.node_options = dense_options.node_options;
    lean_options.scheme = dense_options.scheme;

    const qmomkit::TestDenseCqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.1, velocity, zero_source, dense_options);
    const qmomkit::TestLeanCqmomTransportResult lean = qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.1, velocity, zero_source, lean_options);

    assert(dense.ok());
    assert(lean.ok());
    assert(dense.quadrature.has_value());
    assert(lean.quadrature.has_value());
    assert(dense.diagnostics.first_continuous_diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(lean.diagnostics.first_continuous_diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(dense.diagnostics.first_continuous_diagnostics.step_count > 1);
    assert(lean.diagnostics.first_continuous_diagnostics.step_count == dense.diagnostics.first_continuous_diagnostics.step_count);
    assert_cqmom_quadrature_near(*lean.quadrature, *dense.quadrature, 1.0e-11);
}

void test_explicit_dense_and_lean_source_options_are_equivalent()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestLeanCqmomTransportState lean_state = lean_state_from_dense(moments, hierarchy);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& values)
    { return qmomkit::CqmomMomentRate{std::vector<double>(values.size(), 0.0)}; };

    qmomkit::TestDenseCqmomTransportOptions dense_options;
    dense_options.scheme = qmomkit::TransportScheme::dc;
    dense_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    dense_options.source_options.ssp_rk3_options.max_retries = 2;
    dense_options.source_options.ssp_rk3_options.accept = [](const std::vector<double>&) { return false; };
    qmomkit::TestLeanCqmomTransportOptions lean_options;
    lean_options.scheme = dense_options.scheme;
    lean_options.node_options = dense_options.node_options;
    lean_options.source_options = dense_options.source_options;

    const qmomkit::TestDenseCqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.1, zero_velocity, zero_source, dense_options);
    const qmomkit::TestLeanCqmomTransportResult lean =
        qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.1, zero_velocity, zero_source, lean_options);

    assert(!dense.ok());
    assert(!lean.ok());
    assert(dense.status == qmomkit::TransportStatus::source_failed);
    assert(lean.status == qmomkit::TransportStatus::source_failed);
    assert(dense.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
    assert(lean.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::rejected_by_predicate);
    assert(dense.diagnostics.source_retries == 2);
    assert(lean.diagnostics.source_retries == 2);
}

void test_explicit_dense_and_lean_accept_same_cqmom_native_source()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestLeanCqmomTransportState lean_state = lean_state_from_dense(moments, hierarchy);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    std::size_t inspected_context_count = 0;

    const qmomkit::CqmomMomentSource source = [&](const double time, const qmomkit::CqmomSourceContext& context)
    {
        assert_near(time, context.time(), 1.0e-15);
        assert(context.dimension_count() == hierarchy.number_of_dimensions());
        if (context.indices() == lean_state.moments.layout.required_indices)
        {
            ++inspected_context_count;
            assert(context.moment_count() == lean_state.moments.values.size());
            assert(context.indices().size() == context.values().size());
            for (std::size_t entry = 0; entry < context.size(); ++entry)
            {
                const qmomkit::MultiIndex& index = context.index(entry);
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

    qmomkit::TestDenseCqmomTransportOptions dense_options;
    dense_options.scheme = qmomkit::TransportScheme::dc;
    dense_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    qmomkit::TestLeanCqmomTransportOptions lean_options;
    lean_options.scheme = dense_options.scheme;
    lean_options.node_options = dense_options.node_options;

    const qmomkit::TestDenseCqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.25, 0.02, zero_velocity, source, dense_options);
    const qmomkit::TestLeanCqmomTransportResult lean = qmomkit::advance_transport_cqmom(lean_state, 0.25, 0.02, zero_velocity, source, lean_options);

    assert(dense.ok());
    assert(lean.ok());
    assert(inspected_context_count > 0);
    assert(dense.diagnostics.source_status == lean.diagnostics.source_status);
    assert_vector_near(lean.state.moments.values, qmomkit::extract_cqmom_moments(dense.moments, lean_state.moments.layout).values, 1.0e-11);
}

void test_explicit_dense_and_lean_accept_same_cqmom_native_boundary()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestLeanCqmomTransportState lean_state = lean_state_from_dense(moments, hierarchy);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomMomentSource source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };
    const qmomkit::CqmomMomentBoundaryContribution boundary = [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            rate[entry] = 0.02 * context.value(entry);
        }
        return qmomkit::CqmomMomentRate{rate};
    };

    qmomkit::TestDenseCqmomTransportOptions dense_options;
    dense_options.scheme = qmomkit::TransportScheme::dc;
    dense_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    qmomkit::TestLeanCqmomTransportOptions lean_options;
    lean_options.scheme = dense_options.scheme;
    lean_options.node_options = dense_options.node_options;

    const qmomkit::TestDenseCqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.01, zero_velocity, source, boundary, dense_options);
    const qmomkit::TestLeanCqmomTransportResult lean =
        qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.01, zero_velocity, source, boundary, lean_options);

    assert(dense.ok());
    assert(lean.ok());
    assert(dense.diagnostics.source_status == lean.diagnostics.source_status);
    assert_vector_near(lean.state.moments.values, qmomkit::extract_cqmom_moments(dense.moments, lean_state.moments.layout).values, 1.0e-11);
}

void test_explicit_dense_and_lean_cqmom_native_source_failures_are_equivalent()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestLeanCqmomTransportState lean_state = lean_state_from_dense(moments, hierarchy);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomMomentSource bad_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size() + 1, 0.0)}; };

    qmomkit::TestDenseCqmomTransportOptions dense_options;
    dense_options.scheme = qmomkit::TransportScheme::dc;
    dense_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    qmomkit::TestLeanCqmomTransportOptions lean_options;
    lean_options.scheme = dense_options.scheme;
    lean_options.node_options = dense_options.node_options;

    const qmomkit::TestDenseCqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.1, zero_velocity, bad_source, dense_options);
    const qmomkit::TestLeanCqmomTransportResult lean =
        qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.1, zero_velocity, bad_source, lean_options);

    assert(!dense.ok());
    assert(!lean.ok());
    assert(dense.status == qmomkit::TransportStatus::source_failed);
    assert(lean.status == qmomkit::TransportStatus::source_failed);
    assert(dense.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(lean.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(dense.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(lean.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(dense.message == lean.message);
}

void test_explicit_dense_and_lean_cqmom_native_source_invalid_values_are_equivalent()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestLeanCqmomTransportState lean_state = lean_state_from_dense(moments, hierarchy);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    qmomkit::TestDenseCqmomTransportOptions dense_options;
    dense_options.scheme = qmomkit::TransportScheme::dc;
    dense_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    qmomkit::TestLeanCqmomTransportOptions lean_options;
    lean_options.scheme = dense_options.scheme;
    lean_options.node_options = dense_options.node_options;

    const qmomkit::CqmomMomentSource nonfinite_source = [](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        rate.front() = std::numeric_limits<double>::infinity();
        return qmomkit::CqmomMomentRate{rate};
    };
    const qmomkit::TestDenseCqmomTransportResult dense_nonfinite =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.1, zero_velocity, nonfinite_source, dense_options);
    const qmomkit::TestLeanCqmomTransportResult lean_nonfinite =
        qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.1, zero_velocity, nonfinite_source, lean_options);
    assert(!dense_nonfinite.ok());
    assert(!lean_nonfinite.ok());
    assert(dense_nonfinite.status == qmomkit::TransportStatus::source_failed);
    assert(lean_nonfinite.status == qmomkit::TransportStatus::source_failed);
    assert(dense_nonfinite.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(lean_nonfinite.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(dense_nonfinite.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(lean_nonfinite.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);

    const qmomkit::CqmomMomentSource throwing_source = [](double, const qmomkit::CqmomSourceContext&) -> qmomkit::CqmomMomentRate
    { throw std::runtime_error("native CQMOM source failure"); };
    const qmomkit::TestDenseCqmomTransportResult dense_throw =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.1, zero_velocity, throwing_source, dense_options);
    const qmomkit::TestLeanCqmomTransportResult lean_throw =
        qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.1, zero_velocity, throwing_source, lean_options);
    assert(!dense_throw.ok());
    assert(!lean_throw.ok());
    assert(dense_throw.status == qmomkit::TransportStatus::source_failed);
    assert(lean_throw.status == qmomkit::TransportStatus::source_failed);
    assert(dense_throw.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(lean_throw.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(dense_throw.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(lean_throw.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(dense_throw.message.find("native CQMOM source failure") != std::string::npos);
    assert(lean_throw.message.find("native CQMOM source failure") != std::string::npos);
}

void test_explicit_dense_and_lean_source_growth_is_automatic()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = collapsed_branch_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const qmomkit::TestLeanCqmomTransportState lean_state{initial_required};
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    qmomkit::TestDenseCqmomTransportOptions dense_growth_options;
    dense_growth_options.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    dense_growth_options.adaptive_branch_counts = true;
    dense_growth_options.scheme = qmomkit::TransportScheme::dc;
    dense_growth_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    qmomkit::TestLeanCqmomTransportOptions lean_growth_options;
    lean_growth_options.qmom_options = dense_growth_options.qmom_options;
    lean_growth_options.adaptive_branch_counts = dense_growth_options.adaptive_branch_counts;
    lean_growth_options.scheme = dense_growth_options.scheme;
    lean_growth_options.node_options = dense_growth_options.node_options;

    const qmomkit::TestDenseCqmomTransportResult dense_growth = qmomkit::advance_transport_cqmom(
        initial_moments, hierarchy, 0.0, 1.0, zero_velocity, dense_source_from_to(initial_moments, target_moments), dense_growth_options);
    const qmomkit::TestLeanCqmomTransportResult lean_growth =
        qmomkit::advance_transport_cqmom(lean_state, 0.0, 1.0, zero_velocity, source_from_to(initial_required, target_required), lean_growth_options);

    assert(dense_growth.ok());
    assert(lean_growth.ok());
    assert(dense_growth.quadrature.has_value());
    assert(lean_growth.quadrature.has_value());
    assert(dense_growth.quadrature->total_nodes() == 3);
    assert(lean_growth.quadrature->total_nodes() == 3);
    assert_cqmom_quadrature_near(*lean_growth.quadrature, *dense_growth.quadrature, 1.0e-12);
}

void test_public_qmom_continuous_stage_advances_canonical_state()
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = true;
    qmom_options.rmin = 1.0e-10;
    qmom_options.eabs = 1.0e-12;
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.75, 1.25}, {0.4, 0.6}, 6), 3, qmom_options);
    const auto velocity = [](double, const std::vector<double>& xi) { return std::vector<double>{0.2 * xi[0] + 0.1}; };
    qmomkit::QmomTransportOptions options;
    options.qmom_options = qmom_options;
    options.node_options = tight_node_options();

    const qmomkit::QmomTransportResult stage = qmomkit::advance_qmom_continuous(state, 0.0, 0.2, velocity, options);
    const qmomkit::QmomNodeAdvectionResult direct = qmomkit::advance_nodes(state.quadrature(), 0.0, 0.2, velocity, options.node_options);

    assert(stage.ok());
    assert(direct.ok());
    assert(stage.final_time == 0.2);
    assert(stage.state.active_node_count() == state.active_node_count());
    assert(stage.state.active_node_count() == 2);
    assert_qmom_state_canonical(stage.state, 5.0e-10);
    assert_qmom_quadrature_near(stage.state.quadrature(), qmomkit::QmomResult{direct.quadrature.abscissas, direct.quadrature.weights}, 5.0e-10);
    assert_vector_near(stage.state.moments(), qmomkit::reconstruct_raw_moments(direct.quadrature, state.moments().size()), 5.0e-10);
    assert(stage.diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::success);
    assert(stage.diagnostics.first_reconstruction_diagnostics.success);
}

void test_public_qmom_continuous_stage_failure_rolls_back()
{
    const qmomkit::QmomState state = qmom_state_from_moments({1.0, 2.0}, 1);
    const auto bad_velocity = [](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0}; };

    const qmomkit::QmomTransportResult result = qmomkit::advance_qmom_continuous(state, 0.0, 0.1, bad_velocity);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::continuous_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::first_continuous_half_step);
    assert_vector_near(result.state.moments(), state.moments(), 0.0);
    assert_qmom_quadrature_near(result.state.quadrature(), state.quadrature(), 0.0);
}

void test_public_qmom_discontinuous_stage_fixed_mode_is_strict()
{
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.0}, {0.4, 0.6}, 4), 2);
    const std::vector<double> one_node_moments = raw_moments_1d({0.0}, {1.0}, 4);
    qmomkit::QmomTransportOptions options;
    const qmomkit::QmomTransportResult result =
        qmomkit::advance_qmom_discontinuous(state, 2, 0.0, 1.0, qmom_source_from_to(state.moments(), one_node_moments), options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::post_source_inversion);
    assert(result.diagnostics.post_source_qmom_status != qmomkit::QmomStatus::success);
    assert_vector_near(result.state.moments(), state.moments(), 0.0);
    assert_qmom_quadrature_near(result.state.quadrature(), state.quadrature(), 0.0);
}

void test_public_qmom_discontinuous_stage_adaptive_policy()
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = true;
    qmom_options.rmin = 1.0e-10;
    qmom_options.eabs = 1.0e-12;
    qmomkit::QmomTransportOptions options;
    options.qmom_options = qmom_options;

    const qmomkit::QmomState two_node_state = qmom_state_from_moments(raw_moments_1d({-0.75, 1.25}, {0.4, 0.6}, 6), 3, qmom_options);
    const qmomkit::QmomTransportResult retained = qmomkit::advance_qmom_discontinuous(
        two_node_state, 3, 0.0, 0.2, [](double, const std::vector<double>& moments) { return std::vector<double>(moments.size(), 0.0); }, options);
    assert(retained.ok());
    assert(retained.state.active_node_count() == two_node_state.active_node_count());
    assert_qmom_state_canonical(retained.state, 1.0e-12);

    const std::vector<double> one_node_moments = raw_moments_1d({0.5}, {1.0}, 6);
    const qmomkit::QmomTransportResult reduction =
        qmomkit::advance_qmom_discontinuous(two_node_state, 3, 0.0, 1.0, qmom_source_from_to(two_node_state.moments(), one_node_moments), options);
    assert(reduction.ok());
    assert(reduction.state.active_node_count() == 1);
    assert_qmom_state_canonical(reduction.state, 1.0e-12);

    const qmomkit::QmomState one_node_state = qmom_state_from_moments(one_node_moments, 3, qmom_options);
    const std::vector<double> two_node_moments = raw_moments_1d({-0.25, 1.5}, {0.35, 0.65}, 6);
    const qmomkit::QmomTransportResult growth =
        qmomkit::advance_qmom_discontinuous(one_node_state, 3, 0.0, 1.0, qmom_source_from_to(one_node_state.moments(), two_node_moments), options);
    assert(growth.ok());
    assert(growth.state.active_node_count() == 2);
    assert_qmom_state_canonical(growth.state, 1.0e-12);

    const std::vector<double> three_node_moments = raw_moments_1d({-1.0, 0.25, 1.5}, {0.2, 0.3, 0.5}, 6);
    const qmomkit::QmomTransportResult capped_growth =
        qmomkit::advance_qmom_discontinuous(one_node_state, 3, 0.0, 1.0, qmom_source_from_to(one_node_state.moments(), three_node_moments), options);
    assert(capped_growth.ok());
    assert(capped_growth.state.active_node_count() == 2);
    assert_qmom_state_canonical(capped_growth.state, 1.0e-12);
}

void test_public_qmom_discontinuous_stage_failures_roll_back()
{
    const qmomkit::QmomState state = qmom_state_from_moments({1.0, 2.0}, 1);
    const auto bad_source = [](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0, 3.0}; };
    const qmomkit::QmomTransportResult source_failure = qmomkit::advance_qmom_discontinuous(state, 1, 0.0, 0.1, bad_source);
    assert(!source_failure.ok());
    assert(source_failure.status == qmomkit::TransportStatus::source_failed);
    assert(source_failure.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert_vector_near(source_failure.state.moments(), state.moments(), 0.0);

    const auto nonrealizable_source = [](double, const std::vector<double>& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        rate.front() = -20.0;
        return rate;
    };
    const qmomkit::QmomTransportResult inversion_failure = qmomkit::advance_qmom_discontinuous(state, 1, 0.0, 0.1, nonrealizable_source);
    assert(!inversion_failure.ok());
    assert(inversion_failure.status == qmomkit::TransportStatus::inversion_failed);
    assert(inversion_failure.diagnostics.failed_stage == qmomkit::TransportStage::post_source_inversion);
    assert_vector_near(inversion_failure.state.moments(), state.moments(), 0.0);
    assert_qmom_quadrature_near(inversion_failure.state.quadrature(), state.quadrature(), 0.0);
}

void test_public_qmom_discontinuous_stage_boundary_matches_source_update()
{
    const qmomkit::QmomState state = qmom_state_from_moments({1.0, 2.0}, 1);
    const qmomkit::MomentSource source = [](double, const std::vector<double>& moments)
    {
        assert(moments.size() == 2);
        return std::vector<double>{0.0, 1.0};
    };
    const qmomkit::MomentBoundaryContribution boundary = [](double, const std::vector<double>& moments)
    {
        assert(moments.size() == 2);
        return qmomkit::MomentBoundaryRate{std::vector<double>{0.25, -0.5}};
    };

    const qmomkit::QmomTransportResult stage = qmomkit::advance_qmom_discontinuous(state, 1, 0.125, 0.25, source, boundary, {});
    const qmomkit::MomentSourceUpdateResult direct = qmomkit::advance_moments_with_source(state.moments(), 0.125, 0.25, source, boundary, {});

    assert(stage.ok());
    assert(direct.ok());
    assert_vector_near(stage.state.moments(), direct.moments, 1.0e-13);
    assert_qmom_state_canonical(stage.state, 1.0e-13);
}

void assert_qmom_transport_results_match(const qmomkit::QmomTransportResult& actual, const qmomkit::QmomTransportResult& expected,
                                         const double tolerance)
{
    assert(actual.ok());
    assert(expected.ok());
    assert_near(actual.final_time, expected.final_time, tolerance);
    assert(actual.state.active_node_count() == expected.state.active_node_count());
    assert_vector_near(actual.state.moments(), expected.state.moments(), tolerance);
    assert_qmom_quadrature_near(actual.state.quadrature(), expected.state.quadrature(), tolerance);
}

void test_public_qmom_stage_composition_matches_strang_transport()
{
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4), 2);
    const auto velocity = [](const double time, const std::vector<double>& xi) { return std::vector<double>{0.35 + 0.1 * time - 0.05 * xi[0]}; };
    const auto source = [](const double time, const std::vector<double>& moments)
    {
        std::vector<double> rhs(moments.size(), 0.0);
        for (std::size_t order = 1; order < moments.size(); ++order)
        {
            rhs[order] = -(0.2 + 0.05 * time) * static_cast<double>(order) * moments[order];
        }
        return rhs;
    };
    qmomkit::QmomTransportOptions options;
    options.node_options = tight_node_options();
    options.scheme = qmomkit::TransportScheme::strang;
    constexpr double dt = 0.2;

    const qmomkit::QmomTransportResult high_level = qmomkit::advance_transport_qmom(state, 2, 0.0, dt, velocity, source, options);
    const qmomkit::QmomTransportResult first = qmomkit::advance_qmom_continuous(state, 0.0, 0.5 * dt, velocity, options);
    assert(first.ok());
    const qmomkit::QmomTransportResult discontinuous = qmomkit::advance_qmom_discontinuous(first.state, 2, 0.5 * dt, dt, source, options);
    assert(discontinuous.ok());
    const qmomkit::QmomTransportResult second = qmomkit::advance_qmom_continuous(discontinuous.state, 0.5 * dt, 0.5 * dt, velocity, options);

    assert_qmom_transport_results_match(second, high_level, 5.0e-10);
    assert(high_level.diagnostics.first_continuous_diagnostics.step_count == first.diagnostics.first_continuous_diagnostics.step_count);
    assert(high_level.diagnostics.source_status == discontinuous.diagnostics.source_status);
    assert(high_level.diagnostics.second_continuous_diagnostics.step_count == second.diagnostics.first_continuous_diagnostics.step_count);
}

void test_public_qmom_stage_composition_matches_lie_transport()
{
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4), 2);
    const auto velocity = [](const double time, const std::vector<double>& xi) { return std::vector<double>{0.2 + time - 0.1 * xi[0]}; };
    const auto source = [](const double time, const std::vector<double>& moments)
    {
        std::vector<double> rhs(moments.size(), 0.0);
        for (std::size_t order = 1; order < moments.size(); ++order)
        {
            rhs[order] = -(0.1 + 0.02 * time) * static_cast<double>(order) * moments[order];
        }
        return rhs;
    };
    constexpr double dt = 0.2;

    qmomkit::QmomTransportOptions cd_options;
    cd_options.node_options = tight_node_options();
    cd_options.scheme = qmomkit::TransportScheme::cd;
    const qmomkit::QmomTransportResult cd_high = qmomkit::advance_transport_qmom(state, 2, 0.0, dt, velocity, source, cd_options);
    const qmomkit::QmomTransportResult cd_continuous = qmomkit::advance_qmom_continuous(state, 0.0, dt, velocity, cd_options);
    assert(cd_continuous.ok());
    const qmomkit::QmomTransportResult cd_explicit = qmomkit::advance_qmom_discontinuous(cd_continuous.state, 2, 0.0, dt, source, cd_options);
    assert_qmom_transport_results_match(cd_explicit, cd_high, 5.0e-10);

    qmomkit::QmomTransportOptions dc_options = cd_options;
    dc_options.scheme = qmomkit::TransportScheme::dc;
    const qmomkit::QmomTransportResult dc_high = qmomkit::advance_transport_qmom(state, 2, 0.0, dt, velocity, source, dc_options);
    const qmomkit::QmomTransportResult dc_discontinuous = qmomkit::advance_qmom_discontinuous(state, 2, 0.0, dt, source, dc_options);
    assert(dc_discontinuous.ok());
    const qmomkit::QmomTransportResult dc_explicit = qmomkit::advance_qmom_continuous(dc_discontinuous.state, 0.0, dt, velocity, dc_options);
    assert_qmom_transport_results_match(dc_explicit, dc_high, 5.0e-10);
}

void test_qmom_transport_failure_after_first_stage_rolls_back_to_macro_input()
{
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4), 2);
    const auto velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.5}; };
    const auto bad_source = [](double, const std::vector<double>& moments) { return std::vector<double>(moments.size() + 1, 0.0); };
    qmomkit::QmomTransportOptions options;
    options.node_options = tight_node_options();
    options.scheme = qmomkit::TransportScheme::strang;

    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(state, 2, 0.0, 0.2, velocity, bad_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(result.diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::success);
    assert_vector_near(result.state.moments(), state.moments(), 0.0);
    assert_qmom_quadrature_near(result.state.quadrature(), state.quadrature(), 0.0);
}

void test_cqmom_initial_state_stores_accepted_quadrature_fixed_dense_and_lean()
{
    assert_cqmom_initial_state_stores_accepted_quadrature(fixed_cqmom_case_with_storage("dense"), 5.0e-10);
    assert_cqmom_initial_state_stores_accepted_quadrature(fixed_cqmom_case_with_storage("lean"), 5.0e-10);
}

void test_cqmom_initial_state_stores_accepted_quadrature_adaptive_dense_and_lean()
{
    const qmomkit::RawMomentTensor moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::CqmomCase dense_case = adaptive_cqmom_case_from_moments_with_storage(moments, "dense");
    const qmomkit::CqmomCase lean_case = adaptive_cqmom_case_from_moments_with_storage(moments, "lean");

    assert_cqmom_initial_state_stores_accepted_quadrature(dense_case, 5.0e-10);
    assert_cqmom_initial_state_stores_accepted_quadrature(lean_case, 5.0e-10);

    const qmomkit::CqmomState dense_state = dense_case.initial_state();
    const qmomkit::CqmomState lean_state = lean_case.initial_state();
    const std::optional<qmomkit::CqmomQuadrature>& dense_quadrature = qmomkit::internal::cqmom_state_accepted_quadrature(dense_state);
    const std::optional<qmomkit::CqmomQuadrature>& lean_quadrature = qmomkit::internal::cqmom_state_accepted_quadrature(lean_state);

    assert(dense_quadrature.has_value());
    assert(lean_quadrature.has_value());
    assert(dense_quadrature->terminal_paths() == lean_quadrature->terminal_paths());
    assert(active_counts_equal(dense_state.active_branch_counts(), lean_state.active_branch_counts()));
    assert(active_counts_equal(dense_state.active_branch_counts(), qmomkit::internal::active_branch_counts_from_quadrature(*dense_quadrature)));
}

void test_cqmom_first_continuous_stage_uses_initial_accepted_quadrature()
{
    const qmomkit::CqmomState state = adaptive_cqmom_case_from_moments(branch_specific_ragged_cqmom_moments()).initial_state();
    const std::optional<qmomkit::CqmomQuadrature>& initial_quadrature = qmomkit::internal::cqmom_state_accepted_quadrature(state);
    assert(initial_quadrature.has_value());
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    qmomkit::CqmomTransportOptions options;
    options.node_options = tight_node_options();
    const qmomkit::CqmomTransportResult stage = qmomkit::advance_cqmom_continuous(state, 0.0, 0.1, zero_velocity, options);

    assert(stage.ok());
    assert(stage.quadrature.has_value());
    assert(stage.diagnostics.initial_cqmom_diagnostics.empty());
    assert(stage.quadrature->terminal_paths() == initial_quadrature->terminal_paths());
    assert_cqmom_quadrature_near(*stage.quadrature, *initial_quadrature, 5.0e-10);
    assert(active_counts_equal(stage.state.active_branch_counts(), state.active_branch_counts()));
}

void assert_cqmom_transport_results_match(const qmomkit::CqmomTransportResult& actual, const qmomkit::CqmomTransportResult& expected,
                                          const double tolerance)
{
    assert(actual.ok());
    assert(expected.ok());
    assert_near(actual.final_time, expected.final_time, tolerance);
    assert(actual.quadrature.has_value());
    assert(expected.quadrature.has_value());
    assert_cqmom_quadrature_near(*actual.quadrature, *expected.quadrature, tolerance);
    assert(active_counts_equal(actual.state.active_branch_counts(), expected.state.active_branch_counts()));
    assert_common_cqmom_moments_near(actual.state, expected.state, tolerance);
}

void test_public_cqmom_continuous_stage_preserves_topology()
{
    const qmomkit::CqmomCase case_ = adaptive_cqmom_case_from_moments(branch_specific_ragged_cqmom_moments());
    const qmomkit::CqmomState state = case_.initial_state();
    const qmomkit::CqmomResult initial = qmomkit::compute_cqmom(case_);
    assert(initial.ok());
    const auto velocity = [](double, const std::vector<double>& xi)
    {
        assert(xi.size() == 2);
        return std::vector<double>{0.1 + 0.05 * xi[0], -0.08 * xi[1]};
    };

    qmomkit::CqmomTransportOptions options;
    options.node_options = tight_node_options();
    const qmomkit::CqmomTransportResult stage = qmomkit::advance_cqmom_continuous(state, 0.0, 0.2, velocity, options);

    assert(stage.ok());
    assert(stage.quadrature.has_value());
    assert(stage.final_time == 0.2);
    assert(stage.quadrature->hierarchy().node_counts() == initial.quadrature.hierarchy().node_counts());
    assert(stage.quadrature->hierarchy().conditioning_order() == initial.quadrature.hierarchy().conditioning_order());
    assert(stage.quadrature->terminal_paths() == initial.quadrature.terminal_paths());
    assert_vector_near(stage.quadrature->weights(), initial.quadrature.weights(), 0.0);
    assert(active_counts_equal(stage.state.active_branch_counts(), state.active_branch_counts()));
    assert(stage.diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::success);
    assert(stage.diagnostics.first_reconstruction_diagnostics.success);
}

void test_public_cqmom_continuous_stage_reconstructs_dense_and_lean_moments()
{
    const qmomkit::CqmomCase dense_case = fixed_cqmom_case_with_storage("dense");
    const qmomkit::CqmomCase lean_case = fixed_cqmom_case_with_storage("lean");
    const auto velocity = [](double, const std::vector<double>& xi)
    {
        assert(xi.size() == 2);
        return std::vector<double>{0.15 + 0.02 * xi[0], -0.05 + 0.03 * xi[1]};
    };
    qmomkit::CqmomTransportOptions options;
    options.node_options = tight_node_options();

    const qmomkit::CqmomTransportResult dense = qmomkit::advance_cqmom_continuous(dense_case.initial_state(), 0.0, 0.1, velocity, options);
    const qmomkit::CqmomTransportResult lean = qmomkit::advance_cqmom_continuous(lean_case.initial_state(), 0.0, 0.1, velocity, options);

    assert(dense.ok());
    assert(lean.ok());
    assert(dense.quadrature.has_value());
    assert(lean.quadrature.has_value());
    assert_cqmom_quadrature_near(*dense.quadrature, *lean.quadrature, 5.0e-10);
    assert_common_state_moments_near(dense.state, lean.state, 5.0e-10);
    assert_dense_state_matches_reconstructed_quadrature(dense, dense.state.hierarchy().moment_extents(), 5.0e-10);
    assert_cqmom_state_moments_match_quadrature(lean.state, *lean.quadrature, 5.0e-10);
}

void test_public_cqmom_discontinuous_stage_adaptive_policy()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor retained_moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData retained_required = qmomkit::extract_cqmom_moments(retained_moments, layout);
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomTransportOptions options;
    const qmomkit::CqmomState retained_state = adaptive_cqmom_case_from_moments(retained_moments).initial_state();
    const qmomkit::CqmomTransportResult retained = qmomkit::advance_cqmom_discontinuous(retained_state, 0.0, 0.1, zero_source, options);
    assert(retained.ok());
    assert(retained.quadrature.has_value());
    assert(active_counts_equal(retained.state.active_branch_counts(), retained_state.active_branch_counts()));
    assert_cqmom_state_moments_match_quadrature(retained.state, *retained.quadrature, 5.0e-10);

    const qmomkit::RawMomentTensor reduced_moments = collapsed_branch_cqmom_moments();
    const qmomkit::CqmomMomentData reduced_required = qmomkit::extract_cqmom_moments(reduced_moments, layout);
    const qmomkit::CqmomTransportResult reduced =
        qmomkit::advance_cqmom_discontinuous(retained_state, 0.0, 1.0, source_from_to(retained_required, reduced_required), options);
    assert(reduced.ok());
    assert(reduced.quadrature.has_value());
    assert(reduced.quadrature->total_nodes() == 2);
    assert_cqmom_state_moments_match_quadrature(reduced.state, *reduced.quadrature, 5.0e-10);
    assert(transitions_contain(reduced.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{1}}, 2, 1, qmomkit::CqmomNodeCountTransitionKind::reduced));

    const qmomkit::CqmomState collapsed_state = adaptive_cqmom_case_from_moments(reduced_moments).initial_state();
    const qmomkit::CqmomTransportResult growth =
        qmomkit::advance_cqmom_discontinuous(collapsed_state, 0.0, 1.0, source_from_to(reduced_required, retained_required), options);
    assert(growth.ok());
    assert(growth.quadrature.has_value());
    assert(growth.quadrature->total_nodes() == 3);
    assert_cqmom_state_moments_match_quadrature(growth.state, *growth.quadrature, 5.0e-10);
    assert(transitions_contain(growth.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{1}}, 1, 2, qmomkit::CqmomNodeCountTransitionKind::grew));
}

void test_public_cqmom_discontinuous_stage_growth_is_limited_to_one_node()
{
    constexpr std::size_t configured_node_count = 3;
    const qmomkit::RawMomentTensor initial_moments = one_dimensional_cqmom_moments(configured_node_count, {0.0}, {1.0});
    const qmomkit::RawMomentTensor target_moments = one_dimensional_cqmom_moments(configured_node_count, {-1.0, 0.25, 2.0}, {0.2, 0.3, 0.5});
    const qmomkit::PhaseSpaceHierarchy hierarchy({configured_node_count});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);

    const qmomkit::CqmomTransportResult stage = qmomkit::advance_cqmom_discontinuous(one_dimensional_adaptive_cqmom_state(initial_required), 0.0, 1.0,
                                                                                     source_from_to(initial_required, target_required));

    assert(stage.ok());
    assert(stage.quadrature.has_value());
    assert(stage.quadrature->total_nodes() == 2);
    assert(active_counts_contain(stage.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, configured_node_count));
    assert_cqmom_state_moments_match_quadrature(stage.state, *stage.quadrature, 5.0e-10);
}

void test_lean_cqmom_discontinuous_stage_canonicalizes_adaptive_reduction()
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

    const qmomkit::CqmomTransportResult stage = qmomkit::advance_cqmom_discontinuous(one_dimensional_adaptive_cqmom_state(initial_required), 0.0, 1.0,
                                                                                     source_from_to(initial_required, source_candidate_required));

    assert(stage.ok());
    assert(stage.quadrature.has_value());
    assert(stage.state.moment_storage() == qmomkit::CqmomMomentStorage::lean);
    assert(stage.quadrature->total_nodes() == 1);
    assert(active_counts_contain(stage.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 1, 0, configured_node_count));
    assert(transitions_contain(stage.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 0,
                               qmomkit::ConditioningPath{}, 2, 1, qmomkit::CqmomNodeCountTransitionKind::reduced));
    assert_cqmom_state_moments_match_quadrature(stage.state, *stage.quadrature, 5.0e-12);

    const qmomkit::RawMomentTensor reconstructed = qmomkit::reconstruct_raw_moments(*stage.quadrature, stage.state.hierarchy().moment_extents());
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

void assert_public_cqmom_stage_composition_matches_transport(const qmomkit::TransportScheme scheme)
{
    const qmomkit::CqmomState state = fixed_cqmom_case().initial_state();
    const auto velocity = [](const double time, const std::vector<double>& xi)
    {
        assert(xi.size() == 2);
        return std::vector<double>{0.2 + 0.03 * time - 0.04 * xi[0], -0.1 + 0.02 * time + 0.05 * xi[1]};
    };
    const qmomkit::CqmomMomentSource source = [](const double time, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < context.size(); ++i)
        {
            rate[i] = -(0.18 + 0.04 * time) * static_cast<double>(context.index(i).size()) * context.value(i);
        }
        return qmomkit::CqmomMomentRate{rate};
    };
    qmomkit::CqmomTransportOptions options;
    options.node_options = tight_node_options();
    options.scheme = scheme;
    constexpr double dt = 0.2;

    const qmomkit::CqmomTransportResult high_level = qmomkit::advance_transport_cqmom(state, 0.0, dt, velocity, source, options);
    qmomkit::CqmomTransportResult explicit_result = high_level;
    if (scheme == qmomkit::TransportScheme::strang)
    {
        const qmomkit::CqmomTransportResult first = qmomkit::advance_cqmom_continuous(state, 0.0, 0.5 * dt, velocity, options);
        assert(first.ok());
        const qmomkit::CqmomTransportResult discontinuous = qmomkit::advance_cqmom_discontinuous(first.state, 0.5 * dt, dt, source, options);
        assert(discontinuous.ok());
        explicit_result = qmomkit::advance_cqmom_continuous(discontinuous.state, 0.5 * dt, 0.5 * dt, velocity, options);
    }
    else if (scheme == qmomkit::TransportScheme::cd)
    {
        const qmomkit::CqmomTransportResult continuous = qmomkit::advance_cqmom_continuous(state, 0.0, dt, velocity, options);
        assert(continuous.ok());
        explicit_result = qmomkit::advance_cqmom_discontinuous(continuous.state, 0.0, dt, source, options);
    }
    else
    {
        const qmomkit::CqmomTransportResult discontinuous = qmomkit::advance_cqmom_discontinuous(state, 0.0, dt, source, options);
        assert(discontinuous.ok());
        explicit_result = qmomkit::advance_cqmom_continuous(discontinuous.state, 0.0, dt, velocity, options);
    }

    assert_cqmom_transport_results_match(explicit_result, high_level, 5.0e-9);
}

void test_public_cqmom_stage_composition_matches_transport()
{
    assert_public_cqmom_stage_composition_matches_transport(qmomkit::TransportScheme::strang);
    assert_public_cqmom_stage_composition_matches_transport(qmomkit::TransportScheme::cd);
    assert_public_cqmom_stage_composition_matches_transport(qmomkit::TransportScheme::dc);
}

void test_public_cqmom_stage_failures_roll_back()
{
    const qmomkit::CqmomState state = fixed_cqmom_case().initial_state();
    const auto bad_velocity = [](double, const std::vector<double>&) { return std::vector<double>{1.0}; };
    const qmomkit::CqmomTransportResult continuous = qmomkit::advance_cqmom_continuous(state, 0.0, 0.1, bad_velocity);
    assert(!continuous.ok());
    assert(continuous.status == qmomkit::TransportStatus::continuous_failed);
    assert_vector_near(continuous.state.moment_values(), state.moment_values(), 0.0);
    assert(active_counts_equal(continuous.state.active_branch_counts(), state.active_branch_counts()));

    const auto bad_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size() + 1, 0.0)}; };
    const qmomkit::CqmomTransportResult discontinuous = qmomkit::advance_cqmom_discontinuous(state, 0.0, 0.1, bad_source);
    assert(!discontinuous.ok());
    assert(discontinuous.status == qmomkit::TransportStatus::source_failed);
    assert_vector_near(discontinuous.state.moment_values(), state.moment_values(), 0.0);
    assert(active_counts_equal(discontinuous.state.active_branch_counts(), state.active_branch_counts()));
}

void test_cqmom_transport_failure_after_first_stage_rolls_back_to_macro_input()
{
    const qmomkit::CqmomState state = fixed_cqmom_case().initial_state();
    const auto velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.2); };
    const auto bad_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size() + 1, 0.0)}; };
    qmomkit::CqmomTransportOptions options;
    options.node_options = tight_node_options();
    options.scheme = qmomkit::TransportScheme::strang;

    const qmomkit::CqmomTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, 0.2, velocity, bad_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(result.diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::success);
    assert_vector_near(result.state.moment_values(), state.moment_values(), 0.0);
    assert(active_counts_equal(result.state.active_branch_counts(), state.active_branch_counts()));
}

void test_public_cqmom_boundary_discontinuous_stage_matches_transport()
{
    const qmomkit::CqmomState state = fixed_cqmom_case().initial_state();
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
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
    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::CqmomTransportResult high_level = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, source, boundary, options);
    const qmomkit::CqmomTransportResult stage = qmomkit::advance_cqmom_discontinuous(state, 0.0, 0.1, source, boundary, options);

    assert_cqmom_transport_results_match(stage, high_level, 5.0e-10);
    assert(stage.diagnostics.source_status == high_level.diagnostics.source_status);
}

void test_public_cqmom_stage_irrelevant_options_are_ignored()
{
    const qmomkit::CqmomState state = fixed_cqmom_case().initial_state();
    const auto velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    qmomkit::CqmomTransportOptions continuous_options;
    continuous_options.source_options.ssp_rk3_options.min_time_step = 1.0;
    continuous_options.node_options.adaptive_options.initial_step = 0.1;
    const qmomkit::CqmomTransportResult continuous = qmomkit::advance_cqmom_continuous(state, 0.0, 0.1, velocity, continuous_options);
    assert(continuous.ok());

    qmomkit::CqmomTransportOptions discontinuous_options;
    discontinuous_options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    discontinuous_options.node_options.adaptive_options.max_steps = 0;
    const qmomkit::CqmomTransportResult discontinuous = qmomkit::advance_cqmom_discontinuous(state, 0.0, 0.1, zero_source, discontinuous_options);
    assert(discontinuous.ok());
}

void test_continuous_only_matches_direct_operator_qmom()
{
    const std::vector<double> moments = raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4);
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 2);

    const auto velocity = [](double, const std::vector<double>& xi) { return std::vector<double>{0.5 * xi[0]}; };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };

    qmomkit::QmomTransportOptions options;
    options.node_options = tight_node_options();
    options.scheme = qmomkit::TransportScheme::strang;
    const qmomkit::QmomTransportResult transport = qmomkit::advance_transport_qmom(state, 2, 0.0, 0.3, velocity, zero_source, options);
    const qmomkit::QmomNodeAdvectionResult direct = qmomkit::advance_nodes(state.quadrature(), 0.0, 0.3, velocity, options.node_options);

    assert(transport.ok());
    assert(direct.ok());
    const std::vector<double> expected = qmomkit::reconstruct_raw_moments(direct.quadrature, moments.size());
    assert_vector_near(transport.state.moments(), expected, 5.0e-9);
    assert_qmom_state_canonical(transport.state, 5.0e-9);
}

void test_source_only_matches_direct_operator_qmom()
{
    const std::vector<double> moments{1.0, 2.0};
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto source = [](double, const std::vector<double>& state)
    {
        assert(state.size() == 2);
        return std::vector<double>{0.0, 3.0};
    };

    qmomkit::QmomTransportOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::QmomTransportResult transport = qmomkit::advance_transport_qmom(state, 1, 0.0, 0.25, zero_velocity, source, options);
    const qmomkit::MomentSourceUpdateResult direct = qmomkit::advance_moments_with_source(moments, 0.125, 0.25, source);

    assert(transport.ok());
    assert(direct.ok());
    assert_vector_near(transport.state.moments(), direct.moments, 1.0e-13);
    assert_qmom_state_canonical(transport.state, 1.0e-13);
}

void test_source_plus_boundary_matches_direct_operator_qmom()
{
    const std::vector<double> moments{1.0, 2.0};
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
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

    qmomkit::QmomTransportOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::QmomTransportResult transport = qmomkit::advance_transport_qmom(state, 1, 0.0, 0.25, zero_velocity, source, boundary, options);
    const qmomkit::MomentSourceUpdateResult direct = qmomkit::advance_moments_with_source(moments, 0.125, 0.25, source, boundary, {});

    assert(transport.ok());
    assert(direct.ok());
    assert_vector_near(transport.state.moments(), direct.moments, 1.0e-13);
    assert_qmom_state_canonical(transport.state, 1.0e-13);
}

void test_adaptive_branch_specific_cqmom_case()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = adaptive_cqmom_moments();
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& state)
    { return qmomkit::CqmomMomentRate{std::vector<double>(state.size(), 0.0)}; };

    qmomkit::TestDenseCqmomTransportOptions options;
    options.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    options.adaptive_branch_counts = true;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::TestDenseCqmomTransportResult result =
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
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& state_values)
    { return qmomkit::CqmomMomentRate{std::vector<double>(state_values.size(), 0.0)}; };

    qmomkit::TestDenseCqmomTransportOptions dense_options;
    dense_options.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    dense_options.adaptive_branch_counts = true;
    dense_options.node_options = tight_node_options();
    dense_options.scheme = qmomkit::TransportScheme::cd;
    const qmomkit::TestDenseCqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.2, velocity, zero_source, dense_options);

    qmomkit::CqmomTransportOptions state_options;
    state_options.node_options = dense_options.node_options;
    state_options.scheme = dense_options.scheme;
    const qmomkit::CqmomTransportResult transported = qmomkit::advance_transport_cqmom(state, 0.0, 0.2, velocity, zero_source, state_options);

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

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult transported = qmomkit::advance_transport_cqmom(state, 0.0, 1.0, zero_velocity, source, options);

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

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    qmomkit::CqmomTransportResult first = qmomkit::advance_transport_cqmom(case_.initial_state(), 0.0, 0.1, zero_velocity, zero_source, options);
    assert(first.ok());
    assert(first.quadrature.has_value());
    assert(first.quadrature->total_nodes() == 3);
    assert(active_counts_contain(first.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1));
    assert(active_counts_contain(first.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2));
    assert(first.diagnostics.cqmom_node_count_transitions.empty());
    assert_cqmom_state_moments_match_quadrature(first.state, *first.quadrature, 5.0e-10);

    const qmomkit::CqmomTransportResult second =
        qmomkit::advance_transport_cqmom(first.state, first.final_time, 0.1, zero_velocity, zero_source, options);
    assert(second.ok());
    assert(second.quadrature.has_value());
    assert(second.quadrature->terminal_paths() == first.quadrature->terminal_paths());
    assert(active_counts_contain(second.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1));
    assert(active_counts_contain(second.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2));
    assert_cqmom_state_moments_match_quadrature(second.state, *second.quadrature, 5.0e-10);
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

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult transported =
        qmomkit::advance_transport_cqmom(adaptive_cqmom_case_from_moments(initial_moments).initial_state(), 0.0, 1.0, zero_velocity,
                                         source_from_to(initial_required, target_required), options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->total_nodes() == 2);
    assert(active_counts_contain(transported.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 1));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{1}}, 2, 1, qmomkit::CqmomNodeCountTransitionKind::reduced));
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

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult transported =
        qmomkit::advance_transport_cqmom(adaptive_cqmom_case_from_moments(initial_moments).initial_state(), 0.0, 1.0, zero_velocity,
                                         source_from_to(initial_required, target_required), options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->total_nodes() == 3);
    assert(active_counts_contain(transported.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1));
    assert(active_counts_contain(transported.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{1}}, 1, 2, qmomkit::CqmomNodeCountTransitionKind::grew));
}

void test_adaptive_state_source_growth_forces_adaptive_wheeler_without_mutating_qmom_options()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = collapsed_branch_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = branch_specific_ragged_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const qmomkit::CqmomCase cqmom_case = adaptive_cqmom_case_from_moments_with_fixed_qmom(initial_moments);
    const qmomkit::CqmomState state = cqmom_case.initial_state();
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    assert(cqmom_case.adaptive_branch_counts());
    assert(!cqmom_case.qmom_options().adaptive);
    assert(state.adaptive_branch_counts());
    assert(!state.qmom_options().adaptive);
    assert(active_counts_contain(state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, 2));
    assert(active_counts_contain(state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1, 1, 2));
    assert(active_counts_contain(state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 1, 1, 2));

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult transported =
        qmomkit::advance_transport_cqmom(state, 0.0, 1.0, zero_velocity, source_from_to(initial_required, target_required), options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->adaptive_branch_counts_enabled());
    assert(!transported.quadrature->has_rectangular_storage());
    assert(transported.quadrature->total_nodes() == 3);
    assert(!transported.state.qmom_options().adaptive);
    assert(active_counts_contain(transported.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, 2));
    assert(active_counts_contain(transported.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 1, 1, 2));
    assert(active_counts_contain(transported.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2, 1, 2));
    assert(diagnostics_contain_branch_count(transported.diagnostics.post_source_cqmom_diagnostics, 1, qmomkit::ConditioningPath{{1}}, 2, 2));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{1}}, 1, 2, qmomkit::CqmomNodeCountTransitionKind::grew));
    assert(transported.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(transported.diagnostics.post_source_cqmom_status == qmomkit::CqmomStatus::success);

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
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    qmomkit::CqmomTransportOptions source_last_options;
    source_last_options.scheme = qmomkit::TransportScheme::cd;
    source_last_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult first =
        qmomkit::advance_transport_cqmom(one_dimensional_adaptive_cqmom_state(initial_required), 0.0, 1.0, zero_velocity,
                                         source_from_to(initial_required, target_required), source_last_options);

    assert(first.ok());
    assert(first.quadrature.has_value());
    assert(first.quadrature->total_nodes() == 2);
    assert(active_counts_contain(first.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, configured_node_count));
    assert(diagnostics_contain_branch_count(first.diagnostics.post_source_cqmom_diagnostics, 0, qmomkit::ConditioningPath{}, 2, 2));
    assert_cqmom_state_moments_match_quadrature(first.state, *first.quadrature, 5.0e-10);

    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    qmomkit::CqmomTransportOptions growth_options;
    growth_options.scheme = qmomkit::TransportScheme::dc;
    growth_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult second =
        qmomkit::advance_transport_cqmom(first.state, first.final_time, 0.1, zero_velocity, zero_source, growth_options);

    assert(second.ok());
    assert(second.quadrature.has_value());
    assert(second.quadrature->total_nodes() == 2);
    assert(active_counts_contain(second.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, configured_node_count));
    assert(second.diagnostics.cqmom_node_count_transitions.empty());
    assert_cqmom_state_moments_match_quadrature(second.state, *second.quadrature, 5.0e-10);
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
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult transported = qmomkit::advance_transport_cqmom(
        one_dimensional_adaptive_cqmom_state(initial_required), 0.0, 1.0, zero_velocity, source_from_to(initial_required, target_required), options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->total_nodes() == 2);
    assert(active_counts_contain(transported.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, configured_node_count));
    assert(diagnostics_contain_branch_count(transported.diagnostics.post_source_cqmom_diagnostics, 0, qmomkit::ConditioningPath{}, 2, 2));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::initial_inversion, 0,
                               qmomkit::ConditioningPath{}, std::nullopt, 2, qmomkit::CqmomNodeCountTransitionKind::created));
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
    qmomkit::CqmomCoordinateSupport support(2);
    support[0] = qmomkit::BoundedInterval{0.0, 1.0};

    qmomkit::CqmomLeanCase lean_case{initial_required, {}, {true, 1.0e-10, 1.0e-10, 1.0e-12}, true};
    lean_case.coordinate_support = support;
    const qmomkit::CqmomState state = qmomkit::internal::create_cqmom_state(lean_case);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::CqmomTransportResult transported =
        qmomkit::advance_transport_cqmom(state, 0.0, 1.0, zero_velocity, source_from_to(initial_required, target_required), options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->total_nodes() == 2);
    assert(transported.state.coordinate_support().size() == 2);
    assert(transported.state.coordinate_support()[0].has_value());
    assert(active_counts_contain(transported.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, 3));
    assert(!transported.diagnostics.post_source_cqmom_diagnostics.empty());
    const qmomkit::CqmomDiagnostics& diagnostic = transported.diagnostics.post_source_cqmom_diagnostics.front();
    assert(diagnostic.bounded_support_checked);
    assert(diagnostic.bounded_candidate_attempt_type == "growth");
    assert(diagnostic.bounded_candidate_node_count == 3);
    assert(diagnostic.bounded_accepted_fallback_node_count == 2);
    assert(!diagnostic.bounded_support_determinants.empty());
}

void test_dense_and_lean_bounded_post_source_reduction_are_equivalent()
{
    constexpr std::size_t configured_node_count = 3;
    const qmomkit::PhaseSpaceHierarchy hierarchy({configured_node_count});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::RawMomentTensor initial_moments = one_dimensional_cqmom_moments(configured_node_count, {0.2, 0.5, 0.8}, {0.2, 0.3, 0.5});
    qmomkit::RawMomentTensor source_candidate_moments{{2 * configured_node_count}, std::vector<double>(2 * configured_node_count, 0.0)};
    const std::vector<double> feasible_prefix = raw_moments_1d({0.25, 0.75}, {0.45, 0.55}, 4);
    for (std::size_t order = 0; order < feasible_prefix.size(); ++order)
    {
        source_candidate_moments.values[order] = feasible_prefix[order];
    }
    source_candidate_moments.values[4] = 10.0;
    source_candidate_moments.values[5] = 10.0;

    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData source_candidate_required = qmomkit::extract_cqmom_moments(source_candidate_moments, layout);
    qmomkit::CqmomCoordinateSupport support(1);
    support[0] = qmomkit::BoundedInterval{0.0, 1.0};
    const qmomkit::QmomOptions adaptive_options{true, 1.0e-10, 1.0e-10, 1.0e-12};

    qmomkit::CqmomDenseCase dense_case{hierarchy, initial_moments, {}, adaptive_options, true};
    dense_case.coordinate_support = support;
    qmomkit::CqmomLeanCase lean_case{initial_required, {}, adaptive_options, true};
    lean_case.coordinate_support = support;

    const qmomkit::CqmomMomentSource source = source_from_to(initial_required, source_candidate_required);
    const qmomkit::CqmomTransportResult dense =
        qmomkit::advance_cqmom_discontinuous(qmomkit::internal::create_cqmom_state(dense_case), 0.0, 1.0, source);
    const qmomkit::CqmomTransportResult lean =
        qmomkit::advance_cqmom_discontinuous(qmomkit::internal::create_cqmom_state(lean_case), 0.0, 1.0, source);

    assert(dense.ok());
    assert(lean.ok());
    assert(dense.quadrature.has_value());
    assert(lean.quadrature.has_value());
    assert(dense.quadrature->total_nodes() == 2);
    assert(lean.quadrature->total_nodes() == 2);
    assert_cqmom_quadrature_near(*dense.quadrature, *lean.quadrature, 5.0e-10);
    assert(active_counts_equal(dense.state.active_branch_counts(), lean.state.active_branch_counts()));
    assert(active_counts_contain(dense.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2, 0, configured_node_count));
    assert(transitions_contain(dense.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 0,
                               qmomkit::ConditioningPath{}, 3, 2, qmomkit::CqmomNodeCountTransitionKind::reduced));
    assert_transition_diagnostics_equal(dense.diagnostics.cqmom_node_count_transitions, lean.diagnostics.cqmom_node_count_transitions);

    assert(dense.state.coordinate_support().size() == 1);
    assert(lean.state.coordinate_support().size() == 1);
    assert(dense.state.coordinate_support()[0].has_value());
    assert(lean.state.coordinate_support()[0].has_value());
    assert_near(dense.state.coordinate_support()[0]->lower, 0.0, 0.0);
    assert_near(lean.state.coordinate_support()[0]->upper, 1.0, 0.0);
    assert_cqmom_state_moments_match_quadrature(dense.state, *dense.quadrature, 5.0e-10);
    assert_cqmom_state_moments_match_quadrature(lean.state, *lean.quadrature, 5.0e-10);
    assert_common_state_moments_near(dense.state, lean.state, 5.0e-10);

    const qmomkit::RawMomentTensor dense_canonical = qmomkit::reconstruct_raw_moments(*dense.quadrature, hierarchy.moment_extents());
    const qmomkit::CqmomMomentData lean_canonical = qmomkit::extract_cqmom_moments(dense_canonical, layout);
    assert_vector_near(dense.state.moment_values(), dense_canonical.values, 5.0e-10);
    assert_vector_near(lean.state.moment_values(), lean_canonical.values, 5.0e-10);
    assert(std::abs(dense.state.moment_values()[4] - source_candidate_moments.values[4]) > 1.0);
    assert(std::abs(lean.state.moment_values()[4] - source_candidate_required.values[4]) > 1.0);
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
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult transported = qmomkit::advance_transport_cqmom(
        one_dimensional_adaptive_cqmom_state(initial_required), 0.0, 1.0, zero_velocity, source_from_to(initial_required, target_required), options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->total_nodes() == 1);
    assert(active_counts_contain(transported.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 1, 0, configured_node_count));
    assert(diagnostics_contain_branch_count(transported.diagnostics.post_source_cqmom_diagnostics, 0, qmomkit::ConditioningPath{}, 2, 1));
}

void test_adaptive_state_parent_growth_creates_unconstrained_descendants()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor initial_moments = root_collapsed_child_spread_cqmom_moments();
    const qmomkit::RawMomentTensor target_moments = fixed_cqmom_moments();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
    const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult transported =
        qmomkit::advance_transport_cqmom(adaptive_cqmom_case_from_moments(initial_moments).initial_state(), 0.0, 1.0, zero_velocity,
                                         source_from_to(initial_required, target_required), options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->total_nodes() == 4);
    assert(active_counts_contain(transported.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 2));
    assert(active_counts_contain(transported.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 2));
    assert(active_counts_contain(transported.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{1}}, 2));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 0,
                               qmomkit::ConditioningPath{}, 1, 2, qmomkit::CqmomNodeCountTransitionKind::grew));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{1}}, std::nullopt, 2, qmomkit::CqmomNodeCountTransitionKind::created));
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

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult transported =
        qmomkit::advance_transport_cqmom(adaptive_cqmom_case_from_moments(initial_moments).initial_state(), 0.0, 1.0, zero_velocity,
                                         source_from_to(initial_required, target_required), options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.quadrature->total_nodes() == 2);
    assert(active_counts_contain(transported.state.active_branch_counts(), 0, qmomkit::ConditioningPath{}, 1));
    assert(active_counts_contain(transported.state.active_branch_counts(), 1, qmomkit::ConditioningPath{{0}}, 2));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 0,
                               qmomkit::ConditioningPath{}, 2, 1, qmomkit::CqmomNodeCountTransitionKind::reduced));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{0}}, 1, std::nullopt, qmomkit::CqmomNodeCountTransitionKind::removed));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 1,
                               qmomkit::ConditioningPath{{0}}, std::nullopt, 2, qmomkit::CqmomNodeCountTransitionKind::created));
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

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult transported = qmomkit::advance_transport_cqmom(
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
                               qmomkit::ConditioningPath{{0}}, 2, 1, qmomkit::CqmomNodeCountTransitionKind::reduced));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 2,
                               qmomkit::ConditioningPath{{0, 0}}, 1, std::nullopt, qmomkit::CqmomNodeCountTransitionKind::removed));
    assert(transitions_contain(transported.diagnostics.cqmom_node_count_transitions, qmomkit::TransportStage::post_source_inversion, 2,
                               qmomkit::ConditioningPath{{0, 0}}, std::nullopt, 2, qmomkit::CqmomNodeCountTransitionKind::created));
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

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult transported = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, zero_source, options);

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
    assert(good_state.moment_storage() == qmomkit::CqmomMomentStorage::lean);
    const qmomkit::CqmomLeanCase* lean_case = qmomkit::internal::cqmom_state_lean_case(good_state);
    assert(lean_case != nullptr);
    std::vector<qmomkit::CqmomActiveBranchCount> bad_counts = good_state.active_branch_counts();
    assert(!bad_counts.empty());
    bad_counts.push_back(bad_counts.front());
    const qmomkit::CqmomState bad_state = qmomkit::internal::create_cqmom_state(*lean_case, std::move(bad_counts));
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult result = qmomkit::advance_transport_cqmom(bad_state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::initial_inversion);
    assert(result.diagnostics.initial_cqmom_status == qmomkit::CqmomStatus::invalid_input);
}

void test_adaptive_state_rejects_missing_active_count_metadata()
{
    const qmomkit::CqmomState good_state = adaptive_cqmom_case().initial_state();
    assert(good_state.moment_storage() == qmomkit::CqmomMomentStorage::lean);
    const qmomkit::CqmomLeanCase* lean_case = qmomkit::internal::cqmom_state_lean_case(good_state);
    assert(lean_case != nullptr);
    std::vector<qmomkit::CqmomActiveBranchCount> bad_counts = good_state.active_branch_counts();
    assert(!bad_counts.empty());
    bad_counts.pop_back();
    const qmomkit::CqmomState bad_state = qmomkit::internal::create_cqmom_state(*lean_case, std::move(bad_counts));
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult result = qmomkit::advance_transport_cqmom(bad_state, 0.0, 0.1, zero_velocity, zero_source, options);

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

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult failed = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, bad_source, options);
    assert(!failed.ok());
    assert(active_counts_equal(state.active_branch_counts(), original_counts));
    assert(active_counts_equal(failed.state.active_branch_counts(), original_counts));

    const qmomkit::CqmomTransportResult retry = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, zero_source, options);
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

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult failed = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, bad_realizability_source, options);
    assert(!failed.ok());
    assert(failed.status == qmomkit::TransportStatus::inversion_failed);
    assert(active_counts_equal(state.active_branch_counts(), original_counts));
    assert(active_counts_equal(failed.state.active_branch_counts(), original_counts));

    const qmomkit::CqmomTransportResult retry = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, zero_source, options);
    assert(retry.ok());
    assert(active_counts_equal(retry.state.active_branch_counts(), original_counts));
}

void assert_dense_and_lean_adaptive_results_match(const qmomkit::CqmomTransportResult& dense, const qmomkit::CqmomTransportResult& lean,
                                                  const double tolerance)
{
    assert(dense.ok());
    assert(lean.ok());
    assert(dense.state.moment_storage() == qmomkit::CqmomMomentStorage::dense);
    assert(lean.state.moment_storage() == qmomkit::CqmomMomentStorage::lean);
    assert(dense.quadrature.has_value());
    assert(lean.quadrature.has_value());
    assert_cqmom_quadrature_near(*dense.quadrature, *lean.quadrature, tolerance);
    assert(active_counts_equal(dense.state.active_branch_counts(), lean.state.active_branch_counts()));
    assert_common_state_moments_near(dense.state, lean.state, tolerance);
    assert_cqmom_state_moments_match_quadrature(dense.state, *dense.quadrature, tolerance);
    assert_cqmom_state_moments_match_quadrature(lean.state, *lean.quadrature, tolerance);
    assert(dense.diagnostics.initial_cqmom_status == lean.diagnostics.initial_cqmom_status);
    assert(dense.diagnostics.post_source_cqmom_status == lean.diagnostics.post_source_cqmom_status);
    assert(dense.diagnostics.source_status == lean.diagnostics.source_status);
    assert_transition_diagnostics_equal(dense.diagnostics.cqmom_node_count_transitions, lean.diagnostics.cqmom_node_count_transitions);
}

void test_dense_and_lean_adaptive_state_active_count_persistence_are_equivalent()
{
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };

    {
        const qmomkit::CqmomCase dense_case = adaptive_cqmom_case_from_moments_with_storage(branch_specific_ragged_cqmom_moments(), "dense");
        const qmomkit::CqmomCase lean_case = adaptive_cqmom_case_from_moments_with_storage(branch_specific_ragged_cqmom_moments(), "lean");
        const qmomkit::CqmomState dense_initial = dense_case.initial_state();
        const qmomkit::CqmomState lean_initial = lean_case.initial_state();
        assert(dense_initial.moment_storage() == qmomkit::CqmomMomentStorage::dense);
        assert(lean_initial.moment_storage() == qmomkit::CqmomMomentStorage::lean);
        assert(active_counts_equal(dense_initial.active_branch_counts(), lean_initial.active_branch_counts()));

        qmomkit::CqmomTransportOptions options;
        options.scheme = qmomkit::TransportScheme::dc;
        options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
        const qmomkit::CqmomTransportResult dense_first =
            qmomkit::advance_transport_cqmom(dense_initial, 0.0, 0.1, zero_velocity, zero_source, options);
        const qmomkit::CqmomTransportResult lean_first =
            qmomkit::advance_transport_cqmom(lean_initial, 0.0, 0.1, zero_velocity, zero_source, options);
        assert_dense_and_lean_adaptive_results_match(dense_first, lean_first, 5.0e-10);
        assert(dense_first.diagnostics.cqmom_node_count_transitions.empty());
        assert(lean_first.diagnostics.cqmom_node_count_transitions.empty());

        const qmomkit::CqmomTransportResult dense_second =
            qmomkit::advance_transport_cqmom(dense_first.state, dense_first.final_time, 0.1, zero_velocity, zero_source, options);
        const qmomkit::CqmomTransportResult lean_second =
            qmomkit::advance_transport_cqmom(lean_first.state, lean_first.final_time, 0.1, zero_velocity, zero_source, options);
        assert_dense_and_lean_adaptive_results_match(dense_second, lean_second, 5.0e-10);
    }

    const auto run_source_case = [&](const qmomkit::RawMomentTensor& initial_moments, const qmomkit::RawMomentTensor& target_moments,
                                     const qmomkit::PhaseSpaceHierarchy& hierarchy, const qmomkit::CqmomTransportOptions& options)
    {
        const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
        const qmomkit::CqmomMomentData initial_required = qmomkit::extract_cqmom_moments(initial_moments, layout);
        const qmomkit::CqmomMomentData target_required = qmomkit::extract_cqmom_moments(target_moments, layout);
        const qmomkit::CqmomCase dense_case = adaptive_cqmom_case_from_moments_with_hierarchy_and_storage(initial_moments, hierarchy.node_counts(),
                                                                                                          hierarchy.conditioning_order(), "dense");
        const qmomkit::CqmomCase lean_case = adaptive_cqmom_case_from_moments_with_hierarchy_and_storage(initial_moments, hierarchy.node_counts(),
                                                                                                         hierarchy.conditioning_order(), "lean");
        const qmomkit::CqmomMomentSource source = source_from_to(initial_required, target_required);
        const qmomkit::CqmomTransportResult dense =
            qmomkit::advance_transport_cqmom(dense_case.initial_state(), 0.0, 1.0, zero_velocity, source, options);
        const qmomkit::CqmomTransportResult lean =
            qmomkit::advance_transport_cqmom(lean_case.initial_state(), 0.0, 1.0, zero_velocity, source, options);
        assert_dense_and_lean_adaptive_results_match(dense, lean, 5.0e-9);
        assert_dense_state_matches_reconstructed_quadrature(dense, initial_moments.extents, 5.0e-9);
    };

    qmomkit::CqmomTransportOptions preserve_options;
    preserve_options.scheme = qmomkit::TransportScheme::dc;
    preserve_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    run_source_case(branch_specific_ragged_cqmom_moments(), collapsed_branch_cqmom_moments(), qmomkit::PhaseSpaceHierarchy({2, 2}), preserve_options);
    run_source_case(branch_specific_ragged_cqmom_moments(), root_collapsed_child_spread_cqmom_moments(), qmomkit::PhaseSpaceHierarchy({2, 2}),
                    preserve_options);
    run_source_case(second_level_split_3d_cqmom_moments(), second_level_collapsed_grandchild_spread_3d_cqmom_moments(),
                    qmomkit::PhaseSpaceHierarchy({2, 2, 2}, {0, 1, 2}), preserve_options);

    qmomkit::CqmomTransportOptions growth_options;
    growth_options.scheme = qmomkit::TransportScheme::dc;
    growth_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    run_source_case(collapsed_branch_cqmom_moments(), branch_specific_ragged_cqmom_moments(), qmomkit::PhaseSpaceHierarchy({2, 2}), growth_options);

    {
        const qmomkit::CqmomCase dense_case = adaptive_cqmom_case_from_moments_with_storage(adaptive_cqmom_moments(), "dense");
        const qmomkit::CqmomCase lean_case = adaptive_cqmom_case_from_moments_with_storage(adaptive_cqmom_moments(), "lean");
        const qmomkit::CqmomState dense_state = dense_case.initial_state();
        const qmomkit::CqmomState lean_state = lean_case.initial_state();
        const std::vector<qmomkit::CqmomActiveBranchCount> dense_counts = dense_state.active_branch_counts();
        const std::vector<qmomkit::CqmomActiveBranchCount> lean_counts = lean_state.active_branch_counts();
        const auto bad_source = [](double, const qmomkit::CqmomSourceContext& moments)
        { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; };
        const qmomkit::CqmomTransportResult dense_failed =
            qmomkit::advance_transport_cqmom(dense_state, 0.0, 0.1, zero_velocity, bad_source, preserve_options);
        const qmomkit::CqmomTransportResult lean_failed =
            qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.1, zero_velocity, bad_source, preserve_options);
        assert(!dense_failed.ok());
        assert(!lean_failed.ok());
        assert(dense_failed.status == lean_failed.status);
        assert(dense_failed.diagnostics.failed_stage == lean_failed.diagnostics.failed_stage);
        assert(active_counts_equal(dense_failed.state.active_branch_counts(), dense_counts));
        assert(active_counts_equal(lean_failed.state.active_branch_counts(), lean_counts));
        assert(active_counts_equal(dense_counts, lean_counts));
        assert_vector_near(dense_failed.state.moment_values(), dense_state.moment_values(), 0.0);
        assert_vector_near(lean_failed.state.moment_values(), lean_state.moment_values(), 0.0);

        const auto bad_realizability_source = [](double, const qmomkit::CqmomSourceContext& moments)
        {
            std::vector<double> rate(moments.size(), 0.0);
            rate.front() = -20.0;
            return qmomkit::CqmomMomentRate{rate};
        };
        const qmomkit::CqmomTransportResult dense_post_source_failed =
            qmomkit::advance_transport_cqmom(dense_state, 0.0, 0.1, zero_velocity, bad_realizability_source, preserve_options);
        const qmomkit::CqmomTransportResult lean_post_source_failed =
            qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.1, zero_velocity, bad_realizability_source, preserve_options);
        assert(!dense_post_source_failed.ok());
        assert(!lean_post_source_failed.ok());
        assert(dense_post_source_failed.status == lean_post_source_failed.status);
        assert(dense_post_source_failed.diagnostics.failed_stage == lean_post_source_failed.diagnostics.failed_stage);
        assert(active_counts_equal(dense_post_source_failed.state.active_branch_counts(), dense_counts));
        assert(active_counts_equal(lean_post_source_failed.state.active_branch_counts(), lean_counts));
        assert_vector_near(dense_post_source_failed.state.moment_values(), dense_state.moment_values(), 0.0);
        assert_vector_near(lean_post_source_failed.state.moment_values(), lean_state.moment_values(), 0.0);
    }
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
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& state_values)
    { return qmomkit::CqmomMomentRate{std::vector<double>(state_values.size(), 0.0)}; };

    qmomkit::TestDenseCqmomTransportOptions dense_options;
    dense_options.qmom_options = {true, 1.0e-10, 1.0e-10, 1.0e-12};
    dense_options.adaptive_branch_counts = true;
    dense_options.node_options = tight_node_options();
    dense_options.scheme = qmomkit::TransportScheme::strang;
    const qmomkit::TestDenseCqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(moments, hierarchy, 0.0, 0.2, velocity, zero_source, dense_options);

    qmomkit::CqmomTransportOptions state_options;
    state_options.node_options = dense_options.node_options;
    state_options.scheme = dense_options.scheme;
    const qmomkit::CqmomTransportResult transported = qmomkit::advance_transport_cqmom(state, 0.0, 0.2, velocity, zero_source, state_options);

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
        static_cast<void>(parsed.value().initial_state());
    }
    catch (const std::runtime_error& error)
    {
        threw = true;
        assert(std::string(error.what()).find("quadrature reconstruction failed") != std::string::npos);
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
    const qmomkit::TestLeanCqmomTransportState lean_state = lean_state_from_dense(moments, hierarchy);
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

    qmomkit::CqmomTransportOptions state_options;
    state_options.node_options = tight_node_options();
    state_options.scheme = qmomkit::TransportScheme::strang;
    const qmomkit::CqmomTransportResult state_result = qmomkit::advance_transport_cqmom(state, 0.0, 0.2, velocity, source, state_options);

    qmomkit::TestLeanCqmomTransportOptions lean_options;
    lean_options.node_options = state_options.node_options;
    lean_options.scheme = state_options.scheme;
    const qmomkit::TestLeanCqmomTransportResult lean_result = qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.2, velocity, source, lean_options);

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

void test_configured_dense_and_lean_state_transport_are_equivalent()
{
    const qmomkit::CqmomCase dense_case = fixed_cqmom_case_with_storage("dense");
    const qmomkit::CqmomCase lean_case = fixed_cqmom_case_with_storage("lean");
    assert(dense_case.moment_storage() == qmomkit::CqmomMomentStorage::dense);
    assert(lean_case.moment_storage() == qmomkit::CqmomMomentStorage::lean);

    const auto velocity = [](double, const std::vector<double>& xi)
    {
        assert(xi.size() == 2);
        return std::vector<double>{0.15 + 0.03 * xi[0], -0.08 + 0.02 * xi[1]};
    };
    const qmomkit::CqmomMomentSource source = [](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        for (std::size_t entry = 0; entry < moments.size(); ++entry)
        {
            rate[entry] = -0.04 * moments.value(entry);
        }
        return qmomkit::CqmomMomentRate{rate};
    };

    qmomkit::CqmomTransportOptions options;
    options.node_options = tight_node_options();
    options.scheme = qmomkit::TransportScheme::cd;

    const qmomkit::CqmomTransportResult dense = qmomkit::advance_transport_cqmom(dense_case.initial_state(), 0.0, 0.1, velocity, source, options);
    const qmomkit::CqmomTransportResult lean = qmomkit::advance_transport_cqmom(lean_case.initial_state(), 0.0, 0.1, velocity, source, options);

    assert(dense.ok());
    assert(lean.ok());
    assert(dense.quadrature.has_value());
    assert(lean.quadrature.has_value());
    assert(dense.diagnostics.cqmom_moment_storage == qmomkit::CqmomMomentStorage::dense);
    assert(lean.diagnostics.cqmom_moment_storage == qmomkit::CqmomMomentStorage::lean);
    assert_cqmom_quadrature_near(*dense.quadrature, *lean.quadrature, 5.0e-9);

    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(dense_case.hierarchy());
    const qmomkit::RawMomentTensor dense_moments{dense_case.hierarchy().moment_extents(), dense.state.moment_values()};
    const qmomkit::CqmomMomentData dense_required = qmomkit::extract_cqmom_moments(dense_moments, layout);
    assert_vector_near(dense_required.values, lean.state.moment_values(), 5.0e-9);
}

void test_dense_and_lean_state_source_contexts_are_identical()
{
    const qmomkit::CqmomCase dense_case = fixed_cqmom_case_with_storage("dense");
    const qmomkit::CqmomCase lean_case = fixed_cqmom_case_with_storage("lean");
    SourceContextSnapshot dense_snapshot;
    SourceContextSnapshot lean_snapshot;

    const auto make_source = [](SourceContextSnapshot& snapshot)
    {
        return [&snapshot](const double time, const qmomkit::CqmomSourceContext& context)
        {
            capture_source_context_once(snapshot, context);
            const qmomkit::MultiIndex beta(context.dimension_count(), 0);
            const double coupled_moment = context.value(beta).value_or(0.0);
            std::vector<double> rate(context.size(), 0.0);
            for (std::size_t entry = 0; entry < context.size(); ++entry)
            {
                const qmomkit::MultiIndex& alpha = context.index(entry);
                double order_sum = 0.0;
                for (const std::size_t order : alpha)
                {
                    order_sum += static_cast<double>(order);
                }
                const double coefficient = -0.015 * (1.0 + order_sum);
                rate[entry] = coefficient * context.value(entry) + 0.0025 * coupled_moment + 0.001 * time * (1.0 + order_sum);
            }
            return qmomkit::CqmomMomentRate{rate};
        };
    };

    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    constexpr double start_time = 0.125;
    constexpr double dt = 0.02;

    const qmomkit::CqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(dense_case.initial_state(), start_time, dt, zero_velocity, make_source(dense_snapshot), options);
    const qmomkit::CqmomTransportResult lean =
        qmomkit::advance_transport_cqmom(lean_case.initial_state(), start_time, dt, zero_velocity, make_source(lean_snapshot), options);

    assert(dense.ok());
    assert(lean.ok());
    assert(dense.quadrature.has_value());
    assert(lean.quadrature.has_value());
    assert_source_context_snapshots_near(dense_snapshot, lean_snapshot, 5.0e-13);
    assert_cqmom_quadrature_near(*dense.quadrature, *lean.quadrature, 5.0e-9);
    assert_common_state_moments_near(dense.state, lean.state, 5.0e-9);
    assert(dense.diagnostics.source_status == lean.diagnostics.source_status);
    assert(dense.diagnostics.post_source_cqmom_status == lean.diagnostics.post_source_cqmom_status);
}

qmomkit::CqmomMomentSource proportional_cqmom_source(const double coefficient)
{
    return [coefficient](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        for (std::size_t entry = 0; entry < moments.size(); ++entry)
        {
            rate[entry] = coefficient * moments.value(entry);
        }
        return qmomkit::CqmomMomentRate{rate};
    };
}

qmomkit::CqmomCase fixed_cqmom_case_with_hierarchy_and_storage(const qmomkit::RawMomentTensor& moments,
                                                               const std::vector<std::size_t>& conditioning_order, const char* const storage)
{
    std::string text = cqmom_case_text_from_moments_with_hierarchy(moments, {2, 2}, conditioning_order, false, false);
    text += "\n[cqmom]\nmoment_storage = ";
    text += storage;
    text += "\n";
    const auto parsed = qmomkit::config::parse_cqmom(text, "fixed-cqmom-case-with-hierarchy-and-storage");
    assert(parsed.ok());
    return parsed.value();
}

void assert_dense_and_lean_source_results_match(const qmomkit::CqmomTransportResult& dense, const qmomkit::CqmomTransportResult& lean,
                                                const std::vector<std::size_t>& dense_extents, const double tolerance)
{
    assert(dense.ok());
    assert(lean.ok());
    assert(dense.quadrature.has_value());
    assert(lean.quadrature.has_value());
    assert_cqmom_quadrature_near(*dense.quadrature, *lean.quadrature, tolerance);
    assert_common_state_moments_near(dense.state, lean.state, tolerance);
    assert(active_counts_equal(dense.state.active_branch_counts(), lean.state.active_branch_counts()));
    assert(dense.status == lean.status);
    assert(dense.diagnostics.source_status == lean.diagnostics.source_status);
    assert(dense.diagnostics.post_source_cqmom_status == lean.diagnostics.post_source_cqmom_status);
    assert_dense_state_matches_reconstructed_quadrature(dense, dense_extents, tolerance);
}

void test_dense_state_reconstructs_dense_only_moments_after_discontinuous_source()
{
    const qmomkit::RawMomentTensor initial = fixed_cqmom_moments();
    const qmomkit::CqmomCase dense_case = fixed_cqmom_case_with_storage("dense");
    const qmomkit::CqmomCase lean_case = fixed_cqmom_case_with_storage("lean");
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomMomentSource source = proportional_cqmom_source(0.35);
    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::cd;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    options.source_options.ssp_rk3_options.estimate_error = true;
    const std::size_t expected_source_state_size = qmomkit::create_cqmom_moment_layout(dense_case.hierarchy()).required_indices.size();
    std::size_t stage_accept_calls = 0;
    std::size_t accept_calls = 0;
    options.source_options.ssp_rk3_options.stage_accept = [&](qmomkit::SspRk3Stage, const std::vector<double>& state)
    {
        ++stage_accept_calls;
        assert(state.size() == expected_source_state_size);
        return true;
    };
    options.source_options.ssp_rk3_options.accept = [&](const std::vector<double>& state)
    {
        ++accept_calls;
        assert(state.size() == expected_source_state_size);
        return true;
    };

    const qmomkit::CqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(dense_case.initial_state(), 0.0, 0.1, zero_velocity, source, options);
    const qmomkit::CqmomTransportResult lean = qmomkit::advance_transport_cqmom(lean_case.initial_state(), 0.0, 0.1, zero_velocity, source, options);

    assert_dense_and_lean_source_results_match(dense, lean, initial.extents, 5.0e-10);
    assert(dense_only_moment_changed(initial, dense.state, dense_case.hierarchy(), 1.0e-8));
    assert(stage_accept_calls > 0);
    assert(accept_calls > 0);
}

void test_dense_state_zero_source_keeps_full_dense_closure_consistent()
{
    const qmomkit::RawMomentTensor initial = fixed_cqmom_moments();
    const qmomkit::CqmomCase dense_case = fixed_cqmom_case_with_storage("dense");
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size(), 0.0)}; };
    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::cd;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::CqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(dense_case.initial_state(), 0.0, 0.1, zero_velocity, zero_source, options);

    assert_dense_state_matches_reconstructed_quadrature(dense, initial.extents, 5.0e-10);
    assert_vector_near(dense.state.moment_values(), initial.values, 5.0e-10);
}

void test_dense_state_reconstructs_dense_only_moments_after_source_and_boundary()
{
    const qmomkit::RawMomentTensor initial = fixed_cqmom_moments();
    const qmomkit::CqmomCase dense_case = fixed_cqmom_case_with_storage("dense");
    const qmomkit::CqmomCase lean_case = fixed_cqmom_case_with_storage("lean");
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomMomentSource source = proportional_cqmom_source(0.10);
    const qmomkit::CqmomMomentBoundaryContribution boundary = proportional_cqmom_source(0.15);
    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::cd;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::CqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(dense_case.initial_state(), 0.0, 0.1, zero_velocity, source, boundary, options);
    const qmomkit::CqmomTransportResult lean =
        qmomkit::advance_transport_cqmom(lean_case.initial_state(), 0.0, 0.1, zero_velocity, source, boundary, options);

    assert_dense_and_lean_source_results_match(dense, lean, initial.extents, 5.0e-10);
    assert(dense_only_moment_changed(initial, dense.state, dense_case.hierarchy(), 1.0e-8));
}

void test_dense_state_reconstructs_dense_only_moments_after_repeated_source_steps()
{
    const qmomkit::RawMomentTensor initial = fixed_cqmom_moments();
    const qmomkit::CqmomCase dense_case = fixed_cqmom_case_with_storage("dense");
    const qmomkit::CqmomCase lean_case = fixed_cqmom_case_with_storage("lean");
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomMomentSource source = proportional_cqmom_source(0.20);
    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::cd;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    qmomkit::CqmomTransportResult dense = qmomkit::advance_transport_cqmom(dense_case.initial_state(), 0.0, 0.05, zero_velocity, source, options);
    qmomkit::CqmomTransportResult lean = qmomkit::advance_transport_cqmom(lean_case.initial_state(), 0.0, 0.05, zero_velocity, source, options);
    assert_dense_and_lean_source_results_match(dense, lean, initial.extents, 5.0e-10);

    dense = qmomkit::advance_transport_cqmom(dense.state, dense.final_time, 0.05, zero_velocity, source, options);
    lean = qmomkit::advance_transport_cqmom(lean.state, lean.final_time, 0.05, zero_velocity, source, options);
    assert_dense_and_lean_source_results_match(dense, lean, initial.extents, 5.0e-10);
    assert(dense_only_moment_changed(initial, dense.state, dense_case.hierarchy(), 1.0e-8));
}

void test_dense_state_reconstructs_dense_only_moments_with_nonidentity_conditioning_order()
{
    const qmomkit::RawMomentTensor initial = fixed_cqmom_moments();
    const qmomkit::CqmomCase dense_case = fixed_cqmom_case_with_hierarchy_and_storage(initial, {1, 0}, "dense");
    const qmomkit::CqmomCase lean_case = fixed_cqmom_case_with_hierarchy_and_storage(initial, {1, 0}, "lean");
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::CqmomMomentSource source = proportional_cqmom_source(0.25);
    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::cd;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::CqmomTransportResult dense =
        qmomkit::advance_transport_cqmom(dense_case.initial_state(), 0.0, 0.1, zero_velocity, source, options);
    const qmomkit::CqmomTransportResult lean = qmomkit::advance_transport_cqmom(lean_case.initial_state(), 0.0, 0.1, zero_velocity, source, options);

    assert_dense_and_lean_source_results_match(dense, lean, initial.extents, 5.0e-10);
    assert(dense_only_moment_changed(initial, dense.state, dense_case.hierarchy(), 1.0e-8));
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

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult transported = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, zero_source, options);

    assert(transported.ok());
    assert(transported.quadrature.has_value());
    assert(transported.state.vandermonde_solver() == qmomkit::VandermondeSolver::dense);
    assert(transported.state.vandermonde_solver() == state.vandermonde_solver());
    assert_cqmom_quadrature_near(*transported.quadrature, initial.quadrature, 5.0e-10);
}

qmomkit::RawMomentTensor transport_vandermonde_policy_moments()
{
    return fixed_cqmom_moments();
}

bool cqmom_diagnostics_saw_dense_fallback(const std::vector<qmomkit::CqmomDiagnostics>& diagnostics)
{
    for (const qmomkit::CqmomDiagnostics& diagnostic : diagnostics)
    {
        if (diagnostic.dense_vandermonde_fallback_occurred)
        {
            return true;
        }
    }
    return false;
}

void test_cqmom_transport_vandermonde_policy_independent_of_diagnostic_level()
{
    const qmomkit::RawMomentTensor moments = transport_vandermonde_policy_moments();
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    std::optional<qmomkit::CqmomTransportResult> reference_dense;
    std::optional<qmomkit::CqmomTransportResult> reference_lean;
    for (const qmomkit::CqmomDiagnosticLevel level :
         {qmomkit::CqmomDiagnosticLevel::none, qmomkit::CqmomDiagnosticLevel::basic, qmomkit::CqmomDiagnosticLevel::detailed})
    {
        const qmomkit::CqmomDenseCase dense_case{
            hierarchy, moments, {}, {}, false, level, qmomkit::VandermondeSolver::bjorck_pereyra_with_dense_fallback};
        const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
        const qmomkit::CqmomLeanCase lean_case{
            qmomkit::extract_cqmom_moments(moments, layout), {}, {}, false, level, qmomkit::VandermondeSolver::bjorck_pereyra_with_dense_fallback};

        const qmomkit::CqmomTransportResult dense =
            qmomkit::advance_transport_cqmom(qmomkit::internal::create_cqmom_state(dense_case), 0.0, 0.1, zero_velocity, zero_source, options);
        const qmomkit::CqmomTransportResult lean =
            qmomkit::advance_transport_cqmom(qmomkit::internal::create_cqmom_state(lean_case), 0.0, 0.1, zero_velocity, zero_source, options);

        assert(dense.status == lean.status);
        assert(dense.ok() == lean.ok());
        assert(dense.diagnostics.failed_stage == lean.diagnostics.failed_stage);
        assert(dense.diagnostics.initial_cqmom_status == lean.diagnostics.initial_cqmom_status);
        assert(dense.diagnostics.post_source_cqmom_status == lean.diagnostics.post_source_cqmom_status);
        assert(dense.quadrature.has_value() == lean.quadrature.has_value());
        if (dense.quadrature.has_value())
        {
            assert_cqmom_quadrature_near(*dense.quadrature, *lean.quadrature, 1.0e-5);
        }
        if (dense.ok())
        {
            assert_common_cqmom_moments_near(dense.state, lean.state, 1.0e-5);
        }
        assert(active_counts_equal(dense.state.active_branch_counts(), lean.state.active_branch_counts()));

        if (level == qmomkit::CqmomDiagnosticLevel::none)
        {
            assert(dense.diagnostics.initial_cqmom_diagnostics.empty());
            assert(lean.diagnostics.initial_cqmom_diagnostics.empty());
            assert(dense.diagnostics.post_source_cqmom_diagnostics.empty());
            assert(lean.diagnostics.post_source_cqmom_diagnostics.empty());
        }
        else
        {
            assert(!cqmom_diagnostics_saw_dense_fallback(dense.diagnostics.initial_cqmom_diagnostics));
            assert(!cqmom_diagnostics_saw_dense_fallback(lean.diagnostics.initial_cqmom_diagnostics));
            assert(!cqmom_diagnostics_saw_dense_fallback(dense.diagnostics.post_source_cqmom_diagnostics));
            assert(!cqmom_diagnostics_saw_dense_fallback(lean.diagnostics.post_source_cqmom_diagnostics));
        }

        if (!reference_dense.has_value())
        {
            reference_dense = dense;
            reference_lean = lean;
            continue;
        }

        assert(dense.status == reference_dense->status);
        assert(lean.status == reference_lean->status);
        assert(dense.ok() == reference_dense->ok());
        assert(lean.ok() == reference_lean->ok());
        assert(dense.diagnostics.failed_stage == reference_dense->diagnostics.failed_stage);
        assert(lean.diagnostics.failed_stage == reference_lean->diagnostics.failed_stage);
        assert(dense.quadrature.has_value() == reference_dense->quadrature.has_value());
        assert(lean.quadrature.has_value() == reference_lean->quadrature.has_value());
        if (dense.quadrature.has_value())
        {
            assert_cqmom_quadrature_near(*dense.quadrature, *reference_dense->quadrature, 1.0e-5);
        }
        if (lean.quadrature.has_value())
        {
            assert_cqmom_quadrature_near(*lean.quadrature, *reference_lean->quadrature, 1.0e-5);
        }
        if (dense.ok())
        {
            assert_common_cqmom_moments_near(dense.state, reference_dense->state, 1.0e-5);
        }
        if (lean.ok())
        {
            assert_common_cqmom_moments_near(lean.state, reference_lean->state, 1.0e-5);
        }
        assert(active_counts_equal(dense.state.active_branch_counts(), reference_dense->state.active_branch_counts()));
        assert(active_counts_equal(lean.state.active_branch_counts(), reference_lean->state.active_branch_counts()));
    }
}

void test_state_cqmom_source_evolution_matches_lean()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::TestLeanCqmomTransportState lean_state = lean_state_from_dense(fixed_cqmom_moments(), hierarchy);
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

    qmomkit::CqmomTransportOptions state_options;
    state_options.scheme = qmomkit::TransportScheme::dc;
    state_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult state_result = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, source, state_options);

    qmomkit::TestLeanCqmomTransportOptions lean_options;
    lean_options.scheme = state_options.scheme;
    lean_options.node_options = state_options.node_options;
    const qmomkit::TestLeanCqmomTransportResult lean_result =
        qmomkit::advance_transport_cqmom(lean_state, 0.0, 0.1, zero_velocity, source, lean_options);

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

    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, bad_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.message.find("Lean CQMOM") == std::string::npos);
    assert(result.diagnostics.message.find("Lean CQMOM") == std::string::npos);
}

void test_state_cqmom_transport_api_does_not_expose_moments()
{
    static_assert(!HasPublicMomentsMember<qmomkit::CqmomTransportOptions>::value);
    static_assert(!HasPublicMomentsMember<qmomkit::CqmomTransportResult>::value);
    static_assert(std::is_same_v<decltype(std::declval<qmomkit::CqmomTransportResult&>().state), qmomkit::CqmomState>);
}

void test_state_qmom_transport_api_does_not_expose_moments()
{
    static_assert(!HasPublicMomentsMember<qmomkit::QmomState>::value);
    static_assert(!HasPublicQuadratureMember<qmomkit::QmomState>::value);
    static_assert(!HasPublicMomentsMember<qmomkit::QmomTransportResult>::value);
    static_assert(std::is_same_v<decltype(std::declval<qmomkit::QmomTransportResult&>().state), qmomkit::QmomState>);
    static_assert(std::is_copy_constructible_v<qmomkit::QmomState>);
    static_assert(std::is_move_constructible_v<qmomkit::QmomState>);
    static_assert(std::is_copy_assignable_v<qmomkit::QmomState>);
    static_assert(std::is_move_assignable_v<qmomkit::QmomState>);
    static_assert(!std::is_default_constructible_v<qmomkit::QmomState>);
    static_assert(!std::is_aggregate_v<qmomkit::QmomState>);
    static_assert(!std::is_default_constructible_v<qmomkit::QmomTransportResult>);
    static_assert(std::is_copy_constructible_v<qmomkit::QmomTransportResult>);
    static_assert(std::is_move_constructible_v<qmomkit::QmomTransportResult>);
}

void test_qmom_state_copy_move_workflow()
{
    qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.0}, {0.4, 0.6}, 4), 2);
    qmomkit::QmomState copied = state;
    qmomkit::QmomState moved = std::move(copied);
    assert_qmom_state_canonical(moved, 1.0e-12);

    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const std::vector<double>& moments) { return std::vector<double>(moments.size(), 0.0); };
    qmomkit::QmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::cd;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(moved, 2, 0.0, 0.1, zero_velocity, zero_source, options);
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

    assert(state.active_node_count() == 2);
    assert_qmom_state_canonical(state, 1.0e-12);
}

void test_qmom_adaptive_transport_preserves_active_count_with_zero_source()
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = true;
    qmom_options.rmin = 1.0e-10;
    qmom_options.eabs = 1.0e-12;
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.75, 1.25}, {0.4, 0.6}, 6), 3, qmom_options);
    assert(state.active_node_count() == 2);

    const auto velocity = [](double, const std::vector<double>& xi) { return std::vector<double>{0.2 * xi[0]}; };
    const auto zero_source = [](double, const std::vector<double>& moments) { return std::vector<double>(moments.size(), 0.0); };
    qmomkit::QmomTransportOptions options;
    options.qmom_options = qmom_options;
    options.scheme = qmomkit::TransportScheme::cd;
    options.node_options = tight_node_options();

    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(state, 3, 0.0, 0.2, velocity, zero_source, options);

    assert(result.ok());
    assert(result.state.active_node_count() == 2);
    assert_qmom_state_canonical(result.state, 5.0e-10);
}

void test_qmom_adaptive_post_source_reduction_and_growth_policy()
{
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = true;
    qmom_options.rmin = 1.0e-10;
    qmom_options.eabs = 1.0e-12;
    qmomkit::QmomTransportOptions options;
    options.qmom_options = qmom_options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    const auto source_from_to = [](const std::vector<double>& initial, const std::vector<double>& target)
    {
        return [initial, target](double, const std::vector<double>& state)
        {
            assert(state.size() == initial.size());
            std::vector<double> rate(state.size(), 0.0);
            for (std::size_t i = 0; i < rate.size(); ++i)
            {
                rate[i] = target[i] - initial[i];
            }
            return rate;
        };
    };

    const qmomkit::QmomState two_node_state = qmom_state_from_moments(raw_moments_1d({-0.75, 1.25}, {0.4, 0.6}, 6), 3, qmom_options);
    const std::vector<double> one_node_moments = raw_moments_1d({0.5}, {1.0}, 6);
    const qmomkit::QmomTransportResult reduction = qmomkit::advance_transport_qmom(
        two_node_state, 3, 0.0, 1.0, zero_velocity, source_from_to(two_node_state.moments(), one_node_moments), options);
    assert(reduction.ok());
    assert(reduction.state.active_node_count() == 1);
    assert_qmom_state_canonical(reduction.state, 1.0e-12);

    const qmomkit::QmomState one_node_state = qmom_state_from_moments(one_node_moments, 3, qmom_options);
    const std::vector<double> two_node_moments = raw_moments_1d({-0.25, 1.5}, {0.35, 0.65}, 6);
    const qmomkit::QmomTransportResult growth = qmomkit::advance_transport_qmom(one_node_state, 3, 0.0, 1.0, zero_velocity,
                                                                                source_from_to(one_node_state.moments(), two_node_moments), options);
    assert(growth.ok());
    assert(growth.state.active_node_count() == 2);
    assert_qmom_state_canonical(growth.state, 1.0e-12);

    const std::vector<double> three_node_moments = raw_moments_1d({-1.0, 0.25, 1.5}, {0.2, 0.3, 0.5}, 6);
    const qmomkit::QmomTransportResult capped_growth = qmomkit::advance_transport_qmom(
        one_node_state, 3, 0.0, 1.0, zero_velocity, source_from_to(one_node_state.moments(), three_node_moments), options);
    assert(capped_growth.ok());
    assert(capped_growth.state.active_node_count() == 2);
    assert_qmom_state_canonical(capped_growth.state, 1.0e-12);

    const std::vector<double> failed_growth_moments{1.0, 0.8, 2.8, 4.4, 0.0, 0.0};
    const qmomkit::QmomTransportResult failed_growth = qmomkit::advance_transport_qmom(
        two_node_state, 3, 0.0, 1.0, zero_velocity, source_from_to(two_node_state.moments(), failed_growth_moments), options);
    assert(failed_growth.ok());
    assert(failed_growth.state.active_node_count() == 2);
    assert_qmom_state_canonical(failed_growth.state, 1.0e-12);
}

void test_qmom_fixed_post_source_requires_configured_count()
{
    const qmomkit::QmomState state = qmom_state_from_moments(raw_moments_1d({-0.5, 1.0}, {0.4, 0.6}, 4), 2);
    const std::vector<double> one_node_moments = raw_moments_1d({0.0}, {1.0}, 4);
    const auto source = [initial = state.moments(), one_node_moments](double, const std::vector<double>& moments)
    {
        assert(moments.size() == initial.size());
        std::vector<double> rate(moments.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            rate[i] = one_node_moments[i] - initial[i];
        }
        return rate;
    };
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    qmomkit::QmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(state, 2, 0.0, 1.0, zero_velocity, source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::post_source_inversion);
    assert(result.diagnostics.post_source_qmom_status != qmomkit::QmomStatus::success);
    assert_vector_near(result.state.moments(), state.moments(), 0.0);
    assert_vector_near(result.state.quadrature().weights, state.quadrature().weights, 0.0);
    assert_vector_near(result.state.quadrature().abscissas, state.quadrature().abscissas, 0.0);
}

void test_first_order_schemes_differ_for_noncommuting_qmom()
{
    const std::vector<double> moments = raw_moments_1d({-0.5, 1.25}, {0.4, 0.6}, 4);
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 2);
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
    const qmomkit::QmomTransportResult cd = qmomkit::advance_transport_qmom(state, 2, 0.0, 0.2, velocity, source, cd_options);

    qmomkit::QmomTransportOptions dc_options = cd_options;
    dc_options.scheme = qmomkit::TransportScheme::dc;
    const qmomkit::QmomTransportResult dc = qmomkit::advance_transport_qmom(state, 2, 0.0, 0.2, velocity, source, dc_options);

    assert(cd.ok());
    assert(dc.ok());
    assert(cd.diagnostics.scheme == qmomkit::TransportScheme::cd);
    assert(dc.diagnostics.scheme == qmomkit::TransportScheme::dc);
    assert(std::abs(cd.state.moments()[1] - dc.state.moments()[1]) > 1.0e-4);
}

void test_diagnostics_are_aggregated()
{
    const std::vector<double> moments = raw_moments_1d({1.0}, {2.0}, 2);
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    const auto velocity = [](double, const std::vector<double>& xi) { return std::vector<double>{xi[0]}; };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };

    qmomkit::QmomTransportOptions options;
    options.node_options = tight_node_options();
    options.source_options.ssp_rk3_options.estimate_error = true;
    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(state, 1, 0.0, 0.5, velocity, zero_source, options);

    assert(result.ok());
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
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 2);
    const auto velocity = [](double, const std::vector<double>&) { return std::vector<double>{0.5}; };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };

    qmomkit::QmomTransportOptions ssp_options;
    ssp_options.scheme = qmomkit::TransportScheme::cd;
    ssp_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::QmomTransportResult ssp = qmomkit::advance_transport_qmom(state, 2, 0.0, 0.2, velocity, zero_source, ssp_options);
    assert(ssp.ok());
    assert(ssp.diagnostics.first_continuous_diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    assert(ssp.diagnostics.first_continuous_diagnostics.step_count == 1);
    assert(ssp.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(ssp.diagnostics.total_node_advection_rejected_steps == 0);

    qmomkit::QmomTransportOptions rk_options = ssp_options;
    rk_options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    rk_options.node_options.adaptive_options.relative_tolerance = 1.0e-10;
    rk_options.node_options.adaptive_options.absolute_tolerance = 1.0e-12;
    rk_options.node_options.adaptive_options.initial_step = 0.01;
    const qmomkit::QmomTransportResult rk = qmomkit::advance_transport_qmom(state, 2, 0.0, 0.2, velocity, zero_source, rk_options);
    assert(rk.ok());
    assert(rk.diagnostics.first_continuous_diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(rk.diagnostics.first_continuous_diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::success);
    assert(rk.diagnostics.total_node_advection_steps > 0);

    assert_vector_near(ssp.state.moments(), rk.state.moments(), 1.0e-10);

    qmomkit::QmomTransportOptions rk45_options = rk_options;
    rk45_options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk45;
    const qmomkit::QmomTransportResult rk45 = qmomkit::advance_transport_qmom(state, 2, 0.0, 0.2, velocity, zero_source, rk45_options);
    assert(rk45.ok());
    assert(rk45.diagnostics.first_continuous_diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk45);
    assert(rk45.diagnostics.first_continuous_diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::success);
    assert_vector_near(rk.state.moments(), rk45.state.moments(), 1.0e-10);
}

void test_cqmom_state_transport_node_advection_scheme_selection()
{
    const qmomkit::CqmomCase cqmom_case = fixed_cqmom_case();
    qmomkit::CqmomState state = cqmom_case.initial_state();
    const auto velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.25); };
    const auto zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    qmomkit::CqmomTransportOptions ssp_options;
    ssp_options.scheme = qmomkit::TransportScheme::cd;
    ssp_options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::CqmomTransportResult ssp = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, velocity, zero_source, ssp_options);
    assert(ssp.ok());
    assert(ssp.quadrature.has_value());
    assert(ssp.diagnostics.first_continuous_diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    assert(ssp.diagnostics.first_continuous_diagnostics.step_count == 1);
    assert(ssp.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::success);

    qmomkit::CqmomTransportOptions rk_options = ssp_options;
    rk_options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk23;
    rk_options.node_options.adaptive_options.relative_tolerance = 1.0e-10;
    rk_options.node_options.adaptive_options.absolute_tolerance = 1.0e-12;
    rk_options.node_options.adaptive_options.initial_step = 0.01;
    const qmomkit::CqmomTransportResult rk = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, velocity, zero_source, rk_options);
    assert(rk.ok());
    assert(rk.quadrature.has_value());
    assert(rk.diagnostics.first_continuous_diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(rk.diagnostics.first_continuous_diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::success);

    qmomkit::CqmomTransportOptions rk45_options = rk_options;
    rk45_options.node_options.scheme = qmomkit::NodeAdvectionScheme::rk45;
    const qmomkit::CqmomTransportResult rk45 = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, velocity, zero_source, rk45_options);
    assert(rk45.ok());
    assert(rk45.quadrature.has_value());
    assert(rk45.diagnostics.first_continuous_diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk45);
    assert(rk45.diagnostics.first_continuous_diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::success);
    assert_cqmom_quadrature_near(*rk.quadrature, *rk45.quadrature, 1.0e-10);
}

void test_failure_propagation_from_continuous_operator()
{
    const std::vector<double> moments{1.0, 2.0};
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    const auto bad_velocity = [](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0}; };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };
    qmomkit::QmomTransportOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(state, 1, 0.0, 0.1, bad_velocity, zero_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::continuous_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::first_continuous_half_step);
    assert(result.diagnostics.first_continuous_status == qmomkit::NodeAdvectionStatus::velocity_failed);
    assert_vector_near(result.state.moments(), state.moments(), 0.0);
}

void test_failure_propagation_from_source_operator()
{
    const std::vector<double> moments{1.0, 2.0};
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto bad_source = [](double, const std::vector<double>&) { return std::vector<double>{1.0, 2.0, 3.0}; };
    qmomkit::QmomTransportOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(state, 1, 0.0, 0.1, zero_velocity, bad_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert_vector_near(result.state.moments(), state.moments(), 0.0);
}

void test_failure_propagation_from_boundary_contribution()
{
    const std::vector<double> moments{1.0, 2.0};
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };
    const qmomkit::MomentBoundaryContribution bad_boundary = [](double, const std::vector<double>& state)
    { return qmomkit::MomentBoundaryRate{std::vector<double>(state.size() + 1, 0.0)}; };
    qmomkit::QmomTransportOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::QmomTransportResult result =
        qmomkit::advance_transport_qmom(state, 1, 0.0, 0.1, zero_velocity, zero_source, bad_boundary, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
    assert(result.message.find("boundary contribution") != std::string::npos);
    assert_vector_near(result.state.moments(), state.moments(), 0.0);
}

void test_boundary_contribution_post_source_inversion_failure()
{
    const std::vector<double> moments{1.0, 0.0};
    const qmomkit::QmomState state = qmom_state_from_moments(moments, 1);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };
    const qmomkit::MomentBoundaryContribution removal = [](double, const std::vector<double>& state)
    {
        std::vector<double> rate(state.size(), 0.0);
        rate.front() = -20.0;
        return qmomkit::MomentBoundaryRate{rate};
    };
    qmomkit::QmomTransportOptions options;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(state, 1, 0.0, 0.1, zero_velocity, zero_source, removal, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::post_source_inversion);
    assert(result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::success);
    assert(result.diagnostics.post_source_qmom_status != qmomkit::QmomStatus::success);
    assert_vector_near(result.state.moments(), state.moments(), 0.0);
}

void test_context_invalid_qmom_state_is_rejected_before_transport()
{
    const std::vector<double> moments{1.0, 0.0, -1.0, 0.0};
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto zero_source = [](double, const std::vector<double>& state) { return std::vector<double>(state.size(), 0.0); };
    const qmomkit::QmomState invalid_state = qmomkit::QmomCase{qmomkit::MomentSet(raw_moments_1d({0.0}, {1.0}, moments.size())), 1}.initial_state();

    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(invalid_state, 2, 0.0, 0.1, zero_velocity, zero_source);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::invalid_input);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::none);
    assert(result.message.find("configured node count") != std::string::npos);
    assert_vector_near(result.state.moments(), invalid_state.moments(), 0.0);
}

void test_lean_cqmom_failure_propagation_from_inversion()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::TestLeanCqmomTransportState state{qmomkit::CqmomMomentData{layout, std::vector<double>(layout.required_indices.size(), 0.0)}};
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    const qmomkit::TestLeanCqmomTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::inversion_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::initial_inversion);
    assert(result.diagnostics.initial_cqmom_status != qmomkit::CqmomStatus::success);
}

void test_lean_cqmom_transport_preserves_invalid_state_messages()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    qmomkit::TestLeanCqmomTransportState missing_values = lean_state_from_dense(moments, hierarchy);
    missing_values.moments.values.pop_back();
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };

    const qmomkit::TestLeanCqmomTransportResult missing = qmomkit::advance_transport_cqmom(missing_values, 0.0, 0.1, zero_velocity);
    assert(!missing.ok());
    assert(missing.status == qmomkit::TransportStatus::invalid_input);
    assert(missing.message == "CQMOM transport moments must match the fixed layout size.");
    assert(missing.diagnostics.message == missing.message);

    qmomkit::TestLeanCqmomTransportState wrong_layout = lean_state_from_dense(moments, hierarchy);
    wrong_layout.moments.layout.required_indices.pop_back();
    const qmomkit::TestLeanCqmomTransportResult layout = qmomkit::advance_transport_cqmom(wrong_layout, 0.0, 0.1, zero_velocity);
    assert(!layout.ok());
    assert(layout.status == qmomkit::TransportStatus::invalid_input);
    assert(layout.message == "CQMOM transport moments must match the fixed layout size.");
    assert(layout.diagnostics.message == layout.message);
}

void test_lean_cqmom_failure_propagation_from_reconstruction()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestLeanCqmomTransportState state = lean_state_from_dense(moments, hierarchy);
    const auto huge_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 1.0e200); };

    qmomkit::TestLeanCqmomTransportOptions options;
    options.node_options.adaptive_options.initial_step = 1.0;
    options.node_options.adaptive_options.max_step = 1.0;
    options.scheme = qmomkit::TransportScheme::cd;
    const qmomkit::TestLeanCqmomTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, 1.0, huge_velocity, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::moment_reconstruction_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::first_moment_reconstruction);
    assert(!result.diagnostics.first_reconstruction_diagnostics.success);
}

void test_lean_cqmom_failure_propagation_from_source_operator()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestLeanCqmomTransportState state = lean_state_from_dense(moments, hierarchy);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto bad_source = [](double, const qmomkit::CqmomSourceContext& moments)
    { return qmomkit::CqmomMomentRate{std::vector<double>(moments.size() + 1, 0.0)}; };

    qmomkit::TestLeanCqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::TestLeanCqmomTransportResult result = qmomkit::advance_transport_cqmom(state, 0.0, 0.1, zero_velocity, bad_source, options);

    assert(!result.ok());
    assert(result.status == qmomkit::TransportStatus::source_failed);
    assert(result.diagnostics.failed_stage == qmomkit::TransportStage::discontinuous_source_step);
    assert(result.diagnostics.source_status == qmomkit::MomentSourceUpdateStatus::source_failed);
}

void test_lean_cqmom_failure_propagation_from_post_source_inversion()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::RawMomentTensor moments = fixed_cqmom_moments();
    const qmomkit::TestLeanCqmomTransportState state = lean_state_from_dense(moments, hierarchy);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const auto bad_realizability_source = [](double, const qmomkit::CqmomSourceContext& moments)
    {
        std::vector<double> rate(moments.size(), 0.0);
        rate.front() = -20.0;
        return qmomkit::CqmomMomentRate{rate};
    };

    qmomkit::TestLeanCqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    options.node_options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;
    const qmomkit::TestLeanCqmomTransportResult result =
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
    test_public_qmom_continuous_stage_advances_canonical_state();
    test_public_qmom_continuous_stage_failure_rolls_back();
    test_public_qmom_discontinuous_stage_fixed_mode_is_strict();
    test_public_qmom_discontinuous_stage_adaptive_policy();
    test_public_qmom_discontinuous_stage_failures_roll_back();
    test_public_qmom_discontinuous_stage_boundary_matches_source_update();
    test_public_qmom_stage_composition_matches_strang_transport();
    test_public_qmom_stage_composition_matches_lie_transport();
    test_qmom_transport_failure_after_first_stage_rolls_back_to_macro_input();
    test_cqmom_initial_state_stores_accepted_quadrature_fixed_dense_and_lean();
    test_cqmom_initial_state_stores_accepted_quadrature_adaptive_dense_and_lean();
    test_cqmom_first_continuous_stage_uses_initial_accepted_quadrature();
    test_public_cqmom_continuous_stage_preserves_topology();
    test_public_cqmom_continuous_stage_reconstructs_dense_and_lean_moments();
    test_public_cqmom_discontinuous_stage_adaptive_policy();
    test_public_cqmom_discontinuous_stage_growth_is_limited_to_one_node();
    test_lean_cqmom_discontinuous_stage_canonicalizes_adaptive_reduction();
    test_public_cqmom_stage_composition_matches_transport();
    test_public_cqmom_stage_failures_roll_back();
    test_cqmom_transport_failure_after_first_stage_rolls_back_to_macro_input();
    test_public_cqmom_boundary_discontinuous_stage_matches_transport();
    test_public_cqmom_stage_irrelevant_options_are_ignored();
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
    test_adaptive_state_source_growth_forces_adaptive_wheeler_without_mutating_qmom_options();
    test_adaptive_state_source_growth_canonicalizes_between_steps();
    test_adaptive_state_failed_growth_trial_keeps_current_count();
    test_adaptive_state_bounded_support_rejects_post_source_growth();
    test_dense_and_lean_bounded_post_source_reduction_are_equivalent();
    test_adaptive_state_reduction_prevents_source_growth_trial();
    test_adaptive_state_parent_growth_creates_unconstrained_descendants();
    test_adaptive_state_preserved_source_invalidates_descendants_after_root_collapse();
    test_adaptive_state_preserved_source_invalidates_grandchildren_after_second_level_collapse();
    test_adaptive_state_initial_counts_preserve_nonidentity_conditioning_order();
    test_fixed_state_transport_does_not_store_adaptive_branch_counts();
    test_internal_cqmom_state_backend_accessors_are_explicit_and_guarded();
    test_adaptive_state_rejects_duplicate_active_count_metadata();
    test_adaptive_state_rejects_missing_active_count_metadata();
    test_adaptive_state_failure_preserves_active_counts_for_retry();
    test_adaptive_state_post_source_failure_preserves_active_counts_for_retry();
    test_dense_and_lean_adaptive_state_active_count_persistence_are_equivalent();
    test_adaptive_state_strang_transport_matches_dense();
    test_adaptive_state_transport_reports_initial_failure();
    test_lean_cqmom_strang_matches_dense_layout_closed_source();
    test_lean_cqmom_cd_matches_dense_layout_closed_source();
    test_lean_cqmom_dc_matches_dense_layout_closed_source();
    test_explicit_dense_and_lean_adaptive_branch_options_are_equivalent();
    test_explicit_dense_and_lean_diagnostic_level_none_is_propagated();
    test_explicit_dense_and_lean_adaptive_options_are_equivalent();
    test_explicit_dense_and_lean_source_options_are_equivalent();
    test_explicit_dense_and_lean_accept_same_cqmom_native_source();
    test_explicit_dense_and_lean_accept_same_cqmom_native_boundary();
    test_explicit_dense_and_lean_cqmom_native_source_failures_are_equivalent();
    test_explicit_dense_and_lean_cqmom_native_source_invalid_values_are_equivalent();
    test_explicit_dense_and_lean_source_growth_is_automatic();
    test_state_cqmom_transport_matches_lean_fixed_layout();
    test_configured_dense_and_lean_state_transport_are_equivalent();
    test_dense_and_lean_state_source_contexts_are_identical();
    test_dense_state_reconstructs_dense_only_moments_after_discontinuous_source();
    test_dense_state_zero_source_keeps_full_dense_closure_consistent();
    test_dense_state_reconstructs_dense_only_moments_after_source_and_boundary();
    test_dense_state_reconstructs_dense_only_moments_after_repeated_source_steps();
    test_dense_state_reconstructs_dense_only_moments_with_nonidentity_conditioning_order();
    test_state_cqmom_transport_preserves_vandermonde_solver_policy();
    test_cqmom_transport_vandermonde_policy_independent_of_diagnostic_level();
    test_state_cqmom_source_evolution_matches_lean();
    test_state_cqmom_failure_propagation_from_source_operator();
    test_state_cqmom_transport_api_does_not_expose_moments();
    test_state_qmom_transport_api_does_not_expose_moments();
    test_qmom_state_copy_move_workflow();
    test_qmom_initial_state_canonicalizes_adaptive_moments();
    test_qmom_adaptive_transport_preserves_active_count_with_zero_source();
    test_qmom_adaptive_post_source_reduction_and_growth_policy();
    test_qmom_fixed_post_source_requires_configured_count();
    test_first_order_schemes_differ_for_noncommuting_qmom();
    test_diagnostics_are_aggregated();
    test_qmom_transport_node_advection_scheme_selection();
    test_cqmom_state_transport_node_advection_scheme_selection();
    test_failure_propagation_from_continuous_operator();
    test_failure_propagation_from_source_operator();
    test_failure_propagation_from_boundary_contribution();
    test_boundary_contribution_post_source_inversion_failure();
    test_context_invalid_qmom_state_is_rejected_before_transport();
    test_lean_cqmom_failure_propagation_from_inversion();
    test_lean_cqmom_transport_preserves_invalid_state_messages();
    test_lean_cqmom_failure_propagation_from_reconstruction();
    test_lean_cqmom_failure_propagation_from_source_operator();
    test_lean_cqmom_failure_propagation_from_post_source_inversion();
    return 0;
}
