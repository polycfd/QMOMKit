#include "internal/cqmom_case.hpp"
#include "internal/moment_source.hpp"
#include "internal/performance.hpp"
#include "internal/vandermonde.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/qmomkit.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

double raw_moment_1d(const std::vector<double>& nodes, const std::vector<double>& weights, const std::size_t order)
{
    double value = 0.0;
    for (std::size_t node = 0; node < nodes.size(); ++node)
    {
        value += weights[node] * std::pow(nodes[node], static_cast<int>(order));
    }
    return value;
}

qmomkit::CqmomMomentData independent_product_moments(const qmomkit::CqmomMomentLayout& layout)
{
    const std::vector<std::vector<double>> nodes_by_dimension{{-0.5, 1.0}, {0.25, 1.5}, {-1.25, 0.75}};
    const std::vector<std::vector<double>> weights_by_dimension{{0.4, 0.6}, {0.7, 0.3}, {0.2, 0.8}};
    assert(layout.hierarchy.number_of_dimensions() <= nodes_by_dimension.size());

    std::vector<double> values;
    values.reserve(layout.required_indices.size());
    for (const qmomkit::MultiIndex& index : layout.required_indices)
    {
        double value = 1.0;
        for (std::size_t dimension = 0; dimension < layout.hierarchy.number_of_dimensions(); ++dimension)
        {
            value *= raw_moment_1d(nodes_by_dimension[dimension], weights_by_dimension[dimension], index[dimension]);
        }
        values.push_back(value);
    }
    return qmomkit::CqmomMomentData{layout, std::move(values)};
}

void assert_near(const double actual, const double expected)
{
    assert(std::abs(actual - expected) < 1.0e-12);
}

void assert_vector_near(const std::vector<double>& actual, const std::vector<double>& expected)
{
    assert(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        assert_near(actual[i], expected[i]);
    }
}

qmomkit::RawMomentTensor independent_product_dense_moments(const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    const std::vector<std::vector<double>> nodes_by_dimension{{-0.5, 1.0}, {0.25, 1.5}, {-1.25, 0.75}};
    const std::vector<std::vector<double>> weights_by_dimension{{0.4, 0.6}, {0.7, 0.3}, {0.2, 0.8}};
    const std::vector<std::size_t> extents = hierarchy.moment_extents();
    std::vector<double> values(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_index = 0; flat_index < values.size(); ++flat_index)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat_index, extents);
        double value = 1.0;
        for (std::size_t dimension = 0; dimension < hierarchy.number_of_dimensions(); ++dimension)
        {
            value *= raw_moment_1d(nodes_by_dimension[dimension], weights_by_dimension[dimension], index[dimension]);
        }
        values[flat_index] = value;
    }
    return qmomkit::RawMomentTensor{extents, std::move(values)};
}

void assert_same_cqmom_result(const qmomkit::CqmomResult& actual, const qmomkit::CqmomResult& expected)
{
    assert(actual.status == expected.status);
    assert(actual.ok() == expected.ok());
    assert(actual.quadrature.weights().size() == expected.quadrature.weights().size());
    assert(actual.quadrature.abscissas().size() == expected.quadrature.abscissas().size());
    for (std::size_t node = 0; node < actual.quadrature.weights().size(); ++node)
    {
        assert_near(actual.quadrature.weights()[node], expected.quadrature.weights()[node]);
    }
    for (std::size_t dimension = 0; dimension < actual.quadrature.abscissas().size(); ++dimension)
    {
        assert(actual.quadrature.abscissas()[dimension].size() == expected.quadrature.abscissas()[dimension].size());
        for (std::size_t node = 0; node < actual.quadrature.abscissas()[dimension].size(); ++node)
        {
            assert_near(actual.quadrature.abscissas()[dimension][node], expected.quadrature.abscissas()[dimension][node]);
        }
    }
}

void assert_active_count_counters(const qmomkit::internal::performance::PerformanceCounters& counters, const std::size_t expected_derivations,
                                  const std::size_t expected_validations, const std::size_t expected_transition_builds)
{
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.cqmom_topology.active_count_derivations == expected_derivations);
    assert(counters.cqmom_topology.active_count_validations == expected_validations);
    assert(counters.cqmom_topology.active_count_transition_builds == expected_transition_builds);
#else
    (void)expected_derivations;
    (void)expected_validations;
    (void)expected_transition_builds;
    assert(counters.cqmom_topology.active_count_derivations == 0);
    assert(counters.cqmom_topology.active_count_validations == 0);
    assert(counters.cqmom_topology.active_count_transition_builds == 0);
#endif
}

qmomkit::RawMomentTensor adaptive_ragged_dense_moments()
{
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.2, 0.3, 0.5};
    const std::vector<std::vector<double>> abscissas{{-1.0, 1.0, 1.0}, {0.0, -0.5, 0.5}};
    std::vector<double> values(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_index = 0; flat_index < values.size(); ++flat_index)
    {
        const qmomkit::MultiIndex index = qmomkit::unflatten_multi_index(flat_index, extents);
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][node], static_cast<int>(index[dimension]));
            }
            values[flat_index] += contribution;
        }
    }
    return {extents, std::move(values)};
}

qmomkit::CqmomState adaptive_ragged_lean_state()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = true;
    qmom_options.rmin = 1.0e-10;
    qmom_options.eabs = 1.0e-12;
    qmom_options.cutoff = 1.0e-12;
    const qmomkit::CqmomLeanCase lean_case{qmomkit::extract_cqmom_moments(adaptive_ragged_dense_moments(), layout), {}, qmom_options, true};
    return qmomkit::internal::CqmomConstructionAccess::create_case(lean_case).initial_state();
}

qmomkit::CqmomState adaptive_ragged_dense_state()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    qmomkit::QmomOptions qmom_options;
    qmom_options.adaptive = true;
    qmom_options.rmin = 1.0e-10;
    qmom_options.eabs = 1.0e-12;
    qmom_options.cutoff = 1.0e-12;
    const qmomkit::CqmomDenseCase dense_case{hierarchy, adaptive_ragged_dense_moments(), {}, qmom_options, true};
    return qmomkit::internal::CqmomConstructionAccess::create_case(dense_case).initial_state();
}

qmomkit::CqmomState adaptive_ragged_state_without_cache(const qmomkit::CqmomState& state)
{
    if (const qmomkit::CqmomDenseCase* dense_case = qmomkit::internal::CqmomConstructionAccess::state_dense_case(state))
    {
        return qmomkit::internal::CqmomConstructionAccess::create_state(*dense_case, state.active_branch_counts());
    }
    const qmomkit::CqmomLeanCase* lean_case = qmomkit::internal::CqmomConstructionAccess::state_lean_case(state);
    assert(lean_case != nullptr);
    return qmomkit::internal::CqmomConstructionAccess::create_state(*lean_case, state.active_branch_counts());
}

void test_instrumentation_mode_matches_build_option()
{
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(qmomkit::internal::performance::enabled());
#else
    assert(!qmomkit::internal::performance::enabled());
#endif
}

void test_counter_scope_records_direct_events_when_enabled()
{
    qmomkit::internal::performance::CounterSession session;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        qmomkit::internal::performance::record_qmom_reconstruction();
        qmomkit::internal::performance::record_wheeler_call();
        qmomkit::internal::performance::record_vandermonde_solve(qmomkit::internal::performance::VandermondeSolveKind::bjorck_pereyra);
        qmomkit::internal::performance::record_vandermonde_solve(qmomkit::internal::performance::VandermondeSolveKind::pivoted_qr);
        qmomkit::internal::performance::record_vandermonde_solve(qmomkit::internal::performance::VandermondeSolveKind::svd);
        qmomkit::internal::performance::record_active_count_derivation();
        qmomkit::internal::performance::record_active_count_validation();
        qmomkit::internal::performance::record_active_count_transition_build();
        qmomkit::internal::performance::record_embedded_rk_step(true);
        qmomkit::internal::performance::record_embedded_rk_step(false);
    }

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.reconstruction.qmom_reconstructions == 1);
    assert(counters.reconstruction.wheeler_calls == 1);
    assert(counters.linear_algebra.vandermonde_solves == 3);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 1);
    assert(counters.linear_algebra.pivoted_qr_solves == 1);
    assert(counters.linear_algebra.svd_solves == 1);
    assert_active_count_counters(counters, 1, 1, 1);
    assert(counters.time_integration.embedded_rk_accepted_steps == 1);
    assert(counters.time_integration.embedded_rk_rejected_steps == 1);

    session.reset();
    const qmomkit::internal::performance::PerformanceCounters reset_counters = session.snapshot();
    assert(reset_counters.linear_algebra.vandermonde_solves == 0);
    assert(reset_counters.linear_algebra.svd_solves == 0);
#else
    assert(counters.reconstruction.qmom_reconstructions == 0);
    assert(counters.reconstruction.wheeler_calls == 0);
    assert(counters.linear_algebra.vandermonde_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.svd_solves == 0);
    assert_active_count_counters(counters, 0, 0, 0);
    assert(counters.time_integration.embedded_rk_accepted_steps == 0);
    assert(counters.time_integration.embedded_rk_rejected_steps == 0);
#endif
}

void test_adaptive_cqmom_transport_records_active_count_work_when_enabled()
{
    const auto velocity = [](double, const std::vector<double>& coordinate)
    { return std::vector<double>{0.05 + 0.01 * coordinate[0], -0.02 + 0.03 * coordinate[1]}; };
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    const auto exercise_adaptive_storage = [&](const qmomkit::CqmomState& initialized_state)
    {
        qmomkit::CqmomTransportOptions options;
        qmomkit::internal::performance::CounterSession continuous_hit_session;
        {
            qmomkit::internal::performance::CounterScope scope(continuous_hit_session);
            const qmomkit::CqmomTransportResult result = qmomkit::advance_cqmom_continuous(initialized_state, 0.0, 0.1, velocity, options);
            assert(result.ok());
        }
        assert_active_count_counters(continuous_hit_session.snapshot(), 0, 0, 0);

        const qmomkit::CqmomState cache_miss_state = adaptive_ragged_state_without_cache(initialized_state);
        qmomkit::internal::performance::CounterSession continuous_miss_session;
        {
            qmomkit::internal::performance::CounterScope scope(continuous_miss_session);
            const qmomkit::CqmomTransportResult result = qmomkit::advance_cqmom_continuous(cache_miss_state, 0.0, 0.1, velocity, options);
            assert(result.ok());
        }
        assert_active_count_counters(continuous_miss_session.snapshot(), 1, 1, 1);

        qmomkit::internal::performance::CounterSession discontinuous_hit_session;
        {
            qmomkit::internal::performance::CounterScope scope(discontinuous_hit_session);
            const qmomkit::CqmomTransportResult result = qmomkit::advance_cqmom_discontinuous(initialized_state, 0.0, 0.1, zero_source, options);
            assert(result.ok());
        }
        assert_active_count_counters(discontinuous_hit_session.snapshot(), 1, 1, 1);

        const qmomkit::CqmomState discontinuous_cache_miss_state = adaptive_ragged_state_without_cache(initialized_state);
        qmomkit::internal::performance::CounterSession discontinuous_miss_session;
        {
            qmomkit::internal::performance::CounterScope scope(discontinuous_miss_session);
            const qmomkit::CqmomTransportResult result =
                qmomkit::advance_cqmom_discontinuous(discontinuous_cache_miss_state, 0.0, 0.1, zero_source, options);
            assert(result.ok());
        }
        assert_active_count_counters(discontinuous_miss_session.snapshot(), 2, 2, 2);

        for (const qmomkit::TransportScheme scheme : {qmomkit::TransportScheme::cd, qmomkit::TransportScheme::dc, qmomkit::TransportScheme::strang})
        {
            options.scheme = scheme;
            qmomkit::internal::performance::CounterSession hit_session;
            {
                qmomkit::internal::performance::CounterScope scope(hit_session);
                const qmomkit::CqmomTransportResult result =
                    qmomkit::advance_transport_cqmom(initialized_state, 0.0, 0.1, velocity, zero_source, options);
                assert(result.ok());
            }
            assert_active_count_counters(hit_session.snapshot(), 1, 1, 1);

            const qmomkit::CqmomState cache_miss_scheme_state = adaptive_ragged_state_without_cache(initialized_state);
            qmomkit::internal::performance::CounterSession miss_session;
            {
                qmomkit::internal::performance::CounterScope scope(miss_session);
                const qmomkit::CqmomTransportResult result =
                    qmomkit::advance_transport_cqmom(cache_miss_scheme_state, 0.0, 0.1, velocity, zero_source, options);
                assert(result.ok());
            }
            assert_active_count_counters(miss_session.snapshot(), 2, 2, 2);
        }
    };

    qmomkit::internal::performance::CounterSession lean_initialization_session;
    std::optional<qmomkit::CqmomState> lean_state;
    {
        qmomkit::internal::performance::CounterScope scope(lean_initialization_session);
        lean_state.emplace(adaptive_ragged_lean_state());
    }
    assert_active_count_counters(lean_initialization_session.snapshot(), 1, 0, 0);
    assert(lean_state.has_value());
    exercise_adaptive_storage(*lean_state);

    qmomkit::internal::performance::CounterSession dense_initialization_session;
    std::optional<qmomkit::CqmomState> dense_state;
    {
        qmomkit::internal::performance::CounterScope scope(dense_initialization_session);
        dense_state.emplace(adaptive_ragged_dense_state());
    }
    assert_active_count_counters(dense_initialization_session.snapshot(), 1, 0, 0);
    assert(dense_state.has_value());
    exercise_adaptive_storage(*dense_state);

    const qmomkit::PhaseSpaceHierarchy fixed_hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout fixed_layout = qmomkit::create_cqmom_moment_layout(fixed_hierarchy);
    const qmomkit::CqmomLeanCase fixed_case{qmomkit::extract_cqmom_moments(independent_product_dense_moments(fixed_hierarchy), fixed_layout)};
    const qmomkit::CqmomState fixed_state = qmomkit::internal::CqmomConstructionAccess::create_case(fixed_case).initial_state();
    const qmomkit::CqmomState fixed_cache_miss_state = qmomkit::internal::CqmomConstructionAccess::create_state(
        fixed_case, qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(fixed_state));
    const auto& fixed_metadata = qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(fixed_state);
    assert(fixed_metadata != nullptr);
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(fixed_cache_miss_state).get() == fixed_metadata.get());
    qmomkit::internal::performance::CounterSession fixed_session;
    {
        qmomkit::internal::performance::CounterScope scope(fixed_session);
        qmomkit::CqmomTransportOptions options;
        const qmomkit::CqmomTransportResult continuous = qmomkit::advance_cqmom_continuous(fixed_cache_miss_state, 0.0, 0.1, velocity, options);
        const qmomkit::CqmomTransportResult discontinuous = qmomkit::advance_cqmom_discontinuous(fixed_state, 0.0, 0.1, zero_source, options);
        assert(continuous.ok());
        assert(discontinuous.ok());
        assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(continuous.state).get() == fixed_metadata.get());
        assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(discontinuous.state).get() == fixed_metadata.get());
    }
    assert_active_count_counters(fixed_session.snapshot(), 0, 0, 0);
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(fixed_session.snapshot().reconstruction.cqmom_prepared_metadata_uses == 2);
    assert(fixed_session.snapshot().reconstruction.cqmom_dynamic_metadata_builds == 0);
#endif
}

qmomkit::internal::DenseMatrix make_vandermonde_rhs(const std::vector<double>& nodes, const std::vector<double>& weights,
                                                    const qmomkit::internal::DenseMatrix& solution)
{
    const qmomkit::internal::DenseMatrix system = qmomkit::internal::make_conditional_vandermonde_matrix(nodes, weights);
    qmomkit::internal::DenseMatrix rhs(system.rows(), solution.cols());
    for (std::size_t row = 0; row < system.rows(); ++row)
    {
        for (std::size_t col = 0; col < solution.cols(); ++col)
        {
            double value = 0.0;
            for (std::size_t k = 0; k < system.cols(); ++k)
            {
                value += system(row, k) * solution(k, col);
            }
            rhs(row, col) = value;
        }
    }
    return rhs;
}

void test_explicit_svd_vandermonde_records_only_svd_counter_when_enabled()
{
    const std::vector<double> nodes{-0.5, 0.25, 1.5};
    const std::vector<double> weights{0.75, -1.25, 0.5};
    qmomkit::internal::DenseMatrix solution(3, 2);
    solution(0, 0) = 1.0;
    solution(1, 0) = -0.5;
    solution(2, 0) = 2.0;
    solution(0, 1) = -3.0;
    solution(1, 1) = 0.25;
    solution(2, 1) = 1.5;

    qmomkit::internal::performance::CounterSession session;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::svd, 1.0e-14, true, 1.0e-10};
        const qmomkit::internal::VandermondeSolveResult result = qmomkit::internal::solve_conditional_vandermonde_system(
            nodes, weights, make_vandermonde_rhs(nodes, weights, solution), options, "singular");
        assert(result.success);
    }

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.linear_algebra.vandermonde_solves == 1);
    assert(counters.linear_algebra.svd_solves == 1);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.dense_vandermonde_fallbacks == 0);
    assert(counters.linear_algebra.residual_checks == 1);
#else
    assert(counters.linear_algebra.vandermonde_solves == 0);
    assert(counters.linear_algebra.svd_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.dense_vandermonde_fallbacks == 0);
    assert(counters.linear_algebra.residual_checks == 0);
#endif
}

void test_explicit_svd_vandermonde_rank_deficiency_still_records_attempt_when_enabled()
{
    const std::vector<double> nodes{1.0, 1.0};
    const std::vector<double> weights{0.5, 0.5};
    const qmomkit::internal::DenseMatrix rhs(2, 1);

    qmomkit::internal::performance::CounterSession session;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        const qmomkit::internal::VandermondeSolveOptions options{qmomkit::VandermondeSolver::svd, 1.0e-14, false, 1.0e-10};
        const qmomkit::internal::VandermondeSolveResult result =
            qmomkit::internal::solve_conditional_vandermonde_system(nodes, weights, rhs, options, "singular");
        assert(!result.success);
    }

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.linear_algebra.vandermonde_solves == 1);
    assert(counters.linear_algebra.svd_solves == 1);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.dense_vandermonde_fallbacks == 0);
    assert(counters.linear_algebra.residual_checks == 0);
#else
    assert(counters.linear_algebra.vandermonde_solves == 0);
    assert(counters.linear_algebra.svd_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.dense_linear_solves == 0);
    assert(counters.linear_algebra.dense_vandermonde_fallbacks == 0);
    assert(counters.linear_algebra.residual_checks == 0);
#endif
}

void test_qmom_workflow_records_coarse_reconstruction_events_when_enabled()
{
    qmomkit::internal::performance::CounterSession session;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        const qmomkit::QmomResult result = qmomkit::compute_qmom(qmomkit::MomentSet{1.0, 0.5, 0.5, 0.5}, 2);
        assert(result.ok());
    }

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.reconstruction.qmom_reconstructions == 1);
    assert(counters.reconstruction.wheeler_calls == 1);
#else
    assert(counters.reconstruction.qmom_reconstructions == 0);
    assert(counters.reconstruction.wheeler_calls == 0);
#endif
}

void test_cqmom_source_update_reuses_stage_layout_when_enabled()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    std::vector<double> values(layout.required_indices.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        values[i] = 1.0 + 0.1 * static_cast<double>(i);
    }
    const qmomkit::CqmomMomentData moments{layout, values};

    constexpr double lambda = -0.25;
    constexpr double time_step = 0.1;
    qmomkit::internal::performance::CounterSession session;
    std::optional<qmomkit::CqmomMomentSourceUpdateResult> result;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        result.emplace(qmomkit::advance_cqmom_moments_with_source(moments, 0.0, time_step,
                                                                  [lambda](double, const qmomkit::CqmomSourceContext& context)
                                                                  {
                                                                      std::vector<double> rate(context.size(), 0.0);
                                                                      for (std::size_t i = 0; i < rate.size(); ++i)
                                                                      {
                                                                          rate[i] = lambda * context.value(i);
                                                                      }
                                                                      return qmomkit::CqmomMomentRate{rate};
                                                                  }));
    }

    assert(result.has_value());
    assert(result->ok());
    const double z = lambda * time_step;
    const double factor = 1.0 + z + 0.5 * z * z + (z * z * z) / 6.0;
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        assert(std::abs(result->moments.values[i] - values[i] * factor) < 1.0e-14);
    }

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.reconstruction.moment_layout_builds == 1);
    assert(counters.time_integration.source_rk_stages == 3);
#else
    assert(counters.reconstruction.moment_layout_builds == 0);
    assert(counters.time_integration.source_rk_stages == 0);
#endif
}

void test_trusted_cqmom_source_update_matches_general_path_and_skips_layout_validation()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData moments = independent_product_moments(layout);
    qmomkit::MomentSourceUpdateOptions options;
    options.ssp_rk3_options.estimate_error = true;

    std::vector<double> general_times;
    const qmomkit::CqmomMomentSource general_source = [&general_times](const double time, const qmomkit::CqmomSourceContext& context)
    {
        general_times.push_back(time);
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            rate[i] = -0.25 * context.value(i);
        }
        return qmomkit::CqmomMomentRate{std::move(rate)};
    };

    qmomkit::internal::performance::CounterSession general_session;
    std::optional<qmomkit::CqmomMomentSourceUpdateResult> general_result;
    {
        qmomkit::internal::performance::CounterScope scope(general_session);
        general_result.emplace(qmomkit::advance_cqmom_moments_with_source(moments, 0.2, 0.1, general_source, options));
    }

    std::vector<double> trusted_times;
    const qmomkit::CqmomMomentSource trusted_source = [&trusted_times](const double time, const qmomkit::CqmomSourceContext& context)
    {
        trusted_times.push_back(time);
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            rate[i] = -0.25 * context.value(i);
        }
        return qmomkit::CqmomMomentRate{std::move(rate)};
    };

    qmomkit::internal::performance::CounterSession trusted_session;
    std::optional<qmomkit::CqmomMomentSourceUpdateResult> trusted_result;
    {
        qmomkit::internal::performance::CounterScope scope(trusted_session);
        trusted_result.emplace(qmomkit::internal::advance_trusted_cqmom_moments_with_source(moments, 0.2, 0.1, trusted_source, options));
    }

    assert(general_result.has_value());
    assert(trusted_result.has_value());
    assert(general_result->ok());
    assert(trusted_result->ok());
    assert(general_result->status == trusted_result->status);
    assert(general_result->message == trusted_result->message);
    assert(general_result->diagnostics.ssp_rk3_status == trusted_result->diagnostics.ssp_rk3_status);
    assert(general_result->diagnostics.retries == trusted_result->diagnostics.retries);
    assert_near(general_result->diagnostics.time_step_used, trusted_result->diagnostics.time_step_used);
    assert(general_result->diagnostics.error_estimate.has_value());
    assert(trusted_result->diagnostics.error_estimate.has_value());
    assert_near(*general_result->diagnostics.error_estimate, *trusted_result->diagnostics.error_estimate);
    assert_vector_near(general_result->moments.values, trusted_result->moments.values);
    assert(general_result->moments.layout.required_indices == trusted_result->moments.layout.required_indices);
    assert_vector_near(general_times, trusted_times);
    assert(general_times.size() == 3);

    const qmomkit::internal::performance::PerformanceCounters general_counters = general_session.snapshot();
    const qmomkit::internal::performance::PerformanceCounters trusted_counters = trusted_session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(general_counters.reconstruction.moment_layout_builds == 1);
    assert(trusted_counters.reconstruction.moment_layout_builds == 0);
    assert(general_counters.time_integration.source_rk_stages == 3);
    assert(trusted_counters.time_integration.source_rk_stages == 3);
#else
    assert(general_counters.reconstruction.moment_layout_builds == 0);
    assert(trusted_counters.reconstruction.moment_layout_builds == 0);
    assert(general_counters.time_integration.source_rk_stages == 0);
    assert(trusted_counters.time_integration.source_rk_stages == 0);
#endif
}

void test_trusted_cqmom_source_update_matches_general_boundary_failures()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData moments = independent_product_moments(layout);
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    const qmomkit::CqmomMomentBoundaryContribution throwing_boundary = [](double, const qmomkit::CqmomSourceContext&) -> qmomkit::CqmomMomentRate
    { throw std::runtime_error("boundary failure"); };
    const qmomkit::CqmomMomentSourceUpdateResult general =
        qmomkit::advance_cqmom_moments_with_source(moments, 0.0, 0.1, zero_source, throwing_boundary, {});
    const qmomkit::CqmomMomentSourceUpdateResult trusted =
        qmomkit::internal::advance_trusted_cqmom_moments_with_source(moments, 0.0, 0.1, zero_source, throwing_boundary, {});

    assert(!general.ok());
    assert(!trusted.ok());
    assert(general.status == trusted.status);
    assert(general.message == trusted.message);
    assert_vector_near(general.moments.values, moments.values);
    assert_vector_near(trusted.moments.values, moments.values);

    const qmomkit::CqmomMomentSource wrong_size_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size() + 1, 0.0)}; };
    const qmomkit::CqmomMomentSourceUpdateResult general_wrong_size =
        qmomkit::advance_cqmom_moments_with_source(moments, 0.0, 0.1, wrong_size_source, {});
    const qmomkit::CqmomMomentSourceUpdateResult trusted_wrong_size =
        qmomkit::internal::advance_trusted_cqmom_moments_with_source(moments, 0.0, 0.1, wrong_size_source, {});
    assert(!general_wrong_size.ok());
    assert(!trusted_wrong_size.ok());
    assert(general_wrong_size.status == trusted_wrong_size.status);
    assert(general_wrong_size.message == trusted_wrong_size.message);
    assert_vector_near(general_wrong_size.moments.values, moments.values);
    assert_vector_near(trusted_wrong_size.moments.values, moments.values);
}

void test_state_transport_source_update_skips_redundant_layout_validation_when_enabled()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData lean_moments = independent_product_moments(layout);
    const qmomkit::CqmomLeanCase lean_case{lean_moments};
    const qmomkit::CqmomDenseCase dense_case{hierarchy, independent_product_dense_moments(hierarchy)};
    const qmomkit::CqmomState lean_state = qmomkit::internal::CqmomConstructionAccess::create_case(lean_case).initial_state();
    const qmomkit::CqmomState dense_state = qmomkit::internal::CqmomConstructionAccess::create_case(dense_case).initial_state();
    const qmomkit::CqmomMomentSource source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    qmomkit::internal::performance::CounterSession lean_session;
    std::optional<qmomkit::CqmomTransportResult> lean_result;
    {
        qmomkit::internal::performance::CounterScope scope(lean_session);
        lean_result.emplace(qmomkit::advance_cqmom_discontinuous(lean_state, 0.0, 0.1, source));
    }

    qmomkit::internal::performance::CounterSession dense_session;
    std::optional<qmomkit::CqmomTransportResult> dense_result;
    {
        qmomkit::internal::performance::CounterScope scope(dense_session);
        dense_result.emplace(qmomkit::advance_cqmom_discontinuous(dense_state, 0.0, 0.1, source));
    }

    assert(lean_result.has_value());
    assert(dense_result.has_value());
    assert(lean_result->ok());
    assert(dense_result->ok());
    const qmomkit::internal::performance::PerformanceCounters lean_counters = lean_session.snapshot();
    const qmomkit::internal::performance::PerformanceCounters dense_counters = dense_session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(lean_counters.reconstruction.moment_layout_builds == 0);
    assert(dense_counters.reconstruction.moment_layout_builds == 0);
    assert(lean_counters.reconstruction.cqmom_prepared_metadata_uses == 1);
    assert(dense_counters.reconstruction.cqmom_prepared_metadata_uses == 1);
    assert(lean_counters.reconstruction.cqmom_dynamic_metadata_builds == 0);
    assert(dense_counters.reconstruction.cqmom_dynamic_metadata_builds == 0);
    assert(lean_counters.time_integration.source_rk_stages == 3);
    assert(dense_counters.time_integration.source_rk_stages == 3);
    assert(lean_counters.transport.post_source_inversions == 1);
    assert(dense_counters.transport.post_source_inversions == 1);
    assert(lean_counters.linear_algebra.dense_vandermonde_fallbacks == 0);
    assert(dense_counters.linear_algebra.dense_vandermonde_fallbacks == 0);
#else
    assert(lean_counters.reconstruction.moment_layout_builds == 0);
    assert(dense_counters.reconstruction.moment_layout_builds == 0);
    assert(lean_counters.reconstruction.cqmom_prepared_metadata_uses == 0);
    assert(dense_counters.reconstruction.cqmom_prepared_metadata_uses == 0);
    assert(lean_counters.time_integration.source_rk_stages == 0);
    assert(dense_counters.time_integration.source_rk_stages == 0);
#endif
}

void test_ssp_node_advection_does_not_record_source_rk_stages_when_enabled()
{
    qmomkit::Quadrature1D quadrature;
    quadrature.abscissas = {1.0, 2.0};
    quadrature.weights = {0.4, 0.6};

    qmomkit::NodeAdvectionOptions options;
    options.scheme = qmomkit::NodeAdvectionScheme::ssp_rk3;

    qmomkit::internal::performance::CounterSession session;
    std::optional<qmomkit::QmomNodeAdvectionResult> result;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        result.emplace(qmomkit::advance_nodes(
            quadrature, 0.0, 0.1, [](double, const std::vector<double>& xi) { return std::vector<double>{0.25 * xi[0]}; }, options));
    }

    assert(result.has_value());
    assert(result->ok());
    assert(result->diagnostics.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.time_integration.node_rhs_evaluations == 6);
    assert(counters.time_integration.source_rk_stages == 0);
    assert(counters.time_integration.source_rk_retries == 0);
    assert(counters.time_integration.embedded_rk_accepted_steps == 0);
    assert(counters.time_integration.embedded_rk_rejected_steps == 0);
#else
    assert(counters.time_integration.node_rhs_evaluations == 0);
    assert(counters.time_integration.source_rk_stages == 0);
    assert(counters.time_integration.source_rk_retries == 0);
    assert(counters.time_integration.embedded_rk_accepted_steps == 0);
    assert(counters.time_integration.embedded_rk_rejected_steps == 0);
#endif
}

void test_public_lean_cqmom_reconstruction_reuses_validated_global_layout_when_enabled()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData moments = independent_product_moments(layout);
    const qmomkit::CqmomLeanCase lean_case{moments};

    qmomkit::internal::performance::CounterSession session;
    std::optional<qmomkit::CqmomResult> result;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        result.emplace(qmomkit::compute_cqmom(lean_case));
    }

    assert(result.has_value());
    assert(result->ok());
    assert(result->quadrature.total_nodes() == 4);
    assert_near(result->quadrature.weight({0, 0}), 0.28);
    assert_near(result->quadrature.weight({0, 1}), 0.12);
    assert_near(result->quadrature.weight({1, 0}), 0.42);
    assert_near(result->quadrature.weight({1, 1}), 0.18);
    assert_near(result->quadrature.abscissa(0, {0, 0}), -0.5);
    assert_near(result->quadrature.abscissa(0, {1, 1}), 1.0);
    assert_near(result->quadrature.abscissa(1, {0, 0}), 0.25);
    assert_near(result->quadrature.abscissa(1, {1, 1}), 1.5);

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.reconstruction.moment_layout_builds == 3);
    assert(counters.reconstruction.cqmom_prepared_metadata_builds == 0);
    assert(counters.reconstruction.cqmom_prepared_metadata_uses == 0);
    assert(counters.reconstruction.cqmom_dynamic_metadata_builds == 1);
#else
    assert(counters.reconstruction.moment_layout_builds == 0);
    assert(counters.reconstruction.cqmom_prepared_metadata_builds == 0);
    assert(counters.reconstruction.cqmom_prepared_metadata_uses == 0);
    assert(counters.reconstruction.cqmom_dynamic_metadata_builds == 0);
#endif
}

void assert_trusted_case_reconstruction_uses_prepared_metadata(const std::vector<std::size_t>& node_counts, const std::size_t expected_wheeler_calls,
                                                               const std::size_t expected_vandermonde_solves)
{
    const qmomkit::PhaseSpaceHierarchy hierarchy(node_counts);
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData moments = independent_product_moments(layout);
    const qmomkit::CqmomLeanCase lean_case{moments};
    const qmomkit::CqmomResult public_result = qmomkit::compute_cqmom(lean_case);
    assert(public_result.ok());

    const qmomkit::CqmomCase case_ = qmomkit::internal::CqmomConstructionAccess::create_case(lean_case);
    qmomkit::internal::performance::CounterSession session;
    std::optional<qmomkit::CqmomResult> trusted_result;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        trusted_result.emplace(qmomkit::compute_cqmom(case_));
    }

    assert(trusted_result.has_value());
    assert(trusted_result->ok());
    assert_same_cqmom_result(*trusted_result, public_result);

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.reconstruction.moment_layout_builds == 0);
    assert(counters.reconstruction.cqmom_prepared_metadata_builds == 0);
    assert(counters.reconstruction.cqmom_prepared_metadata_uses == 1);
    assert(counters.reconstruction.cqmom_dynamic_metadata_builds == 0);
    assert(counters.reconstruction.wheeler_calls == expected_wheeler_calls);
    assert(counters.linear_algebra.vandermonde_solves == expected_vandermonde_solves);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.svd_solves == 0);
#else
    assert(counters.reconstruction.moment_layout_builds == 0);
    assert(counters.reconstruction.cqmom_prepared_metadata_builds == 0);
    assert(counters.reconstruction.cqmom_prepared_metadata_uses == 0);
    assert(counters.reconstruction.cqmom_dynamic_metadata_builds == 0);
    assert(counters.reconstruction.wheeler_calls == 0);
    assert(counters.linear_algebra.vandermonde_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.linear_algebra.svd_solves == 0);
#endif
}

void test_trusted_cqmom_case_reconstruction_skips_global_layout_validation_when_enabled()
{
    assert_trusted_case_reconstruction_uses_prepared_metadata({2, 2}, 3, 1);
    assert_trusted_case_reconstruction_uses_prepared_metadata({2, 2, 2}, 7, 3);
}

void test_fixed_and_dynamic_cqmom_metadata_counters()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData moments = independent_product_moments(layout);
    const qmomkit::CqmomLeanCase fixed_lean_case{moments};

    qmomkit::internal::performance::CounterSession fixed_session;
    {
        qmomkit::internal::performance::CounterScope scope(fixed_session);
        const qmomkit::CqmomCase fixed_case = qmomkit::internal::CqmomConstructionAccess::create_case(fixed_lean_case);
        assert(qmomkit::compute_cqmom(fixed_case).ok());
        assert(qmomkit::compute_cqmom(fixed_case).ok());
    }

    const qmomkit::CqmomLeanCase adaptive_lean_case{moments, {}, {true, 1.0e-10, 1.0e-12, 1.0e-12}, true};
    qmomkit::internal::performance::CounterSession adaptive_session;
    {
        qmomkit::internal::performance::CounterScope scope(adaptive_session);
        const qmomkit::CqmomCase adaptive_case = qmomkit::internal::CqmomConstructionAccess::create_case(adaptive_lean_case);
        assert(qmomkit::compute_cqmom(adaptive_case).ok());
    }

    const qmomkit::internal::performance::PerformanceCounters fixed_counters = fixed_session.snapshot();
    const qmomkit::internal::performance::PerformanceCounters adaptive_counters = adaptive_session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(fixed_counters.reconstruction.cqmom_prepared_metadata_builds == 1);
    assert(fixed_counters.reconstruction.cqmom_prepared_metadata_uses == 2);
    assert(fixed_counters.reconstruction.cqmom_dynamic_metadata_builds == 0);
    assert(adaptive_counters.reconstruction.cqmom_prepared_metadata_builds == 0);
    assert(adaptive_counters.reconstruction.cqmom_prepared_metadata_uses == 0);
    assert(adaptive_counters.reconstruction.cqmom_dynamic_metadata_builds == 1);
#else
    assert(fixed_counters.reconstruction.cqmom_prepared_metadata_builds == 0);
    assert(fixed_counters.reconstruction.cqmom_prepared_metadata_uses == 0);
    assert(fixed_counters.reconstruction.cqmom_dynamic_metadata_builds == 0);
    assert(adaptive_counters.reconstruction.cqmom_prepared_metadata_builds == 0);
    assert(adaptive_counters.reconstruction.cqmom_prepared_metadata_uses == 0);
    assert(adaptive_counters.reconstruction.cqmom_dynamic_metadata_builds == 0);
#endif
}

void test_every_fixed_state_factory_prepares_and_reconstructs_without_dynamic_metadata()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomLeanCase lean_case{independent_product_moments(layout)};
    const qmomkit::CqmomDenseCase dense_case{hierarchy, independent_product_dense_moments(hierarchy)};
    const qmomkit::CqmomState dense_initial = qmomkit::internal::CqmomConstructionAccess::create_case(dense_case).initial_state();
    const qmomkit::CqmomState lean_initial = qmomkit::internal::CqmomConstructionAccess::create_case(lean_case).initial_state();
    const qmomkit::CqmomQuadrature dense_quadrature = *qmomkit::internal::CqmomConstructionAccess::state_accepted_quadrature(dense_initial);
    const qmomkit::CqmomQuadrature lean_quadrature = *qmomkit::internal::CqmomConstructionAccess::state_accepted_quadrature(lean_initial);
    const qmomkit::CqmomMomentSource zero_source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    qmomkit::internal::performance::CounterSession session;
    {
        qmomkit::internal::performance::CounterScope scope(session);
        std::vector<qmomkit::CqmomState> states;
        states.reserve(8);
        states.push_back(qmomkit::internal::CqmomConstructionAccess::create_state(dense_case));
        states.push_back(qmomkit::internal::CqmomConstructionAccess::create_state(dense_case, std::vector<qmomkit::CqmomActiveBranchCount>{}));
        states.push_back(qmomkit::internal::CqmomConstructionAccess::create_state(dense_case, {}, std::nullopt));
        states.push_back(qmomkit::internal::CqmomConstructionAccess::create_state(dense_case, {}, qmomkit::CqmomQuadrature{dense_quadrature}));
        states.push_back(qmomkit::internal::CqmomConstructionAccess::create_state(lean_case));
        states.push_back(qmomkit::internal::CqmomConstructionAccess::create_state(lean_case, std::vector<qmomkit::CqmomActiveBranchCount>{}));
        states.push_back(qmomkit::internal::CqmomConstructionAccess::create_state(lean_case, {}, std::nullopt));
        states.push_back(qmomkit::internal::CqmomConstructionAccess::create_state(lean_case, {}, qmomkit::CqmomQuadrature{lean_quadrature}));

        for (const qmomkit::CqmomState& state : states)
        {
            assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(state) != nullptr);
            const qmomkit::CqmomTransportResult result = qmomkit::advance_cqmom_discontinuous(state, 0.0, 0.1, zero_source);
            assert(result.ok());
            assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(result.state).get() ==
                   qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(state).get());
        }
    }

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.reconstruction.cqmom_prepared_metadata_builds == 8);
    assert(counters.reconstruction.cqmom_prepared_metadata_uses == 14);
    assert(counters.reconstruction.cqmom_dynamic_metadata_builds == 0);
#else
    assert(counters.reconstruction.cqmom_prepared_metadata_builds == 0);
    assert(counters.reconstruction.cqmom_prepared_metadata_uses == 0);
    assert(counters.reconstruction.cqmom_dynamic_metadata_builds == 0);
#endif
}

} // namespace

int main()
{
    test_instrumentation_mode_matches_build_option();
    test_counter_scope_records_direct_events_when_enabled();
    test_adaptive_cqmom_transport_records_active_count_work_when_enabled();
    test_explicit_svd_vandermonde_records_only_svd_counter_when_enabled();
    test_explicit_svd_vandermonde_rank_deficiency_still_records_attempt_when_enabled();
    test_qmom_workflow_records_coarse_reconstruction_events_when_enabled();
    test_cqmom_source_update_reuses_stage_layout_when_enabled();
    test_trusted_cqmom_source_update_matches_general_path_and_skips_layout_validation();
    test_trusted_cqmom_source_update_matches_general_boundary_failures();
    test_state_transport_source_update_skips_redundant_layout_validation_when_enabled();
    test_ssp_node_advection_does_not_record_source_rk_stages_when_enabled();
    test_public_lean_cqmom_reconstruction_reuses_validated_global_layout_when_enabled();
    test_trusted_cqmom_case_reconstruction_skips_global_layout_validation_when_enabled();
    test_fixed_and_dynamic_cqmom_metadata_counters();
    test_every_fixed_state_factory_prepares_and_reconstructs_without_dynamic_metadata();

    return 0;
}
