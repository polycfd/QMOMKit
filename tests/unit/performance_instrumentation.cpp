#include "internal/cqmom_case.hpp"
#include "internal/performance.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/qmomkit.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <optional>
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
        qmomkit::internal::performance::record_rk23_step(true);
        qmomkit::internal::performance::record_rk23_step(false);
    }

    const qmomkit::internal::performance::PerformanceCounters counters = session.snapshot();
#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)
    assert(counters.reconstruction.qmom_reconstructions == 1);
    assert(counters.reconstruction.wheeler_calls == 1);
    assert(counters.linear_algebra.vandermonde_solves == 2);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 1);
    assert(counters.linear_algebra.pivoted_qr_solves == 1);
    assert(counters.time_integration.rk23_accepted_steps == 1);
    assert(counters.time_integration.rk23_rejected_steps == 1);
#else
    assert(counters.reconstruction.qmom_reconstructions == 0);
    assert(counters.reconstruction.wheeler_calls == 0);
    assert(counters.linear_algebra.vandermonde_solves == 0);
    assert(counters.linear_algebra.bjorck_pereyra_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
    assert(counters.time_integration.rk23_accepted_steps == 0);
    assert(counters.time_integration.rk23_rejected_steps == 0);
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
    assert(counters.time_integration.rk23_accepted_steps == 0);
    assert(counters.time_integration.rk23_rejected_steps == 0);
#else
    assert(counters.time_integration.node_rhs_evaluations == 0);
    assert(counters.time_integration.source_rk_stages == 0);
    assert(counters.time_integration.source_rk_retries == 0);
    assert(counters.time_integration.rk23_accepted_steps == 0);
    assert(counters.time_integration.rk23_rejected_steps == 0);
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
#else
    assert(counters.reconstruction.moment_layout_builds == 0);
#endif
}

void assert_trusted_case_reconstruction_uses_only_branch_layout_builds(const std::vector<std::size_t>& node_counts,
                                                                       const std::size_t expected_layout_builds,
                                                                       const std::size_t expected_wheeler_calls,
                                                                       const std::size_t expected_vandermonde_solves)
{
    const qmomkit::PhaseSpaceHierarchy hierarchy(node_counts);
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData moments = independent_product_moments(layout);
    const qmomkit::CqmomLeanCase lean_case{moments};
    const qmomkit::CqmomResult public_result = qmomkit::compute_cqmom(lean_case);
    assert(public_result.ok());

    const qmomkit::CqmomCase case_ = qmomkit::internal::create_cqmom_case(lean_case);
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
    assert(counters.reconstruction.moment_layout_builds == expected_layout_builds);
    assert(counters.reconstruction.wheeler_calls == expected_wheeler_calls);
    assert(counters.linear_algebra.vandermonde_solves == expected_vandermonde_solves);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
#else
    assert(counters.reconstruction.moment_layout_builds == 0);
    assert(counters.reconstruction.wheeler_calls == 0);
    assert(counters.linear_algebra.vandermonde_solves == 0);
    assert(counters.linear_algebra.pivoted_qr_solves == 0);
#endif
}

void test_trusted_cqmom_case_reconstruction_skips_global_layout_validation_when_enabled()
{
    assert_trusted_case_reconstruction_uses_only_branch_layout_builds({2, 2}, 2, 3, 1);
    assert_trusted_case_reconstruction_uses_only_branch_layout_builds({2, 2, 2}, 3, 7, 3);
}

} // namespace

int main()
{
    test_instrumentation_mode_matches_build_option();
    test_counter_scope_records_direct_events_when_enabled();
    test_qmom_workflow_records_coarse_reconstruction_events_when_enabled();
    test_cqmom_source_update_reuses_stage_layout_when_enabled();
    test_ssp_node_advection_does_not_record_source_rk_stages_when_enabled();
    test_public_lean_cqmom_reconstruction_reuses_validated_global_layout_when_enabled();
    test_trusted_cqmom_case_reconstruction_skips_global_layout_validation_when_enabled();

    return 0;
}
