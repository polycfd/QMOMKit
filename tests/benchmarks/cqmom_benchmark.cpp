#include "internal/recursive_quadrature.hpp"
#include "qmomkit/qmomkit.hpp"

#include <algorithm>
#include <benchmark/benchmark.h>
#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace
{

std::vector<double> raw_moments_from_terminal_quadrature(const std::vector<std::vector<double>>& abscissas, const std::vector<double>& weights,
                                                         const std::vector<std::size_t>& extents)
{
    std::vector<double> moments(qmomkit::product_extents(extents), 0.0);
    for (std::size_t flat_moment = 0; flat_moment < moments.size(); ++flat_moment)
    {
        const qmomkit::MultiIndex moment_index = qmomkit::unflatten_multi_index(flat_moment, extents);
        double value = 0.0;
        for (std::size_t node = 0; node < weights.size(); ++node)
        {
            double contribution = weights[node];
            for (std::size_t dimension = 0; dimension < abscissas.size(); ++dimension)
            {
                contribution *= std::pow(abscissas[dimension][node], static_cast<int>(moment_index[dimension]));
            }
            value += contribution;
        }
        moments[flat_moment] = value;
    }
    return moments;
}

struct CqmomFixture
{
    qmomkit::PhaseSpaceHierarchy hierarchy;
    std::vector<std::size_t> extents;
    std::vector<double> moments;
    qmomkit::QmomOptions qmom_options;
    bool adaptive_branch_counts = false;
};

struct BranchSolveStats
{
    std::size_t solve_count = 0;
    std::size_t bp_count = 0;
    std::size_t qr_count = 0;
    std::size_t selected_bp_count = 0;
    std::size_t selected_qr_count = 0;
    std::size_t fallback_count = 0;
    double total_node_count = 0.0;
    double total_rhs_width = 0.0;
    std::size_t max_node_count = 0;
    std::size_t max_rhs_width = 0;
};

struct CqmomRepresentationFixture
{
    CqmomFixture dense_fixture;
    qmomkit::CqmomMomentLayout layout;
    qmomkit::CqmomMomentData lean_moments;
    qmomkit::internal::IndexedMomentSet indexed_moments;
};

CqmomFixture make_conditioning_order_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> weights{0.08, 0.12, 0.56, 0.24};
    const std::vector<std::vector<double>> abscissas{{-2.0, 1.0, -0.5, 2.25}, {-1.5, -1.5, 1.0, 1.0}};
    return {hierarchy, extents, raw_moments_from_terminal_quadrature(abscissas, weights, extents), {}, false};
}

CqmomFixture make_adaptive_ragged_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2, 2});
    const std::vector<std::size_t> extents{4, 4, 4};
    const std::vector<double> weights{0.05, 0.05, 0.30, 0.60};
    const std::vector<std::vector<double>> abscissas{{-1.0, -1.0, -1.0, 1.0}, {-2.0, -2.0, 2.0, 3.0}, {-1.0, 1.0, 0.0, 2.0}};
    return {hierarchy, extents, raw_moments_from_terminal_quadrature(abscissas, weights, extents), {true, 1.0e-10, 1.0e-10, 1.0e-12}, true};
}

std::vector<double> make_separated_nodes(const std::size_t node_count, const std::size_t dimension)
{
    std::vector<double> nodes(node_count, 0.0);
    if (node_count == 1)
    {
        nodes.front() = 0.07 * static_cast<double>(dimension + 1);
        return nodes;
    }

    for (std::size_t node = 0; node < node_count; ++node)
    {
        const double fraction = static_cast<double>(node) / static_cast<double>(node_count - 1);
        nodes[node] = -0.85 + 1.70 * fraction + 0.04 * static_cast<double>(dimension);
    }
    return nodes;
}

std::vector<double> make_nonuniform_weights(const std::size_t node_count, const std::size_t dimension)
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

CqmomFixture make_independent_fixed_fixture(const std::vector<std::size_t>& node_counts)
{
    const qmomkit::PhaseSpaceHierarchy hierarchy(node_counts);
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
        moments_by_dimension.push_back(one_dimensional_raw_moments(make_separated_nodes(node_counts[dimension], dimension),
                                                                   make_nonuniform_weights(node_counts[dimension], dimension), extents[dimension]));
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

    return {hierarchy, extents, std::move(moments), {}, false};
}

CqmomFixture make_uniform_independent_fixed_fixture(const std::size_t dimensions, const std::size_t node_count)
{
    return make_independent_fixed_fixture(std::vector<std::size_t>(dimensions, node_count));
}

CqmomRepresentationFixture make_representation_fixture(const std::size_t dimensions, const std::size_t node_count)
{
    CqmomFixture fixture = make_uniform_independent_fixed_fixture(dimensions, node_count);
    qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(fixture.hierarchy);
    qmomkit::CqmomMomentData lean_moments = qmomkit::extract_cqmom_moments({fixture.extents, fixture.moments}, layout);
    qmomkit::internal::IndexedMomentSet indexed_moments = qmomkit::internal::create_indexed_cqmom_moment_set(lean_moments);
    return {std::move(fixture), std::move(layout), std::move(lean_moments), std::move(indexed_moments)};
}

qmomkit::CqmomDenseCase create_dense_case(const CqmomFixture& fixture,
                                          const std::optional<qmomkit::CqmomDiagnosticLevel> diagnostic_level = std::nullopt)
{
    qmomkit::CqmomDenseCase dense_case{
        fixture.hierarchy, {fixture.extents, fixture.moments}, {}, fixture.qmom_options, fixture.adaptive_branch_counts};
    if (diagnostic_level.has_value())
    {
        dense_case.diagnostic_level = *diagnostic_level;
    }
    return dense_case;
}

qmomkit::WheelerOptions make_wheeler_options(const qmomkit::QmomOptions& options)
{
    return {options.adaptive, options.rmin, options.eabs, options.cutoff};
}

qmomkit::internal::RecursiveQuadratureRequest
make_recursive_request(const CqmomFixture& fixture, const qmomkit::CqmomDiagnosticLevel diagnostic_level, const qmomkit::VandermondeSolver solver)
{
    return {fixture.hierarchy,
            {fixture.extents, fixture.moments},
            {},
            make_wheeler_options(fixture.qmom_options),
            fixture.adaptive_branch_counts,
            diagnostic_level,
            solver};
}

qmomkit::internal::RecursiveQuadratureMomentDataRequest make_recursive_moment_data_request(const CqmomRepresentationFixture& fixture,
                                                                                           const qmomkit::CqmomDiagnosticLevel diagnostic_level,
                                                                                           const qmomkit::VandermondeSolver solver)
{
    return {fixture.dense_fixture.hierarchy,
            fixture.lean_moments,
            {},
            make_wheeler_options(fixture.dense_fixture.qmom_options),
            fixture.dense_fixture.adaptive_branch_counts,
            diagnostic_level,
            solver};
}

qmomkit::internal::RecursiveQuadratureMomentViewRequest make_recursive_lean_kernel_request(const CqmomRepresentationFixture& fixture,
                                                                                           const qmomkit::CqmomDiagnosticLevel diagnostic_level,
                                                                                           const qmomkit::VandermondeSolver solver)
{
    return {fixture.dense_fixture.hierarchy,
            qmomkit::internal::CqmomMomentView(fixture.indexed_moments),
            {},
            make_wheeler_options(fixture.dense_fixture.qmom_options),
            fixture.dense_fixture.adaptive_branch_counts,
            diagnostic_level,
            solver};
}

std::size_t rhs_width_for_level(const CqmomFixture& fixture, const std::size_t level)
{
    std::size_t rhs_width = 1;
    const std::vector<std::size_t>& order = fixture.hierarchy.conditioning_order();
    for (std::size_t ordered_axis = level + 1; ordered_axis < order.size(); ++ordered_axis)
    {
        rhs_width *= fixture.extents[order[ordered_axis]];
    }
    return rhs_width;
}

BranchSolveStats collect_branch_solve_stats(const CqmomFixture& fixture, const qmomkit::internal::RecursiveQuadratureResult& result)
{
    BranchSolveStats stats;
    for (const qmomkit::internal::RecursiveQuadratureBranchDiagnostics& diagnostics : result.diagnostics)
    {
        if (diagnostics.conditioning_level + 1 >= fixture.hierarchy.number_of_dimensions())
        {
            continue;
        }

        ++stats.solve_count;
        stats.total_node_count += static_cast<double>(diagnostics.actual_node_count);
        const std::size_t rhs_width = rhs_width_for_level(fixture, diagnostics.conditioning_level);
        stats.total_rhs_width += static_cast<double>(rhs_width);
        stats.max_node_count = std::max(stats.max_node_count, diagnostics.actual_node_count);
        stats.max_rhs_width = std::max(stats.max_rhs_width, rhs_width);
        if (diagnostics.conditional_vandermonde_algorithm_used == qmomkit::internal::VandermondeAlgorithm::bjorck_pereyra)
        {
            ++stats.bp_count;
        }
        if (diagnostics.conditional_vandermonde_algorithm_used == qmomkit::internal::VandermondeAlgorithm::pivoted_qr)
        {
            ++stats.qr_count;
        }
        if (diagnostics.conditional_vandermonde_selected_solver == qmomkit::VandermondeSolver::bjorck_pereyra ||
            diagnostics.conditional_vandermonde_selected_solver == qmomkit::VandermondeSolver::bjorck_pereyra_with_dense_fallback)
        {
            ++stats.selected_bp_count;
        }
        if (diagnostics.conditional_vandermonde_selected_solver == qmomkit::VandermondeSolver::pivoted_qr)
        {
            ++stats.selected_qr_count;
        }
        if (diagnostics.dense_vandermonde_fallback_occurred)
        {
            ++stats.fallback_count;
        }
    }
    return stats;
}

qmomkit::CqmomDiagnosticLevel stats_diagnostic_level(const qmomkit::CqmomDiagnosticLevel diagnostic_level)
{
    return diagnostic_level == qmomkit::CqmomDiagnosticLevel::none ? qmomkit::CqmomDiagnosticLevel::basic : diagnostic_level;
}

void set_fixture_counters(benchmark::State& state, const CqmomFixture& fixture, const std::size_t terminal_nodes)
{
    double total_requested_nodes = 0.0;
    std::size_t max_requested_nodes = 0;
    for (const std::size_t node_count : fixture.hierarchy.node_counts())
    {
        total_requested_nodes += static_cast<double>(node_count);
        max_requested_nodes = std::max(max_requested_nodes, node_count);
    }

    const double dimensions = static_cast<double>(fixture.hierarchy.number_of_dimensions());
    state.counters["dimensions"] = dimensions;
    state.counters["requested_terminal_nodes"] = static_cast<double>(fixture.hierarchy.total_nodes());
    state.counters["terminal_nodes"] = static_cast<double>(terminal_nodes);
    state.counters["moment_count"] = static_cast<double>(fixture.moments.size());
    state.counters["avg_requested_nodes"] = total_requested_nodes / dimensions;
    state.counters["max_requested_nodes"] = static_cast<double>(max_requested_nodes);
}

void set_representation_counters(benchmark::State& state, const CqmomRepresentationFixture& fixture)
{
    const double dense_count = static_cast<double>(fixture.dense_fixture.moments.size());
    const double lean_count = static_cast<double>(fixture.layout.required_indices.size());
    state.counters["dense_moment_count"] = dense_count;
    state.counters["lean_moment_count"] = lean_count;
    state.counters["compression_ratio"] = dense_count / lean_count;
}

void set_branch_solve_counters(benchmark::State& state, const BranchSolveStats& stats)
{
    if (stats.solve_count == 0)
    {
        return;
    }

    const double solve_count = static_cast<double>(stats.solve_count);
    state.counters["branch_solves"] = solve_count;
    state.counters["bp_fraction"] = static_cast<double>(stats.bp_count) / solve_count;
    state.counters["qr_fraction"] = static_cast<double>(stats.qr_count) / solve_count;
    state.counters["selected_bp_fraction"] = static_cast<double>(stats.selected_bp_count) / solve_count;
    state.counters["selected_qr_fraction"] = static_cast<double>(stats.selected_qr_count) / solve_count;
    state.counters["dense_fallback_fraction"] = static_cast<double>(stats.fallback_count) / solve_count;
    state.counters["avg_branch_nodes"] = stats.total_node_count / solve_count;
    state.counters["avg_rhs_width"] = stats.total_rhs_width / solve_count;
    state.counters["max_branch_nodes"] = static_cast<double>(stats.max_node_count);
    state.counters["max_rhs_width"] = static_cast<double>(stats.max_rhs_width);
}

void run_cqmom_benchmark(benchmark::State& state, const CqmomFixture& fixture,
                         const std::optional<qmomkit::CqmomDiagnosticLevel> diagnostic_level = std::nullopt)
{
    const qmomkit::CqmomDenseCase dense_case = create_dense_case(fixture, diagnostic_level);
    const qmomkit::CqmomResult setup_check = qmomkit::compute_cqmom(dense_case);
    if (!setup_check.ok())
    {
        state.SkipWithError("deterministic CQMOM fixture failed setup");
        return;
    }
    set_fixture_counters(state, fixture, setup_check.quadrature.total_nodes());

    const qmomkit::internal::RecursiveQuadratureRequest stats_request =
        make_recursive_request(fixture, stats_diagnostic_level(dense_case.diagnostic_level), qmomkit::VandermondeSolver::automatic);
    const qmomkit::internal::RecursiveQuadratureResult stats_check = qmomkit::internal::compute_recursive_quadrature(stats_request);
    if (stats_check.ok())
    {
        set_branch_solve_counters(state, collect_branch_solve_stats(fixture, stats_check));
    }

    for (auto _ : state)
    {
        const qmomkit::CqmomResult result = qmomkit::compute_cqmom(dense_case);
        benchmark::DoNotOptimize(result.quadrature.weights().data());
        benchmark::DoNotOptimize(result.quadrature.abscissas().data());
    }
}

void run_recursive_cqmom_benchmark(benchmark::State& state, const CqmomFixture& fixture, const qmomkit::CqmomDiagnosticLevel diagnostic_level,
                                   const qmomkit::VandermondeSolver solver)
{
    const qmomkit::internal::RecursiveQuadratureRequest request = make_recursive_request(fixture, diagnostic_level, solver);
    const qmomkit::internal::RecursiveQuadratureResult setup_check = qmomkit::internal::compute_recursive_quadrature(request);
    if (!setup_check.ok())
    {
        state.SkipWithError("recursive CQMOM fixture failed setup");
        return;
    }
    set_fixture_counters(state, fixture, setup_check.quadrature.total_nodes());

    const qmomkit::internal::RecursiveQuadratureRequest stats_request =
        make_recursive_request(fixture, stats_diagnostic_level(diagnostic_level), solver);
    const qmomkit::internal::RecursiveQuadratureResult stats_check = qmomkit::internal::compute_recursive_quadrature(stats_request);
    if (stats_check.ok())
    {
        set_branch_solve_counters(state, collect_branch_solve_stats(fixture, stats_check));
    }

    for (auto _ : state)
    {
        const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature(request);
        benchmark::DoNotOptimize(result.quadrature.weights().data());
        benchmark::DoNotOptimize(result.quadrature.abscissas().data());
    }
}

void run_recursive_cqmom_application_benchmark(benchmark::State& state, const CqmomFixture& fixture)
{
    run_recursive_cqmom_benchmark(state, fixture, qmomkit::CqmomDiagnosticLevel::none, qmomkit::VandermondeSolver::automatic);
}

void run_cqmom_public_dense_representation_benchmark(benchmark::State& state, const std::size_t dimensions, const std::size_t node_count)
{
    const CqmomRepresentationFixture fixture = make_representation_fixture(dimensions, node_count);
    qmomkit::CqmomDenseCase dense_case = create_dense_case(fixture.dense_fixture, qmomkit::CqmomDiagnosticLevel::none);
    const qmomkit::CqmomResult setup_check = qmomkit::compute_cqmom(dense_case);
    if (!setup_check.ok())
    {
        state.SkipWithError("dense public CQMOM representation fixture failed setup");
        return;
    }
    set_fixture_counters(state, fixture.dense_fixture, setup_check.quadrature.total_nodes());
    set_representation_counters(state, fixture);

    for (auto _ : state)
    {
        const qmomkit::CqmomResult result = qmomkit::compute_cqmom(dense_case);
        benchmark::DoNotOptimize(result.quadrature.weights().data());
        benchmark::DoNotOptimize(result.quadrature.abscissas().data());
    }
}

void run_cqmom_public_lean_representation_benchmark(benchmark::State& state, const std::size_t dimensions, const std::size_t node_count)
{
    const CqmomRepresentationFixture fixture = make_representation_fixture(dimensions, node_count);
    qmomkit::CqmomLeanCase lean_case{fixture.lean_moments,
                                     {},
                                     fixture.dense_fixture.qmom_options,
                                     fixture.dense_fixture.adaptive_branch_counts,
                                     qmomkit::CqmomDiagnosticLevel::none};
    const qmomkit::CqmomResult setup_check = qmomkit::compute_cqmom(lean_case);
    if (!setup_check.ok())
    {
        state.SkipWithError("lean public CQMOM representation fixture failed setup");
        return;
    }
    set_fixture_counters(state, fixture.dense_fixture, setup_check.quadrature.total_nodes());
    set_representation_counters(state, fixture);

    for (auto _ : state)
    {
        const qmomkit::CqmomResult result = qmomkit::compute_cqmom(lean_case);
        benchmark::DoNotOptimize(result.quadrature.weights().data());
        benchmark::DoNotOptimize(result.quadrature.abscissas().data());
    }
}

void run_cqmom_recursive_dense_representation_benchmark(benchmark::State& state, const std::size_t dimensions, const std::size_t node_count)
{
    const CqmomRepresentationFixture fixture = make_representation_fixture(dimensions, node_count);
    const qmomkit::internal::RecursiveQuadratureRequest request =
        make_recursive_request(fixture.dense_fixture, qmomkit::CqmomDiagnosticLevel::none, qmomkit::VandermondeSolver::automatic);
    const qmomkit::internal::RecursiveQuadratureResult setup_check = qmomkit::internal::compute_recursive_quadrature(request);
    if (!setup_check.ok())
    {
        state.SkipWithError("dense recursive CQMOM representation fixture failed setup");
        return;
    }
    set_fixture_counters(state, fixture.dense_fixture, setup_check.quadrature.total_nodes());
    set_representation_counters(state, fixture);

    for (auto _ : state)
    {
        const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature(request);
        benchmark::DoNotOptimize(result.quadrature.weights().data());
        benchmark::DoNotOptimize(result.quadrature.abscissas().data());
    }
}

void run_cqmom_recursive_lean_representation_benchmark(benchmark::State& state, const std::size_t dimensions, const std::size_t node_count)
{
    const CqmomRepresentationFixture fixture = make_representation_fixture(dimensions, node_count);
    const qmomkit::internal::RecursiveQuadratureMomentDataRequest moment_data_request =
        make_recursive_moment_data_request(fixture, qmomkit::CqmomDiagnosticLevel::none, qmomkit::VandermondeSolver::automatic);
    const qmomkit::internal::RecursiveQuadratureResult setup_check =
        qmomkit::internal::compute_recursive_quadrature_from_moment_data(moment_data_request);
    if (!setup_check.ok())
    {
        state.SkipWithError("lean recursive CQMOM representation fixture failed setup");
        return;
    }
    set_fixture_counters(state, fixture.dense_fixture, setup_check.quadrature.total_nodes());
    set_representation_counters(state, fixture);

    for (auto _ : state)
    {
        const qmomkit::internal::RecursiveQuadratureResult result =
            qmomkit::internal::compute_recursive_quadrature_from_moment_data(moment_data_request);
        benchmark::DoNotOptimize(result.quadrature.weights().data());
        benchmark::DoNotOptimize(result.quadrature.abscissas().data());
    }
}

void run_cqmom_kernel_dense_representation_benchmark(benchmark::State& state, const std::size_t dimensions, const std::size_t node_count)
{
    const CqmomRepresentationFixture fixture = make_representation_fixture(dimensions, node_count);
    const qmomkit::RawMomentTensor dense_moments{fixture.dense_fixture.extents, fixture.dense_fixture.moments};
    const qmomkit::internal::RecursiveQuadratureMomentViewRequest request{fixture.dense_fixture.hierarchy,
                                                                          qmomkit::internal::CqmomMomentView(dense_moments),
                                                                          {},
                                                                          make_wheeler_options(fixture.dense_fixture.qmom_options),
                                                                          fixture.dense_fixture.adaptive_branch_counts,
                                                                          qmomkit::CqmomDiagnosticLevel::none,
                                                                          qmomkit::VandermondeSolver::automatic};
    const qmomkit::internal::RecursiveQuadratureResult setup_check = qmomkit::internal::compute_recursive_quadrature_from_moment_view(request);
    if (!setup_check.ok())
    {
        state.SkipWithError("dense recursive CQMOM kernel fixture failed setup");
        return;
    }
    set_fixture_counters(state, fixture.dense_fixture, setup_check.quadrature.total_nodes());
    set_representation_counters(state, fixture);

    for (auto _ : state)
    {
        const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_view(request);
        benchmark::DoNotOptimize(result.quadrature.weights().data());
        benchmark::DoNotOptimize(result.quadrature.abscissas().data());
    }
}

void run_cqmom_kernel_lean_representation_benchmark(benchmark::State& state, const std::size_t dimensions, const std::size_t node_count)
{
    const CqmomRepresentationFixture fixture = make_representation_fixture(dimensions, node_count);
    const qmomkit::internal::RecursiveQuadratureMomentViewRequest request =
        make_recursive_lean_kernel_request(fixture, qmomkit::CqmomDiagnosticLevel::none, qmomkit::VandermondeSolver::automatic);
    const qmomkit::internal::RecursiveQuadratureResult setup_check = qmomkit::internal::compute_recursive_quadrature_from_moment_view(request);
    if (!setup_check.ok())
    {
        state.SkipWithError("lean recursive CQMOM kernel fixture failed setup");
        return;
    }
    set_fixture_counters(state, fixture.dense_fixture, setup_check.quadrature.total_nodes());
    set_representation_counters(state, fixture);

    for (auto _ : state)
    {
        const qmomkit::internal::RecursiveQuadratureResult result = qmomkit::internal::compute_recursive_quadrature_from_moment_view(request);
        benchmark::DoNotOptimize(result.quadrature.weights().data());
        benchmark::DoNotOptimize(result.quadrature.abscissas().data());
    }
}

void bm_cqmom_2d_conditioning_order_baseline(benchmark::State& state)
{
    run_cqmom_benchmark(state, make_conditioning_order_fixture(), qmomkit::CqmomDiagnosticLevel::none);
}

void bm_cqmom_3d_adaptive_ragged_baseline(benchmark::State& state)
{
    run_cqmom_benchmark(state, make_adaptive_ragged_fixture(), qmomkit::CqmomDiagnosticLevel::none);
}

void bm_cqmom_4d_moderate_fixed(benchmark::State& state)
{
    run_recursive_cqmom_application_benchmark(state, make_independent_fixed_fixture({4, 2, 2, 2}));
}

void bm_cqmom_6d_low_node_deep_tree(benchmark::State& state)
{
    run_recursive_cqmom_application_benchmark(state, make_independent_fixed_fixture({2, 2, 2, 2, 2, 2}));
}

void bm_cqmom_6d_wide_conditional_reconstruction(benchmark::State& state)
{
    run_recursive_cqmom_application_benchmark(state, make_independent_fixed_fixture({8, 2, 2, 2, 2, 2}));
}

void bm_cqmom_diagnostic_overhead(benchmark::State& state, const qmomkit::CqmomDiagnosticLevel diagnostic_level)
{
    run_recursive_cqmom_benchmark(state, make_independent_fixed_fixture({4, 2, 2, 2}), diagnostic_level, qmomkit::VandermondeSolver::automatic);
}

void bm_cqmom_representation_scaling_public_dense(benchmark::State& state)
{
    run_cqmom_public_dense_representation_benchmark(state, static_cast<std::size_t>(state.range(0)), static_cast<std::size_t>(state.range(1)));
}

void bm_cqmom_representation_scaling_public_lean(benchmark::State& state)
{
    run_cqmom_public_lean_representation_benchmark(state, static_cast<std::size_t>(state.range(0)), static_cast<std::size_t>(state.range(1)));
}

void bm_cqmom_representation_scaling_recursive_dense(benchmark::State& state)
{
    run_cqmom_recursive_dense_representation_benchmark(state, static_cast<std::size_t>(state.range(0)), static_cast<std::size_t>(state.range(1)));
}

void bm_cqmom_representation_scaling_recursive_lean(benchmark::State& state)
{
    run_cqmom_recursive_lean_representation_benchmark(state, static_cast<std::size_t>(state.range(0)), static_cast<std::size_t>(state.range(1)));
}

void bm_cqmom_representation_scaling_kernel_dense(benchmark::State& state)
{
    run_cqmom_kernel_dense_representation_benchmark(state, static_cast<std::size_t>(state.range(0)), static_cast<std::size_t>(state.range(1)));
}

void bm_cqmom_representation_scaling_kernel_lean(benchmark::State& state)
{
    run_cqmom_kernel_lean_representation_benchmark(state, static_cast<std::size_t>(state.range(0)), static_cast<std::size_t>(state.range(1)));
}

void apply_representation_scaling_args(benchmark::Benchmark* benchmark)
{
    for (const std::size_t dimensions : {std::size_t{2}, std::size_t{3}, std::size_t{4}, std::size_t{5}, std::size_t{6}})
    {
        benchmark->Args({static_cast<int>(dimensions), 2});
        benchmark->Args({static_cast<int>(dimensions), 4});
    }
    for (const std::size_t dimensions : {std::size_t{2}, std::size_t{3}, std::size_t{4}})
    {
        benchmark->Args({static_cast<int>(dimensions), 8});
    }
}

BENCHMARK(bm_cqmom_2d_conditioning_order_baseline);
BENCHMARK(bm_cqmom_3d_adaptive_ragged_baseline);
BENCHMARK(bm_cqmom_4d_moderate_fixed);
BENCHMARK(bm_cqmom_6d_low_node_deep_tree);
BENCHMARK(bm_cqmom_6d_wide_conditional_reconstruction);

BENCHMARK_CAPTURE(bm_cqmom_diagnostic_overhead, none, qmomkit::CqmomDiagnosticLevel::none);
BENCHMARK_CAPTURE(bm_cqmom_diagnostic_overhead, basic, qmomkit::CqmomDiagnosticLevel::basic);
BENCHMARK_CAPTURE(bm_cqmom_diagnostic_overhead, detailed, qmomkit::CqmomDiagnosticLevel::detailed);

BENCHMARK(bm_cqmom_representation_scaling_public_dense)->Apply(apply_representation_scaling_args);
BENCHMARK(bm_cqmom_representation_scaling_public_lean)->Apply(apply_representation_scaling_args);
BENCHMARK(bm_cqmom_representation_scaling_recursive_dense)->Apply(apply_representation_scaling_args);
BENCHMARK(bm_cqmom_representation_scaling_recursive_lean)->Apply(apply_representation_scaling_args);
BENCHMARK(bm_cqmom_representation_scaling_kernel_dense)->Apply(apply_representation_scaling_args);
BENCHMARK(bm_cqmom_representation_scaling_kernel_lean)->Apply(apply_representation_scaling_args);

} // namespace
