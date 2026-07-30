#include "internal/cqmom_active_counts.hpp"
#include "internal/cqmom_case.hpp"
#include "qmomkit/qmomkit.hpp"

#include <benchmark/benchmark.h>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

enum class CqmomStorage
{
    dense,
    lean
};

struct QmomEventFixture
{
    std::size_t maximum_node_count = 0;
    std::size_t active_node_count_before = 0;
    std::size_t active_node_count_after = 0;
    qmomkit::QmomOptions qmom_options;
    qmomkit::QmomState initial_state;
    std::vector<double> target_moments;
};

struct CqmomExpectedBranchCount
{
    std::size_t conditioning_level = 0;
    qmomkit::ConditioningPath parent_path;
    std::size_t active_node_count = 0;
};

struct CqmomEventFixture
{
    qmomkit::PhaseSpaceHierarchy hierarchy;
    std::vector<std::size_t> extents;
    std::vector<double> initial_moments;
    std::vector<double> target_moments;
    qmomkit::QmomOptions qmom_options;
    std::vector<CqmomExpectedBranchCount> expected_after_counts;
    std::size_t expected_terminal_nodes_after = 0;
};

struct QmomFallbackFixture
{
    std::vector<double> moments;
    std::size_t direct_ceiling = 0;
    std::size_t fallback_ceiling = 0;
    std::size_t accepted_node_count = 0;
};

struct CqmomFallbackFixture
{
    qmomkit::PhaseSpaceHierarchy hierarchy;
    std::vector<std::size_t> extents;
    std::vector<double> moments;
    std::vector<qmomkit::CqmomActiveBranchCount> direct_active_counts;
    std::vector<qmomkit::CqmomActiveBranchCount> fallback_active_counts;
    std::size_t direct_ceiling = 0;
    std::size_t fallback_ceiling = 0;
    std::size_t accepted_count = 0;
    std::size_t expected_terminal_nodes = 0;
};

std::vector<double> raw_moments_from_1d_quadrature(const std::vector<double>& abscissas, const std::vector<double>& weights,
                                                   const std::size_t moment_count)
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

std::vector<double> raw_moments_from_terminal_quadrature(const std::vector<std::vector<double>>& abscissas, const std::vector<double>& weights,
                                                         const std::vector<std::size_t>& extents)
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

std::vector<double> canonical_moments_from_qmom_result(const qmomkit::QmomResult& result, const std::size_t moment_count)
{
    return raw_moments_from_1d_quadrature(result.abscissas, result.weights, moment_count);
}

qmomkit::QmomOptions adaptive_qmom_options()
{
    qmomkit::QmomOptions options;
    options.adaptive = true;
    options.rmin = 1.0e-10;
    options.eabs = 1.0e-12;
    options.cutoff = 1.0e-12;
    return options;
}

qmomkit::QmomState make_qmom_state(const std::vector<double>& moments, const std::size_t maximum_node_count, const qmomkit::QmomOptions& options)
{
    return qmomkit::QmomCase{qmomkit::MomentSet(moments), maximum_node_count, options}.initial_state();
}

QmomEventFixture make_qmom_event_fixture(std::vector<double> initial_moments, std::vector<double> target_moments,
                                         const std::size_t maximum_node_count, const std::size_t active_after)
{
    const qmomkit::QmomOptions options = adaptive_qmom_options();
    qmomkit::QmomState state = make_qmom_state(initial_moments, maximum_node_count, options);
    return {maximum_node_count, state.active_node_count(), active_after, options, std::move(state), std::move(target_moments)};
}

QmomEventFixture qmom_retain_max_fixture(const std::size_t maximum_node_count)
{
    std::vector<double> weights;
    std::vector<double> abscissas;
    weights.reserve(maximum_node_count);
    abscissas.reserve(maximum_node_count);
    double weight_sum = 0.0;
    for (std::size_t node = 0; node < maximum_node_count; ++node)
    {
        weights.push_back(1.0 + static_cast<double>(node));
        weight_sum += weights.back();
        const double fraction = maximum_node_count == 1 ? 0.0 : static_cast<double>(node) / static_cast<double>(maximum_node_count - 1);
        abscissas.push_back(-0.9 + 1.8 * fraction);
    }
    for (double& weight : weights)
    {
        weight /= weight_sum;
    }
    const std::vector<double> moments = raw_moments_from_1d_quadrature(abscissas, weights, 2 * maximum_node_count);
    return make_qmom_event_fixture(moments, moments, maximum_node_count, maximum_node_count);
}

QmomEventFixture qmom_rejected_growth_fixture(const std::size_t maximum_node_count)
{
    const std::size_t active = maximum_node_count - 1;
    std::vector<double> active_abscissas;
    std::vector<double> active_weights;
    active_abscissas.reserve(active);
    active_weights.reserve(active);
    double weight_sum = 0.0;
    for (std::size_t node = 0; node < active; ++node)
    {
        const double fraction = active == 1 ? 0.0 : static_cast<double>(node) / static_cast<double>(active - 1);
        active_abscissas.push_back(-0.75 + 2.0 * fraction);
        active_weights.push_back(1.0 + static_cast<double>(node));
        weight_sum += active_weights.back();
    }
    for (double& weight : active_weights)
    {
        weight /= weight_sum;
    }
    const std::vector<double> moments = raw_moments_from_1d_quadrature(active_abscissas, active_weights, 2 * maximum_node_count);
    return make_qmom_event_fixture(moments, moments, maximum_node_count, active);
}

QmomEventFixture qmom_retain_n2_max2_rejected_growth_pair_fixture()
{
    const QmomEventFixture rejected_growth = qmom_rejected_growth_fixture(3);
    return make_qmom_event_fixture(rejected_growth.target_moments, rejected_growth.target_moments, 2, 2);
}

QmomEventFixture qmom_reduce_fixture(const std::size_t maximum_node_count, const std::size_t active_after)
{
    const QmomEventFixture retained = qmom_retain_max_fixture(maximum_node_count);
    std::vector<double> abscissas;
    std::vector<double> weights;
    abscissas.reserve(active_after);
    weights.reserve(active_after);
    double weight_sum = 0.0;
    for (std::size_t node = 0; node < active_after; ++node)
    {
        const double fraction = active_after == 1 ? 0.0 : static_cast<double>(node) / static_cast<double>(active_after - 1);
        abscissas.push_back(0.2 + 1.4 * fraction);
        weights.push_back(2.0 + static_cast<double>(node));
        weight_sum += weights.back();
    }
    for (double& weight : weights)
    {
        weight /= weight_sum;
    }
    return make_qmom_event_fixture(retained.initial_state.moments(), raw_moments_from_1d_quadrature(abscissas, weights, 2 * maximum_node_count),
                                   maximum_node_count, active_after);
}

QmomEventFixture qmom_growth_fixture(const std::size_t maximum_node_count, const std::size_t active_before)
{
    std::vector<double> initial_abscissas;
    std::vector<double> initial_weights;
    initial_abscissas.reserve(active_before);
    initial_weights.reserve(active_before);
    double initial_weight_sum = 0.0;
    for (std::size_t node = 0; node < active_before; ++node)
    {
        const double fraction = active_before == 1 ? 0.0 : static_cast<double>(node) / static_cast<double>(active_before - 1);
        initial_abscissas.push_back(-0.6 + 1.3 * fraction);
        initial_weights.push_back(1.0 + static_cast<double>(node));
        initial_weight_sum += initial_weights.back();
    }
    for (double& weight : initial_weights)
    {
        weight /= initial_weight_sum;
    }

    const std::size_t active_after = active_before + 1;
    std::vector<double> target_abscissas;
    std::vector<double> target_weights;
    target_abscissas.reserve(active_after);
    target_weights.reserve(active_after);
    double target_weight_sum = 0.0;
    for (std::size_t node = 0; node < active_after; ++node)
    {
        const double fraction = active_after == 1 ? 0.0 : static_cast<double>(node) / static_cast<double>(active_after - 1);
        target_abscissas.push_back(-0.8 + 1.7 * fraction);
        target_weights.push_back(1.0 + static_cast<double>((node + 1) * 2));
        target_weight_sum += target_weights.back();
    }
    for (double& weight : target_weights)
    {
        weight /= target_weight_sum;
    }
    return make_qmom_event_fixture(raw_moments_from_1d_quadrature(initial_abscissas, initial_weights, 2 * maximum_node_count),
                                   raw_moments_from_1d_quadrature(target_abscissas, target_weights, 2 * maximum_node_count), maximum_node_count,
                                   active_after);
}

QmomFallbackFixture qmom_direct_fallback_to_one_fixture(const std::size_t fallback_ceiling)
{
    return {raw_moments_from_1d_quadrature({0.35}, {1.0}, 12), 1, fallback_ceiling, 1};
}

QmomFallbackFixture qmom_direct_fallback_to_two_fixture(const std::size_t fallback_ceiling)
{
    return {raw_moments_from_1d_quadrature({-0.55, 1.15}, {0.42, 0.58}, 12), 2, fallback_ceiling, 2};
}

qmomkit::QmomResult select_qmom_event_quadrature(const QmomEventFixture& fixture)
{
    qmomkit::QmomResult current =
        qmomkit::compute_qmom(qmomkit::MomentSet(fixture.target_moments), fixture.active_node_count_before, fixture.qmom_options);
    if (!current.ok())
    {
        return current;
    }
    if (current.weights.size() < fixture.active_node_count_before || fixture.active_node_count_before >= fixture.maximum_node_count)
    {
        return current;
    }

    qmomkit::QmomResult growth =
        qmomkit::compute_qmom(qmomkit::MomentSet(fixture.target_moments), fixture.active_node_count_before + 1, fixture.qmom_options);
    if (growth.ok() && growth.weights.size() == fixture.active_node_count_before + 1)
    {
        return growth;
    }
    return current;
}

void validate_qmom_fixture(const QmomEventFixture& fixture)
{
    if (fixture.initial_state.active_node_count() != fixture.active_node_count_before)
    {
        throw std::runtime_error("QMOM fixture did not construct the expected initial active node count.");
    }
    const qmomkit::QmomResult selected = select_qmom_event_quadrature(fixture);
    if (!selected.ok() || selected.weights.size() != fixture.active_node_count_after)
    {
        throw std::runtime_error("QMOM fixture did not produce the expected post-source node-count event.");
    }
    (void)canonical_moments_from_qmom_result(selected, fixture.target_moments.size());
}

qmomkit::MomentSource qmom_source_to_target(const QmomEventFixture& fixture)
{
    const std::vector<double> initial = fixture.initial_state.moments();
    const std::vector<double> target = fixture.target_moments;
    return [initial, target](double, const std::vector<double>& state)
    {
        if (state.size() != initial.size())
        {
            throw std::runtime_error("QMOM source fixture received an unexpected moment layout.");
        }
        std::vector<double> rate(state.size(), 0.0);
        for (std::size_t i = 0; i < rate.size(); ++i)
        {
            rate[i] = target[i] - initial[i];
        }
        return rate;
    };
}

void validate_qmom_step_fixture(const QmomEventFixture& fixture)
{
    validate_qmom_fixture(fixture);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    qmomkit::QmomTransportOptions options;
    options.qmom_options = fixture.qmom_options;
    options.scheme = qmomkit::TransportScheme::dc;
    const qmomkit::QmomTransportResult result = qmomkit::advance_transport_qmom(fixture.initial_state, fixture.maximum_node_count, 0.0, 1.0,
                                                                                zero_velocity, qmom_source_to_target(fixture), options);
    if (!result.ok() || result.state.active_node_count() != fixture.active_node_count_after)
    {
        throw std::runtime_error("QMOM full-step fixture did not produce the expected accepted event.");
    }
    (void)canonical_moments_from_qmom_result(result.state.quadrature(), result.state.moments().size());
}

void set_qmom_counters(benchmark::State& state, const QmomEventFixture& fixture)
{
    state.counters["max_nodes"] = static_cast<double>(fixture.maximum_node_count);
    state.counters["active_before"] = static_cast<double>(fixture.active_node_count_before);
    state.counters["active_after"] = static_cast<double>(fixture.active_node_count_after);
    state.counters["moment_count"] = static_cast<double>(fixture.target_moments.size());
}

qmomkit::QmomResult compute_qmom_fallback_fixture(const QmomFallbackFixture& fixture, const std::size_t ceiling)
{
    return qmomkit::compute_qmom(qmomkit::MomentSet(fixture.moments), ceiling, adaptive_qmom_options());
}

void validate_qmom_fallback_fixture(const QmomFallbackFixture& fixture)
{
    const qmomkit::QmomResult direct = compute_qmom_fallback_fixture(fixture, fixture.direct_ceiling);
    const qmomkit::QmomResult fallback = compute_qmom_fallback_fixture(fixture, fixture.fallback_ceiling);
    if (!direct.ok() || !fallback.ok())
    {
        throw std::runtime_error("QMOM direct/fallback fixture reconstruction failed.");
    }
    if (direct.weights.size() != fixture.accepted_node_count || fallback.weights.size() != fixture.accepted_node_count)
    {
        throw std::runtime_error("QMOM direct/fallback fixture did not return the expected accepted count.");
    }
    for (std::size_t node = 0; node < fixture.accepted_node_count; ++node)
    {
        if (std::abs(direct.weights[node] - fallback.weights[node]) > 1.0e-10 ||
            std::abs(direct.abscissas[node] - fallback.abscissas[node]) > 1.0e-10)
        {
            throw std::runtime_error("QMOM direct/fallback quadratures are not equivalent.");
        }
    }
    const std::vector<double> direct_moments = canonical_moments_from_qmom_result(direct, fixture.moments.size());
    const std::vector<double> fallback_moments = canonical_moments_from_qmom_result(fallback, fixture.moments.size());
    for (std::size_t i = 0; i < direct_moments.size(); ++i)
    {
        if (std::abs(direct_moments[i] - fallback_moments[i]) > 1.0e-10)
        {
            throw std::runtime_error("QMOM direct/fallback canonical moments are not equivalent.");
        }
    }
}

void set_qmom_fallback_counters(benchmark::State& state, const QmomFallbackFixture& fixture, const std::size_t ceiling)
{
    state.counters["ceiling"] = static_cast<double>(ceiling);
    state.counters["direct_ceiling"] = static_cast<double>(fixture.direct_ceiling);
    state.counters["fallback_ceiling"] = static_cast<double>(fixture.fallback_ceiling);
    state.counters["accepted_count"] = static_cast<double>(fixture.accepted_node_count);
    state.counters["moment_count"] = static_cast<double>(fixture.moments.size());
}

void run_qmom_direct_fallback_reconstruction(benchmark::State& state, const QmomFallbackFixture& fixture, const std::size_t ceiling)
{
    try
    {
        validate_qmom_fallback_fixture(fixture);
    }
    catch (const std::exception& error)
    {
        state.SkipWithError(error.what());
        return;
    }
    set_qmom_fallback_counters(state, fixture, ceiling);

    for (auto _ : state)
    {
        const qmomkit::QmomResult result = compute_qmom_fallback_fixture(fixture, ceiling);
        const std::vector<double> canonical = canonical_moments_from_qmom_result(result, fixture.moments.size());
        benchmark::DoNotOptimize(result.weights.data());
        benchmark::DoNotOptimize(result.abscissas.data());
        benchmark::DoNotOptimize(canonical.data());
    }
}

void run_qmom_reconstruction_event(benchmark::State& state, const QmomEventFixture& fixture)
{
    try
    {
        validate_qmom_fixture(fixture);
    }
    catch (const std::exception& error)
    {
        state.SkipWithError(error.what());
        return;
    }
    set_qmom_counters(state, fixture);

    for (auto _ : state)
    {
        qmomkit::QmomResult selected = select_qmom_event_quadrature(fixture);
        std::vector<double> canonical = canonical_moments_from_qmom_result(selected, fixture.target_moments.size());
        benchmark::DoNotOptimize(selected.weights.data());
        benchmark::DoNotOptimize(selected.abscissas.data());
        benchmark::DoNotOptimize(canonical.data());
    }
}

void run_qmom_full_step(benchmark::State& state, const QmomEventFixture& fixture)
{
    try
    {
        validate_qmom_step_fixture(fixture);
    }
    catch (const std::exception& error)
    {
        state.SkipWithError(error.what());
        return;
    }
    set_qmom_counters(state, fixture);

    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    const qmomkit::MomentSource source = qmom_source_to_target(fixture);
    qmomkit::QmomTransportOptions options;
    options.qmom_options = fixture.qmom_options;
    options.scheme = qmomkit::TransportScheme::dc;

    for (auto _ : state)
    {
        state.PauseTiming();
        const qmomkit::QmomState iteration_state = fixture.initial_state;
        state.ResumeTiming();
        const qmomkit::QmomTransportResult result =
            qmomkit::advance_transport_qmom(iteration_state, fixture.maximum_node_count, 0.0, 1.0, zero_velocity, source, options);
        benchmark::DoNotOptimize(result.state.moments().data());
        benchmark::DoNotOptimize(result.state.quadrature().weights.data());
    }
}

struct TwoDimensionalBranches
{
    std::vector<double> root_abscissas;
    std::vector<double> root_weights;
    std::vector<std::vector<double>> child_abscissas;
    std::vector<std::vector<double>> child_weights;
};

std::vector<double> moments_from_2d_branches(const TwoDimensionalBranches& branches, const std::vector<std::size_t>& extents)
{
    std::vector<double> weights;
    std::vector<std::vector<double>> abscissas(2);
    for (std::size_t root = 0; root < branches.root_abscissas.size(); ++root)
    {
        for (std::size_t child = 0; child < branches.child_abscissas[root].size(); ++child)
        {
            weights.push_back(branches.root_weights[root] * branches.child_weights[root][child]);
            abscissas[0].push_back(branches.root_abscissas[root]);
            abscissas[1].push_back(branches.child_abscissas[root][child]);
        }
    }
    return raw_moments_from_terminal_quadrature(abscissas, weights, extents);
}

TwoDimensionalBranches make_2d_branches(const std::vector<std::size_t>& child_counts)
{
    TwoDimensionalBranches branches;
    branches.root_abscissas = {-0.7, 1.1};
    branches.root_weights = {0.45, 0.55};
    branches.child_abscissas.reserve(child_counts.size());
    branches.child_weights.reserve(child_counts.size());
    for (std::size_t root = 0; root < child_counts.size(); ++root)
    {
        std::vector<double> nodes;
        std::vector<double> weights;
        double weight_sum = 0.0;
        for (std::size_t child = 0; child < child_counts[root]; ++child)
        {
            const double fraction = child_counts[root] == 1 ? 0.0 : static_cast<double>(child) / static_cast<double>(child_counts[root] - 1);
            nodes.push_back(-1.2 + 2.1 * fraction + 0.2 * static_cast<double>(root));
            weights.push_back(1.0 + static_cast<double>(child + root + 1));
            weight_sum += weights.back();
        }
        for (double& weight : weights)
        {
            weight /= weight_sum;
        }
        branches.child_abscissas.push_back(std::move(nodes));
        branches.child_weights.push_back(std::move(weights));
    }
    return branches;
}

struct ThreeDimensionalBranches
{
    TwoDimensionalBranches xy;
    std::vector<std::vector<std::vector<double>>> z_abscissas;
    std::vector<std::vector<std::vector<double>>> z_weights;
};

ThreeDimensionalBranches make_3d_branches(const std::vector<std::vector<std::size_t>>& z_counts, const std::vector<std::size_t>& y_counts = {2, 3})
{
    ThreeDimensionalBranches branches;
    branches.xy = make_2d_branches(y_counts);
    branches.z_abscissas.resize(y_counts.size());
    branches.z_weights.resize(y_counts.size());
    for (std::size_t root = 0; root < y_counts.size(); ++root)
    {
        branches.z_abscissas[root].resize(y_counts[root]);
        branches.z_weights[root].resize(y_counts[root]);
        for (std::size_t child = 0; child < y_counts[root]; ++child)
        {
            const std::size_t count = z_counts[root][child];
            double weight_sum = 0.0;
            for (std::size_t z = 0; z < count; ++z)
            {
                const double fraction = count == 1 ? 0.0 : static_cast<double>(z) / static_cast<double>(count - 1);
                branches.z_abscissas[root][child].push_back(-0.55 + fraction + 0.15 * static_cast<double>(root) + 0.05 * static_cast<double>(child));
                branches.z_weights[root][child].push_back(1.0 + static_cast<double>(z + child + 1));
                weight_sum += branches.z_weights[root][child].back();
            }
            for (double& weight : branches.z_weights[root][child])
            {
                weight /= weight_sum;
            }
        }
    }
    return branches;
}

std::vector<double> moments_from_3d_branches(const ThreeDimensionalBranches& branches, const std::vector<std::size_t>& extents)
{
    std::vector<double> weights;
    std::vector<std::vector<double>> abscissas(3);
    for (std::size_t root = 0; root < branches.xy.root_abscissas.size(); ++root)
    {
        for (std::size_t child = 0; child < branches.xy.child_abscissas[root].size(); ++child)
        {
            for (std::size_t z = 0; z < branches.z_abscissas[root][child].size(); ++z)
            {
                weights.push_back(branches.xy.root_weights[root] * branches.xy.child_weights[root][child] * branches.z_weights[root][child][z]);
                abscissas[0].push_back(branches.xy.root_abscissas[root]);
                abscissas[1].push_back(branches.xy.child_abscissas[root][child]);
                abscissas[2].push_back(branches.z_abscissas[root][child][z]);
            }
        }
    }
    return raw_moments_from_terminal_quadrature(abscissas, weights, extents);
}

CqmomEventFixture cqmom_2d_conditional_retain_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3});
    const std::vector<std::size_t> extents{4, 6};
    const std::vector<double> moments = moments_from_2d_branches(make_2d_branches({3, 3}), extents);
    return {hierarchy, extents, moments, moments, adaptive_qmom_options(), {{1, qmomkit::ConditioningPath({0}), 3}}, 6};
}

CqmomEventFixture cqmom_2d_conditional_rejected_growth_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3});
    const std::vector<std::size_t> extents{4, 6};
    const std::vector<double> moments = moments_from_2d_branches(make_2d_branches({2, 2}), extents);
    return {hierarchy, extents, moments, moments, adaptive_qmom_options(), {{1, qmomkit::ConditioningPath({0}), 2}}, 4};
}

CqmomEventFixture cqmom_2d_conditional_retain_n2_max2_rejected_growth_pair_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    const std::vector<double> moments = moments_from_2d_branches(make_2d_branches({2, 2}), extents);
    return {hierarchy, extents, moments, moments, adaptive_qmom_options(), {{1, qmomkit::ConditioningPath({0}), 2}}, 4};
}

CqmomEventFixture cqmom_2d_conditional_reduce_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3});
    const std::vector<std::size_t> extents{4, 6};
    return {hierarchy,
            extents,
            moments_from_2d_branches(make_2d_branches({3, 3}), extents),
            moments_from_2d_branches(make_2d_branches({2, 3}), extents),
            adaptive_qmom_options(),
            {{1, qmomkit::ConditioningPath({0}), 2}},
            5};
}

CqmomEventFixture cqmom_2d_conditional_growth_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3});
    const std::vector<std::size_t> extents{4, 6};
    return {hierarchy,
            extents,
            moments_from_2d_branches(make_2d_branches({2, 2}), extents),
            moments_from_2d_branches(make_2d_branches({3, 2}), extents),
            adaptive_qmom_options(),
            {{1, qmomkit::ConditioningPath({0}), 3}},
            5};
}

CqmomEventFixture cqmom_2d_root_reduce_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({3, 2});
    const std::vector<std::size_t> extents{6, 4};
    const TwoDimensionalBranches initial{
        {-1.0, 0.25, 1.25}, {0.25, 0.35, 0.40}, {{-0.8, 0.7}, {-0.6, 1.0}, {-0.4, 1.2}}, {{0.45, 0.55}, {0.5, 0.5}, {0.35, 0.65}}};
    const TwoDimensionalBranches target{{-1.0, 1.25}, {0.45, 0.55}, {{-0.8, 0.7}, {-0.4, 1.2}}, {{0.45, 0.55}, {0.35, 0.65}}};
    return {hierarchy,
            extents,
            moments_from_2d_branches(initial, extents),
            moments_from_2d_branches(target, extents),
            adaptive_qmom_options(),
            {{0, qmomkit::ConditioningPath(), 2}},
            4};
}

CqmomEventFixture cqmom_2d_root_growth_nonidentity_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({3, 2}, {1, 0});
    const std::vector<std::size_t> extents{6, 4};
    const std::vector<double> initial = moments_from_2d_branches(make_2d_branches({2, 2}), extents);
    const TwoDimensionalBranches target{
        {-0.9, 0.2, 1.2}, {0.25, 0.35, 0.40}, {{-1.0, 0.8}, {-0.7, 1.0}, {-0.5, 1.2}}, {{0.4, 0.6}, {0.45, 0.55}, {0.35, 0.65}}};
    return {hierarchy, extents, initial, moments_from_2d_branches(target, extents), adaptive_qmom_options(), {{1, qmomkit::ConditioningPath({0}), 3}},
            6};
}

CqmomEventFixture cqmom_3d_ragged_retain_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 2});
    const std::vector<std::size_t> extents{4, 6, 4};
    const std::vector<double> moments = moments_from_3d_branches(make_3d_branches({{1, 2}, {2, 1, 2}}), extents);
    return {hierarchy, extents, moments, moments, adaptive_qmom_options(), {{2, qmomkit::ConditioningPath({0, 0}), 1}}, 8};
}

CqmomEventFixture cqmom_3d_terminal_growth_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 2});
    const std::vector<std::size_t> extents{4, 6, 4};
    return {hierarchy,
            extents,
            moments_from_3d_branches(make_3d_branches({{1, 2}, {2, 1, 2}}), extents),
            moments_from_3d_branches(make_3d_branches({{2, 2}, {2, 1, 2}}), extents),
            adaptive_qmom_options(),
            {{2, qmomkit::ConditioningPath({0, 0}), 2}},
            9};
}

CqmomEventFixture cqmom_3d_nonterminal_reduce_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3, 2});
    const std::vector<std::size_t> extents{4, 6, 4};
    return {hierarchy,
            extents,
            moments_from_3d_branches(make_3d_branches({{2, 2}, {2, 2, 2}}, {2, 3}), extents),
            moments_from_3d_branches(make_3d_branches({{2, 2}, {2, 2}}, {2, 2}), extents),
            adaptive_qmom_options(),
            {{1, qmomkit::ConditioningPath({1}), 2}},
            8};
}

qmomkit::CqmomActiveBranchCount make_active_branch_count(const qmomkit::PhaseSpaceHierarchy& hierarchy, const std::size_t conditioning_level,
                                                         qmomkit::ConditioningPath parent_path, const std::size_t active_node_count)
{
    const std::size_t phase_space_dimension = hierarchy.conditioning_order()[conditioning_level];
    return {conditioning_level, std::move(parent_path), phase_space_dimension, hierarchy.node_count(phase_space_dimension), active_node_count};
}

std::vector<qmomkit::CqmomActiveBranchCount> make_2d_conditional_active_counts(const qmomkit::PhaseSpaceHierarchy& hierarchy,
                                                                               const std::size_t child_active_count)
{
    return {make_active_branch_count(hierarchy, 0, qmomkit::ConditioningPath(), 2),
            make_active_branch_count(hierarchy, 1, qmomkit::ConditioningPath({0}), child_active_count),
            make_active_branch_count(hierarchy, 1, qmomkit::ConditioningPath({1}), child_active_count)};
}

CqmomFallbackFixture cqmom_2d_conditional_n3_to_n2_fallback_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 3});
    const std::vector<std::size_t> extents{4, 6};
    return {hierarchy,
            extents,
            moments_from_2d_branches(make_2d_branches({2, 2}), extents),
            make_2d_conditional_active_counts(hierarchy, 2),
            make_2d_conditional_active_counts(hierarchy, 3),
            2,
            3,
            2,
            4};
}

CqmomFallbackFixture cqmom_2d_conditional_n2_to_n1_fallback_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2});
    const std::vector<std::size_t> extents{4, 4};
    return {hierarchy,
            extents,
            moments_from_2d_branches(make_2d_branches({1, 1}), extents),
            make_2d_conditional_active_counts(hierarchy, 1),
            make_2d_conditional_active_counts(hierarchy, 2),
            1,
            2,
            1,
            2};
}

CqmomFallbackFixture cqmom_2d_parent_n3_to_n2_fallback_fixture()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy({3, 2});
    const std::vector<std::size_t> extents{6, 4};
    const TwoDimensionalBranches target{{-1.0, 1.25}, {0.45, 0.55}, {{-0.8, 0.7}, {-0.4, 1.2}}, {{0.45, 0.55}, {0.35, 0.65}}};
    return {hierarchy,
            extents,
            moments_from_2d_branches(target, extents),
            {make_active_branch_count(hierarchy, 0, qmomkit::ConditioningPath(), 2),
             make_active_branch_count(hierarchy, 1, qmomkit::ConditioningPath({0}), 2),
             make_active_branch_count(hierarchy, 1, qmomkit::ConditioningPath({1}), 2)},
            {make_active_branch_count(hierarchy, 0, qmomkit::ConditioningPath(), 3),
             make_active_branch_count(hierarchy, 1, qmomkit::ConditioningPath({0}), 2),
             make_active_branch_count(hierarchy, 1, qmomkit::ConditioningPath({1}), 2),
             make_active_branch_count(hierarchy, 1, qmomkit::ConditioningPath({2}), 2)},
            2,
            3,
            2,
            4};
}

qmomkit::CqmomDenseCase make_dense_case(const CqmomEventFixture& fixture, const std::vector<double>& moments)
{
    qmomkit::CqmomDenseCase dense_case{fixture.hierarchy, {fixture.extents, moments}, {}, fixture.qmom_options, true};
    dense_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::none;
    dense_case.vandermonde_solver = qmomkit::VandermondeSolver::automatic;
    return dense_case;
}

qmomkit::CqmomLeanCase make_lean_case(const CqmomEventFixture& fixture, const std::vector<double>& moments)
{
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(fixture.hierarchy);
    qmomkit::CqmomLeanCase lean_case{qmomkit::extract_cqmom_moments({fixture.extents, moments}, layout),
                                     {},
                                     fixture.qmom_options,
                                     true,
                                     qmomkit::CqmomDiagnosticLevel::none,
                                     qmomkit::VandermondeSolver::automatic};
    return lean_case;
}

qmomkit::CqmomCase make_cqmom_case(const CqmomEventFixture& fixture, const std::vector<double>& moments, const CqmomStorage storage)
{
    if (storage == CqmomStorage::dense)
    {
        return qmomkit::internal::create_cqmom_case(make_dense_case(fixture, moments));
    }
    return qmomkit::internal::create_cqmom_case(make_lean_case(fixture, moments));
}

qmomkit::CqmomDenseCase make_dense_case(const CqmomFallbackFixture& fixture)
{
    qmomkit::CqmomDenseCase dense_case{fixture.hierarchy, {fixture.extents, fixture.moments}, {}, adaptive_qmom_options(), true};
    dense_case.diagnostic_level = qmomkit::CqmomDiagnosticLevel::none;
    dense_case.vandermonde_solver = qmomkit::VandermondeSolver::automatic;
    return dense_case;
}

qmomkit::CqmomLeanCase make_lean_case(const CqmomFallbackFixture& fixture)
{
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(fixture.hierarchy);
    return {qmomkit::extract_cqmom_moments({fixture.extents, fixture.moments}, layout),
            {},
            adaptive_qmom_options(),
            true,
            qmomkit::CqmomDiagnosticLevel::none,
            qmomkit::VandermondeSolver::automatic};
}

qmomkit::CqmomResult compute_cqmom_fallback_fixture(const CqmomFallbackFixture& fixture,
                                                    const std::vector<qmomkit::CqmomActiveBranchCount>& active_counts, const CqmomStorage storage)
{
    if (storage == CqmomStorage::dense)
    {
        return qmomkit::internal::compute_trusted_cqmom(make_dense_case(fixture), active_counts, false);
    }
    return qmomkit::internal::compute_trusted_cqmom(make_lean_case(fixture), active_counts, false);
}

void compare_cqmom_active_counts(const std::vector<qmomkit::CqmomActiveBranchCount>& direct,
                                 const std::vector<qmomkit::CqmomActiveBranchCount>& fallback)
{
    if (direct.size() != fallback.size())
    {
        throw std::runtime_error("CQMOM direct/fallback active branch-count sizes differ.");
    }
    for (std::size_t i = 0; i < direct.size(); ++i)
    {
        if (direct[i].conditioning_level != fallback[i].conditioning_level || direct[i].parent_path != fallback[i].parent_path ||
            direct[i].phase_space_dimension != fallback[i].phase_space_dimension ||
            direct[i].configured_node_count != fallback[i].configured_node_count || direct[i].active_node_count != fallback[i].active_node_count)
        {
            throw std::runtime_error("CQMOM direct/fallback active branch counts differ.");
        }
    }
}

void compare_cqmom_quadratures(const qmomkit::CqmomQuadrature& direct, const qmomkit::CqmomQuadrature& fallback)
{
    if (direct.total_nodes() != fallback.total_nodes() || direct.terminal_paths() != fallback.terminal_paths() ||
        direct.abscissas().size() != fallback.abscissas().size())
    {
        throw std::runtime_error("CQMOM direct/fallback accepted hierarchies differ.");
    }
    for (std::size_t node = 0; node < direct.weights().size(); ++node)
    {
        if (std::abs(direct.weights()[node] - fallback.weights()[node]) > 1.0e-10)
        {
            throw std::runtime_error("CQMOM direct/fallback weights differ.");
        }
    }
    for (std::size_t dimension = 0; dimension < direct.abscissas().size(); ++dimension)
    {
        if (direct.abscissas()[dimension].size() != fallback.abscissas()[dimension].size())
        {
            throw std::runtime_error("CQMOM direct/fallback abscissa storage differs.");
        }
        for (std::size_t node = 0; node < direct.abscissas()[dimension].size(); ++node)
        {
            if (std::abs(direct.abscissas()[dimension][node] - fallback.abscissas()[dimension][node]) > 1.0e-10)
            {
                throw std::runtime_error("CQMOM direct/fallback abscissas differ.");
            }
        }
    }
}

void compare_cqmom_reconstructed_moments(const CqmomFallbackFixture& fixture, const qmomkit::CqmomQuadrature& direct,
                                         const qmomkit::CqmomQuadrature& fallback)
{
    const qmomkit::RawMomentTensor direct_moments = qmomkit::reconstruct_raw_moments(direct, fixture.extents);
    const qmomkit::RawMomentTensor fallback_moments = qmomkit::reconstruct_raw_moments(fallback, fixture.extents);
    if (direct_moments.values.size() != fallback_moments.values.size())
    {
        throw std::runtime_error("CQMOM direct/fallback reconstructed moment sizes differ.");
    }
    for (std::size_t i = 0; i < direct_moments.values.size(); ++i)
    {
        if (std::abs(direct_moments.values[i] - fallback_moments.values[i]) > 1.0e-10)
        {
            throw std::runtime_error("CQMOM direct/fallback reconstructed moments differ.");
        }
    }
}

void validate_cqmom_fallback_fixture(const CqmomFallbackFixture& fixture, const CqmomStorage storage)
{
    const qmomkit::CqmomResult direct = compute_cqmom_fallback_fixture(fixture, fixture.direct_active_counts, storage);
    const qmomkit::CqmomResult fallback = compute_cqmom_fallback_fixture(fixture, fixture.fallback_active_counts, storage);
    if (!direct.ok() || !fallback.ok())
    {
        throw std::runtime_error("CQMOM direct/fallback fixture reconstruction failed.");
    }
    if (direct.quadrature.total_nodes() != fixture.expected_terminal_nodes || fallback.quadrature.total_nodes() != fixture.expected_terminal_nodes)
    {
        throw std::runtime_error("CQMOM direct/fallback fixture returned an unexpected terminal-node count.");
    }
    compare_cqmom_quadratures(direct.quadrature, fallback.quadrature);
    compare_cqmom_active_counts(qmomkit::internal::active_branch_counts_from_quadrature(direct.quadrature),
                                qmomkit::internal::active_branch_counts_from_quadrature(fallback.quadrature));
    compare_cqmom_reconstructed_moments(fixture, direct.quadrature, fallback.quadrature);
}

void set_cqmom_fallback_counters(benchmark::State& state, const CqmomFallbackFixture& fixture, const CqmomStorage storage, const std::size_t ceiling)
{
    state.counters["dimensions"] = static_cast<double>(fixture.hierarchy.number_of_dimensions());
    state.counters["configured_terminal_nodes"] = static_cast<double>(fixture.hierarchy.total_nodes());
    state.counters["terminal_nodes"] = static_cast<double>(fixture.expected_terminal_nodes);
    state.counters["ceiling"] = static_cast<double>(ceiling);
    state.counters["direct_ceiling"] = static_cast<double>(fixture.direct_ceiling);
    state.counters["fallback_ceiling"] = static_cast<double>(fixture.fallback_ceiling);
    state.counters["accepted_count"] = static_cast<double>(fixture.accepted_count);
    state.counters["moment_count"] = static_cast<double>(fixture.moments.size());
    state.counters["storage"] = storage == CqmomStorage::dense ? 0.0 : 1.0;
}

void run_cqmom_direct_fallback_reconstruction(benchmark::State& state, const CqmomFallbackFixture& fixture, const CqmomStorage storage,
                                              const bool fallback)
{
    try
    {
        validate_cqmom_fallback_fixture(fixture, storage);
    }
    catch (const std::exception& error)
    {
        state.SkipWithError(error.what());
        return;
    }

    const std::vector<qmomkit::CqmomActiveBranchCount>& active_counts = fallback ? fixture.fallback_active_counts : fixture.direct_active_counts;
    set_cqmom_fallback_counters(state, fixture, storage, fallback ? fixture.fallback_ceiling : fixture.direct_ceiling);
    for (auto _ : state)
    {
        const qmomkit::CqmomResult result = compute_cqmom_fallback_fixture(fixture, active_counts, storage);
        if (storage == CqmomStorage::dense)
        {
            const qmomkit::RawMomentTensor canonical = qmomkit::reconstruct_raw_moments(result.quadrature, fixture.extents);
            benchmark::DoNotOptimize(canonical.values.data());
        }
        else
        {
            const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(fixture.hierarchy);
            const qmomkit::CqmomMomentData canonical = qmomkit::reconstruct_cqmom_moments(result.quadrature, layout);
            benchmark::DoNotOptimize(canonical.values.data());
        }
        benchmark::DoNotOptimize(result.quadrature.weights().data());
        benchmark::DoNotOptimize(result.quadrature.abscissas().data());
    }
}

std::optional<qmomkit::CqmomActiveBranchCount> find_active_count(const std::vector<qmomkit::CqmomActiveBranchCount>& counts,
                                                                 const std::size_t conditioning_level, const qmomkit::ConditioningPath& parent_path)
{
    for (const qmomkit::CqmomActiveBranchCount& count : counts)
    {
        if (count.conditioning_level == conditioning_level && count.parent_path == parent_path)
        {
            return count;
        }
    }
    return std::nullopt;
}

void validate_expected_cqmom_counts(const CqmomEventFixture& fixture, const qmomkit::CqmomQuadrature& quadrature)
{
    if (quadrature.total_nodes() != fixture.expected_terminal_nodes_after)
    {
        throw std::runtime_error("CQMOM fixture produced an unexpected terminal-node count.");
    }
    const std::vector<qmomkit::CqmomActiveBranchCount> after_counts = qmomkit::internal::active_branch_counts_from_quadrature(quadrature);
    for (const CqmomExpectedBranchCount& expected : fixture.expected_after_counts)
    {
        const std::optional<qmomkit::CqmomActiveBranchCount> count =
            find_active_count(after_counts, expected.conditioning_level, expected.parent_path);
        if (!count.has_value() || count->active_node_count != expected.active_node_count)
        {
            throw std::runtime_error("CQMOM fixture did not produce the expected branch-local active count.");
        }
    }
}

qmomkit::CqmomResult run_cqmom_event_reconstruction(const CqmomEventFixture& fixture, const CqmomStorage storage,
                                                    const std::vector<qmomkit::CqmomActiveBranchCount>& active_counts)
{
    if (storage == CqmomStorage::dense)
    {
        return qmomkit::internal::compute_trusted_cqmom(make_dense_case(fixture, fixture.target_moments), active_counts, true);
    }
    return qmomkit::internal::compute_trusted_cqmom(make_lean_case(fixture, fixture.target_moments), active_counts, true);
}

void canonicalize_cqmom_result(const CqmomEventFixture& fixture, const qmomkit::CqmomResult& result, const CqmomStorage storage)
{
    if (storage == CqmomStorage::dense)
    {
        const qmomkit::RawMomentTensor canonical = qmomkit::reconstruct_raw_moments(result.quadrature, fixture.extents);
        benchmark::DoNotOptimize(canonical.values.data());
        return;
    }
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(fixture.hierarchy);
    const qmomkit::CqmomMomentData canonical = qmomkit::reconstruct_cqmom_moments(result.quadrature, layout);
    benchmark::DoNotOptimize(canonical.values.data());
}

void validate_cqmom_fixture(const CqmomEventFixture& fixture, const CqmomStorage storage)
{
    const qmomkit::CqmomCase initial_case = make_cqmom_case(fixture, fixture.initial_moments, storage);
    const qmomkit::CqmomState initial_state = initial_case.initial_state();
    const qmomkit::CqmomResult selected = run_cqmom_event_reconstruction(fixture, storage, initial_state.active_branch_counts());
    if (!selected.ok())
    {
        throw std::runtime_error("CQMOM fixture reconstruction failed: " + selected.message);
    }
    validate_expected_cqmom_counts(fixture, selected.quadrature);
    canonicalize_cqmom_result(fixture, selected, storage);
}

qmomkit::CqmomMomentSource cqmom_source_to_target(const qmomkit::CqmomState& initial_state, const qmomkit::CqmomMomentData& target)
{
    return [initial_state, target](double, const qmomkit::CqmomSourceContext& context)
    {
        std::vector<double> rate(context.size(), 0.0);
        for (std::size_t entry = 0; entry < context.size(); ++entry)
        {
            const std::optional<double> initial_value = initial_state.moment(context.index(entry));
            if (!initial_value.has_value())
            {
                throw std::runtime_error("CQMOM source fixture received an unexpected moment index.");
            }
            rate[entry] = target.moment(context.index(entry)) - *initial_value;
        }
        return qmomkit::CqmomMomentRate{rate};
    };
}

void validate_cqmom_state_canonical(const qmomkit::CqmomState& state, const qmomkit::CqmomQuadrature& quadrature)
{
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(state.hierarchy());
    const qmomkit::CqmomMomentData canonical = qmomkit::reconstruct_cqmom_moments(quadrature, layout);
    for (const qmomkit::MultiIndex& index : layout.required_indices)
    {
        const std::optional<double> state_value = state.moment(index);
        if (!state_value.has_value())
        {
            throw std::runtime_error("CQMOM state is missing a canonical moment index.");
        }
        if (std::abs(*state_value - canonical.moment(index)) > 1.0e-8)
        {
            throw std::runtime_error("CQMOM full-step result is not canonical with the accepted quadrature.");
        }
    }
}

void validate_cqmom_step_fixture(const CqmomEventFixture& fixture, const CqmomStorage storage)
{
    validate_cqmom_fixture(fixture, storage);
    const qmomkit::CqmomCase initial_case = make_cqmom_case(fixture, fixture.initial_moments, storage);
    const qmomkit::CqmomState initial_state = initial_case.initial_state();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(fixture.hierarchy);
    const qmomkit::CqmomMomentData target = qmomkit::extract_cqmom_moments({fixture.extents, fixture.target_moments}, layout);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    const qmomkit::CqmomTransportResult result =
        qmomkit::advance_transport_cqmom(initial_state, 0.0, 1.0, zero_velocity, cqmom_source_to_target(initial_state, target), options);
    if (!result.ok() || !result.quadrature.has_value())
    {
        throw std::runtime_error("CQMOM full-step fixture did not produce an accepted transport step.");
    }
    validate_expected_cqmom_counts(fixture, *result.quadrature);
    validate_cqmom_state_canonical(result.state, *result.quadrature);
}

void set_cqmom_counters(benchmark::State& state, const CqmomEventFixture& fixture, const CqmomStorage storage,
                        const std::size_t active_branch_count_before)
{
    state.counters["dimensions"] = static_cast<double>(fixture.hierarchy.number_of_dimensions());
    state.counters["configured_terminal_nodes"] = static_cast<double>(fixture.hierarchy.total_nodes());
    state.counters["terminal_nodes_after"] = static_cast<double>(fixture.expected_terminal_nodes_after);
    state.counters["active_branches_before"] = static_cast<double>(active_branch_count_before);
    state.counters["moment_count"] = static_cast<double>(fixture.target_moments.size());
    state.counters["storage"] = storage == CqmomStorage::dense ? 0.0 : 1.0;
}

void run_cqmom_reconstruction_event(benchmark::State& state, const CqmomEventFixture& fixture, const CqmomStorage storage)
{
    try
    {
        validate_cqmom_fixture(fixture, storage);
    }
    catch (const std::exception& error)
    {
        state.SkipWithError(error.what());
        return;
    }

    const qmomkit::CqmomCase initial_case = make_cqmom_case(fixture, fixture.initial_moments, storage);
    const qmomkit::CqmomState initial_state = initial_case.initial_state();
    const std::vector<qmomkit::CqmomActiveBranchCount> active_counts = initial_state.active_branch_counts();
    set_cqmom_counters(state, fixture, storage, active_counts.size());

    for (auto _ : state)
    {
        const qmomkit::CqmomResult result = run_cqmom_event_reconstruction(fixture, storage, active_counts);
        canonicalize_cqmom_result(fixture, result, storage);
        benchmark::DoNotOptimize(result.quadrature.weights().data());
        benchmark::DoNotOptimize(result.quadrature.abscissas().data());
    }
}

void run_cqmom_full_step(benchmark::State& state, const CqmomEventFixture& fixture, const CqmomStorage storage)
{
    try
    {
        validate_cqmom_step_fixture(fixture, storage);
    }
    catch (const std::exception& error)
    {
        state.SkipWithError(error.what());
        return;
    }

    const qmomkit::CqmomCase initial_case = make_cqmom_case(fixture, fixture.initial_moments, storage);
    const qmomkit::CqmomState initial_state = initial_case.initial_state();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(fixture.hierarchy);
    const qmomkit::CqmomMomentData target = qmomkit::extract_cqmom_moments({fixture.extents, fixture.target_moments}, layout);
    const qmomkit::CqmomMomentSource source = cqmom_source_to_target(initial_state, target);
    const auto zero_velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.0); };
    qmomkit::CqmomTransportOptions options;
    options.scheme = qmomkit::TransportScheme::dc;
    set_cqmom_counters(state, fixture, storage, initial_state.active_branch_counts().size());

    for (auto _ : state)
    {
        state.PauseTiming();
        const qmomkit::CqmomState iteration_state = initial_state;
        state.ResumeTiming();
        const qmomkit::CqmomTransportResult result = qmomkit::advance_transport_cqmom(iteration_state, 0.0, 1.0, zero_velocity, source, options);
        benchmark::DoNotOptimize(result.state.moment_values().data());
        if (result.quadrature.has_value())
        {
            benchmark::DoNotOptimize(result.quadrature->weights().data());
        }
    }
}

void register_qmom_benchmarks()
{
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Direct/N1", [](benchmark::State& state)
                                 { run_qmom_direct_fallback_reconstruction(state, qmom_direct_fallback_to_one_fixture(2), 1); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Fallback/N2_to_N1", [](benchmark::State& state)
                                 { run_qmom_direct_fallback_reconstruction(state, qmom_direct_fallback_to_one_fixture(2), 2); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Fallback/N3_to_N1", [](benchmark::State& state)
                                 { run_qmom_direct_fallback_reconstruction(state, qmom_direct_fallback_to_one_fixture(3), 3); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Fallback/N4_to_N1", [](benchmark::State& state)
                                 { run_qmom_direct_fallback_reconstruction(state, qmom_direct_fallback_to_one_fixture(4), 4); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Fallback/N6_to_N1", [](benchmark::State& state)
                                 { run_qmom_direct_fallback_reconstruction(state, qmom_direct_fallback_to_one_fixture(6), 6); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Direct/N2", [](benchmark::State& state)
                                 { run_qmom_direct_fallback_reconstruction(state, qmom_direct_fallback_to_two_fixture(3), 2); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Fallback/N3_to_N2", [](benchmark::State& state)
                                 { run_qmom_direct_fallback_reconstruction(state, qmom_direct_fallback_to_two_fixture(3), 3); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Fallback/N4_to_N2", [](benchmark::State& state)
                                 { run_qmom_direct_fallback_reconstruction(state, qmom_direct_fallback_to_two_fixture(4), 4); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Fallback/N6_to_N2", [](benchmark::State& state)
                                 { run_qmom_direct_fallback_reconstruction(state, qmom_direct_fallback_to_two_fixture(6), 6); });

    benchmark::RegisterBenchmark("Qmom/Reconstruction/RetainMax/N2",
                                 [](benchmark::State& state) { run_qmom_reconstruction_event(state, qmom_retain_max_fixture(2)); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/RetainMax/N3",
                                 [](benchmark::State& state) { run_qmom_reconstruction_event(state, qmom_retain_max_fixture(3)); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/RetainMax/N4",
                                 [](benchmark::State& state) { run_qmom_reconstruction_event(state, qmom_retain_max_fixture(4)); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/RetainMax/N6",
                                 [](benchmark::State& state) { run_qmom_reconstruction_event(state, qmom_retain_max_fixture(6)); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Retain/N2_Max2", [](benchmark::State& state)
                                 { run_qmom_reconstruction_event(state, qmom_retain_n2_max2_rejected_growth_pair_fixture()); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/RejectedGrowth/N2_Max3",
                                 [](benchmark::State& state) { run_qmom_reconstruction_event(state, qmom_rejected_growth_fixture(3)); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/RejectedGrowth/N3_Max4",
                                 [](benchmark::State& state) { run_qmom_reconstruction_event(state, qmom_rejected_growth_fixture(4)); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Reduce/N4_to_N3",
                                 [](benchmark::State& state) { run_qmom_reconstruction_event(state, qmom_reduce_fixture(4, 3)); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Reduce/N4_to_N2",
                                 [](benchmark::State& state) { run_qmom_reconstruction_event(state, qmom_reduce_fixture(4, 2)); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Reduce/N6_to_N2",
                                 [](benchmark::State& state) { run_qmom_reconstruction_event(state, qmom_reduce_fixture(6, 2)); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Reduce/N3_to_N1",
                                 [](benchmark::State& state) { run_qmom_reconstruction_event(state, qmom_reduce_fixture(3, 1)); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Grow/N1_to_N2_Max3",
                                 [](benchmark::State& state) { run_qmom_reconstruction_event(state, qmom_growth_fixture(3, 1)); });
    benchmark::RegisterBenchmark("Qmom/Reconstruction/Grow/N2_to_N3_Max4",
                                 [](benchmark::State& state) { run_qmom_reconstruction_event(state, qmom_growth_fixture(4, 2)); });

    benchmark::RegisterBenchmark("Qmom/FullStep/RetainMax/N3",
                                 [](benchmark::State& state) { run_qmom_full_step(state, qmom_retain_max_fixture(3)); });
    benchmark::RegisterBenchmark("Qmom/FullStep/Retain/N2_Max2",
                                 [](benchmark::State& state) { run_qmom_full_step(state, qmom_retain_n2_max2_rejected_growth_pair_fixture()); });
    benchmark::RegisterBenchmark("Qmom/FullStep/RejectedGrowth/N2_Max3",
                                 [](benchmark::State& state) { run_qmom_full_step(state, qmom_rejected_growth_fixture(3)); });
    benchmark::RegisterBenchmark("Qmom/FullStep/Reduce/N3_to_N1",
                                 [](benchmark::State& state) { run_qmom_full_step(state, qmom_reduce_fixture(3, 1)); });
    benchmark::RegisterBenchmark("Qmom/FullStep/Grow/N1_to_N2_Max3",
                                 [](benchmark::State& state) { run_qmom_full_step(state, qmom_growth_fixture(3, 1)); });
}

void register_cqmom_benchmarks()
{
    benchmark::RegisterBenchmark(
        "Cqmom2D/Reconstruction/Lean/Direct/ConditionalN2", [](benchmark::State& state)
        { run_cqmom_direct_fallback_reconstruction(state, cqmom_2d_conditional_n3_to_n2_fallback_fixture(), CqmomStorage::lean, false); });
    benchmark::RegisterBenchmark(
        "Cqmom2D/Reconstruction/Lean/Fallback/ConditionalN3_to_N2", [](benchmark::State& state)
        { run_cqmom_direct_fallback_reconstruction(state, cqmom_2d_conditional_n3_to_n2_fallback_fixture(), CqmomStorage::lean, true); });
    benchmark::RegisterBenchmark(
        "Cqmom2D/Reconstruction/Dense/Direct/ConditionalN2", [](benchmark::State& state)
        { run_cqmom_direct_fallback_reconstruction(state, cqmom_2d_conditional_n3_to_n2_fallback_fixture(), CqmomStorage::dense, false); });
    benchmark::RegisterBenchmark(
        "Cqmom2D/Reconstruction/Dense/Fallback/ConditionalN3_to_N2", [](benchmark::State& state)
        { run_cqmom_direct_fallback_reconstruction(state, cqmom_2d_conditional_n3_to_n2_fallback_fixture(), CqmomStorage::dense, true); });
    benchmark::RegisterBenchmark(
        "Cqmom2D/Reconstruction/Lean/Direct/ConditionalN1", [](benchmark::State& state)
        { run_cqmom_direct_fallback_reconstruction(state, cqmom_2d_conditional_n2_to_n1_fallback_fixture(), CqmomStorage::lean, false); });
    benchmark::RegisterBenchmark(
        "Cqmom2D/Reconstruction/Lean/Fallback/ConditionalN2_to_N1", [](benchmark::State& state)
        { run_cqmom_direct_fallback_reconstruction(state, cqmom_2d_conditional_n2_to_n1_fallback_fixture(), CqmomStorage::lean, true); });
    benchmark::RegisterBenchmark(
        "Cqmom2D/Reconstruction/Dense/Direct/ConditionalN1", [](benchmark::State& state)
        { run_cqmom_direct_fallback_reconstruction(state, cqmom_2d_conditional_n2_to_n1_fallback_fixture(), CqmomStorage::dense, false); });
    benchmark::RegisterBenchmark(
        "Cqmom2D/Reconstruction/Dense/Fallback/ConditionalN2_to_N1", [](benchmark::State& state)
        { run_cqmom_direct_fallback_reconstruction(state, cqmom_2d_conditional_n2_to_n1_fallback_fixture(), CqmomStorage::dense, true); });
    benchmark::RegisterBenchmark(
        "Cqmom2D/Reconstruction/Lean/Direct/ParentN2", [](benchmark::State& state)
        { run_cqmom_direct_fallback_reconstruction(state, cqmom_2d_parent_n3_to_n2_fallback_fixture(), CqmomStorage::lean, false); });
    benchmark::RegisterBenchmark(
        "Cqmom2D/Reconstruction/Lean/Fallback/ParentN3_to_N2", [](benchmark::State& state)
        { run_cqmom_direct_fallback_reconstruction(state, cqmom_2d_parent_n3_to_n2_fallback_fixture(), CqmomStorage::lean, true); });
    benchmark::RegisterBenchmark(
        "Cqmom2D/Reconstruction/Dense/Direct/ParentN2", [](benchmark::State& state)
        { run_cqmom_direct_fallback_reconstruction(state, cqmom_2d_parent_n3_to_n2_fallback_fixture(), CqmomStorage::dense, false); });
    benchmark::RegisterBenchmark(
        "Cqmom2D/Reconstruction/Dense/Fallback/ParentN3_to_N2", [](benchmark::State& state)
        { run_cqmom_direct_fallback_reconstruction(state, cqmom_2d_parent_n3_to_n2_fallback_fixture(), CqmomStorage::dense, true); });

    benchmark::RegisterBenchmark("Cqmom2D/Reconstruction/Dense/ConditionalRetain", [](benchmark::State& state)
                                 { run_cqmom_reconstruction_event(state, cqmom_2d_conditional_retain_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom2D/Reconstruction/Dense/ConditionalRejectedGrowth", [](benchmark::State& state)
                                 { run_cqmom_reconstruction_event(state, cqmom_2d_conditional_rejected_growth_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom2D/Reconstruction/Dense/ConditionalReduce", [](benchmark::State& state)
                                 { run_cqmom_reconstruction_event(state, cqmom_2d_conditional_reduce_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom2D/Reconstruction/Dense/ConditionalGrow", [](benchmark::State& state)
                                 { run_cqmom_reconstruction_event(state, cqmom_2d_conditional_growth_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom2D/Reconstruction/Dense/RootReduce", [](benchmark::State& state)
                                 { run_cqmom_reconstruction_event(state, cqmom_2d_root_reduce_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom2D/Reconstruction/Dense/RootGrowNonIdentityOrder", [](benchmark::State& state)
                                 { run_cqmom_reconstruction_event(state, cqmom_2d_root_growth_nonidentity_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom3D/Reconstruction/Dense/RaggedRetain", [](benchmark::State& state)
                                 { run_cqmom_reconstruction_event(state, cqmom_3d_ragged_retain_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom3D/Reconstruction/Dense/TerminalGrow", [](benchmark::State& state)
                                 { run_cqmom_reconstruction_event(state, cqmom_3d_terminal_growth_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom3D/Reconstruction/Dense/NonterminalReduce", [](benchmark::State& state)
                                 { run_cqmom_reconstruction_event(state, cqmom_3d_nonterminal_reduce_fixture(), CqmomStorage::dense); });

    benchmark::RegisterBenchmark("Cqmom2D/FullStep/Dense/ConditionalRetain", [](benchmark::State& state)
                                 { run_cqmom_full_step(state, cqmom_2d_conditional_retain_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom2D/FullStep/Lean/ConditionalRetain", [](benchmark::State& state)
                                 { run_cqmom_full_step(state, cqmom_2d_conditional_retain_fixture(), CqmomStorage::lean); });
    benchmark::RegisterBenchmark(
        "Cqmom2D/FullStep/Dense/ConditionalRetainN2_Max2", [](benchmark::State& state)
        { run_cqmom_full_step(state, cqmom_2d_conditional_retain_n2_max2_rejected_growth_pair_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark(
        "Cqmom2D/FullStep/Lean/ConditionalRetainN2_Max2", [](benchmark::State& state)
        { run_cqmom_full_step(state, cqmom_2d_conditional_retain_n2_max2_rejected_growth_pair_fixture(), CqmomStorage::lean); });
    benchmark::RegisterBenchmark("Cqmom2D/FullStep/Dense/ConditionalRejectedGrowth", [](benchmark::State& state)
                                 { run_cqmom_full_step(state, cqmom_2d_conditional_rejected_growth_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom2D/FullStep/Lean/ConditionalRejectedGrowth", [](benchmark::State& state)
                                 { run_cqmom_full_step(state, cqmom_2d_conditional_rejected_growth_fixture(), CqmomStorage::lean); });
    benchmark::RegisterBenchmark("Cqmom2D/FullStep/Dense/ConditionalReduce", [](benchmark::State& state)
                                 { run_cqmom_full_step(state, cqmom_2d_conditional_reduce_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom2D/FullStep/Lean/ConditionalReduce", [](benchmark::State& state)
                                 { run_cqmom_full_step(state, cqmom_2d_conditional_reduce_fixture(), CqmomStorage::lean); });
    benchmark::RegisterBenchmark("Cqmom2D/FullStep/Dense/ConditionalGrow", [](benchmark::State& state)
                                 { run_cqmom_full_step(state, cqmom_2d_conditional_growth_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom2D/FullStep/Lean/ConditionalGrow", [](benchmark::State& state)
                                 { run_cqmom_full_step(state, cqmom_2d_conditional_growth_fixture(), CqmomStorage::lean); });
    benchmark::RegisterBenchmark("Cqmom2D/FullStep/Dense/RootReduce",
                                 [](benchmark::State& state) { run_cqmom_full_step(state, cqmom_2d_root_reduce_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom2D/FullStep/Lean/RootReduce",
                                 [](benchmark::State& state) { run_cqmom_full_step(state, cqmom_2d_root_reduce_fixture(), CqmomStorage::lean); });
    benchmark::RegisterBenchmark("Cqmom3D/FullStep/Dense/RaggedRetain",
                                 [](benchmark::State& state) { run_cqmom_full_step(state, cqmom_3d_ragged_retain_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom3D/FullStep/Lean/RaggedRetain",
                                 [](benchmark::State& state) { run_cqmom_full_step(state, cqmom_3d_ragged_retain_fixture(), CqmomStorage::lean); });
    benchmark::RegisterBenchmark("Cqmom3D/FullStep/Dense/TerminalGrow", [](benchmark::State& state)
                                 { run_cqmom_full_step(state, cqmom_3d_terminal_growth_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom3D/FullStep/Lean/TerminalGrow",
                                 [](benchmark::State& state) { run_cqmom_full_step(state, cqmom_3d_terminal_growth_fixture(), CqmomStorage::lean); });
    benchmark::RegisterBenchmark("Cqmom3D/FullStep/Dense/NonterminalReduce", [](benchmark::State& state)
                                 { run_cqmom_full_step(state, cqmom_3d_nonterminal_reduce_fixture(), CqmomStorage::dense); });
    benchmark::RegisterBenchmark("Cqmom3D/FullStep/Lean/NonterminalReduce", [](benchmark::State& state)
                                 { run_cqmom_full_step(state, cqmom_3d_nonterminal_reduce_fixture(), CqmomStorage::lean); });
}

const bool registered = []
{
    register_qmom_benchmarks();
    register_cqmom_benchmarks();
    return true;
}();

} // namespace
