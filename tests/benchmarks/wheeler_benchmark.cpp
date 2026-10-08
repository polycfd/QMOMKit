#include "qmomkit/qmomkit.hpp"

#include <benchmark/benchmark.h>
#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

std::vector<double> raw_moments_from_1d_quadrature(const std::vector<double>& weights, const std::vector<double>& abscissas,
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

void bm_wheeler_three_node(benchmark::State& state)
{
    const std::vector<double> moments = raw_moments_from_1d_quadrature({0.15, 0.35, 0.50}, {-1.25, 0.25, 2.0}, 6);
    const qmomkit::MomentSet moment_set(moments);
    qmomkit::WheelerOptions options;
    options.adaptive = false;
    const qmomkit::Quadrature1D setup_check = qmomkit::compute_quadrature(moment_set, 3, options);
    if (!setup_check.ok())
    {
        state.SkipWithError("deterministic Wheeler fixture is not realizable");
        return;
    }

    for (auto _ : state)
    {
        const qmomkit::Quadrature1D result = qmomkit::compute_quadrature(moment_set, 3, options);
        benchmark::DoNotOptimize(result.weights.data());
        benchmark::DoNotOptimize(result.abscissas.data());
    }
}

void bm_adaptive_wheeler_reduced_two_node(benchmark::State& state)
{
    const std::vector<double> moments = raw_moments_from_1d_quadrature({1.0}, {0.75}, 4);
    const qmomkit::MomentSet moment_set(moments);
    const qmomkit::WheelerOptions options{true, 1.0e-10, 1.0e-10, 1.0e-12};
    const qmomkit::Quadrature1D setup_check = qmomkit::compute_quadrature(moment_set, 2, options);
    if (!setup_check.ok())
    {
        state.SkipWithError("deterministic adaptive Wheeler fixture is not realizable");
        return;
    }

    for (auto _ : state)
    {
        const qmomkit::Quadrature1D result = qmomkit::compute_quadrature(moment_set, 2, options);
        benchmark::DoNotOptimize(result.weights.data());
        benchmark::DoNotOptimize(result.abscissas.data());
    }
}

BENCHMARK(bm_wheeler_three_node);
BENCHMARK(bm_adaptive_wheeler_reduced_two_node);

} // namespace
