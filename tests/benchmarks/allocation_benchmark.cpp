#include "qmomkit/qmomkit.hpp"

#include <benchmark/benchmark.h>
#include <cstddef>
#include <cstdlib>
#include <new>
#include <optional>
#include <utility>
#include <vector>

// Heap-allocation counting benchmark.
//
// Profiling the CQMOM phase-space evolution hot path attributes roughly a third
// of wall-clock samples to `operator new`. Wall-clock benchmarks cannot
// separate that cost from arithmetic, so this executable counts allocations
// directly: it replaces global `operator new`/`operator delete` and reports
// allocations and bytes per measured operation as Google Benchmark counters.
//
// The replacement is process-wide, so counting is gated on a flag that is
// enabled only around the measured call. Allocations performed by Google
// Benchmark itself are therefore excluded.

namespace
{

std::size_t g_allocation_count = 0;
std::size_t g_allocated_bytes = 0;
bool g_counting_enabled = false;

class AllocationScope
{
public:
    AllocationScope() noexcept
    {
        g_counting_enabled = true;
    }

    AllocationScope(const AllocationScope&) = delete;
    AllocationScope& operator=(const AllocationScope&) = delete;

    ~AllocationScope()
    {
        g_counting_enabled = false;
    }
};

void reset_allocation_counters() noexcept
{
    g_allocation_count = 0;
    g_allocated_bytes = 0;
}

void report_allocation_counters(benchmark::State& state)
{
    const double iterations = static_cast<double>(state.iterations());
    if (iterations <= 0.0)
    {
        return;
    }
    state.counters["allocations_per_op"] = static_cast<double>(g_allocation_count) / iterations;
    state.counters["bytes_per_op"] = static_cast<double>(g_allocated_bytes) / iterations;
}

} // namespace

void* operator new(std::size_t size)
{
    if (g_counting_enabled)
    {
        ++g_allocation_count;
        g_allocated_bytes += size;
    }
    void* const memory = std::malloc(size == 0 ? 1 : size);
    if (memory == nullptr)
    {
        throw std::bad_alloc();
    }
    return memory;
}

void operator delete(void* memory) noexcept
{
    std::free(memory);
}

void operator delete(void* memory, std::size_t) noexcept
{
    std::free(memory);
}

namespace
{

// Raw moments of a separable product quadrature with two nodes per dimension.
qmomkit::RawMomentTensor product_moment_tensor(const std::size_t dimensions)
{
    const std::vector<double> weights{0.5, 0.5};
    const std::vector<std::vector<double>> abscissas{{1.0, 2.0}, {0.5, 1.5}, {0.75, 1.25}};

    std::vector<std::size_t> extents(dimensions, 4);
    std::size_t total = 1;
    for (const std::size_t extent : extents)
    {
        total *= extent;
    }

    qmomkit::RawMomentTensor tensor{extents, std::vector<double>(total, 0.0)};
    for (std::size_t flat = 0; flat < total; ++flat)
    {
        const qmomkit::MultiIndex order = qmomkit::unflatten_multi_index(flat, extents);
        double moment = 0.0;
        for (std::size_t a = 0; a < 2; ++a)
        {
            for (std::size_t b = 0; b < (dimensions > 1 ? 2u : 1u); ++b)
            {
                for (std::size_t c = 0; c < (dimensions > 2 ? 2u : 1u); ++c)
                {
                    const std::size_t node[3] = {a, b, c};
                    double term = 1.0;
                    for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                    {
                        term *= weights[node[dimension]];
                        for (std::size_t power = 0; power < order[dimension]; ++power)
                        {
                            term *= abscissas[dimension][node[dimension]];
                        }
                    }
                    moment += term;
                }
            }
        }
        tensor.values[flat] = moment;
    }
    return tensor;
}

qmomkit::CqmomPreparedCase make_cqmom_prepared_case(const std::vector<std::size_t>& node_counts)
{
    const qmomkit::PhaseSpaceHierarchy hierarchy(node_counts);
    const qmomkit::RawMomentTensor dense = product_moment_tensor(node_counts.size());
    qmomkit::CqmomCase lean{qmomkit::extract_cqmom_moments(dense, qmomkit::create_cqmom_moment_layout(hierarchy))};
    return qmomkit::CqmomPreparedCase{std::move(lean)};
}

qmomkit::PhaseSpaceVelocity make_velocity()
{
    return qmomkit::make_phase_space_velocity(
        [](double, const std::vector<double>& coordinate)
        {
            std::vector<double> velocity(coordinate.size(), 0.0);
            for (std::size_t dimension = 0; dimension < coordinate.size(); ++dimension)
            {
                velocity[dimension] = 0.1 * coordinate[dimension];
            }
            return velocity;
        });
}

qmomkit::PhaseSpaceVelocity make_directly_written_velocity()
{
    return [](double, const std::vector<double>& coordinates, std::vector<double>& velocities, std::size_t)
    {
        for (std::size_t entry = 0; entry < coordinates.size(); ++entry)
        {
            velocities[entry] = 0.1 * coordinates[entry];
        }
    };
}

enum class VelocityForm
{
    adapted_per_node,
    packed
};

[[nodiscard]] qmomkit::PhaseSpaceVelocity velocity_for(const VelocityForm form)
{
    return form == VelocityForm::packed ? make_directly_written_velocity() : make_velocity();
}

void bm_allocations_cqmom_phase_space_step(benchmark::State& state, const std::vector<std::size_t>& node_counts,
                                           const VelocityForm form = VelocityForm::adapted_per_node)
{
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = make_cqmom_prepared_case(node_counts);
    const qmomkit::PhaseSpaceVelocity velocity = velocity_for(form);
    const qmomkit::CqmomMomentSource source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;

    qmomkit::CqmomState initial = cqmom_prepared_case.initial_state();
    const qmomkit::CqmomPhaseSpaceEvolutionResult warmup = qmomkit::advance_cqmom_phase_space(initial, 0.0, 1.0e-3, velocity, source, options);
    if (!warmup.ok())
    {
        state.SkipWithError("deterministic CQMOM phase-space fixture did not advance");
        return;
    }

    qmomkit::CqmomState current = initial;
    double time = 0.0;
    reset_allocation_counters();
    for (auto _ : state)
    {
        std::optional<qmomkit::CqmomPhaseSpaceEvolutionResult> result;
        {
            const AllocationScope scope;
            result.emplace(qmomkit::advance_cqmom_phase_space(current, time, 1.0e-3, velocity, source, options));
        }
        if (!result->ok())
        {
            state.SkipWithError("CQMOM phase-space step failed");
            return;
        }
        current = std::move(result->state);
        time += 1.0e-3;
    }
    report_allocation_counters(state);
}

void bm_allocations_qmom_phase_space_step(benchmark::State& state, const VelocityForm form = VelocityForm::adapted_per_node)
{
    const qmomkit::MomentSet moments({1.0, 2.0, 14.0 / 3.0, 12.0, 98.0 / 3.0, 92.0});
    const qmomkit::QmomCase qmom_case{moments, 3, {}};
    const qmomkit::PhaseSpaceVelocity velocity = velocity_for(form);
    const qmomkit::QmomMomentSource source = [](double, const qmomkit::QmomSourceContext& context)
    { return qmomkit::QmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    qmomkit::QmomPhaseSpaceEvolutionOptions options;
    options.scheme = qmomkit::PhaseSpaceEvolutionScheme::strang;

    const qmomkit::QmomState initial = qmom_case.initial_state();
    const qmomkit::QmomPhaseSpaceEvolutionResult warmup = qmomkit::advance_qmom_phase_space(initial, 0.0, 1.0e-3, velocity, source, options);
    if (!warmup.ok())
    {
        state.SkipWithError("deterministic QMOM phase-space fixture did not advance");
        return;
    }

    qmomkit::QmomState current = initial;
    double time = 0.0;
    reset_allocation_counters();
    for (auto _ : state)
    {
        std::optional<qmomkit::QmomPhaseSpaceEvolutionResult> result;
        {
            const AllocationScope scope;
            result.emplace(qmomkit::advance_qmom_phase_space(current, time, 1.0e-3, velocity, source, options));
        }
        if (!result->ok())
        {
            state.SkipWithError("QMOM phase-space step failed");
            return;
        }
        current = std::move(result->state);
        time += 1.0e-3;
    }
    report_allocation_counters(state);
}

void bm_allocations_cqmom_reconstruction(benchmark::State& state, const std::vector<std::size_t>& node_counts)
{
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = make_cqmom_prepared_case(node_counts);
    const qmomkit::CqmomResult setup_check = qmomkit::compute_cqmom(cqmom_prepared_case);
    if (!setup_check.ok())
    {
        state.SkipWithError("deterministic CQMOM reconstruction fixture is not realizable");
        return;
    }

    reset_allocation_counters();
    for (auto _ : state)
    {
        const AllocationScope scope;
        const qmomkit::CqmomResult result = qmomkit::compute_cqmom(cqmom_prepared_case);
        benchmark::DoNotOptimize(result.quadrature.weights().data());
    }
    report_allocation_counters(state);
}

void bm_allocations_qmom_reconstruction(benchmark::State& state)
{
    const qmomkit::MomentSet moments({1.0, 2.0, 14.0 / 3.0, 12.0, 98.0 / 3.0, 92.0});
    const qmomkit::QmomResult setup_check = qmomkit::compute_qmom(moments, 3);
    if (!setup_check.ok())
    {
        state.SkipWithError("deterministic QMOM reconstruction fixture is not realizable");
        return;
    }

    reset_allocation_counters();
    for (auto _ : state)
    {
        const AllocationScope scope;
        const qmomkit::QmomResult result = qmomkit::compute_qmom(moments, 3);
        benchmark::DoNotOptimize(result.weights.data());
    }
    report_allocation_counters(state);
}

BENCHMARK_CAPTURE(bm_allocations_cqmom_phase_space_step, cqmom_2d_2x2, std::vector<std::size_t>{2, 2})
    ->Name("allocations_cqmom_phase_space_step_2d_2x2");
BENCHMARK_CAPTURE(bm_allocations_cqmom_phase_space_step, cqmom_3d_2x2x2, std::vector<std::size_t>{2, 2, 2})
    ->Name("allocations_cqmom_phase_space_step_3d_2x2x2");
BENCHMARK_CAPTURE(bm_allocations_cqmom_phase_space_step, cqmom_2d_2x2_packed, std::vector<std::size_t>{2, 2}, VelocityForm::packed)
    ->Name("allocations_cqmom_phase_space_step_2d_2x2_packed_velocity");
BENCHMARK_CAPTURE(bm_allocations_cqmom_phase_space_step, cqmom_3d_2x2x2_packed, std::vector<std::size_t>{2, 2, 2}, VelocityForm::packed)
    ->Name("allocations_cqmom_phase_space_step_3d_2x2x2_packed_velocity");
BENCHMARK_CAPTURE(bm_allocations_qmom_phase_space_step, qmom_n3, VelocityForm::adapted_per_node)->Name("allocations_qmom_phase_space_step_n3");
BENCHMARK_CAPTURE(bm_allocations_qmom_phase_space_step, qmom_n3_packed, VelocityForm::packed)
    ->Name("allocations_qmom_phase_space_step_n3_packed_velocity");
BENCHMARK_CAPTURE(bm_allocations_cqmom_reconstruction, cqmom_2d_2x2, std::vector<std::size_t>{2, 2})->Name("allocations_cqmom_reconstruction_2d_2x2");
BENCHMARK(bm_allocations_qmom_reconstruction)->Name("allocations_qmom_reconstruction_n3");

} // namespace
