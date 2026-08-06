#pragma once

#include <cstddef>

namespace qmomkit::internal::performance
{

enum class VandermondeSolveKind
{
    dense,
    pivoted_qr,
    svd,
    bjorck_pereyra
};

struct ReconstructionCounters
{
    std::size_t qmom_reconstructions = 0;
    std::size_t cqmom_reconstructions = 0;
    std::size_t cqmom_branch_solves = 0;
    std::size_t wheeler_calls = 0;
    std::size_t wheeler_adaptive_reductions = 0;
    std::size_t terminal_nodes_assembled = 0;
    std::size_t moment_layout_builds = 0;
    std::size_t cqmom_prepared_metadata_builds = 0;
    std::size_t cqmom_prepared_metadata_uses = 0;
    std::size_t cqmom_dynamic_metadata_builds = 0;
};

struct LinearAlgebraCounters
{
    std::size_t vandermonde_solves = 0;
    std::size_t bjorck_pereyra_solves = 0;
    std::size_t pivoted_qr_solves = 0;
    std::size_t svd_solves = 0;
    std::size_t dense_linear_solves = 0;
    std::size_t dense_vandermonde_fallbacks = 0;
    std::size_t residual_checks = 0;
    std::size_t condition_estimates = 0;
};

struct PhaseSpaceEvolutionCounters
{
    std::size_t phase_space_evolution_macro_steps = 0;
    std::size_t advection_operator_calls = 0;
    std::size_t moment_source_operator_calls = 0;
    std::size_t moment_reconstruction_calls = 0;
    std::size_t post_source_inversions = 0;
};

struct CqmomTopologyCounters
{
    std::size_t active_count_derivations = 0;
    std::size_t active_count_validations = 0;
    std::size_t active_count_transition_builds = 0;
};

struct TimeIntegrationCounters
{
    std::size_t source_rk_stages = 0;
    std::size_t source_rk_retries = 0;
    std::size_t node_rhs_evaluations = 0;
    std::size_t embedded_rk_accepted_steps = 0;
    std::size_t embedded_rk_rejected_steps = 0;
};

struct PerformanceCounters
{
    ReconstructionCounters reconstruction;
    LinearAlgebraCounters linear_algebra;
    PhaseSpaceEvolutionCounters phase_space_evolution;
    CqmomTopologyCounters cqmom_topology;
    TimeIntegrationCounters time_integration;
};

#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)

class CounterSession
{
public:
    [[nodiscard]] PerformanceCounters snapshot() const noexcept
    {
        return counters_;
    }

    void reset() noexcept
    {
        counters_ = {};
    }

private:
    friend class CounterScope;
    friend PerformanceCounters* active_counters() noexcept;

    PerformanceCounters counters_;
};

inline thread_local CounterSession* active_session = nullptr;

[[nodiscard]] inline PerformanceCounters* active_counters() noexcept
{
    return active_session == nullptr ? nullptr : &active_session->counters_;
}

class CounterScope
{
public:
    explicit CounterScope(CounterSession& session) noexcept : previous_(active_session)
    {
        active_session = &session;
    }

    CounterScope(const CounterScope&) = delete;
    CounterScope& operator=(const CounterScope&) = delete;

    ~CounterScope()
    {
        active_session = previous_;
    }

private:
    CounterSession* previous_;
};

[[nodiscard]] inline bool enabled() noexcept
{
    return true;
}

inline void record_qmom_reconstruction() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->reconstruction.qmom_reconstructions;
    }
}

inline void record_cqmom_reconstruction() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->reconstruction.cqmom_reconstructions;
    }
}

inline void record_cqmom_branch_solve() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->reconstruction.cqmom_branch_solves;
    }
}

inline void record_wheeler_call() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->reconstruction.wheeler_calls;
    }
}

inline void record_wheeler_adaptive_reduction() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->reconstruction.wheeler_adaptive_reductions;
    }
}

inline void record_terminal_node_assembled() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->reconstruction.terminal_nodes_assembled;
    }
}

inline void record_moment_layout_build() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->reconstruction.moment_layout_builds;
    }
}

inline void record_cqmom_prepared_metadata_build() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->reconstruction.cqmom_prepared_metadata_builds;
    }
}

inline void record_cqmom_prepared_metadata_use() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->reconstruction.cqmom_prepared_metadata_uses;
    }
}

inline void record_cqmom_dynamic_metadata_build() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->reconstruction.cqmom_dynamic_metadata_builds;
    }
}

inline void record_vandermonde_solve(const VandermondeSolveKind algorithm) noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->linear_algebra.vandermonde_solves;
        if (algorithm == VandermondeSolveKind::bjorck_pereyra)
        {
            ++counters->linear_algebra.bjorck_pereyra_solves;
        }
        if (algorithm == VandermondeSolveKind::pivoted_qr)
        {
            ++counters->linear_algebra.pivoted_qr_solves;
        }
        if (algorithm == VandermondeSolveKind::svd)
        {
            ++counters->linear_algebra.svd_solves;
        }
    }
}

inline void record_dense_linear_solve() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->linear_algebra.dense_linear_solves;
    }
}

inline void record_vandermonde_fallback() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->linear_algebra.dense_vandermonde_fallbacks;
    }
}

inline void record_residual_check() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->linear_algebra.residual_checks;
    }
}

inline void record_condition_estimate() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->linear_algebra.condition_estimates;
    }
}

inline void record_phase_space_evolution_macro_step() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->phase_space_evolution.phase_space_evolution_macro_steps;
    }
}

inline void record_advection_operator_call() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->phase_space_evolution.advection_operator_calls;
    }
}

inline void record_moment_source_operator_call() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->phase_space_evolution.moment_source_operator_calls;
    }
}

inline void record_moment_reconstruction_call() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->phase_space_evolution.moment_reconstruction_calls;
    }
}

inline void record_post_source_inversion() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->phase_space_evolution.post_source_inversions;
    }
}

inline void record_active_count_derivation() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->cqmom_topology.active_count_derivations;
    }
}

inline void record_active_count_validation() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->cqmom_topology.active_count_validations;
    }
}

inline void record_active_count_transition_build() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->cqmom_topology.active_count_transition_builds;
    }
}

inline void record_source_rk_stage() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->time_integration.source_rk_stages;
    }
}

inline void record_source_retry() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->time_integration.source_rk_retries;
    }
}

inline void record_node_rhs_evaluation() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->time_integration.node_rhs_evaluations;
    }
}

inline void record_embedded_rk_step(const bool accepted) noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        if (accepted)
        {
            ++counters->time_integration.embedded_rk_accepted_steps;
        }
        else
        {
            ++counters->time_integration.embedded_rk_rejected_steps;
        }
    }
}

#else

class CounterSession
{
public:
    [[nodiscard]] PerformanceCounters snapshot() const noexcept
    {
        return {};
    }

    void reset() noexcept {}
};

class CounterScope
{
public:
    explicit CounterScope(CounterSession&) noexcept {}
};

[[nodiscard]] inline bool enabled() noexcept
{
    return false;
}

inline void record_qmom_reconstruction() noexcept {}
inline void record_cqmom_reconstruction() noexcept {}
inline void record_cqmom_branch_solve() noexcept {}
inline void record_wheeler_call() noexcept {}
inline void record_wheeler_adaptive_reduction() noexcept {}
inline void record_terminal_node_assembled() noexcept {}
inline void record_moment_layout_build() noexcept {}
inline void record_cqmom_prepared_metadata_build() noexcept {}
inline void record_cqmom_prepared_metadata_use() noexcept {}
inline void record_cqmom_dynamic_metadata_build() noexcept {}
inline void record_vandermonde_solve(VandermondeSolveKind) noexcept {}
inline void record_dense_linear_solve() noexcept {}
inline void record_vandermonde_fallback() noexcept {}
inline void record_residual_check() noexcept {}
inline void record_condition_estimate() noexcept {}
inline void record_phase_space_evolution_macro_step() noexcept {}
inline void record_advection_operator_call() noexcept {}
inline void record_moment_source_operator_call() noexcept {}
inline void record_moment_reconstruction_call() noexcept {}
inline void record_post_source_inversion() noexcept {}
inline void record_active_count_derivation() noexcept {}
inline void record_active_count_validation() noexcept {}
inline void record_active_count_transition_build() noexcept {}
inline void record_source_rk_stage() noexcept {}
inline void record_source_retry() noexcept {}
inline void record_node_rhs_evaluation() noexcept {}
inline void record_embedded_rk_step(bool) noexcept {}

#endif

} // namespace qmomkit::internal::performance
