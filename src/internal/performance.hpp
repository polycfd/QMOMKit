#pragma once

#include <cstddef>

namespace qmomkit::internal::performance
{

enum class VandermondeSolveKind
{
    dense,
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
};

struct LinearAlgebraCounters
{
    std::size_t vandermonde_solves = 0;
    std::size_t bjorck_pereyra_solves = 0;
    std::size_t dense_linear_solves = 0;
    std::size_t dense_vandermonde_fallbacks = 0;
    std::size_t residual_checks = 0;
    std::size_t condition_estimates = 0;
};

struct TransportCounters
{
    std::size_t transport_macro_steps = 0;
    std::size_t continuous_operator_calls = 0;
    std::size_t discontinuous_operator_calls = 0;
    std::size_t moment_reconstruction_calls = 0;
    std::size_t post_source_inversions = 0;
};

struct TimeIntegrationCounters
{
    std::size_t source_rk_stages = 0;
    std::size_t source_rk_retries = 0;
    std::size_t node_rhs_evaluations = 0;
    std::size_t rk23_accepted_steps = 0;
    std::size_t rk23_rejected_steps = 0;
};

struct PerformanceCounters
{
    ReconstructionCounters reconstruction;
    LinearAlgebraCounters linear_algebra;
    TransportCounters transport;
    TimeIntegrationCounters time_integration;
};

#if defined(QMOMKIT_ENABLE_PERFORMANCE_INSTRUMENTATION)

class CounterSession
{
public:
    [[nodiscard]] const PerformanceCounters& counters() const noexcept
    {
        return counters_;
    }

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

inline void record_vandermonde_solve(const VandermondeSolveKind algorithm) noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->linear_algebra.vandermonde_solves;
        if (algorithm == VandermondeSolveKind::bjorck_pereyra)
        {
            ++counters->linear_algebra.bjorck_pereyra_solves;
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

inline void record_transport_macro_step() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->transport.transport_macro_steps;
    }
}

inline void record_continuous_operator_call() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->transport.continuous_operator_calls;
    }
}

inline void record_discontinuous_operator_call() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->transport.discontinuous_operator_calls;
    }
}

inline void record_moment_reconstruction_call() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->transport.moment_reconstruction_calls;
    }
}

inline void record_post_source_inversion() noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        ++counters->transport.post_source_inversions;
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

inline void record_rk23_step(const bool accepted) noexcept
{
    if (PerformanceCounters* counters = active_counters())
    {
        if (accepted)
        {
            ++counters->time_integration.rk23_accepted_steps;
        }
        else
        {
            ++counters->time_integration.rk23_rejected_steps;
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
inline void record_vandermonde_solve(VandermondeSolveKind) noexcept {}
inline void record_dense_linear_solve() noexcept {}
inline void record_vandermonde_fallback() noexcept {}
inline void record_residual_check() noexcept {}
inline void record_condition_estimate() noexcept {}
inline void record_transport_macro_step() noexcept {}
inline void record_continuous_operator_call() noexcept {}
inline void record_discontinuous_operator_call() noexcept {}
inline void record_moment_reconstruction_call() noexcept {}
inline void record_post_source_inversion() noexcept {}
inline void record_source_rk_stage() noexcept {}
inline void record_source_retry() noexcept {}
inline void record_node_rhs_evaluation() noexcept {}
inline void record_rk23_step(bool) noexcept {}

#endif

} // namespace qmomkit::internal::performance
