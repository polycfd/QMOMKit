#include "internal/cqmom_active_counts.hpp"
#include "internal/cqmom_case.hpp"
#include "internal/cqmom_validation.hpp"
#include "internal/performance.hpp"
#include "internal/recursive_quadrature.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/wheeler.hpp"

#include <exception>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace qmomkit
{
namespace
{
[[nodiscard]] WheelerOptions to_wheeler_options(const QmomOptions& options)
{
    return WheelerOptions{options.adaptive, options.rmin, options.eabs, options.cutoff};
}

[[nodiscard]] QmomStatus to_qmom_status(const QuadratureStatus status)
{
    switch (status)
    {
    case QuadratureStatus::success:
        return QmomStatus::success;
    case QuadratureStatus::invalid_input:
        return QmomStatus::invalid_input;
    case QuadratureStatus::not_realizable:
        return QmomStatus::not_realizable;
    case QuadratureStatus::no_adaptive_quadrature:
        return QmomStatus::no_adaptive_quadrature;
    }
    return QmomStatus::invalid_input;
}
} // namespace

QmomResult compute_qmom(const MomentSet& moments, const std::size_t node_count, const QmomOptions& options)
{
    internal::performance::record_qmom_reconstruction();
    const Quadrature1D result = compute_quadrature(moments, node_count, to_wheeler_options(options));
    return QmomResult{result.abscissas, result.weights, to_qmom_status(result.status), result.message};
}

QmomResult compute_qmom(const QmomCase& case_)
{
    return compute_qmom(case_.moments, case_.node_count, case_.options);
}

struct CqmomCase::Data
{
    CqmomLeanCase lean_case;
};

struct CqmomState::Data
{
    CqmomLeanCase lean_case;
    std::vector<CqmomActiveBranchCount> active_branch_counts;
};

CqmomCase::CqmomCase(std::unique_ptr<Data> data) : data_(std::move(data)) {}

CqmomCase::CqmomCase(const CqmomCase& other) : data_(other.data_ == nullptr ? nullptr : std::make_unique<Data>(*other.data_)) {}

CqmomCase& CqmomCase::operator=(const CqmomCase& other)
{
    if (this != &other)
    {
        data_ = other.data_ == nullptr ? nullptr : std::make_unique<Data>(*other.data_);
    }
    return *this;
}

CqmomCase::CqmomCase(CqmomCase&&) noexcept = default;

CqmomCase& CqmomCase::operator=(CqmomCase&&) noexcept = default;

CqmomCase::~CqmomCase() = default;

const PhaseSpaceHierarchy& CqmomCase::hierarchy() const noexcept
{
    return data_->lean_case.moments.layout.hierarchy;
}

const QmomOptions& CqmomCase::qmom_options() const noexcept
{
    return data_->lean_case.qmom_options;
}

const ConditionalMomentOptions& CqmomCase::conditional_options() const noexcept
{
    return data_->lean_case.conditional_options;
}

bool CqmomCase::adaptive_branch_counts() const noexcept
{
    return data_->lean_case.adaptive_branch_counts;
}

CqmomDiagnosticLevel CqmomCase::diagnostic_level() const noexcept
{
    return data_->lean_case.diagnostic_level;
}

VandermondeSolver CqmomCase::vandermonde_solver() const noexcept
{
    return data_->lean_case.vandermonde_solver;
}

CqmomState CqmomCase::initial_state() const
{
    if (!data_->lean_case.adaptive_branch_counts)
    {
        return internal::create_cqmom_state(data_->lean_case);
    }

    const CqmomResult seeded = internal::compute_trusted_cqmom(data_->lean_case);
    if (!seeded.ok())
    {
        throw std::runtime_error("CQMOM initial state active-count seeding failed: " + seeded.message);
    }
    return internal::create_cqmom_state(data_->lean_case, internal::active_branch_counts_from_quadrature(seeded.quadrature));
}

CqmomState::CqmomState(std::unique_ptr<Data> data) : data_(std::move(data)) {}

CqmomState::CqmomState(const CqmomState& other) : data_(other.data_ == nullptr ? nullptr : std::make_unique<Data>(*other.data_)) {}

CqmomState& CqmomState::operator=(const CqmomState& other)
{
    if (this != &other)
    {
        data_ = other.data_ == nullptr ? nullptr : std::make_unique<Data>(*other.data_);
    }
    return *this;
}

CqmomState::CqmomState(CqmomState&&) noexcept = default;

CqmomState& CqmomState::operator=(CqmomState&&) noexcept = default;

CqmomState::~CqmomState() = default;

const PhaseSpaceHierarchy& CqmomState::hierarchy() const noexcept
{
    return data_->lean_case.moments.layout.hierarchy;
}

const QmomOptions& CqmomState::qmom_options() const noexcept
{
    return data_->lean_case.qmom_options;
}

const ConditionalMomentOptions& CqmomState::conditional_options() const noexcept
{
    return data_->lean_case.conditional_options;
}

bool CqmomState::adaptive_branch_counts() const noexcept
{
    return data_->lean_case.adaptive_branch_counts;
}

CqmomDiagnosticLevel CqmomState::diagnostic_level() const noexcept
{
    return data_->lean_case.diagnostic_level;
}

VandermondeSolver CqmomState::vandermonde_solver() const noexcept
{
    return data_->lean_case.vandermonde_solver;
}

const std::vector<CqmomActiveBranchCount>& CqmomState::active_branch_counts() const noexcept
{
    return data_->active_branch_counts;
}

namespace internal
{
CqmomCase CqmomConstructionAccess::create_case(CqmomLeanCase lean_case)
{
    return CqmomCase(std::make_unique<CqmomCase::Data>(CqmomCase::Data{std::move(lean_case)}));
}

CqmomState CqmomConstructionAccess::create_state(CqmomLeanCase lean_case)
{
    return CqmomState(std::make_unique<CqmomState::Data>(CqmomState::Data{std::move(lean_case), {}}));
}

CqmomState CqmomConstructionAccess::create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts)
{
    return CqmomState(std::make_unique<CqmomState::Data>(CqmomState::Data{std::move(lean_case), std::move(active_branch_counts)}));
}

const CqmomLeanCase& CqmomConstructionAccess::state_case(const CqmomState& state)
{
    return state.data_->lean_case;
}

const std::vector<CqmomActiveBranchCount>& CqmomConstructionAccess::state_active_branch_counts(const CqmomState& state)
{
    return state.data_->active_branch_counts;
}

CqmomQuadrature CqmomConstructionAccess::create_quadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights,
                                                           std::vector<std::vector<double>> abscissas, const bool adaptive_branch_counts_enabled,
                                                           const bool rectangular_storage, std::vector<ConditioningPath> terminal_paths)
{
    return CqmomQuadrature(std::move(hierarchy), std::move(weights), std::move(abscissas), adaptive_branch_counts_enabled, rectangular_storage,
                           std::move(terminal_paths));
}

CqmomCase create_cqmom_case(CqmomLeanCase lean_case)
{
    return CqmomConstructionAccess::create_case(std::move(lean_case));
}

CqmomState create_cqmom_state(CqmomLeanCase lean_case)
{
    return CqmomConstructionAccess::create_state(std::move(lean_case));
}

CqmomState create_cqmom_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts)
{
    return CqmomConstructionAccess::create_state(std::move(lean_case), std::move(active_branch_counts));
}

const CqmomLeanCase& cqmom_state_case(const CqmomState& state)
{
    return CqmomConstructionAccess::state_case(state);
}

const std::vector<CqmomActiveBranchCount>& cqmom_state_active_branch_counts(const CqmomState& state)
{
    return CqmomConstructionAccess::state_active_branch_counts(state);
}

CqmomQuadrature create_cqmom_quadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights, std::vector<std::vector<double>> abscissas,
                                        const bool adaptive_branch_counts_enabled, const bool rectangular_storage,
                                        std::vector<ConditioningPath> terminal_paths)
{
    return CqmomConstructionAccess::create_quadrature(std::move(hierarchy), std::move(weights), std::move(abscissas), adaptive_branch_counts_enabled,
                                                      rectangular_storage, std::move(terminal_paths));
}
} // namespace internal

CqmomQuadrature::CqmomQuadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights, std::vector<std::vector<double>> abscissas,
                                 const bool adaptive_branch_counts_enabled, const bool rectangular_storage,
                                 std::vector<ConditioningPath> terminal_paths)
    : hierarchy_(std::move(hierarchy)), weights_(std::move(weights)), abscissas_(std::move(abscissas)), terminal_paths_(std::move(terminal_paths)),
      adaptive_branch_counts_enabled_(adaptive_branch_counts_enabled), rectangular_storage_(rectangular_storage)
{
}

const PhaseSpaceHierarchy& CqmomQuadrature::hierarchy() const noexcept
{
    return hierarchy_;
}

std::size_t CqmomQuadrature::number_of_dimensions() const noexcept
{
    return hierarchy_.number_of_dimensions();
}

std::size_t CqmomQuadrature::total_nodes() const noexcept
{
    return weights_.size();
}

const std::vector<double>& CqmomQuadrature::weights() const noexcept
{
    return weights_;
}

const std::vector<std::vector<double>>& CqmomQuadrature::abscissas() const noexcept
{
    return abscissas_;
}

bool CqmomQuadrature::adaptive_branch_counts_enabled() const noexcept
{
    return adaptive_branch_counts_enabled_;
}

bool CqmomQuadrature::has_rectangular_storage() const noexcept
{
    return rectangular_storage_;
}

const std::vector<ConditioningPath>& CqmomQuadrature::terminal_paths() const noexcept
{
    return terminal_paths_;
}

double CqmomQuadrature::weight(const MultiIndex& node_index) const
{
    if (!rectangular_storage_)
    {
        throw std::logic_error("Ragged CQMOM quadrature cannot be addressed by rectangular node multi-index.");
    }
    return weights_.at(hierarchy_.flatten_node_index(node_index));
}

double CqmomQuadrature::abscissa(const std::size_t dimension, const MultiIndex& node_index) const
{
    if (!rectangular_storage_)
    {
        throw std::logic_error("Ragged CQMOM quadrature cannot be addressed by rectangular node multi-index.");
    }
    if (dimension >= hierarchy_.number_of_dimensions())
    {
        throw std::out_of_range("Phase-space dimension is outside the CQMOM quadrature.");
    }
    return abscissas_.at(dimension).at(hierarchy_.flatten_node_index(node_index));
}

std::vector<double> CqmomQuadrature::node(const MultiIndex& node_index) const
{
    std::vector<double> values;
    values.reserve(hierarchy_.number_of_dimensions());
    for (std::size_t dimension = 0; dimension < hierarchy_.number_of_dimensions(); ++dimension)
    {
        values.push_back(abscissa(dimension, node_index));
    }
    return values;
}

namespace
{
[[nodiscard]] CqmomStatus to_cqmom_status(const internal::RecursiveQuadratureStatus status)
{
    switch (status)
    {
    case internal::RecursiveQuadratureStatus::success:
        return CqmomStatus::success;
    case internal::RecursiveQuadratureStatus::invalid_input:
        return CqmomStatus::invalid_input;
    case internal::RecursiveQuadratureStatus::wheeler_failed:
        return CqmomStatus::qmom_failed;
    case internal::RecursiveQuadratureStatus::conditional_moment_failed:
        return CqmomStatus::conditional_moment_failed;
    }
    return CqmomStatus::invalid_input;
}

[[nodiscard]] CqmomDiagnostics to_cqmom_diagnostics(const internal::RecursiveQuadratureBranchDiagnostics& diagnostics)
{
    return CqmomDiagnostics{diagnostics.conditioning_level,
                            diagnostics.path,
                            diagnostics.phase_space_dimension,
                            diagnostics.requested_node_count,
                            diagnostics.actual_node_count,
                            diagnostics.adaptive_branch_counts_enabled,
                            diagnostics.adaptive_node_count_reduced,
                            to_qmom_status(diagnostics.wheeler_status),
                            diagnostics.conditional_moment_status,
                            diagnostics.conditional_moment_diagnostics,
                            diagnostics.conditional_vandermonde_requested_solver,
                            diagnostics.conditional_vandermonde_selected_solver,
                            diagnostics.dense_vandermonde_fallback_occurred,
                            diagnostics.message};
}

[[nodiscard]] std::vector<CqmomDiagnostics>
to_cqmom_diagnostics(const std::vector<internal::RecursiveQuadratureBranchDiagnostics>& branch_diagnostics)
{
    std::vector<CqmomDiagnostics> diagnostics;
    diagnostics.reserve(branch_diagnostics.size());
    for (const internal::RecursiveQuadratureBranchDiagnostics& diagnostic : branch_diagnostics)
    {
        diagnostics.push_back(to_cqmom_diagnostics(diagnostic));
    }
    return diagnostics;
}

[[nodiscard]] ConditioningPath phase_space_index_to_conditioning_path(const PhaseSpaceHierarchy& hierarchy, const MultiIndex& node_index)
{
    std::vector<std::size_t> path(hierarchy.number_of_dimensions(), 0);
    for (std::size_t level = 0; level < hierarchy.number_of_dimensions(); ++level)
    {
        path[level] = node_index[hierarchy.conditioning_order()[level]];
    }
    return ConditioningPath(std::move(path));
}

[[nodiscard]] std::vector<ConditioningPath> terminal_paths_for_public_quadrature(const internal::RecursiveQuadrature& quadrature)
{
    if (!quadrature.terminal_paths().empty() || !quadrature.has_rectangular_storage())
    {
        return quadrature.terminal_paths();
    }

    std::vector<ConditioningPath> paths;
    paths.reserve(quadrature.weights().size());
    for (const MultiIndex& node_index : enumerate_quadrature_node_multi_indices(quadrature.hierarchy()))
    {
        paths.push_back(phase_space_index_to_conditioning_path(quadrature.hierarchy(), node_index));
    }
    return paths;
}

[[nodiscard]] std::size_t failed_phase_space_dimension(const internal::RecursiveQuadratureResult& result, const PhaseSpaceHierarchy& hierarchy)
{
    std::size_t failed_dimension = 0;
    if (result.failed_conditioning_level < hierarchy.conditioning_order().size())
    {
        failed_dimension = hierarchy.conditioning_order()[result.failed_conditioning_level];
    }
    return failed_dimension;
}

[[nodiscard]] CqmomResult to_cqmom_result(const internal::RecursiveQuadratureResult& result, const PhaseSpaceHierarchy& hierarchy)
{
    return CqmomResult{
        internal::create_cqmom_quadrature(result.quadrature.hierarchy(), result.quadrature.weights(), result.quadrature.abscissas(),
                                          result.quadrature.adaptive_branch_counts_enabled(), result.quadrature.has_rectangular_storage(),
                                          terminal_paths_for_public_quadrature(result.quadrature)),
        to_cqmom_status(result.status), result.message, to_cqmom_diagnostics(result.diagnostics), failed_phase_space_dimension(result, hierarchy)};
}
} // namespace

CqmomResult compute_cqmom(const CqmomDenseCase& dense_case)
{
    return internal::compute_trusted_cqmom(dense_case);
}

CqmomResult compute_cqmom(const CqmomLeanCase& lean_case)
{
    const PhaseSpaceHierarchy& hierarchy = lean_case.moments.layout.hierarchy;
    const internal::RecursiveQuadratureResult result =
        internal::compute_recursive_quadrature_from_moment_data({hierarchy,
                                                                 lean_case.moments,
                                                                 lean_case.conditional_options,
                                                                 to_wheeler_options(lean_case.qmom_options),
                                                                 lean_case.adaptive_branch_counts,
                                                                 lean_case.diagnostic_level,
                                                                 lean_case.vandermonde_solver,
                                                                 {}});

    return to_cqmom_result(result, hierarchy);
}

namespace internal
{

CqmomResult compute_trusted_cqmom(const CqmomDenseCase& dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts)
{
    RecursiveQuadratureResult result{RecursiveQuadrature(dense_case.hierarchy), RecursiveQuadratureStatus::success, "", {}, {}, 0};
    result.quadrature.set_adaptive_branch_counts_enabled(dense_case.adaptive_branch_counts);
    std::optional<CqmomMomentData> lean_moments;

    try
    {
        validate_dense_cqmom_reconstruction_tensor(dense_case.moments, dense_case.hierarchy, "Raw moment extents must match the hierarchy dimension.",
                                                   "Raw moment data size does not match its extents.", "Raw moments must be finite.",
                                                   "Raw moment extents must contain at least 2*N orders per dimension.");
        const CqmomMomentLayout layout = create_cqmom_moment_layout(dense_case.hierarchy);
        lean_moments = extract_cqmom_moments(dense_case.moments, layout);
    }
    catch (const std::exception& ex)
    {
        result.status = RecursiveQuadratureStatus::invalid_input;
        result.message = ex.what();
    }
    if (lean_moments.has_value())
    {
        result = compute_recursive_quadrature_from_moment_data(
            {dense_case.hierarchy, *lean_moments, dense_case.conditional_options, to_wheeler_options(dense_case.qmom_options),
             dense_case.adaptive_branch_counts, dense_case.diagnostic_level, dense_case.vandermonde_solver, std::move(active_branch_counts)});
    }

    return to_cqmom_result(result, dense_case.hierarchy);
}

CqmomResult compute_trusted_cqmom(const CqmomLeanCase& lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts)
{
    const PhaseSpaceHierarchy& hierarchy = lean_case.moments.layout.hierarchy;
    const RecursiveQuadratureResult result = compute_recursive_quadrature_from_trusted_moment_data(
        {hierarchy, lean_case.moments, lean_case.conditional_options, to_wheeler_options(lean_case.qmom_options), lean_case.adaptive_branch_counts,
         lean_case.diagnostic_level, lean_case.vandermonde_solver, std::move(active_branch_counts)});

    return to_cqmom_result(result, hierarchy);
}

} // namespace internal

CqmomResult compute_cqmom(const CqmomCase& case_)
{
    return internal::compute_trusted_cqmom(case_.data_->lean_case);
}

} // namespace qmomkit
