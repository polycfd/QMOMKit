#include "internal/bounded_quadrature.hpp"
#include "internal/cqmom_case.hpp"
#include "internal/cqmom_validation.hpp"
#include "internal/moment_reconstruction.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/performance.hpp"
#include "internal/qmom_state_access.hpp"
#include "internal/realizability.hpp"
#include "internal/recursive_quadrature.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/realizability_policy.hpp"
#include "qmomkit/wheeler.hpp"

#include <algorithm>
#include <exception>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
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

[[nodiscard]] QmomStatus to_qmom_status(QuadratureStatus status);

[[nodiscard]] QmomResult compute_qmom_from_moments(const MomentSet& moments, const std::size_t node_count, const QmomOptions& options,
                                                   const std::optional<AbscissaBounds> abscissa_bounds)
{
    internal::performance::record_qmom_reconstruction();
    if (!abscissa_bounds.has_value())
    {
        const Quadrature1D result = compute_quadrature(moments, node_count, to_wheeler_options(options));
        return QmomResult{result.abscissas, result.weights, to_qmom_status(result.status), result.message};
    }

    const internal::BoundedQuadratureSelectionRequest request{node_count, node_count, std::nullopt, options.adaptive, false};
    const internal::BoundedQuadratureSelection selection =
        internal::select_bounded_quadrature(moments.values(), request, *abscissa_bounds, to_wheeler_options(options));
    return QmomResult{selection.quadrature.abscissas, selection.quadrature.weights, to_qmom_status(selection.quadrature.status),
                      selection.quadrature.message};
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
    }
    return QmomStatus::invalid_input;
}

// A case holds exactly one moment representation, so these read its settings
// directly. The hierarchy is reached through the moment layout rather than
// stored separately.
[[nodiscard]] const PhaseSpaceHierarchy& hierarchy(const CqmomCase& cqmom_case)
{
    return cqmom_case.moments.layout.hierarchy;
}

[[nodiscard]] const QmomOptions& qmom_options(const CqmomCase& cqmom_case)
{
    return cqmom_case.qmom_options;
}

[[nodiscard]] const ConditionalMomentOptions& conditional_options(const CqmomCase& cqmom_case)
{
    return cqmom_case.conditional_options;
}

[[nodiscard]] CqmomDiagnosticLevel diagnostic_level(const CqmomCase& cqmom_case)
{
    return cqmom_case.diagnostic_level;
}

[[nodiscard]] VandermondeSolver vandermonde_solver(const CqmomCase& cqmom_case)
{
    return cqmom_case.vandermonde_solver;
}

[[nodiscard]] const CqmomAbscissaBounds& coordinate_support(const CqmomCase& cqmom_case)
{
    return cqmom_case.abscissa_bounds;
}

[[nodiscard]] std::string diagnostic_level_name(const CqmomDiagnosticLevel level)
{
    switch (level)
    {
    case CqmomDiagnosticLevel::none:
        return "none";
    case CqmomDiagnosticLevel::basic:
        return "basic";
    case CqmomDiagnosticLevel::detailed:
        return "detailed";
    }
    return "unknown";
}

[[nodiscard]] std::string vandermonde_solver_name(const VandermondeSolver solver)
{
    switch (solver)
    {
    case VandermondeSolver::automatic:
        return "automatic";
    case VandermondeSolver::dense:
        return "dense";
    case VandermondeSolver::pivoted_qr:
        return "pivoted_qr";
    case VandermondeSolver::svd:
        return "svd";
    case VandermondeSolver::moore_penrose:
        return "moore_penrose";
    case VandermondeSolver::bjorck_pereyra:
        return "bjorck_pereyra";
    }
    return "unknown";
}
} // namespace

QmomResult compute_qmom(const MomentSet& moments, const std::size_t node_count, const QmomOptions& options)
{
    return compute_qmom_from_moments(moments, node_count, options, std::nullopt);
}

QmomResult compute_qmom(const MomentSet& moments, const std::size_t node_count, const QmomOptions& options,
                        const std::optional<AbscissaBounds> abscissa_bounds)
{
    return compute_qmom_from_moments(moments, node_count, options, abscissa_bounds);
}

QmomResult compute_qmom(const QmomCase& qmom_case)
{
    return compute_qmom(qmom_case.moments, qmom_case.node_count, qmom_case.qmom_options, qmom_case.abscissa_bounds);
}

QmomState::QmomState(std::vector<double> moments, QmomResult quadrature, const std::size_t maximum_node_count, QmomOptions qmom_options,
                     std::optional<AbscissaBounds> abscissa_bounds)
    : moments_(std::move(moments)), quadrature_(std::move(quadrature)), maximum_node_count_(maximum_node_count),
      qmom_options_(std::move(qmom_options)), abscissa_bounds_(std::move(abscissa_bounds))
{
}

namespace
{
// Whether moments hold no number density to reconstruct from, as in an empty
// cell: the reconstruction refuses an m0 below its cutoff.
[[nodiscard]] bool below_number_density_cutoff(const double m0, const QmomOptions& options)
{
    return m0 >= 0.0 && m0 < options.cutoff;
}

[[nodiscard]] RealizabilityResult policy_verdict(const MomentSet& moments, const RealizabilityPolicy& policy)
{
    return (policy.support == MomentSupport::bounded_interval && policy.interval.has_value())
               ? check_realizability(moments, *policy.interval, policy.options)
               : check_realizability(moments, policy.support, policy.options);
}

// Moments below the number-density cutoff, such as the zero or rounding-level
// moments of an empty cell, are accepted whatever the verdict, so that the
// reconstruction reports them as empty, as it does without a check.
[[nodiscard]] bool policy_accepts(const RealizabilityResult& verdict, const RealizabilityPolicy& policy, const double m0, const QmomOptions& options)
{
    return accepts_realizability_status(policy, verdict.status) || below_number_density_cutoff(m0, options);
}

// The check is optional, and the quadrature it is followed by is the one the
// state is built from, so the moments are inverted once. Moments too few for
// the node count are left to the reconstruction, which reports them as invalid
// input whether or not a check was requested.
[[nodiscard]] QmomStateCandidate qmom_state_candidate(const QmomCase& qmom_case, const MomentSet& moments, const RealizabilityPolicy* const policy)
{
    QmomStateCandidate candidate;
    if (policy != nullptr)
    {
        if (const char* error = internal::realizability_policy_error(*policy); error != nullptr)
        {
            candidate.message = std::string("QMOM realizability policy: ") + error;
            return candidate;
        }
    }
    const std::size_t node_count = qmom_case.node_count;
    const bool enough_moments = node_count > 0 && node_count <= std::numeric_limits<std::size_t>::max() / 2 && moments.size() >= 2 * node_count;
    if (policy != nullptr && policy->check_mode != RealizabilityCheckMode::None && enough_moments)
    {
        RealizabilityResult realizability = policy_verdict(moments, *policy);
        if (!policy_accepts(realizability, *policy, moments[0], qmom_case.qmom_options))
        {
            candidate.status = realizability.status == RealizabilityStatus::invalid_input ? StateCandidateStatus::invalid_input
                                                                                          : StateCandidateStatus::realizability_rejected;
            candidate.message = realizability.message;
            candidate.realizability = std::move(realizability);
            return candidate;
        }
        candidate.realizability = std::move(realizability);
    }

    QmomResult quadrature = compute_qmom(moments, node_count, qmom_case.qmom_options, qmom_case.abscissa_bounds);
    candidate.reconstruction_status = quadrature.status;
    if (!quadrature.ok())
    {
        if (quadrature.status == QmomStatus::invalid_input)
        {
            candidate.status = StateCandidateStatus::invalid_input;
        }
        else
        {
            candidate.status = below_number_density_cutoff(moments[0], qmom_case.qmom_options) ? StateCandidateStatus::empty
                                                                                               : StateCandidateStatus::reconstruction_failed;
        }
        candidate.message = std::move(quadrature.message);
        return candidate;
    }
    try
    {
        candidate.state = internal::QmomStateAccess::from_quadrature(std::move(quadrature), moments.size(), node_count, qmom_case.qmom_options,
                                                                     qmom_case.abscissa_bounds);
    }
    catch (const std::overflow_error& error)
    {
        // The stored moments, those of the accepted quadrature, can exceed
        // the range of a double although the given moments do not.
        candidate.status = StateCandidateStatus::reconstruction_failed;
        candidate.message = error.what();
        return candidate;
    }
    candidate.status = StateCandidateStatus::success;
    return candidate;
}

// The given moments must be finite throughout, including those above order
// 2n - 1, which the reconstruction does not read and the state replaces with
// the moments of its quadrature.
[[nodiscard]] QmomStateCandidate checked_qmom_state_candidate(const QmomCase& qmom_case, const MomentSet& moments,
                                                              const RealizabilityPolicy* const policy)
{
    if (!internal::all_finite(moments.values()))
    {
        QmomStateCandidate candidate;
        candidate.message = "Moment sequence contains non-finite values.";
        return candidate;
    }
    return qmom_state_candidate(qmom_case, moments, policy);
}
} // namespace

QmomStateCandidate QmomCase::state_candidate(const MomentSet& candidate_moments) const
{
    return checked_qmom_state_candidate(*this, candidate_moments, nullptr);
}

QmomStateCandidate QmomCase::state_candidate(const MomentSet& candidate_moments, const RealizabilityPolicy& policy) const
{
    return checked_qmom_state_candidate(*this, candidate_moments, &policy);
}

QmomState QmomCase::initial_state() const
{
    QmomStateCandidate candidate = checked_qmom_state_candidate(*this, moments, nullptr);
    if (!candidate.ok())
    {
        throw std::runtime_error("QMOM initial state reconstruction failed: " + candidate.message);
    }
    return std::move(*candidate.state);
}

struct CqmomPreparedCase::Data
{
    CqmomCase cqmom_case;
    internal::PreparedCqmomMetadataPtr prepared_metadata;
};

struct CqmomState::Data
{
    CqmomCase cqmom_case;
    internal::PreparedCqmomMetadataPtr prepared_metadata;
    std::vector<CqmomActiveBranchCount> active_branch_counts;
    CqmomQuadrature accepted_quadrature;
};

std::string format_cqmom_settings(const bool adaptive, const VandermondeSolver vandermonde_solver, const CqmomDiagnosticLevel diagnostic_level)
{
    std::ostringstream output;
    output << "CQMOM settings:\n"
           << "  adaptive = " << (adaptive ? "enabled" : "disabled") << '\n'
           << "  Vandermonde solver = " << vandermonde_solver_name(vandermonde_solver) << '\n'
           << "  diagnostic level = " << diagnostic_level_name(diagnostic_level) << '\n';
    return output.str();
}

CqmomPreparedCase::CqmomPreparedCase(std::unique_ptr<Data> data) : data_(std::move(data)) {}

CqmomPreparedCase::CqmomPreparedCase(CqmomCase cqmom_case)
    : CqmomPreparedCase(internal::CqmomConstructionAccess::create_prepared_case(std::move(cqmom_case)))
{
}

CqmomPreparedCase::CqmomPreparedCase(const CqmomPreparedCase& other) : data_(other.data_ == nullptr ? nullptr : std::make_unique<Data>(*other.data_))
{
}

CqmomPreparedCase& CqmomPreparedCase::operator=(const CqmomPreparedCase& other)
{
    if (this != &other)
    {
        data_ = other.data_ == nullptr ? nullptr : std::make_unique<Data>(*other.data_);
    }
    return *this;
}

CqmomPreparedCase::CqmomPreparedCase(CqmomPreparedCase&&) noexcept = default;

CqmomPreparedCase& CqmomPreparedCase::operator=(CqmomPreparedCase&&) noexcept = default;

CqmomPreparedCase::~CqmomPreparedCase() = default;

const PhaseSpaceHierarchy& CqmomPreparedCase::hierarchy() const noexcept
{
    return qmomkit::hierarchy(data_->cqmom_case);
}

const QmomOptions& CqmomPreparedCase::qmom_options() const noexcept
{
    return qmomkit::qmom_options(data_->cqmom_case);
}

const ConditionalMomentOptions& CqmomPreparedCase::conditional_options() const noexcept
{
    return qmomkit::conditional_options(data_->cqmom_case);
}

CqmomDiagnosticLevel CqmomPreparedCase::diagnostic_level() const noexcept
{
    return qmomkit::diagnostic_level(data_->cqmom_case);
}

VandermondeSolver CqmomPreparedCase::vandermonde_solver() const noexcept
{
    return qmomkit::vandermonde_solver(data_->cqmom_case);
}

const CqmomAbscissaBounds& CqmomPreparedCase::abscissa_bounds() const noexcept
{
    return qmomkit::coordinate_support(data_->cqmom_case);
}

namespace
{
// Both cases build their states here. The metadata handle is passed on to the
// new state, so its later steps reuse the branch structure; the raw pointer is
// what the reconstruction call takes. The moment data must already be valid
// for the hierarchy. The checks are optional, and the quadrature they are
// followed by is the one the state is built from, so the moments are
// reconstructed once.
[[nodiscard]] CqmomStateCandidate cqmom_state_candidate(const CqmomCase& cqmom_case,
                                                        const internal::PreparedCqmomMetadataPtr& prepared_metadata_handle,
                                                        const std::vector<std::optional<RealizabilityPolicy>>* const coordinate_policies)
{
    CqmomStateCandidate candidate;
    if (coordinate_policies != nullptr)
    {
        // As in the moment-source update, the policies are validated before
        // any coordinate is checked.
        const CqmomMomentLayout& layout = cqmom_case.moments.layout;
        const std::size_t dimensions = layout.hierarchy.number_of_dimensions();
        if (std::string error = internal::coordinate_realizability_policies_error(*coordinate_policies, dimensions); !error.empty())
        {
            candidate.message = std::move(error);
            return candidate;
        }
        for (std::size_t coordinate = 0; coordinate < std::min(dimensions, coordinate_policies->size()); ++coordinate)
        {
            const std::optional<RealizabilityPolicy>& policy = (*coordinate_policies)[coordinate];
            if (!policy.has_value() || policy->check_mode == RealizabilityCheckMode::None)
            {
                continue;
            }
            const MomentSet marginal(internal::extract_coordinate_marginal(cqmom_case.moments.values, layout, coordinate));
            RealizabilityResult realizability = policy_verdict(marginal, *policy);
            if (!policy_accepts(realizability, *policy, marginal[0], cqmom_case.qmom_options))
            {
                candidate.status = realizability.status == RealizabilityStatus::invalid_input ? StateCandidateStatus::invalid_input
                                                                                              : StateCandidateStatus::realizability_rejected;
                candidate.message = "Coordinate " + std::to_string(coordinate) + ": " + realizability.message;
                candidate.realizability = std::move(realizability);
                candidate.rejected_coordinate = coordinate;
                return candidate;
            }
        }
    }

    std::vector<CqmomActiveBranchCount> active_counts;
    CqmomResult accepted = internal::compute_cqmom_internal(cqmom_case, *prepared_metadata_handle, {}, false, &active_counts);
    candidate.reconstruction_status = accepted.status;
    if (!accepted.ok())
    {
        if (accepted.status == CqmomStatus::invalid_input)
        {
            candidate.status = StateCandidateStatus::invalid_input;
        }
        else
        {
            const double m0 = cqmom_case.moments.moment(MultiIndex(cqmom_case.moments.layout.hierarchy.number_of_dimensions(), 0));
            candidate.status =
                below_number_density_cutoff(m0, cqmom_case.qmom_options) ? StateCandidateStatus::empty : StateCandidateStatus::reconstruction_failed;
        }
        candidate.message = std::move(accepted.message);
        return candidate;
    }
    try
    {
        CqmomQuadrature accepted_quadrature = std::move(accepted.quadrature);
        CqmomCase canonical_cqmom_case =
            internal::with_moments(cqmom_case, internal::reconstruct_trusted_cqmom_moments(accepted_quadrature, cqmom_case.moments.layout));
        candidate.state = internal::CqmomConstructionAccess::create_state(std::move(canonical_cqmom_case), std::move(active_counts),
                                                                          std::move(accepted_quadrature), prepared_metadata_handle);
    }
    catch (const std::overflow_error& error)
    {
        // The stored moments, those of the accepted quadrature, can exceed
        // the range of a double although the given moments do not.
        candidate.status = StateCandidateStatus::reconstruction_failed;
        candidate.message = error.what();
        return candidate;
    }
    candidate.status = StateCandidateStatus::success;
    return candidate;
}

[[nodiscard]] CqmomState initial_cqmom_state(const CqmomCase& cqmom_case, const internal::PreparedCqmomMetadataPtr& prepared_metadata_handle)
{
    CqmomStateCandidate candidate = cqmom_state_candidate(cqmom_case, prepared_metadata_handle, nullptr);
    if (!candidate.ok())
    {
        throw std::runtime_error("CQMOM initial state quadrature reconstruction failed: " + candidate.message);
    }
    return std::move(*candidate.state);
}

// A candidate's moments come from outside the case, so they are validated
// against its layout before the reconstruction, which trusts them. The layout
// of a prepared case is already the required one of its hierarchy.
[[nodiscard]] CqmomStateCandidate prepared_cqmom_state_candidate(const CqmomCase& base_cqmom_case,
                                                                 const internal::PreparedCqmomMetadataPtr& prepared_metadata_handle,
                                                                 const CqmomMomentData& candidate_moments,
                                                                 const std::vector<std::optional<RealizabilityPolicy>>* const coordinate_policies)
{
    try
    {
        const internal::CqmomMomentDataValidationMessages messages{
            "CQMOM moment data layout must match the required CQMOM moment layout of the base hierarchy.",
            "CQMOM moment data value count does not match its required CQMOM moment layout.", "CQMOM moment data values must be finite."};
        internal::validate_cqmom_moment_layout_matches_expected(candidate_moments.layout, base_cqmom_case.moments.layout, messages.fixed_layout);
        internal::validate_cqmom_moment_values_match_layout(candidate_moments, messages);
    }
    catch (const std::invalid_argument& error)
    {
        CqmomStateCandidate candidate;
        candidate.message = error.what();
        return candidate;
    }
    return cqmom_state_candidate(internal::with_moments(base_cqmom_case, candidate_moments), prepared_metadata_handle, coordinate_policies);
}
} // namespace

CqmomState CqmomCase::initial_state() const
{
    // The branch structure is built directly rather than through a
    // CqmomPreparedCase, which would first copy this case.
    return initial_cqmom_state(*this, internal::create_prepared_cqmom_metadata(moments.layout, qmom_options.adaptive));
}

CqmomStateCandidate CqmomPreparedCase::state_candidate(const CqmomMomentData& candidate_moments) const
{
    return prepared_cqmom_state_candidate(data_->cqmom_case, data_->prepared_metadata, candidate_moments, nullptr);
}

CqmomStateCandidate CqmomPreparedCase::state_candidate(const CqmomMomentData& candidate_moments,
                                                       const std::vector<std::optional<RealizabilityPolicy>>& coordinate_policies) const
{
    return prepared_cqmom_state_candidate(data_->cqmom_case, data_->prepared_metadata, candidate_moments, &coordinate_policies);
}

CqmomState CqmomPreparedCase::initial_state() const
{
    return initial_cqmom_state(data_->cqmom_case, data_->prepared_metadata);
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
    return qmomkit::hierarchy(data_->cqmom_case);
}

const CqmomQuadrature& CqmomState::quadrature() const noexcept
{
    return data_->accepted_quadrature;
}

const QmomOptions& CqmomState::qmom_options() const noexcept
{
    return qmomkit::qmom_options(data_->cqmom_case);
}

const ConditionalMomentOptions& CqmomState::conditional_options() const noexcept
{
    return qmomkit::conditional_options(data_->cqmom_case);
}

CqmomDiagnosticLevel CqmomState::diagnostic_level() const noexcept
{
    return qmomkit::diagnostic_level(data_->cqmom_case);
}

VandermondeSolver CqmomState::vandermonde_solver() const noexcept
{
    return qmomkit::vandermonde_solver(data_->cqmom_case);
}

const CqmomAbscissaBounds& CqmomState::abscissa_bounds() const noexcept
{
    return qmomkit::coordinate_support(data_->cqmom_case);
}

const std::vector<CqmomActiveBranchCount>& CqmomState::active_branch_counts() const noexcept
{
    return data_->active_branch_counts;
}

std::vector<MultiIndex> CqmomState::moment_indices() const
{
    return data_->cqmom_case.moments.layout.required_indices.to_multi_indices();
}

const std::vector<double>& CqmomState::moment_values() const noexcept
{
    return data_->cqmom_case.moments.values;
}

std::optional<double> CqmomState::moment(const MultiIndexView index) const
{
    if (!data_->cqmom_case.moments.contains_required_index(index))
    {
        return std::nullopt;
    }
    return data_->cqmom_case.moments.moment(index);
}

namespace internal
{
namespace
{
void validate_prepared_metadata(const PreparedCqmomMetadata& metadata, const CqmomCase& cqmom_case)
{
    if (!prepared_cqmom_metadata_is_compatible(metadata, cqmom_case.moments.layout.hierarchy, cqmom_case.qmom_options.adaptive,
                                               &cqmom_case.moments.layout))
    {
        throw std::invalid_argument("Prepared CQMOM metadata is incompatible with the target CQMOM case.");
    }
}

[[nodiscard]] PreparedCqmomMetadataPtr make_prepared_metadata(const CqmomCase& cqmom_case, PreparedCqmomMetadataPtr prepared_metadata)
{
    if (prepared_metadata != nullptr)
    {
        validate_prepared_metadata(*prepared_metadata, cqmom_case);
        return prepared_metadata;
    }
    return create_prepared_cqmom_metadata(cqmom_case.moments.layout, cqmom_case.qmom_options.adaptive);
}
} // namespace

CqmomPreparedCase CqmomConstructionAccess::create_prepared_case(CqmomCase cqmom_case, PreparedCqmomMetadataPtr prepared_metadata)
{
    prepared_metadata = make_prepared_metadata(cqmom_case, std::move(prepared_metadata));
    return CqmomPreparedCase(std::make_unique<CqmomPreparedCase::Data>(CqmomPreparedCase::Data{std::move(cqmom_case), std::move(prepared_metadata)}));
}

CqmomState CqmomConstructionAccess::create_state(CqmomCase cqmom_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 CqmomQuadrature accepted_quadrature, PreparedCqmomMetadataPtr prepared_metadata)
{
    prepared_metadata = make_prepared_metadata(cqmom_case, std::move(prepared_metadata));
    return CqmomState(std::make_unique<CqmomState::Data>(
        CqmomState::Data{std::move(cqmom_case), std::move(prepared_metadata), std::move(active_branch_counts), std::move(accepted_quadrature)}));
}

const CqmomCase& CqmomConstructionAccess::case_definition(const CqmomPreparedCase& cqmom_prepared_case)
{
    return cqmom_prepared_case.data_->cqmom_case;
}

const CqmomCase& CqmomConstructionAccess::state_definition(const CqmomState& state)
{
    return state.data_->cqmom_case;
}

const std::vector<CqmomActiveBranchCount>& CqmomConstructionAccess::state_active_branch_counts(const CqmomState& state)
{
    return state.data_->active_branch_counts;
}

const PreparedCqmomMetadataPtr& CqmomConstructionAccess::case_prepared_metadata(const CqmomPreparedCase& cqmom_prepared_case)
{
    return cqmom_prepared_case.data_->prepared_metadata;
}

const PreparedCqmomMetadataPtr& CqmomConstructionAccess::state_prepared_metadata(const CqmomState& state)
{
    return state.data_->prepared_metadata;
}

CqmomQuadrature CqmomConstructionAccess::create_quadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights,
                                                           std::vector<std::vector<double>> abscissas, const bool adaptive_branch_counts_enabled,
                                                           const bool rectangular_storage, ConditioningPathTable terminal_paths)
{
    return CqmomQuadrature(std::move(hierarchy), std::move(weights), std::move(abscissas), adaptive_branch_counts_enabled, rectangular_storage,
                           std::move(terminal_paths));
}

} // namespace internal

CqmomQuadrature::CqmomQuadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights, std::vector<std::vector<double>> abscissas,
                                 const bool adaptive_branch_counts_enabled, const bool rectangular_storage, ConditioningPathTable terminal_paths)
    : hierarchy_(std::move(hierarchy)), weights_(std::move(weights)), abscissas_(std::move(abscissas)), terminal_paths_(std::move(terminal_paths)),
      adaptive_branch_counts_enabled_(adaptive_branch_counts_enabled), rectangular_storage_(rectangular_storage)
{
}

double CqmomQuadrature::weight(const MultiIndexView node_index) const
{
    if (!rectangular_storage_)
    {
        throw std::logic_error("Ragged CQMOM quadrature cannot be addressed by rectangular node multi-index.");
    }
    return weights_.at(hierarchy_.flatten_node_index(node_index));
}

double CqmomQuadrature::abscissa(const std::size_t dimension, const MultiIndexView node_index) const
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

std::vector<double> CqmomQuadrature::node(const MultiIndexView node_index) const
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
    CqmomDiagnostics converted{};
    converted.conditioning_level = diagnostics.conditioning_level;
    converted.parent_path = diagnostics.path;
    converted.phase_space_dimension = diagnostics.phase_space_dimension;
    converted.requested_node_count = diagnostics.requested_node_count;
    converted.actual_node_count = diagnostics.actual_node_count;
    converted.adaptive_branch_counts_enabled = diagnostics.adaptive_branch_counts_enabled;
    converted.adaptive_node_count_reduced = diagnostics.adaptive_node_count_reduced;
    converted.qmom_status = to_qmom_status(diagnostics.wheeler_status);
    converted.conditional_moment_status = diagnostics.conditional_moment_status;
    converted.conditional_moment_diagnostics = diagnostics.conditional_moment_diagnostics;
    converted.conditional_vandermonde_requested_solver = diagnostics.conditional_vandermonde_requested_solver;
    converted.conditional_vandermonde_selected_solver = diagnostics.conditional_vandermonde_selected_solver;
    converted.message = diagnostics.message;
    converted.bounded_support_checked = diagnostics.bounded_support_checked;
    converted.bounded_support_interval_present = diagnostics.bounded_support_interval_present;
    converted.bounded_support_interval = diagnostics.bounded_support_interval;
    converted.bounded_previous_active_node_count = diagnostics.bounded_previous_active_node_count;
    converted.bounded_candidate_node_count = diagnostics.bounded_candidate_node_count;
    converted.bounded_candidate_attempt_type = diagnostics.bounded_candidate_attempt_type;
    converted.bounded_accepted_fallback_node_count = diagnostics.bounded_accepted_fallback_node_count;
    converted.bounded_returned_nodes = diagnostics.bounded_returned_nodes;
    converted.bounded_returned_weights = diagnostics.bounded_returned_weights;
    return converted;
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

[[nodiscard]] ConditioningPathTable terminal_paths_for_public_quadrature(const internal::RecursiveQuadrature& quadrature)
{
    if (!quadrature.terminal_paths().empty() || !quadrature.has_rectangular_storage())
    {
        return quadrature.terminal_paths();
    }

    // Built flat: one table allocation rather than one vector per path.
    const PhaseSpaceHierarchy& hierarchy = quadrature.hierarchy();
    const std::size_t levels = hierarchy.number_of_dimensions();
    std::vector<std::size_t> node_indices;
    node_indices.reserve(quadrature.weights().size() * levels);
    for (const MultiIndex& node_index : enumerate_multi_indices(hierarchy.node_counts()))
    {
        for (std::size_t level = 0; level < levels; ++level)
        {
            node_indices.push_back(node_index[hierarchy.conditioning_order()[level]]);
        }
    }
    return ConditioningPathTable(levels, std::move(node_indices));
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
    return CqmomResult{internal::CqmomConstructionAccess::create_quadrature(
                           result.quadrature.hierarchy(), result.quadrature.weights(), result.quadrature.abscissas(),
                           result.quadrature.adaptive_branch_counts_enabled(), result.quadrature.has_rectangular_storage(),
                           terminal_paths_for_public_quadrature(result.quadrature)),
                       to_cqmom_status(result.status), result.message, to_cqmom_diagnostics(result.diagnostics),
                       failed_phase_space_dimension(result, hierarchy)};
}

[[nodiscard]] CqmomResult invalid_cqmom_hierarchy_result(const PhaseSpaceHierarchy& hierarchy, const bool adaptive)
{
    return CqmomResult{internal::CqmomConstructionAccess::create_quadrature(hierarchy, {}, {}, adaptive, !adaptive, {}),
                       CqmomStatus::invalid_input,
                       "CQMOM node count is too large to derive the required moment extent.",
                       {},
                       0};
}
} // namespace

CqmomResult compute_cqmom(const CqmomCase& cqmom_case)
{
    const PhaseSpaceHierarchy& hierarchy = cqmom_case.moments.layout.hierarchy;
    if (!internal::cqmom_required_moment_extents_representable(hierarchy))
    {
        return invalid_cqmom_hierarchy_result(hierarchy, cqmom_case.qmom_options.adaptive);
    }
    const internal::RecursiveQuadratureResult result =
        internal::compute_recursive_quadrature_from_moment_data({hierarchy,
                                                                 cqmom_case.moments,
                                                                 cqmom_case.conditional_options,
                                                                 to_wheeler_options(cqmom_case.qmom_options),
                                                                 cqmom_case.diagnostic_level,
                                                                 cqmom_case.vandermonde_solver,
                                                                 {},
                                                                 false,
                                                                 cqmom_case.abscissa_bounds});

    return to_cqmom_result(result, hierarchy);
}

namespace internal
{

CqmomResult compute_cqmom_internal(const CqmomCase& cqmom_case, const PreparedCqmomMetadata& prepared_metadata,
                                   std::vector<CqmomActiveBranchCount> active_branch_counts, const bool allow_active_branch_count_growth,
                                   std::vector<CqmomActiveBranchCount>* const accepted_active_branch_counts)
{
    const PhaseSpaceHierarchy& hierarchy = cqmom_case.moments.layout.hierarchy;
    if (!cqmom_required_moment_extents_representable(hierarchy))
    {
        return invalid_cqmom_hierarchy_result(hierarchy, cqmom_case.qmom_options.adaptive);
    }
    // The traversal already knows the node count it accepts for each branch,
    // so a state takes its active counts from the reconstruction.
    const bool record_active_branch_counts = accepted_active_branch_counts != nullptr;
    const auto finish = [&](RecursiveQuadratureResult& result)
    {
        CqmomResult converted = to_cqmom_result(result, hierarchy);
        if (record_active_branch_counts && converted.ok())
        {
            *accepted_active_branch_counts = std::move(result.active_branch_counts);
            if (cqmom_case.qmom_options.adaptive)
            {
                performance::record_active_count_derivation();
            }
        }
        return converted;
    };
    RecursiveQuadratureResult result = compute_recursive_quadrature_from_trusted_moment_data(
        {hierarchy, cqmom_case.moments, cqmom_case.conditional_options, to_wheeler_options(cqmom_case.qmom_options), cqmom_case.diagnostic_level,
         cqmom_case.vandermonde_solver, std::move(active_branch_counts), allow_active_branch_count_growth, cqmom_case.abscissa_bounds,
         &prepared_metadata, record_active_branch_counts});

    return finish(result);
}

} // namespace internal

CqmomResult compute_cqmom(const CqmomPreparedCase& cqmom_prepared_case)
{
    return internal::compute_cqmom_internal(cqmom_prepared_case.data_->cqmom_case, *cqmom_prepared_case.data_->prepared_metadata);
}

} // namespace qmomkit
