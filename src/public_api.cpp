#include "internal/bounded_quadrature.hpp"
#include "internal/cqmom_active_counts.hpp"
#include "internal/cqmom_case.hpp"
#include "internal/cqmom_validation.hpp"
#include "internal/moment_reconstruction.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/performance.hpp"
#include "internal/qmom_state_access.hpp"
#include "internal/reconstruction_scaling.hpp"
#include "internal/recursive_quadrature.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/wheeler.hpp"

#include <exception>
#include <optional>
#include <sstream>
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

[[nodiscard]] QmomStatus to_qmom_status(QuadratureStatus status);

[[nodiscard]] QmomResult compute_qmom_unscaled(const MomentSet& moments, const std::size_t node_count, const QmomOptions& options,
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

[[nodiscard]] QmomResult compute_qmom_with_automatic_scaling(const MomentSet& moments, const std::size_t node_count, const QmomOptions& options,
                                                             const std::optional<AbscissaBounds> abscissa_bounds)
{
    if (!options.automatic_scaling)
    {
        return compute_qmom_unscaled(moments, node_count, options, abscissa_bounds);
    }

    const std::optional<internal::QmomScalingContext> scaling = internal::make_qmom_scaling_context(moments.values(), abscissa_bounds);
    if (!scaling.has_value())
    {
        return compute_qmom_unscaled(moments, node_count, options, abscissa_bounds);
    }

    QmomResult scaled_result = compute_qmom_unscaled(MomentSet{scaling->moments}, node_count, options, scaling->abscissa_bounds);
    if (!internal::restore_abscissas(scaled_result.abscissas, scaling->scale))
    {
        return compute_qmom_unscaled(moments, node_count, options, abscissa_bounds);
    }
    return scaled_result;
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

[[nodiscard]] bool adaptive_branch_counts(const CqmomCase& cqmom_case)
{
    return cqmom_case.adaptive_branch_counts;
}

[[nodiscard]] CqmomDiagnosticLevel diagnostic_level(const CqmomCase& cqmom_case)
{
    return cqmom_case.diagnostic_level;
}

[[nodiscard]] VandermondeSolver vandermonde_solver(const CqmomCase& cqmom_case)
{
    return cqmom_case.vandermonde_solver;
}

[[nodiscard]] const CqmomCoordinateSupport& coordinate_support(const CqmomCase& cqmom_case)
{
    return cqmom_case.coordinate_support;
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
    return compute_qmom_with_automatic_scaling(moments, node_count, options, std::nullopt);
}

QmomResult compute_qmom(const MomentSet& moments, const std::size_t node_count, const QmomOptions& options,
                        const std::optional<AbscissaBounds> abscissa_bounds)
{
    return compute_qmom_with_automatic_scaling(moments, node_count, options, abscissa_bounds);
}

QmomResult compute_qmom(const QmomCase& qmom_case)
{
    return compute_qmom(qmom_case.moments, qmom_case.node_count, qmom_case.options, qmom_case.abscissa_bounds);
}

QmomState::QmomState(std::vector<double> moments, QmomResult quadrature, const std::size_t maximum_node_count, QmomOptions qmom_options,
                     std::optional<AbscissaBounds> abscissa_bounds)
    : moments_(std::move(moments)), quadrature_(std::move(quadrature)), maximum_node_count_(maximum_node_count),
      qmom_options_(std::move(qmom_options)), abscissa_bounds_(std::move(abscissa_bounds))
{
}

QmomState QmomCase::initial_state() const
{
    QmomResult quadrature = compute_qmom(*this);
    if (!quadrature.ok())
    {
        throw std::runtime_error("QMOM initial state reconstruction failed: " + quadrature.message);
    }
    return internal::QmomStateAccess::from_quadrature(std::move(quadrature), moments.size(), node_count, options, abscissa_bounds);
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
    std::optional<CqmomQuadrature> accepted_quadrature;
};

std::string format_cqmom_settings(const bool adaptive_branch_counts, const VandermondeSolver vandermonde_solver,
                                  const CqmomDiagnosticLevel diagnostic_level)
{
    std::ostringstream output;
    output << "CQMOM settings:\n"
           << "  adaptive branch counts = " << (adaptive_branch_counts ? "enabled" : "disabled") << '\n'
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

bool CqmomPreparedCase::adaptive_branch_counts() const noexcept
{
    return qmomkit::adaptive_branch_counts(data_->cqmom_case);
}

CqmomDiagnosticLevel CqmomPreparedCase::diagnostic_level() const noexcept
{
    return qmomkit::diagnostic_level(data_->cqmom_case);
}

VandermondeSolver CqmomPreparedCase::vandermonde_solver() const noexcept
{
    return qmomkit::vandermonde_solver(data_->cqmom_case);
}

const CqmomCoordinateSupport& CqmomPreparedCase::coordinate_support() const noexcept
{
    return qmomkit::coordinate_support(data_->cqmom_case);
}

namespace
{
// Both cases build their initial state here. The metadata handle is passed on
// to the new state, so its later steps reuse the branch structure; the raw
// pointer is what the reconstruction call takes.
CqmomState initial_cqmom_state(const CqmomCase& cqmom_case, const internal::PreparedCqmomMetadataPtr& prepared_metadata_handle)
{
    const internal::PreparedCqmomMetadata* prepared_metadata = prepared_metadata_handle.get();
    CqmomResult accepted = internal::compute_cqmom_internal(cqmom_case, {}, false, prepared_metadata, true);
    if (!accepted.ok())
    {
        throw std::runtime_error("CQMOM initial state quadrature reconstruction failed: " + accepted.message);
    }
    CqmomQuadrature accepted_quadrature = std::move(accepted.quadrature);
    std::vector<CqmomActiveBranchCount> active_counts = internal::active_branch_counts_from_quadrature(accepted_quadrature);
    CqmomCase canonical_cqmom_case =
        internal::with_moments(cqmom_case, internal::reconstruct_trusted_cqmom_moments(accepted_quadrature, cqmom_case.moments.layout));
    return internal::CqmomConstructionAccess::create_state(std::move(canonical_cqmom_case), std::move(active_counts), std::move(accepted_quadrature),
                                                           prepared_metadata_handle);
}
} // namespace

CqmomState CqmomCase::initial_state() const
{
    // The branch structure is built directly rather than through a
    // CqmomPreparedCase, which would first copy this case.
    return initial_cqmom_state(*this, internal::create_prepared_cqmom_metadata(moments.layout, adaptive_branch_counts));
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

const CqmomQuadrature& CqmomState::quadrature() const
{
    if (!data_->accepted_quadrature.has_value())
    {
        throw std::logic_error("CQMOM state has no accepted quadrature.");
    }
    return *data_->accepted_quadrature;
}

const QmomOptions& CqmomState::qmom_options() const noexcept
{
    return qmomkit::qmom_options(data_->cqmom_case);
}

const ConditionalMomentOptions& CqmomState::conditional_options() const noexcept
{
    return qmomkit::conditional_options(data_->cqmom_case);
}

bool CqmomState::adaptive_branch_counts() const noexcept
{
    return qmomkit::adaptive_branch_counts(data_->cqmom_case);
}

CqmomDiagnosticLevel CqmomState::diagnostic_level() const noexcept
{
    return qmomkit::diagnostic_level(data_->cqmom_case);
}

VandermondeSolver CqmomState::vandermonde_solver() const noexcept
{
    return qmomkit::vandermonde_solver(data_->cqmom_case);
}

const CqmomCoordinateSupport& CqmomState::coordinate_support() const noexcept
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
    if (!prepared_cqmom_metadata_is_compatible(metadata, cqmom_case.moments.layout.hierarchy, cqmom_case.adaptive_branch_counts,
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
    return create_prepared_cqmom_metadata(cqmom_case.moments.layout, cqmom_case.adaptive_branch_counts);
}
} // namespace

CqmomPreparedCase CqmomConstructionAccess::create_prepared_case(CqmomCase cqmom_case, PreparedCqmomMetadataPtr prepared_metadata)
{
    prepared_metadata = make_prepared_metadata(cqmom_case, std::move(prepared_metadata));
    return CqmomPreparedCase(std::make_unique<CqmomPreparedCase::Data>(CqmomPreparedCase::Data{std::move(cqmom_case), std::move(prepared_metadata)}));
}

CqmomState CqmomConstructionAccess::create_state(CqmomCase cqmom_case, PreparedCqmomMetadataPtr prepared_metadata)
{
    return create_state(std::move(cqmom_case), {}, std::nullopt, std::move(prepared_metadata));
}

CqmomState CqmomConstructionAccess::create_state(CqmomCase cqmom_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 PreparedCqmomMetadataPtr prepared_metadata)
{
    return create_state(std::move(cqmom_case), std::move(active_branch_counts), std::nullopt, std::move(prepared_metadata));
}

CqmomState CqmomConstructionAccess::create_state(CqmomCase cqmom_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 std::optional<CqmomQuadrature> accepted_quadrature, PreparedCqmomMetadataPtr prepared_metadata)
{
    prepared_metadata = make_prepared_metadata(cqmom_case, std::move(prepared_metadata));
    return CqmomState(std::make_unique<CqmomState::Data>(
        CqmomState::Data{std::move(cqmom_case), std::move(prepared_metadata), std::move(active_branch_counts), std::move(accepted_quadrature)}));
}

CqmomState CqmomConstructionAccess::create_state(CqmomCase cqmom_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 CqmomQuadrature&& accepted_quadrature, PreparedCqmomMetadataPtr prepared_metadata)
{
    return create_state(std::move(cqmom_case), std::move(active_branch_counts), std::optional<CqmomQuadrature>{std::move(accepted_quadrature)},
                        std::move(prepared_metadata));
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

const std::optional<CqmomQuadrature>& CqmomConstructionAccess::state_accepted_quadrature(const CqmomState& state)
{
    return state.data_->accepted_quadrature;
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

[[nodiscard]] CqmomDiagnostics to_cqmom_diagnostics(const internal::RecursiveQuadratureBranchDiagnostics& diagnostics,
                                                    const std::vector<double>* const scales = nullptr)
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
    if (scales != nullptr && converted.phase_space_dimension < scales->size())
    {
        const double scale = scales->at(converted.phase_space_dimension);
        if (converted.bounded_support_interval_present && !internal::restore_bounds(converted.bounded_support_interval, scale))
        {
            converted.bounded_support_interval_present = false;
        }
        if (!internal::restore_abscissas(converted.bounded_returned_nodes, scale))
        {
            converted.bounded_returned_nodes.clear();
        }
    }
    return converted;
}

[[nodiscard]] std::vector<CqmomDiagnostics>
to_cqmom_diagnostics(const std::vector<internal::RecursiveQuadratureBranchDiagnostics>& branch_diagnostics,
                     const std::vector<double>* const scales = nullptr)
{
    std::vector<CqmomDiagnostics> diagnostics;
    diagnostics.reserve(branch_diagnostics.size());
    for (const internal::RecursiveQuadratureBranchDiagnostics& diagnostic : branch_diagnostics)
    {
        diagnostics.push_back(to_cqmom_diagnostics(diagnostic, scales));
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

[[nodiscard]] CqmomResult to_cqmom_result(const internal::RecursiveQuadratureResult& result, const PhaseSpaceHierarchy& hierarchy,
                                          const std::vector<double>* const scales = nullptr)
{
    std::vector<std::vector<double>> abscissas = result.quadrature.abscissas();
    if (scales != nullptr)
    {
        (void)internal::restore_cqmom_abscissas(abscissas, *scales);
    }
    return CqmomResult{internal::CqmomConstructionAccess::create_quadrature(result.quadrature.hierarchy(), result.quadrature.weights(),
                                                                            std::move(abscissas), result.quadrature.adaptive_branch_counts_enabled(),
                                                                            result.quadrature.has_rectangular_storage(),
                                                                            terminal_paths_for_public_quadrature(result.quadrature)),
                       to_cqmom_status(result.status), result.message, to_cqmom_diagnostics(result.diagnostics, scales),
                       failed_phase_space_dimension(result, hierarchy)};
}

[[nodiscard]] CqmomResult invalid_cqmom_hierarchy_result(const PhaseSpaceHierarchy& hierarchy, const bool adaptive_branch_counts)
{
    return CqmomResult{internal::CqmomConstructionAccess::create_quadrature(hierarchy, {}, {}, adaptive_branch_counts, !adaptive_branch_counts, {}),
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
        return invalid_cqmom_hierarchy_result(hierarchy, cqmom_case.adaptive_branch_counts);
    }
    if (cqmom_case.qmom_options.automatic_scaling)
    {
        const std::optional<internal::CqmomScalingContext> scaling =
            internal::make_cqmom_scaling_context(cqmom_case.moments, cqmom_case.coordinate_support);
        if (scaling.has_value())
        {
            const internal::RecursiveQuadratureResult result =
                internal::compute_recursive_quadrature_from_moment_data({hierarchy,
                                                                         scaling->moments,
                                                                         cqmom_case.conditional_options,
                                                                         to_wheeler_options(cqmom_case.qmom_options),
                                                                         cqmom_case.adaptive_branch_counts,
                                                                         cqmom_case.diagnostic_level,
                                                                         cqmom_case.vandermonde_solver,
                                                                         {},
                                                                         false,
                                                                         scaling->coordinate_support});
            if (internal::can_restore_cqmom_result(result, scaling->scales))
            {
                return to_cqmom_result(result, hierarchy, &scaling->scales);
            }
        }
    }
    const internal::RecursiveQuadratureResult result =
        internal::compute_recursive_quadrature_from_moment_data({hierarchy,
                                                                 cqmom_case.moments,
                                                                 cqmom_case.conditional_options,
                                                                 to_wheeler_options(cqmom_case.qmom_options),
                                                                 cqmom_case.adaptive_branch_counts,
                                                                 cqmom_case.diagnostic_level,
                                                                 cqmom_case.vandermonde_solver,
                                                                 {},
                                                                 false,
                                                                 cqmom_case.coordinate_support});

    return to_cqmom_result(result, hierarchy);
}

namespace internal
{

CqmomResult compute_cqmom_internal(const CqmomCase& cqmom_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                   const bool allow_active_branch_count_growth, const PreparedCqmomMetadata* const prepared_metadata,
                                   const bool require_prepared_metadata)
{
    const PhaseSpaceHierarchy& hierarchy = cqmom_case.moments.layout.hierarchy;
    if (!cqmom_required_moment_extents_representable(hierarchy))
    {
        return invalid_cqmom_hierarchy_result(hierarchy, cqmom_case.adaptive_branch_counts);
    }
    if (require_prepared_metadata && !cqmom_case.adaptive_branch_counts && prepared_metadata == nullptr)
    {
        throw std::logic_error("Fixed normal CQMOM reconstruction requires prepared metadata.");
    }
    if (cqmom_case.qmom_options.automatic_scaling)
    {
        const std::optional<CqmomScalingContext> scaling = make_cqmom_scaling_context(cqmom_case.moments, cqmom_case.coordinate_support);
        if (scaling.has_value())
        {
            const RecursiveQuadratureResult result = compute_recursive_quadrature_from_trusted_moment_data(
                {hierarchy, scaling->moments, cqmom_case.conditional_options, to_wheeler_options(cqmom_case.qmom_options),
                 cqmom_case.adaptive_branch_counts, cqmom_case.diagnostic_level, cqmom_case.vandermonde_solver, active_branch_counts,
                 allow_active_branch_count_growth, scaling->coordinate_support, prepared_metadata, require_prepared_metadata});
            if (can_restore_cqmom_result(result, scaling->scales))
            {
                return to_cqmom_result(result, hierarchy, &scaling->scales);
            }
        }
    }
    const RecursiveQuadratureResult result = compute_recursive_quadrature_from_trusted_moment_data(
        {hierarchy, cqmom_case.moments, cqmom_case.conditional_options, to_wheeler_options(cqmom_case.qmom_options),
         cqmom_case.adaptive_branch_counts, cqmom_case.diagnostic_level, cqmom_case.vandermonde_solver, std::move(active_branch_counts),
         allow_active_branch_count_growth, cqmom_case.coordinate_support, prepared_metadata, require_prepared_metadata});

    return to_cqmom_result(result, hierarchy);
}

} // namespace internal

CqmomResult compute_cqmom(const CqmomPreparedCase& cqmom_prepared_case)
{
    const internal::PreparedCqmomMetadata* prepared_metadata = cqmom_prepared_case.data_->prepared_metadata.get();
    return internal::compute_cqmom_internal(cqmom_prepared_case.data_->cqmom_case, {}, false, prepared_metadata, true);
}

} // namespace qmomkit
