#include "internal/cqmom_active_counts.hpp"
#include "internal/cqmom_case.hpp"
#include "internal/cqmom_validation.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/performance.hpp"
#include "internal/qmom_state_access.hpp"
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
#include <variant>
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
    }
    return QmomStatus::invalid_input;
}

template <class... Ts> struct Overloaded : Ts...
{
    using Ts::operator()...;
};
template <class... Ts> Overloaded(Ts...) -> Overloaded<Ts...>;

using CqmomCaseVariant = std::variant<CqmomDenseCase, CqmomLeanCase>;

[[nodiscard]] const PhaseSpaceHierarchy& hierarchy(const CqmomCaseVariant& case_)
{
    return std::visit(Overloaded{[](const CqmomDenseCase& dense_case) -> const PhaseSpaceHierarchy& { return dense_case.hierarchy; },
                                 [](const CqmomLeanCase& lean_case) -> const PhaseSpaceHierarchy& { return lean_case.moments.layout.hierarchy; }},
                      case_);
}

[[nodiscard]] const QmomOptions& qmom_options(const CqmomCaseVariant& case_)
{
    return std::visit(Overloaded{[](const CqmomDenseCase& dense_case) -> const QmomOptions& { return dense_case.qmom_options; },
                                 [](const CqmomLeanCase& lean_case) -> const QmomOptions& { return lean_case.qmom_options; }},
                      case_);
}

[[nodiscard]] const ConditionalMomentOptions& conditional_options(const CqmomCaseVariant& case_)
{
    return std::visit(Overloaded{[](const CqmomDenseCase& dense_case) -> const ConditionalMomentOptions& { return dense_case.conditional_options; },
                                 [](const CqmomLeanCase& lean_case) -> const ConditionalMomentOptions& { return lean_case.conditional_options; }},
                      case_);
}

[[nodiscard]] bool adaptive_branch_counts(const CqmomCaseVariant& case_)
{
    return std::visit(Overloaded{[](const CqmomDenseCase& dense_case) { return dense_case.adaptive_branch_counts; },
                                 [](const CqmomLeanCase& lean_case) { return lean_case.adaptive_branch_counts; }},
                      case_);
}

[[nodiscard]] CqmomDiagnosticLevel diagnostic_level(const CqmomCaseVariant& case_)
{
    return std::visit(Overloaded{[](const CqmomDenseCase& dense_case) { return dense_case.diagnostic_level; },
                                 [](const CqmomLeanCase& lean_case) { return lean_case.diagnostic_level; }},
                      case_);
}

[[nodiscard]] VandermondeSolver vandermonde_solver(const CqmomCaseVariant& case_)
{
    return std::visit(Overloaded{[](const CqmomDenseCase& dense_case) { return dense_case.vandermonde_solver; },
                                 [](const CqmomLeanCase& lean_case) { return lean_case.vandermonde_solver; }},
                      case_);
}

[[nodiscard]] const CqmomCoordinateSupport& coordinate_support(const CqmomCaseVariant& case_)
{
    return std::visit(Overloaded{[](const CqmomDenseCase& dense_case) -> const CqmomCoordinateSupport& { return dense_case.coordinate_support; },
                                 [](const CqmomLeanCase& lean_case) -> const CqmomCoordinateSupport& { return lean_case.coordinate_support; }},
                      case_);
}

[[nodiscard]] CqmomMomentStorage moment_storage(const CqmomCaseVariant& case_)
{
    return std::holds_alternative<CqmomDenseCase>(case_) ? CqmomMomentStorage::dense : CqmomMomentStorage::lean;
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
    case VandermondeSolver::bjorck_pereyra:
        return "bjorck_pereyra";
    case VandermondeSolver::bjorck_pereyra_with_dense_fallback:
        return "bjorck_pereyra_with_dense_fallback";
    }
    return "unknown";
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

QmomState::QmomState(std::vector<double> moments, QmomResult quadrature) : moments_(std::move(moments)), quadrature_(std::move(quadrature)) {}

QmomState QmomCase::initial_state() const
{
    QmomResult quadrature = compute_qmom(*this);
    if (!quadrature.ok())
    {
        throw std::runtime_error("QMOM initial state reconstruction failed: " + quadrature.message);
    }
    return internal::QmomStateAccess::from_quadrature(std::move(quadrature), moments.size());
}

struct CqmomCase::Data
{
    CqmomCaseVariant case_;
    internal::PreparedFixedCqmomMetadataPtr prepared_fixed_metadata;
};

struct CqmomState::Data
{
    CqmomCaseVariant case_;
    internal::PreparedFixedCqmomMetadataPtr prepared_fixed_metadata;
    std::vector<CqmomActiveBranchCount> active_branch_counts;
    std::optional<CqmomQuadrature> accepted_quadrature;
};

std::string to_string(const CqmomMomentStorage storage)
{
    switch (storage)
    {
    case CqmomMomentStorage::lean:
        return "lean";
    case CqmomMomentStorage::dense:
        return "dense";
    }
    return "unknown";
}

std::string format_cqmom_settings(const CqmomMomentStorage storage, const bool adaptive_branch_counts, const VandermondeSolver vandermonde_solver,
                                  const CqmomDiagnosticLevel diagnostic_level)
{
    std::ostringstream output;
    output << "CQMOM settings:\n"
           << "  moment storage = " << to_string(storage) << '\n'
           << "  adaptive branch counts = " << (adaptive_branch_counts ? "enabled" : "disabled") << '\n'
           << "  Vandermonde solver = " << vandermonde_solver_name(vandermonde_solver) << '\n'
           << "  diagnostic level = " << diagnostic_level_name(diagnostic_level) << '\n';
    return output.str();
}

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
    return qmomkit::hierarchy(data_->case_);
}

const QmomOptions& CqmomCase::qmom_options() const noexcept
{
    return qmomkit::qmom_options(data_->case_);
}

const ConditionalMomentOptions& CqmomCase::conditional_options() const noexcept
{
    return qmomkit::conditional_options(data_->case_);
}

bool CqmomCase::adaptive_branch_counts() const noexcept
{
    return qmomkit::adaptive_branch_counts(data_->case_);
}

CqmomDiagnosticLevel CqmomCase::diagnostic_level() const noexcept
{
    return qmomkit::diagnostic_level(data_->case_);
}

VandermondeSolver CqmomCase::vandermonde_solver() const noexcept
{
    return qmomkit::vandermonde_solver(data_->case_);
}

const CqmomCoordinateSupport& CqmomCase::coordinate_support() const noexcept
{
    return qmomkit::coordinate_support(data_->case_);
}

CqmomMomentStorage CqmomCase::moment_storage() const noexcept
{
    return qmomkit::moment_storage(data_->case_);
}

CqmomState CqmomCase::initial_state() const
{
    const internal::PreparedFixedCqmomMetadataPtr& prepared_fixed_metadata = data_->prepared_fixed_metadata;
    const internal::PreparedFixedCqmomMetadata* prepared_metadata = prepared_fixed_metadata.get();
    CqmomResult accepted = std::visit(
        Overloaded{[prepared_metadata](const CqmomDenseCase& dense_case)
                   { return internal::compute_cqmom_internal(dense_case, {}, false, prepared_metadata, !dense_case.adaptive_branch_counts); },
                   [prepared_metadata](const CqmomLeanCase& lean_case)
                   { return internal::compute_cqmom_internal(lean_case, {}, false, prepared_metadata, !lean_case.adaptive_branch_counts); }},
        data_->case_);
    if (!accepted.ok())
    {
        throw std::runtime_error("CQMOM initial state quadrature reconstruction failed: " + accepted.message);
    }
    CqmomQuadrature accepted_quadrature = std::move(accepted.quadrature);
    std::vector<CqmomActiveBranchCount> active_counts = internal::active_branch_counts_from_quadrature(accepted_quadrature);
    return std::visit(Overloaded{[&active_counts, &accepted_quadrature, &prepared_fixed_metadata](const CqmomDenseCase& dense_case)
                                 {
                                     return internal::CqmomConstructionAccess::create_state(dense_case, std::move(active_counts),
                                                                                            std::move(accepted_quadrature), prepared_fixed_metadata);
                                 },
                                 [&active_counts, &accepted_quadrature, &prepared_fixed_metadata](const CqmomLeanCase& lean_case)
                                 {
                                     return internal::CqmomConstructionAccess::create_state(lean_case, std::move(active_counts),
                                                                                            std::move(accepted_quadrature), prepared_fixed_metadata);
                                 }},
                      data_->case_);
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
    return qmomkit::hierarchy(data_->case_);
}

const QmomOptions& CqmomState::qmom_options() const noexcept
{
    return qmomkit::qmom_options(data_->case_);
}

const ConditionalMomentOptions& CqmomState::conditional_options() const noexcept
{
    return qmomkit::conditional_options(data_->case_);
}

bool CqmomState::adaptive_branch_counts() const noexcept
{
    return qmomkit::adaptive_branch_counts(data_->case_);
}

CqmomDiagnosticLevel CqmomState::diagnostic_level() const noexcept
{
    return qmomkit::diagnostic_level(data_->case_);
}

VandermondeSolver CqmomState::vandermonde_solver() const noexcept
{
    return qmomkit::vandermonde_solver(data_->case_);
}

const CqmomCoordinateSupport& CqmomState::coordinate_support() const noexcept
{
    return qmomkit::coordinate_support(data_->case_);
}

CqmomMomentStorage CqmomState::moment_storage() const noexcept
{
    return qmomkit::moment_storage(data_->case_);
}

const std::vector<CqmomActiveBranchCount>& CqmomState::active_branch_counts() const noexcept
{
    return data_->active_branch_counts;
}

std::vector<MultiIndex> CqmomState::moment_indices() const
{
    return std::visit(Overloaded{[](const CqmomDenseCase& dense_case) { return enumerate_multi_indices(dense_case.moments.extents); },
                                 [](const CqmomLeanCase& lean_case) { return lean_case.moments.layout.required_indices; }},
                      data_->case_);
}

const std::vector<double>& CqmomState::moment_values() const noexcept
{
    return std::visit(Overloaded{[](const CqmomDenseCase& dense_case) -> const std::vector<double>& { return dense_case.moments.values; },
                                 [](const CqmomLeanCase& lean_case) -> const std::vector<double>& { return lean_case.moments.values; }},
                      data_->case_);
}

std::optional<double> CqmomState::moment(const MultiIndex& index) const
{
    return std::visit(Overloaded{[&index](const CqmomDenseCase& dense_case) -> std::optional<double>
                                 {
                                     if (index.size() != dense_case.moments.extents.size())
                                     {
                                         return std::nullopt;
                                     }
                                     for (std::size_t dimension = 0; dimension < index.size(); ++dimension)
                                     {
                                         if (index[dimension] >= dense_case.moments.extents[dimension])
                                         {
                                             return std::nullopt;
                                         }
                                     }
                                     return dense_case.moments.values[flatten_multi_index(index, dense_case.moments.extents)];
                                 },
                                 [&index](const CqmomLeanCase& lean_case) -> std::optional<double>
                                 {
                                     if (!lean_case.moments.contains_required_index(index))
                                     {
                                         return std::nullopt;
                                     }
                                     return lean_case.moments.moment(index);
                                 }},
                      data_->case_);
}

namespace internal
{
namespace
{
void validate_prepared_fixed_metadata(const PreparedFixedCqmomMetadata& metadata, const CqmomDenseCase& dense_case)
{
    // Fixed prepared metadata is created only from the canonical layout. For a
    // dense target, that canonical ordering is uniquely determined by the
    // hierarchy; dense extents are checked against the hierarchy at every
    // reconstruction and are not part of the cached graph.
    if (!prepared_fixed_cqmom_metadata_is_compatible(metadata, dense_case.hierarchy))
    {
        throw std::invalid_argument("Prepared fixed CQMOM metadata is incompatible with the target dense case.");
    }
}

void validate_prepared_fixed_metadata(const PreparedFixedCqmomMetadata& metadata, const CqmomLeanCase& lean_case)
{
    if (!prepared_fixed_cqmom_metadata_is_compatible(metadata, lean_case.moments.layout.hierarchy, &lean_case.moments.layout))
    {
        throw std::invalid_argument("Prepared fixed CQMOM metadata is incompatible with the target lean case.");
    }
}

[[nodiscard]] PreparedFixedCqmomMetadataPtr make_prepared_fixed_metadata(const CqmomDenseCase& dense_case,
                                                                         PreparedFixedCqmomMetadataPtr prepared_fixed_metadata)
{
    if (dense_case.adaptive_branch_counts)
    {
        if (prepared_fixed_metadata != nullptr)
        {
            throw std::invalid_argument("Adaptive CQMOM cases must not carry fixed prepared metadata.");
        }
        return nullptr;
    }
    if (prepared_fixed_metadata != nullptr)
    {
        validate_prepared_fixed_metadata(*prepared_fixed_metadata, dense_case);
        return prepared_fixed_metadata;
    }
    return create_prepared_fixed_cqmom_metadata(create_cqmom_moment_layout(dense_case.hierarchy));
}

[[nodiscard]] PreparedFixedCqmomMetadataPtr make_prepared_fixed_metadata(const CqmomLeanCase& lean_case,
                                                                         PreparedFixedCqmomMetadataPtr prepared_fixed_metadata)
{
    if (lean_case.adaptive_branch_counts)
    {
        if (prepared_fixed_metadata != nullptr)
        {
            throw std::invalid_argument("Adaptive CQMOM cases must not carry fixed prepared metadata.");
        }
        return nullptr;
    }
    if (prepared_fixed_metadata != nullptr)
    {
        validate_prepared_fixed_metadata(*prepared_fixed_metadata, lean_case);
        return prepared_fixed_metadata;
    }
    return create_prepared_fixed_cqmom_metadata(lean_case.moments.layout);
}
} // namespace

CqmomCase CqmomConstructionAccess::create_case(CqmomDenseCase dense_case, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata)
{
    prepared_fixed_metadata = make_prepared_fixed_metadata(dense_case, std::move(prepared_fixed_metadata));
    return CqmomCase(std::make_unique<CqmomCase::Data>(CqmomCase::Data{std::move(dense_case), std::move(prepared_fixed_metadata)}));
}

CqmomCase CqmomConstructionAccess::create_case(CqmomLeanCase lean_case, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata)
{
    prepared_fixed_metadata = make_prepared_fixed_metadata(lean_case, std::move(prepared_fixed_metadata));
    return CqmomCase(std::make_unique<CqmomCase::Data>(CqmomCase::Data{std::move(lean_case), std::move(prepared_fixed_metadata)}));
}

CqmomState CqmomConstructionAccess::create_state(CqmomDenseCase dense_case, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata)
{
    return create_state(std::move(dense_case), {}, std::nullopt, std::move(prepared_fixed_metadata));
}

CqmomState CqmomConstructionAccess::create_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 PreparedFixedCqmomMetadataPtr prepared_fixed_metadata)
{
    return create_state(std::move(dense_case), std::move(active_branch_counts), std::nullopt, std::move(prepared_fixed_metadata));
}

CqmomState CqmomConstructionAccess::create_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 std::optional<CqmomQuadrature> accepted_quadrature,
                                                 PreparedFixedCqmomMetadataPtr prepared_fixed_metadata)
{
    prepared_fixed_metadata = make_prepared_fixed_metadata(dense_case, std::move(prepared_fixed_metadata));
    return CqmomState(std::make_unique<CqmomState::Data>(CqmomState::Data{std::move(dense_case), std::move(prepared_fixed_metadata),
                                                                          std::move(active_branch_counts), std::move(accepted_quadrature)}));
}

CqmomState CqmomConstructionAccess::create_state(CqmomDenseCase dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 CqmomQuadrature&& accepted_quadrature, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata)
{
    return create_state(std::move(dense_case), std::move(active_branch_counts), std::optional<CqmomQuadrature>{std::move(accepted_quadrature)},
                        std::move(prepared_fixed_metadata));
}

CqmomState CqmomConstructionAccess::create_state(CqmomLeanCase lean_case, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata)
{
    return create_state(std::move(lean_case), {}, std::nullopt, std::move(prepared_fixed_metadata));
}

CqmomState CqmomConstructionAccess::create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 PreparedFixedCqmomMetadataPtr prepared_fixed_metadata)
{
    return create_state(std::move(lean_case), std::move(active_branch_counts), std::nullopt, std::move(prepared_fixed_metadata));
}

CqmomState CqmomConstructionAccess::create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 std::optional<CqmomQuadrature> accepted_quadrature,
                                                 PreparedFixedCqmomMetadataPtr prepared_fixed_metadata)
{
    prepared_fixed_metadata = make_prepared_fixed_metadata(lean_case, std::move(prepared_fixed_metadata));
    return CqmomState(std::make_unique<CqmomState::Data>(
        CqmomState::Data{std::move(lean_case), std::move(prepared_fixed_metadata), std::move(active_branch_counts), std::move(accepted_quadrature)}));
}

CqmomState CqmomConstructionAccess::create_state(CqmomLeanCase lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                                 CqmomQuadrature&& accepted_quadrature, PreparedFixedCqmomMetadataPtr prepared_fixed_metadata)
{
    return create_state(std::move(lean_case), std::move(active_branch_counts), std::optional<CqmomQuadrature>{std::move(accepted_quadrature)},
                        std::move(prepared_fixed_metadata));
}

const CqmomDenseCase* CqmomConstructionAccess::case_dense_case(const CqmomCase& case_)
{
    return std::get_if<CqmomDenseCase>(&case_.data_->case_);
}

const CqmomLeanCase* CqmomConstructionAccess::case_lean_case(const CqmomCase& case_)
{
    return std::get_if<CqmomLeanCase>(&case_.data_->case_);
}

const CqmomDenseCase* CqmomConstructionAccess::state_dense_case(const CqmomState& state)
{
    return std::get_if<CqmomDenseCase>(&state.data_->case_);
}

const CqmomLeanCase* CqmomConstructionAccess::state_lean_case(const CqmomState& state)
{
    return std::get_if<CqmomLeanCase>(&state.data_->case_);
}

const std::vector<CqmomActiveBranchCount>& CqmomConstructionAccess::state_active_branch_counts(const CqmomState& state)
{
    return state.data_->active_branch_counts;
}

const std::optional<CqmomQuadrature>& CqmomConstructionAccess::state_accepted_quadrature(const CqmomState& state)
{
    return state.data_->accepted_quadrature;
}

const PreparedFixedCqmomMetadataPtr& CqmomConstructionAccess::case_prepared_fixed_metadata(const CqmomCase& case_)
{
    return case_.data_->prepared_fixed_metadata;
}

const PreparedFixedCqmomMetadataPtr& CqmomConstructionAccess::state_prepared_fixed_metadata(const CqmomState& state)
{
    return state.data_->prepared_fixed_metadata;
}

CqmomQuadrature CqmomConstructionAccess::create_quadrature(PhaseSpaceHierarchy hierarchy, std::vector<double> weights,
                                                           std::vector<std::vector<double>> abscissas, const bool adaptive_branch_counts_enabled,
                                                           const bool rectangular_storage, std::vector<ConditioningPath> terminal_paths)
{
    return CqmomQuadrature(std::move(hierarchy), std::move(weights), std::move(abscissas), adaptive_branch_counts_enabled, rectangular_storage,
                           std::move(terminal_paths));
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
    CqmomDiagnostics converted;
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
    converted.dense_vandermonde_fallback_occurred = diagnostics.dense_vandermonde_fallback_occurred;
    converted.message = diagnostics.message;
    converted.bounded_support_checked = diagnostics.bounded_support_checked;
    converted.bounded_support_interval_present = diagnostics.bounded_support_interval_present;
    converted.bounded_support_interval = diagnostics.bounded_support_interval;
    converted.bounded_previous_active_node_count = diagnostics.bounded_previous_active_node_count;
    converted.bounded_candidate_node_count = diagnostics.bounded_candidate_node_count;
    converted.bounded_candidate_attempt_type = diagnostics.bounded_candidate_attempt_type;
    converted.bounded_candidate_moments = diagnostics.bounded_candidate_moments;
    converted.bounded_realizability_status = diagnostics.bounded_realizability_status;
    converted.bounded_effective_tolerance = diagnostics.bounded_effective_tolerance;
    converted.bounded_failed_matrix = diagnostics.bounded_failed_matrix;
    converted.bounded_minimum_margin = diagnostics.bounded_minimum_margin;
    converted.bounded_accepted_fallback_node_count = diagnostics.bounded_accepted_fallback_node_count;
    converted.bounded_returned_nodes = diagnostics.bounded_returned_nodes;
    converted.bounded_returned_weights = diagnostics.bounded_returned_weights;
    converted.bounded_support_determinants = diagnostics.bounded_support_determinants;
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
    return CqmomResult{internal::CqmomConstructionAccess::create_quadrature(
                           result.quadrature.hierarchy(), result.quadrature.weights(), result.quadrature.abscissas(),
                           result.quadrature.adaptive_branch_counts_enabled(), result.quadrature.has_rectangular_storage(),
                           terminal_paths_for_public_quadrature(result.quadrature)),
                       to_cqmom_status(result.status), result.message, to_cqmom_diagnostics(result.diagnostics),
                       failed_phase_space_dimension(result, hierarchy)};
}
} // namespace

CqmomResult compute_cqmom(const CqmomDenseCase& dense_case)
{
    return internal::compute_cqmom_internal(dense_case);
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
                                                                 {},
                                                                 false,
                                                                 lean_case.coordinate_support});

    return to_cqmom_result(result, hierarchy);
}

namespace internal
{

CqmomResult compute_cqmom_internal(const CqmomDenseCase& dense_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                   const bool allow_active_branch_count_growth, const PreparedFixedCqmomMetadata* const prepared_fixed_metadata,
                                   const bool require_prepared_fixed_metadata)
{
    if (require_prepared_fixed_metadata && !dense_case.adaptive_branch_counts && prepared_fixed_metadata == nullptr)
    {
        throw std::logic_error("Fixed normal CQMOM reconstruction requires prepared metadata.");
    }
    RecursiveQuadratureResult result{RecursiveQuadrature(dense_case.hierarchy), RecursiveQuadratureStatus::success, "", {}, {}, 0};
    result.quadrature.set_adaptive_branch_counts_enabled(dense_case.adaptive_branch_counts);
    std::optional<CqmomMomentData> lean_moments;

    try
    {
        validate_dense_cqmom_reconstruction_tensor(dense_case.moments, dense_case.hierarchy, "Raw moment extents must match the hierarchy dimension.",
                                                   "Raw moment data size does not match its extents.", "Raw moments must be finite.",
                                                   "Raw moment extents must contain at least 2*N orders per dimension.");
        std::optional<CqmomMomentLayout> generated_layout;
        const CqmomMomentLayout* layout = prepared_fixed_metadata == nullptr ? nullptr : &prepared_fixed_metadata->moment_layout;
        if (layout == nullptr)
        {
            generated_layout.emplace(create_cqmom_moment_layout(dense_case.hierarchy));
            layout = &*generated_layout;
        }
        lean_moments = extract_cqmom_moments(dense_case.moments, *layout);
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
             dense_case.adaptive_branch_counts, dense_case.diagnostic_level, dense_case.vandermonde_solver, std::move(active_branch_counts),
             allow_active_branch_count_growth, dense_case.coordinate_support, prepared_fixed_metadata, require_prepared_fixed_metadata});
    }

    return to_cqmom_result(result, dense_case.hierarchy);
}

CqmomResult compute_cqmom_internal(const CqmomLeanCase& lean_case, std::vector<CqmomActiveBranchCount> active_branch_counts,
                                   const bool allow_active_branch_count_growth, const PreparedFixedCqmomMetadata* const prepared_fixed_metadata,
                                   const bool require_prepared_fixed_metadata)
{
    if (require_prepared_fixed_metadata && !lean_case.adaptive_branch_counts && prepared_fixed_metadata == nullptr)
    {
        throw std::logic_error("Fixed normal CQMOM reconstruction requires prepared metadata.");
    }
    const PhaseSpaceHierarchy& hierarchy = lean_case.moments.layout.hierarchy;
    const RecursiveQuadratureResult result = compute_recursive_quadrature_from_trusted_moment_data(
        {hierarchy, lean_case.moments, lean_case.conditional_options, to_wheeler_options(lean_case.qmom_options), lean_case.adaptive_branch_counts,
         lean_case.diagnostic_level, lean_case.vandermonde_solver, std::move(active_branch_counts), allow_active_branch_count_growth,
         lean_case.coordinate_support, prepared_fixed_metadata, require_prepared_fixed_metadata});

    return to_cqmom_result(result, hierarchy);
}

} // namespace internal

CqmomResult compute_cqmom(const CqmomCase& case_)
{
    const internal::PreparedFixedCqmomMetadata* prepared_fixed_metadata = case_.data_->prepared_fixed_metadata.get();
    return std::visit(
        Overloaded{[prepared_fixed_metadata](const CqmomDenseCase& dense_case)
                   { return internal::compute_cqmom_internal(dense_case, {}, false, prepared_fixed_metadata, !dense_case.adaptive_branch_counts); },
                   [prepared_fixed_metadata](const CqmomLeanCase& lean_case)
                   { return internal::compute_cqmom_internal(lean_case, {}, false, prepared_fixed_metadata, !lean_case.adaptive_branch_counts); }},
        case_.data_->case_);
}

} // namespace qmomkit
