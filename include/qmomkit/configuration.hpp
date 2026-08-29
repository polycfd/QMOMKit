#pragma once

#include "qmomkit/cqmom.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/output.hpp"
#include "qmomkit/phase_space_evolution.hpp"
#include "qmomkit/qmom.hpp"
#include "qmomkit/realizability_policy.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace qmomkit::config
{

enum class ConfigurationErrorCode
{
    file_error,
    syntax_error,
    duplicate_entry,
    unknown_entry,
    missing_entry,
    invalid_value,
    invalid_configuration
};

struct ConfigurationError
{
    ConfigurationErrorCode code = ConfigurationErrorCode::syntax_error;
    std::string source_name;
    std::size_t line = 0;
    std::string message;
};

template <class T> class ConfigurationResult
{
public:
    ConfigurationResult() = default;

    ConfigurationResult(T value) : value_(std::move(value)) {}

    [[nodiscard]] static ConfigurationResult failure(std::string error)
    {
        ConfigurationResult result;
        result.error_ = std::move(error);
        return result;
    }

    [[nodiscard]] bool ok() const noexcept
    {
        return value_.has_value() && error_.empty();
    }

    [[nodiscard]] explicit operator bool() const noexcept
    {
        return ok();
    }

    [[nodiscard]] T& value() &
    {
        return value_.value();
    }

    [[nodiscard]] const T& value() const&
    {
        return value_.value();
    }

    [[nodiscard]] T value() &&
    {
        return std::move(value_).value();
    }

    [[nodiscard]] const std::string& error() const noexcept
    {
        return error_;
    }

private:
    std::optional<T> value_;
    std::string error_;
};

struct ExplicitMomentInitialization
{
    MomentSet moments;
};

struct ExplicitTensorInitialization
{
    RawMomentTensor moments;
};

struct DistributionInitialization
{
    DistributionSpecification distribution;
    std::size_t maximum_order = 0;
    double number_density = 1.0;
};

struct MixtureInitialization
{
    DistributionMixture mixture;
    std::size_t maximum_order = 0;
};

struct IndependentProductInitialization
{
    IndependentProductDistribution distribution;
    std::optional<std::vector<std::size_t>> maximum_orders;
};

using InitialConditionSpecification = std::variant<ExplicitMomentInitialization, ExplicitTensorInitialization, DistributionInitialization,
                                                   MixtureInitialization, IndependentProductInitialization>;

/**
 * Generate one-dimensional raw moments from an initialization specification.
 *
 * Explicit one-dimensional moments are returned unchanged. Analytical
 * distribution and mixture specifications are evaluated through the
 * raw-moment generation API.
 */
[[nodiscard]] ConfigurationResult<MomentSet> initialize_moments(const InitialConditionSpecification& specification,
                                                                std::string source_name = "<initialization>");

/**
 * Generate dense multidimensional raw moments from an initialization specification.
 *
 * Explicit tensors are returned unchanged. Independent product specifications
 * are evaluated through the tensor raw-moment generation API and require dense
 * tensor maximum-order metadata.
 */
[[nodiscard]] ConfigurationResult<RawMomentTensor> initialize_cqmom_dense_moments(const InitialConditionSpecification& specification,
                                                                                  std::string source_name = "<initialization>");

/**
 * Generate the multidimensional raw moments required by a CQMOM hierarchy.
 *
 * Explicit tensor specifications are adapted by extracting the required CQMOM
 * moments. Independent product specifications are evaluated through selected
 * raw-moment generation, avoiding construction of a dense `RawMomentTensor`.
 */
[[nodiscard]] ConfigurationResult<CqmomMomentData> initialize_cqmom_moment_data(const InitialConditionSpecification& specification,
                                                                                const PhaseSpaceHierarchy& hierarchy,
                                                                                std::string source_name = "<initialization>");

/** Parse current QMOM initialization into explicit or generated moments. */
[[nodiscard]] ConfigurationResult<InitialConditionSpecification> parse_qmom_initial_condition(std::string_view text,
                                                                                              std::string source_name = "<string>");

/** Parse current CQMOM initialization into explicit or generated moments. */
[[nodiscard]] ConfigurationResult<InitialConditionSpecification> parse_cqmom_initial_condition(std::string_view text,
                                                                                               std::string source_name = "<string>");

/** Create a QMOM case from an already parsed initialization specification. */
[[nodiscard]] ConfigurationResult<QmomCase> create_qmom_case(const InitialConditionSpecification& initial_condition, std::size_t node_count,
                                                             const QmomOptions& options = {}, std::string source_name = "<initialization>",
                                                             std::optional<AbscissaBounds> abscissa_bounds = {});

/**
 * Advanced: create an explicit dense CQMOM case from parsed initialization.
 *
 * Coordinate support is indexed by phase-space dimension. Omit it, or leave it
 * empty, for unconstrained real-line reconstruction.
 */
[[nodiscard]] ConfigurationResult<CqmomDenseCase> create_cqmom_dense_case(const InitialConditionSpecification& initial_condition,
                                                                          const PhaseSpaceHierarchy& hierarchy, const QmomOptions& qmom_options = {},
                                                                          bool adaptive_branch_counts = false,
                                                                          std::string source_name = "<initialization>",
                                                                          CqmomCoordinateSupport coordinate_support = {});

/**
 * Advanced: create an explicit lean CQMOM case from parsed initialization.
 *
 * Coordinate support has the same phase-space-indexed semantics as the dense
 * factory and is carried into the resulting case/state.
 */
[[nodiscard]] ConfigurationResult<CqmomLeanCase> create_cqmom_lean_case(const InitialConditionSpecification& initial_condition,
                                                                        const PhaseSpaceHierarchy& hierarchy, const QmomOptions& qmom_options = {},
                                                                        bool adaptive_branch_counts = false,
                                                                        std::string source_name = "<initialization>",
                                                                        CqmomCoordinateSupport coordinate_support = {});

/**
 * Create a normal CQMOM case with replacement explicit initial raw moments.
 *
 * The returned immutable case preserves the base case hierarchy, storage
 * choice, coordinate support, and all reconstruction policies. Dense cases
 * retain the complete supplied tensor; lean cases extract their fixed required
 * CQMOM moment layout from it. This operation validates tensor structure only;
 * CQMOM realizability and reconstruction remain the responsibility of
 * `compute_cqmom(...)` and `CqmomCase::initial_state()`.
 */
[[nodiscard]] ConfigurationResult<CqmomCase> with_initial_moments(const CqmomCase& base_case, RawMomentTensor explicit_moments);

/** Parse a current QMOM configuration string. */
[[nodiscard]] ConfigurationResult<QmomCase> parse_qmom(std::string_view text, std::string source_name = "<string>");

/** Advanced: parse configuration into an explicit dense CQMOM case. */
[[nodiscard]] ConfigurationResult<CqmomDenseCase> parse_cqmom_dense(std::string_view text, std::string source_name = "<string>");

/** Advanced: parse configuration into an explicit lean CQMOM case. */
[[nodiscard]] ConfigurationResult<CqmomLeanCase> parse_cqmom_lean(std::string_view text, std::string source_name = "<string>");

/** Parse a current CQMOM configuration string into a user-facing CQMOM case. */
[[nodiscard]] ConfigurationResult<CqmomCase> parse_cqmom(std::string_view text, std::string source_name = "<string>");

/** Parse optional runtime realizability policy configuration. */
[[nodiscard]] ConfigurationResult<RealizabilityPolicy> parse_realizability_policy(std::string_view text, std::string source_name = "<string>");

/** Parse sparse, phase-space-coordinate-indexed runtime realizability policies. */
[[nodiscard]] ConfigurationResult<std::vector<std::optional<RealizabilityPolicy>>>
parse_coordinate_realizability_policies(std::string_view text, std::string source_name = "<string>");

/** Parse an optional phase-space evolution scheme configuration. */
[[nodiscard]] ConfigurationResult<PhaseSpaceEvolutionScheme> parse_phase_space_evolution_scheme(std::string_view text,
                                                                                                std::string source_name = "<string>");

/** Parse optional node-advection configuration, including embedded-RK controls. */
[[nodiscard]] ConfigurationResult<NodeAdvectionOptions> parse_node_advection_options(std::string_view text, std::string source_name = "<string>");

/** Parse QMOM phase-space evolution options from optional runtime configuration. */
[[nodiscard]] ConfigurationResult<QmomPhaseSpaceEvolutionOptions> parse_qmom_phase_space_evolution_options(std::string_view text,
                                                                                                           std::string source_name = "<string>");

/** Parse CQMOM phase-space evolution options from optional runtime configuration. */
[[nodiscard]] ConfigurationResult<CqmomPhaseSpaceEvolutionOptions> parse_cqmom_phase_space_evolution_options(std::string_view text,
                                                                                                             std::string source_name = "<string>");

/** Parse optional output configuration. */
[[nodiscard]] ConfigurationResult<OutputOptions> parse_output_options(std::string_view text, std::string source_name = "<string>");

/** Load and parse a QMOM configuration from exactly the requested file. */
[[nodiscard]] ConfigurationResult<QmomCase> load_qmom(const std::filesystem::path& filename);

/** Advanced: load an explicit dense CQMOM case from exactly the requested file. */
[[nodiscard]] ConfigurationResult<CqmomDenseCase> load_cqmom_dense(const std::filesystem::path& filename);

/** Advanced: load an explicit lean CQMOM case from exactly the requested file. */
[[nodiscard]] ConfigurationResult<CqmomLeanCase> load_cqmom_lean(const std::filesystem::path& filename);

/** Load and parse a CQMOM configuration into a user-facing CQMOM case from exactly the requested file. */
[[nodiscard]] ConfigurationResult<CqmomCase> load_cqmom(const std::filesystem::path& filename);

/** Load and parse output configuration from exactly the requested file. */
[[nodiscard]] ConfigurationResult<OutputOptions> load_output_options(const std::filesystem::path& filename);

/** Load and parse node-advection configuration from exactly the requested file. */
[[nodiscard]] ConfigurationResult<NodeAdvectionOptions> load_node_advection_options(const std::filesystem::path& filename);

/** Load and parse runtime realizability policy configuration from exactly the requested file. */
[[nodiscard]] ConfigurationResult<RealizabilityPolicy> load_realizability_policy(const std::filesystem::path& filename);

/** Load and parse QMOM phase-space evolution options from exactly the requested file. */
[[nodiscard]] ConfigurationResult<QmomPhaseSpaceEvolutionOptions> load_qmom_phase_space_evolution_options(const std::filesystem::path& filename);

/** Load and parse CQMOM phase-space evolution options from exactly the requested file. */
[[nodiscard]] ConfigurationResult<CqmomPhaseSpaceEvolutionOptions> load_cqmom_phase_space_evolution_options(const std::filesystem::path& filename);

} // namespace qmomkit::config
