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
 * Create a CQMOM case definition from an already parsed initialization
 * specification.
 *
 * Coordinate support is indexed by phase-space dimension and is carried into
 * the resulting case and any state created from it. Construct a
 * `CqmomPreparedCase` from the result for the normal CQMOM workflow.
 */
[[nodiscard]] ConfigurationResult<CqmomCase> create_cqmom_case(const InitialConditionSpecification& initial_condition,
                                                               const PhaseSpaceHierarchy& hierarchy, const QmomOptions& qmom_options = {},
                                                               bool adaptive_branch_counts = false, std::string source_name = "<initialization>",
                                                               CqmomCoordinateSupport coordinate_support = {});

/**
 * Create a prepared CQMOM case with replacement initial moments, taken from a
 * full raw moment tensor.
 *
 * This is the entry point for moment data produced outside QMOMKit, such as the
 * raw moments the Monte Carlo companion evaluates from a sample population: the
 * caller supplies the whole tensor and QMOMKit selects the moments the
 * hierarchy reconstructs from. Prefer this overload over extracting those
 * moments yourself when the tensor comes from outside, because a tensor that is
 * the wrong rank, the wrong size, or too small to contain a required moment is
 * reported here as a configuration error you can inspect, whereas
 * `extract_cqmom_moments(...)` signals the same problems by throwing.
 *
 * The returned immutable case keeps the hierarchy, coordinate support,
 * reconstruction policies, and prepared metadata of the base case, which is
 * left unchanged. Cases derived this way therefore share one recursive branch
 * structure. This operation validates tensor structure only; CQMOM
 * realizability and reconstruction remain the responsibility of
 * `compute_cqmom(...)` and `CqmomPreparedCase::initial_state()`.
 */
[[nodiscard]] ConfigurationResult<CqmomPreparedCase> with_initial_moments(const CqmomPreparedCase& base_cqmom_prepared_case,
                                                                          RawMomentTensor explicit_moments);

/**
 * Create a prepared CQMOM case with replacement initial moments, supplied as
 * the required CQMOM moment set.
 *
 * Use this overload when the caller already holds exactly the moments the
 * hierarchy reconstructs from, so that no full tensor has to be built for
 * entries the case never reads.
 *
 * The returned immutable case keeps the hierarchy, coordinate support,
 * reconstruction policies, and prepared metadata of the base case, which is
 * left unchanged. The supplied layout must match the required CQMOM moment
 * layout of the base hierarchy. This operation validates moment-set structure
 * only; CQMOM realizability and reconstruction remain the responsibility of
 * `compute_cqmom(...)` and `CqmomPreparedCase::initial_state()`.
 */
[[nodiscard]] ConfigurationResult<CqmomPreparedCase> with_initial_moments(const CqmomPreparedCase& base_cqmom_prepared_case,
                                                                          CqmomMomentData explicit_moments);

/**
 * Create a CQMOM case definition with replacement initial moments, taken from a
 * full raw moment tensor.
 *
 * The tensor is validated and its moments selected exactly as by the
 * `CqmomPreparedCase` overload, with the same errors. The returned case keeps
 * the settings of the base case, which is left unchanged, but holds no prepared
 * metadata.
 */
[[nodiscard]] ConfigurationResult<CqmomCase> with_initial_moments(const CqmomCase& base_cqmom_case, RawMomentTensor explicit_moments);

/**
 * Create a CQMOM case definition with replacement initial moments, supplied as
 * the required CQMOM moment set.
 *
 * The moment set is validated exactly as by the `CqmomPreparedCase` overload,
 * with the same errors. The returned case keeps the settings of the base case,
 * which is left unchanged, but holds no prepared metadata.
 */
[[nodiscard]] ConfigurationResult<CqmomCase> with_initial_moments(const CqmomCase& base_cqmom_case, CqmomMomentData explicit_moments);

/** Parse a current QMOM configuration string. */
[[nodiscard]] ConfigurationResult<QmomCase> parse_qmom(std::string_view text, std::string source_name = "<string>");

/** Parse a CQMOM configuration string into a CQMOM case definition. */
[[nodiscard]] ConfigurationResult<CqmomCase> parse_cqmom(std::string_view text, std::string source_name = "<string>");

/** Parse a CQMOM configuration string into a prepared CQMOM case, the normal CQMOM workflow. */
[[nodiscard]] ConfigurationResult<CqmomPreparedCase> parse_cqmom_prepared(std::string_view text, std::string source_name = "<string>");

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

/** Load and parse a CQMOM configuration into a CQMOM case definition from exactly the requested file. */
[[nodiscard]] ConfigurationResult<CqmomCase> load_cqmom(const std::filesystem::path& filename);

/** Load and parse a CQMOM configuration into a prepared CQMOM case, the normal CQMOM workflow, from exactly the requested file. */
[[nodiscard]] ConfigurationResult<CqmomPreparedCase> load_cqmom_prepared(const std::filesystem::path& filename);

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
