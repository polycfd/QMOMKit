#include "qmomkit/configuration.hpp"

#include "internal/cqmom_case.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/moment_reconstruction.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace
{

const char* qmom_text = R"(
    # leading comment
    [phase_space]
    names = size
    node_counts = 2

    [qmom]
    adaptive = false

    [moments]
    values = 1.0, 0.5, 0.5, 0.5
)";

const char* cqmom_text = R"(
[phase_space]
names = size, velocity
node_counts = 1, 1
conditioning_order = 1, 0
adaptive_branch_counts = true

[qmom]
adaptive = true

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 1.5
M[1,0] = 0.5
M[1,1] = 0.75
)";

std::string qmom_case_text_with_runtime(const std::string& runtime)
{
    return R"(
[phase_space]
node_counts = 1

[moments]
values = 1.0, 0.0
)" + runtime;
}

std::string cqmom_case_text_with_runtime(const std::string& runtime)
{
    return R"(
[phase_space]
node_counts = 1, 1

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.0
M[1,0] = 0.0
M[1,1] = 0.0
)" + runtime;
}

template <class T> bool has_failure(const qmomkit::config::ConfigurationResult<T>& result)
{
    return !result && !result.error().empty();
}

template <class T> bool has_message_fragment(const qmomkit::config::ConfigurationResult<T>& result, const std::string& fragment)
{
    return result.error().find(fragment) != std::string::npos;
}

template <class T, class = void> struct has_public_value_member : std::false_type
{
};

template <class T> struct has_public_value_member<T, std::void_t<decltype(std::declval<T&>().value)>> : std::true_type
{
};

void assert_near(const double actual, const double expected)
{
    assert(std::abs(actual - expected) <= 1.0e-13);
}

void assert_moments_near(const qmomkit::MomentSet& actual, const qmomkit::MomentSet& expected)
{
    assert(actual.size() == expected.size());
    for (std::size_t order = 0; order < actual.size(); ++order)
    {
        assert_near(actual[order], expected[order]);
    }
}

void assert_tensor_near(const qmomkit::RawMomentTensor& actual, const qmomkit::RawMomentTensor& expected)
{
    assert(actual.extents == expected.extents);
    assert(actual.values.size() == expected.values.size());
    for (std::size_t flat_index = 0; flat_index < actual.values.size(); ++flat_index)
    {
        assert_near(actual.values[flat_index], expected.values[flat_index]);
    }
}

std::vector<double> raw_moments_from_terminals_2d(const std::vector<double>& weights, const std::vector<std::vector<double>>& abscissas)
{
    std::vector<double> moments(16, 0.0);
    for (std::size_t ix = 0; ix < 4; ++ix)
    {
        for (std::size_t iy = 0; iy < 4; ++iy)
        {
            double moment = 0.0;
            for (std::size_t node = 0; node < weights.size(); ++node)
            {
                moment += weights[node] * std::pow(abscissas[0][node], static_cast<int>(ix)) * std::pow(abscissas[1][node], static_cast<int>(iy));
            }
            moments[4 * ix + iy] = moment;
        }
    }
    return moments;
}

std::string cqmom_adaptive_branch_count_text_with_fixed_qmom(const std::vector<double>& weights, const std::vector<std::vector<double>>& abscissas,
                                                             const std::vector<std::size_t>& conditioning_order)
{
    const std::vector<double> moments = raw_moments_from_terminals_2d(weights, abscissas);
    std::ostringstream text;
    text.precision(17);
    text << "[phase_space]\n"
         << "node_counts = 2, 2\n"
         << "conditioning_order = " << conditioning_order[0] << ", " << conditioning_order[1] << "\n"
         << "adaptive_branch_counts = true\n\n"
         << "[qmom]\n"
         << "adaptive = false\n\n"
         << "[cqmom]\n"
         << "diagnostic_level = basic\n\n"
         << "[moments]\n"
         << "extents = 4, 4\n";
    for (std::size_t ix = 0; ix < 4; ++ix)
    {
        for (std::size_t iy = 0; iy < 4; ++iy)
        {
            text << "M[" << ix << "," << iy << "] = " << moments[4 * ix + iy] << '\n';
        }
    }
    return text.str();
}

void assert_cqmom_moments_near(const qmomkit::CqmomMomentData& actual, const qmomkit::CqmomMomentData& expected)
{
    assert(actual.layout.hierarchy.node_counts() == expected.layout.hierarchy.node_counts());
    assert(actual.layout.hierarchy.conditioning_order() == expected.layout.hierarchy.conditioning_order());
    assert(actual.layout.required_indices == expected.layout.required_indices);
    assert(actual.values.size() == expected.values.size());
    for (std::size_t index = 0; index < actual.values.size(); ++index)
    {
        assert_near(actual.values[index], expected.values[index]);
    }
}

void assert_cqmom_quadratures_near(const qmomkit::CqmomQuadrature& actual, const qmomkit::CqmomQuadrature& expected)
{
    assert(actual.hierarchy().node_counts() == expected.hierarchy().node_counts());
    assert(actual.hierarchy().conditioning_order() == expected.hierarchy().conditioning_order());
    assert(actual.weights().size() == expected.weights().size());
    assert(actual.abscissas().size() == expected.abscissas().size());
    for (std::size_t node = 0; node < actual.weights().size(); ++node)
    {
        assert_near(actual.weights()[node], expected.weights()[node]);
    }
    for (std::size_t dimension = 0; dimension < actual.abscissas().size(); ++dimension)
    {
        assert(actual.abscissas()[dimension].size() == expected.abscissas()[dimension].size());
        for (std::size_t node = 0; node < actual.abscissas()[dimension].size(); ++node)
        {
            assert_near(actual.abscissas()[dimension][node], expected.abscissas()[dimension][node]);
        }
    }
}

void assert_coordinate_support_matches_shifted_interval(const qmomkit::CqmomCoordinateSupport& support)
{
    assert(support.size() == 2);
    assert(!support[0].has_value());
    assert(support[1].has_value());
    assert_near(support[1]->lower, 2.0);
    assert_near(support[1]->upper, 5.0);
}

void assert_realizability_policies_match(const qmomkit::RealizabilityPolicy& actual, const qmomkit::RealizabilityPolicy& expected)
{
    assert(actual.check_mode == expected.check_mode);
    assert(actual.support == expected.support);
    assert(actual.accept_boundary == expected.accept_boundary);
    assert_near(actual.options.tolerance, expected.options.tolerance);
    assert(actual.options.relative_tolerance == expected.options.relative_tolerance);
    assert(actual.interval.has_value() == expected.interval.has_value());
    if (actual.interval.has_value())
    {
        assert_near(actual.interval->lower, expected.interval->lower);
        assert_near(actual.interval->upper, expected.interval->upper);
    }
}

struct MoveOnlyValue
{
    explicit MoveOnlyValue(const int input) : value(std::make_unique<int>(input)) {}
    MoveOnlyValue(MoveOnlyValue&&) noexcept = default;
    MoveOnlyValue& operator=(MoveOnlyValue&&) noexcept = default;
    MoveOnlyValue(const MoveOnlyValue&) = delete;
    MoveOnlyValue& operator=(const MoveOnlyValue&) = delete;

    std::unique_ptr<int> value;
};

void test_configuration_result_api()
{
    qmomkit::config::ConfigurationResult<std::string> success{"moments"};
    static_assert(!has_public_value_member<decltype(success)>::value, "ConfigurationResult must not expose a public value field.");
    assert(success);
    assert(success.ok());
    assert(success.error().empty());
    success.value() = "quadrature";
    const auto& const_success = success;
    assert(const_success.value() == "quadrature");
    assert(std::move(success).value() == "quadrature");

    qmomkit::config::ConfigurationResult<MoveOnlyValue> move_only{MoveOnlyValue{7}};
    MoveOnlyValue moved = std::move(move_only).value();
    assert(*moved.value == 7);

    const auto failed = qmomkit::config::ConfigurationResult<int>::failure("configuration failed");
    assert(!failed);
    assert(!failed.ok());
    assert(failed.error() == "configuration failed");
}

void test_valid_qmom_matches_programmatic_case()
{
    const auto configured = qmomkit::config::parse_qmom(qmom_text, "qmom-test");
    assert(configured.ok());
    const qmomkit::QmomCase direct{qmomkit::MomentSet({1.0, 0.5, 0.5, 0.5}), 2, {}};
    assert(configured.value().node_count == direct.node_count);
    assert(configured.value().moments.values() == direct.moments.values());
    assert(!configured.value().options.adaptive);
    assert(!configured.value().options.automatic_scaling);

    const auto configured_result = qmomkit::compute_qmom(configured.value());
    const auto direct_result = qmomkit::compute_qmom(direct);
    assert(configured_result.ok());
    assert(direct_result.ok());
    assert(configured_result.weights == direct_result.weights);
    assert(configured_result.abscissas == direct_result.abscissas);

    const auto adaptive = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[qmom]\nadaptive=true\n[moments]\nvalues=1,0\n");
    assert(adaptive.ok());
    assert(adaptive.value().options.adaptive);
}

void test_qmom_case_parses_abscissa_bounds()
{
    const auto unbounded = qmomkit::config::parse_qmom(qmom_text, "qmom-unbounded");
    assert(unbounded.ok());
    assert(!unbounded.value().abscissa_bounds.has_value());

    const auto bounded = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 2
support_lower = 0.0
support_upper = 1.0

[moments]
values = 1.0, 0.5, 0.3125, 0.21875
)",
                                                     "qmom-bounded");
    assert(bounded.ok());
    assert(bounded.value().abscissa_bounds.has_value());
    assert_near(bounded.value().abscissa_bounds->lower, 0.0);
    assert_near(bounded.value().abscissa_bounds->upper, 1.0);

    const auto explicitly_unbounded = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1
support_lower = none
support_upper = none

[moments]
values = 1.0, 0.5
)");
    assert(explicitly_unbounded.ok());
    assert(!explicitly_unbounded.value().abscissa_bounds.has_value());

    const auto rejected = [](const char* support_entries)
    {
        const auto result = qmomkit::config::parse_qmom(std::string(R"(
[phase_space]
node_counts = 1
)") + support_entries + R"(
[moments]
values = 1.0, 0.5
)");
        assert(!result.ok());
    };
    rejected("support_lower = 0.0\n");
    rejected("support_lower = none\nsupport_upper = 1.0\n");
    rejected("support_lower = 0.0, 1.0\nsupport_upper = 1.0, 2.0\n");
    rejected("support_lower = 1.0\nsupport_upper = 1.0\n");
}

void test_qmom_distribution_initialization_matches_direct_api()
{
    const std::string omitted_order_text = R"(
[phase_space]
names = size
node_counts = 3

[initial]
type = distribution
family = lognormal
number_density = 2.5
mu_log = -0.2
sigma_log = 0.4
)";
    const std::string explicit_order_text = R"(
[phase_space]
names = size
node_counts = 3

[initial]
type = distribution
family = lognormal
number_density = 2.5
maximum_order = 5
mu_log = -0.2
sigma_log = 0.4
)";
    const auto omitted_order = qmomkit::config::parse_qmom(omitted_order_text);
    const auto explicit_order = qmomkit::config::parse_qmom(explicit_order_text);
    assert(omitted_order.ok());
    assert(explicit_order.ok());

    const qmomkit::RawMomentResult direct = qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{-0.2, 0.4}, 5, 2.5);
    assert(direct.ok());
    assert(omitted_order.value().node_count == 3);
    assert(omitted_order.value().node_count == explicit_order.value().node_count);
    assert(omitted_order.value().moments.size() == 6);
    assert(omitted_order.value().moments.size() == explicit_order.value().moments.size());
    assert_near(omitted_order.value().moments[0], 2.5);
    assert_moments_near(omitted_order.value().moments, explicit_order.value().moments);
    assert_moments_near(omitted_order.value().moments, direct.moments);

    const auto initial = qmomkit::config::parse_qmom_initial_condition(R"(
[phase_space]
node_counts = 2

[initial]
type = distribution
family = gamma
number_density = 4.0
maximum_order = 3
shape = 2.0
scale = 3.0
)");
    assert(initial.ok());
    const auto moments = qmomkit::config::initialize_moments(initial.value());
    assert(moments.ok());
    const qmomkit::RawMomentResult expected = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{2.0, 3.0}, 3, 4.0);
    assert(expected.ok());
    assert_moments_near(moments.value(), expected.moments);

    const auto qmom_case = qmomkit::config::create_qmom_case(initial.value(), 2);
    assert(qmom_case.ok());
    assert(qmom_case.value().moments.size() == 4);
    assert_near(qmom_case.value().moments[0], 4.0);
}

void test_all_supported_distribution_families_parse()
{
    const std::vector<std::string> texts{"family=normal\nmaximum_order=1\nmean=1\nstandard_deviation=2\n",
                                         "family=beta\nmaximum_order=1\nalpha=2\nbeta=5\n", "family=chi\nmaximum_order=1\ndegrees_of_freedom=3\n",
                                         "family=chi-squared\nmaximum_order=1\ndegrees_of_freedom=4\n",
                                         "family=weibull\nmaximum_order=1\nshape=2\nscale=3\n"};

    for (const std::string& fields : texts)
    {
        const auto configured =
            qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nnumber_density=2\n" + fields);
        assert(configured.ok());
        assert(configured.value().moments.size() == 2);
        assert_near(configured.value().moments[0], 2.0);
    }
}

void test_qmom_mixture_initialization_matches_direct_api()
{
    const auto configured = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 4

[initial]
type = mixture
maximum_order = 7

[initial.mixture.component.0]
family = gamma
number_density = 1.5
shape = 2.0
scale = 0.4

[initial.mixture.component.1]
family = lognormal
number_density = 0.8
mu_log = -0.2
sigma_log = 0.5
)");
    assert(configured.ok());

    const qmomkit::DistributionMixture direct{{{qmomkit::GammaDistribution{2.0, 0.4}, 1.5}, {qmomkit::LogNormalDistribution{-0.2, 0.5}, 0.8}}};
    const qmomkit::RawMomentResult expected = qmomkit::compute_raw_moments(direct, 7);
    assert(expected.ok());
    assert(configured.value().node_count == 4);
    assert(configured.value().moments.size() == 8);
    assert_near(configured.value().moments[0], 2.3);
    assert_moments_near(configured.value().moments, expected.moments);
}

void test_mixture_zero_density_component_is_skipped()
{
    const auto configured = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 2

[initial]
type = mixture
maximum_order = 3

[initial.mixture.component.0]
family = normal
number_density = 0.0
mean = 0.0
standard_deviation = 0.0

[initial.mixture.component.1]
family = gamma
number_density = 3.0
shape = 1.0
scale = 2.0
)");
    assert(configured.ok());
    assert_near(configured.value().moments[0], 3.0);
    assert_near(configured.value().moments[1], 6.0);
    assert_near(configured.value().moments[2], 24.0);
    assert_near(configured.value().moments[3], 144.0);
}

void test_valid_cqmom_matches_programmatic_case()
{
    const auto configured = qmomkit::config::parse_cqmom_dense(cqmom_text, "cqmom-test");
    assert(configured.ok());
    const qmomkit::CqmomDenseCase direct{qmomkit::PhaseSpaceHierarchy({1, 1}, {1, 0}), {{2, 2}, {1.0, 1.5, 0.5, 0.75}}, {}, {true}, true};
    assert(configured.value().hierarchy.node_counts() == direct.hierarchy.node_counts());
    assert(configured.value().hierarchy.conditioning_order() == direct.hierarchy.conditioning_order());
    assert(configured.value().moments.extents == direct.moments.extents);
    assert(configured.value().moments.values == direct.moments.values);
    assert(configured.value().adaptive_branch_counts);
    assert(configured.value().qmom_options.adaptive);
    assert(!configured.value().qmom_options.automatic_scaling);

    const auto configured_result = qmomkit::compute_cqmom(configured.value());
    const auto direct_result = qmomkit::compute_cqmom(direct);
    assert(configured_result.ok());
    assert(direct_result.ok());
    assert(configured_result.quadrature.weights() == direct_result.quadrature.weights());
    for (std::size_t dimension = 0; dimension < 2; ++dimension)
    {
        assert_near(configured_result.quadrature.abscissas()[dimension][0], direct_result.quadrature.abscissas()[dimension][0]);
    }
}

void test_phase_space_automatic_scaling_parses_for_qmom_and_cqmom()
{
    const auto qmom = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1
automatic_scaling = true

[moments]
values = 1.0, 0.0
)");
    assert(qmom.ok());
    assert(qmom.value().options.automatic_scaling);

    const auto cqmom = qmomkit::config::parse_cqmom_lean(R"(
[phase_space]
node_counts = 1, 1
automatic_scaling = true

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.0
M[1,0] = 0.0
M[1,1] = 0.0
)");
    assert(cqmom.ok());
    assert(cqmom.value().qmom_options.automatic_scaling);

    const auto invalid = qmomkit::config::parse_qmom("[phase_space]\nnode_counts = 1\nautomatic_scaling = yes\n[moments]\nvalues = 1.0, 0.0\n");
    assert(!invalid.ok());
    assert(has_message_fragment(invalid, "automatic_scaling must be true or false"));
}

void test_qmom_automatic_scaling_handoff_to_evolution_options()
{
    // One .qmom file defines the reconstruction policy for both initial
    // inversion and subsequent QMOM post-source reconstructions.
    const std::string scaled_case = R"(
[phase_space]
node_counts = 1
automatic_scaling = true

[moments]
values = 1.0, 0.0
)";
    const auto qmom_case = qmomkit::config::parse_qmom(scaled_case);
    const auto evolution = qmomkit::config::parse_qmom_phase_space_evolution_options(scaled_case);
    assert(qmom_case.ok());
    assert(evolution.ok());
    assert(qmom_case.value().options.automatic_scaling);

    // Explicit false and an omitted key retain the documented default.
    const auto explicitly_unscaled =
        qmomkit::config::parse_qmom_phase_space_evolution_options("[phase_space]\nnode_counts = 1\nautomatic_scaling = false\n");
    const auto omitted = qmomkit::config::parse_qmom_phase_space_evolution_options("[phase_space]\nnode_counts = 1\n");
    assert(explicitly_unscaled.ok());
    assert(omitted.ok());
    assert(!explicitly_unscaled.value().qmom_options.automatic_scaling);
    assert(!omitted.value().qmom_options.automatic_scaling);
    assert(evolution.value().qmom_options.automatic_scaling);
}

void test_cqmom_independent_product_initialization_matches_direct_api()
{
    const auto configured = qmomkit::config::parse_cqmom_dense(R"(
[phase_space]
names = size, velocity
node_counts = 2, 2
conditioning_order = 1, 0
adaptive_branch_counts = true

[qmom]
adaptive = true

[initial]
type = independent_product
number_density = 2.5
maximum_orders = 3, 3

[initial.dimension.0]
family = gamma
shape = 3.0
scale = 0.2

[initial.dimension.1]
family = normal
mean = -0.5
standard_deviation = 1.25
)");
    assert(configured.ok());

    const qmomkit::IndependentProductDistribution direct{{qmomkit::GammaDistribution{3.0, 0.2}, qmomkit::NormalDistribution{-0.5, 1.25}}, 2.5};
    const qmomkit::TensorMomentResult expected = qmomkit::compute_raw_moments(direct, {3, 3});
    assert(expected.ok());
    assert(configured.value().hierarchy.node_counts() == std::vector<std::size_t>({2, 2}));
    assert(configured.value().hierarchy.conditioning_order() == std::vector<std::size_t>({1, 0}));
    assert(configured.value().adaptive_branch_counts);
    assert(configured.value().qmom_options.adaptive);
    assert_tensor_near(configured.value().moments, expected.moments);
    assert_near(configured.value().moments.values[qmomkit::flatten_multi_index({0, 0}, configured.value().moments.extents)], 2.5);
    assert(configured.value().moments.extents == std::vector<std::size_t>({4, 4}));

    const auto initial = qmomkit::config::parse_cqmom_initial_condition(R"(
[phase_space]
node_counts = 1, 1

[initial]
type = independent_product
number_density = 4.0
maximum_orders = 1, 1

[initial.dimension.0]
family = beta
alpha = 2.0
beta = 5.0

[initial.dimension.1]
family = lognormal
mu_log = 0.1
sigma_log = 0.2
)");
    assert(initial.ok());
    const auto tensor = qmomkit::config::initialize_cqmom_dense_moments(initial.value());
    const qmomkit::TensorMomentResult expected_tensor = qmomkit::compute_raw_moments(
        qmomkit::IndependentProductDistribution{{qmomkit::BetaDistribution{2.0, 5.0}, qmomkit::LogNormalDistribution{0.1, 0.2}}, 4.0}, {1, 1});
    assert(expected_tensor.ok());
    assert(tensor.ok());
    assert_tensor_near(tensor.value(), expected_tensor.moments);

    const auto dense_case = qmomkit::config::create_cqmom_dense_case(initial.value(), qmomkit::PhaseSpaceHierarchy({1, 1}));
    assert(dense_case.ok());
    assert_tensor_near(dense_case.value().moments, expected_tensor.moments);
}

void test_cqmom_independent_product_mixture_marginal_initialization()
{
    const auto initial = qmomkit::config::parse_cqmom_initial_condition(R"(
[phase_space]
node_counts = 1, 1

[initial]
type = independent_product
number_density = 2.75
maximum_orders = 1, 2

[initial.dimension.0]
family = normal
mean = -0.4
standard_deviation = 0.7

[initial.dimension.1]
type = mixture

[initial.dimension.1.mixture.component.0]
family = lognormal
number_density = 1.5
mu_log = -0.3
sigma_log = 0.25

[initial.dimension.1.mixture.component.1]
family = lognormal
number_density = 3.25
mu_log = 0.4
sigma_log = 0.45
)");
    assert(initial.ok());

    const auto* independent_product = std::get_if<qmomkit::config::IndependentProductInitialization>(&initial.value());
    assert(independent_product != nullptr);
    assert(independent_product->distribution.number_density == 2.75);
    assert(independent_product->distribution.marginals.size() == 2);
    const auto* normal = std::get_if<qmomkit::DistributionSpecification>(&independent_product->distribution.marginals[0]);
    const auto* mixture = std::get_if<qmomkit::DistributionMixture>(&independent_product->distribution.marginals[1]);
    assert(normal != nullptr);
    assert(mixture != nullptr);
    assert(std::get<qmomkit::NormalDistribution>(*normal).mean == -0.4);
    assert(std::get<qmomkit::NormalDistribution>(*normal).standard_deviation == 0.7);
    assert(mixture->components.size() == 2);
    assert(mixture->components[0].number_density == 1.5);
    assert(mixture->components[1].number_density == 3.25);
    assert(std::get<qmomkit::LogNormalDistribution>(mixture->components[0].distribution).normal_mean == -0.3);
    assert(std::get<qmomkit::LogNormalDistribution>(mixture->components[1].distribution).normal_standard_deviation == 0.45);

    const qmomkit::IndependentProductDistribution direct{{qmomkit::DistributionSpecification{qmomkit::NormalDistribution{-0.4, 0.7}}, *mixture},
                                                         2.75};
    const qmomkit::TensorMomentResult expected = qmomkit::compute_raw_moments(direct, {1, 2});
    const auto initialized = qmomkit::config::initialize_cqmom_dense_moments(initial.value());
    assert(expected.ok());
    assert(initialized.ok());
    assert_tensor_near(initialized.value(), expected.moments);
}

void test_lean_cqmom_initialization_matches_dense_subselection()
{
    const auto initial = qmomkit::config::parse_cqmom_initial_condition(R"(
[phase_space]
node_counts = 2, 2
conditioning_order = 1, 0

[initial]
type = independent_product
number_density = 2.5
maximum_orders = 3, 3

[initial.dimension.0]
family = gamma
shape = 3.0
scale = 0.2

[initial.dimension.1]
family = normal
mean = -0.5
standard_deviation = 1.25
)");
    assert(initial.ok());

    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {1, 0});
    const auto dense = qmomkit::config::initialize_cqmom_dense_moments(initial.value());
    const auto lean = qmomkit::config::initialize_cqmom_moment_data(initial.value(), hierarchy);
    assert(dense.ok());
    assert(lean.ok());

    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments(dense.value(), layout);
    assert_cqmom_moments_near(lean.value(), expected);
    assert_near(lean.value().moment({0, 0}), 2.5);
}

void test_lean_cqmom_initialization_matches_example_flow()
{
    const auto initial = qmomkit::config::parse_cqmom_initial_condition(R"(
[phase_space]
node_counts = 2, 2
conditioning_order = 0, 1

[initial]
type = independent_product
number_density = 2.5
maximum_orders = 3, 3

[initial.dimension.0]
family = gamma
shape = 3.0
scale = 0.2

[initial.dimension.1]
family = normal
mean = -0.5
standard_deviation = 1.25
)");
    assert(initial.ok());

    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {0, 1});
    const auto configured = qmomkit::config::initialize_cqmom_moment_data(initial.value(), hierarchy);
    assert(configured.ok());

    const qmomkit::IndependentProductDistribution distribution{{qmomkit::GammaDistribution{3.0, 0.2}, qmomkit::NormalDistribution{-0.5, 1.25}}, 2.5};
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::RawMomentRequest moment_request{layout.required_indices};
    const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(distribution, moment_request);
    assert(selected.ok());
    const qmomkit::CqmomMomentData expected = qmomkit::create_cqmom_moment_data(layout, selected);

    assert_cqmom_moments_near(configured.value(), expected);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(qmomkit::CqmomLeanCase{configured.value()});
    assert(result.ok());
}

void test_lean_cqmom_initialization_omits_maximum_orders()
{
    const std::string text = R"(
[phase_space]
names = size, velocity
node_counts = 2, 2
conditioning_order = 0, 1

[initial]
type = independent_product
number_density = 2.5

[initial.dimension.0]
family = gamma
shape = 3.0
scale = 0.2

[initial.dimension.1]
family = normal
mean = -0.5
standard_deviation = 1.25
)";

    const auto initial = qmomkit::config::parse_cqmom_initial_condition(text);
    assert(initial.ok());
    const auto* independent_product = std::get_if<qmomkit::config::IndependentProductInitialization>(&initial.value());
    assert(independent_product != nullptr);
    assert(!independent_product->maximum_orders.has_value());

    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {0, 1});
    const auto initialized = qmomkit::config::initialize_cqmom_moment_data(initial.value(), hierarchy);
    assert(initialized.ok());
    assert(initialized.value().layout.hierarchy.node_counts() == std::vector<std::size_t>({2, 2}));
    assert_near(initialized.value().moment({0, 0}), 2.5);

    const auto lean = qmomkit::config::parse_cqmom_lean(text);
    const auto case_result = qmomkit::config::parse_cqmom(text);
    assert(lean.ok());
    assert(case_result.ok());
    const qmomkit::CqmomResult lean_result = qmomkit::compute_cqmom(lean.value());
    const qmomkit::CqmomResult case_compute = qmomkit::compute_cqmom(case_result.value());
    assert(lean_result.ok());
    assert(case_compute.ok());
    assert_cqmom_quadratures_near(case_compute.quadrature, lean_result.quadrature);
}

void test_dense_cqmom_initialization_omits_maximum_orders()
{
    const std::string omitted_order_text = R"(
[phase_space]
node_counts = 1, 1
conditioning_order = 1, 0
support_lower = none, -2.0
support_upper = none, 2.0

[initial]
type = independent_product
number_density = 2.5

[initial.dimension.0]
family = beta
alpha = 2.0
beta = 5.0

[initial.dimension.1]
family = lognormal
mu_log = 0.1
sigma_log = 0.2
)";
    const std::string explicit_order_text = R"(
[phase_space]
node_counts = 1, 1
conditioning_order = 1, 0
support_lower = none, -2.0
support_upper = none, 2.0

[initial]
type = independent_product
number_density = 2.5
maximum_orders = 1, 1

[initial.dimension.0]
family = beta
alpha = 2.0
beta = 5.0

[initial.dimension.1]
family = lognormal
mu_log = 0.1
sigma_log = 0.2
)";

    const auto initial = qmomkit::config::parse_cqmom_initial_condition(omitted_order_text);
    assert(initial.ok());

    const qmomkit::PhaseSpaceHierarchy hierarchy({1, 1}, {1, 0});
    const std::vector<std::size_t> extents = hierarchy.moment_extents();
    const std::vector<std::size_t> inferred_orders{1, 1};
    assert(inferred_orders.size() == extents.size());
    for (std::size_t dimension = 0; dimension < extents.size(); ++dimension)
    {
        assert(inferred_orders[dimension] == extents[dimension] - 1);
    }

    const auto dense_case = qmomkit::config::create_cqmom_dense_case(initial.value(), hierarchy);
    assert(dense_case.ok());

    const qmomkit::IndependentProductDistribution direct{{qmomkit::BetaDistribution{2.0, 5.0}, qmomkit::LogNormalDistribution{0.1, 0.2}}, 2.5};
    const qmomkit::TensorMomentResult expected = qmomkit::compute_raw_moments(direct, inferred_orders);
    assert(expected.ok());
    assert_tensor_near(dense_case.value().moments, expected.moments);

    const auto omitted_order = qmomkit::config::parse_cqmom_dense(omitted_order_text);
    const auto explicit_order = qmomkit::config::parse_cqmom_dense(explicit_order_text);
    assert(omitted_order.ok());
    assert(explicit_order.ok());
    assert(omitted_order.value().hierarchy.node_counts() == explicit_order.value().hierarchy.node_counts());
    assert(omitted_order.value().hierarchy.conditioning_order() == explicit_order.value().hierarchy.conditioning_order());
    assert(omitted_order.value().coordinate_support.size() == explicit_order.value().coordinate_support.size());
    assert(!omitted_order.value().coordinate_support[0].has_value());
    assert(!explicit_order.value().coordinate_support[0].has_value());
    assert(omitted_order.value().coordinate_support[1].has_value());
    assert(explicit_order.value().coordinate_support[1].has_value());
    assert_near(omitted_order.value().coordinate_support[1]->lower, explicit_order.value().coordinate_support[1]->lower);
    assert_near(omitted_order.value().coordinate_support[1]->upper, explicit_order.value().coordinate_support[1]->upper);
    assert(omitted_order.value().moments.extents == explicit_order.value().moments.extents);
    assert(omitted_order.value().moments.values.size() == explicit_order.value().moments.values.size());
    assert_tensor_near(omitted_order.value().moments, explicit_order.value().moments);
    assert_tensor_near(omitted_order.value().moments, expected.moments);
}

void test_lean_cqmom_initialization_selects_explicit_dense_tensor()
{
    const auto initial = qmomkit::config::parse_cqmom_initial_condition(cqmom_text, "cqmom-test");
    assert(initial.ok());

    const qmomkit::PhaseSpaceHierarchy hierarchy({1, 1}, {1, 0});
    const auto lean = qmomkit::config::initialize_cqmom_moment_data(initial.value(), hierarchy);
    assert(lean.ok());

    const qmomkit::RawMomentTensor dense{{2, 2}, {1.0, 1.5, 0.5, 0.75}};
    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments(dense, qmomkit::create_cqmom_moment_layout(hierarchy));
    assert_cqmom_moments_near(lean.value(), expected);
}

void test_cqmom_lean_case_uses_selected_analytical_moments()
{
    const auto lean = qmomkit::config::parse_cqmom_lean(R"(
[phase_space]
names = size, velocity
node_counts = 2, 2
conditioning_order = 0, 1

[cqmom]
diagnostic_level = basic

[initial]
type = independent_product
number_density = 2.5
maximum_orders = 1, 1

[initial.dimension.0]
family = gamma
shape = 3.0
scale = 0.2

[initial.dimension.1]
family = normal
mean = -0.5
standard_deviation = 1.25
)");
    assert(lean.ok());
    assert(lean.value().moments.layout.hierarchy.node_counts() == std::vector<std::size_t>({2, 2}));
    assert(lean.value().moments.layout.hierarchy.conditioning_order() == std::vector<std::size_t>({0, 1}));
    assert(lean.value().diagnostic_level == qmomkit::CqmomDiagnosticLevel::basic);
    assert_near(lean.value().moments.moment({0, 0}), 2.5);

    const auto dense = qmomkit::config::parse_cqmom_dense(R"(
[phase_space]
node_counts = 2, 2

[initial]
type = independent_product
number_density = 2.5
maximum_orders = 1, 1

[initial.dimension.0]
family = gamma
shape = 3.0
scale = 0.2

[initial.dimension.1]
family = normal
mean = -0.5
standard_deviation = 1.25
)");
    assert(!dense.ok());
    assert(has_failure(dense));

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(lean.value());
    assert(result.ok());
}

void test_cqmom_lean_case_matches_dense_case_path()
{
    const std::string text = R"(
[phase_space]
names = size, velocity
node_counts = 2, 2
conditioning_order = 1, 0

[initial]
type = independent_product
number_density = 2.5
maximum_orders = 3, 3

[initial.dimension.0]
family = gamma
shape = 3.0
scale = 0.2

[initial.dimension.1]
family = normal
mean = -0.5
standard_deviation = 1.25
)";

    const auto dense = qmomkit::config::parse_cqmom_dense(text);
    const auto lean = qmomkit::config::parse_cqmom_lean(text);
    assert(dense.ok());
    assert(lean.ok());

    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(dense.value().hierarchy);
    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments(dense.value().moments, layout);
    assert_cqmom_moments_near(lean.value().moments, expected);

    const qmomkit::CqmomResult dense_result = qmomkit::compute_cqmom(dense.value());
    const qmomkit::CqmomResult lean_result = qmomkit::compute_cqmom(lean.value());
    assert(dense_result.ok());
    assert(lean_result.ok());
    assert(dense_result.quadrature.weights().size() == lean_result.quadrature.weights().size());
    for (std::size_t node = 0; node < dense_result.quadrature.weights().size(); ++node)
    {
        assert_near(lean_result.quadrature.weights()[node], dense_result.quadrature.weights()[node]);
        for (std::size_t dimension = 0; dimension < dense_result.quadrature.number_of_dimensions(); ++dimension)
        {
            assert_near(lean_result.quadrature.abscissas()[dimension][node], dense_result.quadrature.abscissas()[dimension][node]);
        }
    }
}

void test_cqmom_lean_case_selects_explicit_dense_tensor()
{
    const auto lean = qmomkit::config::parse_cqmom_lean(cqmom_text, "cqmom-test");
    assert(lean.ok());
    assert(lean.value().adaptive_branch_counts);
    assert(lean.value().qmom_options.adaptive);

    const qmomkit::PhaseSpaceHierarchy hierarchy({1, 1}, {1, 0});
    const qmomkit::RawMomentTensor dense{{2, 2}, {1.0, 1.5, 0.5, 0.75}};
    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments(dense, qmomkit::create_cqmom_moment_layout(hierarchy));
    assert_cqmom_moments_near(lean.value().moments, expected);

    const auto initial = qmomkit::config::parse_cqmom_initial_condition(cqmom_text, "cqmom-test");
    assert(initial.ok());
    const auto made = qmomkit::config::create_cqmom_lean_case(initial.value(), hierarchy, {true}, true);
    assert(made.ok());
    assert_cqmom_moments_near(made.value().moments, expected);
    assert(made.value().adaptive_branch_counts);
    assert(made.value().qmom_options.adaptive);
}

void test_cqmom_case_uses_selected_analytical_moments()
{
    const std::string text = R"(
[phase_space]
names = size, velocity
node_counts = 2, 2
conditioning_order = 1, 0

[cqmom]
diagnostic_level = basic

[initial]
type = independent_product
number_density = 2.5
maximum_orders = 1, 1

[initial.dimension.0]
family = gamma
shape = 3.0
scale = 0.2

[initial.dimension.1]
family = normal
mean = -0.5
standard_deviation = 1.25
)";

    const auto case_result = qmomkit::config::parse_cqmom(text);
    const auto lean = qmomkit::config::parse_cqmom_lean(text);
    assert(case_result.ok());
    assert(lean.ok());

    const qmomkit::CqmomCase cqmom_case = case_result.value();
    assert(cqmom_case.hierarchy().node_counts() == std::vector<std::size_t>({2, 2}));
    assert(cqmom_case.hierarchy().conditioning_order() == std::vector<std::size_t>({1, 0}));
    assert(cqmom_case.diagnostic_level() == qmomkit::CqmomDiagnosticLevel::basic);

    const qmomkit::CqmomResult case_compute = qmomkit::compute_cqmom(cqmom_case);
    const qmomkit::CqmomResult lean_compute = qmomkit::compute_cqmom(lean.value());
    assert(case_compute.ok());
    assert(lean_compute.ok());
    assert_cqmom_quadratures_near(case_compute.quadrature, lean_compute.quadrature);
}

void test_cqmom_case_accepts_explicit_dense_tensor()
{
    const auto case_result = qmomkit::config::parse_cqmom(cqmom_text, "cqmom-test");
    const auto dense = qmomkit::config::parse_cqmom_dense(cqmom_text, "cqmom-test");
    const auto lean = qmomkit::config::parse_cqmom_lean(cqmom_text, "cqmom-test");
    assert(case_result.ok());
    assert(dense.ok());
    assert(lean.ok());

    const qmomkit::CqmomCase cqmom_case = case_result.value();
    assert(cqmom_case.moment_storage() == qmomkit::CqmomMomentStorage::lean);
    assert(cqmom_case.hierarchy().node_counts() == std::vector<std::size_t>({1, 1}));
    assert(cqmom_case.hierarchy().conditioning_order() == std::vector<std::size_t>({1, 0}));
    assert(cqmom_case.adaptive_branch_counts());
    assert(cqmom_case.qmom_options().adaptive);

    const qmomkit::CqmomResult case_compute = qmomkit::compute_cqmom(cqmom_case);
    const qmomkit::CqmomResult dense_compute = qmomkit::compute_cqmom(dense.value());
    const qmomkit::CqmomResult lean_compute = qmomkit::compute_cqmom(lean.value());
    assert(case_compute.ok());
    assert(dense_compute.ok());
    assert(lean_compute.ok());
    assert_cqmom_quadratures_near(case_compute.quadrature, dense_compute.quadrature);
    assert_cqmom_quadratures_near(case_compute.quadrature, lean_compute.quadrature);
}

void test_cqmom_case_parses_coordinate_support()
{
    const std::string text = std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, 1.0
support_upper = none, 2.0

[cqmom]
diagnostic_level = basic
)";
    const auto case_result = qmomkit::config::parse_cqmom(text, "cqmom-support");
    const auto dense = qmomkit::config::parse_cqmom_dense(text, "cqmom-support");
    const auto lean = qmomkit::config::parse_cqmom_lean(text, "cqmom-support");
    assert(case_result.ok());
    assert(dense.ok());
    assert(lean.ok());

    const qmomkit::CqmomCase cqmom_case = case_result.value();
    assert(cqmom_case.coordinate_support().size() == 2);
    assert(!cqmom_case.coordinate_support()[0].has_value());
    assert(cqmom_case.coordinate_support()[1].has_value());
    assert_near(cqmom_case.coordinate_support()[1]->lower, 1.0);
    assert_near(cqmom_case.coordinate_support()[1]->upper, 2.0);
    assert(dense.value().coordinate_support.size() == cqmom_case.coordinate_support().size());
    assert(lean.value().coordinate_support.size() == cqmom_case.coordinate_support().size());

    const qmomkit::CqmomResult computed = qmomkit::compute_cqmom(cqmom_case);
    assert(computed.ok());
    assert(!computed.diagnostics.empty());
    assert(computed.diagnostics.front().phase_space_dimension == 1);
    assert(computed.diagnostics.front().bounded_support_checked);

    const qmomkit::CqmomState state = cqmom_case.initial_state();
    assert(state.coordinate_support().size() == 2);
    assert(state.coordinate_support()[1].has_value());
    assert_near(state.coordinate_support()[1]->lower, 1.0);
    assert_near(state.coordinate_support()[1]->upper, 2.0);
}

void test_cqmom_factories_preserve_coordinate_support()
{
    const auto initial = qmomkit::config::parse_cqmom_initial_condition(cqmom_text, "cqmom-factory-support");
    assert(initial.ok());

    const qmomkit::PhaseSpaceHierarchy hierarchy({1, 1}, {1, 0});
    const auto dense_default = qmomkit::config::create_cqmom_dense_case(initial.value(), hierarchy);
    const auto lean_default = qmomkit::config::create_cqmom_lean_case(initial.value(), hierarchy);
    assert(dense_default.ok());
    assert(lean_default.ok());
    assert(dense_default.value().coordinate_support.empty());
    assert(lean_default.value().coordinate_support.empty());

    qmomkit::CqmomCoordinateSupport support(2);
    support[1] = qmomkit::AbscissaBounds{2.0, 5.0};

    const auto dense = qmomkit::config::create_cqmom_dense_case(initial.value(), hierarchy, qmomkit::QmomOptions{}, false, "dense-factory", support);
    const auto lean = qmomkit::config::create_cqmom_lean_case(initial.value(), hierarchy, qmomkit::QmomOptions{}, false, "lean-factory", support);
    assert(dense.ok());
    assert(lean.ok());
    assert_coordinate_support_matches_shifted_interval(dense.value().coordinate_support);
    assert_coordinate_support_matches_shifted_interval(lean.value().coordinate_support);
}

void test_cqmom_case_rejects_invalid_coordinate_support()
{
    const auto missing_upper = qmomkit::config::parse_cqmom(std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, 1.0
)",
                                                            "cqmom-support");
    assert(!missing_upper.ok());
    assert(has_message_fragment(missing_upper, "support_lower and phase_space.support_upper must be configured together"));

    const auto wrong_size = qmomkit::config::parse_cqmom(std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, 1.0, 2.0
support_upper = none, 2.0, 3.0
)",
                                                         "cqmom-support");
    assert(!wrong_size.ok());
    assert(has_message_fragment(wrong_size, "one entry per phase-space dimension"));

    const auto invalid_interval = qmomkit::config::parse_cqmom(std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, 2.0
support_upper = none, 1.0
)",
                                                               "cqmom-support");
    assert(!invalid_interval.ok());
    assert(has_message_fragment(invalid_interval, "lower < upper"));

    const auto lower_nan = qmomkit::config::parse_cqmom(std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, nan
support_upper = none, 2.0
)",
                                                        "cqmom-support");
    assert(!lower_nan.ok());
    assert(has_message_fragment(lower_nan, "support_lower must contain finite numbers or none entries"));

    const auto upper_inf = qmomkit::config::parse_cqmom(std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, 1.0
support_upper = none, inf
)",
                                                        "cqmom-support");
    assert(!upper_inf.ok());
    assert(has_message_fragment(upper_inf, "support_upper must contain finite numbers or none entries"));

    const auto lower_negative_inf = qmomkit::config::parse_cqmom(std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, -inf
support_upper = none, 2.0
)",
                                                                 "cqmom-support");
    assert(!lower_negative_inf.ok());
    assert(has_message_fragment(lower_negative_inf, "support_lower must contain finite numbers or none entries"));
}

void test_cqmom_case_moment_storage_selection()
{
    std::string dense_text = cqmom_text;
    dense_text += "\n[cqmom]\nmoment_storage = dense\n";
    const auto dense_case_result = qmomkit::config::parse_cqmom(dense_text, "cqmom-dense-storage");
    const auto explicit_dense = qmomkit::config::parse_cqmom_dense(dense_text, "cqmom-dense-storage");
    const auto conflicting_lean = qmomkit::config::parse_cqmom_lean(dense_text, "cqmom-dense-storage");
    assert(dense_case_result.ok());
    assert(explicit_dense.ok());
    assert(!conflicting_lean.ok());
    assert(has_message_fragment(conflicting_lean, "parse_cqmom_lean requires cqmom.moment_storage = lean"));

    const qmomkit::CqmomCase dense_case = dense_case_result.value();
    assert(dense_case.moment_storage() == qmomkit::CqmomMomentStorage::dense);
    assert(dense_case.initial_state().moment_storage() == qmomkit::CqmomMomentStorage::dense);

    std::string lean_text = cqmom_text;
    lean_text += "\n[cqmom]\nmoment_storage = lean\n";
    const auto lean_case_result = qmomkit::config::parse_cqmom(lean_text, "cqmom-lean-storage");
    const auto explicit_lean = qmomkit::config::parse_cqmom_lean(lean_text, "cqmom-lean-storage");
    const auto conflicting_dense = qmomkit::config::parse_cqmom_dense(lean_text, "cqmom-lean-storage");
    assert(lean_case_result.ok());
    assert(explicit_lean.ok());
    assert(!conflicting_dense.ok());
    assert(has_message_fragment(conflicting_dense, "parse_cqmom_dense requires cqmom.moment_storage = dense"));
    assert(lean_case_result.value().moment_storage() == qmomkit::CqmomMomentStorage::lean);
    assert(lean_case_result.value().initial_state().moment_storage() == qmomkit::CqmomMomentStorage::lean);

    std::string invalid_text = cqmom_text;
    invalid_text += "\n[cqmom]\nmoment_storage = compact\n";
    const auto invalid = qmomkit::config::parse_cqmom(invalid_text, "cqmom-invalid-storage");
    assert(!invalid.ok());
    assert(has_message_fragment(invalid, "cqmom.moment_storage must be lean or dense"));
}

void test_cqmom_case_initial_state_preserves_case_settings()
{
    const auto case_result = qmomkit::config::parse_cqmom(cqmom_text, "cqmom-test");
    assert(case_result.ok());

    const qmomkit::CqmomCase cqmom_case = case_result.value();
    const qmomkit::CqmomState state = cqmom_case.initial_state();

    assert(state.hierarchy().node_counts() == cqmom_case.hierarchy().node_counts());
    assert(state.hierarchy().conditioning_order() == cqmom_case.hierarchy().conditioning_order());
    assert(state.qmom_options().adaptive == cqmom_case.qmom_options().adaptive);
    assert_near(state.qmom_options().rmin, cqmom_case.qmom_options().rmin);
    assert_near(state.qmom_options().eabs, cqmom_case.qmom_options().eabs);
    assert_near(state.qmom_options().cutoff, cqmom_case.qmom_options().cutoff);
    assert_near(state.conditional_options().singular_tolerance, cqmom_case.conditional_options().singular_tolerance);
    assert_near(state.conditional_options().conditioning_warning_threshold, cqmom_case.conditional_options().conditioning_warning_threshold);
    assert(state.adaptive_branch_counts() == cqmom_case.adaptive_branch_counts());
    assert(state.diagnostic_level() == cqmom_case.diagnostic_level());
    assert(state.vandermonde_solver() == cqmom_case.vandermonde_solver());
    assert(state.moment_storage() == cqmom_case.moment_storage());

    const qmomkit::CqmomState copied = state;
    assert(copied.hierarchy().node_counts() == state.hierarchy().node_counts());
    assert(copied.hierarchy().conditioning_order() == state.hierarchy().conditioning_order());
    assert(copied.qmom_options().adaptive == state.qmom_options().adaptive);
    assert(copied.vandermonde_solver() == state.vandermonde_solver());
}

qmomkit::RawMomentTensor cqmom_replacement_tensor(const std::vector<double>& x_values)
{
    return {{4, 4}, raw_moments_from_terminals_2d({0.2, 0.3, 0.1, 0.4}, {x_values, {0.25, 1.5, 0.25, 1.5}})};
}

qmomkit::CqmomCoordinateSupport cqmom_replacement_support()
{
    qmomkit::CqmomCoordinateSupport support(2);
    support[1] = qmomkit::AbscissaBounds{0.0, 2.0};
    return support;
}

void assert_cqmom_case_policies_equal(const qmomkit::CqmomCase& actual, const qmomkit::CqmomCase& expected)
{
    assert(actual.hierarchy().node_counts() == expected.hierarchy().node_counts());
    assert(actual.hierarchy().conditioning_order() == expected.hierarchy().conditioning_order());
    assert(actual.qmom_options().adaptive == expected.qmom_options().adaptive);
    assert_near(actual.qmom_options().rmin, expected.qmom_options().rmin);
    assert_near(actual.qmom_options().eabs, expected.qmom_options().eabs);
    assert_near(actual.qmom_options().cutoff, expected.qmom_options().cutoff);
    assert_near(actual.conditional_options().singular_tolerance, expected.conditional_options().singular_tolerance);
    assert_near(actual.conditional_options().conditioning_warning_threshold, expected.conditional_options().conditioning_warning_threshold);
    assert(actual.adaptive_branch_counts() == expected.adaptive_branch_counts());
    assert(actual.diagnostic_level() == expected.diagnostic_level());
    assert(actual.vandermonde_solver() == expected.vandermonde_solver());
    assert(actual.coordinate_support().size() == expected.coordinate_support().size());
    for (std::size_t dimension = 0; dimension < actual.coordinate_support().size(); ++dimension)
    {
        assert(actual.coordinate_support()[dimension].has_value() == expected.coordinate_support()[dimension].has_value());
        if (actual.coordinate_support()[dimension].has_value())
        {
            assert_near(actual.coordinate_support()[dimension]->lower, expected.coordinate_support()[dimension]->lower);
            assert_near(actual.coordinate_support()[dimension]->upper, expected.coordinate_support()[dimension]->upper);
        }
    }
}

qmomkit::CqmomDenseCase cqmom_dense_replacement_base()
{
    return {qmomkit::PhaseSpaceHierarchy({2, 2}, {1, 0}),
            cqmom_replacement_tensor({0.75, 0.75, 1.5, 1.5}),
            {1.0e-12, 1.0e9},
            {false, 1.0e-9, 1.0e-11, 1.0e-28},
            true,
            qmomkit::CqmomDiagnosticLevel::basic,
            qmomkit::VandermondeSolver::dense,
            cqmom_replacement_support()};
}

void test_cqmom_case_replaces_dense_initial_moments_immutably()
{
    const qmomkit::CqmomDenseCase dense_base = cqmom_dense_replacement_base();
    const qmomkit::CqmomCase base_case = qmomkit::internal::CqmomConstructionAccess::create_case(dense_base);
    qmomkit::RawMomentTensor replacement = cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0});
    const qmomkit::RawMomentTensor expected_replacement = replacement;

    const auto transformed = qmomkit::config::with_initial_moments(base_case, replacement);
    assert(transformed.ok());
    const qmomkit::CqmomCase replaced_case = transformed.value();
    assert(replaced_case.moment_storage() == qmomkit::CqmomMomentStorage::dense);
    assert_cqmom_case_policies_equal(replaced_case, base_case);

    replacement.values[0] = -999.0;
    const qmomkit::CqmomState base_state = base_case.initial_state();
    const qmomkit::CqmomState replaced_state = replaced_case.initial_state();
    assert(base_state.moment_storage() == qmomkit::CqmomMomentStorage::dense);
    assert(replaced_state.moment_storage() == qmomkit::CqmomMomentStorage::dense);
    assert(replaced_state.moment_indices() == qmomkit::enumerate_multi_indices(expected_replacement.extents));

    const qmomkit::CqmomResult base_result = qmomkit::compute_cqmom(base_case);
    const qmomkit::CqmomResult replaced_result = qmomkit::compute_cqmom(replaced_case);
    const qmomkit::CqmomResult expected_base_result = qmomkit::compute_cqmom(dense_base);
    const qmomkit::CqmomResult expected_result = qmomkit::compute_cqmom(qmomkit::CqmomDenseCase{
        dense_base.hierarchy, expected_replacement, dense_base.conditional_options, dense_base.qmom_options, dense_base.adaptive_branch_counts,
        dense_base.diagnostic_level, dense_base.vandermonde_solver, dense_base.coordinate_support});
    assert(base_result.ok());
    assert(replaced_result.ok());
    assert(expected_base_result.ok());
    assert(expected_result.ok());
    assert_cqmom_quadratures_near(base_result.quadrature, expected_base_result.quadrature);
    assert_cqmom_quadratures_near(replaced_result.quadrature, expected_result.quadrature);

    const qmomkit::RawMomentTensor canonical_base = qmomkit::reconstruct_raw_moments(base_result.quadrature, dense_base.moments.extents);
    const qmomkit::RawMomentTensor canonical_replacement = qmomkit::reconstruct_raw_moments(replaced_result.quadrature, expected_replacement.extents);
    assert_tensor_near({dense_base.moments.extents, base_state.moment_values()}, canonical_base);
    assert_tensor_near({expected_replacement.extents, replaced_state.moment_values()}, canonical_replacement);
}

void test_cqmom_case_replaces_lean_initial_moments_immutably()
{
    const qmomkit::CqmomDenseCase dense_base = cqmom_dense_replacement_base();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(dense_base.hierarchy);
    const qmomkit::CqmomLeanCase lean_base{qmomkit::extract_cqmom_moments(dense_base.moments, layout),
                                           dense_base.conditional_options,
                                           dense_base.qmom_options,
                                           dense_base.adaptive_branch_counts,
                                           dense_base.diagnostic_level,
                                           dense_base.vandermonde_solver,
                                           dense_base.coordinate_support};
    const qmomkit::CqmomCase base_case = qmomkit::internal::CqmomConstructionAccess::create_case(lean_base);
    qmomkit::RawMomentTensor replacement = cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0});
    const qmomkit::CqmomMomentData expected_moments = qmomkit::extract_cqmom_moments(replacement, layout);

    const auto transformed = qmomkit::config::with_initial_moments(base_case, replacement);
    assert(transformed.ok());
    const qmomkit::CqmomCase replaced_case = transformed.value();
    assert(replaced_case.moment_storage() == qmomkit::CqmomMomentStorage::lean);
    assert_cqmom_case_policies_equal(replaced_case, base_case);

    replacement.values[0] = -999.0;
    const qmomkit::CqmomState base_state = base_case.initial_state();
    const qmomkit::CqmomState replaced_state = replaced_case.initial_state();
    assert(base_state.moment_storage() == qmomkit::CqmomMomentStorage::lean);
    assert(replaced_state.moment_storage() == qmomkit::CqmomMomentStorage::lean);
    assert(base_state.moment_indices() == layout.required_indices);
    assert(replaced_state.moment_indices() == expected_moments.layout.required_indices);

    const qmomkit::CqmomResult base_result = qmomkit::compute_cqmom(base_case);
    const qmomkit::CqmomResult replaced_result = qmomkit::compute_cqmom(replaced_case);
    assert(base_result.ok());
    assert(replaced_result.ok());
    const qmomkit::CqmomMomentData canonical_base =
        qmomkit::extract_cqmom_moments(qmomkit::reconstruct_raw_moments(base_result.quadrature, dense_base.moments.extents), layout);
    const qmomkit::CqmomMomentData canonical_replacement =
        qmomkit::extract_cqmom_moments(qmomkit::reconstruct_raw_moments(replaced_result.quadrature, replacement.extents), layout);
    assert(base_state.moment_values() == canonical_base.values);
    assert(replaced_state.moment_values() == canonical_replacement.values);
}

void test_fixed_cqmom_prepared_metadata_ownership()
{
    qmomkit::CqmomDenseCase dense_case = cqmom_dense_replacement_base();
    dense_case.adaptive_branch_counts = false;
    dense_case.qmom_options.adaptive = false;
    const qmomkit::CqmomCase fixed_case = qmomkit::internal::CqmomConstructionAccess::create_case(dense_case);
    const auto& fixed_metadata = qmomkit::internal::CqmomConstructionAccess::case_prepared_fixed_metadata(fixed_case);
    assert(fixed_metadata != nullptr);

    const qmomkit::CqmomCase copied_case = fixed_case;
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_fixed_metadata(copied_case).get() == fixed_metadata.get());

    qmomkit::CqmomCase movable_case = fixed_case;
    qmomkit::CqmomCase moved_case = std::move(movable_case);
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_fixed_metadata(moved_case).get() == fixed_metadata.get());

    const qmomkit::CqmomCase distinct_case = qmomkit::internal::CqmomConstructionAccess::create_case(dense_case);
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_fixed_metadata(distinct_case).get() != fixed_metadata.get());

    const auto replaced = qmomkit::config::with_initial_moments(fixed_case, dense_case.moments);
    assert(replaced.ok());
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_fixed_metadata(replaced.value()).get() == fixed_metadata.get());

    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(dense_case.hierarchy);
    const qmomkit::CqmomCase lean_case = qmomkit::internal::CqmomConstructionAccess::create_case(qmomkit::CqmomLeanCase{
        qmomkit::extract_cqmom_moments(dense_case.moments, layout), dense_case.conditional_options, dense_case.qmom_options,
        dense_case.adaptive_branch_counts, dense_case.diagnostic_level, dense_case.vandermonde_solver, dense_case.coordinate_support});
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_fixed_metadata(lean_case) != nullptr);
    const qmomkit::CqmomResult dense_result = qmomkit::compute_cqmom(fixed_case);
    const qmomkit::CqmomResult lean_result = qmomkit::compute_cqmom(lean_case);
    assert(dense_result.ok());
    assert(lean_result.ok());
    assert_cqmom_quadratures_near(lean_result.quadrature, dense_result.quadrature);

    qmomkit::CqmomDenseCase adaptive_dense_case = dense_case;
    adaptive_dense_case.adaptive_branch_counts = true;
    adaptive_dense_case.qmom_options.adaptive = true;
    const qmomkit::CqmomCase adaptive_case = qmomkit::internal::CqmomConstructionAccess::create_case(std::move(adaptive_dense_case));
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_fixed_metadata(adaptive_case) == nullptr);

    const qmomkit::CqmomState initial_state = fixed_case.initial_state();
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(initial_state).get() == fixed_metadata.get());
    const qmomkit::CqmomState copied_state = initial_state;
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(copied_state).get() == fixed_metadata.get());
    qmomkit::CqmomState movable_state = initial_state;
    qmomkit::CqmomState moved_state = std::move(movable_state);
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(moved_state).get() == fixed_metadata.get());

    std::optional<qmomkit::CqmomState> surviving_state;
    const qmomkit::internal::PreparedFixedCqmomMetadata* surviving_metadata = nullptr;
    {
        const qmomkit::CqmomCase temporary_case = qmomkit::internal::CqmomConstructionAccess::create_case(dense_case);
        surviving_metadata = qmomkit::internal::CqmomConstructionAccess::case_prepared_fixed_metadata(temporary_case).get();
        surviving_state.emplace(temporary_case.initial_state());
    }
    assert(surviving_state.has_value());
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(*surviving_state).get() == surviving_metadata);
    assert(surviving_state->hierarchy().node_counts() == dense_case.hierarchy.node_counts());
    assert(surviving_state->hierarchy().conditioning_order() == dense_case.hierarchy.conditioning_order());
}

void test_fixed_cqmom_prepared_metadata_is_concurrently_readable()
{
    qmomkit::CqmomDenseCase dense_case = cqmom_dense_replacement_base();
    dense_case.adaptive_branch_counts = false;
    dense_case.qmom_options.adaptive = false;
    const qmomkit::CqmomCase fixed_case = qmomkit::internal::CqmomConstructionAccess::create_case(std::move(dense_case));
    const qmomkit::CqmomResult reference = qmomkit::compute_cqmom(fixed_case);
    assert(reference.ok());

    std::vector<int> succeeded(4, 0);
    std::vector<std::thread> workers;
    workers.reserve(succeeded.size());
    for (std::size_t worker = 0; worker < succeeded.size(); ++worker)
    {
        workers.emplace_back(
            [&fixed_case, &reference, &succeeded, worker]()
            {
                for (std::size_t iteration = 0; iteration < 16; ++iteration)
                {
                    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(fixed_case);
                    assert(result.ok());
                    assert_cqmom_quadratures_near(result.quadrature, reference.quadrature);
                }
                succeeded[worker] = 1;
            });
    }
    for (std::thread& worker : workers)
    {
        worker.join();
    }
    assert(std::all_of(succeeded.begin(), succeeded.end(), [](const int value) { return value == 1; }));
}

void test_fixed_cqmom_state_prepared_metadata_invariant_and_compatibility()
{
    qmomkit::CqmomDenseCase dense_case = cqmom_dense_replacement_base();
    dense_case.adaptive_branch_counts = false;
    dense_case.qmom_options.adaptive = false;
    const qmomkit::CqmomCase fixed_case = qmomkit::internal::CqmomConstructionAccess::create_case(dense_case);
    const auto& prepared_fixed_metadata = qmomkit::internal::CqmomConstructionAccess::case_prepared_fixed_metadata(fixed_case);
    assert(prepared_fixed_metadata != nullptr);

    const qmomkit::CqmomState new_lineage_state = qmomkit::internal::CqmomConstructionAccess::create_state(dense_case);
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(new_lineage_state) != nullptr);
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(new_lineage_state).get() != prepared_fixed_metadata.get());

    const qmomkit::CqmomState descendant_state = qmomkit::internal::CqmomConstructionAccess::create_state(dense_case, prepared_fixed_metadata);
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(descendant_state).get() == prepared_fixed_metadata.get());

    const auto assert_rejected = [&prepared_fixed_metadata](const qmomkit::CqmomDenseCase& incompatible_case)
    {
        bool rejected = false;
        try
        {
            static_cast<void>(qmomkit::internal::CqmomConstructionAccess::create_state(incompatible_case, prepared_fixed_metadata));
        }
        catch (const std::invalid_argument& error)
        {
            rejected = true;
            assert(std::string(error.what()).find("Prepared fixed CQMOM metadata is incompatible") != std::string::npos);
        }
        assert(rejected);
    };

    qmomkit::CqmomDenseCase different_node_counts = dense_case;
    different_node_counts.hierarchy = qmomkit::PhaseSpaceHierarchy({1, 2});
    assert_rejected(different_node_counts);

    qmomkit::CqmomDenseCase different_conditioning_order = dense_case;
    different_conditioning_order.hierarchy = qmomkit::PhaseSpaceHierarchy({2, 2}, {0, 1});
    assert_rejected(different_conditioning_order);

    qmomkit::CqmomDenseCase different_dimension = dense_case;
    different_dimension.hierarchy = qmomkit::PhaseSpaceHierarchy({1, 1, 1});
    assert_rejected(different_dimension);

    const qmomkit::CqmomMomentLayout fixed_layout = qmomkit::create_cqmom_moment_layout(dense_case.hierarchy);
    qmomkit::CqmomMomentLayout reordered_layout = fixed_layout;
    std::reverse(reordered_layout.required_indices.begin(), reordered_layout.required_indices.end());
    const qmomkit::CqmomLeanCase reordered_lean_case{
        qmomkit::CqmomMomentData{std::move(reordered_layout), std::vector<double>(fixed_layout.required_indices.size(), 1.0)}};
    bool lean_rejected = false;
    try
    {
        static_cast<void>(qmomkit::internal::CqmomConstructionAccess::create_state(reordered_lean_case, prepared_fixed_metadata));
    }
    catch (const std::invalid_argument& error)
    {
        lean_rejected = true;
        assert(std::string(error.what()).find("Prepared fixed CQMOM metadata is incompatible") != std::string::npos);
    }
    assert(lean_rejected);

    bool malformed_new_lineage_rejected = false;
    try
    {
        static_cast<void>(qmomkit::internal::CqmomConstructionAccess::create_state(reordered_lean_case));
    }
    catch (const std::invalid_argument& error)
    {
        malformed_new_lineage_rejected = true;
        assert(std::string(error.what()).find("canonical CQMOM moment layout") != std::string::npos);
    }
    assert(malformed_new_lineage_rejected);

    qmomkit::CqmomDenseCase adaptive_case = dense_case;
    adaptive_case.adaptive_branch_counts = true;
    adaptive_case.qmom_options.adaptive = true;
    bool adaptive_rejected = false;
    try
    {
        static_cast<void>(qmomkit::internal::CqmomConstructionAccess::create_state(adaptive_case, prepared_fixed_metadata));
    }
    catch (const std::invalid_argument& error)
    {
        adaptive_rejected = true;
        assert(std::string(error.what()).find("Adaptive CQMOM cases must not carry") != std::string::npos);
    }
    assert(adaptive_rejected);
    const qmomkit::CqmomState adaptive_state = qmomkit::internal::CqmomConstructionAccess::create_state(std::move(adaptive_case));
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_fixed_metadata(adaptive_state) == nullptr);
}

void test_cqmom_case_replacement_rejects_structural_tensor_errors()
{
    const qmomkit::CqmomDenseCase dense_base = cqmom_dense_replacement_base();
    const qmomkit::CqmomCase dense_case = qmomkit::internal::CqmomConstructionAccess::create_case(dense_base);
    const qmomkit::CqmomLeanCase lean_base{
        qmomkit::extract_cqmom_moments(dense_base.moments, qmomkit::create_cqmom_moment_layout(dense_base.hierarchy)),
        dense_base.conditional_options,
        dense_base.qmom_options,
        dense_base.adaptive_branch_counts,
        dense_base.diagnostic_level,
        dense_base.vandermonde_solver,
        dense_base.coordinate_support};
    const qmomkit::CqmomCase lean_case = qmomkit::internal::CqmomConstructionAccess::create_case(lean_base);
    const std::vector<double> base_values = dense_case.initial_state().moment_values();
    const std::vector<double> lean_base_values = lean_case.initial_state().moment_values();

    const auto rank = qmomkit::config::with_initial_moments(dense_case, {{4}, std::vector<double>(4, 1.0)});
    assert(!rank.ok());
    assert(has_message_fragment(rank, "rank"));

    const auto invalid_extents = qmomkit::config::with_initial_moments(dense_case, {{4, 0}, {}});
    assert(!invalid_extents.ok());
    assert(has_message_fragment(invalid_extents, "extents"));

    const auto value_count = qmomkit::config::with_initial_moments(dense_case, {{4, 4}, std::vector<double>(15, 1.0)});
    assert(!value_count.ok());
    assert(has_message_fragment(value_count, "value count"));

    const auto overflow = qmomkit::config::with_initial_moments(dense_case, {{std::numeric_limits<std::size_t>::max(), 2}, std::vector<double>{}});
    assert(!overflow.ok());
    assert(has_message_fragment(overflow, "overflow"));

    qmomkit::RawMomentTensor nonfinite = cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0});
    nonfinite.values[0] = std::numeric_limits<double>::quiet_NaN();
    const auto finite = qmomkit::config::with_initial_moments(dense_case, nonfinite);
    assert(!finite.ok());
    assert(has_message_fragment(finite, "finite"));

    const auto dense_extent = qmomkit::config::with_initial_moments(dense_case, {{3, 4}, std::vector<double>(12, 1.0)});
    assert(!dense_extent.ok());
    assert(has_message_fragment(dense_extent, "2*N"));

    const auto lean_index = qmomkit::config::with_initial_moments(lean_case, {{3, 4}, std::vector<double>(12, 1.0)});
    assert(!lean_index.ok());
    assert(has_message_fragment(lean_index, "required CQMOM moment index"));

    assert(dense_case.initial_state().moment_values() == base_values);
    assert(lean_case.initial_state().moment_values() == lean_base_values);
}

void test_cqmom_case_replacement_defers_realizability()
{
    const qmomkit::CqmomCase base_case = qmomkit::internal::CqmomConstructionAccess::create_case(cqmom_dense_replacement_base());
    qmomkit::RawMomentTensor nonrealizable = cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0});
    nonrealizable.values[0] = -1.0;

    const auto transformed = qmomkit::config::with_initial_moments(base_case, std::move(nonrealizable));
    assert(transformed.ok());
    const qmomkit::CqmomResult reconstruction = qmomkit::compute_cqmom(transformed.value());
    assert(!reconstruction.ok());
    assert(reconstruction.status == qmomkit::CqmomStatus::qmom_failed);

    bool initial_state_failed = false;
    try
    {
        static_cast<void>(transformed.value().initial_state());
    }
    catch (const std::runtime_error& error)
    {
        initial_state_failed = std::string(error.what()).find("CQMOM initial state quadrature reconstruction failed") != std::string::npos;
    }
    assert(initial_state_failed);
}

void test_cqmom_case_file_loading()
{
    const std::filesystem::path filename = std::filesystem::temp_directory_path() / "qmomkit_cqmom_case_configuration_test.mk";
    {
        std::ofstream output(filename);
        output << cqmom_text;
    }
    const auto loaded = qmomkit::config::load_cqmom(filename);
    assert(loaded.ok());
    assert(loaded.value().qmom_options().adaptive);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(loaded.value());
    assert(result.ok());
    std::filesystem::remove(filename);

    const auto missing = qmomkit::config::load_cqmom(filename);
    assert(!missing.ok());
    assert(has_failure(missing));
    assert(has_message_fragment(missing, filename.string()));

    {
        std::ofstream output(filename);
        output << R"(
[phase_space]
node_counts = 2, 2
conditioning_order = 0, 1

[initial]
type = independent_product
number_density = 2.5

[initial.dimension.0]
family = gamma
shape = 3.0
scale = 0.2

[initial.dimension.1]
family = normal
mean = -0.5
standard_deviation = 1.25
)";
    }
    const auto lean_loaded = qmomkit::config::load_cqmom(filename);
    assert(lean_loaded.ok());
    const qmomkit::CqmomResult lean_result = qmomkit::compute_cqmom(lean_loaded.value());
    assert(lean_result.ok());
    std::filesystem::remove(filename);
}

void test_cqmom_independent_product_zero_density_skips_marginal_validation()
{
    const auto configured = qmomkit::config::parse_cqmom_dense(R"(
[phase_space]
node_counts = 1, 1

[initial]
type = independent_product
number_density = 0.0
maximum_orders = 1, 1

[initial.dimension.0]
family = normal
mean = 0.0
standard_deviation = 0.0

[initial.dimension.1]
family = gamma
shape = 0.0
scale = 1.0
)");
    assert(configured.ok());
    assert(configured.value().moments.extents == std::vector<std::size_t>({2, 2}));
    for (const double value : configured.value().moments.values)
    {
        assert_near(value, 0.0);
    }
}

std::string cqmom_adaptive_option_text(const bool adaptive, const bool adaptive_branch_counts)
{
    std::string text = R"(
[phase_space]
node_counts = 1, 1
conditioning_order = 0, 1
adaptive_branch_counts = )";
    text += adaptive_branch_counts ? "true\n" : "false\n";
    text += R"(
[qmom]
adaptive = )";
    text += adaptive ? "true\n" : "false\n";
    text += R"(
[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 2.0
M[1,0] = 3.0
M[1,1] = 6.0
)";
    return text;
}

void assert_cqmom_adaptive_options(const bool adaptive, const bool adaptive_branch_counts)
{
    const auto configured = qmomkit::config::parse_cqmom_dense(cqmom_adaptive_option_text(adaptive, adaptive_branch_counts));
    assert(configured.ok());
    assert(configured.value().qmom_options.adaptive == adaptive);
    assert(configured.value().adaptive_branch_counts == adaptive_branch_counts);
}

void test_cqmom_adaptive_options_validate_compatibility()
{
    assert_cqmom_adaptive_options(false, false);
    assert_cqmom_adaptive_options(false, true);
    assert_cqmom_adaptive_options(true, true);

    const auto incompatible = qmomkit::config::parse_cqmom_dense(cqmom_adaptive_option_text(true, false));
    assert(!incompatible.ok());
    assert(has_failure(incompatible));
    assert(has_message_fragment(incompatible, "qmom.adaptive = true is incompatible with fixed-count CQMOM"));
    assert(has_message_fragment(incompatible, "phase_space.adaptive_branch_counts = true"));

    const auto omitted_qmom_adaptive = qmomkit::config::parse_cqmom_dense(R"(
[phase_space]
node_counts = 1, 1
adaptive_branch_counts = true

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 2.0
M[1,0] = 3.0
M[1,1] = 6.0
)");
    assert(omitted_qmom_adaptive.ok());
    assert(!omitted_qmom_adaptive.value().qmom_options.adaptive);
    assert(omitted_qmom_adaptive.value().adaptive_branch_counts);
}

void test_cqmom_case_parsing_validates_coordinate_realizability_rank()
{
    const auto valid = qmomkit::config::parse_cqmom(cqmom_case_text_with_runtime(R"(

[time_integration.realizability.1]
check_mode = final_only
support = nonnegative
)"));
    assert(valid.ok());

    const auto invalid = qmomkit::config::parse_cqmom(cqmom_case_text_with_runtime(R"(

[time_integration.realizability.2]
check_mode = final_only
support = nonnegative
)"));
    assert(!invalid.ok());
    assert(has_failure(invalid));
    assert(has_message_fragment(invalid, "coordinate index 2"));
    assert(has_message_fragment(invalid, "phase-space rank 2"));
}

void test_adaptive_cqmom_forces_adaptive_wheeler_without_changing_standalone_qmom()
{
    const auto conditional_reduction = qmomkit::config::parse_cqmom_dense(
        cqmom_adaptive_branch_count_text_with_fixed_qmom({0.2, 0.3, 0.5}, {{-1.0, 1.0, 1.0}, {0.0, -0.5, 0.5}}, {0, 1}));
    assert(conditional_reduction.ok());
    assert(!conditional_reduction.value().qmom_options.adaptive);
    assert(conditional_reduction.value().adaptive_branch_counts);

    const qmomkit::CqmomResult conditional = qmomkit::compute_cqmom(conditional_reduction.value());
    assert(conditional.ok());
    assert(conditional.quadrature.adaptive_branch_counts_enabled());
    assert(conditional.quadrature.total_nodes() == 3);
    assert(std::any_of(conditional.diagnostics.begin(), conditional.diagnostics.end(), [](const qmomkit::CqmomDiagnostics& item)
                       { return item.conditioning_level == 1 && item.actual_node_count < item.requested_node_count; }));

    const auto marginal_reduction =
        qmomkit::config::parse_cqmom_dense(cqmom_adaptive_branch_count_text_with_fixed_qmom({0.5, 0.5}, {{0.0, 0.0}, {-1.0, 1.0}}, {0, 1}));
    assert(marginal_reduction.ok());
    assert(!marginal_reduction.value().qmom_options.adaptive);

    const qmomkit::CqmomResult marginal = qmomkit::compute_cqmom(marginal_reduction.value());
    assert(marginal.ok());
    assert(marginal.quadrature.total_nodes() == 2);
    assert(std::any_of(marginal.diagnostics.begin(), marginal.diagnostics.end(), [](const qmomkit::CqmomDiagnostics& item)
                       { return item.conditioning_level == 0 && item.actual_node_count < item.requested_node_count; }));

    const qmomkit::QmomResult fixed_qmom = qmomkit::compute_qmom(qmomkit::MomentSet({1.0, 0.0, 0.0, 0.0}), 2, qmomkit::QmomOptions{false});
    const qmomkit::QmomResult adaptive_qmom = qmomkit::compute_qmom(qmomkit::MomentSet({1.0, 0.0, 0.0, 0.0}), 2, qmomkit::QmomOptions{true});
    assert(!fixed_qmom.ok());
    assert(adaptive_qmom.ok());
    assert(adaptive_qmom.weights.size() == 1);
}

std::string cqmom_diagnostic_level_text(const std::string& diagnostic_level)
{
    std::string text = R"(
[phase_space]
node_counts = 1, 1

)";
    if (!diagnostic_level.empty())
    {
        text += "[cqmom]\ndiagnostic_level = " + diagnostic_level + "\n\n";
    }
    text += R"(
[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.0
M[1,0] = 0.0
M[1,1] = 0.0
)";
    return text;
}

void test_cqmom_diagnostic_level_parses_values()
{
    const auto omitted = qmomkit::config::parse_cqmom_dense(cqmom_diagnostic_level_text(""));
    assert(omitted.ok());
    assert(omitted.value().diagnostic_level == qmomkit::CqmomDiagnosticLevel::none);

    const auto detailed = qmomkit::config::parse_cqmom_dense(cqmom_diagnostic_level_text("detailed"));
    assert(detailed.ok());
    assert(detailed.value().diagnostic_level == qmomkit::CqmomDiagnosticLevel::detailed);

    const auto basic = qmomkit::config::parse_cqmom_dense(cqmom_diagnostic_level_text("basic"));
    assert(basic.ok());
    assert(basic.value().diagnostic_level == qmomkit::CqmomDiagnosticLevel::basic);

    const auto none = qmomkit::config::parse_cqmom_dense(cqmom_diagnostic_level_text("none"));
    assert(none.ok());
    assert(none.value().diagnostic_level == qmomkit::CqmomDiagnosticLevel::none);

    const auto invalid = qmomkit::config::parse_cqmom_dense(cqmom_diagnostic_level_text("verbose"));
    assert(!invalid.ok());
    assert(has_failure(invalid));
    assert(has_message_fragment(invalid, "diagnostic_level"));
    assert(has_message_fragment(invalid, "detailed, basic, or none"));
}

std::string cqmom_vandermonde_solver_text(const std::string& solver)
{
    std::string text = R"(
[phase_space]
node_counts = 1, 1

)";
    if (!solver.empty())
    {
        text += "[vandermonde]\nsolver = " + solver + "\n\n";
    }
    text += R"(
[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.0
M[1,0] = 0.0
M[1,1] = 0.0
)";
    return text;
}

void test_cqmom_vandermonde_solver_parses_values()
{
    const auto omitted = qmomkit::config::parse_cqmom_dense(cqmom_vandermonde_solver_text(""));
    assert(omitted.ok());
    assert(omitted.value().vandermonde_solver == qmomkit::VandermondeSolver::automatic);

    const auto dense = qmomkit::config::parse_cqmom_dense(cqmom_vandermonde_solver_text("dense"));
    assert(dense.ok());
    assert(dense.value().vandermonde_solver == qmomkit::VandermondeSolver::dense);

    const auto pivoted_qr = qmomkit::config::parse_cqmom_dense(cqmom_vandermonde_solver_text("pivoted_qr"));
    assert(pivoted_qr.ok());
    assert(pivoted_qr.value().vandermonde_solver == qmomkit::VandermondeSolver::pivoted_qr);

    const auto svd = qmomkit::config::parse_cqmom_dense(cqmom_vandermonde_solver_text("svd"));
    assert(svd.ok());
    assert(svd.value().vandermonde_solver == qmomkit::VandermondeSolver::svd);

    const auto bp = qmomkit::config::parse_cqmom_lean(cqmom_vandermonde_solver_text("bjorck_pereyra"));
    assert(bp.ok());
    assert(bp.value().vandermonde_solver == qmomkit::VandermondeSolver::bjorck_pereyra);

    const auto automatic = qmomkit::config::parse_cqmom_lean(cqmom_vandermonde_solver_text("automatic"));
    assert(automatic.ok());
    assert(automatic.value().vandermonde_solver == qmomkit::VandermondeSolver::automatic);

    const auto invalid = qmomkit::config::parse_cqmom_dense(cqmom_vandermonde_solver_text("structured"));
    assert(!invalid.ok());
    assert(has_failure(invalid));
    assert(has_message_fragment(invalid, "vandermonde.solver"));
    assert(has_message_fragment(invalid, "svd"));
    assert(has_message_fragment(invalid, "pivoted_qr"));
}

void assert_cqmom_case_reports_vandermonde_solver(const std::string& configured_solver, const qmomkit::VandermondeSolver expected_solver)
{
    const auto parsed = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text(configured_solver));
    assert(parsed.ok());

    const qmomkit::CqmomCase cqmom_case = parsed.value();
    assert(cqmom_case.vandermonde_solver() == expected_solver);

    const qmomkit::CqmomState state = cqmom_case.initial_state();
    assert(state.vandermonde_solver() == expected_solver);
    assert(state.vandermonde_solver() == cqmom_case.vandermonde_solver());
}

void test_cqmom_case_and_state_report_vandermonde_solver()
{
    assert_cqmom_case_reports_vandermonde_solver("", qmomkit::VandermondeSolver::automatic);
    assert_cqmom_case_reports_vandermonde_solver("dense", qmomkit::VandermondeSolver::dense);
    assert_cqmom_case_reports_vandermonde_solver("pivoted_qr", qmomkit::VandermondeSolver::pivoted_qr);
    assert_cqmom_case_reports_vandermonde_solver("svd", qmomkit::VandermondeSolver::svd);
    assert_cqmom_case_reports_vandermonde_solver("bjorck_pereyra", qmomkit::VandermondeSolver::bjorck_pereyra);
}

void test_realizability_policy_defaults_to_disabled()
{
    const auto policy = qmomkit::config::parse_realizability_policy(qmom_text);
    assert(policy.ok());
    assert(policy.value().check_mode == qmomkit::RealizabilityCheckMode::None);

    const auto empty_section = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
)");
    assert(empty_section.ok());
    assert(empty_section.value().check_mode == qmomkit::RealizabilityCheckMode::None);
}

void test_realizability_policy_parses_check_modes()
{
    const auto none = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = none
)");
    assert(none.ok());
    assert(none.value().check_mode == qmomkit::RealizabilityCheckMode::None);

    const auto every_stage = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = every_stage
support = nonnegative
)");
    assert(every_stage.ok());
    assert(every_stage.value().check_mode == qmomkit::RealizabilityCheckMode::EveryStage);
    assert(every_stage.value().support == qmomkit::MomentSupport::nonnegative);

    const auto final_only = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = nonnegative
accept_boundary = false
tolerance = 1e-10
relative_tolerance = false
)");
    assert(final_only.ok());
    assert(final_only.value().check_mode == qmomkit::RealizabilityCheckMode::FinalOnly);
    assert(final_only.value().support == qmomkit::MomentSupport::nonnegative);
    assert(!final_only.value().accept_boundary);
    assert_near(final_only.value().options.tolerance, 1.0e-10);
    assert(!final_only.value().options.relative_tolerance);
}

void test_realizability_policy_parses_support_and_bounded_interval()
{
    const auto real_line = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = real_line
)");
    assert(real_line.ok());
    assert(real_line.value().support == qmomkit::MomentSupport::real_line);

    const auto bounded = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = bounded_interval
lower = -1.0
upper = 2.0
)");
    assert(bounded.ok());
    assert(bounded.value().support == qmomkit::MomentSupport::bounded_interval);
    assert(bounded.value().interval.has_value());
    assert_near(bounded.value().interval->lower, -1.0);
    assert_near(bounded.value().interval->upper, 2.0);
}

void test_coordinate_realizability_policies_parse_sparse_numbered_sections()
{
    const auto real_line = qmomkit::config::parse_coordinate_realizability_policies(R"(
[time_integration.realizability.0]
check_mode = final_only
support = real_line
)");
    assert(real_line.ok());
    assert(real_line.value().size() == 1);
    assert(real_line.value()[0].has_value());
    assert(real_line.value()[0]->support == qmomkit::MomentSupport::real_line);

    const auto nonnegative = qmomkit::config::parse_coordinate_realizability_policies(R"(
[time_integration.realizability.1]
check_mode = every_stage
support = nonnegative
)");
    assert(nonnegative.ok());
    assert(nonnegative.value().size() == 2);
    assert(!nonnegative.value()[0].has_value());
    assert(nonnegative.value()[1].has_value());
    assert(nonnegative.value()[1]->support == qmomkit::MomentSupport::nonnegative);

    const auto bounded_interval = qmomkit::config::parse_coordinate_realizability_policies(R"(
[time_integration.realizability.2]
check_mode = final_only
support = bounded_interval
lower = 0.0
upper = 1.0
)");
    assert(bounded_interval.ok());
    assert(bounded_interval.value().size() == 3);
    assert(bounded_interval.value()[2].has_value());
    assert(bounded_interval.value()[2]->support == qmomkit::MomentSupport::bounded_interval);
    assert(bounded_interval.value()[2]->interval.has_value());
    assert_near(bounded_interval.value()[2]->interval->lower, 0.0);
    assert_near(bounded_interval.value()[2]->interval->upper, 1.0);

    const std::string sparse_text = R"(
[time_integration.realizability.0]
check_mode = final_only
support = real_line

[time_integration.realizability.3]
check_mode = every_stage
support = nonnegative
)";
    const auto sparse = qmomkit::config::parse_coordinate_realizability_policies(sparse_text);
    assert(sparse.ok());
    assert(sparse.value().size() == 4);
    assert(sparse.value()[0].has_value());
    assert(!sparse.value()[1].has_value());
    assert(!sparse.value()[2].has_value());
    assert(sparse.value()[3].has_value());

    const auto malformed = qmomkit::config::parse_coordinate_realizability_policies(R"(
[time_integration.realizability.3]
check_mode = final_only
support = bounded_interval
lower = 0.0
)");
    assert(!malformed.ok());
    assert(has_failure(malformed));

    const auto qmom_case = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1

[moments]
values = 1.0, 0.0

[time_integration.realizability.0]
check_mode = final_only
support = real_line
)");
    assert(qmom_case.ok());
}

void test_quadrature_bounds_do_not_alter_coordinate_realizability_policy()
{
    const std::string first = R"(
[phase_space]
node_counts = 1, 1
support_lower = -2.0, 10.0
support_upper = 2.0, 20.0

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.0
M[1,0] = 0.0
M[1,1] = 0.0

[time_integration.realizability.1]
check_mode = final_only
support = bounded_interval
lower = 0.0
upper = 1.0
tolerance = 1e-10
relative_tolerance = false
accept_boundary = false
)";
    const std::string second = R"(
[phase_space]
node_counts = 1, 1
support_lower = -7.0, -5.0
support_upper = 7.0, 5.0

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.0
M[1,0] = 0.0
M[1,1] = 0.0

[time_integration.realizability.1]
check_mode = final_only
support = bounded_interval
lower = 0.0
upper = 1.0
tolerance = 1e-10
relative_tolerance = false
accept_boundary = false
)";

    assert(qmomkit::config::parse_cqmom(first).ok());
    assert(qmomkit::config::parse_cqmom(second).ok());
    const auto first_policies = qmomkit::config::parse_coordinate_realizability_policies(first);
    const auto second_policies = qmomkit::config::parse_coordinate_realizability_policies(second);
    assert(first_policies.ok());
    assert(second_policies.ok());
    assert(first_policies.value().size() == 2);
    assert(second_policies.value().size() == 2);
    assert(first_policies.value()[1].has_value());
    assert(second_policies.value()[1].has_value());
    assert_realizability_policies_match(*first_policies.value()[1], *second_policies.value()[1]);
}

void test_coordinate_realizability_policy_does_not_alter_quadrature_bounds()
{
    const std::string without_policy = R"(
[phase_space]
node_counts = 1, 1
support_lower = none, 2.0
support_upper = none, 5.0

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.0
M[1,0] = 0.0
M[1,1] = 0.0
)";
    const std::string with_policy = without_policy + R"(

[time_integration.realizability.0]
check_mode = every_stage
support = bounded_interval
lower = -1.0
upper = 1.0
)";

    const auto without_policy_case = qmomkit::config::parse_cqmom(without_policy);
    const auto with_policy_case = qmomkit::config::parse_cqmom(with_policy);
    assert(without_policy_case.ok());
    assert(with_policy_case.ok());
    const qmomkit::CqmomCoordinateSupport& expected = without_policy_case.value().coordinate_support();
    const qmomkit::CqmomCoordinateSupport& actual = with_policy_case.value().coordinate_support();
    assert(actual.size() == expected.size());
    for (std::size_t dimension = 0; dimension < expected.size(); ++dimension)
    {
        assert(actual[dimension].has_value() == expected[dimension].has_value());
        if (expected[dimension].has_value())
        {
            assert_near(actual[dimension]->lower, expected[dimension]->lower);
            assert_near(actual[dimension]->upper, expected[dimension]->upper);
        }
    }
}

void test_realizability_policy_rejects_invalid_values()
{
    const auto unknown_check_mode = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = accepted_step
support = nonnegative
)");
    assert(!unknown_check_mode.ok());
    assert(has_failure(unknown_check_mode));
    assert(has_message_fragment(unknown_check_mode, "check_mode"));

    const auto missing_support = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
)");
    assert(!missing_support.ok());
    assert(has_failure(missing_support));

    const auto unknown_support = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = positive
)");
    assert(!unknown_support.ok());
    assert(has_failure(unknown_support));
    assert(has_message_fragment(unknown_support, "support"));

    const auto invalid_tolerance = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = nonnegative
tolerance = 0.0
)");
    assert(!invalid_tolerance.ok());
    assert(has_failure(invalid_tolerance));

    const auto nonfinite_tolerance = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = nonnegative
tolerance = nan
)");
    assert(!nonfinite_tolerance.ok());
    assert(has_failure(nonfinite_tolerance));

    const auto malformed_bool = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = nonnegative
accept_boundary = yes
)");
    assert(!malformed_bool.ok());
    assert(has_failure(malformed_bool));
}

void test_realizability_policy_rejects_invalid_bounded_interval()
{
    const auto missing_bound = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = bounded_interval
lower = 0.0
)");
    assert(!missing_bound.ok());
    assert(has_failure(missing_bound));

    const auto invalid_interval = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = bounded_interval
lower = 2.0
upper = 1.0
)");
    assert(!invalid_interval.ok());
    assert(has_failure(invalid_interval));

    const auto nonfinite_bound = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = bounded_interval
lower = -inf
upper = 1.0
)");
    assert(!nonfinite_bound.ok());
    assert(has_failure(nonfinite_bound));
}

void test_realizability_policy_rejects_unknown_entries()
{
    const auto unknown_key = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = nonnegative
enabled = true
)");
    assert(!unknown_key.ok());
    assert(has_failure(unknown_key));

    const auto unknown_section = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability.extra]
check_mode = final_only
)");
    assert(!unknown_section.ok());
    assert(has_failure(unknown_section));
}

void test_qmom_case_parsing_tolerates_realizability_policy_without_installing_it()
{
    const auto qmom_case = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1

[moments]
values = 1.0, 0.0

[time_integration.realizability]
check_mode = final_only
support = real_line
)");
    assert(qmom_case.ok());
    assert(qmom_case.value().moments.values() == std::vector<double>({1.0, 0.0}));

    const auto policy = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = real_line
)");
    assert(policy.ok());
    qmomkit::SspRk3Options options;
    options.accept =
        qmomkit::make_realizability_acceptance_predicate(policy.value(), [](const std::vector<double>& state) { return qmomkit::MomentSet(state); });
    assert(options.accept);
    assert(options.accept({1.0, 0.0, 1.0}));
}

void test_case_parsing_validates_realizability_policy_without_owning_it()
{
    const std::string runtime = R"(

[time_integration]
scheme = srcadv

[time_integration.node_advection]
scheme = rk23

[time_integration.node_advection.adaptive]
relative_tolerance = 1.0e-7
absolute_tolerance = 1.0e-10
initial_step = 1.0e-4
max_step = 1.0e-2
min_step = 1.0e-12
max_steps = 200000

[time_integration.realizability]
check_mode = final_only
support = bounded_interval
lower = -1.0
upper = 1.0
tolerance = 1.0e-12
relative_tolerance = false
accept_boundary = false
)";

    const auto qmom_case = qmomkit::config::parse_qmom(qmom_case_text_with_runtime(runtime));
    assert(qmom_case.ok());
    assert(qmom_case.value().node_count == 1);
    assert(!qmom_case.value().options.adaptive);

    const auto cqmom_case = qmomkit::config::parse_cqmom(cqmom_case_text_with_runtime(runtime));
    assert(cqmom_case.ok());
    assert(cqmom_case.value().moment_storage() == qmomkit::CqmomMomentStorage::lean);
    assert(cqmom_case.value().diagnostic_level() == qmomkit::CqmomDiagnosticLevel::none);
    assert(cqmom_case.value().initial_state().moment_storage() == qmomkit::CqmomMomentStorage::lean);

    const auto policy = qmomkit::config::parse_realizability_policy(runtime);
    assert(policy.ok());
    assert(policy.value().check_mode == qmomkit::RealizabilityCheckMode::FinalOnly);
    assert(policy.value().support == qmomkit::MomentSupport::bounded_interval);

    const auto scheme = qmomkit::config::parse_phase_space_evolution_scheme(runtime);
    assert(scheme.ok());
    assert(scheme.value() == qmomkit::PhaseSpaceEvolutionScheme::srcadv);

    const auto node_options = qmomkit::config::parse_node_advection_options(runtime);
    assert(node_options.ok());
    assert(node_options.value().scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert_near(node_options.value().adaptive_options.relative_tolerance, 1.0e-7);
    assert_near(node_options.value().adaptive_options.max_step, 1.0e-2);
}

void test_case_parsing_rejects_invalid_realizability_policy()
{
    const std::vector<std::string> invalid_blocks = {{
        R"(

[time_integration.realizability]
check_mode = accepted_step
support = nonnegative
)",
        R"(

[time_integration.realizability]
check_mode = final_only
support = nonnegative
tolerance = 0.0
)",
        R"(

[time_integration.realizability]
check_mode = final_only
support = bounded_interval
lower = 2.0
upper = 1.0
)",
        R"(

[time_integration.realizability]
check_mode = final_only
)",
        R"(

[time_integration.realizability]
check_mode = final_only
support = nonnegative
enabled = true
)",
    }};

    for (const auto& block : invalid_blocks)
    {
        const auto qmom_case = qmomkit::config::parse_qmom(qmom_case_text_with_runtime(block));
        assert(!qmom_case.ok());
        assert(has_failure(qmom_case));

        const auto cqmom_case = qmomkit::config::parse_cqmom(cqmom_case_text_with_runtime(block));
        assert(!cqmom_case.ok());
        assert(has_failure(cqmom_case));

        const auto policy = qmomkit::config::parse_realizability_policy(block);
        assert(!policy.ok());
        assert(has_failure(policy));
    }
}

void test_phase_space_evolution_scheme_defaults_and_parses_values()
{
    const auto omitted = qmomkit::config::parse_phase_space_evolution_scheme(qmom_text);
    assert(omitted.ok());
    assert(omitted.value() == qmomkit::PhaseSpaceEvolutionScheme::strang);

    const auto strang = qmomkit::config::parse_phase_space_evolution_scheme(R"(
[time_integration]
scheme = strang
)");
    assert(strang.ok());
    assert(strang.value() == qmomkit::PhaseSpaceEvolutionScheme::strang);

    const auto advsrc = qmomkit::config::parse_phase_space_evolution_scheme(R"(
[time_integration]
scheme = advsrc
)");
    assert(advsrc.ok());
    assert(advsrc.value() == qmomkit::PhaseSpaceEvolutionScheme::advsrc);

    const auto srcadv = qmomkit::config::parse_phase_space_evolution_scheme(R"(
[time_integration]
scheme = srcadv
)");
    assert(srcadv.ok());
    assert(srcadv.value() == qmomkit::PhaseSpaceEvolutionScheme::srcadv);
}

void test_phase_space_evolution_scheme_rejects_invalid_values()
{
    const auto invalid = qmomkit::config::parse_phase_space_evolution_scheme(R"(
[time_integration]
scheme = leapfrog
)");
    assert(!invalid.ok());
    assert(has_failure(invalid));
    assert(has_message_fragment(invalid, "strang, advsrc, or srcadv"));

    for (const char* obsolete_scheme : {"cd", "dc"})
    {
        const auto obsolete =
            qmomkit::config::parse_phase_space_evolution_scheme(std::string("[time_integration]\n") + "scheme = " + obsolete_scheme + "\n");
        assert(!obsolete.ok());
        assert(has_failure(obsolete));
        assert(has_message_fragment(obsolete, "strang, advsrc, or srcadv"));
    }

    const auto unknown_time_key = qmomkit::config::parse_phase_space_evolution_scheme(R"(
[time_integration]
stage_count = 3
)");
    assert(!unknown_time_key.ok());
    assert(has_failure(unknown_time_key));
}

void test_node_advection_options_defaults_and_parses_schemes()
{
    qmomkit::NodeAdvectionOptions default_options;
    assert(default_options.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(default_options.adaptive_options.initial_step == qmomkit::default_adaptive_initial_step);
    assert(default_options.adaptive_options.initial_step == 1.0e-6);

    const auto omitted_options = qmomkit::config::parse_node_advection_options(qmom_text);
    assert(omitted_options.ok());
    assert(omitted_options.value().scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert_near(omitted_options.value().adaptive_options.relative_tolerance, default_options.adaptive_options.relative_tolerance);
    assert_near(omitted_options.value().adaptive_options.absolute_tolerance, default_options.adaptive_options.absolute_tolerance);
    assert(omitted_options.value().adaptive_options.initial_step == qmomkit::default_adaptive_initial_step);
    assert_near(omitted_options.value().adaptive_options.max_step, default_options.adaptive_options.max_step);
    assert_near(omitted_options.value().adaptive_options.min_step, default_options.adaptive_options.min_step);
    assert(omitted_options.value().adaptive_options.max_steps == default_options.adaptive_options.max_steps);

    const auto omitted = qmomkit::config::parse_node_advection_options(qmom_text);
    assert(omitted.ok());
    assert(omitted.value().scheme == qmomkit::NodeAdvectionScheme::rk23);

    const auto ssp_rk3 = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection]
scheme = ssp_rk3
)");
    assert(ssp_rk3.ok());
    assert(ssp_rk3.value().scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);

    const auto rk23 = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection]
scheme = rk23
)");
    assert(rk23.ok());
    assert(rk23.value().scheme == qmomkit::NodeAdvectionScheme::rk23);

    const auto rk45 = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection]
scheme = rk45
)");
    assert(rk45.ok());
    assert(rk45.value().scheme == qmomkit::NodeAdvectionScheme::rk45);

    const auto alongside_phase_space_evolution_scheme = qmomkit::config::parse_node_advection_options(R"(
[time_integration]
scheme = advsrc

[time_integration.node_advection]
scheme = rk23
)");
    assert(alongside_phase_space_evolution_scheme.ok());
    assert(alongside_phase_space_evolution_scheme.value().scheme == qmomkit::NodeAdvectionScheme::rk23);
}

void test_node_advection_options_parse_rk23_controls()
{
    const auto parsed = qmomkit::config::parse_node_advection_options(R"(
[time_integration]
scheme = advsrc

[time_integration.node_advection]
scheme = rk23

[time_integration.node_advection.adaptive]
relative_tolerance = 1.0e-8
absolute_tolerance = 1.0e-11
initial_step = 0.02
max_step = 0.04
min_step = 0.001
max_steps = 123
)");
    assert(parsed.ok());
    assert(parsed.value().scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert_near(parsed.value().adaptive_options.relative_tolerance, 1.0e-8);
    assert_near(parsed.value().adaptive_options.absolute_tolerance, 1.0e-11);
    assert_near(parsed.value().adaptive_options.initial_step, 0.02);
    assert_near(parsed.value().adaptive_options.max_step, 0.04);
    assert_near(parsed.value().adaptive_options.min_step, 0.001);
    assert(parsed.value().adaptive_options.max_steps == 123);

    const auto rk45 = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection]
scheme = rk45

[time_integration.node_advection.adaptive]
initial_step = 0.025
)");
    assert(rk45.ok());
    assert(rk45.value().scheme == qmomkit::NodeAdvectionScheme::rk45);
    assert_near(rk45.value().adaptive_options.initial_step, 0.025);

    const auto rk23_without_initial_step = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection]
scheme = rk23
)");
    assert(rk23_without_initial_step.ok());
    assert(rk23_without_initial_step.value().adaptive_options.initial_step == qmomkit::default_adaptive_initial_step);

    const auto rk45_without_initial_step = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection]
scheme = rk45
)");
    assert(rk45_without_initial_step.ok());
    assert(rk45_without_initial_step.value().adaptive_options.initial_step == qmomkit::default_adaptive_initial_step);

    const auto ssp_without_adaptive_step = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection]
scheme = ssp_rk3
)");
    assert(ssp_without_adaptive_step.ok());
    assert(ssp_without_adaptive_step.value().scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);
    assert(ssp_without_adaptive_step.value().adaptive_options.initial_step == qmomkit::default_adaptive_initial_step);

    const auto scheme_only = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection]
scheme = rk23

[time_integration.node_advection.adaptive]
max_step = 0.1
)");
    assert(scheme_only.ok());
    assert(scheme_only.value().scheme == qmomkit::NodeAdvectionScheme::rk23);
}

void test_node_advection_options_reject_invalid_scheme_values()
{
    const auto invalid = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection]
scheme = rk4
)");
    assert(!invalid.ok());
    assert(has_failure(invalid));
    assert(has_message_fragment(invalid, "rk23, rk45, or ssp_rk3"));

    const auto unknown_time_key = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection]
stage_count = 3
)");
    assert(!unknown_time_key.ok());
    assert(has_failure(unknown_time_key));

    const auto obsolete_scheme = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection]
scheme = adaptive_rk23
)");
    assert(!obsolete_scheme.ok());
    assert(has_failure(obsolete_scheme));

    const auto obsolete_subsection = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection.rk23]
initial_step = 0.01
)");
    assert(!obsolete_subsection.ok());
    assert(has_failure(obsolete_subsection));
}

void test_node_advection_options_reject_invalid_rk23_controls()
{
    const auto zero_relative = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection.adaptive]
relative_tolerance = 0.0
)");
    assert(!zero_relative.ok());
    assert(has_failure(zero_relative));
    assert(has_message_fragment(zero_relative, "relative_tolerance"));
    assert(has_message_fragment(zero_relative, "positive finite"));

    const auto nonfinite_absolute = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection.adaptive]
absolute_tolerance = nan
)");
    assert(!nonfinite_absolute.ok());
    assert(has_failure(nonfinite_absolute));
    assert(has_message_fragment(nonfinite_absolute, "absolute_tolerance"));

    const auto zero_initial_step = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection.adaptive]
initial_step = 0.0
)");
    assert(!zero_initial_step.ok());
    assert(has_failure(zero_initial_step));
    assert(has_message_fragment(zero_initial_step, "initial_step"));

    for (const std::string invalid_value : {"-0.1", "nan", "inf", "-inf"})
    {
        const auto invalid_initial_step =
            qmomkit::config::parse_node_advection_options("[time_integration.node_advection.adaptive]\ninitial_step = " + invalid_value + "\n");
        assert(!invalid_initial_step.ok());
        assert(has_failure(invalid_initial_step));
        assert(has_message_fragment(invalid_initial_step, "initial_step"));
    }

    const auto zero_max_steps = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection.adaptive]
max_steps = 0
)");
    assert(!zero_max_steps.ok());
    assert(has_failure(zero_max_steps));
    assert(has_message_fragment(zero_max_steps, "max_steps"));

    const auto min_exceeds_max = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection.adaptive]
min_step = 0.2
max_step = 0.1
)");
    assert(!min_exceeds_max.ok());
    assert(has_failure(min_exceeds_max));
    assert(has_message_fragment(min_exceeds_max, "min_step must not exceed max_step"));

    const auto unknown_rk23_key = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection.adaptive]
growth_limit = 2.0
)");
    assert(!unknown_rk23_key.ok());
    assert(has_failure(unknown_rk23_key));
}

void test_removed_cqmom_source_policy_section_is_unknown()
{
    const std::string removed_key = std::string("[time_integration.cqmom]\nsource_node_count_") + "evolution = preserve_active_" + "counts\n";
    const auto qmom = qmomkit::config::parse_qmom(qmom_case_text_with_runtime(removed_key));
    assert(!qmom.ok());
    assert(has_failure(qmom));

    const auto cqmom = qmomkit::config::parse_cqmom(cqmom_case_text_with_runtime(removed_key));
    assert(!cqmom.ok());
    assert(has_failure(cqmom));

    const auto phase_space_evolution_options = qmomkit::config::parse_cqmom_phase_space_evolution_options(removed_key);
    assert(!phase_space_evolution_options.ok());
    assert(has_failure(phase_space_evolution_options));
}

void test_phase_space_evolution_options_loaders_parse_defaults_explicit_values_and_callbacks()
{
    const auto qmom_defaults = qmomkit::config::parse_qmom_phase_space_evolution_options("");
    assert(qmom_defaults.ok());
    assert(qmom_defaults.value().scheme == qmomkit::PhaseSpaceEvolutionScheme::strang);
    assert(qmom_defaults.value().node_options.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(!qmom_defaults.value().qmom_options.adaptive);
    assert(qmom_defaults.value().realizability_policy.check_mode == qmomkit::RealizabilityCheckMode::None);
    assert(!qmom_defaults.value().source_options.ssp_rk3_options.accept);
    assert(!qmom_defaults.value().source_options.ssp_rk3_options.stage_accept);

    const std::string qmom_runtime = R"(
[qmom]
adaptive = true

[time_integration]
scheme = srcadv

[time_integration.node_advection]
scheme = rk23

[time_integration.node_advection.adaptive]
relative_tolerance = 1.0e-7
max_step = 0.05

[time_integration.realizability]
check_mode = final_only
support = nonnegative
)";
    const auto qmom_options = qmomkit::config::parse_qmom_phase_space_evolution_options(qmom_runtime);
    assert(qmom_options.ok());
    assert(qmom_options.value().scheme == qmomkit::PhaseSpaceEvolutionScheme::srcadv);
    assert(qmom_options.value().node_options.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert_near(qmom_options.value().node_options.adaptive_options.relative_tolerance, 1.0e-7);
    assert_near(qmom_options.value().node_options.adaptive_options.max_step, 0.05);
    assert(qmom_options.value().qmom_options.adaptive);
    assert(qmom_options.value().realizability_policy.check_mode == qmomkit::RealizabilityCheckMode::FinalOnly);
    assert(qmom_options.value().realizability_policy.support == qmomkit::MomentSupport::nonnegative);
    assert(!qmom_options.value().source_options.ssp_rk3_options.accept);
    assert(!qmom_options.value().source_options.ssp_rk3_options.stage_accept);

    const auto cqmom_defaults = qmomkit::config::parse_cqmom_phase_space_evolution_options("");
    assert(cqmom_defaults.ok());
    assert(cqmom_defaults.value().scheme == qmomkit::PhaseSpaceEvolutionScheme::strang);
    assert(cqmom_defaults.value().node_options.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(!cqmom_defaults.value().source_options.ssp_rk3_options.accept);
    assert(!cqmom_defaults.value().source_options.ssp_rk3_options.stage_accept);

    const std::string cqmom_runtime = R"(
[time_integration]
scheme = advsrc

[time_integration.node_advection]
scheme = rk23

[time_integration.node_advection.adaptive]
absolute_tolerance = 1.0e-12
max_step = 0.025

[time_integration.realizability.1]
check_mode = every_stage
support = nonnegative
)";
    const auto cqmom_options = qmomkit::config::parse_cqmom_phase_space_evolution_options(cqmom_runtime);
    assert(cqmom_options.ok());
    assert(cqmom_options.value().scheme == qmomkit::PhaseSpaceEvolutionScheme::advsrc);
    assert(cqmom_options.value().node_options.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert_near(cqmom_options.value().node_options.adaptive_options.absolute_tolerance, 1.0e-12);
    assert_near(cqmom_options.value().node_options.adaptive_options.max_step, 0.025);
    assert(!cqmom_options.value().source_options.ssp_rk3_options.accept);
    assert(!cqmom_options.value().source_options.ssp_rk3_options.stage_accept);
    assert(cqmom_options.value().source_options.coordinate_realizability_policies.size() == 2);
    assert(!cqmom_options.value().source_options.coordinate_realizability_policies[0].has_value());
    assert(cqmom_options.value().source_options.coordinate_realizability_policies[1].has_value());
    assert(cqmom_options.value().source_options.coordinate_realizability_policies[1]->check_mode == qmomkit::RealizabilityCheckMode::EveryStage);
    assert(cqmom_options.value().source_options.coordinate_realizability_policies[1]->support == qmomkit::MomentSupport::nonnegative);

    const std::filesystem::path filename = std::filesystem::temp_directory_path() / "qmomkit_phase_space_evolution_options_configuration_test.qmom";
    {
        std::ofstream output(filename);
        output << cqmom_runtime;
    }
    const auto loaded = qmomkit::config::load_cqmom_phase_space_evolution_options(filename);
    assert(loaded.ok());
    assert(loaded.value().scheme == cqmom_options.value().scheme);
    assert(loaded.value().node_options.scheme == cqmom_options.value().node_options.scheme);
    assert_near(loaded.value().node_options.adaptive_options.max_step, cqmom_options.value().node_options.adaptive_options.max_step);

    {
        std::ofstream output(filename);
        output << qmom_runtime;
    }
    const auto loaded_qmom = qmomkit::config::load_qmom_phase_space_evolution_options(filename);
    assert(loaded_qmom.ok());
    assert(loaded_qmom.value().scheme == qmom_options.value().scheme);
    assert(loaded_qmom.value().node_options.scheme == qmom_options.value().node_options.scheme);
    assert(loaded_qmom.value().qmom_options.adaptive == qmom_options.value().qmom_options.adaptive);
    std::filesystem::remove(filename);
}

void test_phase_space_evolution_options_loaders_reject_invalid_values()
{
    const auto invalid_scheme = qmomkit::config::parse_qmom_phase_space_evolution_options(R"(
[time_integration]
scheme = leapfrog
)");
    assert(!invalid_scheme.ok());
    assert(has_failure(invalid_scheme));

    const auto invalid_qmom = qmomkit::config::parse_qmom_phase_space_evolution_options(R"(
[qmom]
adaptive = sometimes
)");
    assert(!invalid_qmom.ok());
    assert(has_failure(invalid_qmom));
    assert(has_message_fragment(invalid_qmom, "qmom.adaptive"));

    const auto unknown_qmom = qmomkit::config::parse_qmom_phase_space_evolution_options(R"(
[qmom]
diagnostics = true
)");
    assert(!unknown_qmom.ok());
    assert(has_failure(unknown_qmom));

    const auto invalid_node = qmomkit::config::parse_cqmom_phase_space_evolution_options(R"(
[time_integration.node_advection.adaptive]
relative_tolerance = 0.0
)");
    assert(!invalid_node.ok());
    assert(has_failure(invalid_node));

    const auto removed_source_count = qmomkit::config::parse_cqmom_phase_space_evolution_options(
        std::string("[time_integration.cqmom]\nsource_node_count_") + "evolution = freeze\n");
    assert(!removed_source_count.ok());
    assert(has_failure(removed_source_count));
}

void test_cqmom_phase_space_evolution_options_loader_result_drives_phase_space_evolution_path()
{
    const auto case_result = qmomkit::config::parse_cqmom(cqmom_case_text_with_runtime(""));
    assert(case_result.ok());

    const auto configured_options = qmomkit::config::parse_cqmom_phase_space_evolution_options(R"(
[time_integration]
scheme = advsrc

[time_integration.node_advection]
scheme = rk23

[time_integration.node_advection.adaptive]
relative_tolerance = 1.0e-10
absolute_tolerance = 1.0e-12
initial_step = 0.02
max_step = 0.02

)");
    assert(configured_options.ok());

    const auto velocity = [](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.25); };
    const auto source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved =
        qmomkit::advance_cqmom_phase_space(case_result.value().initial_state(), 0.0, 0.1, velocity, source, configured_options.value());
    assert(phase_space_evolved.ok());
    assert(phase_space_evolved.diagnostics.scheme == qmomkit::PhaseSpaceEvolutionScheme::advsrc);
    assert(phase_space_evolved.diagnostics.first_advection_diagnostics.scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert(phase_space_evolved.diagnostics.first_advection_diagnostics.adaptive_status == qmomkit::AdaptiveRkStatus::success);
}

void test_output_options_default_to_disabled()
{
    const auto options = qmomkit::config::parse_output_options(qmom_text);
    assert(options.ok());
    assert(options.value().directory == std::filesystem::path("."));
    assert(!options.value().write_moments);
    assert(!options.value().write_quadrature);
    assert(!options.value().interval.has_value());

    const auto empty_section = qmomkit::config::parse_output_options(R"(
[output]
)");
    assert(empty_section.ok());
    assert(empty_section.value().directory == std::filesystem::path("."));
    assert(!empty_section.value().write_moments);
    assert(!empty_section.value().write_quadrature);
    assert(!empty_section.value().interval.has_value());
}

void test_output_options_parse_enabled_fields()
{
    const auto moments = qmomkit::config::parse_output_options(R"(
[output]
write_moments = true
)");
    assert(moments.ok());
    assert(moments.value().write_moments);
    assert(!moments.value().write_quadrature);

    const auto quadrature = qmomkit::config::parse_output_options(R"(
[output]
write_quadrature = true
)");
    assert(quadrature.ok());
    assert(!quadrature.value().write_moments);
    assert(quadrature.value().write_quadrature);

    const auto all_fields = qmomkit::config::parse_output_options(R"(
[output]
directory = output
write_moments = true
write_quadrature = false
interval = 0.1
)");
    assert(all_fields.ok());
    assert(all_fields.value().directory == std::filesystem::path("output"));
    assert(all_fields.value().write_moments);
    assert(!all_fields.value().write_quadrature);
    assert(all_fields.value().interval.has_value());
    assert_near(*all_fields.value().interval, 0.1);
}

void test_output_options_reject_invalid_values()
{
    const auto zero_interval = qmomkit::config::parse_output_options(R"(
[output]
interval = 0.0
)");
    assert(!zero_interval.ok());
    assert(has_failure(zero_interval));
    assert(has_message_fragment(zero_interval, "interval"));

    const auto nonfinite_interval = qmomkit::config::parse_output_options(R"(
[output]
interval = nan
)");
    assert(!nonfinite_interval.ok());
    assert(has_failure(nonfinite_interval));

    const auto malformed_bool = qmomkit::config::parse_output_options(R"(
[output]
write_moments = yes
)");
    assert(!malformed_bool.ok());
    assert(has_failure(malformed_bool));
    assert(has_message_fragment(malformed_bool, "write_moments"));

    const auto unknown_key = qmomkit::config::parse_output_options(R"(
[output]
write_nodes = true
)");
    assert(!unknown_key.ok());
    assert(has_failure(unknown_key));

    const auto unknown_section = qmomkit::config::parse_output_options(R"(
[output.moments]
M[0] = true
)");
    assert(!unknown_section.ok());
    assert(has_failure(unknown_section));

    const auto empty_directory = qmomkit::config::parse_output_options(R"(
[output]
directory =
)");
    assert(!empty_directory.ok());
    assert(has_failure(empty_directory));
    assert(has_message_fragment(empty_directory, "must not be empty"));

    const auto invalid_syntax = qmomkit::config::parse_output_options(R"(
[output]
directory
)");
    assert(!invalid_syntax.ok());
    assert(has_failure(invalid_syntax));
}

void test_qmom_case_parsing_tolerates_scheme()
{
    const auto qmom = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1

[time_integration]
scheme = srcadv

[time_integration.node_advection]
scheme = ssp_rk3

[time_integration.node_advection.adaptive]
max_step = 0.1

[moments]
values = 1.0, 0.0
)");
    assert(qmom.ok());

    const auto cqmom = qmomkit::config::parse_cqmom_dense(R"(
[phase_space]
node_counts = 1, 1

[time_integration]
scheme = advsrc

[time_integration.node_advection]
scheme = rk23

[time_integration.node_advection.adaptive]
relative_tolerance = 1.0e-8
absolute_tolerance = 1.0e-10

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.0
M[1,0] = 0.0
M[1,1] = 0.0
)");
    assert(cqmom.ok());

    const auto invalid = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1

[time_integration]
scheme = invalid

[moments]
values = 1.0, 0.0
)");
    assert(!invalid.ok());
    assert(has_failure(invalid));

    const auto invalid_node_advection = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1

[time_integration.node_advection]
scheme = rk4

[moments]
values = 1.0, 0.0
)");
    assert(!invalid_node_advection.ok());
    assert(has_failure(invalid_node_advection));
}

void test_case_parsing_rejects_invalid_runtime_phase_space_evolution_settings()
{
    const std::vector<std::string> invalid_shared_blocks = {{
        R"(

[time_integration]
scheme = invalid
)",
        R"(

[time_integration.node_advection]
scheme = rk4
)",
        R"(

[time_integration.node_advection.adaptive]
relative_tolerance = 0.0
)",
        R"(

[time_integration.node_advection.adaptive]
min_step = 2.0
max_step = 1.0
)",
        R"(

[time_integration.node_advection.adaptive]
relative_tolerence = 1.0e-7
)",
    }};

    for (const auto& block : invalid_shared_blocks)
    {
        const auto qmom_case = qmomkit::config::parse_qmom(qmom_case_text_with_runtime(block));
        assert(!qmom_case.ok());
        assert(has_failure(qmom_case));

        const auto cqmom_case = qmomkit::config::parse_cqmom(cqmom_case_text_with_runtime(block));
        assert(!cqmom_case.ok());
        assert(has_failure(cqmom_case));
    }

    const auto invalid_scheme = qmomkit::config::parse_phase_space_evolution_scheme(R"(
[time_integration]
scheme = invalid
)");
    assert(!invalid_scheme.ok());
    assert(has_failure(invalid_scheme));

    const auto invalid_node_options = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection.adaptive]
min_step = 2.0
max_step = 1.0
)");
    assert(!invalid_node_options.ok());
    assert(has_failure(invalid_node_options));

    const auto unknown_node_key = qmomkit::config::parse_node_advection_options(R"(
[time_integration.node_advection.adaptive]
relative_tolerence = 1.0e-7
)");
    assert(!unknown_node_key.ok());
    assert(has_failure(unknown_node_key));

    const std::string source_policy = std::string("\n[time_integration.cqmom]\nsource_node_count_") + "evolution = freeze\n";
    const auto qmom_case = qmomkit::config::parse_qmom(qmom_case_text_with_runtime(source_policy));
    assert(!qmom_case.ok());
    assert(has_failure(qmom_case));

    const auto cqmom_case = qmomkit::config::parse_cqmom(cqmom_case_text_with_runtime(source_policy));
    assert(!cqmom_case.ok());
    assert(has_failure(cqmom_case));
}

void test_qmom_case_parsing_tolerates_output_options_without_installing_them()
{
    const auto qmom = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1

[output]
directory = output
write_moments = true
write_quadrature = true
interval = 0.1

[moments]
values = 1.0, 0.0
)");
    assert(qmom.ok());
    assert(qmom.value().moments.values() == std::vector<double>({1.0, 0.0}));

    const auto cqmom = qmomkit::config::parse_cqmom_dense(R"(
[phase_space]
node_counts = 1, 1

[output]
directory = output
write_moments = true

[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 0.0
M[1,0] = 0.0
M[1,1] = 0.0
)");
    assert(cqmom.ok());

    const auto options = qmomkit::config::parse_output_options(R"(
[output]
directory = output
write_moments = true
write_quadrature = true
interval = 0.1
)");
    assert(options.ok());
    assert(options.value().write_moments);
    assert(options.value().write_quadrature);
}

void test_duplicate_and_missing_entries()
{
    const auto duplicate_key =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[qmom]\nadaptive=false\nadaptive=false\n[moments]\nvalues=1,0\n");
    assert(!duplicate_key.ok());
    assert(has_failure(duplicate_key));

    const auto duplicate_moment = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[moments]\nextents=2,2\nM[0,0]=1\nM[0,0]=1\nM[0,1]=0\nM[1,0]=0\nM[1,1]=0\n");
    assert(!duplicate_moment.ok());
    assert(has_failure(duplicate_moment));

    const auto missing_moment = qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1,1\n[moments]\nextents=2,2\nM[0,0]=1\n");
    assert(!missing_moment.ok());
    assert(has_failure(missing_moment));

    const auto missing_initial_type = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\nfamily=normal\n");
    assert(!missing_initial_type.ok());
    assert(has_failure(missing_initial_type));

    const auto missing_family = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nmaximum_order=1\n");
    assert(!missing_family.ok());
    assert(has_failure(missing_family));

    const auto missing_parameter =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily=normal\nmaximum_order=1\nmean=0\n");
    assert(!missing_parameter.ok());
    assert(has_failure(missing_parameter));

    const auto omitted_maximum_order =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily=gamma\nshape=1\nscale=2\n");
    assert(omitted_maximum_order.ok());
    assert(omitted_maximum_order.value().moments.size() == 2);

    const auto missing_mixture_components = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n");
    assert(!missing_mixture_components.ok());
    assert(has_failure(missing_mixture_components));

    const auto missing_component_family = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                                                      "[initial.mixture.component.0]\nnumber_density=1\n");
    assert(!missing_component_family.ok());
    assert(has_failure(missing_component_family));

    const auto missing_component_density = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                                                       "[initial.mixture.component.0]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!missing_component_density.ok());
    assert(has_failure(missing_component_density));

    const auto missing_cqmom_initial_type = qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1,1\n[initial]\nnumber_density=1\n");
    assert(!missing_cqmom_initial_type.ok());
    assert(has_failure(missing_cqmom_initial_type));

    const auto missing_product_dimensions = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n");
    assert(!missing_product_dimensions.ok());
    assert(has_failure(missing_product_dimensions));

    const auto missing_product_family = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nmean=0\nstandard_deviation=1\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!missing_product_family.ok());
    assert(has_failure(missing_product_family));
}

void test_malformed_nonfinite_invalid_order_and_unknown_key()
{
    const auto malformed = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=two\n[moments]\nvalues=1,0\n", "malformed");
    assert(!malformed.ok());
    assert(has_failure(malformed));
    assert(has_message_fragment(malformed, "malformed:2"));

    const auto nonfinite = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[moments]\nvalues=1,nan\n");
    assert(!nonfinite.ok());
    assert(has_failure(nonfinite));

    const auto malformed_adaptive = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[qmom]\nadaptive=yes\n[moments]\nvalues=1,0\n");
    assert(!malformed_adaptive.ok());
    assert(has_failure(malformed_adaptive));

    const auto invalid_order = qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1,1\nconditioning_order=0,0\n[moments]\nextents=2,2\n"
                                                                  "M[0,0]=1\nM[0,1]=1\nM[1,0]=1\nM[1,1]=1\n");
    assert(!invalid_order.ok());
    assert(has_failure(invalid_order));

    const auto cqmom_initial_section = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[moments]\nextents=2,2\nM[0,0]=1\nM[0,1]=0\nM[1,0]=0\nM[1,1]=0\n[initial]\n");
    assert(cqmom_initial_section.ok());

    const auto unknown = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\nmystery=4\n[moments]\nvalues=1,0\n");
    assert(!unknown.ok());
    assert(has_failure(unknown));

    const auto unsupported_initial_type = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=independent_product\n");
    assert(!unsupported_initial_type.ok());
    assert(has_failure(unsupported_initial_type));

    const auto unsupported_family =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily=poisson\nmaximum_order=1\n");
    assert(!unsupported_family.ok());
    assert(has_failure(unsupported_family));

    const auto invalid_density = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily="
                                                             "normal\nnumber_density=-1\nmaximum_order=1\nmean=0\nstandard_deviation=1\n");
    assert(!invalid_density.ok());
    assert(has_failure(invalid_density));

    const auto nonfinite_density = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily="
                                                               "normal\nnumber_density=nan\nmaximum_order=1\nmean=0\nstandard_deviation=1\n");
    assert(!nonfinite_density.ok());
    assert(has_failure(nonfinite_density));

    const auto invalid_maximum_order = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily="
                                                                   "normal\nmaximum_order=one\nmean=0\nstandard_deviation=1\n");
    assert(!invalid_maximum_order.ok());
    assert(has_failure(invalid_maximum_order));

    const auto invalid_mixture_maximum_order =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=one\n"
                                    "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=1\nshape=1\nscale=2\n");
    assert(!invalid_mixture_maximum_order.ok());
    assert(has_failure(invalid_mixture_maximum_order));

    const auto invalid_distribution_parameter = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily="
                                                                            "normal\nmaximum_order=1\nmean=0\nstandard_deviation=0\n");
    assert(!invalid_distribution_parameter.ok());
    assert(has_failure(invalid_distribution_parameter));

    const auto invalid_component_density =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                    "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=-1\nshape=1\nscale=2\n");
    assert(!invalid_component_density.ok());
    assert(has_failure(invalid_component_density));

    const auto nonfinite_component_density =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                    "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=nan\nshape=1\nscale=2\n");
    assert(!nonfinite_component_density.ok());
    assert(has_failure(nonfinite_component_density));

    const auto invalid_component_parameter =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                    "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=1\nshape=0\nscale=2\n");
    assert(!invalid_component_parameter.ok());
    assert(has_failure(invalid_component_parameter));
    assert(has_message_fragment(invalid_component_parameter, "Component index: 0"));

    const auto invalid_second_component_parameter =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                    "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=1\nshape=1\nscale=2\n"
                                    "[initial.mixture.component.1]\nfamily=normal\nnumber_density=1\nmean=0\nstandard_deviation=0\n");
    assert(!invalid_second_component_parameter.ok());
    assert(has_failure(invalid_second_component_parameter));
    assert(has_message_fragment(invalid_second_component_parameter, "Component index: 1"));

    const auto all_zero_mixture = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                                              "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=0\nshape=1\nscale=2\n");
    assert(!all_zero_mixture.ok());
    assert(has_failure(all_zero_mixture));

    const auto noncontiguous_mixture =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                    "[initial.mixture.component.1]\nfamily=gamma\nnumber_density=1\nshape=1\nscale=2\n");
    assert(!noncontiguous_mixture.ok());
    assert(has_failure(noncontiguous_mixture));

    const auto unknown_component_key =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                    "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=1\nshape=1\nscale=2\nweight_fraction=0.5\n");
    assert(!unknown_component_key.ok());
    assert(has_failure(unknown_component_key));

    const auto component_with_distribution_type =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily=normal\nmaximum_order=1\nmean=0\n"
                                    "standard_deviation=1\n[initial.mixture.component.0]\nfamily=gamma\nnumber_density=1\nshape=1\nscale=2\n");
    assert(!component_with_distribution_type.ok());
    assert(has_failure(component_with_distribution_type));

    const auto unsupported_cqmom_initial_type =
        qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1,1\n[initial]\ntype=distribution\nmaximum_order=1\nfamily=normal\nmean=0\n"
                                           "standard_deviation=1\n");
    assert(!unsupported_cqmom_initial_type.ok());
    assert(has_failure(unsupported_cqmom_initial_type));

    const auto invalid_product_order = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=one\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!invalid_product_order.ok());
    assert(has_failure(invalid_product_order));

    const auto product_tensor_size_overflow =
        qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\n"
                                           "maximum_orders=4294967295,4294967295\n"
                                           "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
                                           "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!product_tensor_size_overflow.ok());
    assert(has_failure(product_tensor_size_overflow));

    const auto invalid_product_density = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=-1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!invalid_product_density.ok());
    assert(has_failure(invalid_product_density));

    const auto noncontiguous_product_dimensions = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!noncontiguous_product_dimensions.ok());
    assert(has_failure(noncontiguous_product_dimensions));

    const auto extra_product_dimension = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
        "[initial.dimension.1]\nfamily=normal\nmean=0\nstandard_deviation=1\n"
        "[initial.dimension.2]\nfamily=beta\nalpha=1\nbeta=2\n");
    assert(!extra_product_dimension.ok());
    assert(has_failure(extra_product_dimension));

    const auto duplicate_product_dimension = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!duplicate_product_dimension.ok());
    assert(has_failure(duplicate_product_dimension));

    const auto invalid_marginal_parameter = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=0\nscale=2\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!invalid_marginal_parameter.ok());
    assert(has_failure(invalid_marginal_parameter));

    const auto product_dimension_mismatch =
        qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1\n"
                                           "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!product_dimension_mismatch.ok());
    assert(has_failure(product_dimension_mismatch));

    const auto product_moment_extent_too_small = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=2,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!product_moment_extent_too_small.ok());
    assert(has_failure(product_moment_extent_too_small));

    const auto product_rejects_marginal_density = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nfamily=gamma\nnumber_density=2\nshape=1\nscale=2\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!product_rejects_marginal_density.ok());
    assert(has_failure(product_rejects_marginal_density));
}

void test_explicit_moments_take_precedence_over_initial_distribution()
{
    const auto configured = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 2

[moments]
values = 3.0, 1.0, 2.0, 4.0

[initial]
type = distribution
family = normal
number_density = -1.0
maximum_order = 3
mean = 0.0
standard_deviation = 0.0
)");
    assert(configured.ok());
    assert(configured.value().moments.values() == std::vector<double>({3.0, 1.0, 2.0, 4.0}));

    const auto initial = qmomkit::config::parse_qmom_initial_condition(R"(
[phase_space]
node_counts = 2

[moments]
values = 3.0, 1.0, 2.0, 4.0

[initial]
type = mixture
)");
    assert(initial.ok());
    const auto moments = qmomkit::config::initialize_moments(initial.value());
    assert(moments.ok());
    assert(moments.value().values() == std::vector<double>({3.0, 1.0, 2.0, 4.0}));

    const auto mixture = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 2

[moments]
values = 3.0, 1.0, 2.0, 4.0

[initial]
type = mixture
maximum_order = 3

[initial.mixture.component.0]
family = gamma
number_density = -1.0
shape = 0.0
scale = 0.0
)");
    assert(mixture.ok());
    assert(mixture.value().moments.values() == std::vector<double>({3.0, 1.0, 2.0, 4.0}));

    const auto unknown_initial_key = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 2

[moments]
values = 3.0, 1.0, 2.0, 4.0

[initial]
type = mixture
maximum_order = 3
weight_fraction = 0.5
)");
    assert(unknown_initial_key.ok());
    assert(unknown_initial_key.value().moments.values() == std::vector<double>({3.0, 1.0, 2.0, 4.0}));

    const auto malformed_ignored_component_section = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 2

[moments]
values = 3.0, 1.0, 2.0, 4.0

[initial.mixture.component.alpha]
family = gamma
number_density = -1.0
weight_fraction = 0.5
)");
    assert(malformed_ignored_component_section.ok());
    assert(malformed_ignored_component_section.value().moments.values() == std::vector<double>({3.0, 1.0, 2.0, 4.0}));

    const auto unrelated_unknown_section = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 2

[moments]
values = 3.0, 1.0, 2.0, 4.0

[unrelated]
value = 1
)");
    assert(!unrelated_unknown_section.ok());
    assert(has_failure(unrelated_unknown_section));

    const auto cqmom = qmomkit::config::parse_cqmom_dense(R"(
[phase_space]
node_counts = 1, 1

[moments]
extents = 2, 2
M[0,0] = 3.0
M[0,1] = 1.0
M[1,0] = 2.0
M[1,1] = 4.0

[initial]
type = independent_product
number_density = -1.0
maximum_orders = 1, 1

[initial.dimension.0]
family = gamma
shape = 0.0
scale = 0.0
)");
    assert(cqmom.ok());
    assert(cqmom.value().moments.values == std::vector<double>({3.0, 1.0, 2.0, 4.0}));
}

void test_configuration_format_infers_formulation_from_phase_space()
{
    const auto qmom = qmomkit::config::parse_qmom(R"(
[phase_space]
names = d
node_counts = 3

[moments]
values = 1.0, 0.5, 0.4, 0.35, 0.32, 0.3
)");
    assert(qmom.ok());
    assert(qmom.value().node_count == 3);

    const auto cqmom = qmomkit::config::parse_cqmom_lean(R"(
[phase_space]
names = d, u, T
node_counts = 1, 1, 1

[initial]
type = independent_product
number_density = 2.0

[initial.dimension.0]
family = gamma
shape = 2.0
scale = 0.5

[initial.dimension.1]
family = normal
mean = 0.0
standard_deviation = 1.0

[initial.dimension.2]
family = beta
alpha = 2.0
beta = 3.0
)");
    assert(cqmom.ok());
    assert(cqmom.value().moments.layout.hierarchy.node_counts() == std::vector<std::size_t>({1, 1, 1}));
    assert(cqmom.value().moments.layout.hierarchy.conditioning_order() == std::vector<std::size_t>({0, 1, 2}));

    const auto explicit_order = qmomkit::config::parse_cqmom_lean(R"(
[phase_space]
names = d, u, T
node_counts = 1, 1, 1
conditioning_order = 2, 0, 1

[initial]
type = independent_product
number_density = 2.0

[initial.dimension.0]
family = gamma
shape = 2.0
scale = 0.5

[initial.dimension.1]
family = normal
mean = 0.0
standard_deviation = 1.0

[initial.dimension.2]
family = beta
alpha = 2.0
beta = 3.0
)");
    assert(explicit_order.ok());
    assert(explicit_order.value().moments.layout.hierarchy.conditioning_order() == std::vector<std::size_t>({2, 0, 1}));
}

void test_configuration_rejects_unsupported_or_ambiguous_format_fields()
{
    const auto unsupported_version = qmomkit::config::parse_qmom("format_version=2\n[phase_space]\nnode_counts=1\n[moments]\nvalues=1,0\n");
    assert(!unsupported_version.ok());
    assert(has_failure(unsupported_version));
    assert(has_message_fragment(unsupported_version, "format_version must be 1"));

    const auto qmom_with_multivariate_phase_space = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1,1\n[moments]\nvalues=1,0\n");
    assert(!qmom_with_multivariate_phase_space.ok());
    assert(has_failure(qmom_with_multivariate_phase_space));
    assert(has_message_fragment(qmom_with_multivariate_phase_space, "exactly one positive integer"));

    const auto cqmom_with_univariate_phase_space =
        qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1\n[moments]\nextents=2\nM[0]=1\nM[1]=0\n");
    assert(!cqmom_with_univariate_phase_space.ok());
    assert(has_failure(cqmom_with_univariate_phase_space));
    assert(has_message_fragment(cqmom_with_univariate_phase_space, "one-dimensional phase space is interpreted as QMOM"));
}

void test_file_loading_and_file_failure()
{
    const std::filesystem::path filename = std::filesystem::temp_directory_path() / "qmomkit_configuration_test.mk";
    {
        std::ofstream output(filename);
        output << qmom_text;
    }
    const auto loaded = qmomkit::config::load_qmom(filename);
    assert(loaded.ok());
    std::filesystem::remove(filename);

    const auto missing = qmomkit::config::load_qmom(filename);
    assert(!missing.ok());
    assert(has_failure(missing));
    assert(has_message_fragment(missing, filename.string()));

    {
        std::ofstream output(filename);
        output << cqmom_text;
    }
    const auto loaded_cqmom = qmomkit::config::load_cqmom_lean(filename);
    assert(loaded_cqmom.ok());
    assert(loaded_cqmom.value().moments.layout.hierarchy.conditioning_order() == std::vector<std::size_t>({1, 0}));
    assert_near(loaded_cqmom.value().moments.moment({0, 0}), 1.0);
    std::filesystem::remove(filename);

    const auto missing_cqmom = qmomkit::config::load_cqmom_lean(filename);
    assert(!missing_cqmom.ok());
    assert(has_failure(missing_cqmom));
    assert(has_message_fragment(missing_cqmom, filename.string()));

    {
        std::ofstream output(filename);
        output << "[output]\ndirectory = output\nwrite_moments = true\ninterval = 0.25\n";
    }
    const auto output_options = qmomkit::config::load_output_options(filename);
    assert(output_options.ok());
    assert(output_options.value().directory == std::filesystem::path("output"));
    assert(output_options.value().write_moments);
    assert(output_options.value().interval.has_value());
    assert_near(*output_options.value().interval, 0.25);
    std::filesystem::remove(filename);

    const auto missing_output = qmomkit::config::load_output_options(filename);
    assert(!missing_output.ok());
    assert(has_failure(missing_output));

    {
        std::ofstream output(filename);
        output << "[time_integration.node_advection]\nscheme = rk23\n"
                  "[time_integration.node_advection.adaptive]\nmax_step = 0.05\n";
    }
    const auto node_advection_options = qmomkit::config::load_node_advection_options(filename);
    assert(node_advection_options.ok());
    assert(node_advection_options.value().scheme == qmomkit::NodeAdvectionScheme::rk23);
    assert_near(node_advection_options.value().adaptive_options.max_step, 0.05);

    const auto node_advection_scheme = qmomkit::config::load_node_advection_options(filename);
    assert(node_advection_scheme.ok());
    assert(node_advection_scheme.value().scheme == qmomkit::NodeAdvectionScheme::rk23);

    const auto phase_space_evolution_options = qmomkit::config::load_cqmom_phase_space_evolution_options(filename);
    assert(phase_space_evolution_options.ok());
    assert(phase_space_evolution_options.value().node_options.scheme == qmomkit::NodeAdvectionScheme::rk23);
    std::filesystem::remove(filename);

    const auto missing_node_advection = qmomkit::config::load_node_advection_options(filename);
    assert(!missing_node_advection.ok());
    assert(has_failure(missing_node_advection));

    const auto missing_phase_space_evolution_options = qmomkit::config::load_cqmom_phase_space_evolution_options(filename);
    assert(!missing_phase_space_evolution_options.ok());
    assert(has_failure(missing_phase_space_evolution_options));

    {
        std::ofstream output(filename);
        output << "[time_integration.realizability]\ncheck_mode = final_only\nsupport = nonnegative\n";
    }
    const auto loaded_policy = qmomkit::config::load_realizability_policy(filename);
    assert(loaded_policy.ok());
    assert(loaded_policy.value().check_mode == qmomkit::RealizabilityCheckMode::FinalOnly);
    assert(loaded_policy.value().support == qmomkit::MomentSupport::nonnegative);
    std::filesystem::remove(filename);

    const auto missing_policy = qmomkit::config::load_realizability_policy(filename);
    assert(!missing_policy.ok());
    assert(has_failure(missing_policy));

    {
        std::ofstream output(filename);
        output << qmom_case_text_with_runtime(R"(

[time_integration.realizability]
check_mode = final_only
)");
    }
    const auto invalid_qmom_runtime = qmomkit::config::load_qmom(filename);
    assert(!invalid_qmom_runtime.ok());
    assert(has_failure(invalid_qmom_runtime));
    std::filesystem::remove(filename);

    {
        std::ofstream output(filename);
        output << cqmom_case_text_with_runtime(std::string("\n[time_integration.cqmom]\nsource_node_count_") + "evolution = freeze\n");
    }
    const auto invalid_cqmom_runtime = qmomkit::config::load_cqmom(filename);
    assert(!invalid_cqmom_runtime.ok());
    assert(has_failure(invalid_cqmom_runtime));
    std::filesystem::remove(filename);
}

} // namespace

void run_all_configuration_tests()
{
    test_configuration_result_api();
    test_valid_qmom_matches_programmatic_case();
    test_phase_space_automatic_scaling_parses_for_qmom_and_cqmom();
    test_qmom_case_parses_abscissa_bounds();
    test_qmom_distribution_initialization_matches_direct_api();
    test_all_supported_distribution_families_parse();
    test_qmom_mixture_initialization_matches_direct_api();
    test_mixture_zero_density_component_is_skipped();
    test_valid_cqmom_matches_programmatic_case();
    test_cqmom_independent_product_initialization_matches_direct_api();
    test_cqmom_independent_product_mixture_marginal_initialization();
    test_lean_cqmom_initialization_matches_dense_subselection();
    test_lean_cqmom_initialization_matches_example_flow();
    test_lean_cqmom_initialization_omits_maximum_orders();
    test_dense_cqmom_initialization_omits_maximum_orders();
    test_lean_cqmom_initialization_selects_explicit_dense_tensor();
    test_cqmom_lean_case_uses_selected_analytical_moments();
    test_cqmom_lean_case_matches_dense_case_path();
    test_cqmom_lean_case_selects_explicit_dense_tensor();
    test_cqmom_case_uses_selected_analytical_moments();
    test_cqmom_case_accepts_explicit_dense_tensor();
    test_cqmom_case_parses_coordinate_support();
    test_cqmom_factories_preserve_coordinate_support();
    test_cqmom_case_rejects_invalid_coordinate_support();
    test_cqmom_case_moment_storage_selection();
    test_cqmom_case_initial_state_preserves_case_settings();
    test_cqmom_case_replaces_dense_initial_moments_immutably();
    test_cqmom_case_replaces_lean_initial_moments_immutably();
    test_fixed_cqmom_prepared_metadata_ownership();
    test_fixed_cqmom_prepared_metadata_is_concurrently_readable();
    test_fixed_cqmom_state_prepared_metadata_invariant_and_compatibility();
    test_cqmom_case_replacement_rejects_structural_tensor_errors();
    test_cqmom_case_replacement_defers_realizability();
    test_cqmom_case_file_loading();
    test_cqmom_independent_product_zero_density_skips_marginal_validation();
    test_cqmom_adaptive_options_validate_compatibility();
    test_adaptive_cqmom_forces_adaptive_wheeler_without_changing_standalone_qmom();
    test_cqmom_diagnostic_level_parses_values();
    test_cqmom_vandermonde_solver_parses_values();
    test_cqmom_case_and_state_report_vandermonde_solver();
    test_realizability_policy_defaults_to_disabled();
    test_realizability_policy_parses_check_modes();
    test_realizability_policy_parses_support_and_bounded_interval();
    test_coordinate_realizability_policies_parse_sparse_numbered_sections();
    test_quadrature_bounds_do_not_alter_coordinate_realizability_policy();
    test_coordinate_realizability_policy_does_not_alter_quadrature_bounds();
    test_cqmom_case_parsing_validates_coordinate_realizability_rank();
    test_realizability_policy_rejects_invalid_values();
    test_realizability_policy_rejects_invalid_bounded_interval();
    test_realizability_policy_rejects_unknown_entries();
    test_qmom_case_parsing_tolerates_realizability_policy_without_installing_it();
    test_case_parsing_validates_realizability_policy_without_owning_it();
    test_case_parsing_rejects_invalid_realizability_policy();
    test_phase_space_evolution_scheme_defaults_and_parses_values();
    test_phase_space_evolution_scheme_rejects_invalid_values();
    test_node_advection_options_defaults_and_parses_schemes();
    test_node_advection_options_parse_rk23_controls();
    test_node_advection_options_reject_invalid_scheme_values();
    test_node_advection_options_reject_invalid_rk23_controls();
    test_removed_cqmom_source_policy_section_is_unknown();
    test_phase_space_evolution_options_loaders_parse_defaults_explicit_values_and_callbacks();
    test_qmom_automatic_scaling_handoff_to_evolution_options();
    test_phase_space_evolution_options_loaders_reject_invalid_values();
    test_cqmom_phase_space_evolution_options_loader_result_drives_phase_space_evolution_path();
    test_output_options_default_to_disabled();
    test_output_options_parse_enabled_fields();
    test_output_options_reject_invalid_values();
    test_qmom_case_parsing_tolerates_scheme();
    test_case_parsing_rejects_invalid_runtime_phase_space_evolution_settings();
    test_qmom_case_parsing_tolerates_output_options_without_installing_them();
    test_duplicate_and_missing_entries();
    test_malformed_nonfinite_invalid_order_and_unknown_key();
    test_explicit_moments_take_precedence_over_initial_distribution();
    test_configuration_format_infers_formulation_from_phase_space();
    test_configuration_rejects_unsupported_or_ambiguous_format_fields();
    test_file_loading_and_file_failure();
}

void run_initialization_tests()
{
    test_configuration_result_api();
    test_valid_qmom_matches_programmatic_case();
    test_qmom_case_parses_abscissa_bounds();
    test_qmom_distribution_initialization_matches_direct_api();
    test_all_supported_distribution_families_parse();
    test_qmom_mixture_initialization_matches_direct_api();
    test_mixture_zero_density_component_is_skipped();
    test_valid_cqmom_matches_programmatic_case();
    test_cqmom_independent_product_initialization_matches_direct_api();
    test_cqmom_independent_product_mixture_marginal_initialization();
    test_lean_cqmom_initialization_matches_dense_subselection();
    test_lean_cqmom_initialization_matches_example_flow();
    test_lean_cqmom_initialization_omits_maximum_orders();
    test_dense_cqmom_initialization_omits_maximum_orders();
    test_lean_cqmom_initialization_selects_explicit_dense_tensor();
}

void run_cqmom_case_state_tests()
{
    test_cqmom_lean_case_uses_selected_analytical_moments();
    test_cqmom_lean_case_matches_dense_case_path();
    test_cqmom_lean_case_selects_explicit_dense_tensor();
    test_cqmom_case_uses_selected_analytical_moments();
    test_cqmom_case_accepts_explicit_dense_tensor();
    test_cqmom_case_parses_coordinate_support();
    test_cqmom_factories_preserve_coordinate_support();
    test_cqmom_case_rejects_invalid_coordinate_support();
    test_cqmom_case_moment_storage_selection();
    test_cqmom_case_initial_state_preserves_case_settings();
    test_cqmom_case_replaces_dense_initial_moments_immutably();
    test_cqmom_case_replaces_lean_initial_moments_immutably();
    test_fixed_cqmom_prepared_metadata_ownership();
    test_fixed_cqmom_prepared_metadata_is_concurrently_readable();
    test_fixed_cqmom_state_prepared_metadata_invariant_and_compatibility();
    test_cqmom_case_replacement_rejects_structural_tensor_errors();
    test_cqmom_case_replacement_defers_realizability();
    test_cqmom_case_file_loading();
    test_cqmom_independent_product_zero_density_skips_marginal_validation();
}

void run_runtime_option_tests()
{
    test_cqmom_adaptive_options_validate_compatibility();
    test_adaptive_cqmom_forces_adaptive_wheeler_without_changing_standalone_qmom();
    test_cqmom_diagnostic_level_parses_values();
    test_cqmom_vandermonde_solver_parses_values();
    test_cqmom_case_and_state_report_vandermonde_solver();
    test_realizability_policy_defaults_to_disabled();
    test_realizability_policy_parses_check_modes();
    test_realizability_policy_parses_support_and_bounded_interval();
    test_coordinate_realizability_policies_parse_sparse_numbered_sections();
    test_quadrature_bounds_do_not_alter_coordinate_realizability_policy();
    test_coordinate_realizability_policy_does_not_alter_quadrature_bounds();
    test_cqmom_case_parsing_validates_coordinate_realizability_rank();
    test_realizability_policy_rejects_invalid_values();
    test_realizability_policy_rejects_invalid_bounded_interval();
    test_realizability_policy_rejects_unknown_entries();
    test_qmom_case_parsing_tolerates_realizability_policy_without_installing_it();
    test_case_parsing_validates_realizability_policy_without_owning_it();
    test_case_parsing_rejects_invalid_realizability_policy();
    test_phase_space_evolution_scheme_defaults_and_parses_values();
    test_phase_space_evolution_scheme_rejects_invalid_values();
    test_node_advection_options_defaults_and_parses_schemes();
    test_node_advection_options_parse_rk23_controls();
    test_node_advection_options_reject_invalid_scheme_values();
    test_node_advection_options_reject_invalid_rk23_controls();
    test_removed_cqmom_source_policy_section_is_unknown();
    test_phase_space_evolution_options_loaders_parse_defaults_explicit_values_and_callbacks();
    test_qmom_automatic_scaling_handoff_to_evolution_options();
    test_phase_space_evolution_options_loaders_reject_invalid_values();
    test_cqmom_phase_space_evolution_options_loader_result_drives_phase_space_evolution_path();
    test_output_options_default_to_disabled();
    test_output_options_parse_enabled_fields();
    test_output_options_reject_invalid_values();
    test_qmom_case_parsing_tolerates_scheme();
    test_case_parsing_rejects_invalid_runtime_phase_space_evolution_settings();
    test_qmom_case_parsing_tolerates_output_options_without_installing_them();
}

void run_failure_path_tests()
{
    test_duplicate_and_missing_entries();
    test_malformed_nonfinite_invalid_order_and_unknown_key();
    test_explicit_moments_take_precedence_over_initial_distribution();
    test_configuration_format_infers_formulation_from_phase_space();
    test_configuration_rejects_unsupported_or_ambiguous_format_fields();
    test_file_loading_and_file_failure();
}

int main(const int argc, char** argv)
{
    if (argc > 2)
    {
        std::cerr << "Usage: test_configuration [initialization|cqmom_case_state|runtime_options|failure_paths]\n";
        return 2;
    }

    const std::string selector = argc == 2 ? argv[1] : "";
    if (selector.empty())
    {
        run_all_configuration_tests();
        return 0;
    }
    if (selector == "initialization")
    {
        run_initialization_tests();
        return 0;
    }
    if (selector == "cqmom_case_state")
    {
        run_cqmom_case_state_tests();
        return 0;
    }
    if (selector == "runtime_options")
    {
        run_runtime_option_tests();
        return 0;
    }
    if (selector == "failure_paths")
    {
        run_failure_path_tests();
        return 0;
    }

    std::cerr << "Unknown configuration selector: " << selector << '\n';
    return 2;
}
