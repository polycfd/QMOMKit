#include "qmomkit/configuration.hpp"

#include "internal/cqmom_case.hpp"
#include "internal/recursive_quadrature.hpp"
#include "qmomkit/distributions.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "test_case_dispatch.hpp"

#include <algorithm>
#include <cassert>
#include <clocale>
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

std::string cqmom_adaptive_text(const std::vector<double>& weights, const std::vector<std::vector<double>>& abscissas,
                                const std::vector<std::size_t>& conditioning_order)
{
    const std::vector<double> moments = raw_moments_from_terminals_2d(weights, abscissas);
    std::ostringstream text;
    text.precision(17);
    text << "[phase_space]\n"
         << "node_counts = 2, 2\n"
         << "conditioning_order = " << conditioning_order[0] << ", " << conditioning_order[1] << "\n\n"
         << "[qmom]\n"
         << "adaptive = true\n\n"
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

// Compare CQMOM moments against an analytically computed full tensor. Only the
// moments the hierarchy reconstructs from are held by a case, so the expected
// tensor is reduced the same way before comparison; entries outside that set
// are not part of the case and are deliberately not checked here.
void assert_cqmom_moments_near_tensor(const qmomkit::CqmomMomentData& actual, const qmomkit::RawMomentTensor& expected_tensor)
{
    assert_cqmom_moments_near(actual, qmomkit::extract_cqmom_moments(expected_tensor, actual.layout));
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

void assert_coordinate_support_matches_shifted_interval(const qmomkit::CqmomAbscissaBounds& support)
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
    const qmomkit::QmomCase direct_qmom_case{qmomkit::MomentSet({1.0, 0.5, 0.5, 0.5}), 2, {}};
    assert(configured.value().node_count == direct_qmom_case.node_count);
    assert(configured.value().moments.values() == direct_qmom_case.moments.values());
    assert(!configured.value().qmom_options.adaptive);

    const auto configured_result = qmomkit::compute_qmom(configured.value());
    const auto direct_result = qmomkit::compute_qmom(direct_qmom_case);
    assert(configured_result.ok());
    assert(direct_result.ok());
    assert(configured_result.weights == direct_result.weights);
    assert(configured_result.abscissas == direct_result.abscissas);

    const auto adaptive = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[qmom]\nadaptive=true\n[moments]\nvalues=1,0\n");
    assert(adaptive.ok());
    assert(adaptive.value().qmom_options.adaptive);

    const std::size_t impossible_node_count = std::numeric_limits<std::size_t>::max() / 2 + 1;
    const auto impossible =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=" + std::to_string(impossible_node_count) + "\n[moments]\nvalues=1,0\n");
    assert(!impossible.ok());
    assert(has_message_fragment(impossible, "too large"));
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
    const auto moments = qmomkit::config::initialize_qmom_moments(initial.value());
    assert(moments.ok());
    const qmomkit::RawMomentResult expected = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{2.0, 3.0}, 3, 4.0);
    assert(expected.ok());
    assert_moments_near(moments.value(), expected.moments);

    const auto qmom_case_result = qmomkit::config::create_qmom_case(initial.value(), 2);
    assert(qmom_case_result.ok());
    assert(qmom_case_result.value().moments.size() == 4);
    assert_near(qmom_case_result.value().moments[0], 4.0);
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

// A mixture takes each distribution and number density from its components,
// so the same keys at the mixture's own level would have no effect. They are
// reported for a QMOM mixture and for a CQMOM mixture marginal.
void test_mixture_level_distribution_keys_are_reported()
{
    const std::string components = R"(
[initial.mixture.component.0]
family = gamma
number_density = 1.5
shape = 2.0
scale = 0.4
)";
    for (const std::string& key : {std::string("number_density = 10.0"), std::string("mean = 3.0")})
    {
        const std::string text = "[phase_space]\nnode_counts = 1\n[initial]\ntype = mixture\n" + key + "\n" + components;
        const std::string message = "Key 'initial." + key.substr(0, key.find(' ')) + "' has no effect with 'type = mixture'";
        assert(has_message_fragment(qmomkit::config::parse_qmom(text), message));
        assert(has_message_fragment(qmomkit::config::parse_qmom_initial_condition(text), message));
    }
    assert(qmomkit::config::parse_qmom("[phase_space]\nnode_counts = 1\n[initial]\ntype = mixture\n" + components).ok());

    const std::string cqmom_text = R"(
[phase_space]
node_counts = 1, 1

[initial]
type = independent_product
number_density = 2.75

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
)";
    assert(qmomkit::config::parse_cqmom_initial_condition(cqmom_text).ok());
    std::string with_family = cqmom_text;
    with_family.insert(with_family.find("type = mixture\n") + std::string("type = mixture\n").size(), "family = normal\n");
    assert(has_message_fragment(qmomkit::config::parse_cqmom_initial_condition(with_family),
                                "Key 'initial.dimension.1.family' has no effect with 'type = mixture'"));
}

// cqmom.moment_storage selected between two moment representations. Only one
// remains, so the key is no longer recognized and is rejected like any other
// unknown key, through every CQMOM loader.
void test_cqmom_rejects_the_removed_moment_storage_key()
{
    const std::string with_removed_key = std::string(cqmom_text) + "\n[cqmom]\nmoment_storage = lean\n";

    const auto configured = qmomkit::config::parse_cqmom_prepared(with_removed_key, "cqmom-test");
    assert(!configured.ok());
    assert(has_failure(configured));
    assert(has_message_fragment(configured, "moment_storage"));

    const auto cqmom_case_result = qmomkit::config::parse_cqmom(with_removed_key, "cqmom-test");
    assert(!cqmom_case_result.ok());
    assert(has_failure(cqmom_case_result));

    // The same configuration without the key still loads through both loaders.
    assert(qmomkit::config::parse_cqmom_prepared(cqmom_text, "cqmom-test").ok());
    assert(qmomkit::config::parse_cqmom(cqmom_text, "cqmom-test").ok());
}

void test_valid_cqmom_matches_programmatic_case()
{
    const auto configured = qmomkit::config::parse_cqmom(cqmom_text, "cqmom-test");
    assert(configured.ok());
    const qmomkit::PhaseSpaceHierarchy hierarchy({1, 1}, {1, 0});
    const qmomkit::CqmomCase direct_cqmom_case{
        qmomkit::extract_cqmom_moments({{2, 2}, {1.0, 1.5, 0.5, 0.75}}, qmomkit::create_cqmom_moment_layout(hierarchy)), {}, {true}};
    assert(configured.value().moments.layout.hierarchy.node_counts() == direct_cqmom_case.moments.layout.hierarchy.node_counts());
    assert(configured.value().moments.layout.hierarchy.conditioning_order() == direct_cqmom_case.moments.layout.hierarchy.conditioning_order());
    assert(configured.value().moments.values.size() == direct_cqmom_case.moments.values.size());
    for (std::size_t index = 0; index < direct_cqmom_case.moments.values.size(); ++index)
    {
        assert_near(configured.value().moments.values[index], direct_cqmom_case.moments.values[index]);
    }
    assert(configured.value().qmom_options.adaptive);

    const auto configured_result = qmomkit::compute_cqmom(configured.value());
    const auto direct_result = qmomkit::compute_cqmom(direct_cqmom_case);
    assert(configured_result.ok());
    assert(direct_result.ok());
    assert(configured_result.quadrature.weights() == direct_result.quadrature.weights());
    for (std::size_t dimension = 0; dimension < 2; ++dimension)
    {
        assert_near(configured_result.quadrature.abscissas()[dimension][0], direct_result.quadrature.abscissas()[dimension][0]);
    }
}

void test_cqmom_independent_product_initialization_matches_direct_api()
{
    const auto configured = qmomkit::config::parse_cqmom(R"(
[phase_space]
names = size, velocity
node_counts = 2, 2
conditioning_order = 1, 0

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
    assert(configured.value().moments.layout.hierarchy.node_counts() == std::vector<std::size_t>({2, 2}));
    assert(configured.value().moments.layout.hierarchy.conditioning_order() == std::vector<std::size_t>({1, 0}));
    assert(configured.value().qmom_options.adaptive);
    assert_cqmom_moments_near_tensor(configured.value().moments, expected.moments);
    assert_near(configured.value().moments.moment({0, 0}), 2.5);

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
    const auto tensor = qmomkit::config::initialize_cqmom_moment_tensor(initial.value());
    const qmomkit::TensorMomentResult expected_tensor = qmomkit::compute_raw_moments(
        qmomkit::IndependentProductDistribution{{qmomkit::BetaDistribution{2.0, 5.0}, qmomkit::LogNormalDistribution{0.1, 0.2}}, 4.0}, {1, 1});
    assert(expected_tensor.ok());
    assert(tensor.ok());
    assert_tensor_near(tensor.value(), expected_tensor.moments);

    const auto cqmom_case_result = qmomkit::config::create_cqmom_case(initial.value(), qmomkit::PhaseSpaceHierarchy({1, 1}));
    assert(cqmom_case_result.ok());
    assert_cqmom_moments_near_tensor(cqmom_case_result.value().moments, expected_tensor.moments);
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
    const auto initialized = qmomkit::config::initialize_cqmom_moment_tensor(initial.value());
    assert(expected.ok());
    assert(initialized.ok());
    assert_tensor_near(initialized.value(), expected.moments);
}

void test_cqmom_initialization_matches_example_flow()
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
    const auto configured = qmomkit::config::initialize_cqmom_moments(initial.value(), hierarchy);
    assert(configured.ok());

    const qmomkit::IndependentProductDistribution distribution{{qmomkit::GammaDistribution{3.0, 0.2}, qmomkit::NormalDistribution{-0.5, 1.25}}, 2.5};
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::RawMomentRequest moment_request{layout.required_indices.to_multi_indices()};
    const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(distribution, moment_request);
    assert(selected.ok());
    const qmomkit::CqmomMomentData expected = qmomkit::create_cqmom_moment_data(layout, selected);

    assert_cqmom_moments_near(configured.value(), expected);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(qmomkit::CqmomCase{configured.value()});
    assert(result.ok());
}

void test_cqmom_initialization_omits_maximum_orders()
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
    const auto initialized = qmomkit::config::initialize_cqmom_moments(initial.value(), hierarchy);
    assert(initialized.ok());
    assert(initialized.value().layout.hierarchy.node_counts() == std::vector<std::size_t>({2, 2}));
    assert_near(initialized.value().moment({0, 0}), 2.5);

    const auto cqmom_case_result = qmomkit::config::parse_cqmom(text);
    const auto cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(text);
    assert(cqmom_case_result.ok());
    assert(cqmom_prepared_case_result.ok());
    const qmomkit::CqmomResult cqmom_case_compute = qmomkit::compute_cqmom(cqmom_case_result.value());
    const qmomkit::CqmomResult cqmom_prepared_case_compute = qmomkit::compute_cqmom(cqmom_prepared_case_result.value());
    assert(cqmom_case_compute.ok());
    assert(cqmom_prepared_case_compute.ok());
    assert_cqmom_quadratures_near(cqmom_prepared_case_compute.quadrature, cqmom_case_compute.quadrature);
}

// The full tensor of raw moments has no extents of its own for an analytical
// distribution, so building it requires maximum_orders; the error says so.
void test_cqmom_moment_tensor_initialization_requires_maximum_orders()
{
    const auto initial = qmomkit::config::parse_cqmom_initial_condition(R"(
[initial]
type = independent_product
number_density = 1.0

[initial.dimension.0]
family = normal
mean = 0.0
standard_deviation = 1.0

[initial.dimension.1]
family = normal
mean = 0.0
standard_deviation = 1.0
)",
                                                                        "initialization-test");
    assert(initial.ok());
    const auto moment_tensor = qmomkit::config::initialize_cqmom_moment_tensor(initial.value(), "initialization-test");
    assert(!moment_tensor.ok());
    assert(moment_tensor.error() == "initialization-test: initial.maximum_orders is required to initialize the full tensor of raw moments.");
}

// CQMOM moments and the full tensor of raw moments come from an explicit
// tensor or an independent product; both initializers reject a
// one-dimensional moment sequence with the same message naming both sources.
void test_cqmom_moment_initialization_rejects_one_dimensional_moments()
{
    const auto initial = qmomkit::config::parse_qmom_initial_condition("[moments]\nvalues = 1.0, 0.0, 1.0, 0.0\n", "initialization-test");
    assert(initial.ok());
    const std::string expected =
        "initialization-test: CQMOM initialization requires explicit tensor moments or an independent product specification.";
    const auto moments = qmomkit::config::initialize_cqmom_moments(initial.value(), qmomkit::PhaseSpaceHierarchy({2, 2}), "initialization-test");
    assert(!moments.ok());
    assert(moments.error() == expected);
    const auto moment_tensor = qmomkit::config::initialize_cqmom_moment_tensor(initial.value(), "initialization-test");
    assert(!moment_tensor.ok());
    assert(moment_tensor.error() == expected);
}

// Configured CQMOM loading builds the moment set implied by the hierarchy, so
// supplying initial.maximum_orders or omitting it must produce the same case.
// The orders used below are the hierarchy's own extents, used only to compute
// the analytical reference; the loader does not derive them.
void test_cqmom_initialization_ignores_maximum_orders()
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
    const std::vector<std::size_t> hierarchy_orders{1, 1};
    assert(hierarchy_orders.size() == extents.size());
    for (std::size_t dimension = 0; dimension < extents.size(); ++dimension)
    {
        assert(hierarchy_orders[dimension] == extents[dimension] - 1);
    }

    const auto cqmom_case_result = qmomkit::config::create_cqmom_case(initial.value(), hierarchy);
    assert(cqmom_case_result.ok());

    const qmomkit::IndependentProductDistribution direct{{qmomkit::BetaDistribution{2.0, 5.0}, qmomkit::LogNormalDistribution{0.1, 0.2}}, 2.5};
    const qmomkit::TensorMomentResult expected = qmomkit::compute_raw_moments(direct, hierarchy_orders);
    assert(expected.ok());
    assert_cqmom_moments_near_tensor(cqmom_case_result.value().moments, expected.moments);

    const auto omitted_order = qmomkit::config::parse_cqmom(omitted_order_text);
    const auto explicit_order = qmomkit::config::parse_cqmom(explicit_order_text);
    assert(omitted_order.ok());
    assert(explicit_order.ok());
    assert(omitted_order.value().moments.layout.hierarchy.node_counts() == explicit_order.value().moments.layout.hierarchy.node_counts());
    assert(omitted_order.value().moments.layout.hierarchy.conditioning_order() ==
           explicit_order.value().moments.layout.hierarchy.conditioning_order());
    assert(omitted_order.value().abscissa_bounds.size() == explicit_order.value().abscissa_bounds.size());
    assert(!omitted_order.value().abscissa_bounds[0].has_value());
    assert(!explicit_order.value().abscissa_bounds[0].has_value());
    assert(omitted_order.value().abscissa_bounds[1].has_value());
    assert(explicit_order.value().abscissa_bounds[1].has_value());
    assert_near(omitted_order.value().abscissa_bounds[1]->lower, explicit_order.value().abscissa_bounds[1]->lower);
    assert_near(omitted_order.value().abscissa_bounds[1]->upper, explicit_order.value().abscissa_bounds[1]->upper);
    assert_cqmom_moments_near(omitted_order.value().moments, explicit_order.value().moments);
    assert_cqmom_moments_near_tensor(omitted_order.value().moments, expected.moments);
}

void test_cqmom_initialization_selects_explicit_moment_tensor()
{
    const auto initial = qmomkit::config::parse_cqmom_initial_condition(cqmom_text, "cqmom-test");
    assert(initial.ok());

    const qmomkit::PhaseSpaceHierarchy hierarchy({1, 1}, {1, 0});
    const auto initialized = qmomkit::config::initialize_cqmom_moments(initial.value(), hierarchy);
    assert(initialized.ok());

    const qmomkit::RawMomentTensor moment_tensor{{2, 2}, {1.0, 1.5, 0.5, 0.75}};
    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments(moment_tensor, qmomkit::create_cqmom_moment_layout(hierarchy));
    assert_cqmom_moments_near(initialized.value(), expected);
}

void test_cqmom_case_selects_explicit_moment_tensor()
{
    const auto cqmom_case_result = qmomkit::config::parse_cqmom(cqmom_text, "cqmom-test");
    assert(cqmom_case_result.ok());
    assert(cqmom_case_result.value().qmom_options.adaptive);

    const qmomkit::PhaseSpaceHierarchy hierarchy({1, 1}, {1, 0});
    const qmomkit::RawMomentTensor moment_tensor{{2, 2}, {1.0, 1.5, 0.5, 0.75}};
    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments(moment_tensor, qmomkit::create_cqmom_moment_layout(hierarchy));
    assert_cqmom_moments_near(cqmom_case_result.value().moments, expected);

    const auto initial = qmomkit::config::parse_cqmom_initial_condition(cqmom_text, "cqmom-test");
    assert(initial.ok());
    const auto made = qmomkit::config::create_cqmom_case(initial.value(), hierarchy, {true});
    assert(made.ok());
    assert_cqmom_moments_near(made.value().moments, expected);
    assert(made.value().qmom_options.adaptive);
}

void test_cqmom_prepared_case_uses_selected_analytical_moments()
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

    const auto cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(text);
    const auto cqmom_case_result = qmomkit::config::parse_cqmom(text);
    assert(cqmom_prepared_case_result.ok());
    assert(cqmom_case_result.ok());

    const qmomkit::CqmomPreparedCase cqmom_prepared_case = cqmom_prepared_case_result.value();
    assert(cqmom_prepared_case.hierarchy().node_counts() == std::vector<std::size_t>({2, 2}));
    assert(cqmom_prepared_case.hierarchy().conditioning_order() == std::vector<std::size_t>({1, 0}));
    assert(cqmom_prepared_case.diagnostic_level() == qmomkit::CqmomDiagnosticLevel::basic);

    const qmomkit::CqmomResult cqmom_prepared_case_compute = qmomkit::compute_cqmom(cqmom_prepared_case);
    const qmomkit::CqmomResult cqmom_case_compute = qmomkit::compute_cqmom(cqmom_case_result.value());
    assert(cqmom_prepared_case_compute.ok());
    assert(cqmom_case_compute.ok());
    assert_cqmom_quadratures_near(cqmom_prepared_case_compute.quadrature, cqmom_case_compute.quadrature);
}

void test_cqmom_prepared_case_accepts_explicit_moment_tensor()
{
    const auto cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(cqmom_text, "cqmom-test");
    const auto cqmom_case_result = qmomkit::config::parse_cqmom(cqmom_text, "cqmom-test");
    assert(cqmom_prepared_case_result.ok());
    assert(cqmom_case_result.ok());

    const qmomkit::CqmomPreparedCase cqmom_prepared_case = cqmom_prepared_case_result.value();
    assert(cqmom_prepared_case.hierarchy().node_counts() == std::vector<std::size_t>({1, 1}));
    assert(cqmom_prepared_case.hierarchy().conditioning_order() == std::vector<std::size_t>({1, 0}));
    assert(cqmom_prepared_case.qmom_options().adaptive);

    const qmomkit::CqmomResult cqmom_prepared_case_compute = qmomkit::compute_cqmom(cqmom_prepared_case);
    const qmomkit::CqmomResult cqmom_case_compute = qmomkit::compute_cqmom(cqmom_case_result.value());
    assert(cqmom_prepared_case_compute.ok());
    assert(cqmom_case_compute.ok());
    assert_cqmom_quadratures_near(cqmom_prepared_case_compute.quadrature, cqmom_case_compute.quadrature);
}

void test_cqmom_prepared_case_parses_coordinate_support()
{
    const std::string text = std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, 1.0
support_upper = none, 2.0

[cqmom]
diagnostic_level = basic
)";
    const auto cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(text, "cqmom-support");
    const auto cqmom_case_result = qmomkit::config::parse_cqmom(text, "cqmom-support");
    assert(cqmom_prepared_case_result.ok());
    assert(cqmom_case_result.ok());

    const qmomkit::CqmomPreparedCase cqmom_prepared_case = cqmom_prepared_case_result.value();
    assert(cqmom_prepared_case.abscissa_bounds().size() == 2);
    assert(!cqmom_prepared_case.abscissa_bounds()[0].has_value());
    assert(cqmom_prepared_case.abscissa_bounds()[1].has_value());
    assert_near(cqmom_prepared_case.abscissa_bounds()[1]->lower, 1.0);
    assert_near(cqmom_prepared_case.abscissa_bounds()[1]->upper, 2.0);
    assert(cqmom_case_result.value().abscissa_bounds.size() == cqmom_prepared_case.abscissa_bounds().size());

    const qmomkit::CqmomResult computed = qmomkit::compute_cqmom(cqmom_prepared_case);
    assert(computed.ok());
    assert(!computed.diagnostics.empty());
    assert(computed.diagnostics.front().phase_space_dimension == 1);
    assert(computed.diagnostics.front().bounded_support_checked);

    const qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
    assert(state.abscissa_bounds().size() == 2);
    assert(state.abscissa_bounds()[1].has_value());
    assert_near(state.abscissa_bounds()[1]->lower, 1.0);
    assert_near(state.abscissa_bounds()[1]->upper, 2.0);
}

void test_cqmom_factories_preserve_coordinate_support()
{
    const auto initial = qmomkit::config::parse_cqmom_initial_condition(cqmom_text, "cqmom-factory-support");
    assert(initial.ok());

    const qmomkit::PhaseSpaceHierarchy hierarchy({1, 1}, {1, 0});
    const auto defaulted = qmomkit::config::create_cqmom_case(initial.value(), hierarchy);
    assert(defaulted.ok());
    assert(defaulted.value().abscissa_bounds.empty());

    qmomkit::CqmomAbscissaBounds support(2);
    support[1] = qmomkit::AbscissaBounds{2.0, 5.0};

    const auto configured = qmomkit::config::create_cqmom_case(initial.value(), hierarchy, qmomkit::QmomOptions{false}, "support-factory", support);
    assert(configured.ok());
    assert_coordinate_support_matches_shifted_interval(configured.value().abscissa_bounds);
}

void test_cqmom_prepared_case_rejects_invalid_coordinate_support()
{
    const auto missing_upper = qmomkit::config::parse_cqmom_prepared(std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, 1.0
)",
                                                                     "cqmom-support");
    assert(!missing_upper.ok());
    assert(has_message_fragment(missing_upper, "support_lower and phase_space.support_upper must be configured together"));

    const auto wrong_size = qmomkit::config::parse_cqmom_prepared(std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, 1.0, 2.0
support_upper = none, 2.0, 3.0
)",
                                                                  "cqmom-support");
    assert(!wrong_size.ok());
    assert(has_message_fragment(wrong_size, "one entry per phase-space dimension"));

    const auto invalid_interval = qmomkit::config::parse_cqmom_prepared(std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, 2.0
support_upper = none, 1.0
)",
                                                                        "cqmom-support");
    assert(!invalid_interval.ok());
    assert(has_message_fragment(invalid_interval, "lower < upper"));

    const auto lower_nan = qmomkit::config::parse_cqmom_prepared(std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, nan
support_upper = none, 2.0
)",
                                                                 "cqmom-support");
    assert(!lower_nan.ok());
    assert(has_message_fragment(lower_nan, "support_lower must contain finite numbers or none entries"));

    const auto upper_inf = qmomkit::config::parse_cqmom_prepared(std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, 1.0
support_upper = none, inf
)",
                                                                 "cqmom-support");
    assert(!upper_inf.ok());
    assert(has_message_fragment(upper_inf, "support_upper must contain finite numbers or none entries"));

    const auto lower_negative_inf = qmomkit::config::parse_cqmom_prepared(std::string(cqmom_text) + R"(
[phase_space]
support_lower = none, -inf
support_upper = none, 2.0
)",
                                                                          "cqmom-support");
    assert(!lower_negative_inf.ok());
    assert(has_message_fragment(lower_negative_inf, "support_lower must contain finite numbers or none entries"));
}

void test_cqmom_prepared_case_initial_state_preserves_case_settings()
{
    const auto cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(cqmom_text, "cqmom-test");
    assert(cqmom_prepared_case_result.ok());

    const qmomkit::CqmomPreparedCase cqmom_prepared_case = cqmom_prepared_case_result.value();
    const qmomkit::CqmomState state = cqmom_prepared_case.initial_state();

    assert(state.hierarchy().node_counts() == cqmom_prepared_case.hierarchy().node_counts());
    assert(state.hierarchy().conditioning_order() == cqmom_prepared_case.hierarchy().conditioning_order());
    assert(state.qmom_options().adaptive == cqmom_prepared_case.qmom_options().adaptive);
    assert_near(state.qmom_options().rmin, cqmom_prepared_case.qmom_options().rmin);
    assert_near(state.qmom_options().eabs, cqmom_prepared_case.qmom_options().eabs);
    assert_near(state.qmom_options().cutoff, cqmom_prepared_case.qmom_options().cutoff);
    assert_near(state.conditional_options().singular_tolerance, cqmom_prepared_case.conditional_options().singular_tolerance);
    assert_near(state.conditional_options().conditioning_warning_threshold, cqmom_prepared_case.conditional_options().conditioning_warning_threshold);
    assert(state.qmom_options().adaptive == cqmom_prepared_case.qmom_options().adaptive);
    assert(state.diagnostic_level() == cqmom_prepared_case.diagnostic_level());
    assert(state.vandermonde_solver() == cqmom_prepared_case.vandermonde_solver());

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

qmomkit::CqmomAbscissaBounds cqmom_replacement_support()
{
    qmomkit::CqmomAbscissaBounds support(2);
    support[1] = qmomkit::AbscissaBounds{0.0, 2.0};
    return support;
}

void assert_cqmom_case_policies_equal(const qmomkit::CqmomPreparedCase& actual_cqmom_prepared_case,
                                      const qmomkit::CqmomPreparedCase& expected_cqmom_prepared_case)
{
    assert(actual_cqmom_prepared_case.hierarchy().node_counts() == expected_cqmom_prepared_case.hierarchy().node_counts());
    assert(actual_cqmom_prepared_case.hierarchy().conditioning_order() == expected_cqmom_prepared_case.hierarchy().conditioning_order());
    assert(actual_cqmom_prepared_case.qmom_options().adaptive == expected_cqmom_prepared_case.qmom_options().adaptive);
    assert_near(actual_cqmom_prepared_case.qmom_options().rmin, expected_cqmom_prepared_case.qmom_options().rmin);
    assert_near(actual_cqmom_prepared_case.qmom_options().eabs, expected_cqmom_prepared_case.qmom_options().eabs);
    assert_near(actual_cqmom_prepared_case.qmom_options().cutoff, expected_cqmom_prepared_case.qmom_options().cutoff);
    assert_near(actual_cqmom_prepared_case.conditional_options().singular_tolerance,
                expected_cqmom_prepared_case.conditional_options().singular_tolerance);
    assert_near(actual_cqmom_prepared_case.conditional_options().conditioning_warning_threshold,
                expected_cqmom_prepared_case.conditional_options().conditioning_warning_threshold);
    assert(actual_cqmom_prepared_case.qmom_options().adaptive == expected_cqmom_prepared_case.qmom_options().adaptive);
    assert(actual_cqmom_prepared_case.diagnostic_level() == expected_cqmom_prepared_case.diagnostic_level());
    assert(actual_cqmom_prepared_case.vandermonde_solver() == expected_cqmom_prepared_case.vandermonde_solver());
    assert(actual_cqmom_prepared_case.abscissa_bounds().size() == expected_cqmom_prepared_case.abscissa_bounds().size());
    for (std::size_t dimension = 0; dimension < actual_cqmom_prepared_case.abscissa_bounds().size(); ++dimension)
    {
        assert(actual_cqmom_prepared_case.abscissa_bounds()[dimension].has_value() ==
               expected_cqmom_prepared_case.abscissa_bounds()[dimension].has_value());
        if (actual_cqmom_prepared_case.abscissa_bounds()[dimension].has_value())
        {
            assert_near(actual_cqmom_prepared_case.abscissa_bounds()[dimension]->lower,
                        expected_cqmom_prepared_case.abscissa_bounds()[dimension]->lower);
            assert_near(actual_cqmom_prepared_case.abscissa_bounds()[dimension]->upper,
                        expected_cqmom_prepared_case.abscissa_bounds()[dimension]->upper);
        }
    }
}

qmomkit::PhaseSpaceHierarchy cqmom_replacement_hierarchy()
{
    return qmomkit::PhaseSpaceHierarchy({2, 2}, {1, 0});
}

// A case carrying the moments its hierarchy reconstructs from, with every other
// setting pinned so that replacement and prepared-metadata tests can vary one
// thing at a time.
qmomkit::CqmomCase replacement_base_cqmom_case()
{
    const qmomkit::PhaseSpaceHierarchy hierarchy = cqmom_replacement_hierarchy();
    return {qmomkit::extract_cqmom_moments(cqmom_replacement_tensor({0.75, 0.75, 1.5, 1.5}), qmomkit::create_cqmom_moment_layout(hierarchy)),
            {1.0e-12, 1.0e9},
            {true, 1.0e-9, 1.0e-11, 1.0e-28},
            qmomkit::CqmomDiagnosticLevel::basic,
            qmomkit::VandermondeSolver::dense,
            cqmom_replacement_support()};
}

// The same case rebuilt for a different hierarchy. Prepared metadata is keyed on
// the hierarchy, so this is how an incompatible case is produced now that the
// hierarchy cannot be reassigned on its own.
qmomkit::CqmomCase replacement_cqmom_case_with_hierarchy(const qmomkit::PhaseSpaceHierarchy& hierarchy)
{
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    qmomkit::CqmomCase cqmom_case = replacement_base_cqmom_case();
    cqmom_case.moments = qmomkit::CqmomMomentData{layout, std::vector<double>(layout.required_indices.size(), 1.0)};
    return cqmom_case;
}

void test_cqmom_prepared_case_replaces_initial_moments_immutably()
{
    const qmomkit::CqmomCase base_cqmom_case = replacement_base_cqmom_case();
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(cqmom_replacement_hierarchy());
    const qmomkit::CqmomCase reassembled_base_cqmom_case{
        base_cqmom_case.moments,          base_cqmom_case.conditional_options, base_cqmom_case.qmom_options,
        base_cqmom_case.diagnostic_level, base_cqmom_case.vandermonde_solver,  base_cqmom_case.abscissa_bounds};
    const qmomkit::CqmomPreparedCase base_cqmom_prepared_case =
        qmomkit::internal::CqmomConstructionAccess::create_prepared_case(reassembled_base_cqmom_case);
    qmomkit::RawMomentTensor replacement = cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0});
    const qmomkit::CqmomMomentData expected_moments = qmomkit::extract_cqmom_moments(replacement, layout);

    const auto transformed = qmomkit::config::with_initial_moments(base_cqmom_prepared_case, replacement);
    assert(transformed.ok());
    const qmomkit::CqmomPreparedCase replaced_cqmom_prepared_case = transformed.value();
    assert_cqmom_case_policies_equal(replaced_cqmom_prepared_case, base_cqmom_prepared_case);

    replacement.values[0] = -999.0;
    const qmomkit::CqmomState base_state = base_cqmom_prepared_case.initial_state();
    const qmomkit::CqmomState replaced_state = replaced_cqmom_prepared_case.initial_state();
    assert(base_state.moment_indices() == layout.required_indices.to_multi_indices());
    assert(replaced_state.moment_indices() == expected_moments.layout.required_indices.to_multi_indices());

    const qmomkit::CqmomResult base_result = qmomkit::compute_cqmom(base_cqmom_prepared_case);
    const qmomkit::CqmomResult replaced_result = qmomkit::compute_cqmom(replaced_cqmom_prepared_case);
    assert(base_result.ok());
    assert(replaced_result.ok());
    const qmomkit::CqmomMomentData canonical_base = qmomkit::extract_cqmom_moments(
        qmomkit::reconstruct_raw_moments(base_result.quadrature, cqmom_replacement_hierarchy().moment_extents()), layout);
    const qmomkit::CqmomMomentData canonical_replacement =
        qmomkit::extract_cqmom_moments(qmomkit::reconstruct_raw_moments(replaced_result.quadrature, replacement.extents), layout);
    assert(base_state.moment_values() == canonical_base.values);
    assert(replaced_state.moment_values() == canonical_replacement.values);

    // The replacement case shares the base case's prepared metadata rather than
    // building the branch structure again.
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(replaced_cqmom_prepared_case).get() ==
           qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(base_cqmom_prepared_case).get());
}

void test_cqmom_prepared_case_replaces_initial_moments_from_moment_data()
{
    const qmomkit::CqmomCase base_cqmom_case = replacement_base_cqmom_case();
    const qmomkit::CqmomPreparedCase base_cqmom_prepared_case = qmomkit::internal::CqmomConstructionAccess::create_prepared_case(base_cqmom_case);
    const qmomkit::CqmomMomentLayout& layout = base_cqmom_case.moments.layout;
    const std::vector<double> base_values = base_cqmom_prepared_case.initial_state().moment_values();

    const qmomkit::RawMomentTensor replacement_tensor = cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0});
    qmomkit::CqmomMomentData replacement = qmomkit::extract_cqmom_moments(replacement_tensor, layout);

    const auto transformed = qmomkit::config::with_initial_moments(base_cqmom_prepared_case, replacement);
    assert(transformed.ok());
    const qmomkit::CqmomPreparedCase replaced_cqmom_prepared_case = transformed.value();
    assert_cqmom_case_policies_equal(replaced_cqmom_prepared_case, base_cqmom_prepared_case);

    // The supplied moment data is copied into the returned case, so mutating it
    // afterwards must not reach either case.
    replacement.values[0] = -999.0;
    assert(base_cqmom_prepared_case.initial_state().moment_values() == base_values);

    const qmomkit::CqmomState replaced_state = replaced_cqmom_prepared_case.initial_state();
    assert(replaced_state.moment_indices() == layout.required_indices.to_multi_indices());

    // Replacing the moments from moment data must be indistinguishable from
    // supplying the same moments as a full tensor of raw moments: same
    // quadrature, same canonical state.
    const auto tensor_transformed = qmomkit::config::with_initial_moments(base_cqmom_prepared_case, replacement_tensor);
    assert(tensor_transformed.ok());
    const qmomkit::CqmomResult moment_data_result = qmomkit::compute_cqmom(replaced_cqmom_prepared_case);
    const qmomkit::CqmomResult tensor_result = qmomkit::compute_cqmom(tensor_transformed.value());
    assert(moment_data_result.ok());
    assert(tensor_result.ok());
    assert_cqmom_quadratures_near(moment_data_result.quadrature, tensor_result.quadrature);
    assert(replaced_state.moment_values() == tensor_transformed.value().initial_state().moment_values());

    const qmomkit::CqmomMomentData canonical_replacement =
        qmomkit::extract_cqmom_moments(qmomkit::reconstruct_raw_moments(moment_data_result.quadrature, replacement_tensor.extents), layout);
    assert(replaced_state.moment_values() == canonical_replacement.values);

    // A fixed lineage reuses the immutable prepared metadata rather than
    // rebuilding the branch graph for the replacement case.
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(replaced_cqmom_prepared_case).get() ==
           qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(base_cqmom_prepared_case).get());
}

void test_cqmom_prepared_case_moment_data_replacement_rejects_invalid_input()
{
    const qmomkit::CqmomCase base_cqmom_case = replacement_base_cqmom_case();
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = qmomkit::internal::CqmomConstructionAccess::create_prepared_case(base_cqmom_case);
    const qmomkit::CqmomMomentLayout& layout = base_cqmom_case.moments.layout;
    const std::vector<double> base_values = cqmom_prepared_case.initial_state().moment_values();

    const qmomkit::CqmomMomentData valid = qmomkit::extract_cqmom_moments(cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0}), layout);

    // A layout built for a different conditioning order is not the required
    // layout of the base hierarchy.
    const qmomkit::CqmomMomentLayout other_layout = qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy({2, 2}));
    const qmomkit::CqmomMomentData wrong_layout{other_layout, std::vector<double>(other_layout.required_indices.size(), 1.0)};
    const auto layout_mismatch = qmomkit::config::with_initial_moments(cqmom_prepared_case, wrong_layout);
    assert(!layout_mismatch.ok());
    // "layout" alone also appears in the value-count message, so the fragment
    // has to be specific enough to tell the two rejections apart.
    assert(has_message_fragment(layout_mismatch, "layout must match the required"));

    const qmomkit::CqmomMomentData wrong_value_count{layout, std::vector<double>(valid.values.size() + 1, 1.0)};
    const auto value_count = qmomkit::config::with_initial_moments(cqmom_prepared_case, wrong_value_count);
    assert(!value_count.ok());
    assert(has_message_fragment(value_count, "value count"));

    qmomkit::CqmomMomentData nonfinite = valid;
    nonfinite.values[0] = std::numeric_limits<double>::quiet_NaN();
    const auto finite = qmomkit::config::with_initial_moments(cqmom_prepared_case, nonfinite);
    assert(!finite.ok());
    assert(has_message_fragment(finite, "finite"));

    // Every rejection leaves the base case untouched.
    assert(cqmom_prepared_case.initial_state().moment_values() == base_values);
}

void test_cqmom_prepared_case_moment_data_replacement_defers_realizability()
{
    const qmomkit::CqmomCase base_cqmom_case = replacement_base_cqmom_case();
    const qmomkit::CqmomPreparedCase base_cqmom_prepared_case = qmomkit::internal::CqmomConstructionAccess::create_prepared_case(base_cqmom_case);
    qmomkit::CqmomMomentData nonrealizable =
        qmomkit::extract_cqmom_moments(cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0}), base_cqmom_case.moments.layout);
    // The layout enumerates m_[0,0] first. It is the total number density, so a
    // negative value is structurally valid but not realizable: replacement must
    // accept it and reconstruction must reject it.
    assert(nonrealizable.layout.required_indices[0] == qmomkit::MultiIndex({0, 0}));
    nonrealizable.values[0] = -1.0;

    const auto transformed = qmomkit::config::with_initial_moments(base_cqmom_prepared_case, std::move(nonrealizable));
    assert(transformed.ok());
    const qmomkit::CqmomResult reconstruction = qmomkit::compute_cqmom(transformed.value());
    assert(!reconstruction.ok());
    assert(reconstruction.status == qmomkit::CqmomStatus::qmom_failed);
}

void test_fixed_cqmom_prepared_metadata_ownership()
{
    qmomkit::CqmomCase cqmom_case = replacement_base_cqmom_case();
    cqmom_case.qmom_options.adaptive = false;
    const qmomkit::CqmomPreparedCase fixed_cqmom_prepared_case = qmomkit::internal::CqmomConstructionAccess::create_prepared_case(cqmom_case);
    const auto& fixed_metadata = qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(fixed_cqmom_prepared_case);
    assert(fixed_metadata != nullptr);

    const qmomkit::CqmomPreparedCase copied_cqmom_prepared_case = fixed_cqmom_prepared_case;
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(copied_cqmom_prepared_case).get() == fixed_metadata.get());

    qmomkit::CqmomPreparedCase movable_cqmom_prepared_case = fixed_cqmom_prepared_case;
    qmomkit::CqmomPreparedCase moved_cqmom_prepared_case = std::move(movable_cqmom_prepared_case);
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(moved_cqmom_prepared_case).get() == fixed_metadata.get());

    const qmomkit::CqmomPreparedCase distinct_cqmom_prepared_case = qmomkit::internal::CqmomConstructionAccess::create_prepared_case(cqmom_case);
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(distinct_cqmom_prepared_case).get() != fixed_metadata.get());

    const auto replaced = qmomkit::config::with_initial_moments(fixed_cqmom_prepared_case, cqmom_case.moments);
    assert(replaced.ok());
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(replaced.value()).get() == fixed_metadata.get());

    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(cqmom_replacement_hierarchy());
    const qmomkit::CqmomPreparedCase reassembled_cqmom_prepared_case = qmomkit::internal::CqmomConstructionAccess::create_prepared_case(
        qmomkit::CqmomCase{cqmom_case.moments, cqmom_case.conditional_options, cqmom_case.qmom_options, cqmom_case.diagnostic_level,
                           cqmom_case.vandermonde_solver, cqmom_case.abscissa_bounds});
    assert(qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(reassembled_cqmom_prepared_case) != nullptr);
    // The same case reached two ways: passed on directly, and reassembled field
    // by field. Both must prepare metadata and invert to the same quadrature.
    const qmomkit::CqmomResult direct_result = qmomkit::compute_cqmom(fixed_cqmom_prepared_case);
    const qmomkit::CqmomResult reassembled_result = qmomkit::compute_cqmom(reassembled_cqmom_prepared_case);
    assert(direct_result.ok());
    assert(reassembled_result.ok());
    assert_cqmom_quadratures_near(reassembled_result.quadrature, direct_result.quadrature);

    qmomkit::CqmomCase adaptive_cqmom_case = cqmom_case;
    adaptive_cqmom_case.qmom_options.adaptive = true;
    const qmomkit::CqmomPreparedCase adaptive_cqmom_prepared_case =
        qmomkit::internal::CqmomConstructionAccess::create_prepared_case(std::move(adaptive_cqmom_case));
    // An adaptive case prepares metadata too. Its branch tree carries a child
    // case for every parent node count Wheeler could realize, so it is a pure
    // function of the hierarchy and is built once rather than per reconstruction.
    const auto& adaptive_metadata = qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(adaptive_cqmom_prepared_case);
    assert(adaptive_metadata != nullptr);
    assert(adaptive_metadata->adaptive_branch_counts);

    const qmomkit::CqmomState initial_state = fixed_cqmom_prepared_case.initial_state();
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(initial_state).get() == fixed_metadata.get());
    const qmomkit::CqmomState copied_state = initial_state;
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(copied_state).get() == fixed_metadata.get());
    qmomkit::CqmomState movable_state = initial_state;
    qmomkit::CqmomState moved_state = std::move(movable_state);
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(moved_state).get() == fixed_metadata.get());

    std::optional<qmomkit::CqmomState> surviving_state;
    const qmomkit::internal::PreparedCqmomMetadata* surviving_metadata = nullptr;
    {
        const qmomkit::CqmomPreparedCase temporary_cqmom_prepared_case = qmomkit::internal::CqmomConstructionAccess::create_prepared_case(cqmom_case);
        surviving_metadata = qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(temporary_cqmom_prepared_case).get();
        surviving_state.emplace(temporary_cqmom_prepared_case.initial_state());
    }
    assert(surviving_state.has_value());
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(*surviving_state).get() == surviving_metadata);
    assert(surviving_state->hierarchy().node_counts() == cqmom_case.moments.layout.hierarchy.node_counts());
    assert(surviving_state->hierarchy().conditioning_order() == cqmom_case.moments.layout.hierarchy.conditioning_order());
}

void test_fixed_cqmom_prepared_metadata_is_concurrently_readable()
{
    qmomkit::CqmomCase cqmom_case = replacement_base_cqmom_case();
    cqmom_case.qmom_options.adaptive = false;
    const qmomkit::CqmomPreparedCase fixed_cqmom_prepared_case =
        qmomkit::internal::CqmomConstructionAccess::create_prepared_case(std::move(cqmom_case));
    const qmomkit::CqmomResult reference = qmomkit::compute_cqmom(fixed_cqmom_prepared_case);
    assert(reference.ok());

    std::vector<int> succeeded(4, 0);
    std::vector<std::thread> workers;
    workers.reserve(succeeded.size());
    for (std::size_t worker = 0; worker < succeeded.size(); ++worker)
    {
        workers.emplace_back(
            [&fixed_cqmom_prepared_case, &reference, &succeeded, worker]()
            {
                for (std::size_t iteration = 0; iteration < 16; ++iteration)
                {
                    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(fixed_cqmom_prepared_case);
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
    qmomkit::CqmomCase cqmom_case = replacement_base_cqmom_case();
    cqmom_case.qmom_options.adaptive = false;
    const qmomkit::CqmomPreparedCase fixed_cqmom_prepared_case = qmomkit::internal::CqmomConstructionAccess::create_prepared_case(cqmom_case);
    const auto& prepared_metadata = qmomkit::internal::CqmomConstructionAccess::case_prepared_metadata(fixed_cqmom_prepared_case);
    assert(prepared_metadata != nullptr);

    const qmomkit::CqmomState new_lineage_state = qmomkit::internal::CqmomConstructionAccess::create_state(cqmom_case);
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(new_lineage_state) != nullptr);
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(new_lineage_state).get() != prepared_metadata.get());

    const qmomkit::CqmomState descendant_state = qmomkit::internal::CqmomConstructionAccess::create_state(cqmom_case, prepared_metadata);
    assert(qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(descendant_state).get() == prepared_metadata.get());

    const auto assert_rejected = [&prepared_metadata](const qmomkit::CqmomCase& incompatible_cqmom_case)
    {
        bool rejected = false;
        try
        {
            static_cast<void>(qmomkit::internal::CqmomConstructionAccess::create_state(incompatible_cqmom_case, prepared_metadata));
        }
        catch (const std::invalid_argument& error)
        {
            rejected = true;
            assert(std::string(error.what()) == "Prepared CQMOM metadata is incompatible with the target CQMOM case.");
        }
        assert(rejected);
    };

    assert_rejected(replacement_cqmom_case_with_hierarchy(qmomkit::PhaseSpaceHierarchy({1, 2})));
    assert_rejected(replacement_cqmom_case_with_hierarchy(qmomkit::PhaseSpaceHierarchy({2, 2}, {0, 1})));
    assert_rejected(replacement_cqmom_case_with_hierarchy(qmomkit::PhaseSpaceHierarchy({1, 1, 1})));

    const qmomkit::CqmomMomentLayout fixed_layout = qmomkit::create_cqmom_moment_layout(cqmom_replacement_hierarchy());
    std::vector<qmomkit::MultiIndex> reordered_indices = fixed_layout.required_indices.to_multi_indices();
    std::reverse(reordered_indices.begin(), reordered_indices.end());
    qmomkit::CqmomMomentLayout reordered_layout{fixed_layout.hierarchy, qmomkit::MultiIndexTable(reordered_indices)};
    const qmomkit::CqmomCase reordered_cqmom_case{
        qmomkit::CqmomMomentData{std::move(reordered_layout), std::vector<double>(fixed_layout.required_indices.size(), 1.0)}};
    bool reordered_rejected = false;
    try
    {
        static_cast<void>(qmomkit::internal::CqmomConstructionAccess::create_state(reordered_cqmom_case, prepared_metadata));
    }
    catch (const std::invalid_argument& error)
    {
        reordered_rejected = true;
        assert(std::string(error.what()) == "Prepared CQMOM metadata is incompatible with the target CQMOM case.");
    }
    assert(reordered_rejected);

    bool malformed_new_lineage_rejected = false;
    try
    {
        static_cast<void>(qmomkit::internal::CqmomConstructionAccess::create_state(reordered_cqmom_case));
    }
    catch (const std::invalid_argument& error)
    {
        malformed_new_lineage_rejected = true;
        assert(std::string(error.what()).find("canonical CQMOM moment layout") != std::string::npos);
    }
    assert(malformed_new_lineage_rejected);

    qmomkit::CqmomCase adaptive_cqmom_case = cqmom_case;
    adaptive_cqmom_case.qmom_options.adaptive = true;
    // Metadata built for a fixed case must not be reused by an adaptive one:
    // its branch tree carries only the single configured parent node count,
    // not a case for every count Wheeler could realize. The compatibility check
    // rejects the mismatch.
    bool adaptive_rejected = false;
    try
    {
        static_cast<void>(qmomkit::internal::CqmomConstructionAccess::create_state(adaptive_cqmom_case, prepared_metadata));
    }
    catch (const std::invalid_argument& error)
    {
        adaptive_rejected = true;
        assert(std::string(error.what()).find("incompatible with the target") != std::string::npos);
    }
    assert(adaptive_rejected);
    // An adaptive case prepares its own metadata instead of rebuilding the
    // branch tree at every reconstruction.
    const qmomkit::CqmomState adaptive_state = qmomkit::internal::CqmomConstructionAccess::create_state(std::move(adaptive_cqmom_case));
    const auto& adaptive_metadata = qmomkit::internal::CqmomConstructionAccess::state_prepared_metadata(adaptive_state);
    assert(adaptive_metadata != nullptr);
    assert(adaptive_metadata->adaptive_branch_counts);
}

void test_cqmom_prepared_case_replacement_rejects_structural_tensor_errors()
{
    const qmomkit::CqmomCase base_cqmom_case = replacement_base_cqmom_case();
    const qmomkit::CqmomPreparedCase cqmom_prepared_case = qmomkit::internal::CqmomConstructionAccess::create_prepared_case(base_cqmom_case);
    const qmomkit::CqmomCase reassembled_base_cqmom_case{
        base_cqmom_case.moments,          base_cqmom_case.conditional_options, base_cqmom_case.qmom_options,
        base_cqmom_case.diagnostic_level, base_cqmom_case.vandermonde_solver,  base_cqmom_case.abscissa_bounds};
    const qmomkit::CqmomPreparedCase reassembled_cqmom_prepared_case =
        qmomkit::internal::CqmomConstructionAccess::create_prepared_case(reassembled_base_cqmom_case);
    const std::vector<double> base_values = cqmom_prepared_case.initial_state().moment_values();
    const std::vector<double> reassembled_base_values = reassembled_cqmom_prepared_case.initial_state().moment_values();

    const auto rank = qmomkit::config::with_initial_moments(cqmom_prepared_case, {{4}, std::vector<double>(4, 1.0)});
    assert(!rank.ok());
    assert(has_message_fragment(rank, "rank"));

    const auto invalid_extents = qmomkit::config::with_initial_moments(cqmom_prepared_case, {{4, 0}, {}});
    assert(!invalid_extents.ok());
    assert(has_message_fragment(invalid_extents, "extents"));

    const auto value_count = qmomkit::config::with_initial_moments(cqmom_prepared_case, {{4, 4}, std::vector<double>(15, 1.0)});
    assert(!value_count.ok());
    assert(has_message_fragment(value_count, "value count"));

    const auto overflow =
        qmomkit::config::with_initial_moments(cqmom_prepared_case, {{std::numeric_limits<std::size_t>::max(), 2}, std::vector<double>{}});
    assert(!overflow.ok());
    assert(has_message_fragment(overflow, "overflow"));

    qmomkit::RawMomentTensor nonfinite = cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0});
    nonfinite.values[0] = std::numeric_limits<double>::quiet_NaN();
    const auto finite = qmomkit::config::with_initial_moments(cqmom_prepared_case, nonfinite);
    assert(!finite.ok());
    assert(has_message_fragment(finite, "finite"));

    const auto missing_index = qmomkit::config::with_initial_moments(reassembled_cqmom_prepared_case, {{3, 4}, std::vector<double>(12, 1.0)});
    assert(!missing_index.ok());
    assert(has_message_fragment(missing_index, "required CQMOM moment index"));

    assert(cqmom_prepared_case.initial_state().moment_values() == base_values);
    assert(reassembled_cqmom_prepared_case.initial_state().moment_values() == reassembled_base_values);
}

void test_cqmom_prepared_case_replacement_defers_realizability()
{
    const qmomkit::CqmomPreparedCase base_cqmom_prepared_case =
        qmomkit::internal::CqmomConstructionAccess::create_prepared_case(replacement_base_cqmom_case());
    qmomkit::RawMomentTensor nonrealizable = cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0});
    nonrealizable.values[0] = -1.0;

    const auto transformed = qmomkit::config::with_initial_moments(base_cqmom_prepared_case, std::move(nonrealizable));
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

// Everything in a replaced case except its moments comes from the base case.
void assert_cqmom_case_settings_equal(const qmomkit::CqmomCase& actual_cqmom_case, const qmomkit::CqmomCase& expected_cqmom_case)
{
    assert(actual_cqmom_case.moments.layout.hierarchy.node_counts() == expected_cqmom_case.moments.layout.hierarchy.node_counts());
    assert(actual_cqmom_case.moments.layout.hierarchy.conditioning_order() == expected_cqmom_case.moments.layout.hierarchy.conditioning_order());
    assert(actual_cqmom_case.qmom_options.adaptive == expected_cqmom_case.qmom_options.adaptive);
    assert(actual_cqmom_case.qmom_options.rmin == expected_cqmom_case.qmom_options.rmin);
    assert(actual_cqmom_case.qmom_options.eabs == expected_cqmom_case.qmom_options.eabs);
    assert(actual_cqmom_case.qmom_options.cutoff == expected_cqmom_case.qmom_options.cutoff);
    assert(actual_cqmom_case.conditional_options.singular_tolerance == expected_cqmom_case.conditional_options.singular_tolerance);
    assert(actual_cqmom_case.conditional_options.conditioning_warning_threshold ==
           expected_cqmom_case.conditional_options.conditioning_warning_threshold);
    assert(actual_cqmom_case.conditional_options.moore_penrose_rcond == expected_cqmom_case.conditional_options.moore_penrose_rcond);
    assert(actual_cqmom_case.qmom_options.adaptive == expected_cqmom_case.qmom_options.adaptive);
    assert(actual_cqmom_case.diagnostic_level == expected_cqmom_case.diagnostic_level);
    assert(actual_cqmom_case.vandermonde_solver == expected_cqmom_case.vandermonde_solver);
    assert(actual_cqmom_case.abscissa_bounds.size() == expected_cqmom_case.abscissa_bounds.size());
    for (std::size_t dimension = 0; dimension < actual_cqmom_case.abscissa_bounds.size(); ++dimension)
    {
        assert(actual_cqmom_case.abscissa_bounds[dimension].has_value() == expected_cqmom_case.abscissa_bounds[dimension].has_value());
        if (actual_cqmom_case.abscissa_bounds[dimension].has_value())
        {
            assert(actual_cqmom_case.abscissa_bounds[dimension]->lower == expected_cqmom_case.abscissa_bounds[dimension]->lower);
            assert(actual_cqmom_case.abscissa_bounds[dimension]->upper == expected_cqmom_case.abscissa_bounds[dimension]->upper);
        }
    }
}

void assert_cqmom_moment_data_equal(const qmomkit::CqmomMomentData& actual, const qmomkit::CqmomMomentData& expected)
{
    assert(actual.layout.hierarchy.node_counts() == expected.layout.hierarchy.node_counts());
    assert(actual.layout.hierarchy.conditioning_order() == expected.layout.hierarchy.conditioning_order());
    assert(actual.layout.required_indices == expected.layout.required_indices);
    assert(actual.values == expected.values);
}

// A replacement through the open case must give exactly the case that the
// prepared overload holds, and therefore the same initial state.
template <class Moments>
void assert_cqmom_case_replacement_matches_prepared(const qmomkit::CqmomCase& replaced_cqmom_case, const qmomkit::CqmomCase& base_cqmom_case,
                                                    const Moments& moments)
{
    const auto replaced_cqmom_prepared_case_result = qmomkit::config::with_initial_moments(qmomkit::CqmomPreparedCase(base_cqmom_case), moments);
    assert(replaced_cqmom_prepared_case_result.ok());
    const qmomkit::CqmomCase& expected_cqmom_case =
        qmomkit::internal::CqmomConstructionAccess::case_definition(replaced_cqmom_prepared_case_result.value());
    assert_cqmom_case_settings_equal(replaced_cqmom_case, expected_cqmom_case);
    assert_cqmom_moment_data_equal(replaced_cqmom_case.moments, expected_cqmom_case.moments);
    assert(replaced_cqmom_case.initial_state().moment_values() == replaced_cqmom_prepared_case_result.value().initial_state().moment_values());
}

// Every setting is moved away from its default, so a replacement that dropped
// one would show.
qmomkit::CqmomCase fully_configured_replacement_base_cqmom_case()
{
    qmomkit::CqmomCase base_cqmom_case = replacement_base_cqmom_case();
    base_cqmom_case.qmom_options.adaptive = true;
    base_cqmom_case.conditional_options.moore_penrose_rcond = 1.0e-13;
    return base_cqmom_case;
}

// An open case's moments can be changed after it is built. Replacement reads
// only the new moments, so it must not depend on the old ones being valid.
qmomkit::CqmomCase with_broken_moments(qmomkit::CqmomCase cqmom_case)
{
    cqmom_case.moments.values.assign(2, std::numeric_limits<double>::quiet_NaN());
    return cqmom_case;
}

// An open base case selects the overloads that return an open case.
static_assert(std::is_same_v<decltype(qmomkit::config::with_initial_moments(std::declval<const qmomkit::CqmomCase&>(),
                                                                            std::declval<qmomkit::RawMomentTensor>())),
                             qmomkit::config::ConfigurationResult<qmomkit::CqmomCase>>);
static_assert(std::is_same_v<decltype(qmomkit::config::with_initial_moments(std::declval<const qmomkit::CqmomCase&>(),
                                                                            std::declval<qmomkit::CqmomMomentData>())),
                             qmomkit::config::ConfigurationResult<qmomkit::CqmomCase>>);

void test_cqmom_case_replaces_initial_moments_from_tensor()
{
    const qmomkit::CqmomCase base_cqmom_case = fully_configured_replacement_base_cqmom_case();
    const qmomkit::RawMomentTensor replacement = cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0});

    const auto replaced_cqmom_case_result = qmomkit::config::with_initial_moments(base_cqmom_case, replacement);
    assert(replaced_cqmom_case_result.ok());
    const qmomkit::CqmomCase& replaced_cqmom_case = replaced_cqmom_case_result.value();
    assert_cqmom_case_settings_equal(replaced_cqmom_case, base_cqmom_case);
    assert_cqmom_moment_data_equal(replaced_cqmom_case.moments, qmomkit::extract_cqmom_moments(replacement, base_cqmom_case.moments.layout));
    assert(replaced_cqmom_case.moments.values != base_cqmom_case.moments.values);
    assert_cqmom_case_replacement_matches_prepared(replaced_cqmom_case, base_cqmom_case, replacement);

    // Entries beyond those the hierarchy reads are ignored.
    qmomkit::RawMomentTensor larger{{5, 6}, std::vector<double>(30, 7.0)};
    for (std::size_t first = 0; first < 4; ++first)
    {
        for (std::size_t second = 0; second < 4; ++second)
        {
            larger.values[first * 6 + second] = replacement.values[first * 4 + second];
        }
    }
    const auto larger_cqmom_case_result = qmomkit::config::with_initial_moments(base_cqmom_case, larger);
    assert(larger_cqmom_case_result.ok());
    assert_cqmom_moment_data_equal(larger_cqmom_case_result.value().moments, replaced_cqmom_case.moments);

    const auto from_broken_base_result = qmomkit::config::with_initial_moments(with_broken_moments(base_cqmom_case), replacement);
    assert(from_broken_base_result.ok());
    assert_cqmom_case_settings_equal(from_broken_base_result.value(), replaced_cqmom_case);
    assert_cqmom_moment_data_equal(from_broken_base_result.value().moments, replaced_cqmom_case.moments);
}

void test_cqmom_case_replaces_initial_moments_from_moment_data()
{
    const qmomkit::CqmomCase base_cqmom_case = fully_configured_replacement_base_cqmom_case();
    const qmomkit::CqmomMomentData replacement =
        qmomkit::extract_cqmom_moments(cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0}), base_cqmom_case.moments.layout);

    const auto replaced_cqmom_case_result = qmomkit::config::with_initial_moments(base_cqmom_case, replacement);
    assert(replaced_cqmom_case_result.ok());
    const qmomkit::CqmomCase& replaced_cqmom_case = replaced_cqmom_case_result.value();
    assert_cqmom_case_settings_equal(replaced_cqmom_case, base_cqmom_case);
    assert_cqmom_moment_data_equal(replaced_cqmom_case.moments, replacement);
    assert_cqmom_case_replacement_matches_prepared(replaced_cqmom_case, base_cqmom_case, replacement);

    const auto from_broken_base_result = qmomkit::config::with_initial_moments(with_broken_moments(base_cqmom_case), replacement);
    assert(from_broken_base_result.ok());
    assert_cqmom_case_settings_equal(from_broken_base_result.value(), replaced_cqmom_case);
    assert_cqmom_moment_data_equal(from_broken_base_result.value().moments, replaced_cqmom_case.moments);
}

// The prepared overloads validate through the open-case overloads, so both
// reject a malformed replacement with the same message.
template <class Moments> std::string cqmom_replacement_error(const qmomkit::CqmomCase& base_cqmom_case, const Moments& moments)
{
    const auto replaced_cqmom_case_result = qmomkit::config::with_initial_moments(base_cqmom_case, moments);
    const auto replaced_cqmom_prepared_case_result = qmomkit::config::with_initial_moments(qmomkit::CqmomPreparedCase(base_cqmom_case), moments);
    assert(!replaced_cqmom_case_result.ok());
    assert(!replaced_cqmom_prepared_case_result.ok());
    assert(replaced_cqmom_case_result.error() == replaced_cqmom_prepared_case_result.error());
    return replaced_cqmom_case_result.error();
}

bool contains(const std::string& text, const std::string& fragment)
{
    return text.find(fragment) != std::string::npos;
}

void test_cqmom_case_replacement_rejects_invalid_moments()
{
    const qmomkit::CqmomCase base_cqmom_case = replacement_base_cqmom_case();

    const std::string extents_message = "The extents of the full tensor of raw moments must be positive and their product must not overflow.";
    assert(cqmom_replacement_error(base_cqmom_case, qmomkit::RawMomentTensor{{4}, std::vector<double>(4, 1.0)}) ==
           "The rank of the full tensor of raw moments must match the base CQMOM hierarchy.");
    assert(cqmom_replacement_error(base_cqmom_case, qmomkit::RawMomentTensor{{4, 0}, {}}) == extents_message);
    assert(cqmom_replacement_error(base_cqmom_case, qmomkit::RawMomentTensor{{4, 4}, std::vector<double>(15, 1.0)}) ==
           "The value count of the full tensor of raw moments does not match its extents.");
    assert(cqmom_replacement_error(base_cqmom_case, qmomkit::RawMomentTensor{{std::numeric_limits<std::size_t>::max(), 2}, std::vector<double>{}}) ==
           extents_message);
    qmomkit::RawMomentTensor nonfinite_tensor = cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0});
    nonfinite_tensor.values[5] = std::numeric_limits<double>::quiet_NaN();
    assert(cqmom_replacement_error(base_cqmom_case, nonfinite_tensor) == "The values of the full tensor of raw moments must be finite.");
    // The base hierarchy conditions on the second coordinate first, so each
    // coordinate needs four moments and a tensor short in either one lacks one.
    assert(cqmom_replacement_error(base_cqmom_case, qmomkit::RawMomentTensor{{3, 4}, std::vector<double>(12, 1.0)}) ==
           "The full tensor of raw moments does not contain required CQMOM moment index [3, 0].");
    assert(cqmom_replacement_error(base_cqmom_case, qmomkit::RawMomentTensor{{4, 3}, std::vector<double>(12, 1.0)}) ==
           "The full tensor of raw moments does not contain required CQMOM moment index [0, 3].");

    const std::string prefix = "Unable to replace CQMOM initial moments: ";
    const qmomkit::CqmomMomentLayout& layout = base_cqmom_case.moments.layout;
    const qmomkit::CqmomMomentLayout other_layout = qmomkit::create_cqmom_moment_layout(qmomkit::PhaseSpaceHierarchy({2, 2}));
    assert(cqmom_replacement_error(base_cqmom_case,
                                   qmomkit::CqmomMomentData{other_layout, std::vector<double>(other_layout.required_indices.size(), 1.0)}) ==
           prefix + "CQMOM moment data layout must match the required CQMOM moment layout of the base hierarchy.");
    assert(cqmom_replacement_error(base_cqmom_case, qmomkit::CqmomMomentData{layout, std::vector<double>(layout.required_indices.size() + 1, 1.0)}) ==
           prefix + "CQMOM moment data value count does not match its required CQMOM moment layout.");
    qmomkit::CqmomMomentData nonfinite_data = base_cqmom_case.moments;
    nonfinite_data.values[2] = std::numeric_limits<double>::infinity();
    assert(cqmom_replacement_error(base_cqmom_case, nonfinite_data) == prefix + "CQMOM moment data values must be finite.");
}

// Replacement checks structure only. Moments that are well formed but not
// realizable are accepted, and reconstructing them reports the failure.
void test_cqmom_case_replacement_defers_realizability()
{
    const qmomkit::CqmomCase base_cqmom_case = replacement_base_cqmom_case();
    qmomkit::RawMomentTensor nonrealizable_tensor = cqmom_replacement_tensor({1.0, 1.0, 2.0, 2.0});
    nonrealizable_tensor.values[0] = -1.0;
    qmomkit::CqmomMomentData nonrealizable_data = qmomkit::extract_cqmom_moments(nonrealizable_tensor, base_cqmom_case.moments.layout);

    const auto from_tensor = qmomkit::config::with_initial_moments(base_cqmom_case, nonrealizable_tensor);
    const auto from_moment_data = qmomkit::config::with_initial_moments(base_cqmom_case, nonrealizable_data);
    for (const qmomkit::config::ConfigurationResult<qmomkit::CqmomCase>* replaced_cqmom_case_result : {&from_tensor, &from_moment_data})
    {
        assert(replaced_cqmom_case_result->ok());
        const qmomkit::CqmomResult reconstruction = qmomkit::compute_cqmom(replaced_cqmom_case_result->value());
        assert(!reconstruction.ok());
        assert(reconstruction.status == qmomkit::CqmomStatus::qmom_failed);
        bool initial_state_failed = false;
        try
        {
            static_cast<void>(replaced_cqmom_case_result->value().initial_state());
        }
        catch (const std::runtime_error& error)
        {
            initial_state_failed = contains(error.what(), "CQMOM initial state quadrature reconstruction failed");
        }
        assert(initial_state_failed);
    }
}

void test_cqmom_prepared_case_file_loading()
{
    const std::filesystem::path filename = std::filesystem::temp_directory_path() / "qmomkit_cqmom_case_configuration_test.mk";
    {
        std::ofstream output(filename);
        output << cqmom_text;
    }
    const auto loaded = qmomkit::config::load_cqmom_prepared(filename);
    assert(loaded.ok());
    assert(loaded.value().qmom_options().adaptive);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(loaded.value());
    assert(result.ok());
    std::filesystem::remove(filename);

    const auto missing = qmomkit::config::load_cqmom_prepared(filename);
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
    const auto cqmom_prepared_case_result = qmomkit::config::load_cqmom_prepared(filename);
    assert(cqmom_prepared_case_result.ok());
    const qmomkit::CqmomResult independent_product_result = qmomkit::compute_cqmom(cqmom_prepared_case_result.value());
    assert(independent_product_result.ok());
    std::filesystem::remove(filename);
}

void test_cqmom_independent_product_zero_density_skips_marginal_validation()
{
    const auto configured = qmomkit::config::parse_cqmom(R"(
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
    assert(configured.value().moments.layout.hierarchy.node_counts() == std::vector<std::size_t>({1, 1}));
    for (const double value : configured.value().moments.values)
    {
        assert_near(value, 0.0);
    }
}

std::string cqmom_option_text(const std::string& phase_space_extra, const std::string& qmom_section)
{
    return "[phase_space]\nnode_counts = 1, 1\nconditioning_order = 0, 1\n" + phase_space_extra + "\n" + qmom_section + R"(
[moments]
extents = 2, 2
M[0,0] = 1.0
M[0,1] = 2.0
M[1,0] = 3.0
M[1,1] = 6.0
)";
}

void test_cqmom_adaptive_option_controls_the_reconstruction()
{
    // [qmom] adaptive alone decides whether CQMOM is adaptive, and it is
    // adaptive unless asked otherwise.
    const auto omitted = qmomkit::config::parse_cqmom(cqmom_option_text("", ""));
    assert(omitted.ok());
    assert(omitted.value().qmom_options.adaptive);

    const auto fixed = qmomkit::config::parse_cqmom(cqmom_option_text("", "[qmom]\nadaptive = false\n"));
    assert(fixed.ok());
    assert(!fixed.value().qmom_options.adaptive);
    const qmomkit::CqmomResult fixed_result = qmomkit::compute_cqmom(fixed.value());
    assert(fixed_result.ok());
    assert(fixed_result.quadrature.has_rectangular_storage());

    // The former branch-count key names its replacement, also when only the
    // initial condition is parsed.
    const std::string removed_text = cqmom_option_text("adaptive_branch_counts = false\n", "[qmom]\nadaptive = false\n");
    const auto removed = qmomkit::config::parse_cqmom(removed_text);
    assert(!removed.ok());
    assert(has_message_fragment(removed, "phase_space.adaptive_branch_counts was removed; use [qmom] adaptive."));
    const auto removed_initial = qmomkit::config::parse_cqmom_initial_condition(removed_text);
    assert(!removed_initial.ok());
    assert(has_message_fragment(removed_initial, "phase_space.adaptive_branch_counts was removed; use [qmom] adaptive."));
}

void test_cqmom_prepared_case_parsing_validates_coordinate_realizability_rank()
{
    const auto valid = qmomkit::config::parse_cqmom_prepared(cqmom_case_text_with_runtime(R"(

[time_integration.realizability.1]
check_mode = final_only
support = nonnegative
)"));
    assert(valid.ok());

    const auto invalid = qmomkit::config::parse_cqmom_prepared(cqmom_case_text_with_runtime(R"(

[time_integration.realizability.2]
check_mode = final_only
support = nonnegative
)"));
    assert(!invalid.ok());
    assert(has_failure(invalid));
    assert(has_message_fragment(invalid, "coordinate index 2"));
    assert(has_message_fragment(invalid, "phase-space rank 2"));
}

void test_huge_declared_extents_are_reported_as_incomplete()
{
    // The declared extents fix the size of the tensor, but only the supplied
    // entries can fill it. Extents far beyond them are an incomplete tensor and
    // are reported as one, not allocated: a 1e8 x 1e8 tensor would need 80 PB.
    for (const auto& [extents, missing] : {std::pair{std::string("100000000, 100000000"), std::string("9999999999999999 of its 10000000000000000")},
                                           std::pair{std::string("2, 1000000000000000"), std::string("1999999999999999 of its 2000000000000000")}})
    {
        const std::string text = "[phase_space]\nnode_counts = 1, 1\n[moments]\nextents = " + extents + "\nM[0,0] = 1.0\n";
        const auto open_case = qmomkit::config::parse_cqmom(text);
        const auto prepared_case = qmomkit::config::parse_cqmom_prepared(text);
        assert(!open_case.ok());
        assert(!prepared_case.ok());
        const std::string message = "The full tensor of raw moments is incomplete; it is missing " + missing + " entries.";
        assert(has_message_fragment(open_case, message));
        assert(has_message_fragment(prepared_case, message));
    }
}

void test_declared_rank_bounds_coordinate_realizability_indices()
{
    // Options read from a file that declares its phase space reject a
    // coordinate index beyond the rank, as the case loaders do, without sizing
    // the policy vector by the index. An index of 1e11 would otherwise ask for
    // about 6 TB.
    const std::string policy = "]\ncheck_mode = final_only\nsupport = nonnegative\n";
    for (const std::string index : {"2", "100000000000"})
    {
        const std::string text = "[phase_space]\nnode_counts = 2, 2\n[time_integration.realizability." + index + policy;
        const auto evolution_options = qmomkit::config::parse_cqmom_phase_space_evolution_options(text);
        const auto policies = qmomkit::config::parse_coordinate_realizability_policies(text);
        assert(!evolution_options.ok());
        assert(!policies.ok());
        assert(has_message_fragment(evolution_options, "coordinate index " + index + " is outside the configured phase-space rank 2"));
        assert(has_message_fragment(policies, "coordinate index " + index + " is outside the configured phase-space rank 2"));
    }

    const std::string within_rank_text = "[phase_space]\nnode_counts = 2, 2\n[time_integration.realizability.1" + policy;
    const auto within_rank = qmomkit::config::parse_cqmom_phase_space_evolution_options(within_rank_text);
    const auto within_rank_policies = qmomkit::config::parse_coordinate_realizability_policies(within_rank_text);
    assert(within_rank.ok());
    assert(within_rank_policies.ok());
    assert(within_rank.value().source_options.coordinate_realizability_policies.size() == 2);
    assert(within_rank.value().source_options.coordinate_realizability_policies[1].has_value());
    assert(within_rank_policies.value().size() == 2);
    assert(within_rank_policies.value()[1].has_value());

    // Every section is checked before the indices, as the case loaders do, so
    // all parsers report the same first error.
    const std::string misspelled = "[phase_space]\nnode_counts = 2, 2\n[time_integration.realizability.0]\ncheck_mode = final_only\n"
                                   "support = nonnegativ\n[time_integration.realizability.2" +
                                   policy;
    const std::string support_message = "time_integration.realizability.0.support must be";
    assert(has_message_fragment(qmomkit::config::parse_cqmom_phase_space_evolution_options(misspelled), support_message));
    assert(has_message_fragment(qmomkit::config::parse_coordinate_realizability_policies(misspelled), support_message));
}

void test_adaptive_cqmom_reduces_conditional_and_marginal_branches()
{
    const std::vector<double> conditional_weights{0.2, 0.3, 0.5};
    const std::vector<std::vector<double>> conditional_abscissas{{-1.0, 1.0, 1.0}, {0.0, -0.5, 0.5}};
    const auto conditional_reduction = qmomkit::config::parse_cqmom(cqmom_adaptive_text(conditional_weights, conditional_abscissas, {0, 1}));
    assert(conditional_reduction.ok());
    assert(conditional_reduction.value().qmom_options.adaptive);

    const qmomkit::CqmomResult conditional = qmomkit::compute_cqmom(conditional_reduction.value());
    assert(conditional.ok());
    assert(conditional.quadrature.adaptive_branch_counts_enabled());
    assert(conditional.quadrature.total_nodes() == 3);
    assert(std::any_of(conditional.diagnostics.begin(), conditional.diagnostics.end(), [](const qmomkit::CqmomDiagnostics& item)
                       { return item.conditioning_level == 1 && item.actual_node_count < item.requested_node_count; }));

    const auto marginal_reduction = qmomkit::config::parse_cqmom(cqmom_adaptive_text({0.5, 0.5}, {{0.0, 0.0}, {-1.0, 1.0}}, {0, 1}));
    assert(marginal_reduction.ok());

    const qmomkit::CqmomResult marginal = qmomkit::compute_cqmom(marginal_reduction.value());
    assert(marginal.ok());
    assert(marginal.quadrature.total_nodes() == 2);
    assert(std::any_of(marginal.diagnostics.begin(), marginal.diagnostics.end(), [](const qmomkit::CqmomDiagnostics& item)
                       { return item.conditioning_level == 0 && item.actual_node_count < item.requested_node_count; }));
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
    const auto omitted = qmomkit::config::parse_cqmom(cqmom_diagnostic_level_text(""));
    assert(omitted.ok());
    assert(omitted.value().diagnostic_level == qmomkit::CqmomDiagnosticLevel::none);

    const auto detailed = qmomkit::config::parse_cqmom(cqmom_diagnostic_level_text("detailed"));
    assert(detailed.ok());
    assert(detailed.value().diagnostic_level == qmomkit::CqmomDiagnosticLevel::detailed);

    const auto basic = qmomkit::config::parse_cqmom(cqmom_diagnostic_level_text("basic"));
    assert(basic.ok());
    assert(basic.value().diagnostic_level == qmomkit::CqmomDiagnosticLevel::basic);

    const auto none = qmomkit::config::parse_cqmom(cqmom_diagnostic_level_text("none"));
    assert(none.ok());
    assert(none.value().diagnostic_level == qmomkit::CqmomDiagnosticLevel::none);

    const auto invalid = qmomkit::config::parse_cqmom(cqmom_diagnostic_level_text("verbose"));
    assert(!invalid.ok());
    assert(has_failure(invalid));
    assert(has_message_fragment(invalid, "diagnostic_level"));
    assert(has_message_fragment(invalid, "detailed, basic, or none"));
}

std::string cqmom_vandermonde_solver_text(const std::string& solver, const std::string& rcond = "")
{
    std::string text = R"(
[phase_space]
node_counts = 1, 1

)";
    if (!solver.empty() || !rcond.empty())
    {
        text += "[vandermonde]\n";
        if (!solver.empty())
        {
            text += "solver = " + solver + "\n";
        }
        if (!rcond.empty())
        {
            text += "rcond = " + rcond + "\n";
        }
        text += "\n";
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
    const auto omitted = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text(""));
    assert(omitted.ok());
    assert(omitted.value().vandermonde_solver == qmomkit::VandermondeSolver::automatic);

    const auto dense = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text("dense"));
    assert(dense.ok());
    assert(dense.value().vandermonde_solver == qmomkit::VandermondeSolver::dense);

    const auto pivoted_qr = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text("pivoted_qr"));
    assert(pivoted_qr.ok());
    assert(pivoted_qr.value().vandermonde_solver == qmomkit::VandermondeSolver::pivoted_qr);

    const auto svd = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text("svd"));
    assert(svd.ok());
    assert(svd.value().vandermonde_solver == qmomkit::VandermondeSolver::svd);

    const auto moore_penrose = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text("moore_penrose", "0.25"));
    assert(moore_penrose.ok());
    assert(moore_penrose.value().vandermonde_solver == qmomkit::VandermondeSolver::moore_penrose);
    assert(std::abs(moore_penrose.value().conditional_options.moore_penrose_rcond - 0.25) <= 1.0e-15);

    const auto bp = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text("bjorck_pereyra"));
    assert(bp.ok());
    assert(bp.value().vandermonde_solver == qmomkit::VandermondeSolver::bjorck_pereyra);

    const auto automatic = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text("automatic"));
    assert(automatic.ok());
    assert(automatic.value().vandermonde_solver == qmomkit::VandermondeSolver::automatic);

    const auto invalid = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text("structured"));
    assert(!invalid.ok());
    assert(has_failure(invalid));
    assert(has_message_fragment(invalid, "vandermonde.solver"));
    assert(has_message_fragment(invalid, "svd"));
    assert(has_message_fragment(invalid, "pivoted_qr"));
}

void test_cqmom_vandermonde_rcond_parses_and_validates_for_all_solvers()
{
    const auto omitted = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text(""));
    assert(omitted.ok());
    assert(std::abs(omitted.value().conditional_options.moore_penrose_rcond - 1.0e-14) <= 1.0e-28);

    const auto zero = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text("moore_penrose", "0"));
    assert(zero.ok());
    assert(zero.value().conditional_options.moore_penrose_rcond == 0.0);

    const auto ignored_for_dense = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text("dense", "0.25"));
    assert(ignored_for_dense.ok());
    assert(std::abs(ignored_for_dense.value().conditional_options.moore_penrose_rcond - 0.25) <= 1.0e-15);

    const auto representation_neutral = qmomkit::config::parse_cqmom_prepared(cqmom_vandermonde_solver_text("moore_penrose", "0.25"));
    assert(representation_neutral.ok());
    assert(std::abs(representation_neutral.value().conditional_options().moore_penrose_rcond - 0.25) <= 1.0e-15);
    assert(std::abs(representation_neutral.value().initial_state().conditional_options().moore_penrose_rcond - 0.25) <= 1.0e-15);

    for (const std::string invalid_rcond : {"-0.1", "1", "1.1", "nan", "inf"})
    {
        const auto invalid = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text("dense", invalid_rcond));
        assert(!invalid.ok());
        assert(has_message_fragment(invalid, "vandermonde.rcond"));
    }
}

void assert_cqmom_case_reports_vandermonde_solver(const std::string& configured_solver, const qmomkit::VandermondeSolver expected_solver)
{
    const auto parsed = qmomkit::config::parse_cqmom_prepared(cqmom_vandermonde_solver_text(configured_solver));
    assert(parsed.ok());

    const qmomkit::CqmomPreparedCase cqmom_prepared_case = parsed.value();
    assert(cqmom_prepared_case.vandermonde_solver() == expected_solver);

    const qmomkit::CqmomState state = cqmom_prepared_case.initial_state();
    assert(state.vandermonde_solver() == expected_solver);
    assert(state.vandermonde_solver() == cqmom_prepared_case.vandermonde_solver());
}

void test_cqmom_prepared_case_and_state_report_vandermonde_solver()
{
    assert_cqmom_case_reports_vandermonde_solver("", qmomkit::VandermondeSolver::automatic);
    assert_cqmom_case_reports_vandermonde_solver("dense", qmomkit::VandermondeSolver::dense);
    assert_cqmom_case_reports_vandermonde_solver("pivoted_qr", qmomkit::VandermondeSolver::pivoted_qr);
    assert_cqmom_case_reports_vandermonde_solver("svd", qmomkit::VandermondeSolver::svd);
    assert_cqmom_case_reports_vandermonde_solver("moore_penrose", qmomkit::VandermondeSolver::moore_penrose);
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

    // Per-coordinate policies are the CQMOM setting; a QMOM case reads the
    // unnumbered policy instead and reports a numbered one.
    const auto qmom_case_result = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1

[moments]
values = 1.0, 0.0

[time_integration.realizability.0]
check_mode = final_only
support = real_line
)");
    assert(!qmom_case_result.ok());
    assert(has_message_fragment(qmom_case_result, "Unknown section 'time_integration.realizability.0'"));

    // Before the checked addition, SIZE_MAX + 1 wrapped to zero before an
    // unchecked policies[SIZE_MAX] access. The bounded regression verifies
    // the configuration error without attempting the unsafe old path.
    const auto impossible_index = qmomkit::config::parse_coordinate_realizability_policies("[time_integration.realizability." +
                                                                                           std::to_string(std::numeric_limits<std::size_t>::max()) +
                                                                                           "]\ncheck_mode = final_only\nsupport = real_line\n");
    assert(!impossible_index.ok());
    assert(has_message_fragment(impossible_index, "cannot be represented"));
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

    assert(qmomkit::config::parse_cqmom_prepared(first).ok());
    assert(qmomkit::config::parse_cqmom_prepared(second).ok());
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

    const auto without_policy_cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(without_policy);
    const auto with_policy_cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(with_policy);
    assert(without_policy_cqmom_prepared_case_result.ok());
    assert(with_policy_cqmom_prepared_case_result.ok());
    const qmomkit::CqmomAbscissaBounds& expected = without_policy_cqmom_prepared_case_result.value().abscissa_bounds();
    const qmomkit::CqmomAbscissaBounds& actual = with_policy_cqmom_prepared_case_result.value().abscissa_bounds();
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
    const auto qmom_case_result = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1

[moments]
values = 1.0, 0.0

[time_integration.realizability]
check_mode = final_only
support = real_line
)");
    assert(qmom_case_result.ok());
    assert(qmom_case_result.value().moments.values() == std::vector<double>({1.0, 0.0}));

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

    const auto qmom_case_result = qmomkit::config::parse_qmom(qmom_case_text_with_runtime(runtime));
    assert(qmom_case_result.ok());
    assert(qmom_case_result.value().node_count == 1);
    assert(qmom_case_result.value().qmom_options.adaptive);

    // The unnumbered policy is the QMOM setting; a CQMOM case reads the
    // per-coordinate policies instead and reports the unnumbered one.
    const auto cqmom_with_qmom_policy = qmomkit::config::parse_cqmom_prepared(cqmom_case_text_with_runtime(runtime));
    assert(!cqmom_with_qmom_policy.ok());
    assert(has_message_fragment(cqmom_with_qmom_policy, "Unknown section 'time_integration.realizability'"));

    std::string coordinate_runtime = runtime;
    const std::string unnumbered_header = "[time_integration.realizability]";
    coordinate_runtime.replace(coordinate_runtime.find(unnumbered_header), unnumbered_header.size(), "[time_integration.realizability.0]");
    const auto cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(cqmom_case_text_with_runtime(coordinate_runtime));
    assert(cqmom_prepared_case_result.ok());
    assert(cqmom_prepared_case_result.value().diagnostic_level() == qmomkit::CqmomDiagnosticLevel::none);

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
        const auto qmom_case_result = qmomkit::config::parse_qmom(qmom_case_text_with_runtime(block));
        assert(!qmom_case_result.ok());
        assert(has_failure(qmom_case_result));

        const auto cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(cqmom_case_text_with_runtime(block));
        assert(!cqmom_prepared_case_result.ok());
        assert(has_failure(cqmom_prepared_case_result));

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

    const auto cqmom = qmomkit::config::parse_cqmom_prepared(cqmom_case_text_with_runtime(removed_key));
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
    const auto cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(cqmom_case_text_with_runtime(""));
    assert(cqmom_prepared_case_result.ok());

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

    const qmomkit::PhaseSpaceVelocity velocity =
        qmomkit::make_phase_space_velocity([](double, const std::vector<double>& xi) { return std::vector<double>(xi.size(), 0.25); });
    const auto source = [](double, const qmomkit::CqmomSourceContext& context)
    { return qmomkit::CqmomMomentRate{std::vector<double>(context.size(), 0.0)}; };

    const qmomkit::CqmomPhaseSpaceEvolutionResult phase_space_evolved = qmomkit::advance_cqmom_phase_space(
        cqmom_prepared_case_result.value().initial_state(), 0.0, 0.1, velocity, source, configured_options.value());
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
    assert(options.value().generalized_moment_orders.empty());

    const auto empty_section = qmomkit::config::parse_output_options(R"(
[output]
)");
    assert(empty_section.ok());
    assert(empty_section.value().directory == std::filesystem::path("."));
    assert(!empty_section.value().write_moments);
    assert(!empty_section.value().write_quadrature);
    assert(!empty_section.value().interval.has_value());
    assert(empty_section.value().generalized_moment_orders.empty());
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
generalized_moments = 0.5, 0 | 1, 1.5 | 2, 0
)");
    assert(all_fields.ok());
    assert(all_fields.value().directory == std::filesystem::path("output"));
    assert(all_fields.value().write_moments);
    assert(!all_fields.value().write_quadrature);
    assert(all_fields.value().interval.has_value());
    assert_near(*all_fields.value().interval, 0.1);
    assert((all_fields.value().generalized_moment_orders == std::vector<std::vector<double>>{{0.5, 0.0}, {1.0, 1.5}, {2.0, 0.0}}));

    const auto one_dimensional = qmomkit::config::parse_output_options(R"(
[output]
generalized_moments = 0.5 | 2 | 1.5
)");
    assert(one_dimensional.ok());
    assert((one_dimensional.value().generalized_moment_orders == std::vector<std::vector<double>>{{0.5}, {2.0}, {1.5}}));
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

    const auto negative_order = qmomkit::config::parse_output_options(R"(
[output]
generalized_moments = 1, -0.5
)");
    assert(!negative_order.ok());
    assert(has_message_fragment(negative_order, "generalized_moments"));

    const auto nonfinite_order = qmomkit::config::parse_output_options(R"(
[output]
generalized_moments = 1, nan
)");
    assert(!nonfinite_order.ok());
    assert(has_message_fragment(nonfinite_order, "generalized_moments"));

    const auto empty_order = qmomkit::config::parse_output_options(R"(
[output]
    generalized_moments =
)");
    assert(!empty_order.ok());
    assert(has_failure(empty_order));
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

    const auto cqmom = qmomkit::config::parse_cqmom(R"(
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
        const auto qmom_case_result = qmomkit::config::parse_qmom(qmom_case_text_with_runtime(block));
        assert(!qmom_case_result.ok());
        assert(has_failure(qmom_case_result));

        const auto cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(cqmom_case_text_with_runtime(block));
        assert(!cqmom_prepared_case_result.ok());
        assert(has_failure(cqmom_prepared_case_result));
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
    const auto qmom_case_result = qmomkit::config::parse_qmom(qmom_case_text_with_runtime(source_policy));
    assert(!qmom_case_result.ok());
    assert(has_failure(qmom_case_result));

    const auto cqmom_prepared_case_result = qmomkit::config::parse_cqmom_prepared(cqmom_case_text_with_runtime(source_policy));
    assert(!cqmom_prepared_case_result.ok());
    assert(has_failure(cqmom_prepared_case_result));
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

    const auto cqmom = qmomkit::config::parse_cqmom(R"(
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

    const auto duplicate_moment =
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[moments]\nextents=2,2\nM[0,0]=1\nM[0,0]=1\nM[0,1]=0\nM[1,0]=0\nM[1,1]=0\n");
    assert(!duplicate_moment.ok());
    assert(has_failure(duplicate_moment));

    const auto missing_moment = qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[moments]\nextents=2,2\nM[0,0]=1\n");
    assert(!missing_moment.ok());
    assert(has_failure(missing_moment));

    const auto missing_initial_type = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\nfamily=normal\n");
    assert(!missing_initial_type.ok());
    assert(has_failure(missing_initial_type));

    // The node_counts entry (settings-phase validation) is well-formed here, so
    // only the initial-condition phase detects this error. Pinning the exact
    // source line guards against a future change accidentally handing that
    // phase a stale or unrelated ParsedFile instead of the one derived from
    // this text: a wrong-object substitution could keep the message text but
    // would very likely corrupt or zero out the reported line number.
    const auto invalid_initial_type = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=bogus\n");
    assert(!invalid_initial_type.ok());
    assert(has_message_fragment(invalid_initial_type, "<string>:4: Only '[initial] type = distribution' and '[initial] type = mixture'"));

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

    const auto missing_cqmom_initial_type = qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\nnumber_density=1\n");
    assert(!missing_cqmom_initial_type.ok());
    assert(has_failure(missing_cqmom_initial_type));

    // Same rationale as invalid_initial_type above, exercised through the
    // shared ParsedFile path used by parse_cqmom/parse_cqmom/parse_cqmom_prepared.
    const auto invalid_cqmom_initial_type = qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\ntype=bogus\n");
    assert(!invalid_cqmom_initial_type.ok());
    assert(has_message_fragment(invalid_cqmom_initial_type, "<string>:4: Only '[initial] type = independent_product' is supported"));

    const auto missing_product_dimensions =
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n");
    assert(!missing_product_dimensions.ok());
    assert(has_failure(missing_product_dimensions));

    const auto missing_product_family =
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
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

    const auto empty_section_name = qmomkit::config::parse_qmom("[]\nnode_counts=1\n[moments]\nvalues=1,0\n", "empty-section-name");
    assert(!empty_section_name.ok());
    assert(has_failure(empty_section_name));
    assert(has_message_fragment(empty_section_name, "Malformed section header."));
    assert(has_message_fragment(empty_section_name, "empty-section-name:1"));

    const auto nonfinite = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[moments]\nvalues=1,nan\n");
    assert(!nonfinite.ok());
    assert(has_failure(nonfinite));

    const auto malformed_adaptive = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[qmom]\nadaptive=yes\n[moments]\nvalues=1,0\n");
    assert(!malformed_adaptive.ok());
    assert(has_failure(malformed_adaptive));

    const auto invalid_order = qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\nconditioning_order=0,0\n[moments]\nextents=2,2\n"
                                                            "M[0,0]=1\nM[0,1]=1\nM[1,0]=1\nM[1,1]=1\n");
    assert(!invalid_order.ok());
    assert(has_failure(invalid_order));

    const auto cqmom_initial_section =
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[moments]\nextents=2,2\nM[0,0]=1\nM[0,1]=0\nM[1,0]=0\nM[1,1]=0\n[initial]\n");
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
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\ntype=distribution\nmaximum_order=1\nfamily=normal\nmean=0\n"
                                     "standard_deviation=1\n");
    assert(!unsupported_cqmom_initial_type.ok());
    assert(has_failure(unsupported_cqmom_initial_type));

    const auto invalid_product_order =
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=one\n"
                                     "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
                                     "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!invalid_product_order.ok());
    assert(has_failure(invalid_product_order));

    const auto product_tensor_size_overflow =
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\n"
                                     "maximum_orders=4294967295,4294967295\n"
                                     "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
                                     "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!product_tensor_size_overflow.ok());
    assert(has_failure(product_tensor_size_overflow));
    assert(
        has_message_fragment(product_tensor_size_overflow,
                             "The full tensor of raw moments given by initial.maximum_orders is too large: the product of its extents overflows."));

    const auto invalid_product_density =
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=-1\nmaximum_orders=1,1\n"
                                     "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
                                     "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!invalid_product_density.ok());
    assert(has_failure(invalid_product_density));

    const auto noncontiguous_product_dimensions =
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
                                     "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!noncontiguous_product_dimensions.ok());
    assert(has_failure(noncontiguous_product_dimensions));

    const auto extra_product_dimension =
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
                                     "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
                                     "[initial.dimension.1]\nfamily=normal\nmean=0\nstandard_deviation=1\n"
                                     "[initial.dimension.2]\nfamily=beta\nalpha=1\nbeta=2\n");
    assert(!extra_product_dimension.ok());
    assert(has_failure(extra_product_dimension));

    const auto duplicate_product_dimension =
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
                                     "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
                                     "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!duplicate_product_dimension.ok());
    assert(has_failure(duplicate_product_dimension));

    // Non-overlapping keys across the two blocks isolate the duplicate-section check
    // itself; identical keys would instead (also correctly) fail via the separate
    // duplicate-key check, without exercising the section-level detection.
    const auto duplicate_mixture_component = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                                                         "[initial.mixture.component.0]\nfamily=gamma\n"
                                                                         "[initial.mixture.component.0]\nnumber_density=1\nshape=1\nscale=2\n");
    assert(!duplicate_mixture_component.ok());
    assert(has_failure(duplicate_mixture_component));
    assert(has_message_fragment(duplicate_mixture_component, "Duplicate section 'initial.mixture.component.0'."));

    const auto invalid_marginal_parameter =
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
                                     "[initial.dimension.0]\nfamily=gamma\nshape=0\nscale=2\n"
                                     "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!invalid_marginal_parameter.ok());
    assert(has_failure(invalid_marginal_parameter));

    const auto product_dimension_mismatch =
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1\n"
                                     "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!product_dimension_mismatch.ok());
    assert(has_failure(product_dimension_mismatch));

    const auto product_rejects_marginal_density =
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
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
    const auto moments = qmomkit::config::initialize_qmom_moments(initial.value());
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

    const auto cqmom = qmomkit::config::parse_cqmom(R"(
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
    assert_cqmom_moments_near_tensor(cqmom.value().moments, qmomkit::RawMomentTensor{{2, 2}, {3.0, 1.0, 2.0, 4.0}});
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

    const auto cqmom = qmomkit::config::parse_cqmom(R"(
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

    const auto explicit_order = qmomkit::config::parse_cqmom(R"(
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
        qmomkit::config::parse_cqmom("[phase_space]\nnode_counts=1\n[moments]\nextents=2\nM[0]=1\nM[1]=0\n");
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
    const auto loaded_cqmom = qmomkit::config::load_cqmom(filename);
    assert(loaded_cqmom.ok());
    assert(loaded_cqmom.value().moments.layout.hierarchy.conditioning_order() == std::vector<std::size_t>({1, 0}));
    assert_near(loaded_cqmom.value().moments.moment({0, 0}), 1.0);
    std::filesystem::remove(filename);

    const auto missing_cqmom = qmomkit::config::load_cqmom(filename);
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
    const auto invalid_cqmom_runtime = qmomkit::config::load_cqmom_prepared(filename);
    assert(!invalid_cqmom_runtime.ok());
    assert(has_failure(invalid_cqmom_runtime));
    std::filesystem::remove(filename);
}

const char* const qmom_mixture_head = R"(
[phase_space]
node_counts = 1

[initial]
type = mixture

[initial.mixture.component.0]
family = normal
mean = 0.0
standard_deviation = 1.0
number_density = 1.0
)";

const char* const cqmom_product_head = R"(
[phase_space]
node_counts = 1, 1

[initial]
type = independent_product
number_density = 1.0

[initial.dimension.0]
family = normal
mean = 0.0
standard_deviation = 1.0

[initial.dimension.1]
family = normal
mean = 0.0
standard_deviation = 1.0
)";

// Indexed sections are looked up by their canonical index, so "00" names a
// section that nothing reads. It is reported instead of silently ignored.
void test_indexed_sections_require_canonical_indices()
{
    const auto mixture = qmomkit::config::parse_qmom(std::string(qmom_mixture_head) +
                                                     "\n[initial.mixture.component.00]\nfamily = normal\nmean = 5.0\nstandard_deviation = 1.0\n"
                                                     "number_density = 100.0\n");
    assert(!mixture.ok());
    assert(has_message_fragment(mixture, "Unknown section 'initial.mixture.component.00'"));

    const auto product = qmomkit::config::parse_cqmom_prepared(std::string(cqmom_product_head) +
                                                               "\n[initial.dimension.00]\nfamily = normal\nmean = 0.0\nstandard_deviation = 1.0\n");
    assert(!product.ok());
    assert(has_message_fragment(product, "Unknown section 'initial.dimension.00'"));

    const auto realizability = qmomkit::config::parse_coordinate_realizability_policies(
        "[time_integration.realizability.0]\ncheck_mode = final_only\nsupport = real_line\n"
        "[time_integration.realizability.00]\ncheck_mode = every_stage\nsupport = real_line\n");
    assert(!realizability.ok());
    assert(has_message_fragment(realizability, "Unknown section 'time_integration.realizability.00'"));
}

// Every section and key under an indexed prefix must belong to an indexed
// section the initialization reads; anything else is reported.
void test_indexed_entries_that_nothing_reads_are_reported()
{
    const auto garbage_component = qmomkit::config::parse_qmom(std::string(qmom_mixture_head) + "\n[initial.mixture.component.x]\nbogus = 3\n");
    assert(has_message_fragment(garbage_component, "Unknown section 'initial.mixture.component.x'"));

    const auto component_without_section =
        qmomkit::config::parse_qmom(std::string(qmom_mixture_head) + "\n[initial]\nmixture.component.1.bogus = 1\n");
    assert(has_message_fragment(component_without_section, "Unknown key 'initial.mixture.component.1.bogus'"));

    const auto component_for_distribution = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1

[initial]
type = distribution
family = normal
mean = 0.0
standard_deviation = 1.0

[initial.mixture.component.x]
family = normal
)");
    assert(has_message_fragment(component_for_distribution, "Unknown section 'initial.mixture.component.x'"));

    const auto garbage_dimension = qmomkit::config::parse_cqmom_prepared(std::string(cqmom_product_head) + "\n[initial.dimension.x]\nmean = 1.0\n");
    assert(has_message_fragment(garbage_dimension, "Unknown section 'initial.dimension.x'"));

    const auto garbage_marginal_component = qmomkit::config::parse_cqmom_prepared(R"(
[phase_space]
node_counts = 1, 1

[initial]
type = independent_product
number_density = 1.0

[initial.dimension.0]
type = mixture

[initial.dimension.0.mixture.component.0]
family = normal
mean = 0.0
standard_deviation = 1.0
number_density = 1.0

[initial.dimension.0.mixture.component.x]
whatever = 1

[initial.dimension.1]
family = normal
mean = 0.0
standard_deviation = 1.0
)");
    assert(has_message_fragment(garbage_marginal_component, "Unknown section 'initial.dimension.0.mixture.component.x'"));

    // The well-formed files above still load without the extra entries.
    assert(qmomkit::config::parse_qmom(qmom_mixture_head).ok());
    assert(qmomkit::config::parse_cqmom_prepared(cqmom_product_head).ok());
}

// A case loader validates the output settings of the file it reads, as the
// standalone output parser does.
void test_case_loaders_validate_output_settings()
{
    const std::string invalid_interval = "\n[output]\ninterval = -1\n";
    const auto qmom_case = qmomkit::config::parse_qmom(qmom_case_text_with_runtime(invalid_interval));
    assert(has_message_fragment(qmom_case, "output.interval must be positive."));
    const auto cqmom_case = qmomkit::config::parse_cqmom_prepared(cqmom_case_text_with_runtime(invalid_interval));
    assert(has_message_fragment(cqmom_case, "output.interval must be positive."));
    const auto output = qmomkit::config::parse_output_options(invalid_interval);
    assert(has_message_fragment(output, "output.interval must be positive."));

    // A trailing '|' is an empty exponent vector, as a trailing ',' is an empty list entry.
    const std::string trailing_separator = "\n[output]\ngeneralized_moments = 0.5 |\n";
    assert(has_message_fragment(qmomkit::config::parse_output_options(trailing_separator), "output.generalized_moments"));
    assert(has_message_fragment(qmomkit::config::parse_qmom(qmom_case_text_with_runtime(trailing_separator)), "output.generalized_moments"));

    assert(qmomkit::config::parse_qmom(qmom_case_text_with_runtime("\n[output]\ninterval = 0.5\ngeneralized_moments = 0.5 | 1.5\n")).ok());
}

// QMOM reads the unnumbered realizability policy and CQMOM the per-coordinate
// policies. Each parser reports the other formulation's settings, which would
// otherwise be accepted and have no effect.
void test_realizability_settings_match_the_formulation()
{
    const std::string unnumbered = "[time_integration.realizability]\ncheck_mode = final_only\nsupport = real_line\n";
    const std::string numbered = "[time_integration.realizability.0]\ncheck_mode = final_only\nsupport = real_line\n";

    assert(qmomkit::config::parse_qmom_phase_space_evolution_options(unnumbered).ok());
    assert(has_message_fragment(qmomkit::config::parse_qmom_phase_space_evolution_options(numbered),
                                "Unknown section 'time_integration.realizability.0'"));
    assert(qmomkit::config::parse_cqmom_phase_space_evolution_options(numbered).ok());
    assert(has_message_fragment(qmomkit::config::parse_cqmom_phase_space_evolution_options(unnumbered),
                                "Unknown section 'time_integration.realizability'"));
    assert(has_message_fragment(qmomkit::config::parse_coordinate_realizability_policies(unnumbered),
                                "Unknown section 'time_integration.realizability'"));

    // A per-coordinate key written inside the unnumbered QMOM section.
    const auto policy = qmomkit::config::parse_realizability_policy("[time_integration.realizability]\ncheck_mode = final_only\nsupport = real_line\n"
                                                                    "0.check_mode = every_stage\n");
    assert(has_message_fragment(policy, "Unknown key 'time_integration.realizability.0.check_mode'"));

    // Parsers that read neither setting accept both.
    assert(qmomkit::config::parse_phase_space_evolution_scheme(unnumbered + numbered).ok());
    assert(qmomkit::config::parse_node_advection_options(unnumbered + numbered).ok());
}

// Requests far beyond any node count cannot be allocated. They are reported
// as configuration errors instead of escaping as exceptions.
void test_oversized_requests_are_reported_not_thrown()
{
    const auto maximum_order = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1

[initial]
type = distribution
family = normal
mean = 0.0
standard_deviation = 1.0
maximum_order = 1000000000000000
)");
    assert(has_message_fragment(maximum_order, "moment generation failed"));

    const auto node_count = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1000000000000000

[initial]
type = distribution
family = normal
mean = 0.0
standard_deviation = 1.0
)");
    assert(has_message_fragment(node_count, "moment generation failed"));

    const auto coordinate_index = qmomkit::config::parse_coordinate_realizability_policies(
        "[time_integration.realizability.100000000000000000]\ncheck_mode = final_only\nsupport = real_line\n");
    assert(has_message_fragment(coordinate_index, "Coordinate realizability policy index is too large."));
}

// Numbers are plain decimal and are read in the classic "C" locale, so a host
// application's locale cannot change how a configuration file is read. The
// locale part runs where a decimal-comma locale is installed.
void test_numbers_are_read_independently_of_the_locale()
{
    const std::string text = "[phase_space]\nnode_counts = 1\n[moments]\nvalues = 1.5, 0.25\n";

    const std::string previous_locale = std::setlocale(LC_NUMERIC, nullptr);
    for (const char* const name : {"de_DE.UTF-8", "de_DE.utf8", "fr_FR.UTF-8", "fr_FR.utf8"})
    {
        if (std::setlocale(LC_NUMERIC, name) != nullptr && std::string(std::localeconv()->decimal_point) == ",")
        {
            const auto qmom_case = qmomkit::config::parse_qmom(text);
            assert(qmom_case.ok());
            assert(qmom_case.value().moments[0] == 1.5);
            assert(qmom_case.value().moments[1] == 0.25);
            break;
        }
    }
    std::setlocale(LC_NUMERIC, previous_locale.c_str());

    const auto qmom_case = qmomkit::config::parse_qmom(text);
    assert(qmom_case.ok());
    assert(qmom_case.value().moments[0] == 1.5);
    assert(qmomkit::config::parse_qmom("[phase_space]\nnode_counts = 1\n[moments]\nvalues = +1.0, -2.5e-1\n").ok());
    for (const char* const value : {"0x1p0", "1.5.0", "1.5e", "inf", "1e400"})
    {
        const auto rejected = qmomkit::config::parse_qmom("[phase_space]\nnode_counts = 1\n[moments]\nvalues = 1.0, " + std::string(value) + "\n");
        assert(has_message_fragment(rejected, "QMOM moments must be a comma-separated list of finite numbers."));
    }
}

// Each error names the entry that is actually wrong, on its own line.
void test_errors_name_the_offending_entry()
{
    const auto missing_initialization = qmomkit::config::parse_qmom("[phase_space]\nnode_counts = 1\n");
    assert(has_message_fragment(missing_initialization, "QMOM configuration requires either '[moments] values'"));

    const std::string upper_count_text = cqmom_case_text_with_runtime("") + "[phase_space]\nsupport_lower = 0.0, 0.0\nsupport_upper = 1.0\n";
    const std::size_t upper_line =
        static_cast<std::size_t>(std::count(upper_count_text.begin(), upper_count_text.begin() + upper_count_text.find("support_upper"), '\n')) + 1;
    const auto upper_count = qmomkit::config::parse_cqmom_prepared(upper_count_text);
    assert(has_message_fragment(upper_count, ":" + std::to_string(upper_line) +
                                                 ": phase_space.support_lower and phase_space.support_upper must contain one entry"));

    const auto huge_moment_index = qmomkit::config::parse_cqmom_prepared(R"(
[phase_space]
node_counts = 1, 1

[moments]
extents = 2, 2
M[99999999999999999999,0] = 1.0
)");
    assert(has_message_fragment(huge_moment_index, "Invalid moment index 'M[99999999999999999999,0]'"));

    const auto tensor_rank = qmomkit::config::parse_cqmom_prepared(R"(
[phase_space]
node_counts = 1, 1

[moments]
extents = 2, 2, 2
M[0,0] = 1.0
)");
    assert(has_message_fragment(tensor_rank, "does not match the number of entries in moments.extents."));

    const auto marginal_count = qmomkit::config::parse_cqmom_prepared(R"(
[phase_space]
node_counts = 1, 1, 1

[initial]
type = independent_product
number_density = 1.0

[initial.dimension.0]
family = normal
mean = 0.0
standard_deviation = 1.0

[initial.dimension.1]
family = normal
mean = 0.0
standard_deviation = 1.0
)");
    assert(has_message_fragment(marginal_count, ":3: Independent product initialization must define one '[initial.dimension.N]' section"));
}

} // namespace

int main(const int argc, char** argv)
{
    return qmomkit::test::dispatch_test_case(
        argc, argv,
        {
            {"configuration_result_api", test_configuration_result_api},
            {"valid_qmom_matches_programmatic_case", test_valid_qmom_matches_programmatic_case},
            {"qmom_case_parses_abscissa_bounds", test_qmom_case_parses_abscissa_bounds},
            {"qmom_distribution_initialization_matches_direct_api", test_qmom_distribution_initialization_matches_direct_api},
            {"all_supported_distribution_families_parse", test_all_supported_distribution_families_parse},
            {"qmom_mixture_initialization_matches_direct_api", test_qmom_mixture_initialization_matches_direct_api},
            {"mixture_zero_density_component_is_skipped", test_mixture_zero_density_component_is_skipped},
            {"mixture_level_distribution_keys_are_reported", test_mixture_level_distribution_keys_are_reported},
            {"cqmom_rejects_the_removed_moment_storage_key", test_cqmom_rejects_the_removed_moment_storage_key},
            {"valid_cqmom_matches_programmatic_case", test_valid_cqmom_matches_programmatic_case},
            {"cqmom_independent_product_initialization_matches_direct_api", test_cqmom_independent_product_initialization_matches_direct_api},
            {"cqmom_independent_product_mixture_marginal_initialization", test_cqmom_independent_product_mixture_marginal_initialization},
            {"cqmom_initialization_matches_example_flow", test_cqmom_initialization_matches_example_flow},
            {"cqmom_initialization_omits_maximum_orders", test_cqmom_initialization_omits_maximum_orders},
            {"cqmom_initialization_ignores_maximum_orders", test_cqmom_initialization_ignores_maximum_orders},
            {"cqmom_moment_tensor_initialization_requires_maximum_orders", test_cqmom_moment_tensor_initialization_requires_maximum_orders},
            {"cqmom_moment_initialization_rejects_one_dimensional_moments", test_cqmom_moment_initialization_rejects_one_dimensional_moments},
            {"cqmom_initialization_selects_explicit_moment_tensor", test_cqmom_initialization_selects_explicit_moment_tensor},
            {"cqmom_case_selects_explicit_moment_tensor", test_cqmom_case_selects_explicit_moment_tensor},
            {"cqmom_prepared_case_uses_selected_analytical_moments", test_cqmom_prepared_case_uses_selected_analytical_moments},
            {"cqmom_prepared_case_accepts_explicit_moment_tensor", test_cqmom_prepared_case_accepts_explicit_moment_tensor},
            {"cqmom_prepared_case_parses_coordinate_support", test_cqmom_prepared_case_parses_coordinate_support},
            {"cqmom_factories_preserve_coordinate_support", test_cqmom_factories_preserve_coordinate_support},
            {"cqmom_prepared_case_rejects_invalid_coordinate_support", test_cqmom_prepared_case_rejects_invalid_coordinate_support},
            {"cqmom_prepared_case_initial_state_preserves_case_settings", test_cqmom_prepared_case_initial_state_preserves_case_settings},
            {"cqmom_prepared_case_replaces_initial_moments_immutably", test_cqmom_prepared_case_replaces_initial_moments_immutably},
            {"cqmom_prepared_case_replaces_initial_moments_from_moment_data", test_cqmom_prepared_case_replaces_initial_moments_from_moment_data},
            {"cqmom_prepared_case_moment_data_replacement_rejects_invalid_input",
             test_cqmom_prepared_case_moment_data_replacement_rejects_invalid_input},
            {"cqmom_prepared_case_moment_data_replacement_defers_realizability",
             test_cqmom_prepared_case_moment_data_replacement_defers_realizability},
            {"fixed_cqmom_prepared_metadata_ownership", test_fixed_cqmom_prepared_metadata_ownership},
            {"fixed_cqmom_prepared_metadata_is_concurrently_readable", test_fixed_cqmom_prepared_metadata_is_concurrently_readable},
            {"fixed_cqmom_state_prepared_metadata_invariant_and_compatibility", test_fixed_cqmom_state_prepared_metadata_invariant_and_compatibility},
            {"cqmom_prepared_case_replacement_rejects_structural_tensor_errors",
             test_cqmom_prepared_case_replacement_rejects_structural_tensor_errors},
            {"cqmom_prepared_case_replacement_defers_realizability", test_cqmom_prepared_case_replacement_defers_realizability},
            {"cqmom_case_replaces_initial_moments_from_tensor", test_cqmom_case_replaces_initial_moments_from_tensor},
            {"cqmom_case_replaces_initial_moments_from_moment_data", test_cqmom_case_replaces_initial_moments_from_moment_data},
            {"cqmom_case_replacement_rejects_invalid_moments", test_cqmom_case_replacement_rejects_invalid_moments},
            {"cqmom_case_replacement_defers_realizability", test_cqmom_case_replacement_defers_realizability},
            {"cqmom_prepared_case_file_loading", test_cqmom_prepared_case_file_loading},
            {"cqmom_independent_product_zero_density_skips_marginal_validation",
             test_cqmom_independent_product_zero_density_skips_marginal_validation},
            {"cqmom_adaptive_option_controls_the_reconstruction", test_cqmom_adaptive_option_controls_the_reconstruction},
            {"huge_declared_extents_are_reported_as_incomplete", test_huge_declared_extents_are_reported_as_incomplete},
            {"declared_rank_bounds_coordinate_realizability_indices", test_declared_rank_bounds_coordinate_realizability_indices},
            {"cqmom_prepared_case_parsing_validates_coordinate_realizability_rank",
             test_cqmom_prepared_case_parsing_validates_coordinate_realizability_rank},
            {"adaptive_cqmom_reduces_conditional_and_marginal_branches", test_adaptive_cqmom_reduces_conditional_and_marginal_branches},
            {"cqmom_diagnostic_level_parses_values", test_cqmom_diagnostic_level_parses_values},
            {"cqmom_vandermonde_solver_parses_values", test_cqmom_vandermonde_solver_parses_values},
            {"cqmom_vandermonde_rcond_parses_and_validates_for_all_solvers", test_cqmom_vandermonde_rcond_parses_and_validates_for_all_solvers},
            {"cqmom_prepared_case_and_state_report_vandermonde_solver", test_cqmom_prepared_case_and_state_report_vandermonde_solver},
            {"realizability_policy_defaults_to_disabled", test_realizability_policy_defaults_to_disabled},
            {"realizability_policy_parses_check_modes", test_realizability_policy_parses_check_modes},
            {"realizability_policy_parses_support_and_bounded_interval", test_realizability_policy_parses_support_and_bounded_interval},
            {"coordinate_realizability_policies_parse_sparse_numbered_sections",
             test_coordinate_realizability_policies_parse_sparse_numbered_sections},
            {"quadrature_bounds_do_not_alter_coordinate_realizability_policy", test_quadrature_bounds_do_not_alter_coordinate_realizability_policy},
            {"coordinate_realizability_policy_does_not_alter_quadrature_bounds",
             test_coordinate_realizability_policy_does_not_alter_quadrature_bounds},
            {"realizability_policy_rejects_invalid_values", test_realizability_policy_rejects_invalid_values},
            {"realizability_policy_rejects_invalid_bounded_interval", test_realizability_policy_rejects_invalid_bounded_interval},
            {"realizability_policy_rejects_unknown_entries", test_realizability_policy_rejects_unknown_entries},
            {"qmom_case_parsing_tolerates_realizability_policy_without_installing_it",
             test_qmom_case_parsing_tolerates_realizability_policy_without_installing_it},
            {"case_parsing_validates_realizability_policy_without_owning_it", test_case_parsing_validates_realizability_policy_without_owning_it},
            {"case_parsing_rejects_invalid_realizability_policy", test_case_parsing_rejects_invalid_realizability_policy},
            {"phase_space_evolution_scheme_defaults_and_parses_values", test_phase_space_evolution_scheme_defaults_and_parses_values},
            {"phase_space_evolution_scheme_rejects_invalid_values", test_phase_space_evolution_scheme_rejects_invalid_values},
            {"node_advection_options_defaults_and_parses_schemes", test_node_advection_options_defaults_and_parses_schemes},
            {"node_advection_options_parse_rk23_controls", test_node_advection_options_parse_rk23_controls},
            {"node_advection_options_reject_invalid_scheme_values", test_node_advection_options_reject_invalid_scheme_values},
            {"node_advection_options_reject_invalid_rk23_controls", test_node_advection_options_reject_invalid_rk23_controls},
            {"removed_cqmom_source_policy_section_is_unknown", test_removed_cqmom_source_policy_section_is_unknown},
            {"phase_space_evolution_options_loaders_parse_defaults_explicit_values_and_callbacks",
             test_phase_space_evolution_options_loaders_parse_defaults_explicit_values_and_callbacks},
            {"phase_space_evolution_options_loaders_reject_invalid_values", test_phase_space_evolution_options_loaders_reject_invalid_values},
            {"cqmom_phase_space_evolution_options_loader_result_drives_phase_space_evolution_path",
             test_cqmom_phase_space_evolution_options_loader_result_drives_phase_space_evolution_path},
            {"output_options_default_to_disabled", test_output_options_default_to_disabled},
            {"output_options_parse_enabled_fields", test_output_options_parse_enabled_fields},
            {"output_options_reject_invalid_values", test_output_options_reject_invalid_values},
            {"qmom_case_parsing_tolerates_scheme", test_qmom_case_parsing_tolerates_scheme},
            {"case_parsing_rejects_invalid_runtime_phase_space_evolution_settings",
             test_case_parsing_rejects_invalid_runtime_phase_space_evolution_settings},
            {"qmom_case_parsing_tolerates_output_options_without_installing_them",
             test_qmom_case_parsing_tolerates_output_options_without_installing_them},
            {"duplicate_and_missing_entries", test_duplicate_and_missing_entries},
            {"malformed_nonfinite_invalid_order_and_unknown_key", test_malformed_nonfinite_invalid_order_and_unknown_key},
            {"explicit_moments_take_precedence_over_initial_distribution", test_explicit_moments_take_precedence_over_initial_distribution},
            {"configuration_format_infers_formulation_from_phase_space", test_configuration_format_infers_formulation_from_phase_space},
            {"configuration_rejects_unsupported_or_ambiguous_format_fields", test_configuration_rejects_unsupported_or_ambiguous_format_fields},
            {"file_loading_and_file_failure", test_file_loading_and_file_failure},
            {"indexed_sections_require_canonical_indices", test_indexed_sections_require_canonical_indices},
            {"indexed_entries_that_nothing_reads_are_reported", test_indexed_entries_that_nothing_reads_are_reported},
            {"case_loaders_validate_output_settings", test_case_loaders_validate_output_settings},
            {"realizability_settings_match_the_formulation", test_realizability_settings_match_the_formulation},
            {"oversized_requests_are_reported_not_thrown", test_oversized_requests_are_reported_not_thrown},
            {"numbers_are_read_independently_of_the_locale", test_numbers_are_read_independently_of_the_locale},
            {"errors_name_the_offending_entry", test_errors_name_the_offending_entry},
        });
}
