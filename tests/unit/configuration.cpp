#include "qmomkit/configuration.hpp"

#include "qmomkit/distributions.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <variant>
#include <vector>

namespace
{

using qmomkit::config::ConfigurationErrorCode;

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

template <class T> bool has_error(const qmomkit::config::ConfigurationResult<T>& result, const ConfigurationErrorCode code)
{
    for (const auto& error : result.errors)
    {
        if (error.code == code)
        {
            return true;
        }
    }
    return false;
}

template <class T> bool has_message_fragment(const qmomkit::config::ConfigurationResult<T>& result, const std::string& fragment)
{
    for (const auto& error : result.errors)
    {
        if (error.message.find(fragment) != std::string::npos)
        {
            return true;
        }
    }
    return false;
}

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

void test_valid_qmom_matches_programmatic_case()
{
    const auto configured = qmomkit::config::parse_qmom(qmom_text, "qmom-test");
    assert(configured.ok());
    const qmomkit::QmomCase direct{qmomkit::MomentSet({1.0, 0.5, 0.5, 0.5}), 2, {}};
    assert(configured.value->node_count == direct.node_count);
    assert(configured.value->moments.values() == direct.moments.values());
    assert(!configured.value->options.adaptive);

    const auto configured_result = qmomkit::compute_qmom(*configured.value);
    const auto direct_result = qmomkit::compute_qmom(direct);
    assert(configured_result.ok());
    assert(direct_result.ok());
    assert(configured_result.weights == direct_result.weights);
    assert(configured_result.abscissas == direct_result.abscissas);

    const auto adaptive = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[qmom]\nadaptive=true\n[moments]\nvalues=1,0\n");
    assert(adaptive.ok());
    assert(adaptive.value->options.adaptive);
}

void test_qmom_distribution_initialization_matches_direct_api()
{
    const auto configured = qmomkit::config::parse_qmom(R"(
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
)");
    assert(configured.ok());

    const qmomkit::RawMomentResult direct = qmomkit::compute_raw_moments(qmomkit::LogNormalDistribution{-0.2, 0.4}, 5, 2.5);
    assert(direct.ok());
    assert(configured.value->node_count == 3);
    assert(configured.value->moments.size() == 6);
    assert_near(configured.value->moments[0], 2.5);
    assert_moments_near(configured.value->moments, direct.moments);

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
    const auto moments = qmomkit::config::initialize_moments(*initial.value);
    assert(moments.ok());
    const qmomkit::RawMomentResult expected = qmomkit::compute_raw_moments(qmomkit::GammaDistribution{2.0, 3.0}, 3, 4.0);
    assert(expected.ok());
    assert_moments_near(*moments.value, expected.moments);

    const auto qmom_case = qmomkit::config::create_qmom_case(*initial.value, 2);
    assert(qmom_case.ok());
    assert(qmom_case.value->moments.size() == 4);
    assert_near(qmom_case.value->moments[0], 4.0);
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
        assert(configured.value->moments.size() == 2);
        assert_near(configured.value->moments[0], 2.0);
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
    assert(configured.value->node_count == 4);
    assert(configured.value->moments.size() == 8);
    assert_near(configured.value->moments[0], 2.3);
    assert_moments_near(configured.value->moments, expected.moments);
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
    assert_near(configured.value->moments[0], 3.0);
    assert_near(configured.value->moments[1], 6.0);
    assert_near(configured.value->moments[2], 24.0);
    assert_near(configured.value->moments[3], 144.0);
}

void test_valid_cqmom_matches_programmatic_case()
{
    const auto configured = qmomkit::config::parse_cqmom_dense(cqmom_text, "cqmom-test");
    assert(configured.ok());
    const qmomkit::CqmomDenseCase direct{qmomkit::PhaseSpaceHierarchy({1, 1}, {1, 0}), {{2, 2}, {1.0, 1.5, 0.5, 0.75}}, {}, {true}, true};
    assert(configured.value->hierarchy.node_counts() == direct.hierarchy.node_counts());
    assert(configured.value->hierarchy.conditioning_order() == direct.hierarchy.conditioning_order());
    assert(configured.value->moments.extents == direct.moments.extents);
    assert(configured.value->moments.values == direct.moments.values);
    assert(configured.value->adaptive_branch_counts);
    assert(configured.value->qmom_options.adaptive);

    const auto configured_result = qmomkit::compute_cqmom(*configured.value);
    const auto direct_result = qmomkit::compute_cqmom(direct);
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
    assert(configured.value->hierarchy.node_counts() == std::vector<std::size_t>({2, 2}));
    assert(configured.value->hierarchy.conditioning_order() == std::vector<std::size_t>({1, 0}));
    assert(configured.value->adaptive_branch_counts);
    assert(configured.value->qmom_options.adaptive);
    assert_tensor_near(configured.value->moments, expected.moments);
    assert_near(configured.value->moments.values[qmomkit::flatten_multi_index({0, 0}, configured.value->moments.extents)], 2.5);
    assert(configured.value->moments.extents == std::vector<std::size_t>({4, 4}));

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
    const auto tensor = qmomkit::config::initialize_cqmom_dense_moments(*initial.value);
    const qmomkit::TensorMomentResult expected_tensor = qmomkit::compute_raw_moments(
        qmomkit::IndependentProductDistribution{{qmomkit::BetaDistribution{2.0, 5.0}, qmomkit::LogNormalDistribution{0.1, 0.2}}, 4.0}, {1, 1});
    assert(expected_tensor.ok());
    assert(tensor.ok());
    assert_tensor_near(*tensor.value, expected_tensor.moments);

    const auto dense_case = qmomkit::config::create_cqmom_dense_case(*initial.value, qmomkit::PhaseSpaceHierarchy({1, 1}));
    assert(dense_case.ok());
    assert_tensor_near(dense_case.value->moments, expected_tensor.moments);
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
    const auto dense = qmomkit::config::initialize_cqmom_dense_moments(*initial.value);
    const auto lean = qmomkit::config::initialize_cqmom_moment_data(*initial.value, hierarchy);
    assert(dense.ok());
    assert(lean.ok());

    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments(*dense.value, layout);
    assert_cqmom_moments_near(*lean.value, expected);
    assert_near(lean.value->moment({0, 0}), 2.5);
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
    const auto configured = qmomkit::config::initialize_cqmom_moment_data(*initial.value, hierarchy);
    assert(configured.ok());

    const qmomkit::IndependentProductDistribution distribution{{qmomkit::GammaDistribution{3.0, 0.2}, qmomkit::NormalDistribution{-0.5, 1.25}}, 2.5};
    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(hierarchy);
    const qmomkit::RawMomentRequest moment_request{layout.required_indices};
    const qmomkit::SelectedRawMomentResult selected = qmomkit::compute_raw_moments(distribution, moment_request);
    assert(selected.ok());
    const qmomkit::CqmomMomentData expected = qmomkit::create_cqmom_moment_data(layout, selected);

    assert_cqmom_moments_near(*configured.value, expected);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(qmomkit::CqmomLeanCase{*configured.value});
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
    const auto* independent_product = std::get_if<qmomkit::config::IndependentProductInitialization>(&*initial.value);
    assert(independent_product != nullptr);
    assert(!independent_product->maximum_orders.has_value());

    const qmomkit::PhaseSpaceHierarchy hierarchy({2, 2}, {0, 1});
    const auto initialized = qmomkit::config::initialize_cqmom_moment_data(*initial.value, hierarchy);
    assert(initialized.ok());
    assert(initialized.value->layout.hierarchy.node_counts() == std::vector<std::size_t>({2, 2}));
    assert_near(initialized.value->moment({0, 0}), 2.5);

    const auto lean = qmomkit::config::parse_cqmom_lean(text);
    const auto case_result = qmomkit::config::parse_cqmom(text);
    assert(lean.ok());
    assert(case_result.ok());
    const qmomkit::CqmomResult lean_result = qmomkit::compute_cqmom(*lean.value);
    const qmomkit::CqmomResult case_compute = qmomkit::compute_cqmom(*case_result.value);
    assert(lean_result.ok());
    assert(case_compute.ok());
    assert_cqmom_quadratures_near(case_compute.quadrature, lean_result.quadrature);
}

void test_dense_cqmom_initialization_requires_maximum_orders()
{
    const std::string text = R"(
[phase_space]
node_counts = 1, 1

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

    const auto initial = qmomkit::config::parse_cqmom_initial_condition(text);
    assert(initial.ok());

    const auto tensor = qmomkit::config::initialize_cqmom_dense_moments(*initial.value);
    assert(!tensor.ok());
    assert(has_error(tensor, ConfigurationErrorCode::missing_entry));
    assert(has_message_fragment(tensor, "initial.maximum_orders"));
    assert(has_message_fragment(tensor, "dense/full tensor"));

    const auto dense_case = qmomkit::config::create_cqmom_dense_case(*initial.value, qmomkit::PhaseSpaceHierarchy({1, 1}));
    assert(!dense_case.ok());
    assert(has_error(dense_case, ConfigurationErrorCode::missing_entry));
    assert(has_message_fragment(dense_case, "initial.maximum_orders"));

    const auto parsed_case = qmomkit::config::parse_cqmom_dense(text);
    assert(!parsed_case.ok());
    assert(has_error(parsed_case, ConfigurationErrorCode::missing_entry));
    assert(has_message_fragment(parsed_case, "dense/full tensor"));
}

void test_lean_cqmom_initialization_selects_explicit_dense_tensor()
{
    const auto initial = qmomkit::config::parse_cqmom_initial_condition(cqmom_text, "cqmom-test");
    assert(initial.ok());

    const qmomkit::PhaseSpaceHierarchy hierarchy({1, 1}, {1, 0});
    const auto lean = qmomkit::config::initialize_cqmom_moment_data(*initial.value, hierarchy);
    assert(lean.ok());

    const qmomkit::RawMomentTensor dense{{2, 2}, {1.0, 1.5, 0.5, 0.75}};
    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments(dense, qmomkit::create_cqmom_moment_layout(hierarchy));
    assert_cqmom_moments_near(*lean.value, expected);
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
    assert(lean.value->moments.layout.hierarchy.node_counts() == std::vector<std::size_t>({2, 2}));
    assert(lean.value->moments.layout.hierarchy.conditioning_order() == std::vector<std::size_t>({0, 1}));
    assert(lean.value->diagnostic_level == qmomkit::CqmomDiagnosticLevel::basic);
    assert_near(lean.value->moments.moment({0, 0}), 2.5);

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
    assert(has_error(dense, ConfigurationErrorCode::invalid_configuration));

    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(*lean.value);
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

    const qmomkit::CqmomMomentLayout layout = qmomkit::create_cqmom_moment_layout(dense.value->hierarchy);
    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments(dense.value->moments, layout);
    assert_cqmom_moments_near(lean.value->moments, expected);

    const qmomkit::CqmomResult dense_result = qmomkit::compute_cqmom(*dense.value);
    const qmomkit::CqmomResult lean_result = qmomkit::compute_cqmom(*lean.value);
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
    assert(lean.value->adaptive_branch_counts);
    assert(lean.value->qmom_options.adaptive);

    const qmomkit::PhaseSpaceHierarchy hierarchy({1, 1}, {1, 0});
    const qmomkit::RawMomentTensor dense{{2, 2}, {1.0, 1.5, 0.5, 0.75}};
    const qmomkit::CqmomMomentData expected = qmomkit::extract_cqmom_moments(dense, qmomkit::create_cqmom_moment_layout(hierarchy));
    assert_cqmom_moments_near(lean.value->moments, expected);

    const auto initial = qmomkit::config::parse_cqmom_initial_condition(cqmom_text, "cqmom-test");
    assert(initial.ok());
    const auto made = qmomkit::config::create_cqmom_lean_case(*initial.value, hierarchy, {true}, true);
    assert(made.ok());
    assert_cqmom_moments_near(made.value->moments, expected);
    assert(made.value->adaptive_branch_counts);
    assert(made.value->qmom_options.adaptive);
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

    const qmomkit::CqmomCase cqmom_case = *case_result.value;
    assert(cqmom_case.hierarchy().node_counts() == std::vector<std::size_t>({2, 2}));
    assert(cqmom_case.hierarchy().conditioning_order() == std::vector<std::size_t>({1, 0}));
    assert(cqmom_case.diagnostic_level() == qmomkit::CqmomDiagnosticLevel::basic);

    const qmomkit::CqmomResult case_compute = qmomkit::compute_cqmom(cqmom_case);
    const qmomkit::CqmomResult lean_compute = qmomkit::compute_cqmom(*lean.value);
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

    const qmomkit::CqmomCase cqmom_case = *case_result.value;
    assert(cqmom_case.hierarchy().node_counts() == std::vector<std::size_t>({1, 1}));
    assert(cqmom_case.hierarchy().conditioning_order() == std::vector<std::size_t>({1, 0}));
    assert(cqmom_case.adaptive_branch_counts());
    assert(cqmom_case.qmom_options().adaptive);

    const qmomkit::CqmomResult case_compute = qmomkit::compute_cqmom(cqmom_case);
    const qmomkit::CqmomResult dense_compute = qmomkit::compute_cqmom(*dense.value);
    const qmomkit::CqmomResult lean_compute = qmomkit::compute_cqmom(*lean.value);
    assert(case_compute.ok());
    assert(dense_compute.ok());
    assert(lean_compute.ok());
    assert_cqmom_quadratures_near(case_compute.quadrature, dense_compute.quadrature);
    assert_cqmom_quadratures_near(case_compute.quadrature, lean_compute.quadrature);
}

void test_cqmom_case_initial_state_preserves_case_settings()
{
    const auto case_result = qmomkit::config::parse_cqmom(cqmom_text, "cqmom-test");
    assert(case_result.ok());

    const qmomkit::CqmomCase cqmom_case = *case_result.value;
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

    const qmomkit::CqmomState copied = state;
    assert(copied.hierarchy().node_counts() == state.hierarchy().node_counts());
    assert(copied.hierarchy().conditioning_order() == state.hierarchy().conditioning_order());
    assert(copied.qmom_options().adaptive == state.qmom_options().adaptive);
    assert(copied.vandermonde_solver() == state.vandermonde_solver());
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
    assert(loaded.value->qmom_options().adaptive);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(*loaded.value);
    assert(result.ok());
    std::filesystem::remove(filename);

    const auto missing = qmomkit::config::load_cqmom(filename);
    assert(!missing.ok());
    assert(has_error(missing, ConfigurationErrorCode::file_error));
    assert(missing.errors.front().source_name == filename.string());

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
    const qmomkit::CqmomResult lean_result = qmomkit::compute_cqmom(*lean_loaded.value);
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
    assert(configured.value->moments.extents == std::vector<std::size_t>({2, 2}));
    for (const double value : configured.value->moments.values)
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
    assert(configured.value->qmom_options.adaptive == adaptive);
    assert(configured.value->adaptive_branch_counts == adaptive_branch_counts);
}

void test_cqmom_adaptive_options_are_independent()
{
    assert_cqmom_adaptive_options(false, false);
    assert_cqmom_adaptive_options(true, false);
    assert_cqmom_adaptive_options(false, true);
    assert_cqmom_adaptive_options(true, true);

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
    assert(!omitted_qmom_adaptive.value->qmom_options.adaptive);
    assert(omitted_qmom_adaptive.value->adaptive_branch_counts);
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
    assert(omitted.value->diagnostic_level == qmomkit::CqmomDiagnosticLevel::detailed);

    const auto detailed = qmomkit::config::parse_cqmom_dense(cqmom_diagnostic_level_text("detailed"));
    assert(detailed.ok());
    assert(detailed.value->diagnostic_level == qmomkit::CqmomDiagnosticLevel::detailed);

    const auto basic = qmomkit::config::parse_cqmom_dense(cqmom_diagnostic_level_text("basic"));
    assert(basic.ok());
    assert(basic.value->diagnostic_level == qmomkit::CqmomDiagnosticLevel::basic);

    const auto none = qmomkit::config::parse_cqmom_dense(cqmom_diagnostic_level_text("none"));
    assert(none.ok());
    assert(none.value->diagnostic_level == qmomkit::CqmomDiagnosticLevel::none);

    const auto invalid = qmomkit::config::parse_cqmom_dense(cqmom_diagnostic_level_text("verbose"));
    assert(!invalid.ok());
    assert(has_error(invalid, ConfigurationErrorCode::invalid_value));
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
    assert(omitted.value->vandermonde_solver == qmomkit::VandermondeSolver::automatic);

    const auto dense = qmomkit::config::parse_cqmom_dense(cqmom_vandermonde_solver_text("dense"));
    assert(dense.ok());
    assert(dense.value->vandermonde_solver == qmomkit::VandermondeSolver::dense);

    const auto bp = qmomkit::config::parse_cqmom_lean(cqmom_vandermonde_solver_text("bjorck_pereyra"));
    assert(bp.ok());
    assert(bp.value->vandermonde_solver == qmomkit::VandermondeSolver::bjorck_pereyra);

    const auto fallback = qmomkit::config::parse_cqmom_dense(cqmom_vandermonde_solver_text("bjorck_pereyra_with_dense_fallback"));
    assert(fallback.ok());
    assert(fallback.value->vandermonde_solver == qmomkit::VandermondeSolver::bjorck_pereyra_with_dense_fallback);

    const auto automatic = qmomkit::config::parse_cqmom_lean(cqmom_vandermonde_solver_text("automatic"));
    assert(automatic.ok());
    assert(automatic.value->vandermonde_solver == qmomkit::VandermondeSolver::automatic);

    const auto invalid = qmomkit::config::parse_cqmom_dense(cqmom_vandermonde_solver_text("structured"));
    assert(!invalid.ok());
    assert(has_error(invalid, ConfigurationErrorCode::invalid_value));
    assert(has_message_fragment(invalid, "vandermonde.solver"));
}

void assert_cqmom_case_reports_vandermonde_solver(const std::string& configured_solver, const qmomkit::VandermondeSolver expected_solver)
{
    const auto parsed = qmomkit::config::parse_cqmom(cqmom_vandermonde_solver_text(configured_solver));
    assert(parsed.ok());

    const qmomkit::CqmomCase cqmom_case = *parsed.value;
    assert(cqmom_case.vandermonde_solver() == expected_solver);

    const qmomkit::CqmomState state = cqmom_case.initial_state();
    assert(state.vandermonde_solver() == expected_solver);
    assert(state.vandermonde_solver() == cqmom_case.vandermonde_solver());
}

void test_cqmom_case_and_state_report_vandermonde_solver()
{
    assert_cqmom_case_reports_vandermonde_solver("", qmomkit::VandermondeSolver::automatic);
    assert_cqmom_case_reports_vandermonde_solver("dense", qmomkit::VandermondeSolver::dense);
    assert_cqmom_case_reports_vandermonde_solver("bjorck_pereyra", qmomkit::VandermondeSolver::bjorck_pereyra);
    assert_cqmom_case_reports_vandermonde_solver("bjorck_pereyra_with_dense_fallback",
                                                 qmomkit::VandermondeSolver::bjorck_pereyra_with_dense_fallback);
}

void test_realizability_policy_defaults_to_disabled()
{
    const auto policy = qmomkit::config::parse_realizability_policy(qmom_text);
    assert(policy.ok());
    assert(policy.value->check_mode == qmomkit::RealizabilityCheckMode::None);

    const auto empty_section = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
)");
    assert(empty_section.ok());
    assert(empty_section.value->check_mode == qmomkit::RealizabilityCheckMode::None);
}

void test_realizability_policy_parses_check_modes()
{
    const auto none = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = none
)");
    assert(none.ok());
    assert(none.value->check_mode == qmomkit::RealizabilityCheckMode::None);

    const auto every_stage = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = every_stage
support = nonnegative
)");
    assert(every_stage.ok());
    assert(every_stage.value->check_mode == qmomkit::RealizabilityCheckMode::EveryStage);
    assert(every_stage.value->support == qmomkit::MomentSupport::nonnegative);

    const auto final_only = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = nonnegative
accept_boundary = false
tolerance = 1e-10
relative_tolerance = false
)");
    assert(final_only.ok());
    assert(final_only.value->check_mode == qmomkit::RealizabilityCheckMode::FinalOnly);
    assert(final_only.value->support == qmomkit::MomentSupport::nonnegative);
    assert(!final_only.value->accept_boundary);
    assert_near(final_only.value->options.tolerance, 1.0e-10);
    assert(!final_only.value->options.relative_tolerance);
}

void test_realizability_policy_parses_support_and_bounded_interval()
{
    const auto real_line = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = real_line
)");
    assert(real_line.ok());
    assert(real_line.value->support == qmomkit::MomentSupport::real_line);

    const auto bounded = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = bounded_interval
lower = -1.0
upper = 2.0
)");
    assert(bounded.ok());
    assert(bounded.value->support == qmomkit::MomentSupport::bounded_interval);
    assert(bounded.value->interval.has_value());
    assert_near(bounded.value->interval->lower, -1.0);
    assert_near(bounded.value->interval->upper, 2.0);
}

void test_realizability_policy_rejects_invalid_values()
{
    const auto unknown_check_mode = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = accepted_step
support = nonnegative
)");
    assert(!unknown_check_mode.ok());
    assert(has_error(unknown_check_mode, ConfigurationErrorCode::invalid_value));
    assert(has_message_fragment(unknown_check_mode, "check_mode"));

    const auto missing_support = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
)");
    assert(!missing_support.ok());
    assert(has_error(missing_support, ConfigurationErrorCode::missing_entry));

    const auto unknown_support = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = positive
)");
    assert(!unknown_support.ok());
    assert(has_error(unknown_support, ConfigurationErrorCode::invalid_value));
    assert(has_message_fragment(unknown_support, "support"));

    const auto invalid_tolerance = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = nonnegative
tolerance = 0.0
)");
    assert(!invalid_tolerance.ok());
    assert(has_error(invalid_tolerance, ConfigurationErrorCode::invalid_value));

    const auto nonfinite_tolerance = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = nonnegative
tolerance = nan
)");
    assert(!nonfinite_tolerance.ok());
    assert(has_error(nonfinite_tolerance, ConfigurationErrorCode::invalid_value));

    const auto malformed_bool = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = nonnegative
accept_boundary = yes
)");
    assert(!malformed_bool.ok());
    assert(has_error(malformed_bool, ConfigurationErrorCode::invalid_value));
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
    assert(has_error(missing_bound, ConfigurationErrorCode::missing_entry));

    const auto invalid_interval = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = bounded_interval
lower = 2.0
upper = 1.0
)");
    assert(!invalid_interval.ok());
    assert(has_error(invalid_interval, ConfigurationErrorCode::invalid_value));

    const auto nonfinite_bound = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = bounded_interval
lower = -inf
upper = 1.0
)");
    assert(!nonfinite_bound.ok());
    assert(has_error(nonfinite_bound, ConfigurationErrorCode::invalid_value));
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
    assert(has_error(unknown_key, ConfigurationErrorCode::unknown_entry));

    const auto unknown_section = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability.extra]
check_mode = final_only
)");
    assert(!unknown_section.ok());
    assert(has_error(unknown_section, ConfigurationErrorCode::unknown_entry));
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
    assert(qmom_case.value->moments.values() == std::vector<double>({1.0, 0.0}));

    const auto policy = qmomkit::config::parse_realizability_policy(R"(
[time_integration.realizability]
check_mode = final_only
support = real_line
)");
    assert(policy.ok());
    qmomkit::SspRk3Options options;
    options.accept =
        qmomkit::make_realizability_acceptance_predicate(*policy.value, [](const std::vector<double>& state) { return qmomkit::MomentSet(state); });
    assert(options.accept);
    assert(options.accept({1.0, 0.0, 1.0}));
}

void test_scheme_defaults_and_parses_values()
{
    const auto omitted = qmomkit::config::parse_transport_scheme(qmom_text);
    assert(omitted.ok());
    assert(*omitted.value == qmomkit::TransportScheme::strang);

    const auto strang = qmomkit::config::parse_transport_scheme(R"(
[time_integration]
scheme = strang
)");
    assert(strang.ok());
    assert(*strang.value == qmomkit::TransportScheme::strang);

    const auto cd = qmomkit::config::parse_transport_scheme(R"(
[time_integration]
scheme = cd
)");
    assert(cd.ok());
    assert(*cd.value == qmomkit::TransportScheme::cd);

    const auto dc = qmomkit::config::parse_transport_scheme(R"(
[time_integration]
scheme = dc
)");
    assert(dc.ok());
    assert(*dc.value == qmomkit::TransportScheme::dc);
}

void test_scheme_rejects_invalid_values()
{
    const auto invalid = qmomkit::config::parse_transport_scheme(R"(
[time_integration]
scheme = leapfrog
)");
    assert(!invalid.ok());
    assert(has_error(invalid, ConfigurationErrorCode::invalid_value));
    assert(has_message_fragment(invalid, "strang, cd, or dc"));

    const auto unknown_time_key = qmomkit::config::parse_transport_scheme(R"(
[time_integration]
stage_count = 3
)");
    assert(!unknown_time_key.ok());
    assert(has_error(unknown_time_key, ConfigurationErrorCode::unknown_entry));
}

void test_node_advection_scheme_defaults_and_parses_values()
{
    qmomkit::NodeAdvectionOptions default_options;
    assert(default_options.scheme == qmomkit::NodeAdvectionScheme::ssp_rk3);

    const auto omitted = qmomkit::config::parse_node_advection_scheme(qmom_text);
    assert(omitted.ok());
    assert(*omitted.value == qmomkit::NodeAdvectionScheme::ssp_rk3);

    const auto ssp_rk3 = qmomkit::config::parse_node_advection_scheme(R"(
[time_integration.node_advection]
scheme = ssp_rk3
)");
    assert(ssp_rk3.ok());
    assert(*ssp_rk3.value == qmomkit::NodeAdvectionScheme::ssp_rk3);

    const auto adaptive_rk23 = qmomkit::config::parse_node_advection_scheme(R"(
[time_integration.node_advection]
scheme = adaptive_rk23
)");
    assert(adaptive_rk23.ok());
    assert(*adaptive_rk23.value == qmomkit::NodeAdvectionScheme::adaptive_rk23);

    const auto alongside_transport_scheme = qmomkit::config::parse_node_advection_scheme(R"(
[time_integration]
scheme = cd

[time_integration.node_advection]
scheme = adaptive_rk23
)");
    assert(alongside_transport_scheme.ok());
    assert(*alongside_transport_scheme.value == qmomkit::NodeAdvectionScheme::adaptive_rk23);
}

void test_node_advection_scheme_rejects_invalid_values()
{
    const auto invalid = qmomkit::config::parse_node_advection_scheme(R"(
[time_integration.node_advection]
scheme = rk4
)");
    assert(!invalid.ok());
    assert(has_error(invalid, ConfigurationErrorCode::invalid_value));
    assert(has_message_fragment(invalid, "ssp_rk3 or adaptive_rk23"));

    const auto unknown_time_key = qmomkit::config::parse_node_advection_scheme(R"(
[time_integration.node_advection]
stage_count = 3
)");
    assert(!unknown_time_key.ok());
    assert(has_error(unknown_time_key, ConfigurationErrorCode::unknown_entry));
}

void test_output_options_default_to_disabled()
{
    const auto options = qmomkit::config::parse_output_options(qmom_text);
    assert(options.ok());
    assert(options.value->directory == std::filesystem::path("."));
    assert(!options.value->write_moments);
    assert(!options.value->write_quadrature);
    assert(!options.value->interval.has_value());

    const auto empty_section = qmomkit::config::parse_output_options(R"(
[output]
)");
    assert(empty_section.ok());
    assert(empty_section.value->directory == std::filesystem::path("."));
    assert(!empty_section.value->write_moments);
    assert(!empty_section.value->write_quadrature);
    assert(!empty_section.value->interval.has_value());
}

void test_output_options_parse_enabled_fields()
{
    const auto moments = qmomkit::config::parse_output_options(R"(
[output]
write_moments = true
)");
    assert(moments.ok());
    assert(moments.value->write_moments);
    assert(!moments.value->write_quadrature);

    const auto quadrature = qmomkit::config::parse_output_options(R"(
[output]
write_quadrature = true
)");
    assert(quadrature.ok());
    assert(!quadrature.value->write_moments);
    assert(quadrature.value->write_quadrature);

    const auto all_fields = qmomkit::config::parse_output_options(R"(
[output]
directory = output
write_moments = true
write_quadrature = false
interval = 0.1
)");
    assert(all_fields.ok());
    assert(all_fields.value->directory == std::filesystem::path("output"));
    assert(all_fields.value->write_moments);
    assert(!all_fields.value->write_quadrature);
    assert(all_fields.value->interval.has_value());
    assert_near(*all_fields.value->interval, 0.1);
}

void test_output_options_reject_invalid_values()
{
    const auto zero_interval = qmomkit::config::parse_output_options(R"(
[output]
interval = 0.0
)");
    assert(!zero_interval.ok());
    assert(has_error(zero_interval, ConfigurationErrorCode::invalid_value));
    assert(has_message_fragment(zero_interval, "interval"));

    const auto nonfinite_interval = qmomkit::config::parse_output_options(R"(
[output]
interval = nan
)");
    assert(!nonfinite_interval.ok());
    assert(has_error(nonfinite_interval, ConfigurationErrorCode::invalid_value));

    const auto malformed_bool = qmomkit::config::parse_output_options(R"(
[output]
write_moments = yes
)");
    assert(!malformed_bool.ok());
    assert(has_error(malformed_bool, ConfigurationErrorCode::invalid_value));
    assert(has_message_fragment(malformed_bool, "write_moments"));

    const auto unknown_key = qmomkit::config::parse_output_options(R"(
[output]
write_nodes = true
)");
    assert(!unknown_key.ok());
    assert(has_error(unknown_key, ConfigurationErrorCode::unknown_entry));

    const auto unknown_section = qmomkit::config::parse_output_options(R"(
[output.moments]
M[0] = true
)");
    assert(!unknown_section.ok());
    assert(has_error(unknown_section, ConfigurationErrorCode::unknown_entry));

    const auto empty_directory = qmomkit::config::parse_output_options(R"(
[output]
directory =
)");
    assert(!empty_directory.ok());
    assert(has_error(empty_directory, ConfigurationErrorCode::syntax_error));
    assert(has_message_fragment(empty_directory, "must not be empty"));

    const auto invalid_syntax = qmomkit::config::parse_output_options(R"(
[output]
directory
)");
    assert(!invalid_syntax.ok());
    assert(has_error(invalid_syntax, ConfigurationErrorCode::syntax_error));
}

void test_qmom_case_parsing_tolerates_scheme()
{
    const auto qmom = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1

[time_integration]
scheme = dc

[time_integration.node_advection]
scheme = ssp_rk3

[moments]
values = 1.0, 0.0
)");
    assert(qmom.ok());

    const auto cqmom = qmomkit::config::parse_cqmom_dense(R"(
[phase_space]
node_counts = 1, 1

[time_integration]
scheme = cd

[time_integration.node_advection]
scheme = adaptive_rk23

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
    assert(has_error(invalid, ConfigurationErrorCode::invalid_value));

    const auto invalid_node_advection = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 1

[time_integration.node_advection]
scheme = rk4

[moments]
values = 1.0, 0.0
)");
    assert(!invalid_node_advection.ok());
    assert(has_error(invalid_node_advection, ConfigurationErrorCode::invalid_value));
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
    assert(qmom.value->moments.values() == std::vector<double>({1.0, 0.0}));

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
    assert(options.value->write_moments);
    assert(options.value->write_quadrature);
}

void test_duplicate_and_missing_entries()
{
    const auto duplicate_key =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[qmom]\nadaptive=false\nadaptive=false\n[moments]\nvalues=1,0\n");
    assert(!duplicate_key.ok());
    assert(has_error(duplicate_key, ConfigurationErrorCode::duplicate_entry));

    const auto duplicate_moment = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[moments]\nextents=2,2\nM[0,0]=1\nM[0,0]=1\nM[0,1]=0\nM[1,0]=0\nM[1,1]=0\n");
    assert(!duplicate_moment.ok());
    assert(has_error(duplicate_moment, ConfigurationErrorCode::duplicate_entry));

    const auto missing_moment = qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1,1\n[moments]\nextents=2,2\nM[0,0]=1\n");
    assert(!missing_moment.ok());
    assert(has_error(missing_moment, ConfigurationErrorCode::missing_entry));

    const auto missing_initial_type = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\nfamily=normal\n");
    assert(!missing_initial_type.ok());
    assert(has_error(missing_initial_type, ConfigurationErrorCode::missing_entry));

    const auto missing_family = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nmaximum_order=1\n");
    assert(!missing_family.ok());
    assert(has_error(missing_family, ConfigurationErrorCode::missing_entry));

    const auto missing_parameter =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily=normal\nmaximum_order=1\nmean=0\n");
    assert(!missing_parameter.ok());
    assert(has_error(missing_parameter, ConfigurationErrorCode::missing_entry));

    const auto missing_maximum_order =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!missing_maximum_order.ok());
    assert(has_error(missing_maximum_order, ConfigurationErrorCode::missing_entry));

    const auto missing_mixture_components = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n");
    assert(!missing_mixture_components.ok());
    assert(has_error(missing_mixture_components, ConfigurationErrorCode::missing_entry));

    const auto missing_component_family = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                                                      "[initial.mixture.component.0]\nnumber_density=1\n");
    assert(!missing_component_family.ok());
    assert(has_error(missing_component_family, ConfigurationErrorCode::missing_entry));

    const auto missing_component_density = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                                                       "[initial.mixture.component.0]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!missing_component_density.ok());
    assert(has_error(missing_component_density, ConfigurationErrorCode::missing_entry));

    const auto missing_cqmom_initial_type = qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1,1\n[initial]\nnumber_density=1\n");
    assert(!missing_cqmom_initial_type.ok());
    assert(has_error(missing_cqmom_initial_type, ConfigurationErrorCode::missing_entry));

    const auto missing_product_dimensions = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n");
    assert(!missing_product_dimensions.ok());
    assert(has_error(missing_product_dimensions, ConfigurationErrorCode::missing_entry));

    const auto missing_product_family = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nmean=0\nstandard_deviation=1\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!missing_product_family.ok());
    assert(has_error(missing_product_family, ConfigurationErrorCode::missing_entry));
}

void test_malformed_nonfinite_invalid_order_and_unknown_key()
{
    const auto malformed = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=two\n[moments]\nvalues=1,0\n", "malformed");
    assert(!malformed.ok());
    assert(has_error(malformed, ConfigurationErrorCode::invalid_value));
    assert(!malformed.errors.front().source_name.empty());
    assert(malformed.errors.front().line > 0);

    const auto nonfinite = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[moments]\nvalues=1,nan\n");
    assert(!nonfinite.ok());
    assert(has_error(nonfinite, ConfigurationErrorCode::invalid_value));

    const auto malformed_adaptive = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[qmom]\nadaptive=yes\n[moments]\nvalues=1,0\n");
    assert(!malformed_adaptive.ok());
    assert(has_error(malformed_adaptive, ConfigurationErrorCode::invalid_value));

    const auto invalid_order = qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1,1\nconditioning_order=0,0\n[moments]\nextents=2,2\n"
                                                                  "M[0,0]=1\nM[0,1]=1\nM[1,0]=1\nM[1,1]=1\n");
    assert(!invalid_order.ok());
    assert(has_error(invalid_order, ConfigurationErrorCode::invalid_configuration));

    const auto cqmom_initial_section = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[moments]\nextents=2,2\nM[0,0]=1\nM[0,1]=0\nM[1,0]=0\nM[1,1]=0\n[initial]\n");
    assert(cqmom_initial_section.ok());

    const auto unknown = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\nmystery=4\n[moments]\nvalues=1,0\n");
    assert(!unknown.ok());
    assert(has_error(unknown, ConfigurationErrorCode::unknown_entry));

    const auto unsupported_initial_type = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=independent_product\n");
    assert(!unsupported_initial_type.ok());
    assert(has_error(unsupported_initial_type, ConfigurationErrorCode::invalid_value));

    const auto unsupported_family =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily=poisson\nmaximum_order=1\n");
    assert(!unsupported_family.ok());
    assert(has_error(unsupported_family, ConfigurationErrorCode::invalid_value));

    const auto invalid_density = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily="
                                                             "normal\nnumber_density=-1\nmaximum_order=1\nmean=0\nstandard_deviation=1\n");
    assert(!invalid_density.ok());
    assert(has_error(invalid_density, ConfigurationErrorCode::invalid_value));

    const auto nonfinite_density = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily="
                                                               "normal\nnumber_density=nan\nmaximum_order=1\nmean=0\nstandard_deviation=1\n");
    assert(!nonfinite_density.ok());
    assert(has_error(nonfinite_density, ConfigurationErrorCode::invalid_value));

    const auto invalid_maximum_order = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily="
                                                                   "normal\nmaximum_order=one\nmean=0\nstandard_deviation=1\n");
    assert(!invalid_maximum_order.ok());
    assert(has_error(invalid_maximum_order, ConfigurationErrorCode::invalid_value));

    const auto invalid_mixture_maximum_order =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=one\n"
                                    "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=1\nshape=1\nscale=2\n");
    assert(!invalid_mixture_maximum_order.ok());
    assert(has_error(invalid_mixture_maximum_order, ConfigurationErrorCode::invalid_value));

    const auto invalid_distribution_parameter = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily="
                                                                            "normal\nmaximum_order=1\nmean=0\nstandard_deviation=0\n");
    assert(!invalid_distribution_parameter.ok());
    assert(has_error(invalid_distribution_parameter, ConfigurationErrorCode::invalid_value));

    const auto invalid_component_density =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                    "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=-1\nshape=1\nscale=2\n");
    assert(!invalid_component_density.ok());
    assert(has_error(invalid_component_density, ConfigurationErrorCode::invalid_value));

    const auto nonfinite_component_density =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                    "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=nan\nshape=1\nscale=2\n");
    assert(!nonfinite_component_density.ok());
    assert(has_error(nonfinite_component_density, ConfigurationErrorCode::invalid_value));

    const auto invalid_component_parameter =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                    "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=1\nshape=0\nscale=2\n");
    assert(!invalid_component_parameter.ok());
    assert(has_error(invalid_component_parameter, ConfigurationErrorCode::invalid_value));
    assert(has_message_fragment(invalid_component_parameter, "Component index: 0"));

    const auto invalid_second_component_parameter =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                    "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=1\nshape=1\nscale=2\n"
                                    "[initial.mixture.component.1]\nfamily=normal\nnumber_density=1\nmean=0\nstandard_deviation=0\n");
    assert(!invalid_second_component_parameter.ok());
    assert(has_error(invalid_second_component_parameter, ConfigurationErrorCode::invalid_value));
    assert(has_message_fragment(invalid_second_component_parameter, "Component index: 1"));

    const auto all_zero_mixture = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                                              "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=0\nshape=1\nscale=2\n");
    assert(!all_zero_mixture.ok());
    assert(has_error(all_zero_mixture, ConfigurationErrorCode::invalid_value));

    const auto noncontiguous_mixture =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                    "[initial.mixture.component.1]\nfamily=gamma\nnumber_density=1\nshape=1\nscale=2\n");
    assert(!noncontiguous_mixture.ok());
    assert(has_error(noncontiguous_mixture, ConfigurationErrorCode::invalid_configuration));

    const auto unknown_component_key =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=mixture\nmaximum_order=1\n"
                                    "[initial.mixture.component.0]\nfamily=gamma\nnumber_density=1\nshape=1\nscale=2\nweight_fraction=0.5\n");
    assert(!unknown_component_key.ok());
    assert(has_error(unknown_component_key, ConfigurationErrorCode::unknown_entry));

    const auto component_with_distribution_type =
        qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1\n[initial]\ntype=distribution\nfamily=normal\nmaximum_order=1\nmean=0\n"
                                    "standard_deviation=1\n[initial.mixture.component.0]\nfamily=gamma\nnumber_density=1\nshape=1\nscale=2\n");
    assert(!component_with_distribution_type.ok());
    assert(has_error(component_with_distribution_type, ConfigurationErrorCode::invalid_configuration));

    const auto unsupported_cqmom_initial_type =
        qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1,1\n[initial]\ntype=distribution\nmaximum_order=1\nfamily=normal\nmean=0\n"
                                           "standard_deviation=1\n");
    assert(!unsupported_cqmom_initial_type.ok());
    assert(has_error(unsupported_cqmom_initial_type, ConfigurationErrorCode::unknown_entry) ||
           has_error(unsupported_cqmom_initial_type, ConfigurationErrorCode::invalid_value));

    const auto invalid_product_order = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=one\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!invalid_product_order.ok());
    assert(has_error(invalid_product_order, ConfigurationErrorCode::invalid_value));

    const auto product_tensor_size_overflow =
        qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\n"
                                           "maximum_orders=4294967295,4294967295\n"
                                           "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
                                           "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!product_tensor_size_overflow.ok());
    assert(has_error(product_tensor_size_overflow, ConfigurationErrorCode::invalid_configuration));

    const auto invalid_product_density = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=-1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!invalid_product_density.ok());
    assert(has_error(invalid_product_density, ConfigurationErrorCode::invalid_value));

    const auto noncontiguous_product_dimensions = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!noncontiguous_product_dimensions.ok());
    assert(has_error(noncontiguous_product_dimensions, ConfigurationErrorCode::invalid_configuration));

    const auto extra_product_dimension = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
        "[initial.dimension.1]\nfamily=normal\nmean=0\nstandard_deviation=1\n"
        "[initial.dimension.2]\nfamily=beta\nalpha=1\nbeta=2\n");
    assert(!extra_product_dimension.ok());
    assert(has_error(extra_product_dimension, ConfigurationErrorCode::invalid_configuration));

    const auto duplicate_product_dimension = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!duplicate_product_dimension.ok());
    assert(has_error(duplicate_product_dimension, ConfigurationErrorCode::duplicate_entry));

    const auto invalid_marginal_parameter = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=0\nscale=2\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!invalid_marginal_parameter.ok());
    assert(has_error(invalid_marginal_parameter, ConfigurationErrorCode::invalid_value));

    const auto product_dimension_mismatch =
        qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1\n"
                                           "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!product_dimension_mismatch.ok());
    assert(has_error(product_dimension_mismatch, ConfigurationErrorCode::invalid_configuration));

    const auto product_moment_extent_too_small = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=2,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nfamily=gamma\nshape=1\nscale=2\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!product_moment_extent_too_small.ok());
    assert(has_error(product_moment_extent_too_small, ConfigurationErrorCode::invalid_configuration));

    const auto product_rejects_marginal_density = qmomkit::config::parse_cqmom_dense(
        "[phase_space]\nnode_counts=1,1\n[initial]\ntype=independent_product\nnumber_density=1\nmaximum_orders=1,1\n"
        "[initial.dimension.0]\nfamily=gamma\nnumber_density=2\nshape=1\nscale=2\n"
        "[initial.dimension.1]\nfamily=gamma\nshape=1\nscale=2\n");
    assert(!product_rejects_marginal_density.ok());
    assert(has_error(product_rejects_marginal_density, ConfigurationErrorCode::unknown_entry));
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
    assert(configured.value->moments.values() == std::vector<double>({3.0, 1.0, 2.0, 4.0}));

    const auto initial = qmomkit::config::parse_qmom_initial_condition(R"(
[phase_space]
node_counts = 2

[moments]
values = 3.0, 1.0, 2.0, 4.0

[initial]
type = mixture
)");
    assert(initial.ok());
    const auto moments = qmomkit::config::initialize_moments(*initial.value);
    assert(moments.ok());
    assert(moments.value->values() == std::vector<double>({3.0, 1.0, 2.0, 4.0}));

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
    assert(mixture.value->moments.values() == std::vector<double>({3.0, 1.0, 2.0, 4.0}));

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
    assert(unknown_initial_key.value->moments.values() == std::vector<double>({3.0, 1.0, 2.0, 4.0}));

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
    assert(malformed_ignored_component_section.value->moments.values() == std::vector<double>({3.0, 1.0, 2.0, 4.0}));

    const auto unrelated_unknown_section = qmomkit::config::parse_qmom(R"(
[phase_space]
node_counts = 2

[moments]
values = 3.0, 1.0, 2.0, 4.0

[unrelated]
value = 1
)");
    assert(!unrelated_unknown_section.ok());
    assert(has_error(unrelated_unknown_section, ConfigurationErrorCode::unknown_entry));

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
    assert(cqmom.value->moments.values == std::vector<double>({3.0, 1.0, 2.0, 4.0}));
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
    assert(qmom.value->node_count == 3);

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
    assert(cqmom.value->moments.layout.hierarchy.node_counts() == std::vector<std::size_t>({1, 1, 1}));
    assert(cqmom.value->moments.layout.hierarchy.conditioning_order() == std::vector<std::size_t>({0, 1, 2}));

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
    assert(explicit_order.value->moments.layout.hierarchy.conditioning_order() == std::vector<std::size_t>({2, 0, 1}));
}

void test_configuration_rejects_unsupported_or_ambiguous_format_fields()
{
    const auto unsupported_version = qmomkit::config::parse_qmom("format_version=2\n[phase_space]\nnode_counts=1\n[moments]\nvalues=1,0\n");
    assert(!unsupported_version.ok());
    assert(has_error(unsupported_version, ConfigurationErrorCode::invalid_value));
    assert(has_message_fragment(unsupported_version, "format_version must be 1"));

    const auto qmom_with_multivariate_phase_space = qmomkit::config::parse_qmom("[phase_space]\nnode_counts=1,1\n[moments]\nvalues=1,0\n");
    assert(!qmom_with_multivariate_phase_space.ok());
    assert(has_error(qmom_with_multivariate_phase_space, ConfigurationErrorCode::invalid_value));
    assert(has_message_fragment(qmom_with_multivariate_phase_space, "exactly one positive integer"));

    const auto cqmom_with_univariate_phase_space =
        qmomkit::config::parse_cqmom_dense("[phase_space]\nnode_counts=1\n[moments]\nextents=2\nM[0]=1\nM[1]=0\n");
    assert(!cqmom_with_univariate_phase_space.ok());
    assert(has_error(cqmom_with_univariate_phase_space, ConfigurationErrorCode::invalid_configuration));
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
    assert(has_error(missing, ConfigurationErrorCode::file_error));
    assert(missing.errors.front().source_name == filename.string());

    {
        std::ofstream output(filename);
        output << cqmom_text;
    }
    const auto loaded_cqmom = qmomkit::config::load_cqmom_lean(filename);
    assert(loaded_cqmom.ok());
    assert(loaded_cqmom.value->moments.layout.hierarchy.conditioning_order() == std::vector<std::size_t>({1, 0}));
    assert_near(loaded_cqmom.value->moments.moment({0, 0}), 1.0);
    std::filesystem::remove(filename);

    const auto missing_cqmom = qmomkit::config::load_cqmom_lean(filename);
    assert(!missing_cqmom.ok());
    assert(has_error(missing_cqmom, ConfigurationErrorCode::file_error));
    assert(missing_cqmom.errors.front().source_name == filename.string());

    {
        std::ofstream output(filename);
        output << "[output]\ndirectory = output\nwrite_moments = true\ninterval = 0.25\n";
    }
    const auto output_options = qmomkit::config::load_output_options(filename);
    assert(output_options.ok());
    assert(output_options.value->directory == std::filesystem::path("output"));
    assert(output_options.value->write_moments);
    assert(output_options.value->interval.has_value());
    assert_near(*output_options.value->interval, 0.25);
    std::filesystem::remove(filename);

    const auto missing_output = qmomkit::config::load_output_options(filename);
    assert(!missing_output.ok());
    assert(has_error(missing_output, ConfigurationErrorCode::file_error));

    {
        std::ofstream output(filename);
        output << "[time_integration.node_advection]\nscheme = adaptive_rk23\n";
    }
    const auto node_advection_scheme = qmomkit::config::load_node_advection_scheme(filename);
    assert(node_advection_scheme.ok());
    assert(*node_advection_scheme.value == qmomkit::NodeAdvectionScheme::adaptive_rk23);
    std::filesystem::remove(filename);

    const auto missing_node_advection = qmomkit::config::load_node_advection_scheme(filename);
    assert(!missing_node_advection.ok());
    assert(has_error(missing_node_advection, ConfigurationErrorCode::file_error));
}

} // namespace

int main()
{
    test_valid_qmom_matches_programmatic_case();
    test_qmom_distribution_initialization_matches_direct_api();
    test_all_supported_distribution_families_parse();
    test_qmom_mixture_initialization_matches_direct_api();
    test_mixture_zero_density_component_is_skipped();
    test_valid_cqmom_matches_programmatic_case();
    test_cqmom_independent_product_initialization_matches_direct_api();
    test_lean_cqmom_initialization_matches_dense_subselection();
    test_lean_cqmom_initialization_matches_example_flow();
    test_lean_cqmom_initialization_omits_maximum_orders();
    test_dense_cqmom_initialization_requires_maximum_orders();
    test_lean_cqmom_initialization_selects_explicit_dense_tensor();
    test_cqmom_lean_case_uses_selected_analytical_moments();
    test_cqmom_lean_case_matches_dense_case_path();
    test_cqmom_lean_case_selects_explicit_dense_tensor();
    test_cqmom_case_uses_selected_analytical_moments();
    test_cqmom_case_accepts_explicit_dense_tensor();
    test_cqmom_case_initial_state_preserves_case_settings();
    test_cqmom_case_file_loading();
    test_cqmom_independent_product_zero_density_skips_marginal_validation();
    test_cqmom_adaptive_options_are_independent();
    test_cqmom_diagnostic_level_parses_values();
    test_cqmom_vandermonde_solver_parses_values();
    test_cqmom_case_and_state_report_vandermonde_solver();
    test_realizability_policy_defaults_to_disabled();
    test_realizability_policy_parses_check_modes();
    test_realizability_policy_parses_support_and_bounded_interval();
    test_realizability_policy_rejects_invalid_values();
    test_realizability_policy_rejects_invalid_bounded_interval();
    test_realizability_policy_rejects_unknown_entries();
    test_qmom_case_parsing_tolerates_realizability_policy_without_installing_it();
    test_scheme_defaults_and_parses_values();
    test_scheme_rejects_invalid_values();
    test_node_advection_scheme_defaults_and_parses_values();
    test_node_advection_scheme_rejects_invalid_values();
    test_output_options_default_to_disabled();
    test_output_options_parse_enabled_fields();
    test_output_options_reject_invalid_values();
    test_qmom_case_parsing_tolerates_scheme();
    test_qmom_case_parsing_tolerates_output_options_without_installing_them();
    test_duplicate_and_missing_entries();
    test_malformed_nonfinite_invalid_order_and_unknown_key();
    test_explicit_moments_take_precedence_over_initial_distribution();
    test_configuration_format_infers_formulation_from_phase_space();
    test_configuration_rejects_unsupported_or_ambiguous_format_fields();
    test_file_loading_and_file_failure();
    return 0;
}
